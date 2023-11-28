#ifndef FSM_PARTITIONER_HPP
#define FSM_PARTITIONER_HPP

#include <random>
#include <memory>
#include <unordered_map>

#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"
#include "hep_graph.hpp"

/* Fine-grained SplitMerge Partitioner (FSM) */
class FsmPartitioner : public Partitioner
{
  
private:
    std::unique_ptr<Partitioner> split_partitioner;
    std::string split_method;

    std::string basefilename;

    size_t assigned_edges;
    bid_t p, k;

    // mem_graph_t<vid_eid_t> mem_graph;

    std::vector<eid_t> num_bucket_edges;

    struct BucketInfo {
        dense_bitset is_mirror;
        size_t occupied, old_id, replicas;
        bool is_chosen{false};
        BucketInfo(vid_t num_vertices) {
            is_mirror = dense_bitset(num_vertices);
            old_id = occupied = replicas = 0;
        }
        bool operator < (const BucketInfo& rhs) const {
            if (is_chosen != rhs.is_chosen) return is_chosen > rhs.is_chosen;
            return old_id < rhs.old_id;
        }
    };
    std::vector<BucketInfo> bucket_info;

    edgepart_writer<vid_t, bid_t> writer;

    template<typename EdgeHandler>
    eid_t iterate_edges(EdgeHandler edge_handler)
    {
        LOG(INFO) << "In-memory iterating edges (adjlist in memory)";
        eid_t num_adjlist_edges = 0;
        eid_t offset = 0;
        for (vid_t from = 0; from < num_vertices; ++from) {
            if (degrees[from] > mem_graph.high_degree_threshold) continue;
            for (vid_t i = 0; i < degrees[from]; ++i) {
                auto& to = mem_graph.neighbors[offset + i];
                auto& edge_bucket = reinterpret_cast<bid_t&>(to.bid);
                if (edge_bucket == kInvalidBid) {
                    continue;
                }
                ++num_adjlist_edges;
                edge_handler(edge_bucket, from, to.vid);
            }
            offset += degrees[from];
        }
        LOG(INFO) << "num_adjlist_edges: " << num_adjlist_edges;

        LOG(INFO) << "Streamingly iterating edges (edgelist on disk)";
        eid_t num_edgelist_edges = 0;

        mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::in);
        mem_graph.h2h_file.seekg(0, std::ios::beg);

        std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
        eid_t chunk_size;
        eid_t left_h2h_edges = mem_graph.num_h2h_edges;
        eid_t id_h2h_edges = 0;

        if (left_h2h_edges >= 100000) {
            chunk_size = 100000; // batch read of so many edges
        } else {
            chunk_size = left_h2h_edges;
        }
        stream_edges.resize(chunk_size);

        while (left_h2h_edges > 0) { // edges to be read
            mem_graph.h2h_file.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
            for (eid_t i = 0; i < chunk_size; ++i) {
                auto &edge_bucket = edgelist2bucket[id_h2h_edges++];
                ++num_edgelist_edges;
                edge_handler(edge_bucket, stream_edges[i].first, stream_edges[i].second);
            }

            left_h2h_edges -= chunk_size;
            if (left_h2h_edges < chunk_size) { // adapt chunk size for last batch read
                chunk_size = left_h2h_edges;
            }
        }
        LOG(INFO) << "num_edgelist_edges: " << num_edgelist_edges;

        mem_graph.h2h_file.close();
        return num_adjlist_edges + num_edgelist_edges;
    }
    
    eid_t rearrange_edge(std::vector<edge_t> &e, const std::unordered_map<bid_t, bid_t> &valid_bucket)
    {
        eid_t curr_assigned_edges = 0;
        for (eid_t edge_id = 0; edge_id < e.size(); ++edge_id) {
            auto &edge_bucket = edgelist2bucket[edge_id];
            e[edge_id].recover();
            if (valid_bucket.count(edge_bucket)) {
                edge_bucket = valid_bucket.at(edge_bucket);
                ++num_bucket_edges[edge_bucket] ;
                // e[edge_id].remove();
                ++curr_assigned_edges;
            } else {        // [[unlikely]]  No edge should be left
                LOG(FATAL) << "Edge(id): " << edge_id << ", bucket: " << edge_bucket << ", should not be left in the final round\n";
            }
        }
        return curr_assigned_edges;
    }

    size_t rearrange_edge_hybrid(const std::unordered_map<bid_t, bid_t> &valid_bucket)
    {
        return iterate_edges( 
            [this, &valid_bucket](bid_t& edge_bucket, 
                const vid_t& u, const vid_t& v) {
                if (valid_bucket.count(edge_bucket)) {
                    edge_bucket = valid_bucket.at(edge_bucket);
                    ++num_bucket_edges[edge_bucket];
                } else {
                    LOG(FATAL) << "bucket: " << edge_bucket 
                            << ", should not be left in the final round\n";
                }
            }
        );
    }

    bool check_edge_hybrid()
    {
        std::vector<dense_bitset> dbitsets(p, dense_bitset(num_vertices));
        for (vid_t vid = 0; vid < num_vertices; ++vid) {
            bool assigned_to_a_part = false;
            for (bid_t b = 0; b < p; ++b) {
                if (bucket_info[b].is_mirror.get(vid)) {
                    assigned_to_a_part = true;
                    break;
                }
            }
            if (!assigned_to_a_part) {
                return false;
            }
        }
        iterate_edges(
            [&dbitsets](bid_t& edge_bucket, const vid_t& u, const vid_t& v) {
                dbitsets[edge_bucket].set_bit_unsync(u);
                dbitsets[edge_bucket].set_bit_unsync(v);
            }
        );
        // for (size_t edge_id = 0; edge_id < num_edges; ++edge_id) {
        //     // edges[edge_id].recover();
        //     int16_t edge_bucket = edgelist2bucket[edge_id];
        //     vid_t u = edges[edge_id].first, v = edges[edge_id].second;
        //     dbitsets[edge_bucket].set_bit_unsync(u);
        //     dbitsets[edge_bucket].set_bit_unsync(v);
        // }
        auto equal_dbitset = [&](const dense_bitset &l, const dense_bitset &r) {
            if (l.size() != r.size()) return false;
            for (size_t bit = 0; bit < l.size(); ++bit) {
                if (l.get(bit) != r.get(bit)) { return false; }
            }   
            return true;
        };
        for (bid_t b = 0; b < p; ++b) {
            if (!equal_dbitset(dbitsets[b], bucket_info[b].is_mirror)) { return false; }
        }
        return true;
    }

    void save_edge_hybrid()
    {
        iterate_edges([this](bid_t& edge_bucket, const vid_t& u, const vid_t& v) {
            writer.save_edge(u, v, edge_bucket);
        });
    }

    bool check_edge()
    {
        std::vector<dense_bitset> dbitsets(p, dense_bitset(num_vertices));
        for (vid_t vid = 0; vid < num_vertices; ++vid) {
            bool assigned_to_a_part = false;
            for (bid_t b = 0; b < p; ++b) {
                if (bucket_info[b].is_mirror.get(vid)) {
                    assigned_to_a_part = true;
                    break;
                }
            }
            if (!assigned_to_a_part) {
                return false;
            }
        }
        for (size_t edge_id = 0; edge_id < num_edges; ++edge_id) {
            // edges[edge_id].recover();
            int16_t edge_bucket = edgelist2bucket[edge_id];
            vid_t u = edges[edge_id].first, v = edges[edge_id].second;
            dbitsets[edge_bucket].set_bit_unsync(u);
            dbitsets[edge_bucket].set_bit_unsync(v);
        }
        auto equal_dbitset = [&](const dense_bitset &l, const dense_bitset &r) {
            if (l.size() != r.size()) return false;
            for (size_t bit = 0; bit < l.size(); ++bit) {
                if (l.get(bit) != r.get(bit)) { return false; }
            }   
            return true;
        };
        for (bid_t b = 0; b < p; ++b) {
            if (!equal_dbitset(dbitsets[b], bucket_info[b].is_mirror)) { return false; }
        }
        return true;
    }

    std::unordered_map<bid_t, bid_t> mp;
    bid_t curr_bucket_id;
    bid_t get_final_bucket(bid_t bucket_id)
    {
        if (mp.count(bucket_id)) {
            return mp[bucket_id];
        } else {
            return mp[bucket_id] = curr_bucket_id++;
        }
    }

    void merge();
    void calculate_stats();

    vid_t merge_bucket(bid_t dst, bid_t src, bool &has_intersection);
    std::unordered_map<bid_t, bid_t> fast_merge();
    std::unordered_map<bid_t, bid_t> precise_merge();

public:
    FsmPartitioner(std::string basefilename);
    void split();
};

#endif