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

/* Fine-grained SplitMerge Partitioner (FSM) */
class FsmPartitioner : public Partitioner
{
  
  private:
    std::unique_ptr<Partitioner> split_partitioner;

    std::string basefilename;

    size_t assigned_edges;
    int p, k;

    std::vector<size_t> bucket_edge_cnt;

    struct BucketInfo {
        dense_bitset is_mirror;
        size_t occupied, old_id, replicas;
        bool is_chosen;
        BucketInfo(vid_t num_vertices) {
            is_mirror = dense_bitset(num_vertices);
            old_id = occupied = replicas = 0;
            is_chosen = false;
        }
        bool operator < (const BucketInfo& rhs) const {
            if (is_chosen != rhs.is_chosen) return is_chosen > rhs.is_chosen;
            return old_id < rhs.old_id;
        }
    };
    std::vector<BucketInfo> bucket_info;

    edgepart_writer<vid_t, uint16_t> writer;
    
    size_t rearrange_edge(std::vector<edge_t> &e, const std::unordered_map<int, int> &valid_bucket)
    {
        size_t curr_assigned_edges = 0;
        for (size_t edge_id = 0; edge_id < e.size(); ++ edge_id) {
            auto &edge_bucket = edge2bucket[edge_id];
            e[edge_id].recover();
            if (valid_bucket.count(edge_bucket)) {
                edge_bucket = valid_bucket.at(edge_bucket);
                bucket_edge_cnt[edge_bucket] ++;
                // e[edge_id].remove();
                ++ curr_assigned_edges;
            } else {        // [[unlikely]]  No edge should be left
                LOG(FATAL) << "Edge(id): " << edge_id << ", bucket: " << edge_bucket << ", should not be left in the final round\n";
            }
        }
        return curr_assigned_edges;
    }

    bool check_edge()
    {
        std::vector<dense_bitset> dbitsets(p, dense_bitset(num_vertices));
        for (vid_t vid = 0; vid < num_vertices; ++ vid) {
            bool assigned_to_a_part = false;
            for (int b = 0; b < p; ++ b) {
                if (bucket_info[b].is_mirror.get(vid)) {
                    assigned_to_a_part = true;
                    break;
                }
            }
            if (!assigned_to_a_part) {
                return false;
            }
        }
        for (size_t edge_id = 0; edge_id < num_edges; ++ edge_id) {
            // edges[edge_id].recover();
            int16_t edge_bucket = edge2bucket[edge_id];
            vid_t u = edges[edge_id].first, v = edges[edge_id].second;
            dbitsets[edge_bucket].set_bit_unsync(u);
            dbitsets[edge_bucket].set_bit_unsync(v);
        }
        auto equal_dbitset = [&](const dense_bitset &l, const dense_bitset &r) {
            if (l.size() != r.size()) return false;
            for (size_t bit = 0; bit < l.size(); ++ bit) {
                if (l.get(bit) != r.get(bit)) { return false; }
            }   
            return true;
        };
        for (int b = 0; b < p; ++ b) {
            if (!equal_dbitset(dbitsets[b], bucket_info[b].is_mirror)) { return false; }
        }
        return true;
    }

    std::unordered_map<int, int> mp;
    int curr_bucket_id;
    int get_final_bucket(int bucket_id)
    {
        if (mp.count(bucket_id)) {
            return mp[bucket_id];
        } else {
            return mp[bucket_id] = curr_bucket_id ++;
        }
    }

    void merge();
    void calculate_stats();

    int merge_bucket(int dst, int src, bool &has_intersection);
    std::unordered_map<int, int> fast_merge();
    std::unordered_map<int, int> precise_merge();

  public:
    FsmPartitioner(std::string basefilename);
    void split();
};

#endif