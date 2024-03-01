#ifndef HYBRIDBL_PARTITIONER_HPP
#define HYBRIDBL_PARTITIONER_HPP

#include <queue>
#include <random>
#include <unordered_map>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "ne_graph.hpp"
#include "graph.hpp"

/* Hybrid-BL (Topology refactorization) */
template <typename TAdj>
class HybridBLPartitioner : public AdjListEPartitioner<TAdj>
{
private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    eid_t assigned_edges;
    double average_degree;
    eid_t capacity;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;

    std::vector<eid_t> fission_occupy;
    std::vector<eid_t> fusion_occupy;
    std::unordered_map<vid_t, eid_t> root_assigned;
    std::unordered_map<vid_t, bid_t> root_bucket;

    // degree threshold 
    vid_t degree_threshold = 100;
    // radius threshold
    vid_t gamma = 3;
    eid_t fusion_threshold = 1'000'000;
    std::vector<vid_t> super;

    // V[vid] = 1, if vid is handled
    dense_bitset V;
    std::vector<vid_t> free_vertex;
    // vertex id, dist from super
    std::vector<std::queue<std::tuple<vid_t, vid_t, vid_t>>> Q;

    using AdjListEPartitioner<TAdj>::total_time;
    using AdjListEPartitioner<TAdj>::num_vertices;
    using AdjListEPartitioner<TAdj>::num_edges;
    using AdjListEPartitioner<TAdj>::num_partitions;
    
    using AdjListEPartitioner<TAdj>::occupied;
    using AdjListEPartitioner<TAdj>::is_boundarys;
    using AdjListEPartitioner<TAdj>::edges;
    using AdjListEPartitioner<TAdj>::degrees;
    using AdjListEPartitioner<TAdj>::edgelist2bucket;
    using EdgePartitioner::calculate_stats;

    bool need_k_split;
    edgepart_writer<vid_t, bid_t> writer;

    void assign_edge(bid_t bucket, vid_t from, vid_t to, eid_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        CHECK_EQ(edgelist2bucket[edge_id], kInvalidBid);
        edgelist2bucket[edge_id] = bucket;
        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);
        ++assigned_edges;
        ++occupied[bucket];
    }

    // bool get_free_vertex_by_rand()
    // {
    //     free_vertex = dis(gen);
    //     vid_t vid = free_vertex;
    //     for (vid_t offset = 0; offset < num_vertices; ++offset) {
    //         vid = (free_vertex + offset) % num_vertices;
    //         if (V.get(vid) == 1) continue;
    //         free_vertex = vid;
    //         return true;
    //     }
    //     return false;
    // }

    bool get_free_vertex(bid_t machine)
    {
        while (free_vertex[machine] < num_vertices
            && V.get(free_vertex[machine]) == 1) {
            free_vertex[machine] += num_partitions;
        }
        if (free_vertex[machine] >= num_vertices)
            return false;
        return true;
    }

    void init_fusion(bid_t machine, vid_t vid, vid_t dist);
    void fusion(bid_t machine, vid_t vid, vid_t root, vid_t dist);
    void fission(bid_t machine, vid_t vid);

public:
    HybridBLPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

template class HybridBLPartitioner<adj_t>;

#endif