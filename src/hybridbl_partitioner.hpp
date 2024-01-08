#ifndef HYBRIDBL_PARTITIONER_HPP
#define HYBRIDBL_PARTITIONER_HPP

#include <queue>
#include <random>
#include <unordered_map>

#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* Hybrid-BL (Topology refactorization) */
class HybridBLPartitioner : public Partitioner
{
private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    eid_t assigned_edges;
    bid_t p, bucket;
    double average_degree;
    eid_t capacity;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    // std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<dense_bitset> is_cores;
    // std::vector<dense_bitset> is_boundarys;

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
    vid_t num_visit_vertices;
    eid_t num_visit_edges;

    // V[vid] = 1, if vid is handled
    dense_bitset V;
    std::vector<vid_t> free_vertex;
    // vertex id, dist from super
    std::vector<std::queue<std::tuple<vid_t, vid_t, vid_t>>> Q;

    vid_t rand_vertice_id;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    bool need_k_split;
    edgepart_writer<vid_t, bid_t> writer;

    void assign_edge(int bucket, vid_t from, vid_t to, size_t edge_id)
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
            free_vertex[machine] += p;
        }
        if (free_vertex[machine] >= num_vertices)
            return false;
        return true;
    }

    void calculate_stats();
    void init_fusion(bid_t machine, vid_t vid, vid_t dist);
    void fusion(bid_t machine, vid_t vid, vid_t root, vid_t dist);
    void fission(bid_t machine, vid_t vid);

public:
    HybridBLPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif