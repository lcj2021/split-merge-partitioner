#ifndef HYBRID_PARTITIONER_HPP
#define HYBRID_PARTITIONER_HPP

#include <queue>
#include <random>
#include <unordered_map>

#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* Hybrid (powerlyra) */
class HybridPartitioner : public Partitioner
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

    // degree threshold 
    vid_t degree_threshold = 100;
    vid_t num_visit_vertices;
    eid_t num_visit_edges;

    vid_t free_vertex;
    // vertex id, dist from super

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

    void calculate_stats();

public:
    HybridPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif