#ifndef HYBRID_PARTITIONER_HPP
#define HYBRID_PARTITIONER_HPP

#include <random>
#include <unordered_map>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* Hybrid (powerlyra) */
class HybridPartitioner : public EdgeListEPartitioner
{
private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    eid_t assigned_edges;
    double average_degree;

    /// @ref powerlyra EuroSys2015: degree threshold = 100
    vid_t degree_threshold = 100;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    bool need_k_split;
    edgepart_writer<vid_t, bid_t> writer;

    void assign_edge(bid_t bucket, vid_t from, vid_t to, eid_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        // CHECK_EQ(edgelist2bucket[edge_id], kInvalidBid);
        // edgelist2bucket[edge_id] = bucket;
        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);
        ++assigned_edges;
        ++occupied[bucket];
    }

public:
    HybridPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif