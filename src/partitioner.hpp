#ifndef PARTITIONER_HPP
#define PARTITIONER_HPP

#include "dense_bitset.hpp"
#include "hep_graph.hpp"

class Partitioner
{
protected:
    Timer total_time;

public:
    std::vector<eid_t> occupied;
    std::vector<dense_bitset> is_boundarys;
    vid_t num_vertices;
    eid_t num_edges;
    std::vector<edge_t> edges;
    std::vector<eid_t> degrees;
    std::vector<bid_t> edge2bucket;
    std::vector<bid_t> edgelist2bucket;
    mem_graph_t<vid_eid_t> mem_graph;
    virtual void split() = 0;
};

#endif