#ifndef PARTITIONER_HPP
#define PARTITIONER_HPP

#include "dense_bitset.hpp"
#include "hep_graph.hpp"

class PartitionerBase
{
public:
    vid_t num_vertices;
    eid_t num_edges;
    virtual void split() = 0;        
};

template <typename TAdj>
class AdjListPartitioner: public PartitionerBase
{
protected:
    Timer total_time;

public:
    
    std::vector<eid_t> occupied;
    std::vector<dense_bitset> is_boundarys;
    std::vector<edge_t> edges;
    std::vector<vid_t> degrees;
    std::vector<bid_t> edgelist2bucket;
    mem_graph_t<TAdj> mem_graph;
};

template class AdjListPartitioner<adj_with_bid_t>;
template class AdjListPartitioner<adj_t>;

class EdgeListPartitioner: public PartitionerBase
{
protected:
    Timer total_time;

public:
    
    std::vector<eid_t> occupied;
    std::vector<dense_bitset> is_boundarys;
    std::vector<edge_t> edges;
    std::vector<vid_t> degrees;
    std::vector<bid_t> edgelist2bucket;
};

#endif