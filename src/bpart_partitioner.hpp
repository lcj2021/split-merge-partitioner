#ifndef BPART_PARTITIONER_HPP
#define BPART_PARTITIONER_HPP

#include <random>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "ne_graph.hpp"
#include "partitioner.hpp"

template <typename TAdj>
class BPartPartitioner : public AdjListVPartitioner<TAdj>
{
  private:
    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;

    using AdjListVPartitioner<TAdj>::total_time;
    using AdjListVPartitioner<TAdj>::partition_time;
    using AdjListVPartitioner<TAdj>::num_vertices;
    using AdjListVPartitioner<TAdj>::num_edges;
    using AdjListVPartitioner<TAdj>::num_partitions;
    using AdjListVPartitioner<TAdj>::occupied;
    using AdjListVPartitioner<TAdj>::is_boundarys;
    using AdjListVPartitioner<TAdj>::edges;
    using AdjListVPartitioner<TAdj>::degrees;
    using AdjListVPartitioner<TAdj>::calculate_stats;

    double c_ = 0.5;
    
    std::vector<vid_t> num_bucket_vertices;

    Graph<AdjEntryVid> graph;
    std::vector<double> w_;
    vertexpart_writer<vid_t, bid_t> writer;
    std::vector<bid_t> vertex2bucket;
    double average_degree;
    double capacity;

    // Parameters
    /// @ref: https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L760
    double gamma = 1.5;
    double alpha = 1.5;

    /// https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L1131
    double balance_factor = 0.05;

    double intra_partition_cost(double size);

    std::tuple<vid_t, vid_t> overlap_partition_vertex(vid_t vid, bid_t bucket_id);

    /// @ref: https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L858
    void assign_vertex(bid_t bucket, vid_t vid, eid_t additional_edges);
    
    /// @return: <final bucket, additional edges> whose score is best for vertex v
    std::tuple<bid_t, eid_t> best_scored_partition(vid_t v); 
    std::tuple<double, vid_t> compute_partition_score(vid_t vid, bid_t bucket_id);

  public:
    BPartPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

template class BPartPartitioner<adj_t>;

#endif