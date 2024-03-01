#ifndef FENNEL_PARTITIONER_HPP
#define FENNEL_PARTITIONER_HPP

#include <random>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "ne_graph.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

template <typename TAdj>
class FennelPartitioner : public AdjListVPartitioner<TAdj>
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

    // // use mmap for file input
    // int fin;
    // off_t filesize;
    // char *fin_map, *fin_ptr, *fin_end;

    Graph<AdjEntryVid> graph;

    std::vector<eid_t> w_;
    vertexpart_writer<vid_t, bid_t> writer;
    
    // std::vector<bid_t> vertex2bucket;
    eid_t capacity;
    // Parameters
    double gamma = 1.5;
    double alpha;

    double intra_partition_cost(double sz) {
        return alpha * gamma * pow(sz, gamma - 1.0);
    }

    std::tuple<vid_t, vid_t> overlap_partition_vertex(vid_t vid, bid_t bucket_id);

    void assign_vertex(bid_t bucket, vid_t vid, eid_t additional_edges)
    {
        writer.save_vertex(vid, bucket);
        is_boundarys[bucket].set_bit_unsync(vid);
        // vertex2bucket[vid] = bucket;
        // w_[bucket] += degrees[vid];
        w_[bucket] += 1;
        occupied[bucket] += additional_edges;
    }

    /// @return <final bucket, additional edges> whose score is best for vertex v
    std::tuple<bid_t, eid_t> best_scored_partition(vid_t v); 
    std::tuple<double, eid_t> compute_partition_score(vid_t vid, bid_t bucket_id);

  public:
    FennelPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

template class FennelPartitioner<adj_t>;

#endif