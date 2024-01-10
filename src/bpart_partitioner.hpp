#ifndef BPART_PARTITIONER_HPP
#define BPART_PARTITIONER_HPP

#include <random>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

template <typename TAdj>
class BPartPartitioner : public AdjListPartitioner<TAdj>
{
  private:
    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;

    bid_t p;

    using AdjListPartitioner<TAdj>::total_time;
    using AdjListPartitioner<TAdj>::num_vertices;
    using AdjListPartitioner<TAdj>::num_edges;
    using AdjListPartitioner<TAdj>::occupied;
    using AdjListPartitioner<TAdj>::is_boundarys;
    using AdjListPartitioner<TAdj>::edges;
    using AdjListPartitioner<TAdj>::degrees;
    using AdjListPartitioner<TAdj>::edgelist2bucket;

    double c_ = 0.5;
    
    std::vector<vid_t> num_bucket_vertices;
    std::vector<eid_t> num_bucket_edges;

    graph_t adj_out, adj_in;
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

    std::tuple<bid_t, eid_t> best_scored_partition(vid_t v); // returns <final bucket, additional edges> whose score is best for vertex v
    std::tuple<double, vid_t> compute_partition_score(vid_t vid, bid_t bucket_id);

    void calculate_stats();

    template <typename T>
    double jains_fairness(const std::vector<T>& L);

  public:
    BPartPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

template <typename TAdj>
template <typename T>
double BPartPartitioner<TAdj>::jains_fairness(const std::vector<T>& L)
{
    double a = 0.0;
    for (const auto& x : L) {
        a += static_cast<double>(x);
    }

    double b = 0.0;
    for (const auto& x : L) {
        b += static_cast<double>(x) * static_cast<double>(x);
    }
    b *= static_cast<double>(L.size());

    return (a * a) / b;
}

template class BPartPartitioner<adj_t>;

#endif