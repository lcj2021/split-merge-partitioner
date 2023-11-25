#ifndef BPART_PARTITIONER_HPP
#define BPART_PARTITIONER_HPP

#include <random>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

class BPartPartitioner : public Partitioner
{
  private:
    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;
    // vid_t num_vertices;
    // size_t num_edges;
    int p;

    double c_ = 0.5;
    
    std::vector<vid_t> degrees;
    std::vector<vid_t> num_vertices_bucket;
    std::vector<size_t> num_edges_bucket;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    std::vector<double> w_;
    // std::vector<dense_bitset> is_boundarys;
    vertexpart_writer<vid_t, bid_t> writer;
    std::vector<uint16_t> vertex2bucket;
    int max_degree;
    double average_degree;
    size_t capacity;

    // Parameters
    /// @ref: https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L760
    double gamma = 1.5;
    double alpha = 1.5;

    /// https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L1131
    double balance_factor = 0.05;

    double intra_partition_cost(double sz);

    std::tuple<size_t, size_t> overlap_partition_vertex(vid_t vid, int bucket_id);

    /// @ref: https://github.com/ustcadsl/BPart/blob/master/include/graph.hpp#L858
    void assign_vertex(int bucket, vid_t vid, size_t additional_edges);

    std::tuple<int, size_t> best_scored_partition(vid_t v); // returns <final bucket, additional edges> whose score is best for vertex v
    std::tuple<double, size_t> compute_partition_score(vid_t vid, int bucket_id);

    void calculate_stats();

    template <typename T>
    double jains_fairness(const std::vector<T>& L);

  public:
    BPartPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif