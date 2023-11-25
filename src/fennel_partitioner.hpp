#ifndef FENNEL_PARTITIONER_HPP
#define FENNEL_PARTITIONER_HPP

#include <random>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

class FennelPartitioner : public Partitioner
{
  private:
    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;
    // vid_t num_vertices;
    // size_t num_edges;
    bid_t p;

    // // use mmap for file input
    // int fin;
    // off_t filesize;
    // char *fin_map, *fin_ptr, *fin_end;

    std::vector<vid_t> degrees;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    std::vector<size_t> vcount;
    // std::vector<dense_bitset> is_boundarys;
    vertexpart_writer<vid_t, bid_t> writer;
    std::vector<uint16_t> vertex2bucket;
    int max_degree;
    size_t capacity;
    // Parameters
    double gamma = 1.5;
    double alpha;

    double intra_partition_cost(double sz) {
        return alpha * gamma * pow(sz, gamma - 1.0);
    }

    std::tuple<size_t, size_t> overlap_partition_vertex(vid_t vid, int bucket_id);

    void assign_vertex(int bucket, vid_t vid, size_t additional_edges)
    {
        writer.save_vertex(vid, bucket);
        is_boundarys[bucket].set_bit_unsync(vid);
        vertex2bucket[vid] = bucket;
        // vcount[bucket] += degrees[vid];
        vcount[bucket] += 1;
        occupied[bucket] += additional_edges;
    }

    std::tuple<bid_t, size_t> best_scored_partition(vid_t v); // returns <final bucket, additional edges> whose score is best for vertex v
    std::tuple<double, size_t> compute_partition_score(vid_t vid, int bucket_id);

    void calculate_stats();

  public:
    FennelPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif