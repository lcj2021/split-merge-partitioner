#ifndef HDRF_PARTITIONER_HPP
#define HDRF_PARTITIONER_HPP

#include <numeric>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

class HdrfPartitioner : public Partitioner
{
  private:
    std::string basefilename;

    // vid_t num_vertices;
    // size_t num_edges;
    int p;

    // use mmap for file input
    int fin;
    off_t filesize;
    char *fin_map, *fin_ptr, *fin_end;

    std::vector<vid_t> degrees;

    // std::vector<edge_t> edges;
    std::vector<size_t> vcount;
    // std::vector<dense_bitset> is_boundarys;
    edgepart_writer<vid_t, uint16_t> writer;
    int max_degree;
    size_t capacity;
    size_t min_size = 0; // currently smallest partition
    size_t max_size = 0; // currently largest partition
    double lambda = 1.10;

    void assign_edge(int bucket, vid_t from, vid_t to, size_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        edgelist2bucket[edge_id] = bucket;
        occupied[bucket]++;

        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);

        // ++ degrees[from];
        // ++ degrees[to];
    }

    int best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)
    double compute_partition_score(vid_t u, vid_t v, int bucket_id);

    void calculate_stats();

  public:
    HdrfPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif