#ifndef HDRF_PARTITIONER_HPP
#define HDRF_PARTITIONER_HPP

#include <numeric>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

class HdrfPartitioner : public EdgeListPartitioner
{
  private:
    std::string basefilename;

    bid_t p;

    // use mmap for file input
    int fin;
    off_t filesize;
    char *fin_map, *fin_ptr, *fin_end;

    std::vector<vid_t> degrees;

    edgepart_writer<vid_t, bid_t> writer;
    vid_t max_degree;
    eid_t capacity;
    eid_t min_size = 0; // currently smallest partition
    eid_t max_size = 0; // currently largest partition
    double lambda = 1.10;

    void assign_edge(bid_t bucket, vid_t from, vid_t to, eid_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        // edgelist2bucket[edge_id] = bucket;
        ++occupied[bucket];

        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);
    }

    bid_t best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)
    double compute_partition_score(vid_t u, vid_t v, bid_t bucket_id);

    void calculate_stats();

  public:
    HdrfPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif