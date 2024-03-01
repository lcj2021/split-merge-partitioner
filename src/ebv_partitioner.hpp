#ifndef EBV_PARTITIONER_HPP
#define EBV_PARTITIONER_HPP

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "ne_graph.hpp"
#include "partitioner.hpp"

class EbvPartitioner : public EdgeListEPartitioner
{
  private:
    std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;

    double avg_edge_cnt;

    std::vector<vid_t> num_bucket_vertices;

    edgepart_writer<vid_t, bid_t> writer;
    eid_t num_vertices_all_buckets;

    void assign_edge(bid_t bucket, vid_t from, vid_t to, eid_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        // edgelist2bucket[edge_id] = bucket;
        ++occupied[bucket];
        if (!is_boundarys[bucket].get(from)) {
            ++num_bucket_vertices[bucket];
            ++num_vertices_all_buckets;
            is_boundarys[bucket].set_bit_unsync(from);
        }
        if (!is_boundarys[bucket].get(to)) {
            ++num_bucket_vertices[bucket];
            ++num_vertices_all_buckets;
            is_boundarys[bucket].set_bit_unsync(to);
        }
    }

    // returns bucket id where score is best for edge (u,v)
    bid_t best_scored_partition(vid_t u, vid_t v, eid_t edge_id); 
    double compute_partition_score(vid_t u, vid_t v, bid_t bucket_id, eid_t edge_id);

  public:
    EbvPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif