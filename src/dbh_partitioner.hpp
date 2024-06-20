#ifndef DBH_PARTITIONER_HPP
#define DBH_PARTITIONER_HPP

#include <memory>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "ne_graph.hpp"
#include "partitioner.hpp"

class DbhPartitioner : public EdgeListEPartitioner
{
  private:
    std::string basefilename;

    using EdgeListEPartitioner::partition_time;

    double avg_edge_cnt;

    std::unique_ptr<EdgepartWriterBase<vid_t, bid_t>> writer = nullptr;

    void assign_edge(bid_t bucket, vid_t from, vid_t to, eid_t edge_id)
    {
        writer->save_edge(from, to, bucket);
        // edgelist2bucket[edge_id] = bucket;
        ++occupied[bucket];
        is_boundarys[bucket].set_bit_unsync(from);
        is_boundarys[bucket].set_bit_unsync(to);
    }

  public:
    DbhPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif