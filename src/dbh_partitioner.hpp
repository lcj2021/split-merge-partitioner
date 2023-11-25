#ifndef DBH_PARTITIONER_HPP
#define DBH_PARTITIONER_HPP

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "graph.hpp"
#include "partitioner.hpp"

class DbhPartitioner : public Partitioner
{
  private:
    std::string basefilename;

    // vid_t num_vertices;
    // size_t num_edges;
    double avg_edge_cnt;
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
    size_t all_part_vertice_cnt;

    inline void assign_edge(int bucket, vid_t from, vid_t to, size_t edge_id) noexcept
    {
        writer.save_edge(from, to, bucket);
        edgelist2bucket[edge_id] = bucket;
        occupied[bucket]++;
        if (!is_boundarys[bucket].get(from)) {
            ++ vcount[bucket];
            ++ all_part_vertice_cnt;
            is_boundarys[bucket].set_bit_unsync(from);
        }
        if (!is_boundarys[bucket].get(to)) {
            ++ vcount[bucket];
            ++ all_part_vertice_cnt;
            is_boundarys[bucket].set_bit_unsync(to);
        }
    }

    void calculate_stats();

  public:
    DbhPartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif