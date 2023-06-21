#pragma once

#include <vector>
#include "util.hpp"
#include "dense_bitset.hpp"

class Partitioner
{
  protected:
    Timer total_time;

  public:
    std::vector<size_t> occupied;
    std::vector<dense_bitset> is_boundarys; 
    vid_t num_vertices;
    size_t num_edges;
    std::vector<edge_t> edges;
    std::vector<int16_t> edge2bucket;
    virtual void split() = 0;
};
