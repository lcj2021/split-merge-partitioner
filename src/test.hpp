#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "util.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

class Test : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    size_t assigned_edges;
    int p, bucket;
    double average_degree;
    size_t capacity;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    // std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores;
    // std::vector<dense_bitset> is_boundarys;

  public:
    Test(std::string basefilename);
    void split();
};
