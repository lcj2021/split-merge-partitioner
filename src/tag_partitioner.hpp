#pragma once

#include <random>

#include "util.hpp"
#include "partitioner.hpp"
#include "dense_bitset.hpp"
#include "graph.hpp"
#include "edgepart.hpp"

class TagPartitioner : public Partitioner
{
private:
    /* data */
    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p;

    std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    // MinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores, is_boundarys;


    std::random_device rd;
    std::mt19937 gen;

    edgepart_writer<vid_t, uint16_t> writer;

public:
    TagPartitioner(std::string basefilename);
    void split();
    ~TagPartitioner();
};
