#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

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
    std::vector<size_t> global_tag_distribute;
    std::vector<vid_t> degrees;
    std::vector<vid_t> bfs_src;
    std::vector<int8_t> master;
    std::vector<dense_bitset> vertex2tag;
    std::vector<bool> can_cover;

    std::random_device rd;
    std::mt19937 gen;

    edgepart_writer<vid_t, uint16_t> writer;

    void random_tag(size_t random_cnt);
    void bfs_walk(size_t random_cnt);
    bool seed_check(vid_t uid);

public:
    TagPartitioner(std::string basefilename);
    void split();
    ~TagPartitioner();
};
