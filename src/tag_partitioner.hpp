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
#include "dsu.hpp"

class TagPartitioner : public Partitioner
{
private:
    /* data */
    const double BALANCE_RATIO = 1.05;
    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p;
    size_t threshold;
    size_t capacity;
    double average_degree;

    std::vector<vid_t> degree_1_vertex;
    std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    // MinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> occupied;
    std::vector<size_t> global_tag_distribute;
    std::vector<vid_t> degrees;             // vertex_id [0, n - 1]
    std::unordered_set<vid_t> bfs_src;
    std::vector<dense_bitset> vertex2tag;
    std::vector<bool> can_cover;
    std::vector<bool> edge_covered;
    std::vector<bool> tag_valid;

    std::vector<std::unordered_set<vid_t>> bucket;


    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    edgepart_writer<vid_t, uint16_t> writer;

    DSU dsu;
    int kcore_k;
    std::unordered_set<int> kcore_s;


    void random_tag(size_t random_cnt);
    void bfs_walk(size_t random_cnt);
    bool seed_check(vid_t uid);
    int choose_tag(vid_t uid, bool restrict);
    int choose_tag(vid_t uid, bool restrict, bool debug);
    void assign_tag(vid_t uid, int candidate_tag, bool restrict);
    void union_tag();

    std::vector<bool> seeded;
    bool get_free_vertex(vid_t &vid)
    {
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() == 0 ||
                adj_out[vid].size() >
                    2 * average_degree ||
                seeded[vid])) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    bool get_free_vertex(vid_t &vid, int k)
    {
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() == 0 ||
                adj_out[vid].size() <=
                    k * average_degree ||
                seeded[vid])) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    bool get_free_vertex_kcore(vid_t &vid)
    {
        if (!kcore_s.size()) {
            vid = dis(gen);
            vid_t count = 0;
            while (count < num_vertices &&
                (adj_out[vid].size() == 0 ||
                    adj_out[vid].size() >
                        2 * average_degree ||
                    seeded[vid])) {
                vid = (vid + ++count) % num_vertices;
            }
            if (count == num_vertices)
                return false;
            return true;
        }
        vid = *(kcore_s.begin());
        kcore_s.erase(kcore_s.find(vid));
        return true;
    }
    
    uint64_t opposite(uint64_t edge_id) 
    {
        if (edge_id < num_edges / 2)    return edge_id + num_edges / 2;
        else                            return edge_id - num_edges / 2;
    }

    bool all_neighbor_covered(vid_t vid)
    {
        bool all_neighbor_covered = true;
        for (auto &i : adj_out[vid]) {
            if (edge_covered[i.v] || edge_covered[opposite(i.v)])  continue;

            vid_t nid = edges[i.v].second;
            bool n_covered = false;
            if (degrees[nid] == 1) {
                n_covered = true;
                continue;
            }
            for (int b = 0; b < p; ++ b) {
                if (vertex2tag[nid].get(b) == 1 && vertex2tag[vid].get(b) == 1) {
                    n_covered = true;
                    edge_covered[i.v] = edge_covered[opposite(i.v)] = true;
                    ++ occupied[b];
                }
            }
            all_neighbor_covered &= n_covered;
        }
        return all_neighbor_covered;
    }

    size_t get_uncovered_edge_cnt()
    {
        size_t cnt = 0;
        for (size_t eid = 0; eid < num_edges / 2; ++ eid) {
            bool covered = false;
            vid_t u = edges[eid].first, v = edges[eid].second;
            for (int b = 0; b < p; ++ b) {
                if (vertex2tag[u].get(b) == 1 && vertex2tag[v].get(b) == 1) {
                    covered = true;
                    break;
                }
            }
            if (!covered) ++ cnt;
        }
        return cnt;
    }

public:
    TagPartitioner(std::string basefilename);
    void split();
    ~TagPartitioner();
    std::ofstream data;
};
