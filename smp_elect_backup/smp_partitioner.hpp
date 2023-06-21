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
#include <queue>
#include <unordered_map>

#include "util.hpp"
#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* SplitMerge Partitioner (SMP) */
class SmpPartitioner : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p, bucket, filled_bucket_cnt;
    double average_degree;
    double theta;
    size_t capacity;

    std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> bucket_edge_cnt;
    std::vector<vid_t> degrees;
    std::vector<int16_t> edge2bucket;

    vid_t rand_vertice_id;
    
    struct bucket_info_item {
        dense_bitset is_core, is_boundary;
        size_t occupied, old_id, replicas;
        bool is_chosen;
        bucket_info_item(vid_t num_vertices) {
            is_core = dense_bitset(num_vertices);
            is_boundary = dense_bitset(num_vertices);
            old_id = occupied = replicas = 0;
            is_chosen = false;
        }
        void clear() noexcept {
            is_core.clear();
            is_boundary.clear();
            old_id = occupied = replicas = 0;
            is_chosen = false;
        }
        bool operator < (const bucket_info_item& rhs) {
            if (is_chosen != rhs.is_chosen) return is_chosen > rhs.is_chosen;
            return old_id < rhs.old_id;
        }
    };
    std::vector<bucket_info_item> bucket_info;
    dense_bitset is_in_a_core; 

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    edgepart_writer<vid_t, uint16_t> writer;
    std::array<size_t, 2> valid_threshold;
    std::queue<vid_t> q;

    void assign_edge(int bucket, vid_t from, vid_t to, size_t edge_id)
    {
        edge2bucket[edge_id] = bucket;
        bucket_info[bucket].occupied++;
        degrees[from]--;
        degrees[to]--;
    }

    void add_boundary(vid_t vid)
    {
        auto &is_core = bucket_info[bucket].is_core, &is_boundary = bucket_info[bucket].is_boundary;

        if (is_boundary.get(vid))
            return;
        is_boundary.set_bit_unsync(vid);

        if (!is_core.get(vid)) {
            min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
        }

        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size();) {
                if (edges[neighbors[i].v].valid()) {
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (is_core.get(u)) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid, neighbors[i].v);
                        min_heap.decrease_key(vid);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               bucket_info[bucket].occupied < capacity) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid, neighbors[i].v);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        i++;
                } else {
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    void occupy_vertex(vid_t vid, vid_t d)
    {
        CHECK(!bucket_info[bucket].is_core.get(vid)) << "add " << vid << " to core again";
        bucket_info[bucket].is_core.set_bit_unsync(vid);

        if (d == 0)
            return;

        add_boundary(vid);

        for (auto &i : adj_out[vid])
            if (edges[i.v].valid())
                add_boundary(edges[i.v].second);
        adj_out[vid].clear();

        for (auto &i : adj_in[vid])
            if (edges[i.v].valid())
                add_boundary(edges[i.v].first);
        adj_in[vid].clear();
    }

    void add_boundary_bfs(vid_t vid)
    {
        auto &is_core = bucket_info[bucket].is_core, &is_boundary = bucket_info[bucket].is_boundary;

        if (is_boundary.get(vid))
            return;
        is_boundary.set_bit_unsync(vid);

        if (!is_core.get(vid)) {
            q.push(vid);
        }

        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size();) {
                size_t edge_id = neighbors[i].v;
                if (edges[edge_id].valid()) {
                    vid_t &u = direction ? edges[edge_id].second : edges[edge_id].first;
                    if (is_core.get(u)) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid, edge_id);
                        min_heap.decrease_key(vid);
                        edges[edge_id].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               bucket_info[bucket].occupied < capacity) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid, edge_id);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        edges[edge_id].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else {
                        i++;
                    }
                } else {
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    void occupy_vertex_bfs(vid_t vid)
    {
        CHECK(!bucket_info[bucket].is_core.get(vid)) << "add " << vid << " to core again";
        bucket_info[bucket].is_core.set_bit_unsync(vid);

        add_boundary_bfs(vid);

        for (auto &i : adj_out[vid])
            if (edges[i.v].valid())
                add_boundary_bfs(edges[i.v].second);
        adj_out[vid].clear();

        for (auto &i : adj_in[vid])
            if (edges[i.v].valid())
                add_boundary_bfs(edges[i.v].first);
        adj_in[vid].clear();
    }

    bool get_free_vertex(vid_t &vid)
    {
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                    2 * average_degree ||
                bucket_info[bucket].is_core.get(vid))) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    bool get_free_vertex_by_rand(vid_t &vid)
    {
        rand_vertice_id = dis(gen);
        vid = rand_vertice_id;
        for (size_t offset = 0; offset < num_vertices; ++ offset) {
            vid = (vid + 1) % num_vertices;
            if (bucket_info[bucket].is_core.get(vid)) continue;
            // if (is_in_a_core.get(vid)) continue;
            if (adj_out[vid].size() + adj_in[vid].size() == 0) continue;
            rand_vertice_id = vid;
            return true;
        }
        return false;
    }

    size_t rearrange_edge(std::vector<edge_t> &e, const std::unordered_map<int, int> &valid_bucket, bool is_final_round)
    {
        size_t curr_assigned_edges = 0;
        for (size_t edge_id = 0; edge_id < e.size(); ++ edge_id) {
            int16_t &edge_bucket = edge2bucket[edge_id];
            if (edge_bucket < filled_bucket_cnt) continue;
            e[edge_id].recover();
            if (valid_bucket.count(edge_bucket)) {
                edge_bucket = valid_bucket.at(edge_bucket);
                bucket_edge_cnt[edge_bucket] ++;
                e[edge_id].remove();
                ++ curr_assigned_edges;
            } else {        
                if (is_final_round) {   // [[unlikely]] No edge should be left in the final round
                    std::cerr << "Edge(id): " << edge_id << ", bucket: " << edge_bucket << ", should not be left in the final round\n";
                    exit(1);
                } else {
                    edge_bucket = -1;
                }
            }
        }
        return curr_assigned_edges;
    }

    bool check_edge()
    {
        std::vector<dense_bitset> dbitsets(FLAGS_p, dense_bitset(num_vertices));

        for (size_t edge_id = 0; edge_id < num_edges; ++ edge_id) {
            edges[edge_id].recover();
            int16_t edge_bucket = edge2bucket[edge_id];
            vid_t u = edges[edge_id].first, v = edges[edge_id].second;
            dbitsets[edge_bucket].set_bit_unsync(u);
            dbitsets[edge_bucket].set_bit_unsync(v);
        }
        auto equal_dbitset = [&](const dense_bitset &l, const dense_bitset &r) {
            if (l.size() != r.size()) return false;
            for (size_t bit = 0; bit < l.size(); ++ bit) {
                if (l.get(bit) != r.get(bit)) { return false; }
            }   
            return true;
        };
        for (int b = 0; b < FLAGS_p; ++ b) {
            if (!equal_dbitset(dbitsets[b], bucket_info[b].is_boundary)) {
                return false;
            }
        }
        return true;
    }

    void expand()
    {
        std::cerr << bucket << ", ";
        DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (bucket_info[bucket].occupied < capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) {
                if (!get_free_vertex_by_rand(vid)) {
                    DLOG(INFO) << "partition " << bucket
                               << " stop: no free vertices";
                    break;
                }
                d = adj_out[vid].size() + adj_in[vid].size();
            } else {
                min_heap.remove(vid);
                /* CHECK_EQ(d, adj_out[vid].size() + adj_in[vid].size()); */
            }

            occupy_vertex(vid, d);
        }
        min_heap.clear();
        rep (direction, 2)
            repv (vid, num_vertices) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }

    std::unordered_map<int, int> mp;
    int curr_bucket_id;
    int get_final_bucket(int bucket_id)
    {
        if (mp.count(bucket_id)) {
            return mp[bucket_id];
        } else {
            return mp[bucket_id] = curr_bucket_id ++;
        }
    }

    void assign_remaining(int last_p);
    std::array<size_t, 2> count_mirrors_elect(int iter);   // return {mirrors, filled_bucket_this_iter}
    void count_mirrors_k_split();
    void calculate_stats();
    size_t elect(int iter);     // return the number of filled buckets
    void k_split();

    int merge_bucket(int dst, int src, bool &has_intersection);
    std::unordered_map<int, int> merge_by_size();
    std::unordered_map<int, int> merge_by_overlap();

  public:
    SmpPartitioner(std::string basefilename);
    void split();
};
