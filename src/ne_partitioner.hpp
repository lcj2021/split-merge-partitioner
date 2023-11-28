#ifndef NE_PARTITIONER_HPP
#define NE_PARTITIONER_HPP

#include <random>

#include "min_heap.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "graph.hpp"

/* Neighbor Expansion (NE) */
class NePartitioner : public Partitioner
{
private:
    const double BALANCE_RATIO = 1.00;

    std::string basefilename;

    eid_t assigned_edges;
    bid_t p, bucket;
    double average_degree;
    eid_t capacity;

    // std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    // std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<dense_bitset> is_cores;
    // std::vector<dense_bitset> is_boundarys;

    vid_t rand_vertice_id;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    bool need_k_split;
    edgepart_writer<vid_t, bid_t> writer;

    int check_edge(const edge_t *e)
    {
        for (bid_t b = 0; b < bucket; ++b) {
            auto &is_boundary = is_boundarys[b];
            if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
                occupied[b] < capacity) {
                return b;
            }
        }

        for (bid_t b = 0; b < bucket; ++b) {
            auto &is_core = is_cores[b], &is_boundary = is_boundarys[b];
            if ((is_core.get(e->first) || is_core.get(e->second)) &&
                occupied[b] < capacity) {
                if (is_core.get(e->first) && degrees[e->second] > average_degree)
                    continue;
                if (is_core.get(e->second) && degrees[e->first] > average_degree)
                    continue;
                is_boundary.set_bit(e->first);
                is_boundary.set_bit(e->second);
                return b;
            }
        }

        return p;
    }

    void assign_edge(int bucket, vid_t from, vid_t to, size_t edge_id)
    {
        writer.save_edge(from, to, bucket);
        CHECK_EQ(edgelist2bucket[edge_id], kInvalidBid);
        edgelist2bucket[edge_id] = bucket;
        ++assigned_edges;
        ++occupied[bucket];
        --degrees[from];
        --degrees[to];
    }

    void add_boundary(vid_t vid)
    {
        auto &is_core = is_cores[bucket], &is_boundary = is_boundarys[bucket];

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
                               occupied[bucket] < capacity) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid, neighbors[i].v);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        ++i;
                } else {
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    void occupy_vertex(vid_t vid, vid_t d)
    {
        CHECK(!is_cores[bucket].get(vid)) << "add " << vid << " to core again";
        is_cores[bucket].set_bit_unsync(vid);

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

    bool get_free_vertex_by_rand(vid_t &vid)
    {
        rand_vertice_id = dis(gen);
        vid = rand_vertice_id;
        for (vid_t offset = 0; offset < num_vertices; ++ offset) {
            vid = (vid + 1) % num_vertices;
            if (is_cores[bucket].get(vid)) continue;
            // if (is_in_a_core.get(vid)) continue;
            if (adj_out[vid].size() + adj_in[vid].size() == 0) continue;
            rand_vertice_id = vid;
            return true;
        }
        return false;
    }

    bool get_free_vertex(vid_t &vid)
    {
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                    2 * average_degree ||
                is_cores[bucket].get(vid))) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    void assign_remaining();
    void calculate_stats();

public:
    NePartitioner(std::string basefilename, bool need_k_split);
    void split();
};

#endif