#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "common.hpp"

struct uint40_t {
        uint64_t v:40;
} __attribute__((packed));

class adjlist_t
{
private:
    uint40_t *adj;
    vid_t len;

public:
    adjlist_t() : adj(NULL), len(0) {}
    adjlist_t(uint40_t *adj, vid_t len = 0) : adj(adj), len(len) {}
    uint40_t *begin() { return adj; }
    uint40_t *end() { return adj + len; }
    void increment() { len++; }
    void push_back(size_t data) { adj[len++].v = data; }
    size_t size() const { return len; }
    uint40_t &operator[](size_t idx) { return adj[idx]; };
    const uint40_t &operator[](size_t idx) const { return adj[idx]; };
    uint40_t &back() { return adj[len - 1]; };
    const uint40_t &back() const { return adj[len - 1]; };
    void pop_back() { len--; }
    void clear() { len = 0; }
};

class graph_t
{
  private:
    vid_t num_vertices;
    size_t nedges;
    uint40_t *neighbors;
    std::vector<adjlist_t> vdata;

  public:
    graph_t() : num_vertices(0), nedges(0), neighbors(NULL) {}

    ~graph_t()
    {
        if (neighbors)
            free(neighbors);
    }

    void resize(vid_t _num_vertices)
    {
        num_vertices = _num_vertices;
        vdata.resize(num_vertices);
    }

    size_t num_edges() const { return nedges; }

    void build(const std::vector<edge_t> &edges);

    void build_reverse(const std::vector<edge_t> &edges);

    adjlist_t &operator[](size_t idx) { return vdata[idx]; };
    const adjlist_t &operator[](size_t idx) const { return vdata[idx]; };
};

#endif