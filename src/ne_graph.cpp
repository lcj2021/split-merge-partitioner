#include "ne_graph.hpp"

void graph_t::build(const std::vector<edge_t> &edges)
{
    if (edges.size() > nedges) {
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    }
    
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++) {
        if (edges[i].valid()) {
            count[edges[i].first]++;
        }
    }

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = adjlist_t(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i].valid()) {
            vdata[edges[i].first].push_back(i);
        }
    }
}

void graph_t::build_reverse(const std::vector<edge_t> &edges)
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++) {
        if (edges[i].valid()) {
            count[edges[i].second]++;
        }
    }

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v - 1];
        vdata[v] = adjlist_t(neighbors + count[v - 1]);
    }
    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i].valid()) {
            vdata[edges[i].second].push_back(i);
        }
    }
}
