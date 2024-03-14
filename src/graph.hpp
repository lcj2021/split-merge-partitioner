#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <utility>

#include "common.hpp"
#include "dense_bitset.hpp"

struct AdjEntryVidBid {
    vid_t vid;
    bid_t bid;
    AdjEntryVidBid(vid_t vid, bid_t bid) : vid(vid), bid(bid) {}
    AdjEntryVidBid(vid_t vid) : vid(vid), bid(kInvalidBid) {}
} __attribute__((packed));

struct AdjEntryVid {
    vid_t vid;
    AdjEntryVid(vid_t vid) : vid(vid) {}
};

template <typename TAdj>
class AdjList 
{
public:
	TAdj *adj; // link into the column array
    vid_t len_out;
    vid_t len_in;

public:

    AdjList() : adj(NULL), len_out(0), len_in(0) {}
    AdjList( TAdj *adj) :  adj(adj), len_out(0), len_in(0) {}

    TAdj *begin() { return adj; }
    TAdj *end() { return adj + len_out + len_in; }

    TAdj &operator[](vid_t idx) { return adj[idx]; };

    vid_t size() const { return len_out + len_in; }
    vid_t size_out() const { return len_out; }
    vid_t size_in() const { return len_in; }

    void push_back_out(TAdj data) 
    {
    	adj[len_out++] = data;
    }

    void push_back_in(TAdj data, vid_t offset) 
    {
        // incoming neighbors are written after the outgoing neighbors.
        adj[len_in++ + offset] = data; 
    }


    void erase_out(vid_t to_erase_pos) 
    {
        std::swap(adj[to_erase_pos], adj[size_out() - 1]);
        std::swap(adj[size_out() - 1], adj[size() - 1]);
        
    	pop_back_out();
    }

    void erase_in(vid_t to_erase_pos) 
    {
        std::swap(adj[to_erase_pos], adj[size() - 1]);
        pop_back_in();
    }

    void pop_back_out(){len_out--;}
    void pop_back_in(){len_in--;}
};

template class AdjList<AdjEntryVidBid>;

/// @brief: in-memory graph
template <typename TAdj>
class Graph 
{

public:
    vid_t num_vertices;
    eid_t num_edges;
    TAdj *neighbors;
    std::vector<AdjList<TAdj>> vdata;

public:
    Graph() : num_vertices(0), num_edges(0), neighbors(NULL) {  }

    Graph(Graph&& other) noexcept
        : num_vertices(std::exchange(other.num_vertices, 0)),
          num_edges(std::exchange(other.num_edges, 0)),
          neighbors(std::exchange(other.neighbors, nullptr)),
          vdata(std::move(other.vdata))
    {
    }

    Graph& operator=(Graph&& other) noexcept 
    {
        num_vertices = std::exchange(other.num_vertices, 0);
        num_edges = std::exchange(other.num_edges, 0);
        neighbors = std::exchange(other.neighbors, nullptr);
        vdata = std::move(other.vdata);
        return *this;
    }

    ~Graph()
    {
        if (neighbors) 
            free(neighbors);
    }

    void resize(vid_t _num_vertices)
    {
        num_vertices = _num_vertices;
        vdata.resize(num_vertices);
    }

    eid_t stream_build(std::ifstream &fin, eid_t num_edges, std::vector<vid_t> &count);

    AdjList<TAdj> &operator[](eid_t idx) 
    { 
        return vdata[idx]; 
    }

    const AdjList<TAdj> &operator[](eid_t idx) const 
    {
        return vdata[idx];
    }


};

template class Graph<AdjEntryVidBid>;
template class Graph<AdjEntryVid>;

#endif