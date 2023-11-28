#ifndef HEP_PARTITIONER_HPP
#define HEP_PARTITIONER_HPP

#include <random>

#include "hep_min_heap.hpp"
#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"
#include "hep_graph.hpp"

/* Hybrid Edge Partitioner (HEP) */
template <typename TAdj>
class HepPartitioner : public Partitioner
{
private:
    std::random_device rd;
    std::mt19937 gen;
    const double BALANCE_RATIO = 1.0;
    double lambda;
    bool extended_metrics;

    bool stream_random; // whether to use random streaming for the second phase of HEP

    std::string basefilename;

    // vid_t num_vertices;
    eid_t assigned_edges, num_h2h_edges;
    bid_t p, bucket; 
    double average_degree;
    eid_t capacity;
    eid_t capacity_in_memory; // capacity per partition of in memory partitioning
    eid_t num_invalidated_edges; // number of edges removed in clean up phase overall
    eid_t min_size = 0; // currently smallest partition
    eid_t max_size = 0; // currently largest partition

    bool write_out_partitions = false; // whether the partitions should be written to the out-file or not
    bool write_low_degree_edgelist = false; // whether edges incident to a low-degree vertex should be written out to a file. useful if this sub-graph should be analyzed separately.

    // std::vector<edge_t> edges;
    // mem_graph_t<TAdj> mem_graph; // graph for in-memory processing
    double high_degree_factor;
    HepMinHeap<vid_t, vid_t> min_heap;
    // std::vector<eid_t> occupied;
    // std::vector<dense_bitset> is_boundarys; 
    dense_bitset is_in_a_core;
    dense_bitset is_high_degree;
    dense_bitset has_high_degree_neighbor;


    vid_t search_index_free_vertex = 0;


    std::vector<vid_t> vid_id_not_in_boundary; // degrees of vertices//(num_vertices, 0);
    edgepart_writer<vid_t, bid_t> writer;
    

    void in_memory_clean_up_neighbors(vid_t vid, dense_bitset & is_core, dense_bitset & is_boundary);

    void assign_edge(bid_t cbucket, vid_t from, TAdj& v_e)
    {
        writer.save_edge(from, v_e.vid, cbucket);
        ++assigned_edges;
        ++occupied[cbucket];
        v_e.bid = cbucket;

        is_boundarys[cbucket].set_bit_unsync(from);
        is_boundarys[cbucket].set_bit_unsync(v_e.vid);
    }

    void assign_edge(bid_t cbucket, vid_t from, vid_t to, eid_t edge_id)
    {        
        writer.save_edge(from, to, cbucket);
        ++assigned_edges;
        ++occupied[cbucket];
        edgelist2bucket[edge_id] = cbucket;

        is_boundarys[cbucket].set_bit_unsync(from);
        is_boundarys[cbucket].set_bit_unsync(to);
    }

    void in_memory_add_boundary(vid_t vid)
    {
    	bool bucket_full = false, bucket_full_at_start = false;
    	auto &is_core = is_in_a_core, &is_boundary = is_boundarys[bucket];

        // vid is already in boundary!
		if (is_boundary.get(vid)) {
			return; 
		}

		is_boundary.set_bit_unsync(vid);

        // high degree vertices are treated as if they were in the core
		if (is_high_degree.get(vid)) { 
			is_in_a_core.set_bit_unsync(vid);
            // high degree vertices are ignored in the normal expansion. we do not at all look at their neighbors. we do also not add them to min heap.
			return; 
		}

		bool vid_is_in_core = is_core.get(vid);

		if (!vid_is_in_core) {
            // is not in the core yet: Potential next candidate --> insert in MinHeap
			min_heap.insert(mem_graph[vid].size(), vid); 
		}
		auto &neighbors = mem_graph[vid].adj;
		vid_t count = 0;
		for (; count < mem_graph[vid].size_out(); ++count) //for the adj_out neighbors
		{
			if (occupied[bucket] >= capacity) {
                // full, stop adding vertices to the boundary of this bucket
				bucket_full = true;
			} 
            if (count == 0 && bucket_full) {
                bucket_full_at_start = true;
            }
            
			TAdj &u = neighbors[count];

            // high degree vertices are always considered to be in c
			if (is_high_degree.get(u.vid)) { 
				if (!bucket_full) {
                    // assign edge --> vid is the left vertex
                    assign_edge(bucket, vid, u);
					if (!vid_is_in_core) {
                        // vid has one neighbor less now
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
					}
				} else { 
                    // bucket is full; assign to next bucket
                    assign_edge(bucket + 1, vid, u);
				}
			} else {
                // If the neighbor of vid is in core
				if (is_core.get(u.vid)) { 
					if (!bucket_full) {
                        // assign edge --> vid is the left vertex
                        assign_edge(bucket, vid, u);
						if (!vid_is_in_core) {
                            // vid has one neighbor less now
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
						}
					} else {
						// bucket is full; assign to next bucket
                        assign_edge(bucket + 1, vid, u);
					}
				} else if (is_boundary.get(u.vid)) {
					if (!bucket_full) {
                        assign_edge(bucket, vid, u);
						min_heap.decrease_key(u.vid, 1, mem_graph[u.vid].size());
						if (!vid_is_in_core) {
                            // vid has one neighbor less now
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
						}
					} else {
						// bucket is full; assign to next bucket
                        assign_edge(bucket + 1, vid, u);
					}
				}
			}
		}
        // for the adj_in neighbors
		for (; count < mem_graph[vid].size(); ++count) {
			if (occupied[bucket] >= capacity) {
				bucket_full = true;
			} 
            // full, stop adding vertices to the boundary
            if (count == 0 && bucket_full) {
                bucket_full_at_start = true;
            }

			TAdj &u = neighbors[count];

            // high degree vertices are always considered to be in c
			if (is_high_degree.get(u.vid)) { 
				if (!bucket_full) {
                    // assign edge --> vid is the right vertex
                    assign_edge(bucket, vid, u);
					if (!vid_is_in_core) {
                        // vid has one neighbor less now
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
					}
				} else {
					// bucket is full; assign to next bucket
                    assign_edge(bucket + 1, vid, u);
				}
			} else {
				if (is_core.get(u.vid)) {
					if (!bucket_full) {
                        // vid is on the right side
                        assign_edge(bucket, vid, u);
						if (!vid_is_in_core) {
                            // vid has one neighbor less now
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
						}
					} else {
						// bucket is full; assign to next bucket
                        assign_edge(bucket + 1, vid, u);
					}
				} else if (is_boundary.get(u.vid)) {
					if (!bucket_full) {
                        assign_edge(bucket, vid, u);
						min_heap.decrease_key(u.vid, 1, mem_graph[u.vid].size());
						if (!vid_is_in_core) {
                            // vid has one neighbor less now
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); 
						}
					} else {
						// bucket is full; assign to next bucket
                        assign_edge(bucket + 1, vid, u);
					}
				}
			}
		}
        if (bucket_full_at_start) {
            vid_id_not_in_boundary.emplace_back(vid);
        }
    }


    void in_memory_occupy_vertex(vid_t vid, vid_t d)
    {
    	CHECK(!is_in_a_core.get(vid)) << "add " << vid << " to core again";

    	is_in_a_core.set_bit_unsync(vid);
        // also applies to high degree vertices, as their in-memory degree is 0
    	if (d == 0) { 
    	    return;
    	}
    	in_memory_add_boundary(vid);
        // Set all neighbors of vid to boundary
    	for (vid_t i = 0; i < mem_graph[vid].size(); ++i) { 
    		in_memory_add_boundary(mem_graph[vid].adj[i].vid);
    	}
    }

    bool in_memory_get_free_vertex(vid_t &vid)
    {
        vid = search_index_free_vertex;

      	/*
       	* find a vertex to start expansion with
       	*/
       	while ((mem_graph[vid].size() == 0 || is_in_a_core.get(vid)) && vid < num_vertices) {
       	   	++vid;
       	}

       	search_index_free_vertex = vid;
        // searched to the end, did not find free vertex
        return vid != num_vertices;
    }

    bool check_edge_hybrid()
    {
        std::vector<dense_bitset> dbitsets(p, dense_bitset(num_vertices));
        for (vid_t vid = 0; vid < num_vertices; ++vid) {
            bool assigned_to_a_part = false;
            for (bid_t b = 0; b < p; ++b) {
                if (is_boundarys[b].get(vid)) {
                    assigned_to_a_part = true;
                    break;
                }
            }
            if (!assigned_to_a_part) {
                return false;
            }
        }

        return true;
    }

    void load_in_memory(std::string basefilename, std::ifstream &fin);
    void partition_in_memory();
    void in_memory_assign_remaining();

    double compute_partition_score(vid_t u, vid_t v, bid_t bucket_id); // returns HDRF score for edge (u,v) on partition <bucket_id>
    bid_t best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)

    void compute_stats();

    // void random_streaming();
    void hdrf_streaming();


public:
    HepPartitioner(std::string basefilename, bool need_k_split);
    void split();

    void calculate_stats();
};

template class HepPartitioner<vid_eid_t>;
// template class HepPartitioner<vid_t>;

#endif