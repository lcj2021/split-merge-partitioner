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

#include <set>

#include "util.hpp"
#include "hep_min_heap.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"
#include "hep_graph.hpp"

/* Hybrid Edge Partitioner (HEP) */
class HepPartitioner : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.0;
    double lambda;
    bool extended_metrics;

    bool two_ps; // whether to use 2ps restreaming for the second phase of HEP
    bool stream_random; // whether to use random streaming for the second phase of HEP

    std::string basefilename;

    // vid_t num_vertices;
    size_t assigned_edges, num_h2h_edges;
    int p, bucket; 
    double average_degree;
    size_t capacity;
    size_t capacity_in_memory; // capacity per partition of in memory partitioning
    size_t invalidated_edges_count; // number of edges removed in clean up phase overall
    size_t min_size = 0; // currently smallest partition
    size_t max_size = 0; // currently largest partition

    bool write_out_partitions = false; // whether the partitions should be written to the out-file or not
    bool write_low_degree_edgelist = false; // whether edges incident to a low-degree vertex should be written out to a file. useful if this sub-graph should be analyzed separately.

    // std::vector<edge_t> edges;
    mem_graph_t mem_graph; // graph for in-memory processing
    double high_degree_factor;
    HepMinHeap<vid_t, vid_t> min_heap;
    // std::vector<size_t> occupied;
    // std::vector<dense_bitset> is_boundarys; 
    dense_bitset is_in_a_core;
    dense_bitset is_high_degree;
    dense_bitset has_high_degree_neighbor;
    std::vector<size_t> count; // degrees of vertices//(num_vertices, 0);


    vid_t search_index_free_vertex = 0;


    std::vector<vid_t> vid_id_not_in_boundary; // degrees of vertices//(num_vertices, 0);
    edgepart_writer<vid_t, uint16_t> writer;
    

    void in_memory_clean_up_neighbors(vid_t vid, dense_bitset & is_core, dense_bitset & is_boundary);

    std::vector<bool> assigned;
    void assign_edge(int cbucket, vid_t from, vid_t to, vid_t edge_id)
    {        
        assigned[edge_id] = true;
        writer.save_edge(from, to, cbucket);
        assigned_edges++;
        occupied[cbucket]++;
        CHECK_EQ(edge2bucket[edge_id], -1);
        edge2bucket[edge_id] = cbucket;

        // vid_t u = edges[edge_id].first, v = edges[edge_id].second;
        // // if (!(u == from and v == to) and !(u == to and v == from)) {
        // if (!(u == from and v == to)) {
        //     LOG(FATAL) << "BUG";
        // }

        // if (from == 17340 or to == 17340) {
        //     LOG(INFO) << "17340 -> " << cbucket << ", edge_id: " << edge_id;
        // }

        is_boundarys[cbucket].set_bit_unsync(from);
        is_boundarys[cbucket].set_bit_unsync(to);

    }

    void in_memory_add_boundary(vid_t vid){

    	bool bucket_full = false, bucket_full_at_start = false;
    	auto &is_core = is_in_a_core, &is_boundary = is_boundarys[bucket];

		if (is_boundary.get(vid)){
			return; //vid is already in boundary!
		}

		is_boundary.set_bit_unsync(vid);
        if (vid == 17340) {
            LOG(INFO) << "17340, " << "in_memory_add_boundary " << bucket;
            LOG(INFO) << "mem_graph[vid].size_out(): " << mem_graph[vid].size_out();
            LOG(INFO) << "mem_graph[vid].size(): " << mem_graph[vid].size();
        }

		if (is_high_degree.get(vid)){ // high degree vertices are treated as if they were in the core
			is_in_a_core.set_bit_unsync(vid);
			return; // high degree vertices are ignored in the normal expansion. we do not at all look at their neighbors. we do also not add them to min heap.
		}

		bool vid_is_in_core = is_core.get(vid);

		if (!vid_is_in_core) {
			min_heap.insert(mem_graph[vid].size(), vid); //is not in the core yet: Potential next candidate --> insert in MinHeap
		}
		auto &neighbors = mem_graph[vid].adj;
		vid_t count = 0;
		for(; count < mem_graph[vid].size_out(); count++) //for the adj_out neighbors
		{
			if (occupied[bucket] >= capacity){
				bucket_full = true;
			} // full, stop adding vertices to the boundary of this bucket
            if (count == 0 && bucket_full) {
                bucket_full_at_start = true;
            }
            

			vid_eid_t &u = neighbors[count];

			if (is_high_degree.get(u.vid)){ // high degree vertices are always considered to be in c
				if (!bucket_full){
					assign_edge(bucket, vid , u.vid, u.eid ); //assign edge --> vid is the left vertex
					if (!vid_is_in_core) {
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
					}
				}
				else{ //bucket is full; assign to next bucket
					assign_edge(bucket + 1, vid , u.vid, u.eid );
				}

			}
			else{
				if(is_core.get(u.vid)){ //If the neighbor of vid is in core
					if (!bucket_full){
						assign_edge(bucket, vid , u.vid, u.eid ); //assign edge --> vid is the left vertex
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, vid , u.vid, u.eid );
					}
				}else if (is_boundary.get(u.vid)) {
					if (!bucket_full){
						assign_edge(bucket, vid , u.vid, u.eid );
						min_heap.decrease_key(u.vid, 1, mem_graph[u.vid].size());
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, vid , u.vid, u.eid );
					}
				}
			}

		}
		for(; count < mem_graph[vid].size(); count++) //for the adj_in neighbors
		{
			if (occupied[bucket] >= capacity){
				bucket_full = true;
			} // full, stop adding vertices to the boundary
            if (count == 0 && bucket_full) {
                bucket_full_at_start = true;
            }
            if (count == 0 && vid == 17340) {
                LOG(INFO) << "bucket_full ? " << bucket_full;
                LOG(INFO) << "bucket_full ? " << bucket_full;
            }

			vid_eid_t &u = neighbors[count];

			if (is_high_degree.get(u.vid)){ // high degree vertices are always considered to be in c
				if (!bucket_full){
					assign_edge(bucket, u.vid , vid, u.eid ); //assign edge --> vid is the right vertex
					if (!vid_is_in_core) {
						min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
					}
				}
				else {
					//bucket is full; assign to next bucket
					assign_edge(bucket + 1, u.vid , vid, u.eid );
				}
			}
			else{
				if(is_core.get(u.vid)){
					if (!bucket_full){
						assign_edge(bucket, u.vid, vid, u.eid ); //vid is on the right side
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, u.vid , vid, u.eid );
					}
				}else if (is_boundary.get(u.vid)) {
					if (!bucket_full){
						assign_edge(bucket, u.vid, vid, u.eid );
						min_heap.decrease_key(u.vid, 1, mem_graph[u.vid].size());
						if (!vid_is_in_core) {
							min_heap.decrease_key(vid, 1, mem_graph[vid].size()); //vid has one neighbor less now
						}
					}
					else {
						//bucket is full; assign to next bucket
						assign_edge(bucket + 1, u.vid , vid, u.eid );
					}
				}
			}
		}
        if (bucket_full_at_start) {
		    // is_boundary.set_bit_unsync(vid);
            LOG(INFO) << vid << ", bucket_full_at_start";
            vid_id_not_in_boundary.emplace_back(vid);
        }
    }


    void in_memory_occupy_vertex(vid_t vid, vid_t d){

    	CHECK(!is_in_a_core.get(vid)) << "add " << vid << " to core again";

    	is_in_a_core.set_bit_unsync(vid);
    	if (d == 0){ //also applies to high degree vertices, as their in-memory degree is 0
    	    return;
    	}
    	in_memory_add_boundary(vid);
    	for(vid_t i = 0; i < mem_graph[vid].size(); i++) { //Set all neighbors of vid to boundary
    		in_memory_add_boundary(mem_graph[vid].adj[i].vid);
    	}
    }

    bool in_memory_get_free_vertex(vid_t &vid){

       vid = search_index_free_vertex;

      	/*
       	* find a vertex to start expansion with
       	*/
       	while ((mem_graph[vid].size() == 0 || is_in_a_core.get(vid)) && vid < num_vertices) {
       	   	vid++;
       	}

       	search_index_free_vertex = vid;
       	if (vid == num_vertices){
       		return false;
       	} // searched to the end, did not find free vertex
       	else{
           	return true;
       	}
    }

    // bool check_edge()
    // {
    //     LOG(INFO) << std::accumulate(assigned.begin(), assigned.end(), 0);
    //     std::vector<dense_bitset> dbitsets(p, dense_bitset(num_vertices));
    //     for (vid_t vid = 0; vid < num_vertices; ++ vid) {
    //         bool assigned_to_a_part = false;
    //         for (int b = 0; b < p; ++ b) {
    //             if (is_boundarys[b].get(vid)) {
    //                 assigned_to_a_part = true;
    //                 break;
    //             }
    //         }
    //         if (!assigned_to_a_part) {
    //             LOG(FATAL) << "BUG, vid = " << vid;
    //             return false;
    //         }
    //     }
    //     LOG(INFO) << num_edges;
    //     for (size_t edge_id = 0; edge_id < num_edges; ++ edge_id) {
    //         // edges[edge_id].recover();
    //         int16_t edge_bucket = edge2bucket[edge_id];
    //         if (edge_id == 30202525) {
    //             LOG(INFO) << "eid 30202525 -> " << edge_bucket;
    //         }
    //         CHECK_NE(edge_bucket, -1);
    //         vid_t u = edges[edge_id].first, v = edges[edge_id].second;
    //         dbitsets[edge_bucket].set_bit_unsync(u);
    //         dbitsets[edge_bucket].set_bit_unsync(v);
    //     }
    //     size_t r0 = 0, r1 = 0;
    //     for (int b = 0; b < p; ++ b) {
    //         r0 = dbitsets[b].popcount();
    //         r1 = is_boundarys[b].popcount();
    //         LOG(INFO) << "bucket: " << b << ' ' << r0 << ' ' << r1;
    //         // if (b == 2 or b == 3) {
    //         if (r0 != r1) {
    //             for (vid_t vid = 0; vid < num_vertices; ++vid) {
    //                 if (!dbitsets[b].get(vid) and is_boundarys[b].get(vid)) {
    //                     LOG(INFO) << "vid: " << vid << " not in dbitsets " << b;
    //                 }
    //             }
    //         }
    //     }
    //     auto equal_dbitset = [&](const dense_bitset &l, const dense_bitset &r) {
    //         if (l.size() != r.size()) return false;
    //         for (size_t bit = 0; bit < l.size(); ++ bit) {
    //             if (l.get(bit) != r.get(bit)) {
    //                 if (l.get(bit) and !r.get(bit)) {
    //                     LOG(FATAL) << "LESS";
    //                 } else {
    //                     LOG(FATAL) << "MORE";
    //                 }
    //                 return false; 
    //             }
    //         }   
    //         return true;
    //     };
    //     for (int b = 0; b < p; ++ b) {
    //         if (!equal_dbitset(dbitsets[b], is_boundarys[b])) {
    //             LOG(FATAL) << "BUG, b = " << b;
    //             return false;
    //         }
    //     }
    //     return true;
    // }

    void load_in_memory(std::string basefilename, std::ifstream &fin);
    void partition_in_memory();
    void in_memory_assign_remaining();

    double compute_partition_score(vid_t u, vid_t v, int bucket_id); // returns HDRF score for edge (u,v) on partition <bucket_id>
    int best_scored_partition(vid_t u, vid_t v); // returns bucket id where score is best for edge (u,v)

    size_t count_mirrors();
    void compute_stats();

    // void random_streaming();
    void hdrf_streaming();


  public:
    HepPartitioner(std::string basefilename, bool need_k_split);
    void split();

    void calculate_stats();

    //needed for testing
    std::vector<dense_bitset> get_boundary(){
        return is_boundarys;
    }

	vid_t getNumVertices() const {
		return num_vertices;
	}

	void setNumVertices(vid_t numVertices) {
		num_vertices = numVertices;
	}
};
