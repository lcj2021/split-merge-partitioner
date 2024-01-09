#include "conversions.hpp"
#include "hep_partitioner.hpp"

template <typename TAdj>
HepPartitioner<TAdj>::HepPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename, !need_k_split && FLAGS_write)
{

    Timer convert_timer;
    convert_timer.start();

    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter); //converts the original edgelist in a suitable format for this application
    delete converter;

    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time(); 
    total_time.start();
    LOG(INFO) << "initializing partitioner";
    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate); //are the edges in the format from convert(basefilename, converter)
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);
    //Before the actual edges are coming in the binedgelist-file, there is an entry of the number of vertices and an entry of number of edges
    fin.read((char *)&num_vertices, sizeof(num_vertices)); 
    fin.read((char *)&num_edges, sizeof(num_edges));
    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t), filesize);

    p = FLAGS_p;
    if (need_k_split) {
        p *= FLAGS_k;
    }

    degrees.resize(num_vertices, 0);

    lambda = FLAGS_lambda; //for weighing in balancing score in streaming
    extended_metrics = FLAGS_extended_metrics; // displaying extended metrics
    write_out_partitions = FLAGS_write; // writing out partitions to file
    write_low_degree_edgelist = FLAGS_write_low_degree_edgelist; // writing low degree edgelist to file
    stream_random = FLAGS_random_streaming; // random streaming
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0; //will take track of how many edges are assigned to a bucket so far
    capacity = (double)num_edges * BALANCE_RATIO / p + 1; //will be used to as stopping criterion later
    occupied.assign(p, 0);  //Will count how many edges are in one partition

    is_boundarys.assign(p, dense_bitset(num_vertices)); //shows if a vertex is in S of a bucket
    is_in_a_core = dense_bitset(num_vertices); //Shows if a vertex is in ANY C
    is_high_degree = dense_bitset(num_vertices); // whether a vertex has a high degree and is handled differently
    has_high_degree_neighbor = dense_bitset(num_vertices); // whether the vertex has a high degree neighbor (important in assign_remaining function)

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading and constructing...";

    high_degree_factor = FLAGS_hdf;

    load_in_memory(basefilename, fin);
    capacity_in_memory = ((double)num_edges - num_h2h_edges) * BALANCE_RATIO / p + 1;

    if (std::is_same<TAdj, adj_with_bid_t>::value) {
        edgelist2bucket.resize(num_h2h_edges, kInvalidBid);
    }

    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();
}

// these are the extended stats, including degree distributions etc.
template <typename TAdj>
void HepPartitioner<TAdj>::compute_stats()
{ 
	// average degree of vertices in C and in S\C
	eid_t total_degree_C = 0, total_degree_S = 0;
	vid_t vertex_count_C = 0, vertex_count_S = 0;
	vid_t max_degree = 0;
	for (vid_t i = 0; i < num_vertices; ++i) {
		if (degrees[i] > max_degree) {
			max_degree = degrees[i];
		}
		if (is_in_a_core.get(i)) {
			total_degree_C += degrees[i];
			++vertex_count_C;
		} else { 
            // for those not in core, check whether they are in the boundary of any of the first k-1 partitions (the last partition is not built based on expansion)
			for (bid_t j = 0; j < p - 1; ++j) {
				if (is_boundarys[j].get(i)) {
					total_degree_S += degrees[i];
					++vertex_count_S;
					break;
				}
			}
		}
	}
	double avg_deg_C = static_cast<double>(total_degree_C) / vertex_count_C;
	double avg_deg_S = static_cast<double>(total_degree_S) / vertex_count_S;
	double invalidation_fraction = num_invalidated_edges / (num_edges * 2.0);
	LOG(INFO) << "normalized avg degree C " << avg_deg_C / average_degree << std::endl;
	LOG(INFO) << "normalized avg degree S " << avg_deg_S / average_degree << std::endl;
	LOG(INFO) << "fraction of edges invalidated in clean up phase " << invalidation_fraction << std::endl;

	// computing the distribution of degree to replication factor
	std::vector<eid_t> replication_factor_per_vertex_degree(max_degree + 1, 0);
	std::vector<vid_t> num_vertices_per_vertex_degree(max_degree + 1, 0);
	for (vid_t i = 0; i < num_vertices; ++i) {
		vid_t rep_factor = 0;
		for (dense_bitset &is_boundary : is_boundarys) {
			if (is_boundary.get(i)) {
				++rep_factor;
			}
		}
		vid_t degree = degrees[i];
		replication_factor_per_vertex_degree[degree] += rep_factor;
		++num_vertices_per_vertex_degree[degree];
	}
	std::cout << "Format: degree <whitespace> replication_factor" << std::endl;
	eid_t bucket_min_degree = 1;
	eid_t bucket_max_degree = 10;
	vid_t vertices_in_bucket = 0;
	eid_t replication_in_bucket = 0;
	for (vid_t k = 1; k <= max_degree; ++k) {
		vertices_in_bucket += num_vertices_per_vertex_degree[k];
		replication_in_bucket += replication_factor_per_vertex_degree[k];

		if (k == bucket_max_degree) { // make new bucket
			// finalize it print results
			double result = 0.0;
			if (vertices_in_bucket != 0) {
				result = (double) replication_in_bucket / vertices_in_bucket;
			}
			std::cout << bucket_min_degree << " " << bucket_max_degree << " " << result << std::endl;
			std::cout << bucket_min_degree << " " << bucket_max_degree << " " << (double) vertices_in_bucket / num_vertices << std::endl;
			bucket_min_degree = bucket_max_degree + 1;
			bucket_max_degree = bucket_max_degree * 10;
			vertices_in_bucket = 0;
			replication_in_bucket = 0;
		}
	}
}

template <typename TAdj>
void HepPartitioner<TAdj>::load_in_memory(std::string basefilename, std::ifstream &fin) 
{
	mem_graph.high_degree_factor = high_degree_factor;
	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::out ); // *.h2h_edgelist file
	if (write_low_degree_edgelist) {
		mem_graph.low_degree_file.open(lowedgelist_name(basefilename), std::ios_base::binary | std::ios_base::out ); // *.low_edgelist file;
	}
	mem_graph.resize(num_vertices);
	num_h2h_edges = mem_graph.stream_build(fin, num_edges, is_high_degree, has_high_degree_neighbor, degrees, write_low_degree_edgelist);
	mem_graph.h2h_file.close(); //flushed
	if (write_low_degree_edgelist) {
		mem_graph.low_degree_file.close(); //flushed
	}
}


template <typename TAdj>
void HepPartitioner<TAdj>::in_memory_assign_remaining() 
{

	LOG(INFO) << "Assigned edges before assign_remaining: " << assigned_edges << std::endl;

	for (vid_t vid = 0; vid < num_vertices; ++vid) {
		if (!is_in_a_core.get(vid)) {
			auto &neighbors = mem_graph[vid].adj;
			vid_t i = 0;
			for(; i < mem_graph[vid].size_out(); ++i) {
				bid_t target = best_scored_partition(vid, neighbors[i].vid);
                assign_edge(target, vid, neighbors[i]);
			}

			// in case the vertex has high degree neighbors, the edges from
			// those have not been assigned yet, as the hd vertices were ignored
			// in the expansion. Hence, we have to assign those explicitly.
			if (has_high_degree_neighbor.get(vid)) {
				for(; i < mem_graph[vid].size(); ++i) //for the adj_in neighbors
				{
					if (is_high_degree.get(neighbors[i].vid)) {
						bid_t target = best_scored_partition(neighbors[i].vid, vid);
                        assign_edge(target, vid, neighbors[i]);
					}
				}

			}
		}
	}

	LOG(INFO) << "Assigned edges before streaming: " << assigned_edges << std::endl;
	LOG(INFO) << "Assigning edges between high-degree vertices" << std::endl;

	if (stream_random) {
		// random_streaming();
	}
	else {
		hdrf_streaming();
	}

}



// void HepPartitioner::random_streaming() {

// 	LOG(INFO) << "Streaming randomly." << std::endl;
// 	// assign the edges between two high degree vertices

// 	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::in );
// 	mem_graph.h2h_file.seekg(0, std::ios::beg);

// 	std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
// 	size_t chunk_size;
// 	size_t left_h2h_edges = mem_graph.num_h2h_edges;

// 	if (left_h2h_edges >= 100000) {
// 		chunk_size = 100000; // batch read of so many edges
// 	}
// 	else {
// 		chunk_size = left_h2h_edges;
// 	}
// 	stream_edges.resize(chunk_size);

// 	//	LOG(INFO) << "Chunk size is " << chunk_size << endl;


// 	while (left_h2h_edges > 0) { // edges to be read
// 		mem_graph.h2h_file.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
// 		for (size_t i = 0; i < chunk_size; ++i) {

// 			bucket = std::rand() % p; // random bucket

// 			assign_edge(bucket, stream_edges[i].first, stream_edges[i].second);
// 			is_boundarys[bucket].set_bit_unsync(stream_edges[i].first);
// 			is_boundarys[bucket].set_bit_unsync(stream_edges[i].second);

// 		}

// 		left_h2h_edges -= chunk_size;
// 		if (left_h2h_edges < chunk_size) { // adapt chunk size for last batch read
// 		 chunk_size = left_h2h_edges;
// 		}
// 	}
// }


template <typename TAdj>
void HepPartitioner<TAdj>::hdrf_streaming()
{
	LOG(INFO) << "Streaming using HDRF algorithm." << std::endl;
	// assign the edges between two high degree vertices

	mem_graph.h2h_file.open(h2hedgelist_name(basefilename), std::ios_base::binary | std::ios_base::in);
	mem_graph.h2h_file.seekg(0, std::ios::beg);

	std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
	eid_t chunk_size;
	eid_t left_h2h_edges = mem_graph.num_h2h_edges;
    eid_t id_h2h_edges = 0;

	if (left_h2h_edges >= 100000) {
		chunk_size = 100000; // batch read of so many edges
	} else {
		chunk_size = left_h2h_edges;
	}
	stream_edges.resize(chunk_size);


    /*
    * init min_size
    */
    min_size = *std::min_element(occupied.begin(), occupied.end());

	while (left_h2h_edges > 0) { // edges to be read
		mem_graph.h2h_file.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
		for (eid_t i = 0; i < chunk_size; ++i) {
			bucket = best_scored_partition(stream_edges[i].first, stream_edges[i].second); // according to HDRF scoring
			assign_edge(bucket, stream_edges[i].first, stream_edges[i].second, id_h2h_edges++);

			if (occupied[bucket] > max_size) {
				max_size = occupied[bucket];
			}
			if (occupied[bucket] == min_size) {
				bid_t min_sized_bucket_count = 0;
				for (bid_t b = 0; b < p; ++b) {
					if (occupied[b] == min_size) {
						++min_sized_bucket_count;
					}
				}
				if (min_sized_bucket_count == 1) {
					++min_size;
				}
			}
		}

		left_h2h_edges -= chunk_size;
		if (left_h2h_edges < chunk_size) { // adapt chunk size for last batch read
		    chunk_size = left_h2h_edges;
		}
	}
    mem_graph.h2h_file.close();
}

template <typename TAdj>
void HepPartitioner<TAdj>::in_memory_clean_up_neighbors(vid_t vid, dense_bitset & is_core, dense_bitset & is_boundary) 
{
	mem_adjlist_t<TAdj> &neighbors = mem_graph[vid];

	const vid_t num_neigh_out = neighbors.size_out();
	const vid_t num_neigh_size = neighbors.size();
	vid_t i = 0;
	vid_t j = 0;

	for(; j < num_neigh_out; ++j) {
	    //naive vid_t u = neighbors.get_neighbor(neighbors_file, i);
		vid_t u = neighbors.adj[i].vid;
	    if (is_core.get(u)) { // neighbor u is in core, so edge is removed
	    	++num_invalidated_edges;
	    	neighbors.erase_out(i);
	    } else if (is_boundary.get(u)) { // neighbor u is in boundary, so edge is removed
	        ++num_invalidated_edges;
	    	neighbors.erase_out(i);
	    } else if (is_high_degree.get(u)) {
	    	++num_invalidated_edges;
	    	neighbors.erase_out(i);
	    } else {
	        ++i;
	    }
	 }

	 for(; j < num_neigh_size; ++j) {
		 vid_t u = neighbors.adj[i].vid;
	     if (is_core.get(u)) { // neighbor u is in core, so edge is removed
	    	 ++num_invalidated_edges;
	         neighbors.erase_in(i);
	     } else if (is_boundary.get(u)) { // neighbor u is in boundary, so edge is removed
	    	 ++num_invalidated_edges;
	         neighbors.erase_in(i);
	     } else if (is_high_degree.get(u)) {
			 ++num_invalidated_edges;
			 neighbors.erase_in(i);
		 } else {
        	 ++i;
	     }
	}
}


template <typename TAdj>
void HepPartitioner<TAdj>::partition_in_memory() 
{
	bool expansion_finished = false;

    for (bucket = 0; bucket < p - 1; ++bucket) {
        LOG(INFO) << bucket << ", ";

        //DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (occupied[bucket] < capacity_in_memory) {
            vid_t d, vid = 0;
            if (!min_heap.get_min(d, vid)) {
                if (!in_memory_get_free_vertex(vid)) {
                    LOG(INFO) << "partition " << bucket
                               << " stop: no free vertices";
                    expansion_finished = true;
                    break;
                }
                d= mem_graph[vid].size(); // a high degree vertex will not be chosen by get free vertex. also will not be in min heap.
            } else {
                min_heap.remove(vid);
            }

            in_memory_occupy_vertex(vid, d);

        }

        /*
         * clean up the adjacency lists of vertices from the S set from the last bucket
         */
        auto &is_core = is_in_a_core, &is_boundary = is_boundarys[bucket];

        std::vector<std::pair<vid_t, vid_t>> heap = min_heap.getHeap();
        vid_t size = min_heap.getSize(); // do not iterate over heap more than "size" entries
        /*
         * vid is the vertex in S
         * u is the neighbor of vid in the currently examined edge
         */
        vid_t count = 0;
        for (std::vector<std::pair<vid_t, vid_t>>::iterator it = heap.begin(); it != heap.end(); ++it) {
        	if (count >= size) {
        		break; //stop here
        	}

        	++count;
        	vid_t vid = it->second;

        	in_memory_clean_up_neighbors(vid, is_core, is_boundary);
        }

        min_heap.clear();

        for (vid_t vid : vid_id_not_in_boundary) {
            is_boundary.set_unsync(vid, 0);
        }
        vid_id_not_in_boundary.clear();

        if (expansion_finished) {
        	break;
        }
    }

    in_memory_assign_remaining();
    LOG(INFO) << "Finished partitioning" << std::endl;
    LOG(INFO) << "Core vertice count: " << is_in_a_core.popcount() << ' ' << (double)is_in_a_core.popcount() / num_vertices;
    writer.fout.close();; // flushing out the unwritten edge assignments

}

template <typename TAdj>
double HepPartitioner<TAdj>::compute_partition_score(vid_t u, vid_t v, bid_t bucket_id) 
{
	if (occupied[bucket_id] >= capacity) {
		return -1.0; // partition is full, do not choose it
	}
	vid_t degree_u = degrees[u];
	vid_t degree_v = degrees[v];
	vid_t sum = degree_u + degree_v;
	double gu = 0.0, gv = 0.0;
	if (is_boundarys[bucket_id].get(u)) {
		gu = degree_u;
		gu /= sum;
		gu = 1 + (1 - gu);
	}
	if (is_boundarys[bucket_id].get(v)) {
		 gv = degree_v;
		 gv /= sum;
		 gv = 1 + (1 - gv);
	}

	double bal = (max_size - occupied[bucket_id]) / (1.0 + max_size - min_size);

	double score = gu + gv + lambda * bal;
	return score;
}

template <typename TAdj>
bid_t HepPartitioner<TAdj>::best_scored_partition(vid_t u, vid_t v) 
{
	double best_score = -1.0;
	bid_t best_partition = kInvalidBid;
	for (bid_t i = 0; i < p; ++i) {
		double score = compute_partition_score(u, v, i);
		if (score > best_score) {
			best_score = score;
			best_partition = i;
		}
	}
    if (best_partition == kInvalidBid) {
        best_partition = gen() % p;
    }
	return best_partition;
}

template <typename TAdj>
void HepPartitioner<TAdj>::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();


    partition_in_memory();

    compute_timer.stop();

    calculate_stats();

    CHECK_EQ(assigned_edges, num_edges);
    CHECK_EQ(check_edge_hybrid(), 1);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    /*
     * compute some stats about the partitioned graph (for further analysis)
     * if extended_metrics flag is set
     */
    if (extended_metrics) {
    	compute_stats();
    }
}

template <typename TAdj>
void HepPartitioner<TAdj>::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<vid_t> num_bucket_vertices(p, 0);
    for (bid_t b = 0; b < p; ++b) {
        num_bucket_vertices[b] = is_boundarys[b].popcount();
    }
    vid_t max_part_vertice_cnt = *std::max_element(num_bucket_vertices.begin(), num_bucket_vertices.end());
    vid_t all_part_vertice_cnt = accumulate(num_bucket_vertices.begin(), num_bucket_vertices.end(), (vid_t)0);
    eid_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
    eid_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (eid_t)0);

    for (bid_t b = 0; b < p; ++b) {
        LOG(INFO) << "Bucket_info: " << b
                << ", vertices: " << num_bucket_vertices[b]
                << ", edges: " << occupied[b];
    }
    
    double avg_vertice_cnt = static_cast<double>(all_part_vertice_cnt) / p;
    double avg_edge_cnt = static_cast<double>(all_part_edge_cnt) / p;

    double std_vertice_deviation = 0.0;
    double std_edge_deviation = 0.0;
    for (bid_t b = 0; b < p; ++b) {
        std_vertice_deviation += pow(num_bucket_vertices[b] - avg_vertice_cnt, 2);
        std_edge_deviation += pow(occupied[b] - avg_edge_cnt, 2);
    }
    std_vertice_deviation = sqrt(static_cast<double>(std_vertice_deviation) / p);
    std_edge_deviation = sqrt(static_cast<double>(std_edge_deviation) / p);

    
    LOG(INFO) << std::string(20, '#') << "\tVertice    balance\t" << std::string(20, '#');
    LOG(INFO) << "Max vertice count / avg vertice count: "
              << static_cast<double>(max_part_vertice_cnt) / (num_vertices / p);
    LOG(INFO) << "Max Vertice count: "
              << max_part_vertice_cnt;
    LOG(INFO) << "Avg Vertice count(No replicate): "
              << num_vertices / p;
    LOG(INFO) << "Vertice std_vertice_deviation / avg: "
              << std_vertice_deviation / avg_vertice_cnt;

    LOG(INFO) << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#');
    LOG(INFO) << "Max edge count / avg edge count: "
              << static_cast<double>(max_part_edge_cnt) / avg_edge_cnt;
    LOG(INFO) << "Max Edge count: "
              << max_part_edge_cnt;
    LOG(INFO) << "Avg Edge count: "
              << num_edges / p;
    LOG(INFO) << "Edge std_edge_deviation / avg: "
              << std_edge_deviation / avg_edge_cnt;

    CHECK_EQ(all_part_edge_cnt, num_edges);
    LOG(INFO) << std::string(20, '#') << "\tReplicate    factor\t" << std::string(20, '#');
    LOG(INFO) << "replication factor (final): " << static_cast<double>(all_part_vertice_cnt) / num_vertices;
}