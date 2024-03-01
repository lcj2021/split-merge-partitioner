#include <numeric>

#include "graph.hpp"

#include "conversions.hpp"

// returns number of h2h edges
template <>
eid_t Graph<AdjEntryVidBid>::stream_build(
    std::ifstream &fin, eid_t num_edges, std::vector<vid_t> &degrees
)
{
	// fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);

	// LOG(INFO) << "builder starts...";
	// std::vector<vid_t> offsets(num_vertices, 0); // to put the in-neighbors at the right position when building the column array

	// nedges = num_edges; // num_edges, num_vertices
	// double average_degree = num_edges * 2.0 / num_vertices; // non-rounded average degree
	// high_degree_threshold = average_degree * high_degree_factor; // this is the th, if exceeded, the node is ignored for csr

	// LOG(INFO) << "Average degree: " << average_degree << std::endl;
	// LOG(INFO) << "High degree threshold: " << high_degree_threshold << std::endl;

    // std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
	// eid_t chunk_size;

	// if (num_edges >= 100000) {
	// 	chunk_size = 100000; // batch read of so many edges
	// } else {
	// 	chunk_size = num_edges;
	// }
	// stream_edges.resize(chunk_size);
    // fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);

    // while (num_edges > 0) { // edges to be read
	// 	fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
	//     for (eid_t i = 0; i < chunk_size; ++i) {
	//     	++degrees[stream_edges[i].first];
	//     	++degrees[stream_edges[i].second];
	//     	++offsets[stream_edges[i].first];
	//     }
	//     num_edges -= chunk_size;
	//     if (num_edges < chunk_size) { // adapt chunk size for last batch read
	//     	chunk_size = num_edges;
	//     }
	// }
	// // now, degrees are complete. degrees are the degrees of the vertices.

	// /************************
	//  * build the index array
	//  * **********************
	//  */
	// vid_t h_count = 0; // how many high degree vertices are found
	// if (degrees[0] > high_degree_threshold) {
	// 	is_high_degree.set_bit_unsync(0);
	// 	++h_count;
	// }

	// std::vector<eid_t> index(num_vertices, 0); // for index array

	// for (vid_t v = 1; v < num_vertices; ++v) {
	// 	// we ignore the degrees of vertices that have a high degree; we also ignore them when building the column array
	// 	if (degrees[v - 1] <= high_degree_threshold) {
	// 		index[v] = index[v - 1] + degrees[v - 1];
	// 	} else {
	// 		index[v] = index[v - 1]; // ignoring v - 1, will not use it in CSR
	// 	}

	// 	if (degrees[v] > high_degree_threshold) {
	// 		is_high_degree.set_bit_unsync(v);
	// 		++h_count;
	// 	}
	// }
    // eid_t neighbors_len = index[num_vertices - 1] + degrees[num_vertices - 1];
    // LOG(INFO) << "neighbors_len: " << neighbors_len;

    // neighbors = (AdjEntryVidBid *)realloc(neighbors, sizeof(AdjEntryVidBid) * neighbors_len); // store 2 vids for each edge
    // LOG(INFO) << sizeof(AdjEntryVidBid) << " bytes needed for AdjEntryVidBid";
    // LOG(INFO) << (sizeof(AdjEntryVidBid) * neighbors_len / 1024.0 / 1024 / 1024) << " G bytes needed for neighbors";
	// CHECK(neighbors) << "allocation failed!";

    // #pragma omp parallel for
    // for (eid_t i = 0; i < neighbors_len; ++i) {
    //     neighbors[i].bid = kInvalidBid;
    // }

    // #pragma omp parallel for
    // for (vid_t vid = 0; vid < num_vertices; ++vid) {
	// 	vdata[vid] = mem_adjlist_t(neighbors + index[vid]);
    // }

	// LOG(INFO) << "Number of vertices with high degree " << h_count << std::endl;
	// std::streampos pos(0);
	// h2h_file.seekp(pos);

	// std::streampos pos_low(0);
	// low_degree_file.seekp(pos_low);
	// low_degree_file.write((char *)&num_vertices, sizeof(num_vertices));
	// low_degree_file.write((char *)&num_edges, sizeof(num_edges));

	// /****************************
	//  * build the column array
	//  * **************************
	//  */
    // num_edges = nedges;
	// //resizing the chunk size
	// if (num_edges >= 100000) {
	//    	chunk_size = 100000; // batch read of so many edges
	// } else {
	//    	chunk_size = num_edges;
	// }

	// fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); // start read from beginning

	// eid_t savings = 0;

    // eid_t edge_id = 0;
    // while (num_edges > 0) { // edges to be read
	//    	fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
	//     for (eid_t i = 0; i < chunk_size; ++i, ++edge_id) {
	//    		const auto& [u, v] = stream_edges[i];
	//    		// we do not build column array for high degree vertices
	//    		bool low_degree = false; // needed in case we write a low_degree edge list out to file
	//    		if (degrees[u] <= high_degree_threshold) {
	//    			vdata[u].push_back_out(AdjEntryVidBid(v));
	//    			low_degree = true;
	//    		} else {
	//    			has_high_degree_neighbor.set_bit_unsync(v);
	//    			++savings;
	//    		}
	//    		if (degrees[v] <= high_degree_threshold) {
	//    			vdata[v].push_back_in(AdjEntryVidBid(u), offsets[v]);
	//    			low_degree = true;
	//    		} else {
	//    			has_high_degree_neighbor.set_bit_unsync(u);
	//    			++savings;
    //             // u AND v are both high degree vertices, treat the edge specially
	//    		  	if (degrees[u] > high_degree_threshold) {
	//    		  		edge_t edge = edge_t(u, v);
	//    		  		h2h_file.write((char*)&edge, sizeof(edge_t));
	//    		  		++num_h2h_edges;
	//    			}
	//    		}
	//    		if (write_low_degree_edgelist && low_degree) {
	//    			edge_t edge = edge_t(u, v);
	//    			low_degree_file.write((char*)&edge, sizeof(edge_t));
	//    		}
	//    	}
	//    	num_edges -= chunk_size;
	//     if (num_edges < chunk_size) { // adapt chunk size for last batch read
	//     	chunk_size = num_edges;
	//     }
	// }

	// LOG(INFO) << "Edges to a high-degree vertex: " << savings << std::endl;
	// LOG(INFO) << "Edges between two high-degree vertices: " << num_h2h_edges << std::endl;

	// // write the number of vertices and number of low-degree edges to the low-file
	// eid_t num_low_edges = num_edges - num_h2h_edges;

	// low_degree_file.seekp(pos_low);
	// low_degree_file.write((char *)&num_vertices, sizeof(num_vertices));
	// low_degree_file.write((char *)&num_low_edges, sizeof(num_edges));

	return 0;
}

template <>
eid_t Graph<AdjEntryVid>::stream_build(
    std::ifstream &fin, eid_t num_edges, std::vector<vid_t> &degrees
)
{
	fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);

	LOG(INFO) << "builder starts...";
    // to put the in-neighbors at the right position when building the column array
	std::vector<vid_t> offsets(num_vertices, 0); 

	this->num_edges = num_edges; // num_edges, num_vertices
    this->num_vertices = degrees.size();
	double average_degree = num_edges * 2.0 / num_vertices; // non-rounded average degree

	LOG(INFO) << "Average degree: " << average_degree << std::endl;

    std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
	eid_t chunk_size = std::min(num_edges, (eid_t)100000);

	stream_edges.resize(chunk_size);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);

    eid_t num_remaining_edges = num_edges;
    while (num_remaining_edges > 0) {
		fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
	    for (eid_t i = 0; i < chunk_size; ++i) {
	    	++offsets[stream_edges[i].first];
	    }
	    num_remaining_edges -= chunk_size;
	    if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
	    	chunk_size = num_remaining_edges;
	    }
	}

	std::vector<eid_t> index(num_vertices, 0); // for index array

	for (vid_t v = 1; v < num_vertices; ++v) {
        index[v] = index[v - 1] + degrees[v - 1];
	}
    eid_t neighbors_len = index[num_vertices - 1] + degrees[num_vertices - 1];

    neighbors = (AdjEntryVid *)realloc(neighbors, sizeof(AdjEntryVid) * neighbors_len); 
    LOG(INFO) << sizeof(AdjEntryVid) << " bytes needed for AdjEntryVid";
    LOG(INFO) << (sizeof(AdjEntryVid) * neighbors_len / 1024.0 / 1024 / 1024) << " G bytes needed for neighbors";
	CHECK(neighbors) << "allocation failed!";

    // #pragma omp parallel for
    for (vid_t vid = 0; vid < num_vertices; ++vid) {
		vdata[vid] = AdjList(neighbors + index[vid]);
    }

    // start read from beginning
	fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); 

    eid_t edge_id = 0;
    num_remaining_edges = num_edges;
    chunk_size = std::min(num_edges, (eid_t)100000);
    
    while (num_remaining_edges > 0) { 
	   	fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
	    for (eid_t i = 0; i < chunk_size; ++i, ++edge_id) {
	   		const auto& [u, v] = stream_edges[i];
            vdata[u].push_back_out(AdjEntryVid(v));
            vdata[v].push_back_in(AdjEntryVid(u), offsets[v]);
	   	}
    // LOG(INFO) << "FUCK";
	   	num_remaining_edges -= chunk_size;
	    if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
	    	chunk_size = num_remaining_edges;
	    }
	}
    
	return 0;
}

template<typename EdgeHandler>
eid_t iterate_edges(std::ifstream &fin, EdgeHandler edge_handler)
{
    LOG(INFO) << "Streamingly iterating edges (edgelist on disk)";
    eid_t num_edgelist_edges = 0;

    fin.seekg(0, std::ios::beg);

    vid_t num_vertices;
    eid_t num_edges;
    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); 

    std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
    eid_t chunk_size = std::min(num_edges, (eid_t)100000);
    eid_t num_remaining_edges = num_edges; // number of edges to be read from file

    stream_edges.resize(chunk_size);

    while (num_remaining_edges > 0) { // edges to be read
        fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
        for (eid_t i = 0; i < chunk_size; ++i) {
            ++num_edgelist_edges;
            edge_handler(stream_edges[i].first, stream_edges[i].second);
        }

        num_remaining_edges -= chunk_size;
        if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
            chunk_size = num_remaining_edges;
        }
    }
    LOG(INFO) << "num_edgelist_edges: " << num_edgelist_edges;

    fin.close();
    return num_edgelist_edges;
}