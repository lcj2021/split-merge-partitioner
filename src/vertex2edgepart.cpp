/*
 * vertex2edgePart.cpp
 *
 *  Created on: June 19, 2023
 *      Author: xxx
 *      Base on: Mayerrn
 *      Description: Use to transform METIS edge cut to vertex cut. 
 *              First get adjacent list by method "e2a". 
 *              Then feed the *.adjlist to METIS to get ".part"
 *              Finally use this "v2e" method to get the vertex cut form. 
 *              
 *              If FLAGS_k = 1, this method will output the original result 
 *              of METIS under FLAGS_p partitions. 
 */

#include <queue>
#include <numeric>

#include "vertex2edgepart.hpp"

Vertex2EdgePart::Vertex2EdgePart(std::string basefilename): rd(), gen(rd()) {
	this->basefilename = basefilename;
    CHECK_NE(FLAGS_p, 0); p = FLAGS_p;
    CHECK_NE(FLAGS_k, 0); k = FLAGS_k;
}

Vertex2EdgePart::~Vertex2EdgePart() {
}

void 
Vertex2EdgePart::read_vertexpart() {
	// open the partitioning file
	// assumption: basefilename + ".part"
    std::string partition_method = FLAGS_method.substr(4);
	std::string partfilename = basefilename + ".vertexpart." + partition_method + "." + std::to_string(p * k); //METIS file name format

    LOG(INFO) << partfilename;

    std::string line;
    std::ifstream partfile(partfilename);
    uint32_t current_vertex = 0;
    int result;
    while(std::getline(partfile, line)) {
        result = std::stoi(line);
        current_vertex++; // as in METIS format, first vertex id is 1, not 0!
        vertex2partition[current_vertex] = result;
        bucket_info[result].occupied ++;
    }
}

int 
Vertex2EdgePart::vertex2edgepartID(vid_t u, vid_t v) {
	int target_u = vertex2partition[u];
	int target_v = vertex2partition[v];
	if (u == v) {
		return target_u;
	} else {
		// flip a coin: edge goes to either u's or v's partition
		if (rand() % 2 == 0) {
			return target_u;
		} else {
			return target_v;
		}
	}
}

void
Vertex2EdgePart::merge() {

    size_t max_part_edges_cnt = 0, all_part_edges_cnt = 0; 
    for (int b = 0; b < p * k; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_mirror.popcount();
        max_part_edges_cnt = std::max(max_part_edges_cnt, bucket_info[b].replicas);
        all_part_edges_cnt += bucket_info[b].replicas;
    }

    std::sort(bucket_info.begin(), bucket_info.end(), 
        [&](const bucket_info_item &l, const bucket_info_item &r) {
        return l.replicas > r.replicas;
    });

    for (int b = 0; b < p * k; ++ b) 
        LOG(INFO)   << "Bucket_info " << bucket_info[b].old_id 
                    << " edges: " << bucket_info[b].replicas 
                    << " vertices: " << bucket_info[b].occupied
                    << " is_chosen: " << bucket_info[b].is_chosen << std::endl;
    LOG(INFO) << "Edge cut ratio: " << (double)(all_part_edges_cnt - num_edges) / num_edges;

    curr_bucket_id = 0;
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >

    // valid_bucket = fast_merge();
    valid_bucket = precise_merge();

    // rearrange bucket after heuristic merging
    std::sort(bucket_info.begin(), bucket_info.end());

    for (size_t b = 0; b < bucket_info.size(); ) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        if (is_chosen) {
            old_id = get_final_bucket(old_id);
            std::cerr << "Is chosen, new old_id: " << old_id << ", replicas: " << replicas << '\n';
            ++ b;
        } else {
            std::swap(bucket_info[b], bucket_info.back());
            bucket_info.pop_back();
        }
    }
    bucket_info.shrink_to_fit();
    for (size_t b = 0; b < bucket_info.size(); ++ b) 
        LOG(INFO) << "Bucket_info " << bucket_info[b].old_id 
                    << ", edges: " << bucket_info[b].replicas 
                    << ", vertices: " << bucket_info[b].occupied
                    << ", rank: " << b;

    size_t curr_assigned_vertices = rearrange_vertice(edges, valid_bucket);
    assigned_vertices += curr_assigned_vertices;
}

std::unordered_map<int, int> 
Vertex2EdgePart::fast_merge()
{
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >
    
    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::priority_queue<std::tuple<int, int, int, int>, std::vector<std::tuple<int, int, int, int>>, std::greater<std::tuple<int, int, int, int>>> pq;  
    for (int b = 0; b < p; ++ b) {
        pq.emplace(0, 0, b, b);
    }
    for (int b = 0; b < p * k; ++ b) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = pq.top();
        pq.pop();
        if (partitions_inside == 0) {
            parent_old_id = get_final_bucket(old_id);
            pq.emplace(replicas, partitions_inside + 1, b, parent_old_id);
            is_chosen = true;
        } else {
            bool has_intersection = false;
            size_t new_vertice_cnt = merge_bucket(parent_bucket, b, has_intersection);
            if (!has_intersection) { 
                std::cerr << "No intersection betwenn bucket " << b << " and " << parent_bucket << '\n';
                // exit(1);
            }
            if (partitions_inside + 1 < k) {
                pq.emplace(new_vertice_cnt, partitions_inside + 1, parent_bucket, parent_old_id);
            }
        }
        valid_bucket.emplace(old_id, parent_old_id);
    }
    return valid_bucket;
}

std::unordered_map<int, int> 
Vertex2EdgePart::precise_merge()
{
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >

    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::vector<std::tuple<int, int, int, int>> final_bucket;
    for (int b = 0; b < p; ++ b) {
        final_bucket.emplace_back(0, 0, b, b);
    }
    auto compute_new_bucket_size = [&](int bid_a, int bid_b) -> size_t {
        const auto &is_mirror_a = bucket_info[bid_a].is_mirror, &is_mirror_b = bucket_info[bid_b].is_mirror;
        dense_bitset new_bucket = is_mirror_a | is_mirror_b;
        return new_bucket.popcount();
    };

    for (int b = 0; b < p * k; ++ b) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];

        int best_final_bucket = -1;
        size_t min_size_after_merge = std::numeric_limits<size_t>::max();

        for (size_t fb = 0; fb < final_bucket.size(); ++ fb) {
            auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[fb];
            if (partitions_inside == k) continue;
            if (partitions_inside == 0) {
                best_final_bucket = fb;
                break;
            }
            size_t new_size_after_merge = compute_new_bucket_size(parent_bucket, b);
            if (min_size_after_merge > new_size_after_merge) {
                min_size_after_merge = new_size_after_merge;
                best_final_bucket = fb;
            }
        }
        CHECK_NE(best_final_bucket, -1);

        auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[best_final_bucket];
        if (partitions_inside == 0) {
            parent_old_id = get_final_bucket(old_id);
            final_bucket[best_final_bucket] = {replicas, partitions_inside + 1, b, parent_old_id};
            is_chosen = true;
        } else {
            bool has_intersection = false;
            size_t new_vertice_cnt = merge_bucket(parent_bucket, b, has_intersection);
            if (!has_intersection) { 
                std::cerr << "No intersection betwenn bucket " << b << " and " << parent_bucket << '\n';
                // exit(1);
            }

            final_bucket[best_final_bucket] = {new_vertice_cnt, partitions_inside + 1, parent_bucket, parent_old_id};
        }
        valid_bucket.emplace(old_id, parent_old_id);
    }
    return valid_bucket;
}

int 
Vertex2EdgePart::merge_bucket(int dst, int src, bool &has_intersection)   // dst, src
{
    auto &is_mirror_a = bucket_info[dst].is_mirror, &is_mirror_b = bucket_info[src].is_mirror;
    size_t mirror_cnt = 0;
    for (size_t i = 0; i < is_mirror_a.size(); ++ i) {
        if (is_mirror_a.get(i) == 1 || is_mirror_b.get(i) == 1) {
            is_mirror_a.set_bit_unsync(i);
            ++ mirror_cnt;
        }
        if (!has_intersection && is_mirror_a.get(i) == 1 && is_mirror_b.get(i) == 1) {
            has_intersection = true;
        }
    }
    bucket_info[dst].replicas = mirror_cnt;
    bucket_info[dst].occupied += bucket_info[src].occupied;
    return mirror_cnt;
}

void 
Vertex2EdgePart::split() {
    if (!split_in_edgelist() && !split_in_adjlist()) {
        LOG(FATAL) << "Read graph failed!";
    }

    merge();
    CHECK_EQ(assigned_vertices, num_vertices);

    for (const auto &edge : edges) {
        vid_t from = edge.first, to = edge.second;
        int target_p = vertex2edgepartID(from, to);
        occupied[target_p]++;
        is_boundarys[target_p].set_bit_unsync(from);
        is_boundarys[target_p].set_bit_unsync(to);
        // from --, to --;
    }

	calculate_stats();
}

void 
Vertex2EdgePart::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<size_t> bucket2vcnt(p, 0);
    rep (i, p) {
        bucket2vcnt[i] = is_boundarys[i].popcount();
    }

    size_t max_part_vertice_cnt = *std::max_element(bucket2vcnt.begin(), bucket2vcnt.end());
    size_t all_part_vertice_cnt = accumulate(bucket2vcnt.begin(), bucket2vcnt.end(), (size_t)0);
    size_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
    size_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (size_t)0);

    rep(i, p)
        LOG(INFO) << "Partition id: " << i
                    << ", vertice: " << bucket2vcnt[i]
                    << ", edges: " << occupied[i] << std::endl;
    
    double avg_vertice_cnt = (double)all_part_vertice_cnt / (p);
    double avg_edge_cnt = (double)all_part_edge_cnt / (p);

    double std_vertice_deviation = 0.0;
    double std_edge_deviation = 0.0;
    for (int b = 0; b < p; ++ b) {
        std_vertice_deviation += pow(bucket2vcnt[b] - avg_vertice_cnt, 2);
        std_edge_deviation += pow(occupied[b] - avg_edge_cnt, 2);
    }
    std_vertice_deviation = sqrt((double)std_vertice_deviation / p);
    std_edge_deviation = sqrt((double)std_edge_deviation / p);
    
    LOG(INFO) << std::string(20, '#') << "\tVertice    balance\t" << std::string(20, '#');
    LOG(INFO) << "Max vertice count / avg vertice count: "
              << (double)max_part_vertice_cnt / ((double)num_vertices / (p));
    LOG(INFO) << "Max Vertice count: "
              << max_part_vertice_cnt;
    LOG(INFO) << "Avg Vertice count(No replicate): "
              << num_vertices / p;
    LOG(INFO) << "Vertice std_vertice_deviation / avg: "
              << std_vertice_deviation / avg_vertice_cnt;

    LOG(INFO) << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#');
    LOG(INFO) << "Max edge count / avg edge count: "
              << (double)max_part_edge_cnt / avg_edge_cnt;
    LOG(INFO) << "Max Edge count: "
              << max_part_edge_cnt;
    LOG(INFO) << "Avg Edge count: "
              << num_edges / p;
    LOG(INFO) << "Edge std_edge_deviation / avg: "
              << std_edge_deviation / avg_edge_cnt;

    CHECK_EQ(all_part_edge_cnt, num_edges);
    LOG(INFO) << std::string(20, '#') << "\tReplicate    factor\t" << std::string(20, '#');
    LOG(INFO) << "replication factor (final): " << (double)all_part_vertice_cnt / num_vertices;
}