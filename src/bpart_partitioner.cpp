#include <random>

#include "bpart_partitioner.hpp"
#include "conversions.hpp"
#include <omp.h>

BPartPartitioner::BPartPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename, FLAGS_write) // Vertex partitioner must record basefilename.part 
{
    Timer convert_timer;
    convert_timer.start();
    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter);
    delete converter;
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time();

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);

    p = FLAGS_p;
    if (need_k_split) {
        p *= FLAGS_k;
    }

    // alpha = sqrt(p) * (double)num_edges / pow(num_vertices, 1.5);

    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);


    vertex2bucket.assign(num_vertices, -1);
    // capacity = (double)num_edges * 2 * 1.05 / FLAGS_p + 1; //will be used to as stopping criterion later
    // capacity = 1e18; //will be used to as stopping criterion later

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);

    // edge2bucket.assign(num_edges, -1);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    average_degree = (double)num_edges * 2.0 / num_vertices;
    // max_degree = *std::max_element(degrees.begin(), degrees.end());
}

void BPartPartitioner::split()
{
    vid_t aim_subgraph_vertex = num_vertices / FLAGS_p;
    size_t aim_subgraph_edge = num_edges / FLAGS_p;

    std::vector<vid_t> order(num_vertices);
    std::iota(order.begin(), order.end(), (vid_t)0);
    // std::mt19937 engine(123);
    // std::shuffle(order.begin(), order.end(), engine);

    int num_finished_partition = 0;
    int num_remain_partition = FLAGS_p;
    int iter = 1;
    while (num_finished_partition < FLAGS_p) {
        LOG(INFO) << "Iteration " << iter;
        
        num_remain_partition = FLAGS_p - num_finished_partition;
        int& num_under_partition = p;
        num_under_partition = (1 << (iter)) * num_remain_partition;
        auto& average_balance = capacity;
        average_balance = (double)num_vertices / num_under_partition * 2.0;

        w_.assign(num_under_partition, 0);
        is_boundarys.assign(num_under_partition, dense_bitset(num_vertices));

        // fine-grained 
        num_vertices_bucket.assign(num_under_partition, 0);
        num_edges_bucket.assign(num_under_partition, 0);

        for (vid_t v = 0; v < num_vertices; ++ v) {
            vid_t vid = order[v];
            // already in fininshed partition
            if (vertex2bucket[vid] != (uint16_t)(-1)) {
                continue;
            }
            if (v % 5'000'000 == 0) {
                LOG(INFO) << "Vertex processed " << v;
            }
            // according to bpart scoring
            const auto [bucket, additional_edges] = best_scored_partition(vid); 
            assign_vertex(bucket, vid, additional_edges);
        }

        std::vector<vid_t> num_last_round_vertex(num_vertices_bucket);
        std::vector<size_t> num_last_round_edge(num_edges_bucket);

        std::vector<int> id_last_round(num_under_partition);
        std::iota(id_last_round.begin(), id_last_round.end(), 0);

        std::vector<size_t> parent(num_under_partition);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](auto self, size_t x) -> size_t {
            if (x == parent[x]) {
                return x;
            } else {
                parent[x] = self(self, parent[x]);
                return parent[x];
            }
        };

        auto merge = [&](size_t x, size_t y) -> size_t {
            size_t x_root = find(find, x);
            size_t y_root = find(find, y);
            if (x_root == y_root) {
                return x_root;
            }
            parent[y_root] = x_root;
            return x_root;
        };

        for (int round = 0; round < iter; ++ round) {
            // target bucket num after combination
            int num_target_buckets = num_under_partition / (1 << (round + 1));

            std::vector<vid_t> num_combined_partition_vertex(num_under_partition);
            std::vector<size_t> num_combined_partition_edge(num_under_partition);

            std::vector<int> sorted_buckets_id(id_last_round);
            // std::iota(sorted_buckets_id.begin(), sorted_buckets_id.end(), 0);

            // sort by num_vertices, using the result of last round
            std::sort(sorted_buckets_id.begin(), sorted_buckets_id.end(), [&](const int &a, const int &b) {
                return num_last_round_vertex[a] < num_last_round_vertex[b];
            });
            id_last_round.clear();

            // combine buckets
            std::vector<int> combined_buckets(num_target_buckets);
            std::iota(combined_buckets.begin(), combined_buckets.end(), 0);

            LOG(INFO) << "round: " << round << " num_target_buckets: " << num_target_buckets;
            for (size_t i = 0; i < sorted_buckets_id.size() / 2; ++ i) {
                int bucket_1 = sorted_buckets_id[i];
                int bucket_2 = sorted_buckets_id[sorted_buckets_id.size() - i - 1];

                if (bucket_1 > bucket_2) std::swap(bucket_1, bucket_2);

                num_combined_partition_vertex[bucket_1] = 
                    num_last_round_vertex[bucket_1] + num_last_round_vertex[bucket_2];
                num_combined_partition_edge[bucket_1] = 
                    num_last_round_edge[bucket_1] + num_last_round_edge[bucket_2];

                auto new_root = merge(bucket_1, bucket_2);
                id_last_round.emplace_back(new_root);
            }
            for (auto& b : parent) {
                b = find(find, b);
            }

            std::swap(num_last_round_vertex, num_combined_partition_vertex);
            std::swap(num_last_round_edge, num_combined_partition_edge);
            num_last_round_vertex.shrink_to_fit();
            num_last_round_edge.shrink_to_fit();
        }

        LOG(INFO) << "aim_subgraph_vertex: " << aim_subgraph_vertex << " aim_subgraph_edge: " <<
            aim_subgraph_edge;
        LOG(INFO) << (double)aim_subgraph_vertex * balance_factor << ' ' << (double)aim_subgraph_edge * balance_factor;
        
        for (size_t bucket = 0; bucket < parent.size(); ++bucket) {
            if (bucket != parent[bucket])   continue;
            double delta_vertex = double(num_last_round_vertex[bucket]) - double(aim_subgraph_vertex);
            double delta_edge = (double)num_last_round_edge[bucket] - double(aim_subgraph_edge);

            LOG(INFO) << "bucket: " << bucket << " num_last_round_vertex: " <<
                num_last_round_vertex[bucket] << " delta_vertex: " << delta_vertex << " num_last_round_edge: " << num_last_round_edge[bucket] << " delta_edge: " << delta_edge;
            if ((std::abs(delta_vertex) < (double)aim_subgraph_vertex * balance_factor 
                && std::abs(delta_edge) < (double)aim_subgraph_edge * balance_factor) 
                    or iter == 1 or FLAGS_p - num_finished_partition <= 2) {
                
                #pragma omp parallel for
                for (size_t b = 0; b < parent.size(); ++b) {
                    if (parent[b] != bucket) continue;
                    for (vid_t vid = 0; vid < num_vertices; ++vid) {
                        if (is_boundarys[b].get(vid)) {
                            vertex2bucket[vid] = num_finished_partition;
                        }
                    }
                }

                LOG(INFO) << "bucket: " << bucket << " is ok!";
                ++ num_finished_partition;
            }
        }
        ++ iter;
    }
    

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    
    num_vertices_bucket.assign(FLAGS_p, 0);
    num_edges_bucket.assign(FLAGS_p, 0);

    for (vid_t vid = 0; vid < num_vertices; ++vid) {
        int bucket = vertex2bucket[vid];
        writer.save_vertex(vid, bucket);
        num_vertices_bucket[bucket] += 1;
    }
    size_t total_cut_edges = 0;
    for (const auto &[u, v] : edges) {
        int bu = vertex2bucket[u], bv = vertex2bucket[v];
        if (bu != bv) {
            ++ total_cut_edges;
            num_edges_bucket[bu] += 1;
            num_edges_bucket[bv] += 1;
        } else {
            num_edges_bucket[bu] += 1;
        }
    }
    LOG(INFO) << "edge_cut_ratio: " << (double)(total_cut_edges) / num_edges;

    calculate_stats();
}

// <final bucket, additional edges>
std::tuple<int, size_t> BPartPartitioner::best_scored_partition(vid_t vid) 
{
	double best_score = -1e18;
	int best_partition = -1;
    size_t final_additional_edges = 0;
	for (int i = 0; i < p; ++ i) {
        if (w_[i] >= capacity)  continue;
		const auto &[score, additional_edges] = compute_partition_score(vid, i);
		if (score > best_score) {
			best_score = score;
			best_partition = i;
            final_additional_edges = additional_edges;
		}
	}
    if (best_partition == -1) {
        best_partition = gen() % p;
        const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, best_partition);
        final_additional_edges = neighbors_cnt - overlap;
    }
	return {best_partition, final_additional_edges};
}

std::tuple<size_t, size_t>
BPartPartitioner::overlap_partition_vertex(vid_t vid, int bucket_id) 
{
    size_t overlap = 0, neighbors_cnt = 0;
    rep (direction, 2) {
        auto &neighbors = (direction ? adj_out[vid] : adj_in[vid]);
        neighbors_cnt += neighbors.size();
        for (size_t i = 0; i < neighbors.size(); ++ i) {
            vid_t uid = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
            if (is_boundarys[bucket_id].get(uid)) {
                ++ overlap;
            }
        }
    }
    return {overlap, neighbors_cnt};
}

std::tuple<double, size_t> BPartPartitioner::compute_partition_score(vid_t vid, int bucket_id) 
{
    double intra_partition_score = intra_partition_cost(w_[bucket_id]);
    const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, bucket_id);
    double inter_partition_score = overlap;
	return {inter_partition_score - intra_partition_score, neighbors_cnt - overlap};
}

void BPartPartitioner::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    
    p = FLAGS_p;

    size_t max_part_vertice_cnt = *std::max_element(num_vertices_bucket.begin(), num_vertices_bucket.end());
    size_t all_part_vertice_cnt = accumulate(num_vertices_bucket.begin(), num_vertices_bucket.end(), (size_t)0);
    size_t max_part_edge_cnt = *std::max_element(num_edges_bucket.begin(), num_edges_bucket.end());
    size_t all_part_edge_cnt = accumulate(num_edges_bucket.begin(), num_edges_bucket.end(), (size_t)0);

    double avg_vertice_cnt = (double)all_part_vertice_cnt / (p);
    double avg_edge_cnt = (double)all_part_edge_cnt / (p);

    double std_vertice_deviation = 0.0;
    double std_edge_deviation = 0.0;
    for (int b = 0; b < p; ++ b) {
        std_vertice_deviation += pow(num_vertices_bucket[b] - avg_vertice_cnt, 2);
        std_edge_deviation += pow(num_edges_bucket[b] - avg_edge_cnt, 2);
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
    LOG(INFO) << "Vertice jains_fairness: "
              << jains_fairness(num_vertices_bucket);

    LOG(INFO) << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#');
    LOG(INFO) << "Max edge count / avg edge count: "
              << (double)max_part_edge_cnt / avg_edge_cnt;
    LOG(INFO) << "Max Edge count: "
              << max_part_edge_cnt;
    LOG(INFO) << "Avg Edge count: "
              << num_edges / p;
    LOG(INFO) << "Edge std_edge_deviation / avg: "
              << std_edge_deviation / avg_edge_cnt;
    LOG(INFO) << "Edge jains_fairness: "
              << jains_fairness(num_edges_bucket);


    std::vector<int> sorted_buckets(FLAGS_p);
    std::iota(sorted_buckets.begin(), sorted_buckets.end(), 0);
    std::sort(sorted_buckets.begin(), sorted_buckets.end(), [&](const int &a, const int &b) {
        return num_vertices_bucket[a] < num_vertices_bucket[b];
    });
    LOG(INFO) << "partition |V| w |E|";
    for (const auto& i : sorted_buckets) {
        LOG(INFO) << i << "\t" << num_vertices_bucket[i] << "\t" << w_[i] << "\t" << num_edges_bucket[i];
    }
}