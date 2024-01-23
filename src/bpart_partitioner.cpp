#include <random>

#include "bpart_partitioner.hpp"
#include "conversions.hpp"
#include <omp.h>

template <typename TAdj>
BPartPartitioner<TAdj>::BPartPartitioner(std::string basefilename, bool need_k_split)
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
    CHECK_EQ(sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t), filesize);

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    // alpha = sqrt(num_partitions) * static_cast<double>(num_edges) / pow(num_vertices, 1.5);

    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);


    vertex2bucket.assign(num_vertices, kInvalidBid);
    // capacity = static_cast<double>(num_edges) * 2 * 1.05 / FLAGS_p + 1; //will be used to as stopping criterion later
    // capacity = 1e18; //will be used to as stopping criterion later

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    average_degree = static_cast<double>(num_edges) * 2.0 / num_vertices;
}

template <typename TAdj>
void BPartPartitioner<TAdj>::split()
{
    vid_t aim_subgraph_vertex = (num_vertices + FLAGS_p - 1) / FLAGS_p;
    eid_t aim_subgraph_edge = (num_edges + FLAGS_p - 1) / FLAGS_p;

    // std::vector<vid_t> order(num_vertices);
    // std::iota(order.begin(), order.end(), (vid_t)0);
    // std::mt19937 engine(123);
    // std::shuffle(order.begin(), order.end(), engine);

    int num_finished_partition = 0;
    int num_remain_partition = FLAGS_p;
    int iter = 1;
    while (num_finished_partition < FLAGS_p) {
        LOG(INFO) << "Iteration " << iter;
        
        num_remain_partition = FLAGS_p - num_finished_partition;
        auto& num_under_partition = num_partitions;
        num_under_partition = (1 << (iter)) * num_remain_partition;
        auto& average_balance = capacity;
        average_balance = static_cast<double>(num_vertices) / num_under_partition * 2.0;

        // fine-grained 
        w_.assign(num_under_partition, 0);
        is_boundarys.assign(num_under_partition, dense_bitset(num_vertices));
        num_bucket_vertices.assign(num_under_partition, 0);
        occupied.assign(num_under_partition, 0);

        for (vid_t v = 0; v < num_vertices; ++v) {
            // vid_t vid = order[v];
            vid_t vid = v;
            // already in fininshed partition
            if (vertex2bucket[vid] != kInvalidBid) {
                continue;
            }
            if (v % 5'000'000 == 0) {
                LOG(INFO) << "Vertex processed " << v;
            }
            // according to bpart scoring
            const auto [bucket, additional_edges] = best_scored_partition(vid); 
            assign_vertex(bucket, vid, additional_edges);
        }

        std::vector<vid_t> num_last_round_vertex(num_bucket_vertices);
        std::vector<eid_t> num_last_round_edge(occupied);

        std::vector<uint32_t> id_last_round(num_under_partition);
        std::iota(id_last_round.begin(), id_last_round.end(), 0);

        std::vector<uint32_t> parent(num_under_partition);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](auto self, uint32_t x) -> uint32_t {
            if (x == parent[x]) {
                return x;
            } else {
                parent[x] = self(self, parent[x]);
                return parent[x];
            }
        };

        auto merge = [&](uint32_t x, uint32_t y) -> uint32_t {
            uint32_t x_root = find(find, x);
            uint32_t y_root = find(find, y);
            if (x_root == y_root) {
                return x_root;
            }
            parent[y_root] = x_root;
            return x_root;
        };

        for (int round = 0; round < iter; ++round) {
            // target bucket num after combination
            int num_target_buckets = num_under_partition / (1 << (round + 1));

            std::vector<vid_t> num_combined_partition_vertex(num_under_partition);
            std::vector<eid_t> num_combined_partition_edge(num_under_partition);

            std::vector<uint32_t> sorted_buckets_id(id_last_round);

            // sort by num_vertices, using the result of last round
            std::sort(
                sorted_buckets_id.begin(), sorted_buckets_id.end(), 
                [&](const auto& a, const auto& b) {
                    return num_last_round_vertex[a] < num_last_round_vertex[b];
                }
            );
            id_last_round.clear();

            LOG(INFO) << "round: " << round << " num_target_buckets: " << num_target_buckets;
            for (size_t i = 0; i < sorted_buckets_id.size() / 2; ++i) {
                bid_t bucket_1 = sorted_buckets_id[i];
                bid_t bucket_2 = sorted_buckets_id[sorted_buckets_id.size() - i - 1];

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
        LOG(INFO) << static_cast<double>(aim_subgraph_vertex) * balance_factor << ' ' << static_cast<double>(aim_subgraph_edge) * balance_factor;
        
        for (size_t bucket = 0; bucket < parent.size(); ++bucket) {
            if (bucket != parent[bucket])   continue;
            double delta_vertex = static_cast<double>(num_last_round_vertex[bucket]) - aim_subgraph_vertex;
            double delta_edge = static_cast<double>(num_last_round_edge[bucket]) - aim_subgraph_edge;

            LOG(INFO) << "bucket: " << bucket << " num_last_round_vertex: " <<
                num_last_round_vertex[bucket] << " delta_vertex: " << delta_vertex << " num_last_round_edge: " << num_last_round_edge[bucket] << " delta_edge: " << delta_edge;
            if ((std::abs(delta_vertex) < static_cast<double>(aim_subgraph_vertex) * balance_factor 
                && std::abs(delta_edge) < static_cast<double>(aim_subgraph_edge) * balance_factor) 
                    or iter == 1 or FLAGS_p - num_finished_partition <= 2) {
                
                // #pragma omp parallel for
                for (size_t b = 0; b < parent.size(); ++b) {
                    if (parent[b] != bucket) continue;
                    for (vid_t vid = 0; vid < num_vertices; ++vid) {
                        if (is_boundarys[b].get(vid)) {
                            vertex2bucket[vid] = num_finished_partition;
                        }
                    }
                }

                LOG(INFO) << "bucket: " << bucket << " is ok!";
                ++num_finished_partition;
            }
        }
        ++iter;
    }
    

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    
    num_partitions = FLAGS_p;
    w_.assign(num_partitions, 0);
    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    num_bucket_vertices.assign(num_partitions, 0);
    occupied.assign(num_partitions, 0);

    for (vid_t vid = 0; vid < num_vertices; ++vid) {
        bid_t bucket = vertex2bucket[vid];
        writer.save_vertex(vid, bucket);
        assign_vertex(bucket, vid, 0);
    }
    for (const auto &[u, v] : edges) {
        bid_t bu = vertex2bucket[u], bv = vertex2bucket[v];
        if (bu != bv) {
            occupied[bu] += 1;
            occupied[bv] += 1;
        } else {
            occupied[bu] += 1;
        }
    }

    calculate_stats();
}

// <final bucket, additional edges>
template <typename TAdj>
std::tuple<bid_t, eid_t> 
BPartPartitioner<TAdj>::best_scored_partition(vid_t vid) 
{
	double best_score = -1e18;
	bid_t best_partition = kInvalidBid;
    eid_t final_additional_edges = 0;
	for (bid_t b = 0; b < num_partitions; ++b) {
        if (w_[b] >= capacity)  continue;
		const auto &[score, additional_edges] = compute_partition_score(vid, b);
        if (score > best_score) {
            best_score = score;
            best_partition = b;
            final_additional_edges = additional_edges;
        }
	}
    if (best_partition == kInvalidBid) {
        best_partition = gen() % num_partitions;
        const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, best_partition);
        final_additional_edges = neighbors_cnt - overlap;
    }
	return {best_partition, final_additional_edges};
}

template <typename TAdj>
std::tuple<vid_t, vid_t>
BPartPartitioner<TAdj>::overlap_partition_vertex(vid_t vid, bid_t bucket_id) 
{
    vid_t overlap = 0, neighbors_cnt = 0;
    for (int direction = 0; direction < 2; ++direction) {
        auto &neighbors = (direction ? adj_out[vid] : adj_in[vid]);
        neighbors_cnt += neighbors.size();
        for (size_t i = 0; i < neighbors.size(); ++i) {
            vid_t uid = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
            if (is_boundarys[bucket_id].get(uid)) {
                ++overlap;
            }
        }
    }
    return {overlap, neighbors_cnt};
}

template <typename TAdj>
double BPartPartitioner<TAdj>::intra_partition_cost(double size) 
{
    return alpha * gamma * pow(size, gamma - 1.0);
}

template <typename TAdj>
std::tuple<double, vid_t> 
BPartPartitioner<TAdj>::compute_partition_score(vid_t vid, bid_t bucket_id) 
{
    double intra_partition_score = intra_partition_cost(w_[bucket_id]);
    const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, bucket_id);
    double inter_partition_score = overlap;
	return {inter_partition_score - intra_partition_score, neighbors_cnt - overlap};
}

template <typename TAdj>
void BPartPartitioner<TAdj>::assign_vertex(bid_t bucket, vid_t vid, eid_t additional_edges)
{
    is_boundarys[bucket].set_bit_unsync(vid);
    
    num_bucket_vertices[bucket] += 1;
    occupied[bucket] += 0;
    w_[bucket] = static_cast<double>(num_bucket_vertices[bucket]) + occupied[bucket] * 2.0 / average_degree;
}
