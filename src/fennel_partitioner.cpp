#include <random>

#include "fennel_partitioner.hpp"
#include "conversions.hpp"

FennelPartitioner::FennelPartitioner(std::string basefilename, bool need_k_split)
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

    alpha = sqrt(p) * (double)num_edges / pow(num_vertices, 1.5);
    // alpha = 1.5;

    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_boundarys.assign(p, dense_bitset(num_vertices));
    occupied.assign(p, 0);
    vcount.assign(p, 0);
    vertex2bucket.assign(num_vertices, -1);
    // capacity = (double)num_edges * 2 * 1.05 / p + 1; //will be used to as stopping criterion later
    capacity = (double)num_vertices * 1.1 / p + 1; //will be used to as stopping criterion later
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
    // max_degree = *std::max_element(degrees.begin(), degrees.end());
}

void FennelPartitioner::split()
{
    std::vector<vid_t> order(num_vertices);
    std::iota(order.begin(), order.end(), (vid_t)0);
    // std::mt19937 engine(123);
    // std::shuffle(order.begin(), order.end(), engine);

    for (vid_t v = 0; v < num_vertices; ++ v) {
        vid_t vid = order[v];
        if (v % 5'000'000 == 0) {
            LOG(INFO) << "Vertex processed " << v;
        }
        const auto [bucket, additional_edges] = best_scored_partition(vid); // according to ebv scoring
        assign_vertex(bucket, vid, additional_edges);
    }

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    size_t total_cut_edges = 0;
    for (int i = 0; i < p; ++ i) {
        LOG(INFO) << i << "\t" << is_boundarys[i].popcount() << "\t" << vcount[i] << "\t" << occupied[i];
        total_cut_edges += occupied[i];
    }
    double edge_cut_ratio = (double)(total_cut_edges - num_edges) / num_edges;
    LOG(INFO) << "edge_cut_ratio: " << edge_cut_ratio;
    // calculate_stats();
    total_cut_edges = 0;
    for (const auto &[u, v] : edges) {
        if (vertex2bucket[u] != vertex2bucket[v]) {
            ++ total_cut_edges;
        }
    }
    LOG(INFO) << "edge_cut_ratio: " << (double)(total_cut_edges) / num_edges;
}

// <final bucket, additional edges>
std::tuple<int, size_t> FennelPartitioner::best_scored_partition(vid_t vid) {
	double best_score = -1e18;
	int best_partition = -1;
    size_t final_additional_edges = 0;
	for (int i = 0; i < p; ++ i) {
        if (vcount[i] >= capacity)  continue;
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
FennelPartitioner::overlap_partition_vertex(vid_t vid, int bucket_id) {
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

std::tuple<double, size_t> FennelPartitioner::compute_partition_score(vid_t vid, int bucket_id) {
    double intra_partition_score = intra_partition_cost(vcount[bucket_id]);
    const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, bucket_id);
    double inter_partition_score = overlap;
	return {inter_partition_score - intra_partition_score, neighbors_cnt - overlap};
}

void FennelPartitioner::calculate_stats()
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

    for (int b = 0; b < p; ++ b) 
        LOG(INFO) << "Bucket_info: " << b
                << ", vertices: " << bucket2vcnt[b]
                << ", edges: " << occupied[b];
    
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