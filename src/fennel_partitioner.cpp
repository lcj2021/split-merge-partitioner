#include <random>

#include "fennel_partitioner.hpp"
#include "conversions.hpp"

template <typename TAdj>
FennelPartitioner<TAdj>::FennelPartitioner(std::string basefilename, bool need_k_split)
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
    w_.assign(p, 0);
    vertex2bucket.assign(num_vertices, kInvalidBid);
    // capacity = (double)num_edges * 2 * 1.05 / p + 1; //will be used to as stopping criterion later
    capacity = (double)num_vertices * 1.1 / p + 1; //will be used to as stopping criterion later
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
}

template <typename TAdj>
void FennelPartitioner<TAdj>::split()
{
    // std::vector<vid_t> order(num_vertices);
    // std::iota(order.begin(), order.end(), (vid_t)0);
    // std::mt19937 engine(123);
    // std::shuffle(order.begin(), order.end(), engine);

    for (vid_t v = 0; v < num_vertices; ++v) {
        // vid_t vid = order[v];
        vid_t vid = v;
        if (v % 5'000'000 == 0) {
            LOG(INFO) << "Vertex processed " << v;
        }
        const auto [bucket, additional_edges] = best_scored_partition(vid); // according to ebv scoring
        assign_vertex(bucket, vid, additional_edges);
    }

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
}

// <final bucket, additional edges>
template <typename TAdj>
std::tuple<bid_t, eid_t> FennelPartitioner<TAdj>::best_scored_partition(vid_t vid) 
{
	double best_score = -1e18;
	bid_t best_partition = kInvalidBid;
    eid_t final_additional_edges = 0;
	for (bid_t b = 0; b < p; ++b) {
        if (w_[b] >= capacity)  continue;
		const auto &[score, additional_edges] = compute_partition_score(vid, b);
		if (score > best_score) {
			best_score = score;
			best_partition = b;
            final_additional_edges = additional_edges;
		}
	}
    if (best_partition == kInvalidBid) {
        best_partition = gen() % p;
        const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, best_partition);
        final_additional_edges = neighbors_cnt - overlap;
    }
	return {best_partition, final_additional_edges};
}

template <typename TAdj>
std::tuple<vid_t, vid_t> 
FennelPartitioner<TAdj>::overlap_partition_vertex(vid_t vid, bid_t bucket_id) 
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
std::tuple<double, size_t> 
FennelPartitioner<TAdj>::compute_partition_score(vid_t vid, bid_t bucket_id) 
{
    double intra_partition_score = intra_partition_cost(w_[bucket_id]);
    const auto &[overlap, neighbors_cnt] = overlap_partition_vertex(vid, bucket_id);
    double inter_partition_score = overlap;
	return {inter_partition_score - intra_partition_score, neighbors_cnt - overlap};
}

template <typename TAdj>
void FennelPartitioner<TAdj>::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<vid_t> num_bucket_vertices(p, 0);
    eid_t total_cut_edges = 0;
    for (bid_t b = 0; b < p; ++b) {
        num_bucket_vertices[b] = is_boundarys[b].popcount();
        total_cut_edges += occupied[b];
    }
    total_cut_edges -= num_edges;
    vid_t max_part_vertice_cnt = *std::max_element(num_bucket_vertices.begin(), num_bucket_vertices.end());
    vid_t all_part_vertice_cnt = accumulate(num_bucket_vertices.begin(), num_bucket_vertices.end(), (vid_t)0);
    eid_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
    eid_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (eid_t)0);

    for (bid_t b = 0; b < p; ++b) 
        LOG(INFO) << "Bucket_info: " << b
                << ", vertices: " << num_bucket_vertices[b]
                << ", edges: " << occupied[b];
    
    double avg_vertice_cnt = (double)all_part_vertice_cnt / (p);
    double avg_edge_cnt = (double)all_part_edge_cnt / (p);

    double std_vertice_deviation = 0.0;
    double std_edge_deviation = 0.0;
    for (int b = 0; b < p; ++ b) {
        std_vertice_deviation += pow(num_bucket_vertices[b] - avg_vertice_cnt, 2);
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

    CHECK_EQ(all_part_vertice_cnt, num_vertices);

    LOG(INFO) << std::string(20, '#') << "\tEdge  cut  ratio\t" << std::string(20, '#');
    LOG(INFO) << "Edge cut ratio (final): " << (double)total_cut_edges / num_edges;
}