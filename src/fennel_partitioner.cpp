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

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    alpha = sqrt(num_partitions) * (double)num_edges / pow(num_vertices, 1.5);
    // alpha = 1.5;

    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    occupied.assign(num_partitions, 0);
    w_.assign(num_partitions, 0);
    // vertex2bucket.assign(num_vertices, kInvalidBid);
    // capacity = (double)num_edges * 2 * 1.05 / num_partitions + 1; //will be used to as stopping criterion later
    capacity = (double)num_vertices * 1.1 / num_partitions + 1; // will be used to as stopping criterion later

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
