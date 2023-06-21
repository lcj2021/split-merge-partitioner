#include <utility>

#include "util.hpp"
#include "ebv_partitioner.hpp"
#include "conversions.hpp"

EbvPartitioner::EbvPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), p(p), writer(basefilename, !need_k_split && FLAGS_write)
{
    Timer convert_timer;
    convert_timer.start();
    convert(basefilename, new Converter(basefilename));
    convert_timer.stop();

    LOG(INFO) << "convert time: " << convert_timer.get_time();

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    // In-memory
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

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";

    is_boundarys.assign(p, dense_bitset(num_vertices));
    occupied.assign(p, 0);
    vcount.assign(p, 0);
    avg_edge_cnt = (double)num_edges / FLAGS_p;
    edge2bucket.assign(num_edges, -1);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
}

void EbvPartitioner::split()
{
    // std::shuffle(edges.begin(), edges.end(), rd);
    sort(edges.begin(), edges.end(), [&](const auto &l, const auto &r) {
        vid_t lu = l.first, lv = l.second;
        vid_t ru = r.first, rv = r.second;
        return degrees[lu] + degrees[lv] < degrees[ru] + degrees[rv];
    });

    for (size_t i = 0; i < num_edges; ++ i) {
        edge_t e = edges[i];
        vid_t u = e.first, v = e.second;
        int bucket = best_scored_partition(u, v, i); // according to ebv scoring
        assign_edge(bucket, u, v, i);
        if (i % 50000000 == 0) {
            LOG(INFO) << "Processing edges " << i;
        }
    }

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    for (int i = 0; i < p; ++ i) {
        LOG(INFO) << i << ' ' << is_boundarys[i].popcount() << ' ' << occupied[i];
    }
    calculate_stats();
}

int EbvPartitioner::best_scored_partition(vid_t u, vid_t v, int edge_id) noexcept
{
    double best_score = 1e18;
	int best_partition = 0;
	for (int b = 0; b < p; ++ b) {
		double score = compute_partition_score(u, v, b, edge_id);
		if (score < best_score) {
			best_score = score;
			best_partition = b;
		}
	}
	return best_partition;
}

double EbvPartitioner::compute_partition_score(vid_t u, vid_t v, int bucket_id, int edge_id) noexcept
{
	double su = 0.0, sv = 0.0;
    bool u_is_boundary = is_boundarys[bucket_id].get(u), 
        v_is_boundary = is_boundarys[bucket_id].get(v);
    if (!u_is_boundary) {
        ++ su ;
    } 
    if (!v_is_boundary) {
        ++ sv;
    }

    double imbalance = (double)occupied[bucket_id] / avg_edge_cnt
                    + (double)vcount[bucket_id] / (all_part_vertice_cnt / p);

	double score = (su + sv) + (imbalance);
	return score;
}

void EbvPartitioner::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';

    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    size_t max_part_edge_cnt = 0, all_part_edge_cnt = 0; 
    for (int b = 0; b < p; ++ b) {
        CHECK_EQ(vcount[b], is_boundarys[b].popcount());
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, is_boundarys[b].popcount());
        all_part_vertice_cnt += is_boundarys[b].popcount();
        max_part_edge_cnt = std::max(max_part_edge_cnt, occupied[b]);
        all_part_edge_cnt += occupied[b];
    }
    
    double avg_vertice_cnt = (double)all_part_vertice_cnt / p;
    double avg_edge_cnt = (double)all_part_edge_cnt / p;
    double std_deviation = 0.0;
    for (int b = 0; b < p; ++ b) 
        std_deviation += pow(is_boundarys[b].popcount() - avg_vertice_cnt, 2);
    std_deviation = sqrt((double)std_deviation / p);
    
    LOG(INFO) << "Vertice balance: "
              << (double)max_part_vertice_cnt / ((double)num_vertices / p);
    LOG(INFO) << "Max Vertice count: "
              << max_part_vertice_cnt;
    LOG(INFO) << "Vertice std_deviation / avg: "
              << std_deviation / avg_vertice_cnt;
    LOG(INFO) << "Edge balance: "
              << (double)max_part_edge_cnt / avg_edge_cnt;
    CHECK_EQ(all_part_edge_cnt, num_edges);

    LOG(INFO) << "replication factor (final): " << (double)all_part_vertice_cnt / num_vertices;
}