#include <numeric>
#include <random>

#include "ebv_partitioner.hpp"
#include "conversions.hpp"

EbvPartitioner::EbvPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename, !need_k_split && FLAGS_write)
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
    CHECK_EQ(sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t), filesize);

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";

    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    occupied.assign(num_partitions, 0);
    num_bucket_vertices.assign(num_partitions, 0);
    avg_edge_cnt = (double)num_edges / FLAGS_p;
    // edgelist2bucket.assign(num_edges, kInvalidBid);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
}

void EbvPartitioner::split()
{
    partition_time.start();
    // std::shuffle(edges.begin(), edges.end(), rd);
    sort(edges.begin(), edges.end(), [&](const auto &l, const auto &r) {
        vid_t lu = l.first, lv = l.second;
        vid_t ru = r.first, rv = r.second;
        return degrees[lu] + degrees[lv] < degrees[ru] + degrees[rv];
    });

    for (eid_t eid = 0; eid < num_edges; ++eid) {
        edge_t e = edges[eid];
        vid_t u = e.first, v = e.second;
        bid_t bucket = best_scored_partition(u, v, eid); // according to ebv scoring
        assign_edge(bucket, u, v, eid);
        if (eid % 50000000 == 0) {
            LOG(INFO) << "Processing edges " << eid;
        }
    }
    partition_time.stop();

    total_time.stop();
    LOG(INFO) << "partition time: " << partition_time.get_time();
    calculate_stats();
}

bid_t EbvPartitioner::best_scored_partition(vid_t u, vid_t v, eid_t edge_id)
{
    double best_score = 1e18;
	bid_t best_partition = kInvalidBid;
	for (bid_t b = 0; b < num_partitions; ++b) {
		double score = compute_partition_score(u, v, b, edge_id);
		if (score < best_score) {
			best_score = score;
			best_partition = b;
		}
	}
    if (best_partition == kInvalidBid) {
        best_partition = gen() % num_partitions;
    }
	return best_partition;
}

double EbvPartitioner::compute_partition_score(vid_t u, vid_t v, bid_t bucket_id, eid_t edge_id)
{
	double su = 0.0, sv = 0.0;
    bool u_is_boundary = is_boundarys[bucket_id].get(u), 
        v_is_boundary = is_boundarys[bucket_id].get(v);
    if (!u_is_boundary) {
        ++su;
    } 
    if (!v_is_boundary) {
        ++sv;
    }

    double imbalance = (double)occupied[bucket_id] / avg_edge_cnt
                    + (double)num_bucket_vertices[bucket_id] / (num_vertices_all_buckets / num_partitions);

	double score = (su + sv) + (imbalance);
	return score;
}
