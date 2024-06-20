#include "hdrf_partitioner.hpp"
#include "conversions.hpp"

HdrfPartitioner::HdrfPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename)
{
    if (need_k_split || FLAGS_write == "none") {
        writer = std::make_unique<EdgepartWriterBase<vid_t, bid_t>>(basefilename);
    } else {
        if (FLAGS_write == "onefile") {
            writer = std::make_unique<EdgepartWriterOnefile<vid_t, bid_t>>(basefilename);
        } else if (FLAGS_write == "multifile") {
            writer = std::make_unique<EdgepartWriterMultifile<vid_t, bid_t>>(basefilename);
        }
    }
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

    LOG(INFO) << "constructing...";

    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    occupied.assign(num_partitions, 0);
    capacity = (double)num_edges * 1.0 / num_partitions + 1; //will be used to as stopping criterion later
    // edgelist2bucket.assign(num_edges, kInvalidBid);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    max_degree = *std::max_element(degrees.begin(), degrees.end());
}

void HdrfPartitioner::split()
{
    partition_time.start();
    // std::shuffle(edges.begin(), edges.end(), rd);
    // for (eid_t eid = 0; eid < num_edges; ++eid) {
    //     edge_t e = edges[eid];
    //     if (eid % 50000000 == 0) {
    //         LOG(INFO) << "Processing edges " << eid;
    //     }
    //     vid_t u = e.first, v = e.second;
    //     bid_t bucket = best_scored_partition(u, v); // according to ebv scoring
    //     assign_edge(bucket, u, v, eid);
    //     if (occupied[bucket] > max_size) {
    //         max_size = occupied[bucket];
    //     }
    //     if (occupied[bucket] == min_size) {
    //         int min_sized_bucket_count = 0;
    //         for (bid_t b = 0; b < num_partitions; ++b) {
    //             if (occupied[b] == min_size) {
    //                 ++min_sized_bucket_count;
    //             }
    //         }
    //         if (min_sized_bucket_count == 1) {
    //             ++min_size;
    //         }
    //     }
    // }

    std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
    eid_t chunk_size = std::min(num_edges, (eid_t)100000);
    eid_t num_remaining_edges = num_edges; // number of edges to be read from file

    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); 

    stream_edges.resize(chunk_size);

    eid_t eid = 0; // number of edges read from file
    while (num_remaining_edges > 0) { // edges to be read
        fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
        for (eid_t i = 0; i < chunk_size; ++i, ++eid) {
            const auto& [u, v] = stream_edges[i];
            if (eid % 50000000 == 0) {
                LOG(INFO) << "Processing edges " << eid;
            }
            bid_t bucket = best_scored_partition(u, v); // according to ebv scoring
            assign_edge(bucket, u, v, eid);
            if (occupied[bucket] > max_size) {
                max_size = occupied[bucket];
            }
            if (occupied[bucket] == min_size) {
                int min_sized_bucket_count = 0;
                for (bid_t b = 0; b < num_partitions; ++b) {
                    if (occupied[b] == min_size) {
                        ++min_sized_bucket_count;
                    }
                }
                if (min_sized_bucket_count == 1) {
                    ++min_size;
                }
            }
        }

        num_remaining_edges -= chunk_size;
        if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
            chunk_size = num_remaining_edges;
        }
    }

    partition_time.stop();
    total_time.stop();
    LOG(INFO) << "partition time: " << total_time.get_time();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
}

bid_t HdrfPartitioner::best_scored_partition(vid_t u, vid_t v) 
{
	double best_score = -1.0;
	bid_t best_partition = 0;
	for (bid_t b = 0; b < num_partitions; ++b) {
		double score = compute_partition_score(u, v, b);
		if (score > best_score) {
			best_score = score;
			best_partition = b;
		}
	}
	return best_partition;
}

double HdrfPartitioner::compute_partition_score(vid_t u, vid_t v, bid_t bucket_id) 
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
