#include <utility>

#include "util.hpp"
#include "hdrf_partitioner.hpp"
#include "conversions.hpp"

HdrfPartitioner::HdrfPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), writer(basefilename, !need_k_split && FLAGS_write)
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

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";

    is_boundarys.assign(p, dense_bitset(num_vertices));
    occupied.assign(p, 0);
    vcount.assign(p, 0);
    capacity = (double)num_edges * 1.0 / p + 1; //will be used to as stopping criterion later
    edge2bucket.assign(num_edges, -1);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    max_degree = *std::max_element(degrees.begin(), degrees.end());
}

void HdrfPartitioner::split()
{
    // std::shuffle(edges.begin(), edges.end(), rd);
    for (size_t i = 0; i < num_edges; ++ i) {
        edge_t e = edges[i];
        if (i % 50000000 == 0) {
            LOG(INFO) << "Processing edges " << i;
        }
        vid_t u = e.first, v = e.second;
        int bucket = best_scored_partition(u, v); // according to ebv scoring
        assign_edge(bucket, u, v, i);
        if (occupied[bucket] > max_size){
            max_size = occupied[bucket];
        }
        if (occupied[bucket] == min_size){
            int min_sized_bucket_count = 0;
            for (int i = 0; i < p; i++){
                if (occupied[i] == min_size){
                    min_sized_bucket_count++;
                }
            }
            if (min_sized_bucket_count == 1){
                min_size++;
            }
        }
    }

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    for (int i = 0; i < p; ++ i) {
        LOG(INFO) << i << ' ' << is_boundarys[i].popcount() << ' ' << occupied[i];
    }
    calculate_stats();
}

int HdrfPartitioner::best_scored_partition(vid_t u, vid_t v) {
	double best_score = -1.0;
	int best_partition = 0;
	for (int i = 0; i < p; i++){
		double score = compute_partition_score(u, v, i);
//		cout << "score for partition " << i << " is " << score << endl;
		if (score > best_score){
			best_score = score;
			best_partition = i;
		}
	}
	return best_partition;
}

double HdrfPartitioner::compute_partition_score(vid_t u, vid_t v, int bucket_id) {
	if (occupied[bucket_id] >= capacity){
//		cout << "partition " << bucket_id << " is full with " << occupied[bucket_id] << endl;
		return -1.0; // partition is full, do not choose it
	}
	size_t degree_u = degrees[u];
	size_t degree_v = degrees[v];
	size_t sum = degree_u + degree_v;
	double gu = 0.0, gv = 0.0;
	if (is_boundarys[bucket_id].get(u)){
		gu = degree_u;
		gu /=sum;
		gu = 1 + (1-gu);
	}
	if (is_boundarys[bucket_id].get(v)){
		 gv = degree_v;
		 gv /= sum;
		 gv = 1 + (1-gv);
	}

	double bal = (max_size - occupied[bucket_id]) / (1.0 + max_size - min_size);

	double score = gu + gv + lambda * bal;
	return score;
}

void HdrfPartitioner::calculate_stats()
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