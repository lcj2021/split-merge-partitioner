#include "ne_partitioner.hpp"
#include "conversions.hpp"

NePartitioner::NePartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd())
    , need_k_split(need_k_split), writer(basefilename, !need_k_split && FLAGS_write)
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

    average_degree = num_edges * 2.0 / num_vertices;
    assigned_edges = 0;
    capacity = static_cast<double>(num_edges) * BALANCE_RATIO / p + 1;
    occupied.assign(p, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    edgelist2bucket.assign(num_edges, kInvalidBid);

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";
    LOG(INFO) << sizeof(edge_t) * num_edges / 1024.0 / 1024 / 1024 << " G bytes needed for edges " << __FILE__<<":"<<__LINE__;
    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();

    size_t num_out_edges = 0;
    for (vid_t vid = 0; vid < num_vertices; ++vid) {
        num_out_edges += adj_out[vid].size();
    }
    LOG(INFO) << "num_out_edges: " << num_out_edges;
    LOG(INFO) << "num_edges: " << num_edges;
    CHECK_EQ(num_out_edges, num_edges);
};

void NePartitioner::assign_remaining()
{
    auto &is_boundary = is_boundarys[p - 1], &is_core = is_cores[p - 1];
    for (vid_t u = 0; u < num_vertices; ++u) {
        for (auto &i : adj_out[u]) {
            if (edges[i.v].valid()) {
                assign_edge(p - 1, u, edges[i.v].second, i.v);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }
        }
    }

    for (vid_t i = 0; i < num_vertices; ++i) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            for (bid_t b = 0; b < p - 1; ++b) {
                if (is_cores[b].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
            }
        }
    }
}

void NePartitioner::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<vid_t> num_bucket_vertices(p, 0);
    for (bid_t b = 0; b < p; ++b) {
        num_bucket_vertices[b] = is_boundarys[b].popcount();
    }
    vid_t max_part_vertice_cnt = *std::max_element(num_bucket_vertices.begin(), num_bucket_vertices.end());
    vid_t all_part_vertice_cnt = accumulate(num_bucket_vertices.begin(), num_bucket_vertices.end(), (vid_t)0);
    eid_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
    eid_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (eid_t)0);

    for (bid_t b = 0; b < p; ++b) {
        LOG(INFO) << "Bucket_info: " << b
                << ", vertices: " << num_bucket_vertices[b]
                << ", edges: " << occupied[b];
    }
    
    double avg_vertice_cnt = static_cast<double>(all_part_vertice_cnt) / p;
    double avg_edge_cnt = static_cast<double>(all_part_edge_cnt) / p;

    double std_vertice_deviation = 0.0;
    double std_edge_deviation = 0.0;
    for (bid_t b = 0; b < p; ++b) {
        std_vertice_deviation += pow(num_bucket_vertices[b] - avg_vertice_cnt, 2);
        std_edge_deviation += pow(occupied[b] - avg_edge_cnt, 2);
    }
    std_vertice_deviation = sqrt(static_cast<double>(std_vertice_deviation) / p);
    std_edge_deviation = sqrt(static_cast<double>(std_edge_deviation) / p);

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

void NePartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();
    for (bucket = 0; bucket < p - 1; ++bucket) {
        std::cerr << bucket << ", ";
        DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (occupied[bucket] < capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) {
                if (need_k_split) {
                    if (!get_free_vertex_by_rand(vid)) {
                        DLOG(INFO) << "partition " << bucket
                                << " stop: no free vertices";
                        break;
                    }
                } else {
                    if (!get_free_vertex(vid)) {
                        DLOG(INFO) << "partition " << bucket
                                << " stop: no free vertices";
                        break;
                    }
                }
                
                d = adj_out[vid].size() + adj_in[vid].size();
            } else {
                min_heap.remove(vid);
                /* CHECK_EQ(d, adj_out[vid].size() + adj_in[vid].size()); */
            }

            occupy_vertex(vid, d);
        }
        min_heap.clear();
        for (int direction = 0; direction < 2; ++direction) {
            for (vid_t vid = 0; vid < num_vertices; ++vid) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (vid_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        ++i;
                    } else {
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
        }
    }
    bucket = p - 1;
    std::cerr << bucket << std::endl;
    assign_remaining();
    compute_timer.stop();

    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();

    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
}
