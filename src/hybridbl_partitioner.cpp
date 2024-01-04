#include "hybridbl_partitioner.hpp"
#include "conversions.hpp"

HybridBLPartitioner::HybridBLPartitioner(std::string basefilename, bool need_k_split)
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
    // degree_threshold = average_degree * 1;

    fission_occupy.resize(p);
    fusion_occupy.resize(p);
    LOG(INFO) << "average_degree: " << average_degree;
    assigned_edges = 0;
    capacity = static_cast<double>(num_edges) * BALANCE_RATIO / p + 1;
    occupied.assign(p, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    edgelist2bucket.assign(num_edges, kInvalidBid);

    V = dense_bitset(num_vertices);
    super.resize(num_vertices, kInvalidVid);


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
};

void HybridBLPartitioner::calculate_stats()
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

void HybridBLPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();
    
    // V is not empty
    // for (int iter = 0; iter < 50000; ++iter) {
    while (true) {
        if (Q.empty()) {
            if (!get_free_vertex()) {
                break;
            }
            auto vid = free_vertex;
            // LOG(INFO) << "free_vertex: " << vid;
            // if (degrees[vid] < degree_threshold && super[vid] == kInvalidVid) {
            if (adj_in[vid].size() < degree_threshold && super[vid] == kInvalidVid) {
                // InitFusion(vid)
                // LOG(INFO) << "InitFusion: " << vid;
                init_fusion(vid, 0);
                // LOG(INFO) << "End InitFusion: " << vid;

            // } else if (degrees[vid] >= degree_threshold) {
            } else if (adj_in[vid].size() >= degree_threshold) {
                // Fission(vid)
                // LOG(INFO) << "Fission: " << vid;
                fission(vid);
            }
        } else {
            auto [uid, root, dist] = Q.front();
            Q.pop();
            // LOG(INFO) << "Expand to(uid, dist): " << uid << ' ' << dist 
            //         << ", super: " << super[uid];
            if (dist < gamma && super[uid] == kInvalidVid) {
            // if (root_assigned[root] < fusion_threshold && super[uid] == kInvalidVid) {
                // Fusion(uid)
                // LOG(INFO) << "Fusion(uid, dist): " << uid << ' ' << dist;
                fusion(uid, root, dist);
            }
        }
    }

    compute_timer.stop();

    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    LOG(INFO) << assigned_edges << ' ' << num_edges;
    LOG(INFO) << "num_visit(V, E): " << num_visit_vertices << ' ' << num_visit_edges;
    CHECK_EQ(assigned_edges, num_edges);

    // total_time.stop();
    // LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
    eid_t sum = 0;
    for (bid_t b = 0; b < p; ++b) {
        LOG(INFO) << b << ": " << fusion_occupy[b] << '\t' << fission_occupy[b];
        sum += fusion_occupy[b] + fission_occupy[b];
    }
    LOG(INFO) << "sum|num_edges: " << sum << '|' << num_edges;
    eid_t max_assigned = 0;
    for (const auto& [vid, assigned] : root_assigned) {
        // LOG(INFO) << vid << ' ' << assigned;
        max_assigned = std::max(assigned, max_assigned);
    }
    LOG(INFO) << "max_assigned: " << max_assigned;
}

void HybridBLPartitioner::init_fusion(vid_t vid, vid_t dist)
{
    super[vid] = vid;
    bid_t candidate = kInvalidBid;
    eid_t candidate_size = std::numeric_limits<eid_t>::max();
    for (bid_t b = 0; b < p; ++b) {
        if (occupied[b] < candidate_size) {
            candidate = b;
            candidate_size = occupied[b];
        }
    }
    root_bucket[vid] = candidate;
    fusion(vid, vid, dist);
}

void HybridBLPartitioner::fusion(vid_t vid, vid_t root, vid_t dist)
{
    // LOG(INFO) << "Fussion function begin(vid, dist, ind): " 
    //     << vid << ' ' << dist << ' ' 
    //     << adj_in[vid].size() << ' ' << adj_out[vid].size();
    // if (root_assigned[root] >= fusion_threshold) return;
    ++num_visit_vertices;
    V.set(vid, 1);
    vid_t current_super = super[root];
    assert(super[root] == root);
    super[vid] = current_super;
    bid_t current_super_bid = root_bucket[root];
    for (int direction = 0; direction < 2; ++direction) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        for (vid_t i = 0; i < neighbors.size(); ++i) {
            eid_t eid = neighbors[i].v;
            if (edgelist2bucket[eid] != kInvalidBid) continue;
            num_visit_edges += 1;

            vid_t &uid = direction ? edges[eid].second : edges[eid].first;
            // if (super[uid] != kInvalidVid) continue;

            // super[uid] = current_super;
            assign_edge(current_super_bid, vid, uid, eid);
            ++fusion_occupy[current_super_bid];
            ++root_assigned[root];
            Q.push({uid, root, dist + 1});
        }
    }
    // LOG(INFO) << "Fussion function end(vid, dist, ind)" ;
}

void HybridBLPartitioner::fission(vid_t vid)
{
    ++num_visit_vertices;
    V.set(vid, 1);
    for (int direction = 0; direction < 2; ++direction) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        for (vid_t i = 0; i < neighbors.size(); ++i) {
            eid_t eid = neighbors[i].v;
            if (edgelist2bucket[eid] != kInvalidBid) continue;
            num_visit_edges += 1;

            vid_t &uid = direction ? edges[eid].second : edges[eid].first;
            assign_edge(uid % p, vid, uid, eid);
            ++fission_occupy[uid % p];
        }
    }
}
