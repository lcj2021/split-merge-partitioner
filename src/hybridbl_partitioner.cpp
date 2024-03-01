#include "hybridbl_partitioner.hpp"
#include "conversions.hpp"

template <typename TAdj>
HybridBLPartitioner<TAdj>::HybridBLPartitioner(std::string basefilename, bool need_k_split)
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

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    average_degree = num_edges * 2.0 / num_vertices;
    // degree_threshold = average_degree * 1;

    // fission_occupy.resize(num_partitions);
    // fusion_occupy.resize(num_partitions);
    LOG(INFO) << "average_degree: " << average_degree;
    LOG(INFO) << "degree_threshold: " << degree_threshold;
    LOG(INFO) << "gamma: " << gamma;
    assigned_edges = 0;
    capacity = static_cast<double>(num_edges) * BALANCE_RATIO / num_partitions + 1;
    occupied.assign(num_partitions, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    // graph.resize(num_vertices);
    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    edgelist2bucket.assign(num_edges, kInvalidBid);

    Q.resize(num_partitions);
    V = dense_bitset(num_vertices);
    super.resize(num_vertices, kInvalidVid);
    free_vertex.resize(num_partitions, 0);
    std::iota(free_vertex.begin(), free_vertex.end(), 0);


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
    // graph.stream_build(fin, num_edges, degrees);
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();
};

template <typename TAdj>
void HybridBLPartitioner<TAdj>::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << num_partitions;

    Timer compute_timer;

    LOG(INFO) << "partitioning...";
    compute_timer.start();
    
    // V is not empty
    while (true) {
        bool stop = true;
        for (bid_t m = 0; m < num_partitions; ++m) {
            if (Q[m].empty()) {
                bool has_free_vertex = get_free_vertex(m);
                stop &= !has_free_vertex;
                if (!has_free_vertex) {
                    continue;
                }

                auto vid = free_vertex[m];
                // if (degrees[vid] < degree_threshold) {
                if (adj_in[vid].size() < degree_threshold && super[vid] == kInvalidVid) {
                    // assert(super[vid] == kInvalidVid);
                    init_fusion(m, vid, 0);

                // } else if (degrees[vid] >= degree_threshold) {
                } else if (adj_in[vid].size() >= degree_threshold) {
                    fission(m, vid);
                }
            } else {
                stop = false;
                auto [uid, root, dist] = Q[m].front();
                Q[m].pop();
                if (dist < gamma && super[uid] == kInvalidVid) {
                // if (root_assigned[root] < fusion_threshold && super[uid] == kInvalidVid) {
                    fusion(m, uid, root, dist);
                }
            }
        }
        if (stop) {
            break;
        }
    }

    compute_timer.stop();

    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
    // eid_t sum = 0;
    // for (bid_t b = 0; b < num_partitions; ++b) {
    //     LOG(INFO) << b << ": " << fusion_occupy[b] << '\t' << fission_occupy[b];
    //     sum += fusion_occupy[b] + fission_occupy[b];
    // }
    // eid_t max_assigned = 0;
    // for (const auto& [vid, assigned] : root_assigned) {
    //     // LOG(INFO) << vid << ' ' << assigned;
    //     max_assigned = std::max(assigned, max_assigned);
    // }
    // LOG(INFO) << "max_assigned: " << max_assigned;
}

template <typename TAdj>
void HybridBLPartitioner<TAdj>::init_fusion(bid_t machine, vid_t vid, vid_t dist)
{
    super[vid] = vid;
    bid_t candidate = kInvalidBid;
    eid_t candidate_size = std::numeric_limits<eid_t>::max();
    for (bid_t b = 0; b < num_partitions; ++b) {
        if (occupied[b] < candidate_size) {
            candidate = b;
            candidate_size = occupied[b];
        }
    }
    root_bucket[vid] = candidate;
    fusion(machine, vid, vid, dist);
}

template <typename TAdj>
void HybridBLPartitioner<TAdj>::fusion(bid_t machine, vid_t vid, vid_t root, vid_t dist)
{
    // LOG(INFO) << "Fussion function begin(vid, dist, ind): " 
    //     << vid << ' ' << dist << ' ' 
    //     << adj_in[vid].size() << ' ' << adj_out[vid].size();
    V.set(vid, 1);
    vid_t current_super = super[root];
    assert(super[root] == root);
    super[vid] = current_super;
    bid_t current_super_bid = root_bucket[root];
    for (int direction = 0; direction < 2; ++direction) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        for (vid_t i = 0; i < neighbors.size(); ++i) {
            eid_t eid = neighbors[i].v;
            assert(eid < num_edges && eid >= 0);
            if (edgelist2bucket[eid] != kInvalidBid) continue;

            vid_t &uid = direction ? edges[eid].second : edges[eid].first;

            assign_edge(current_super_bid, direction ? vid : uid, direction ? uid : vid, eid);
            // ++fusion_occupy[current_super_bid];
            // ++root_assigned[root];
            Q[machine].push({uid, root, dist + 1});
        }
    }

    // AdjList<AdjEntryVid>& neighbors = graph[vid];
    // vid_t overlap = 0, neighbors_cnt = neighbors.size();
    // for (vid_t i = 0; i < neighbors.size(); ++i) {
    //     vid_t uid = neighbors[i].vid;
    // }
}

template <typename TAdj>
void HybridBLPartitioner<TAdj>::fission(bid_t machine, vid_t vid)
{
    V.set(vid, 1);
    for (int direction = 1; direction < 2; ++direction) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        for (vid_t i = 0; i < neighbors.size(); ++i) {
            eid_t eid = neighbors[i].v;
            assert(eid < num_edges && eid >= 0);
            if (edgelist2bucket[eid] != kInvalidBid) continue;

            vid_t &uid = direction ? edges[eid].second : edges[eid].first;
            assign_edge(uid % num_partitions, direction ? vid : uid, direction ? uid : vid, eid);
            // ++fission_occupy[uid % num_partitions];
        }
    }
}
