#include "smp_partitioner.hpp"
#include "conversions.hpp"
DECLARE_int32(k);

SmpPartitioner::SmpPartitioner(std::string basefilename)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename)
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
    theta = 0.10;
    LOG(INFO) << "k = " << FLAGS_k;
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0;
    capacity = (double)num_edges * BALANCE_RATIO / p + 1;
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);

    bucket_info.assign(p, bucket_info_item(num_vertices));
    // is_in_a_core = dense_bitset(num_vertices);
    for (int i = 0; i < p; ++ i) bucket_info[i].old_id = i;

    edge2bucket.assign(num_edges, -1);
    bucket_edge_cnt.assign(FLAGS_p, 0);
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";
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

void SmpPartitioner::assign_remaining(int last_p)
{
    auto &is_boundary = bucket_info[last_p].is_boundary, &is_core = bucket_info[last_p].is_core;
    // auto &is_boundary = bucket_info[last_p].is_boundary, &is_core = is_in_a_core;
    repv (u, num_vertices)
        for (auto &i : adj_out[u])
            if (edges[i.v].valid()) {
                assign_edge(last_p, u, edges[i.v].second, i.v);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }

    repv (i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            for (int j = 0; j < last_p; ++ j) {
                if (bucket_info[j].is_core.get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
            }
        }
    }
}

void SmpPartitioner::count_mirrors_k_split()
{
    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    for (int b = filled_bucket_cnt; b < FLAGS_p + (FLAGS_k - 1) * p; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_boundary.popcount();
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
    }

    std::vector<vid_t> boundary_replicate_times(FLAGS_k * FLAGS_p + 1);
    vid_t boundary_vertice_cnt = 0;
    for (vid_t vid = 0; vid < num_vertices; ++ vid) {
        int cnt = 0;
        boundary_vertice_cnt ++;
        for (int b = 0; b < FLAGS_k * FLAGS_p; ++ b) {
            cnt += bucket_info[b].is_boundary.get(vid);
        }
        boundary_replicate_times[cnt] ++;
    }
    for (int i = 0; i <= FLAGS_k * FLAGS_p; ++ i) {
        LOG(INFO) << i << ' ' << boundary_replicate_times[i] << ' ' << (double)boundary_replicate_times[i] / boundary_vertice_cnt;
    }

    std::sort(bucket_info.begin() + filled_bucket_cnt, bucket_info.end(), [&](const bucket_info_item &l, const bucket_info_item &r) {
        return l.replicas > r.replicas;
    });

    for (int b = filled_bucket_cnt; b < FLAGS_p + (FLAGS_k - 1) * p; ++ b) 
        LOG(INFO) << "Bucket_info " << bucket_info[b].old_id << " vertices: " << bucket_info[b].replicas 
            << " edges: " << bucket_info[b].occupied << " is_chosen: " << bucket_info[b].is_chosen << std::endl;

    curr_bucket_id = filled_bucket_cnt;
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >

    // valid_bucket = merge_by_size();
    valid_bucket = merge_by_overlap();

    // rearrange bucket after heuristic merging
    std::sort(bucket_info.begin() + filled_bucket_cnt, bucket_info.end());

    for (size_t b = filled_bucket_cnt; b < bucket_info.size(); ) {
        auto &[is_core, is_boundary, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        if (is_chosen) {
            old_id = get_final_bucket(old_id);
            ++ b;
        } else {
            std::swap(bucket_info[b], bucket_info.back());
            bucket_info.pop_back();
        }
    }
    bucket_info.shrink_to_fit();
    for (size_t b = filled_bucket_cnt; b < bucket_info.size(); ++ b) 
        LOG(INFO) << "Bucket_info " << bucket_info[b].old_id << " vertices: " << bucket_info[b].replicas 
        << " edges: " << bucket_info[b].occupied << " rank: " << b;

    size_t curr_assigned_edges = rearrange_edge(edges, valid_bucket, true);
    assigned_edges += curr_assigned_edges;
}

std::unordered_map<int, int> SmpPartitioner::merge_by_size()
{
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >
    
    // < vertice_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::priority_queue<std::tuple<int, int, int, int>, std::vector<std::tuple<int, int, int, int>>, std::greater<std::tuple<int, int, int, int>>> pq;  
    for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) {
        pq.emplace(0, 0, b, b);
    }
    for (int b = filled_bucket_cnt; b < FLAGS_p + (FLAGS_k - 1) * p; ++ b) {
        auto &[is_core, is_boundary, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        auto [vertice_cnt, partitions_inside, parent_bucket, parent_old_id] = pq.top();
        pq.pop();
        if (partitions_inside == 0) {
            parent_old_id = get_final_bucket(old_id);
            pq.emplace(replicas, partitions_inside + 1, b, parent_old_id);
            is_chosen = true;
        } else {
            bool has_intersection = false;
            size_t new_vertice_cnt = merge_bucket(parent_bucket, b, has_intersection);
            if (!has_intersection) { 
                std::cerr << "No intersection betwenn bucket " << b << " and " << parent_bucket << '\n';
                // exit(1);
            }
            if (partitions_inside + 1 < FLAGS_k) {
                pq.emplace(new_vertice_cnt, partitions_inside + 1, parent_bucket, parent_old_id);
            }
        }
        valid_bucket.emplace(old_id, parent_old_id);
    }
    return valid_bucket;
}

std::unordered_map<int, int> SmpPartitioner::merge_by_overlap()
{
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >

    // < vertice_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::vector<std::tuple<int, int, int, int>> final_bucket;
    for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) {
        final_bucket.emplace_back(0, 0, b, b);
    }
    auto compute_new_bucket_size = [&](int bid_a, int bid_b) -> size_t {
        const auto &is_boundary_a = bucket_info[bid_a].is_boundary, &is_boundary_b = bucket_info[bid_b].is_boundary;
        dense_bitset new_bucket = is_boundary_a | is_boundary_b;
        return new_bucket.popcount();
    };

    for (int b = filled_bucket_cnt; b < FLAGS_p + (FLAGS_k - 1) * p; ++ b) {
        auto &[is_core, is_boundary, occupied, old_id, replicas, is_chosen] = bucket_info[b];

        int best_final_bucket = -1;
        size_t min_size_after_merge = std::numeric_limits<size_t>::max();

        for (int fb = 0; fb < final_bucket.size(); ++ fb) {
            auto [vertice_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[fb];
            if (partitions_inside == FLAGS_k) continue;
            if (partitions_inside == 0) {
                best_final_bucket = fb;
                break;
            }
            size_t new_size_after_merge = compute_new_bucket_size(parent_bucket, b);
            if (min_size_after_merge > new_size_after_merge) {
                min_size_after_merge = new_size_after_merge;
                best_final_bucket = fb;
            }
        }
        CHECK_NE(best_final_bucket, -1);

        auto [vertice_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[best_final_bucket];
        if (partitions_inside == 0) {
            parent_old_id = get_final_bucket(old_id);
            final_bucket[best_final_bucket] = {replicas, partitions_inside + 1, b, parent_old_id};
            is_chosen = true;
        } else {
            bool has_intersection = false;
            size_t new_vertice_cnt = merge_bucket(parent_bucket, b, has_intersection);
            if (!has_intersection) { 
                std::cerr << "No intersection betwenn bucket " << b << " and " << parent_bucket << '\n';
                // exit(1);
            }

            final_bucket[best_final_bucket] = {new_vertice_cnt, partitions_inside + 1, parent_bucket, parent_old_id};
        }
        valid_bucket.emplace(old_id, parent_old_id);
    }
    return valid_bucket;
}

std::array<size_t, 2> SmpPartitioner::count_mirrors_elect(int iter)
{
    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_boundary.popcount();
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
    }

    double avg_vertice_cnt = (double)all_part_vertice_cnt / p;
    LOG(INFO) << "master balance: "
              << (double)max_part_vertice_cnt / avg_vertice_cnt;

    LOG(INFO) << "replication factor: " << (double)all_part_vertice_cnt / num_vertices;
    
    if (iter == 0)  { valid_threshold = {size_t(avg_vertice_cnt * 1.05 + 1), size_t(avg_vertice_cnt * (1.0 + theta) + 1)}; }
    std::unordered_map<int, int> valid_bucket;
    uint8_t valid_bucket_cnt = 0;
    for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) {
        auto &[is_core, is_boundary, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        if (valid_threshold[0] <= replicas && replicas <= valid_threshold[1]) {
            is_chosen = true;
            valid_bucket.emplace(old_id, filled_bucket_cnt + valid_bucket_cnt);
            old_id = valid_bucket[old_id];
            ++ valid_bucket_cnt;
        }
    }

    sort(bucket_info.begin() + filled_bucket_cnt, bucket_info.end());

    for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) 
        LOG(INFO) << "Bucket_info " << bucket_info[b].old_id << " vertices: " << bucket_info[b].replicas 
        << " edges: " << bucket_info[b].occupied << " rank: " << b << " is_chosen: " << bucket_info[b].is_chosen << std::endl;
        
    for (auto [bucket, new_bucket] : valid_bucket) {
        LOG(INFO) << "Valid partition(old/new) " << bucket << " / " << new_bucket << " : " << bucket_info[new_bucket].is_boundary.popcount() << " vertices " << std::endl;
        assert(valid_threshold[0] <= bucket_info[new_bucket].is_boundary.popcount() && bucket_info[new_bucket].is_boundary.popcount() <= valid_threshold[1]);
    }
    size_t curr_assigned_edges = rearrange_edge(edges, valid_bucket, false);

    assigned_edges += curr_assigned_edges;
    // edge_shrink(edges);
    // edges = edge_backup;    // shrink_to_fit ? 
    // edges.shrink_to_fit();
    // CHECK_EQ(edges.size() + assigned_edges, num_edges);
    adj_out.build(edges);
    adj_in.build_reverse(edges);
    p -= valid_bucket_cnt;
    filled_bucket_cnt += valid_bucket_cnt;
    // if (iter < FLAGS_round - 1)
        for (int b = filled_bucket_cnt; b < FLAGS_p; ++ b) { bucket_info[b].clear(); bucket_info[b].old_id = b; }
    LOG(INFO) << valid_threshold[0] << ' ' << valid_threshold[1] << ' ' << assigned_edges << ' ' << FLAGS_p - filled_bucket_cnt;
    return {all_part_vertice_cnt, valid_bucket_cnt};
}

void SmpPartitioner::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    
    sort(bucket_info.begin(), bucket_info.end(), [&](const auto &l, const auto &r) {
        return l.old_id < r.old_id;
    });

    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    size_t max_part_edge_cnt = 0, all_part_edge_cnt = 0; 
    for (int b = 0; b < FLAGS_p; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_boundary.popcount();
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
        max_part_edge_cnt = std::max(max_part_edge_cnt, bucket_edge_cnt[b]);
        all_part_edge_cnt += bucket_edge_cnt[b];
    }

    for (int b = 0; b < FLAGS_p; ++ b) 
        LOG(INFO) << "Bucket_info " << bucket_info[b].old_id << " vertices: " << bucket_info[b].replicas 
                << " edges: " << bucket_info[b].occupied << " rank: " << b << std::endl;
    
    double avg_vertice_cnt = (double)all_part_vertice_cnt / FLAGS_p;
    double avg_edge_cnt = (double)all_part_edge_cnt / FLAGS_p;
    double std_deviation = 0.0;
    for (int b = 0; b < FLAGS_p; ++ b) 
        std_deviation += pow(bucket_info[b].replicas - avg_vertice_cnt, 2);
    std_deviation = sqrt((double)std_deviation / FLAGS_p);
    
    LOG(INFO) << "Vertice balance: "
              << (double)max_part_vertice_cnt / ((double)num_vertices / FLAGS_p);
    LOG(INFO) << "Max Vertice count: "
              << max_part_vertice_cnt;
    LOG(INFO) << "Vertice std_deviation / avg: "
              << std_deviation / avg_vertice_cnt;
    LOG(INFO) << "Edge balance: "
              << (double)max_part_edge_cnt / avg_edge_cnt;
    CHECK_EQ(all_part_edge_cnt, num_edges);

    LOG(INFO) << "replication factor (final): " << (double)all_part_vertice_cnt / num_vertices;
}


int SmpPartitioner::merge_bucket(int dst, int src, bool &has_intersection)   // dst, src
{
    auto &is_boundary_a = bucket_info[dst].is_boundary, &is_boundary_b = bucket_info[src].is_boundary;
    size_t vertice_cnt = 0;
    repv (i, num_vertices) {
        if (is_boundary_a.get(i) == 1 || is_boundary_b.get(i) == 1) {
            is_boundary_a.set_bit_unsync(i);
            ++ vertice_cnt;
        }
        if (!has_intersection && is_boundary_a.get(i) == 1 && is_boundary_b.get(i) == 1) {
            has_intersection = true;
        }
    }
    bucket_info[dst].replicas = vertice_cnt;
    bucket_info[dst].occupied += bucket_info[src].occupied;
    return vertice_cnt;
}

void SmpPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    compute_timer.start();

    k_split();

    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();
    calculate_stats();
    
    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    bool check_res = check_edge();
    CHECK_EQ(check_res, true);
    if (FLAGS_write) {
        LOG(INFO) << "Writing result...";
    }
    for (size_t i = 0; i < edge2bucket.size(); ++ i) {
        const auto [from, to] = edges[i];
        writer.save_edge(from, to, edge2bucket[i]);
    }
}

size_t SmpPartitioner::elect(int iter)
{
    for (bucket = filled_bucket_cnt; bucket < FLAGS_p - 1; bucket++) {
        expand();
    }
    bucket = FLAGS_p - 1;
    std::cerr << bucket << std::endl;
    assign_remaining(bucket);
    LOG(INFO) << "expected edges in each partition: " << num_edges / p;
    for (int b = filled_bucket_cnt; b < filled_bucket_cnt + p; ++ b)
        DLOG(INFO) << "edges in partition " << b << ": " << bucket_info[b].occupied;
    size_t max_occupied = 0;
    rep(b, p) max_occupied = std::max(max_occupied, bucket_info[b].occupied);

    LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
    auto [total_mirrors, filled_bucket_this_iter] = count_mirrors_elect(iter);
    LOG(INFO) << "total mirrors: " << total_mirrors;
    LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;
    std::cerr << std::string(75, '#') << "\n";

    return filled_bucket_this_iter;
}

void SmpPartitioner::k_split()
{
    capacity /= FLAGS_k;
    for (int b = FLAGS_p; b < FLAGS_p + (FLAGS_k - 1) * p; ++ b) {
        bucket_info.emplace_back(bucket_info_item(num_vertices));
        bucket_info[b].old_id = b;
    }
    for (bucket = filled_bucket_cnt; bucket < FLAGS_p + (FLAGS_k - 1) * p - 1; bucket++) {
        expand();
    }
    bucket = FLAGS_p + (FLAGS_k - 1) * p - 1;
    std::cerr << bucket << std::endl;
    assign_remaining(bucket);

    count_mirrors_k_split();
    std::cerr << std::string(75, '#') << "\n";
}