#include "smp_partitioner.hpp"
#include "ne_partitioner.hpp"
#include "hdrf_partitioner.hpp"
#include "hep_partitioner.hpp"
#include "dbh_partitioner.hpp"
#include "ebv_partitioner.hpp"
#include "conversions.hpp"
DECLARE_int32(k);

SmpPartitioner::SmpPartitioner(std::string basefilename)
    : basefilename(basefilename), writer(basefilename, FLAGS_write)
{
    p = FLAGS_p;
    k = FLAGS_k;
    LOG(INFO) << "k = " << k
                << ", p = " << p;

    total_time.start();
    std::string split_method = FLAGS_method == "smp" ? "ne" : FLAGS_method.substr(4);
    if (split_method == "ne") {
        partitioner = new NePartitioner(FLAGS_filename, true);
    } else if (split_method == "dbh") {
        partitioner = new DbhPartitioner(FLAGS_filename, true);
    } else if (split_method == "ebv") {
        partitioner = new EbvPartitioner(FLAGS_filename, true);
    } else if (split_method == "hdrf") {
        partitioner = new HdrfPartitioner(FLAGS_filename, true);
    // } else if (split_method == "hep") {
    //     partitioner = new HepPartitioner(FLAGS_filename, true);
    } else {
        LOG(ERROR) << "Unknown split method!";
    }

    num_vertices = partitioner->num_vertices;
    num_edges = partitioner->num_edges;

    bucket_info.assign(k * p, bucket_info_item(num_vertices));
    for (int i = 0; i < k * p; ++ i) bucket_info[i].old_id = i;

    edge2bucket.assign(num_edges, -1);
    occupied.assign(p, 0);
    bucket_edge_cnt.assign(FLAGS_p, 0);
};

void SmpPartitioner::merge()
{
    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    for (int b = 0; b < p * k; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_mirror.popcount();
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
    }

    // std::vector<vid_t> boundary_replicate_times(FLAGS_k * FLAGS_p + 1);
    // vid_t boundary_vertice_cnt = 0;
    // for (vid_t vid = 0; vid < num_vertices; ++ vid) {
    //     int cnt = 0;
    //     boundary_vertice_cnt ++;
    //     for (int b = 0; b < FLAGS_k * FLAGS_p; ++ b) {
    //         cnt += bucket_info[b].is_mirror.get(vid);
    //     }
    //     boundary_replicate_times[cnt] ++;
    // }
    // for (int i = 0; i <= FLAGS_k * FLAGS_p; ++ i) {
    //     LOG(INFO) << i << ' ' << boundary_replicate_times[i] << ' ' << (double)boundary_replicate_times[i] / boundary_vertice_cnt;
    // }

    std::sort(bucket_info.begin(), bucket_info.end(), 
        [&](const bucket_info_item &l, const bucket_info_item &r) {
        return l.replicas > r.replicas;
    });


    curr_bucket_id = 0;
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >

    // valid_bucket = merge_by_size();
    valid_bucket = merge_by_overlap();
    for (int b = 0; b < p * k; ++ b) 
        DLOG(INFO)   << "Bucket_info " << bucket_info[b].old_id 
                    << " vertices: " << bucket_info[b].replicas 
                    << " edges: " << bucket_info[b].occupied
                    << " is_chosen: " << bucket_info[b].is_chosen << std::endl;

    // rearrange bucket after heuristic merging
    std::sort(bucket_info.begin(), bucket_info.end());

    for (size_t b = 0; b < bucket_info.size(); ) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        if (is_chosen) {
            old_id = get_final_bucket(old_id);
            // std::cerr << "Is chosen, new old_id: " << old_id << ", replicas: " << replicas << '\n';
            ++ b;
        } else {
            std::swap(bucket_info[b], bucket_info.back());
            bucket_info.pop_back();
        }
    }
    bucket_info.shrink_to_fit();

    size_t curr_assigned_edges = rearrange_edge(edges, valid_bucket);
    assigned_edges += curr_assigned_edges;
}

std::unordered_map<int, int> SmpPartitioner::merge_by_size()
{
    std::unordered_map<int, int> valid_bucket;  // < old bucket, new bucket >
    
    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::priority_queue<std::tuple<int, int, int, int>, std::vector<std::tuple<int, int, int, int>>, std::greater<std::tuple<int, int, int, int>>> pq;  
    for (int b = 0; b < FLAGS_p; ++ b) {
        pq.emplace(0, 0, b, b);
    }
    for (int b = 0; b < FLAGS_p * FLAGS_k; ++ b) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = pq.top();
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

    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    std::vector<std::tuple<int, int, int, int>> final_bucket;
    for (int b = 0; b < p; ++ b) {
        final_bucket.emplace_back(0, 0, b, b);
    }
    auto compute_new_bucket_size = [&](int bid_a, int bid_b) -> size_t {
        const auto &is_mirror_a = bucket_info[bid_a].is_mirror, &is_mirror_b = bucket_info[bid_b].is_mirror;
        dense_bitset new_bucket = is_mirror_a | is_mirror_b;
        return new_bucket.popcount();
    };

    for (int b = 0; b < p * k; ++ b) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];

        int best_final_bucket = -1;
        size_t min_size_after_merge = std::numeric_limits<size_t>::max();

        for (size_t fb = 0; fb < final_bucket.size(); ++ fb) {
            auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[fb];
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

        auto [mirror_cnt, partitions_inside, parent_bucket, parent_old_id] = final_bucket[best_final_bucket];
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

void SmpPartitioner::calculate_stats()
{
    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<size_t> bucket2vcnt(p, 0);
    
    sort(bucket_info.begin(), bucket_info.end(), [&](const auto &l, const auto &r) {
        return l.old_id < r.old_id;
    });

    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    size_t max_part_edge_cnt = 0, all_part_edge_cnt = 0; 
    for (int b = 0; b < p; ++ b) {
        bucket_info[b].replicas = bucket_info[b].is_mirror.popcount();
        bucket2vcnt[b] = bucket_info[b].replicas;
        occupied[b] = bucket_info[b].occupied;
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
        max_part_edge_cnt = std::max(max_part_edge_cnt, bucket_info[b].occupied);
        all_part_edge_cnt += bucket_info[b].occupied;
        CHECK_EQ(bucket_info[b].occupied, bucket_edge_cnt[b]);
    }

    for (int b = 0; b < p; ++ b) 
        LOG(INFO) << "Bucket_info: " << bucket_info[b].old_id 
                << ", vertices: " << bucket_info[b].replicas 
                << ", edges: " << bucket_info[b].occupied;
    
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


int SmpPartitioner::merge_bucket(int dst, int src, bool &has_intersection)   // dst, src
{
    auto &is_mirror_a = bucket_info[dst].is_mirror, &is_mirror_b = bucket_info[src].is_mirror;
    size_t mirror_cnt = 0;
    for (size_t i = 0; i < is_mirror_a.size(); ++ i) {
        if (is_mirror_a.get(i) == 1 || is_mirror_b.get(i) == 1) {
            is_mirror_a.set_bit_unsync(i);
            ++ mirror_cnt;
        }
        if (!has_intersection && is_mirror_a.get(i) == 1 && is_mirror_b.get(i) == 1) {
            has_intersection = true;
        }
    }
    bucket_info[dst].replicas = mirror_cnt;
    bucket_info[dst].occupied += bucket_info[src].occupied;
    return mirror_cnt;
}

void SmpPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    LOG(INFO) << "partitioning...";
    compute_timer.start();

    partitioner->split();

    std::swap(edges, partitioner->edges);
    std::swap(edge2bucket, partitioner->edge2bucket);
    for (int bucket = 0; bucket < p * k; bucket++) {
        std::swap(partitioner->is_boundarys[bucket], bucket_info[bucket].is_mirror);
        std::swap(partitioner->occupied[bucket], bucket_info[bucket].occupied);
    }

    std::cerr << "\n" << std::string(25, '#') << " Split phase end, Merge phase start " << std::string(25, '#') << "\n\n";
    merge();
    CHECK_EQ(assigned_edges, num_edges);


    compute_timer.stop();
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();

    calculate_stats();
    

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    bool check_res = check_edge();
    CHECK_EQ(check_res, true);
    
    if (FLAGS_write) 
        LOG(INFO) << "Writing result...";

    for (size_t i = 0; i < edge2bucket.size(); ++ i) 
        writer.save_edge(edges[i].first, edges[i].second, edge2bucket[i]);
}