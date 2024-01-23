#include <queue>

#include "fsm_partitioner.hpp"
#include "ne_partitioner.hpp"
#include "hep_partitioner.hpp"
// #include "hdrf_partitioner.hpp"
// #include "dbh_partitioner.hpp"
// #include "ebv_partitioner.hpp"
#include "conversions.hpp"
DECLARE_int32(k);
DECLARE_bool(fastmerge);

FsmPartitioner::FsmPartitioner(std::string basefilename)
    : basefilename(basefilename), writer(basefilename, FLAGS_write)
{
    num_partitions = FLAGS_p;
    k = FLAGS_k;
    split_partitioner = nullptr;
    LOG(INFO) << "k = " << k
                << ", num_partitions = " << num_partitions;

    split_method = FLAGS_method == "fsm" ? "ne" : FLAGS_method.substr(4);
    if (split_method == "ne") {
        split_partitioner = std::make_unique<NePartitioner<adj_with_bid_t>>(FLAGS_filename, true);
    } else if (split_method == "hep") {
        split_partitioner = std::make_unique<HepPartitioner<adj_with_bid_t>>(FLAGS_filename, true);
    }
    // else if (split_method == "dbh") {
    //     split_partitioner = std::make_unique<DbhPartitioner>(FLAGS_filename, true);
    // } else if (split_method == "ebv") {
    //     split_partitioner = std::make_unique<EbvPartitioner>(FLAGS_filename, true);
    // } else if (split_method == "hdrf") {
    //     split_partitioner = std::make_unique<HdrfPartitioner>(FLAGS_filename, true);
    // } else {
    //     LOG(ERROR) << "Unknown split method!";
    // }

    num_vertices = split_partitioner->num_vertices;
    num_edges = split_partitioner->num_edges;

    bucket_info.assign(k * num_partitions, BucketInfo(num_vertices));
    for (bid_t i = 0; i < k * num_partitions; ++i) bucket_info[i].old_id = i;

    edgelist2bucket.assign(num_edges, kInvalidBid);
    occupied.assign(num_partitions, 0);
};

void FsmPartitioner::merge()
{
    size_t max_part_vertice_cnt = 0, all_part_vertice_cnt = 0; 
    for (bid_t b = 0; b < num_partitions * k; ++b) {
        bucket_info[b].replicas = bucket_info[b].is_mirror.popcount();
        max_part_vertice_cnt = std::max(max_part_vertice_cnt, bucket_info[b].replicas);
        all_part_vertice_cnt += bucket_info[b].replicas;
    }

    // std::vector<vid_t> boundary_replicate_times(FLAGS_k * FLAGS_p + 1);
    // vid_t boundary_vertice_cnt = 0;
    // for (vid_t vid = 0; vid < num_vertices; ++vid) {
    //     int cnt = 0;
    //     boundary_vertice_cnt ++;
    //     for (bid_t b = 0; b < FLAGS_k * FLAGS_p; ++b) {
    //         cnt += bucket_info[b].is_mirror.get(vid);
    //     }
    //     boundary_replicate_times[cnt] ++;
    // }
    // for (int i = 0; i <= FLAGS_k * FLAGS_p; ++i) {
    //     LOG(INFO) << i << ' ' << boundary_replicate_times[i] << ' ' << (double)boundary_replicate_times[i] / boundary_vertice_cnt;
    // }

    std::sort(bucket_info.begin(), bucket_info.end(), 
        [&](const BucketInfo &l, const BucketInfo &r) {
        return l.replicas > r.replicas;
    });


    curr_bucket_id = 0;
    std::unordered_map<bid_t, bid_t> valid_bucket;  // < old bucket, new bucket >

    if (FLAGS_fastmerge) {
        valid_bucket = fast_merge();
    } else {
        valid_bucket = precise_merge();
    }
    for (bid_t b = 0; b < num_partitions * k; ++b) 
        DLOG(INFO)   << "Bucket_info " << bucket_info[b].old_id 
                    << " vertices: " << bucket_info[b].replicas 
                    << " edges: " << bucket_info[b].occupied
                    << " is_chosen: " << bucket_info[b].is_chosen << std::endl;

    // rearrange bucket after heuristic merging
    std::sort(bucket_info.begin(), bucket_info.end());

    for (eid_t b = 0; b < bucket_info.size(); ) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];
        if (is_chosen) {
            old_id = get_final_bucket(old_id);
            // std::cerr << "Is chosen, new old_id: " << old_id << ", replicas: " << replicas << '\n';
            ++b;
        } else {
            std::swap(bucket_info[b], bucket_info.back());
            bucket_info.pop_back();
        }
    }
    bucket_info.shrink_to_fit();

    eid_t curr_assigned_edges;
    if (split_method == "hep") {
        curr_assigned_edges = rearrange_edge_hybrid(valid_bucket);
    } else {
        curr_assigned_edges = rearrange_edge(edges, valid_bucket);
    }
    
    assigned_edges += curr_assigned_edges;
}

std::unordered_map<bid_t, bid_t> FsmPartitioner::fast_merge()
{
    std::unordered_map<bid_t, bid_t> valid_bucket;  // < old bucket, new bucket >
    
    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    using bucket_item = std::tuple<vid_t, bid_t, bid_t, bid_t>;
    std::priority_queue<bucket_item, std::vector<bucket_item>, std::greater<bucket_item>> pq;  
    for (bid_t b = 0; b < FLAGS_p; ++b) {
        pq.emplace(0, 0, b, b);
    }
    for (bid_t b = 0; b < FLAGS_p * FLAGS_k; ++b) {
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

std::unordered_map<bid_t, bid_t> FsmPartitioner::precise_merge()
{
    std::unordered_map<bid_t, bid_t> valid_bucket;  // < old bucket, new bucket >

    // < mirror_cnt, partitions_inside, index_in_bucket_info, old_id >
    using bucket_item = std::tuple<vid_t, bid_t, bid_t, bid_t>;
    std::vector<bucket_item> final_bucket;
    for (bid_t b = 0; b < num_partitions; ++b) {
        final_bucket.emplace_back(0, 0, b, b);
    }
    auto compute_new_bucket_size = [&](bid_t bid_a, bid_t bid_b) {
        const auto &is_mirror_a = bucket_info[bid_a].is_mirror, &is_mirror_b = bucket_info[bid_b].is_mirror;
        dense_bitset new_bucket = is_mirror_a | is_mirror_b;
        return new_bucket.popcount();
    };

    for (bid_t b = 0; b < num_partitions * k; ++b) {
        auto &[is_mirror, occupied, old_id, replicas, is_chosen] = bucket_info[b];

        bid_t best_final_bucket = kInvalidBid;
        size_t min_size_after_merge = std::numeric_limits<size_t>::max();

        for (bid_t fb = 0; fb < final_bucket.size(); ++fb) {
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
        CHECK_NE(best_final_bucket, kInvalidBid);

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

void FsmPartitioner::calculate_stats()
{
    sort(bucket_info.begin(), bucket_info.end(), [&](const auto &l, const auto &r) {
        return l.old_id < r.old_id;
    });

    is_boundarys.resize(num_partitions);
    occupied.resize(num_partitions);
    for (bid_t b = 0; b < num_partitions; ++b) {
        std::swap(bucket_info[b].is_mirror, this->is_boundarys[b]);
        std::swap(bucket_info[b].occupied, this->occupied[b]);
    }
    EdgePartitioner::calculate_stats(true);
    for (bid_t b = 0; b < num_partitions; ++b) {
        std::swap(this->is_boundarys[b], bucket_info[b].is_mirror);
        std::swap(this->occupied[b], bucket_info[b].occupied);
    }
}


vid_t FsmPartitioner::merge_bucket(bid_t dst, bid_t src, bool &has_intersection)   // dst, src
{
    auto &is_mirror_a = bucket_info[dst].is_mirror, &is_mirror_b = bucket_info[src].is_mirror;
    size_t mirror_cnt = 0;
    for (size_t i = 0; i < is_mirror_a.size(); ++i) {
        if (is_mirror_a.get(i) == 1 || is_mirror_b.get(i) == 1) {
            is_mirror_a.set_bit_unsync(i);
            ++mirror_cnt;
        }
        if (!has_intersection && is_mirror_a.get(i) == 1 && is_mirror_b.get(i) == 1) {
            has_intersection = true;
        }
    }
    bucket_info[dst].replicas = mirror_cnt;
    bucket_info[dst].occupied += bucket_info[src].occupied;
    return mirror_cnt;
}

void FsmPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << num_partitions;

    Timer merge_timer, split_timer;

    LOG(INFO) << "partitioning...";

    split_partitioner->split();
    {
        for (bid_t bucket = 0; bucket < num_partitions * k; ++bucket) {
            std::swap(split_partitioner->is_boundarys[bucket], bucket_info[bucket].is_mirror);
            std::swap(split_partitioner->occupied[bucket], bucket_info[bucket].occupied);
        }

        std::swap(split_partitioner->partition_time, partition_time);
        split_timer = partition_time;
        partition_time.start();
        
        if (split_method == "hep") {
            std::swap(split_partitioner->edges, edges);
            std::swap(split_partitioner->degrees, degrees);
            std::swap(split_partitioner->mem_graph, mem_graph);
            std::swap(split_partitioner->edgelist2bucket, edgelist2bucket);
            // std::cerr << "edgelist2bucket.size: " << edgelist2bucket.size() << std::endl;
            
            // if (edges.size() == 0) {
            //     LOG(INFO) << "Loading edges list...";
            //     std::ifstream fin(binedgelist_name(basefilename),
            //             std::ios::binary | std::ios::ate); 
            //     fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
            //     edges.resize(num_edges);
            //     fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);
            // }
        } else {
            std::swap(split_partitioner->edges, edges);
            std::swap(split_partitioner->edgelist2bucket, edgelist2bucket);
            
            if (edges.size() == 0) {
                LOG(INFO) << "Loading edges list...";
                std::ifstream fin(binedgelist_name(basefilename),
                        std::ios::binary | std::ios::ate); 
                fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
                edges.resize(num_edges);
                fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);
            }
        }
        
    }

    std::cerr << "\n" << std::string(25, '#') << " Split phase end, Merge phase start " << std::string(25, '#') << "\n\n";

    merge_timer.start();
    {
        merge();
    }
    merge_timer.stop();
    partition_time.stop();

    CHECK_EQ(assigned_edges, num_edges);

    LOG(INFO) << "spliting time: " << split_timer.get_time();
    LOG(INFO) << "merging time: " << merge_timer.get_time();
    LOG(INFO) << "partitioning time: " << partition_time.get_time();

    calculate_stats();

    if (split_method == "hep") {
        CHECK_EQ(check_edge_hybrid(), true);
    } else {
        CHECK_EQ(check_edge(), true);
    }



    if (FLAGS_write) 
        LOG(INFO) << "Writing result...";

    if (split_method == "hep") {
        save_edge_hybrid();
    } else {
        for (eid_t i = 0; i < edgelist2bucket.size(); ++i) {
            writer.save_edge(edges[i].first, edges[i].second, edgelist2bucket[i]);
        }
    }

    
}