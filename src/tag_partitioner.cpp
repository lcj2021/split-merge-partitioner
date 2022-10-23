#include "tag_partitioner.hpp"
#include "conversions.hpp"
#include "KM.hpp"
using namespace std;
#define all(x) x.begin(), x.end()

TagPartitioner::TagPartitioner(std::string basefilename)
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
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0;
    capacity = (double)num_vertices * BALANCE_RATIO / p + 1;
    occupied.assign(p + 1, 0);
    adj_out.resize(num_vertices + 1);
    adj_in.resize(num_vertices + 1);
    global_tag_distribute.assign(p + 1, 0);
    vertex2tag.assign(num_vertices, dense_bitset(p + 1));
    can_cover.assign(num_vertices, false);
    tag_valid.assign(p + 1, false);
    bucket.resize(p + 1);
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";
    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    num_edges *= 2;
    adj_out.build_bidirection(edges);
    // adj_out.build(edges);
    // adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();

    int sum_d = 0;      for (auto deg : degrees)    sum_d += deg;
    LOG(INFO) << (double) sum_d / num_vertices << endl;
    data = ofstream("1_2", ios::app);

    random_tag(p);
    bfs_walk(p);
    // random_tag(num_vertices);
    // bfs_walk(num_vertices);
    ssize_t c = get_uncovered_edge_cnt();
    LOG(INFO) << "c = " << c << endl;

    int all2 = 0;
    for (int i = 1; i <= p; ++ i)
        cout << i << " : " << global_tag_distribute[i] << endl, all2 += global_tag_distribute[i];
    cout << "VERTEX_CNT : " << num_vertices << endl;
    cout << "ALL_TAG_CNT : " << all2 << endl;
    cout << "RF : " << (double) all2 / num_vertices << endl;
    size_t max_occupied = *std::max_element(global_tag_distribute.begin() + 1, global_tag_distribute.end());
    cout << "BALANCE : " << (double)max_occupied / ((double)num_vertices / p) << endl;

    data << "p : " << p << endl;
    data << "RF : " << (double) all2 / num_vertices << endl;
    data << "BALANCE : " << (double)max_occupied / ((double)num_vertices / p) << endl;

    union_tag();
}

void TagPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer compute_timer;

    LOG(INFO) << "partitioning...";
    compute_timer.start();
    
}

bool TagPartitioner::seed_check(vid_t seed_id)
{
    const int INF = num_vertices, MAX_DIST = 5, degree_threshold = 5;
    if (degrees[seed_id] > degree_threshold) return false;

    vector<bool> vis(num_vertices);
    vector<vid_t> dist(num_vertices, INF);         dist[seed_id] = 0;
    unordered_set<vid_t> seed_set;
    queue<array<vid_t, 2>> q;           // [curr_dist, curr_vertex]

    q.push({0, seed_id});
    while (q.size()) 
    {
        auto [cd, c_id] = q.front();      q.pop();
        if (cd > dist[c_id])              continue;
        seed_set.insert(c_id);
        // if (seed_set.size() > num_vertices / 1000)          goto final_check;
        if (dist[c_id] >= MAX_DIST)                 continue;
        // if (bfs_src.count(c_id))            return false;

        for (auto &i : adj_out[c_id]) {
            vid_t to_id = edges[i.v].second;
            if (degrees[to_id] > degree_threshold) continue;
            if (dist[to_id] > dist[c_id] + 1) {
                dist[to_id] = dist[c_id] + 1;
                q.push({dist[to_id], to_id});
            }
        }
    }

        // if (seed_set.size() < 17)     return false;
        // LOG(INFO) << "Size of all_vertex in dist " << MAX_DIST << " is: " << seed_set.size() << endl;
    for (auto vid : bfs_src)    
        if (seed_set.count(vid))   
            return false;
    return true;
}
void TagPartitioner::random_tag(size_t seed_cnt)
{
    Timer seed_timer;   seed_timer.start();
    vector<vid_t> vertex_seq(num_vertices, 0);
    iota(all(vertex_seq), 0);
    shuffle(all(vertex_seq), gen);

    global_tag_distribute[0] = 1e18;
    
    for (int i = 0; i < num_vertices; ++ i)
    {
        if (bfs_src.size() >= seed_cnt) break;
        int vertex_id = vertex_seq[i];
        
        if (adj_out[vertex_id].size() <= 1)      continue;
        if (!seed_check(vertex_id))                                         continue;

        int tag = (bfs_src.size() % p) + 1;
        bfs_src.insert(vertex_id);
        assign_tag(vertex_id, tag);
        // LOG(INFO) << "vertex, tag " << vertex << ' ' << tag;
    }

    vector<size_t> tag_size(p + 1, num_vertices / p / 2), c_tag_size(p + 1, 0);
    vector<bool> st(num_vertices + 1, false);
    for (int b = 1; b <= num_vertices % p; ++ b)    tag_size[b] ++;
    assert(accumulate(tag_size.begin() + 1, tag_size.end(), 0) == num_vertices);

    // for (vid_t u = 0; u < num_vertices; ++ u)  S.insert(u);

    // for (int b = 1; b <= p - 1; ++ b) {
    //     std::cerr << b << ", " << S.size() << endl;
    //     queue<vid_t> q;
    //     while (occupied[b] < capacity) {
    //         vid_t d, vid;
    //         if (!q.size()) {
    //             if (!get_free_vertex(vid)) {
    //                 DLOG(INFO) << "partition " << b
    //                            << " stop: no free vertices";
    //                 break;
    //             }
    //             assert(S.count(vid));
    //         } else {
    //             vid = q.front();    q.pop();
    //         }

    //         if (!S.count(vid))       continue;
    //         S.erase(vid);
    //         assign_tag(vid, b);
    //         ++ occupied[b];
    //         for (auto &i : adj_out[vid]) {
    //             vid_t nid = edges[i.v].second;
    //             if (!S.count(nid))  continue;
    //             q.push(nid);
    //         }
    //     }
    // }
    // std::cerr << p << ", " << S.size() << endl;
    // for (auto vid : S) {
    //     assign_tag(vid, p);
    //     ++ occupied[p];
    // }
    // cerr << global_tag_distribute[p] << endl;

    // for (int b = 1; b <= p; ++ b) {
    //     queue<vid_t> q;     q.push(*S.begin());
    //     LOG(INFO) << b << ' ' << q.size() << ' ';
    //     while (q.size()) {
    //         auto cid = q.front();   q.pop();

    //         if (!S.count(cid))  continue;
    //         S.erase(cid);

    //         assign_tag(cid, b);
    //         ++ c_tag_size[b];
    //         if (c_tag_size[b] >= tag_size[b]) break;
    //         for (auto &i : adj_out[cid]) {
    //             vid_t nid = edges[i.v].second;
    //             if (!S.count(nid))  continue;
    //             q.push(nid);
    //         }
    //         assert(q.size() != 0);
    //         if (c_tag_size[b] < tag_size[b] && q.size() == 0) {
    //             q.push(*S.begin());
    //             cerr << b << ' ' << *S.begin() << ' ' << c_tag_size[b] << endl;
    //         }
    //     }
    //     cout << c_tag_size[b] << ' ' << S.size() << endl;
    // }
    
    // queue<array<vid_t, 2>> q;   
    // for (auto uid : bfs_src) {
    //     q.push({uid, (vid_t)q.size() + 1});
    // }
    
    // while (q.size()) {
    //     auto [cid, ctag] = q.front();   q.pop();

    //     if (c_tag_size[ctag] >= tag_size[ctag]) {
    //         // q.push({cid, (vid_t)(ctag == p ? 1 : ctag + 1)});
    //         for (int b = 1; b <= p; ++ b) if (b != ctag && c_tag_size[b] < tag_size[b]) {
    //             q.push({cid, (vid_t)b});
    //             break;
    //         }
    //         continue;
    //     }
    //     ++ c_tag_size[ctag];

    //     if (st[cid])        continue;
    //     st[cid] = true;

    //     assign_tag(cid, ctag);
    //     bfs_src.insert(cid);
    //     for (auto &i : adj_out[cid]) {
    //         vid_t nid = edges[i.v].second;
    //         if (st[nid])  continue;
    //         q.push({nid, ctag});
    //     }
    // }
    
    seed_timer.stop();
    LOG(INFO) << "time used for seed generation: " << seed_timer.get_time();
}

void 
TagPartitioner::bfs_walk(size_t seed_cnt)
{
    size_t degree_1_vertex_cnt = 0;
    threshold = num_vertices / p * 10000;
    LOG(INFO) << "threshold " << threshold;
    tag_valid = vector<bool>(p + 5, true);

    vector<vid_t> degree_1_vertex;
    for (vid_t u = 0; u < num_vertices; ++ u)   
        if (degrees[u] == 1)        degree_1_vertex.push_back(u);
    degree_1_vertex_cnt = degree_1_vertex.size();
    LOG(INFO) << "degree_1_vertex_cnt " << degree_1_vertex_cnt;

    queue<vid_t> q;
    for (const auto &v : bfs_src)
        q.push(v);

    vid_t covered_cnt = 0;
    vector<vid_t> next_round_vertex;
    edge_covered = vector<bool>(num_edges, false);

    // for (int round = 1; round <= 20; ++ round)
    for (int round = 1; ; ++ round)
    {
        cout << ("========== ROUND " + to_string(round) + " ==========\n");

        LOG(INFO) << "q.size() " << q.size();
        LOG(INFO) << "covered_cnt " << covered_cnt;
        LOG(INFO) << "num_vertices - degree_1_vertex_cnt " << num_vertices - degree_1_vertex_cnt;
        LOG(INFO) << "next_round_vertex.size() " << next_round_vertex.size();
        // int max_tag = max_element(global_tag_distribute.begin() + 1, global_tag_distribute.end()) - global_tag_distribute.begin();
        // LOG(INFO) << "max_tag.size() " << max_tag << ' ' << global_tag_distribute[max_tag];
        for (int b = 1; b <= p; ++ b) {
            cout << b << ' ' << global_tag_distribute[b] << endl;
        }
        cout << "ALL: " << accumulate(global_tag_distribute.begin() + 1, global_tag_distribute.end(), 0) << endl;

        if (covered_cnt == num_vertices - degree_1_vertex_cnt)  break;

        vector<bool> vis(num_vertices, false);

        if (round != 1)
        {
            // if (next_round_vertex.size() <= 5)    break;
            // Test for bfs order: degree wise & tag count wise
            // sort(all(next_round_vertex), [&](int a, int b) {return degree[a] < degree[b];} );
            // sort(all(next_round_vertex), [&](int a, int b) {return vertex2tag[a].popcount() > vertex2tag[b].popcount();} );
            shuffle(all(next_round_vertex), gen);
            for (const auto & v : next_round_vertex)    q.push(v);
        }
        next_round_vertex.clear();

        while (q.size())
        {
            vid_t top_id = q.front();  q.pop();

            if (vis[top_id])    continue;
            vis[top_id] = 1;

            // current vertex cannot cover neighbors yet
            if (!can_cover[top_id] && degrees[top_id] > 1)
            {
                int candidate_tag = choose_tag(top_id);

                if (vertex2tag[top_id].get(candidate_tag) != 1) 
                    assign_tag(top_id, candidate_tag);
                
                // update neighbor's edge_covered[] 
                bool all_neighbor_covered = true;

                for (auto &i : adj_out[top_id]) {
                    if (edge_covered[i.v] || edge_covered[opposite(i.v)])  continue;

                    vid_t nid = edges[i.v].second;
                    bool n_covered = false;
                    if (degrees[nid] == 1) {
                        n_covered = true;
                        continue;
                    }
                    for (int b = 1; b <= p; ++ b) {
                        if (vertex2tag[nid].get(b) == 1 && vertex2tag[top_id].get(b) == 1) {
                            n_covered = true;
                            edge_covered[i.v] = edge_covered[opposite(i.v)] = true;
                        }
                    }
                    all_neighbor_covered &= n_covered;
                }
                
                if (all_neighbor_covered) 
                    ++ covered_cnt, can_cover[top_id] = 1;
                else
                    next_round_vertex.push_back(top_id);
            }
            
            for (auto &i : adj_out[top_id])
            {
                vid_t to_id = edges[i.v].second;
                if (vis[to_id]) continue;
                q.push(to_id);
            }
        }
    }

    size_t c = 0, cc = 0;
    for (size_t eid = 0; eid < num_edges / 2; ++ eid) {
        bool covered = false;
        vid_t u = edges[eid].first, v = edges[eid].second;
        for (int b = 1; b <= p; ++ b) {
            if (vertex2tag[u].get(b) && vertex2tag[v].get(b)) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            ++ c;
            if (degrees[u] == 1 || degrees[v] == 1) ++ cc;
        }
    }
    LOG(INFO) << "c = " << c << endl;
    LOG(INFO) << "cc = " << cc << endl;

    for (auto u : degree_1_vertex) {
        int candidate_tag = choose_tag(u);
        assign_tag(u, candidate_tag);
    }

    // size_t cnt = 0;
    // for (vid_t vid = 0; vid < num_vertices; ++ vid) if (!all_neighbor_covered(vid)) {
    //     cerr << ++ cnt << ' ' << vid << ' ' << degrees[vid] << endl;
    //     cerr << "vid: " << vid << ": ";
    //     for (int b = 1; b <= p; ++ b)   cerr << vertex2tag[vid].get(b); cerr << endl;
    //     for (auto &i : adj_out[vid]) {
    //         vid_t uid = edges[i.v].second;
    //         cerr << "uid: " << uid << ": ";
    //         for (int b = 1; b <= p; ++ b)   cerr << vertex2tag[uid].get(b); cerr << endl;
    //     }
    //     int candidate = choose_tag(vid);
    //     cerr << "candidate: " << candidate << endl;
    // }
}

inline int
TagPartitioner::choose_tag(vid_t uid)
{
    vector<array<vid_t, 2>> neighbor_tag_cnt(p + 1);
    for (vid_t b = 1; b <= p; ++ b)   neighbor_tag_cnt[b] = {0, b};

    for (auto &i : adj_out[uid]) {
        if (edge_covered[i.v] || edge_covered[opposite(i.v)])            continue;

        vid_t vid = edges[i.v].second;
        const auto & v2t = vertex2tag[vid];
        if (!v2t.empty()) {
            for (int b = 1; b <= p; ++ b)
                neighbor_tag_cnt[b][0] += (v2t.get(b) == 1);
        }
    }
    sort(neighbor_tag_cnt.begin() + 1, neighbor_tag_cnt.end());
    reverse(neighbor_tag_cnt.begin() + 1, neighbor_tag_cnt.end());
    // choose a tag for current vertex

    int candidate_tag = 0, c_tag_cnt = 0;

    for (int b = 1; b <= p; ++ b) {
        auto [tag_cnt, tag] = neighbor_tag_cnt[b];
        if (tag_valid[tag] && vertex2tag[uid].get(tag) == 0) {
            if (candidate_tag == 0) {
                candidate_tag = tag, 
                c_tag_cnt = tag_cnt;
                continue;
            }

            if (tag_cnt >= c_tag_cnt)
            {
                if (tag_cnt == c_tag_cnt
                && global_tag_distribute[b] < global_tag_distribute[candidate_tag]) {
                    candidate_tag = tag;
                    c_tag_cnt = tag_cnt;
                }
                else if (tag_cnt > c_tag_cnt) {
                    candidate_tag = tag;
                    c_tag_cnt = tag_cnt;
                }
            }
        }
    }
    return candidate_tag;
}

inline void 
TagPartitioner::assign_tag(vid_t uid, int candidate_tag)
{
    if (++ global_tag_distribute[candidate_tag] >= threshold) {
        tag_valid[candidate_tag] = false;
        // cout << global_tag_distribute[candidate_tag] << ' ' << threshold << ' ' << tag_valid[candidate_tag] << endl;
    }
    vertex2tag[uid].set_bit(candidate_tag);
    bucket[candidate_tag].insert(uid);
}

void 
TagPartitioner::union_tag()
{   
    vector<int> new_tag(p + 1);
    int tag_sum = 0;

    // // What if merging tags by intersection size?
    // vector<vector<int>> intersection(p + 1, vector<int>(p + 1));
    // for (vid_t uid = 0; uid < num_vertices; ++ uid) {
    //     vector<vid_t> v;
    //     for (int b = 1; b <= p; ++ b) {
    //         if (bucket[b].count(uid))   v.push_back(b);
    //     }
    //     if (v.size() <= 1)      continue;

    //     for (auto b : v)
    //         for (auto c : v) if (b != c) {
    //             intersection[b][c] ++;
    //         }
    // }
    // // for (int b = 1; b <= p; ++ b)
    // //     for (int c = 1; c <= p; ++ c) 
    // //         cout << intersection[b][c] << " \n"[c == p];

    // KM km;  km.n = p / 2;
    // for (int b = 1; b <= p / 2; ++ b) 
    //     for (int c = p / 2 + 1; c <= p; ++ c) 
    //         km.w[b][c - p / 2] = intersection[b][c];
        
    // km.solve();

    // for (int c = p / 2 + 1; c <= p; ++ c) {
    //     int matchc = km.matchb[c - p / 2];
        
    //     LOG(INFO) << "tag_c, matchc: " << c << ' ' << matchc << endl;

    //     for (const auto& uid : bucket[c])
    //         bucket[matchc].insert(uid);
    //     new_tag[c] = matchc;
    // }

    // for (int b = 1; b <= p / 2; ++ b) {
    //     auto tag_size = bucket[b].size(); 
    //     cout << tag_size << ' ' << b << endl;
    //     tag_sum += tag_size;
    // }

    // for (int b = 1; b <= p; ++ b) {
    //     cout << b << ' ' << tag_valid[b] << endl;
    // }


    vector<array<vid_t, 2>> v(1);
    for (int b = 1; b <= p; ++ b) {
        v.push_back({vid_t(bucket[b].size()), vid_t(b)});
    }

    sort(v.begin() + 1, v.end());
    reverse(v.begin() + 1, v.end());
    // [tag_size, tag]      
    priority_queue<array<vid_t, 2>, vector<array<vid_t, 2>>, greater<array<vid_t, 2>>> pq;  

    const int k = 4;
    for (int b = 1; b <= p / k; ++ b) {
        pq.push({ v[b][0], v[b][1] });
    }

    for (int b = p / k + 1; b <= p; ++ b) {          // b -> a
        auto [tag_size_b, tag_b] = v[b];
        auto [tag_size_a, tag_a] = pq.top();
        pq.pop();

        for (const auto& uid : bucket[tag_b])
            bucket[tag_a].insert(uid);
        new_tag[tag_b] = tag_a;

        pq.push({vid_t(bucket[tag_a].size()), vid_t(tag_a)});
    }

    for (int i = 1; i <= p / k; ++ i) {
        auto [tag_size, tag] = pq.top();    
        cout << tag_size << ' ' << tag << endl;
        tag_sum += tag_size;
        pq.pop();
    }
    cout << (double)tag_sum / num_vertices << endl;

    data << "final p : " << p / k << endl;
    data << (double)tag_sum / num_vertices << endl;
    data << endl;
}
