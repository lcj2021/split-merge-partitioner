#include "tag_partitioner.hpp"
#include "conversions.hpp"
#include "kcore.hpp"
#include "dsu.hpp"
using namespace std;
#define all(x) x.begin(), x.end()
#define rall(x) x.rbegin(), x.rend()

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
    occupied.assign(p, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    global_tag_distribute.assign(p, 0);
    vertex2tag.assign(num_vertices, dense_bitset(p));
    can_cover.assign(num_vertices, false);
    tag_valid.assign(p, false);
    bucket.resize(p);
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
    threshold = num_edges / 2 / p * BALANCE_RATIO;
    for (vid_t u = 0; u < num_vertices; ++ u)   
        if (degrees[u] == 1)        degree_1_vertex.push_back(u);

    int sum_d = 0;      for (auto deg : degrees)    sum_d += deg;
    LOG(INFO) << (double) sum_d / num_vertices << endl;
    data = ofstream("EXP_LOG.txt", ios::app);

    random_tag(p);
    bfs_walk(p);
    // size_t c = get_uncovered_edge_cnt();
    size_t cnt = 0, c1 = 0;
        for (size_t eid = 0; eid < num_edges / 2; ++ eid) {
            bool covered = false;
            vid_t u = edges[eid].first, v = edges[eid].second;
            for (int b = 0; b < p; ++ b) {
                if (vertex2tag[u].get(b) == 1 && vertex2tag[v].get(b) == 1) {
                    covered = true;
                    break;
                }
            }
            if (!covered) {
                ++ cnt;
                if (degrees[u] == 1 || degrees[v] == 1)     ++ c1;
                vid_t vid = all_neighbor_covered(u) ? u : v;

                int candidate = choose_tag(vid, true, true);
                cerr << candidate << endl;
                exit(0);
            }
        }
    LOG(INFO) << "cnt = " << cnt << endl;
    data << "cnt = " << cnt << endl;
    LOG(INFO) << "c1 = " << c1 << endl;
    
    size_t all2 = accumulate(all(global_tag_distribute), 0);
    cout << "VERTEX_CNT : " << num_vertices << endl;
    cout << "ALL_TAG_CNT : " << all2 << endl;
    cout << "RF : " << (double) all2 / num_vertices << endl;
    for (int b = 0; b < p; ++ b) {
        cerr << b << ' ' << global_tag_distribute[b] << ' ' << occupied[b] << endl;
    }
    size_t max_occupied = *std::max_element(occupied.begin(), occupied.end());
    cerr << accumulate(all(occupied), 0) << ' ' << num_edges / 2 << endl;
    cout << "BALANCE : " << (double)max_occupied / ((double)num_edges / 2 / p) << endl;

    data << "p : " << p << endl;
    data << "RF : " << (double) all2 / num_vertices << endl;
    data << "BALANCE : " << (double)max_occupied / ((double)num_edges / 2 / p) << endl;

    // union_tag();

    vid_t mirror = 0, cut_vertex = 0;
    map<vid_t, vid_t> deg_cnt;
    for (vid_t vid = 0; vid < num_vertices; ++ vid) {
        int cnt = 0;
        for (int b = 0; b < p; ++ b) {
            cnt += vertex2tag[vid].get(b);
        }
        if (cnt > 1) {
            deg_cnt[degrees[vid]] ++;
            mirror += cnt - 1;
            cut_vertex ++;
        }
    }

    ofstream out(basefilename + ".tag.cut", ios::out);
    out << "Mirror: " << mirror << endl;
    out << "Cut_vertex: " << cut_vertex << endl;
    for (auto [deg, cnt] : deg_cnt) {
        out << deg << '\t' << cnt << endl;
    }


    vid_t cut_kcore_vertex = 0;
    for (auto vid : kcore_s) {
        int cnt = 0;
        for (int b = 0; b < p; ++ b) {
            cnt += vertex2tag[vid].get(b);
        }
        if (cnt > 1) {
            cut_kcore_vertex ++;
        }
    }
    out << "kcore_k: " << kcore_k << '\n';
    out << "Cut_kcore_vertex: " << cut_kcore_vertex << endl;
    out << "Cut_kcore_vertex ratio: " << (double)cut_kcore_vertex / kcore_s.size() << endl;
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

    vector<size_t> tag_size(p, num_vertices / p / 2), c_tag_size(p, 0);
    vector<bool> st(num_vertices, false);
    for (int b = 0; b < num_vertices % p; ++ b)    tag_size[b] ++;

    seeded.assign(num_vertices, false);
    capacity = (double)(num_vertices - degree_1_vertex.size()) / p / 2.25 + 1;

    vector<queue<vid_t>> q(p);
    for (int b = 0; b < p - 1; ++ b) {
        cerr << b << endl;
        while (global_tag_distribute[b] < capacity) {
            vid_t vid;
            if (!q[b].size()) {
                while (!get_free_vertex(vid));
                LOG(INFO) << "partition " << b
                               << " stop: no free vertices";
            } else {
                vid = q[b].front();    q[b].pop();
            }

            if (seeded[vid])            continue;
            if (degrees[vid] == 1)      continue;
            seeded[vid] = true;
            assign_tag(vid, b, true);
            for (auto &i : adj_out[vid]) {
                vid_t nid = edges[i.v].second;
                if (seeded[nid])        continue;
                if (adj_out[nid].size() >
                    2 * average_degree) continue;
                q[b].push(nid);
            }
        }
    }
    cerr << p << endl;
    for (vid_t vid = 0, cnt = 0; vid < num_vertices && cnt < capacity; ++ vid) if (degrees[vid] != 1 && !seeded[vid]) {
        assign_tag(vid, p - 1, true);
        ++ cnt;
    }

    // kcore_t kcore(num_vertices);
    // LOG(INFO) << edges.size();
    // LOG(INFO) << num_vertices;
    // for (auto [u, v] : edges) {
    //     kcore.AddEdge(u, v);
    // }

    // kcore_k = 2;
    // kcore_s = kcore.solve(kcore_k * average_degree);
    // ofstream out("kcore.txt", ios::app);
    // out << "\n==============================\n";
    // out << "kcore_k: " << kcore_k << '\n';
    // out << "kcore_s.size(): " << kcore_s.size() << '\n';
    // out << "kcore_s.DSU:" << '\n';

    // dsu = DSU(num_vertices);
    // for (auto u : kcore_s) {
    //     for (auto &i : adj_out[u]) {
    //         int v = edges[i.v].second;
    //         if (!kcore_s.count(v)) continue;
    //         dsu.merge(u, v);
    //     }
    // }
    // priority_queue<array<int, 2>, vector<array<int, 2>>, greater<array<int, 2>> >  pq;
    // int wcc_id = 0;
    // vector<int> id2dsup;
    // for (auto u : kcore_s) {
    //     if (u == dsu.find(u)) {
    //         out << u << ' ' << dsu.size(u) << '\n';
    //         pq.push({dsu.size(u), wcc_id ++});
    //         id2dsup.push_back(u);
    //     }
    // }
    // DSU new_id(wcc_id);
    // while (pq.size() > p - 1) {
    //     auto [sza, ida] = pq.top(); pq.pop();
    //     auto [szb, idb] = pq.top(); pq.pop();
    //     new_id.merge(ida, idb);
    //     pq.push({sza + szb, new_id.find(ida)});
    // }
    // while (pq.size()) {
    //     auto [sz, id] = pq.top();   pq.pop();
    //     cerr << sz << ' ' << id << endl;
    // }
    // map<int, int> dsup2id;
    // map<int, int> id2qid;
    // int qid = 0;
    // for (int i = 0; i < wcc_id; ++ i) {
    //     // cerr << i << ' ' << id2dsup[i] << ' ' << new_id.find(i) << endl;
    //     if (!id2qid.count(new_id.find(i))) {
    //         id2qid[new_id.find(i)] = qid ++;
    //     }
    //     dsup2id[id2dsup[i]] = id2qid[new_id.find(i)];
    // }

    // for (int i = 0; i < wcc_id; ++ i) {
    //     cerr << i << ' ' << id2dsup[i] << ' ' << new_id.find(i) << ' ' << id2qid[new_id.find(i)] << endl;
    // }


    // vector<queue<vid_t>> q(p);
    // for (auto vid : kcore_s) {
    //     q[id2qid[new_id.find( dsup2id[dsu.find(vid)] )]].push(vid);
    // }

    // for (int b = 0; b < p - 1; ++ b) {
    //     cerr << b << endl;
    //     while (global_tag_distribute[b] < capacity) {
    //         vid_t vid;
    //         if (!q[b].size()) {
    //             while (!get_free_vertex_kcore(vid));
    //             LOG(INFO) << "partition " << b
    //                            << " stop: no free vertices";
    //         } else {
    //             vid = q[b].front();    q[b].pop();
    //         }

    //         if (seeded[vid])            continue;
    //         if (degrees[vid] == 1)      continue;
    //         seeded[vid] = true;
    //         assign_tag(vid, b, true);
    //         for (auto &i : adj_out[vid]) {
    //             vid_t nid = edges[i.v].second;
    //             if (seeded[nid])        continue;
    //             q[b].push(nid);
    //         }
    //     }
    // }
    // cerr << p << endl;
    // for (vid_t vid = 0, cnt = 0; vid < num_vertices && cnt < capacity; ++ vid) if (degrees[vid] != 1 && !seeded[vid]) {
    //     assign_tag(vid, p - 1, true);
    //     ++ cnt;
    // }


    seed_timer.stop();
    LOG(INFO) << "time used for seed generation: " << seed_timer.get_time();
}

void 
TagPartitioner::bfs_walk(size_t seed_cnt)
{
    size_t degree_1_vertex_cnt = 0;
    LOG(INFO) << "threshold " << threshold;
    tag_valid.assign(p, true);

    degree_1_vertex_cnt = degree_1_vertex.size();
    LOG(INFO) << "degree_1_vertex_cnt " << degree_1_vertex_cnt;

    vid_t covered_cnt = 0;
    vector<vid_t> next_round_vertex;
    vector<vid_t> curr_round_vertex(num_vertices);
    for (vid_t vid = 0; vid < num_vertices; ++ vid) if (degrees[vid] > 1) 
        curr_round_vertex.push_back(vid);

    edge_covered = vector<bool>(num_edges, false);

    for (int round = 1; ; ++ round)
    {
        cout << ("========== ROUND " + to_string(round) + " ==========\n");

        LOG(INFO) << "covered_cnt " << covered_cnt;
        LOG(INFO) << "num_vertices - degree_1_vertex_cnt " << num_vertices - degree_1_vertex_cnt;
        LOG(INFO) << "next_round_vertex.size() " << next_round_vertex.size();
        // int max_tag = max_element(global_tag_distribute.begin(), global_tag_distribute.end()) - global_tag_distribute.begin();
        // LOG(INFO) << "max_tag.size() " << max_tag << ' ' << global_tag_distribute[max_tag];
        for (int b = 0; b < p; ++ b) 
            cout << b << ' ' << global_tag_distribute[b] << endl;
        cout << "ALL: " << accumulate(all(global_tag_distribute), 0) << endl;

        if (covered_cnt == num_vertices - degree_1_vertex_cnt)  break;

        vector<bool> vis(num_vertices, false);

        for (auto vid : next_round_vertex)  curr_round_vertex.push_back(vid);
        // shuffle(all(curr_round_vertex), gen);
        // sort(all(curr_round_vertex), [&](vid_t u, vid_t v) {
        //     return degrees[u] > degrees[v];
        // });
        next_round_vertex.clear();

        for (auto vid : curr_round_vertex) {
            if (!can_cover[vid]) {
                if (all_neighbor_covered(vid)) {
                    ++ covered_cnt, can_cover[vid] = 1;
                    continue;
                }
                int candidate_tag = choose_tag(vid, true, true);
                // assert(candidate_tag >= 0 && candidate_tag < p && "Candidate_tag out of range!");
                if (candidate_tag >= p) {
                    assign_tag(vid, candidate_tag, true);
                    assign_tag(vid, candidate_tag, true);
                }
                else
                    assign_tag(vid, candidate_tag, true);
                
                // update neighbor's edge_covered[]
                bool neighbor_covered = all_neighbor_covered(vid);
                
                if (neighbor_covered) 
                    ++ covered_cnt, can_cover[vid] = 1;
                else
                    next_round_vertex.push_back(vid);
            }
        }
    }

    size_t c = get_uncovered_edge_cnt();
    LOG(INFO) << "c = " << c << endl;

    for (auto u : degree_1_vertex) {
        int candidate_tag = choose_tag(u, false);
        assign_tag(u, candidate_tag, false);
        ++ occupied[candidate_tag];
    }

    // size_t cnt = 0;
    // for (vid_t vid = 0; vid < num_vertices; ++ vid) if (!all_neighbor_covered(vid)) {
    //     ++ cnt;
    //     if (cnt >= 10)  continue;
    //     cerr << cnt << ' ' << vid << ' ' << degrees[vid] << endl;
    //     cerr << "vid: " << vid << ": ";
    //     for (int b = 0; b < p; ++ b)   cerr << vertex2tag[vid].get(b); cerr << endl;
    //     for (auto &i : adj_out[vid]) {
    //         vid_t uid = edges[i.v].second;
    //         cerr << "uid: " << uid << ": ";
    //         for (int b = 0; b < p; ++ b)   cerr << vertex2tag[uid].get(b); cerr << endl;
    //     }
    //     int candidate = choose_tag(vid, 1, 1);
    //     cerr << "candidate: " << candidate << endl;
    // }
    // LOG(INFO) << "cnt: " << cnt << endl;
}

inline int
TagPartitioner::choose_tag(vid_t uid, bool restrict)
{
    vector<array<vid_t, 2>> neighbor_tag_cnt(p);
    for (vid_t b = 0; b < p; ++ b)   neighbor_tag_cnt[b] = {0, b};

    for (auto &i : adj_out[uid]) {
        if (edge_covered[i.v] || edge_covered[opposite(i.v)])            continue;

        vid_t vid = edges[i.v].second;
        const auto & v2t = vertex2tag[vid];
        if (!v2t.empty()) {
            for (int b = 0; b < p; ++ b)
                neighbor_tag_cnt[b][0] += (v2t.get(b) == 1);
        }
    }
    sort(rall(neighbor_tag_cnt));
    // choose a tag for current vertex

    int candidate_tag = -1, c_tag_cnt = 0;

    for (int b = 0; b < p; ++ b) {
        auto [tag_cnt, tag] = neighbor_tag_cnt[b];
        if ((!restrict || tag_valid[tag]) && vertex2tag[uid].get(tag) == 0) {
            if (candidate_tag == -1) {
                candidate_tag = tag, 
                c_tag_cnt = tag_cnt;
                continue;
            }

            if (tag_cnt >= c_tag_cnt) {
                if (tag_cnt == c_tag_cnt
                && occupied[tag] < occupied[candidate_tag]) {
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
TagPartitioner::assign_tag(vid_t uid, int candidate_tag, bool restrict)
{
    // if (global_tag_distribute[candidate_tag] >= threshold) {
    //     tag_valid[candidate_tag] = false;
    //     return;
    // }
    if (restrict && occupied[candidate_tag] >= threshold) {
        tag_valid[candidate_tag] = false;
        return;
    }
    ++ global_tag_distribute[candidate_tag];

    vertex2tag[uid].set_bit(candidate_tag);
    bucket[candidate_tag].insert(uid);
}

void 
TagPartitioner::union_tag()
{   
    vector<int> new_tag(p);
    int tag_sum = 0;

    // // What if merging tags by intersection size?
    // vector<vector<int>> intersection(p, vector<int>(p));
    // for (vid_t uid = 0; uid < num_vertices; ++ uid) {
    //     vector<vid_t> v;
    //     for (int b = 0; b < p; ++ b) {
    //         if (bucket[b].count(uid))   v.push_back(b);
    //     }
    //     if (v.size() <= 1)      continue;

    //     for (auto b : v)
    //         for (auto c : v) if (b != c) {
    //             intersection[b][c] ++;
    //         }
    // }
    // // for (int b = 0; b < p; ++ b)
    // //     for (int c = 1; c <= p; ++ c) 
    // //         cout << intersection[b][c] << " \n"[c == p];

    // KM km;  km.n = p / 2;
    // for (int b = 0; b < p / 2; ++ b) 
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

    // for (int b = 0; b < p / 2; ++ b) {
    //     auto tag_size = bucket[b].size(); 
    //     cout << tag_size << ' ' << b << endl;
    //     tag_sum += tag_size;
    // }

    // for (int b = 0; b < p; ++ b) {
    //     cout << b << ' ' << tag_valid[b] << endl;
    // }


    vector<array<vid_t, 2>> v;
    for (int b = 0; b < p; ++ b) {
        v.push_back({vid_t(bucket[b].size()), vid_t(b)});
    }

    sort(rall(v));
    // [tag_size, tag]      
    priority_queue<array<vid_t, 2>, vector<array<vid_t, 2>>, greater<array<vid_t, 2>>> pq;  

    const int k = 4;
    for (int b = 0; b < p / k; ++ b) {
        pq.push({ v[b][0], v[b][1] });
    }

    for (int b = p / k; b < p; ++ b) {          // b -> a
        auto [tag_size_b, tag_b] = v[b];
        auto [tag_size_a, tag_a] = pq.top();
        pq.pop();

        for (const auto& uid : bucket[tag_b])
            bucket[tag_a].insert(uid);
        new_tag[tag_b] = tag_a;

        pq.push({vid_t(bucket[tag_a].size()), vid_t(tag_a)});
    }

    for (int i = 0; i < p / k; ++ i) {
        auto [tag_size, tag] = pq.top();    
        cout << tag_size << ' ' << tag << endl;
        tag_sum += tag_size;
        pq.pop();
    }
    cout << (double)tag_sum / num_vertices << endl;

    data << "final p : " << p / k << endl;
    data << (double)tag_sum / num_vertices << endl;
}

inline int
TagPartitioner::choose_tag(vid_t uid, bool restrict, bool debug)
{
    vector<array<vid_t, 2>> neighbor_tag_cnt(p);
    for (vid_t b = 0; b < p; ++ b)   neighbor_tag_cnt[b] = {0, b};

    for (auto &i : adj_out[uid]) {
        if (edge_covered[i.v] || edge_covered[opposite(i.v)])            continue;

        vid_t vid = edges[i.v].second;
        const auto & v2t = vertex2tag[vid];
        if (!v2t.empty()) {
            for (int b = 0; b < p; ++ b)
                neighbor_tag_cnt[b][0] += (v2t.get(b) == 1);
        }
    }
    sort(rall(neighbor_tag_cnt));
    // choose a tag for current vertex
    // for (auto [tag_c, tag] : neighbor_tag_cnt) {
    //     cerr << tag_c << ' ' << tag << " valid? " << tag_valid[tag] << endl;
    // }
    int candidate_tag = -1, c_tag_cnt = 0;

    for (int b = 0; b < p; ++ b) {
        auto [tag_cnt, tag] = neighbor_tag_cnt[b];
        if (tag_valid[tag] && vertex2tag[uid].get(tag) == 0) {
            if (candidate_tag == -1) {
                candidate_tag = tag, 
                c_tag_cnt = tag_cnt;
                continue;
            }

            if (tag_cnt >= c_tag_cnt)
            {
                if (tag_cnt == c_tag_cnt
                && occupied[tag] < occupied[candidate_tag]) {
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
    if (candidate_tag < 0 || candidate_tag >= p) {
        for (auto [tag_c, tag] : neighbor_tag_cnt) {
            cerr << tag_c << ' ' << tag << " valid? " << tag_valid[tag] << " has? " << vertex2tag[uid].get(tag) << endl;
        }
        cerr << "Candidate tag: " << candidate_tag << endl;
        exit(0);
    }
    return candidate_tag;
}