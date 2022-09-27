#include "tag_partitioner.hpp"
#include "conversions.hpp"
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
    assigned_edges = 0;
    // occupied.assign(p, 0);
    adj_out.resize(num_vertices + 1);
    adj_in.resize(num_vertices + 1);
    global_tag_distribute.assign(p + 1, 0);
    vertex2tag.assign(num_vertices + 1, dense_bitset(p + 1));
    can_cover.assign(num_vertices + 1, false);

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

    // for (auto [u, v] : edges)   cout << u << ' ' << v << endl;  cout << "================" << endl;
    // for (int i = 1; i <= num_vertices; ++ i)
    // {
    //     vid_t v = i - 1;
    //     for (auto &i : adj_out[v])
    //     {
    //         vid_t to = edges[i.v].second + 1;
    //         cout << v + 1 << ' ' << to << endl;
    //     }
    //     cout << endl;
    // }
    int sum_d = 0;      for (auto deg : degrees)    sum_d += deg;
    cout << (double) sum_d / num_vertices << endl;

    random_tag(50* p);
    bfs_walk(50 * p);

    int all2 = 0;
    for (int i = 1; i <= p; ++ i)
        cout << i << " : " << global_tag_distribute[i] << endl, all2 += global_tag_distribute[i];
    cout << "VERTEX_CNT : " << num_vertices << endl;
    cout << "ALL_TAG_CNT : " << all2 << endl;
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
    const int MAX_DEG = 17 / 2, MIN_DEG = 3, INF = num_vertices, MAX_DIST = 15;
    if (degrees[seed_id] > MAX_DEG || degrees[seed_id] < MIN_DEG)             return false;

    vector<bool> vis(num_vertices);
    vector<vid_t> dist(num_vertices, INF);         dist[seed_id] = 0;
    unordered_set<vid_t> all_v;
    queue<array<vid_t, 2>> q;           // [curr_dist, curr_vertex]

    q.push({0, seed_id});
    while (q.size()) 
    {
        auto [cd, c_id] = q.front();      q.pop();
        if (cd > dist[c_id])              continue;
        all_v.insert(c_id);
        // if (all_v.size() > num_vertices / 1000)          goto final_check;
        if (dist[c_id] >= MAX_DIST)                 continue;
        // if (bfs_src.count(c_id))            return false;

        for (auto &i : adj_out[c_id])
        {
            vid_t to_id = edges[i.v].second;
            if (cd <= 3 && degrees[to_id] > 17)                                continue;
            if (cd > 3 && degrees[to_id] > MAX_DEG || degrees[to_id] < MIN_DEG)      continue;
            if (dist[to_id] > dist[c_id] + 1) {
                dist[to_id] = dist[c_id] + 1;
                q.push({dist[to_id], to_id});
            }
        }
    }

    final_check: 
        // if (all_v.size() < 17)     return false;
        // LOG(INFO) << "Size of all_vertex in dist " << MAX_DIST << " is: " << all_v.size() << endl;
        for (auto vid : bfs_src)    
            if (all_v.count(vid))   
                return false;
        return true;
}
void TagPartitioner::random_tag(size_t random_cnt)
{
    // ull tag_cnt[FLAGS_p + 1];
    Timer seed_timer;
    seed_timer.start();

    global_tag_distribute[0] = 1e18;
    while (bfs_src.size() < random_cnt)
    {
        // LOG(INFO) << "bfs_src.size " << bfs_src.size();
        int vertex = (gen() % num_vertices) + 1, vertex_id = vertex - 1;
        // LOG(INFO) << "adj_out[vertex_id].size() " << adj_out[vertex_id].size();
        if (adj_out[vertex_id].size() + adj_in[vertex_id].size() == 0)      continue;

        if (!seed_check(vertex_id))                                         continue;

        LOG(INFO) << "vertex " << vertex;

        int tag = gen() % p + 1;
        vertex2tag[vertex_id].set(tag, 1);
        bfs_src.insert(vertex_id);
        // tag_cnt[tag] ++;

        ++ global_tag_distribute[tag];
    }

    // for (int i = 1; i <= random_cnt; ++ i)
    //     LOG(bfs_src[i]);
    seed_timer.stop();
    LOG(INFO) << "time used for seed generation: " << seed_timer.get_time();
}

void 
TagPartitioner::bfs_walk(size_t random_cnt)
{
    size_t all2 = 0, degree_1_vertex_cnt = 0;
    // const vid_t threshold = round((double)num_vertices / p * 10);
    threshold = num_vertices / p * 2;
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

        // if (q.empty())      break;

        LOG(INFO) << "q.size() " << q.size();
        LOG(INFO) << "covered_cnt " << covered_cnt;
        LOG(INFO) << "num_vertices - degree_1_vertex_cnt " << num_vertices - degree_1_vertex_cnt;
        LOG(INFO) << "next_round_vertex.size() " << next_round_vertex.size();
        int max_tag = max_element(global_tag_distribute.begin() + 1, global_tag_distribute.end()) - global_tag_distribute.begin();
        LOG(INFO) << "max_tag.size() " << max_tag << ' ' << global_tag_distribute[max_tag];

        if (covered_cnt == num_vertices - degree_1_vertex_cnt)  break;

        vector<bool> vis(num_vertices, false);

        if (round != 1)
        {
            // Test for bfs order: degree wise & tag count wise
            // sort(all(next_round_vertex), [&](int a, int b) {return degree[a] < degree[b];} );
            // sort(all(next_round_vertex), [&](int a, int b) {return vertex2tag[a].popcount() > vertex2tag[b].popcount();} );
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
                if (candidate_tag == 0)     exit(0);

                if (vertex2tag[top_id].get(candidate_tag) != 1) 
                    assign_tag(top_id, candidate_tag);
                
                // update neighbor's edge_covered[] 
                bool all_neighbor_covered = true;
                int neighbor_uncovered_cnt = 0;

                for (auto &i : adj_out[top_id])
                {
                    if (edge_covered[i.v])      continue;

                    vid_t to_id = edges[i.v].second;
                    if (degrees[to_id] == 1)    
                    {
                        edge_covered[i.v] = 1;
                        continue;
                    }

                    all_neighbor_covered = false;
                    ++ neighbor_uncovered_cnt;

                    if (vertex2tag[to_id].get(candidate_tag)) {
                        edge_covered[i.v] = 1;
                    }
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

    for (auto u : degree_1_vertex) {
        int candidate_tag = choose_tag(u);
        assign_tag(u, candidate_tag);
    }
}

int
TagPartitioner::choose_tag(vid_t uid)
{
    vector<vid_t> curr_vertex_neighbor_tag_cnt(p + 1, 0);

    for (auto &i : adj_out[uid])
    {
        if (edge_covered[i.v])            continue;

        vid_t to_id = edges[i.v].second;
        const auto & v2t = vertex2tag[to_id];
        if (!v2t.empty())
        {
            for (int b = 1; b <= p; ++ b)
                if (v2t.get(b) == 1)
                    ++ curr_vertex_neighbor_tag_cnt[b];
        }
    }
    
    // choose a tag for current vertex
    int max_tag = max_element(curr_vertex_neighbor_tag_cnt.begin() + 1, curr_vertex_neighbor_tag_cnt.end()) - curr_vertex_neighbor_tag_cnt.begin();
    vid_t max_tag_cnt = curr_vertex_neighbor_tag_cnt[max_tag];

    int candidate_tag = max_tag;
    for (int b = 1; b <= p; ++ b) if (tag_valid[b]) {
        if (curr_vertex_neighbor_tag_cnt[b] == max_tag_cnt)
        {
            if (global_tag_distribute[b] < global_tag_distribute[candidate_tag]) {
                candidate_tag = b;
            }
        }
    }
    return candidate_tag;
}

void 
TagPartitioner::assign_tag(vid_t uid, int candidate_tag)
{
    ++ global_tag_distribute[candidate_tag];
    if (global_tag_distribute[candidate_tag] >= threshold)  
        tag_valid[candidate_tag] = false;
    vertex2tag[uid].set_bit(candidate_tag);
}
