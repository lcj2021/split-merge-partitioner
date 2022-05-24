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

    random_tag(p * 5);
    bfs_walk(p * 5);

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
void TagPartitioner::random_tag(size_t random_cnt)
{
//     bfs_src = new int[random_cnt + 1];
    // ull tag_cnt[FLAGS_p + 1];

    global_tag_distribute[0] = 1e18;

    while (bfs_src.size() < random_cnt)
    {
        // LOG(INFO) << "bfs_src.size " << bfs_src.size();
        int vertex_id = gen() % num_vertices;
        // LOG(INFO) << "adj_out[vertex_id].size() " << adj_out[vertex_id].size();
        if (adj_out[vertex_id].size() + adj_in[vertex_id].size() == 0)     continue;

        LOG(INFO) << "vertex_id " << vertex_id;

        int tag = gen() % p + 1;
        vertex2tag[vertex_id].set(tag, 1);
        bfs_src.push_back(vertex_id);
        // tag_cnt[tag] ++;

        ++ global_tag_distribute[tag];
    }

    // for (int i = 1; i <= random_cnt; ++ i)
    //     LOG(bfs_src[i]);
}

void 
TagPartitioner::bfs_walk(size_t random_cnt)
{
    size_t all2 = 0, degree_1_vertex_cnt = 0;
    // vector<size_t> curr_vertex_neighbor_tag_cnt(p + 1);
    size_t* curr_vertex_neighbor_tag_cnt = new size_t[p + 1];

    for (vid_t u = 0; u < num_vertices; ++ u)   
        if (degrees[u] == 1) 
            ++ degree_1_vertex_cnt;
    LOG(INFO) << "degree_1_vertex_cnt " << degree_1_vertex_cnt;

    queue<vid_t> q;
    for (const auto &v : bfs_src)
        q.push(v);

    vid_t covered_cnt = 0;
    vector<vid_t> next_round_vertex;
    vector<bool> st(num_edges, false);

    // for (int round = 1; round <= 20; ++ round)
    for (int round = 1; ; ++ round)
    {
        cout << ("========== ROUND " + to_string(round) + " ==========\n");

        // if (q.empty())      break;

        LOG(INFO) << "q.size() " << q.size();
        LOG(INFO) << "covered_cnt " << covered_cnt;
        LOG(INFO) << "num_vertices - degree_1_vertex_cnt " << num_vertices - degree_1_vertex_cnt;
        LOG(INFO) << "next_round_vertex.size() " << next_round_vertex.size();

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
            vid_t top = q.front();
            q.pop();

            if (vis[top])    continue;
            vis[top] = 1;

            // LOG(INFO) << "top " << top;

            // current vertex cannot cover neighbors yet
            if (!can_cover[top] && degrees[top] > 1)
            {
                // fill (all(curr_vertex_neighbor_tag_cnt), 0);
                memset(curr_vertex_neighbor_tag_cnt, 0, sizeof(size_t) * (p + 1));

                for (auto &i : adj_out[top])
                {
                    if (st[i.v])            continue;
                    // if (!edges[i.v].valid())      continue;

                    vid_t to = edges[i.v].second;
                    const auto & v2t = vertex2tag[to];
                    if (!v2t.empty())
                    {
                        for (int b = 1; b <= p; ++ b)
                            if (v2t.get(b) == 1)
                                ++ curr_vertex_neighbor_tag_cnt[b];
                    }
                }
                

                // for (int i = 0; i <= p; ++ i) printf("%d, ", curr_vertex_neighbor_tag_cnt[i]); puts("");

                // choose a tag for current vertex
                // vid_t max_tag_cnt = *max_element(curr_vertex_neighbor_tag_cnt.begin() + 1, curr_vertex_neighbor_tag_cnt.end());
                vid_t max_tag_cnt = *max_element(curr_vertex_neighbor_tag_cnt + 1, curr_vertex_neighbor_tag_cnt + 1 + FLAGS_p);
                // LOG(max_tag_cnt);

                vid_t candidate_tag_cnt = 0;
                int candidate_tag = 0;
                for (int b = 1; b <= p; ++ b)
                {
                    if (curr_vertex_neighbor_tag_cnt[b] == max_tag_cnt)
                    {
                        if (global_tag_distribute[b] < global_tag_distribute[candidate_tag])
                        {
                            candidate_tag = b;
                            candidate_tag_cnt = curr_vertex_neighbor_tag_cnt[b];
                        }
                    }
                }

                if (candidate_tag > p || candidate_tag < 1)   
                {
                    // LOG(top);
                }
                else
                {
                    if (vertex2tag[top].get(candidate_tag) != 1)
                    {
                        ++ global_tag_distribute[candidate_tag];
                        vertex2tag[top].set_bit(candidate_tag);
                    }
                }
                

                // update neighbor's st[] 
                bool all_neighbor_covered = true;
                int neighbor_uncovered_cnt = 0;

                for (auto &i : adj_out[top])
                {
                    // if (!edges[i.v].valid())      continue;
                    if (st[i.v])      continue;

                    vid_t to = edges[i.v].second;
                    if (degrees[to] == 1)    
                    {
                        st[i.v] = 1;
                        // edges[i.v].remove();
                        continue;
                    }

                    all_neighbor_covered = false;
                    ++ neighbor_uncovered_cnt;

                    if (vertex2tag[to].get(candidate_tag))
                    {
                        st[i.v] = 1;
                        // edges[i.v].remove();
                        // continue;
                    }
                }
                

                if (all_neighbor_covered)   
                {
                    ++ covered_cnt, can_cover[top] = 1;
                }
                else
                {
                    // if (neighbor_uncovered_cnt == 1)    
                    //     LOG(top);
                    next_round_vertex.push_back(top);
                }

            }
            

            for (auto &i : adj_out[top])
            {
                vid_t to = edges[i.v].second;
                if (vis[to]) continue;
                q.push(to);
            }
            
        }
    }

}