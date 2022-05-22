#include <bits/stdc++.h>
#include "dense_bitset.hpp"
#include "util.hpp"
#include "tag_partitioner.hpp"
using namespace std;
#define endl "\n"
#define all(x) x.begin(), x.end()
typedef unsigned long long ull;
const ull N = 1e7 + 10, M = 3e8 + 10;

DEFINE_string(filetype, "edgelist",
              "the type of input file (supports 'edgelist' and 'adjlist')");
DEFINE_int32(p, 10, "number of parititions");
DEFINE_string(filename, "", "the file name of the input graph");

namespace debugger
{
#ifdef LOCAL
    template <typename T>
    void __print_var(string_view name, const T &x)
    {
        std::cerr << name << " = " << x;
    }
    // template <typename T>
    void __print_var(string_view name, const dense_bitset &x)
    {
        std::cerr << name << " = ";
        bool is_first = true;
        for (int i = 0; i <= FLAGS_p; ++ i)
            std::cerr << (is_first ? (is_first = false, "[") : ", ") << x.get(i);
        std::cerr << "]";
    }
    template <typename T>
    void __print_var(string_view name, const vector<T> &x)
    {
        std::cerr << name << " = ";
        bool is_first = true;
        for (auto &ele : x)
            std::cerr << (is_first ? (is_first = false, "[") : ", ") << ele;
        std::cerr << "]";
    }
    template <typename T>
    void __print_var(string_view name, const set<T> &x)
    {
        std::cerr << name << " = ";
        bool is_first = true;
        for (auto &ele : x)
            std::cerr << (is_first ? (is_first = false, "{") : ", ") << ele;
        std::cerr << "}";
    }
    template <typename K, typename V>
    void __print_var(string_view name, const map<K, V> &x)
    {
        std::cerr << name << " = ";
        bool is_first = true;
        for (auto &[k, v] : x)
            std::cerr << (is_first ? (is_first = false, "{") : ", ") << "(" << k << ": " << v << ")";
        std::cerr << "}";
    }
    template <typename T>
    void __log(string_view name, const T &x)
    {
        __print_var(name, x);
        std::cerr << endl;
    }
    template <typename T, typename... Ts>
    void __log(string_view name, const T &x, const Ts &...others)
    {
        size_t pos = name.find(',');
        __print_var(name.substr(0, pos), x);
        std::cerr << ", ";
        __log(name.substr(pos + 1), others...);
    }

#define LOG(args...)\
    {																   \
        std::cerr << "line " << __LINE__ << ": " << __func__ << "(): "; \
        __log(#args, ##args);										   \
    }
#else
#define LOG(...)
#endif
}
using namespace debugger;

// ull h[N], e[M], ne[M], idx;
// bool st[M], can_cover[N];
// ull degree[N];
vector<ull> h(N), e(M), ne(M);
int idx;
vector<bool> st(M), can_cover(M);
vector<ull> degree(N);
vector<ull> tag2size;

ull degree_1_vertex_cnt;
unordered_set<ull> vertex_set;
dense_bitset vertex2tag[N];
bool vis[N];

ull *global_tag_distribute = NULL;
int *bfs_src = NULL;

void add_edge (int a, int b)
{
    e[idx] = b, ne[idx] = h[a], h[a] = idx ++;
    ++ degree[b];
}

void read_graph()
{
    freopen(FLAGS_filename.data(), "r", stdin);
    // memset(h, - 1, sizeof h);
    fill(h.begin(), h.end(), -1);

    string line;
    int line_cnt = 0;
    while (getline(cin, line))
    {
        if (line[0] == '#')     continue;

        if (++ line_cnt % 100000 == 0)  LOG(line_cnt);

        int from = stoi(line.substr(0, line.find("\t")));
        int to = stoi(line.substr(line.find("\t") + 1));
        // LOG(to_string(from) + " -> " + to_string(to));

        if (from == to)         continue;
        vertex_set.insert(from);
        vertex_set.insert(to);
        add_edge (from, to);
        add_edge (to, from);
    }

    for (auto v : vertex_set)   
        if (degree[v] == 1) 
            ++ degree_1_vertex_cnt;
}

void random_tag(int random_cnt)
{
    bfs_src = new int[random_cnt + 1];
    ull tag_cnt[FLAGS_p + 1];
    std::mt19937 gen(std::random_device{}());

    global_tag_distribute[0] = ULLONG_MAX;

    int curr_cnt = 0;
    while (curr_cnt < random_cnt)
    {
        int vertex_id = gen() % N;
        if (h[vertex_id] == -1)     continue;

        ++ curr_cnt;
        int tag = gen() % FLAGS_p + 1;
        vertex2tag[vertex_id].set(tag, 1);
        bfs_src[curr_cnt] = vertex_id;
        tag_cnt[tag] ++;

        ++ global_tag_distribute[tag];
    }

    // for (int i = 1; i <= random_cnt; ++ i)
    //     LOG(bfs_src[i]);
}

void bfs_walk(int random_cnt)
{
    int all2 = 0;
    size_t* curr_vertex_neighbor_tag_cnt = new size_t[FLAGS_p + 1];

    // LOG(vertex_set.size());

    queue<ull> q;
    for (int i = 1; i <= random_cnt; ++ i)
    {
        q.push(bfs_src[i]);
    }

    ull covered_cnt = 0;
    vector<ull> next_round_vertex;

    // for (int round = 1; round <= 2; ++ round)
    for (int round = 1; ; ++ round)
    {
        cout << ("========== ROUND " + to_string(round) + " ==========\n");

        // if (q.empty())      break;

        LOG(q.size());
        LOG(covered_cnt);
        LOG(vertex_set.size() - degree_1_vertex_cnt);
        LOG(next_round_vertex.size());

        if (covered_cnt == vertex_set.size() - degree_1_vertex_cnt)  break;

        memset(vis, 0, sizeof vis);

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
            ull top = q.front();
            q.pop();

            if (vis[top])    continue;
            vis[top] = 1;

            // current vertex cannot cover neighbors yet
            if (!can_cover[top] && degree[top] > 1)
            {
                memset(curr_vertex_neighbor_tag_cnt, 0, sizeof(size_t) * (FLAGS_p + 1));
                for (int i = h[top]; i != -1; i = ne[i])
                {
                    if (st[i])      continue;

                    int to = e[i];
                    const auto & v2t = vertex2tag[to];
                    if (!v2t.empty())
                    {
                        for (int b = 1; b <= FLAGS_p; ++ b)
                            if (v2t.get(b) == 1)
                                ++ curr_vertex_neighbor_tag_cnt[b];
                    }
                }

                // for (int i = 0; i <= FLAGS_p; ++ i) printf("%d, ", curr_vertex_neighbor_tag_cnt[i]); puts("");

                // choose a tag for current vertex
                ull max_tag_cnt = *max_element(curr_vertex_neighbor_tag_cnt + 1, curr_vertex_neighbor_tag_cnt + 1 + FLAGS_p);
                // LOG(max_tag_cnt);

                ull candidate_tag_cnt = 0;
                int candidate_tag = 0;
                for (int b = 1; b <= FLAGS_p; ++ b)
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

                if (candidate_tag > FLAGS_p || candidate_tag < 1)   
                {
                    LOG(top);
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
                for (int i = h[top]; ~i; i = ne[i])
                {
                    if (st[i])      continue;

                    int to = e[i];
                    if (degree[to] == 1)    
                    {
                        st[i] = 1;  
                        continue;
                    }

                    all_neighbor_covered = false;
                    ++ neighbor_uncovered_cnt;
                    
                    if (vertex2tag[to].get(candidate_tag))
                    {
                        st[i] = 1;
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

            // LOG(top);
            // LOG(vertex2tag[top]);
            
            for (int i = h[top]; i != - 1; i = ne[i])
            {
                int to = e[i];
                if (vis[to]) continue;
                q.push(to);
            }
        }
    }

    delete[] curr_vertex_neighbor_tag_cnt;
}

void degree_1_vertex_assignment()
{
    size_t* curr_vertex_neighbor_tag_cnt = new size_t[FLAGS_p + 1];
    
    for (auto v : vertex_set)
    {
        if (degree[v] > 1)  continue;

        memset(curr_vertex_neighbor_tag_cnt, 0, sizeof(size_t) * (FLAGS_p + 1));

        // choose a tag for current vertex

        int to = e[h[v]];
        const auto & v2t = vertex2tag[to];
        int candidate_tag = 0;

        if (!v2t.empty())
        {
            for (int b = 1; b <= FLAGS_p; ++ b)
                if (v2t.get(b) == 1)
                    if (global_tag_distribute[b] < global_tag_distribute[candidate_tag])
                        candidate_tag = b;
        }

        if (candidate_tag > FLAGS_p || candidate_tag < 1)   
        {
            LOG(v);
        }
        else
        {
            if (vertex2tag[v].get(candidate_tag) != 1)
            {
                ++ global_tag_distribute[candidate_tag];
                vertex2tag[v].set_bit(candidate_tag);
            }
        }
        
    }

    delete[] curr_vertex_neighbor_tag_cnt;
}

void union_tag()
{
	priority_queue<array<int, 2>, vector<array<int, 2>>, greater<array<int, 2>>> pq;

	for (int i = 1; i <= FLAGS_p / 2; i++)
		pq.push({0, i});
	vector<int> block2goal(FLAGS_p + 1, 0);


	sort(global_tag_distribute + 1, global_tag_distribute + 1 + FLAGS_p);
	for (int i = FLAGS_p; i >= 1; -- i)
	{
		const auto [v_cnt, p_id] = pq.top();
        // cout << v_cnt + (int)global_tag_distribute[i] << ' ' << p_id << ' ' << i << endl;
		pq.pop();
		block2goal[i] = p_id;

		//tmp = make_pair(block[i].first, 加入序号)
		pq.push({v_cnt + (int)global_tag_distribute[i], p_id});
	}

    tag2size = vector<ull>(FLAGS_p / 2 + 1, 0);
	while (!pq.empty())
	{
        tag2size[pq.top()[1]] = pq.top()[0];
		pq.pop();
	}
    // LOG(block2goal);
    LOG(tag2size);

    // map<int, int> cnt;
    // int sum = 0;
    // for (int i = 1; i <= FLAGS_p; ++ i)
    // {
    //     cnt[block2goal[i]]++;
    //     sum ++;
    // }
    // LOG(cnt);
    // LOG(sum);
}

// ./main -p 8 -filename ../dataset/mini.txt
// ./main -p 8 -filename ../dataset/com-amazon.ungraph.txt
// ./main -p 16 -filename ../dataset/com-lj.ungraph.txt
// ./main -p 16 -filename ../dataset/com-orkut.ungraph.txt

int main(int argc, char *argv[])
{
    google::ParseCommandLineNonHelpFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    google::HandleCommandLineHelpFlags();

    // LOG(FLAGS_p);
    // LOG(FLAGS_filename);

    global_tag_distribute = new ull[FLAGS_p + 1];
    memset(global_tag_distribute, 0, sizeof (ull) * (FLAGS_p + 1));

    for (int i = 0; i < N; ++ i)   
        vertex2tag[i] = dense_bitset(FLAGS_p + 1);

    read_graph();

    Timer timer;
    timer.start();

    Partitioner *partitioner = NULL;
    partitioner = new TagPartitioner(FLAGS_filename);
    partitioner->split();

    int random_cnt = FLAGS_p * 5;
    random_tag(random_cnt);
    bfs_walk(random_cnt);
    degree_1_vertex_assignment();
    // union_tag();

    int all2 = 0;
    for (int i = 1; i <= FLAGS_p; ++ i)
        cout << i << " : " << global_tag_distribute[i] << endl, all2 += global_tag_distribute[i];
    // for (int i = 1; i <= FLAGS_p / 2; ++ i)
    //     cout << i << " : " << tag2size[i] << endl, all2 += tag2size[i];
    cout << "VERTEX_CNT : " << vertex_set.size() << endl;
    cout << "ALL_TAG_CNT : " << all2 << endl;

    LOG(degree_1_vertex_cnt);
    delete[] global_tag_distribute;
    delete[] bfs_src;

    timer.stop();
    LOG(timer.get_time());
    
    return 0;
}