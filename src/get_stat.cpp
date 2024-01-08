#include <algorithm>
#include <numeric>
#include <cmath>

#include "common.hpp"
#include "dense_bitset.hpp"
using namespace std;
using ull = unsigned long long;

// Removes \n from the end of line
void FIXLINE(char *s)
{
    int len = (int)strlen(s) - 1;
    if (s[len] == '\n')
        s[len] = 0;
}


int main(int argc, char** argv)
{
    string path = string(argv[1]);
    int p = atoi(argv[2]);
    // int num_vertices = atoi(argv[3]);
    dense_bitset v_set;
    v_set.resize(500'000'000);
    vector<dense_bitset> is_boundarys(p);
    is_boundarys.assign(p, dense_bitset(500'000'000));
    vector<size_t> occupied(p);


    FILE *inf = fopen(path.c_str(), "r");
    size_t bytesread = 0;
    size_t linenum = 0;
    if (inf == NULL) {
        cerr << "Could not load:" << path
                   << ", error: " << strerror(errno) << std::endl;
        exit(1);
    }
    ull num_edges = 0;
    cout << "Reading in edge list format!" << std::endl;
    char s[1024];
    while (fgets(s, 1024, inf) != NULL) {
        linenum++;
        if (linenum % 10000000 == 0) {
            cout << "Read " << linenum << " lines, "
                      << bytesread / 1024 / 1024. << " MB" << std::endl;
        }
        FIXLINE(s);
        bytesread += strlen(s);
        if (s[0] == '#')
            continue; // Comment
        if (s[0] == '%')
            continue; // Comment

        char delims[] = "\t, ";
        char *t;
        t = strtok(s, delims);
        if (t == NULL) {
            cerr << "Input file is not in right format. "
                       << "Expecting \"<from>\t<to>\". "
                       << "Current line: \"" << s << "\"\n";
            exit(1);
        }
        ull from = atoi(t);
        t = strtok(NULL, delims);
        if (t == NULL) {
            cerr << "Input file is not in right format. "
                       << "Expecting \"<from>\t<to>\". "
                       << "Current line: \"" << s << "\"\n";
            exit(1);
        }
        ull to = atoi(t);

        t = strtok(NULL, delims);
        if (t == NULL) {
            cerr << "Input file is not in right format. "
                       << "Expecting \"<from>\t<to>\". "
                       << "Current line: \"" << s << "\"\n";
            exit(1);
        }
        ull bid = atoi(t);

        is_boundarys[bid].set_bit_unsync(from);
        is_boundarys[bid].set_bit_unsync(to);
        ++occupied[bid];
        ++num_edges;
        v_set.set_bit_unsync(from);
        v_set.set_bit_unsync(to);
    }
    fclose(inf);

    vid_t num_vertices = v_set.popcount();

    std::cerr << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
    std::vector<size_t> bucket2vcnt(p, 0);
    for (int bid = 0; bid < p; ++bid) {
        bucket2vcnt[bid] = is_boundarys[bid].popcount();
    }
    size_t max_part_vertice_cnt = *std::max_element(bucket2vcnt.begin(), bucket2vcnt.end());
    size_t all_part_vertice_cnt = std::accumulate(bucket2vcnt.begin(), bucket2vcnt.end(), (size_t)0);
    size_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
    size_t all_part_edge_cnt = std::accumulate(occupied.begin(), occupied.end(), (size_t)0);

    for (int b = 0; b < p; ++ b) 
        cout << "Bucket_info: " << b
                << ", vertices: " << bucket2vcnt[b]
                << ", edges: " << occupied[b] << endl;
    
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
    
    cout << std::string(20, '#') << "\tVertice    balance\t" << std::string(20, '#') << endl;
    cout << "Max vertice count / avg vertice count: "
              << (double)max_part_vertice_cnt / ((double)num_vertices / (p)) << endl;
    cout << "Max Vertice count: "
              << max_part_vertice_cnt << endl;
    cout << "Avg Vertice count(No replicate): "
              << num_vertices / p << endl;
    cout << "Vertice std_vertice_deviation / avg: "
              << std_vertice_deviation / avg_vertice_cnt << endl;

    cout << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#') << endl;
    cout << "Max edge count / avg edge count: "
              << (double)max_part_edge_cnt / avg_edge_cnt << endl;
    cout << "Max Edge count: "
              << max_part_edge_cnt << endl;
    cout << "Avg Edge count: "
              << num_edges / p << endl;
    cout << "Edge std_edge_deviation / avg: "
              << std_edge_deviation / avg_edge_cnt << endl;

    cout << std::string(20, '#') << "\tReplicate    factor\t" << std::string(20, '#') << endl;
    cout << "replication factor (final): " << (double)all_part_vertice_cnt / num_vertices << endl;
}
