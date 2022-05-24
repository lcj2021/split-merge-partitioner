#include <bits/stdc++.h>
#include "dense_bitset.hpp"
#include "util.hpp"
#include "tag_partitioner.hpp"
using namespace std;
#define endl "\n"
#define all(x) x.begin(), x.end()
typedef unsigned long long ull;
// const ull N = 1e7 + 10, M = 3e8 + 10;

DEFINE_string(filetype, "edgelist",
              "the type of input file (supports 'edgelist' and 'adjlist')");
DEFINE_int32(p, 10, "number of parititions");
DEFINE_string(filename, "", "the file name of the input graph");

// ./main -p 8 -filename ../dataset/mini.txt
// ./main -p 8 -filename ../dataset/com-amazon.ungraph.txt
// ./main -p 16 -filename ../dataset/com-lj.ungraph.txt
// ./main -p 16 -filename ../dataset/com-orkut.ungraph.txt

int main(int argc, char *argv[])
{
    google::ParseCommandLineNonHelpFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    google::HandleCommandLineHelpFlags();
    FLAGS_logtostderr = 1; // output log to stderr

    Timer timer;
    timer.start();

    Partitioner *partitioner = NULL;
    partitioner = new TagPartitioner(FLAGS_filename);
    partitioner->split();

    timer.stop();
    LOG(INFO) << "Total time used " << timer.get_time();
    
    return 0;
}