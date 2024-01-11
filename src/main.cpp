#include <memory>

#include "fsm_partitioner.hpp"
#include "ne_partitioner.hpp"
#include "hep_partitioner.hpp"
#include "ebv_partitioner.hpp"
#include "dbh_partitioner.hpp"
#include "hdrf_partitioner.hpp"
#include "fennel_partitioner.hpp"
#include "bpart_partitioner.hpp"
#include "hybridbl_partitioner.hpp"
#include "hybrid_partitioner.hpp"
// #include "edgelist2adjlist.hpp"
#include "vertex2edgepart.hpp"
// #include "test.hpp"

DECLARE_bool(help);
DECLARE_bool(helpshort);

DEFINE_int32(p, 32, "number of parititions");
DEFINE_string(filename, "", "the file name of the input graph");
DEFINE_string(filetype, "edgelist",
              "the type of input file (supports 'edgelist' and 'adjlist')");
DEFINE_bool(write, false, "write out partition result");
DEFINE_int32(k, 1, "split factor, i.e. the exact partitions count / remain partitions");
DEFINE_bool(fastmerge, false, "use fast merge?");
DEFINE_string(method, "sne",
              "partition method: ne, sne, random, and dbh");

DEFINE_bool(write_low_degree_edgelist, false, "Should the list of edges incident to a low-degree vertex be written out to a file?");
DEFINE_double(hdf, 100, "High-degree factor: hdf * average_degree = high-degree threshold (hdth). Called \\tau in the paper. Vertices with than hdth neighbors are treated specially in fast NE");
DEFINE_double(lambda, 1.1, "Lambda value to weigh in balancing score in streaming partitioning via HDRF");
DEFINE_bool(extended_metrics, false, "Display extended metrics in the result");
DEFINE_bool(random_streaming, false, "Use random streaming instead of HDRF in the second phase of HEP.");
DEFINE_bool(hybrid_NE, false, "Perform hybrid partitioning in HEP-style, but use NE instead of NE++ for the first phase.");

int main(int argc, char *argv[])
{
    std::string usage = "-filename <path to the input graph> "
                        "[-filetype <edgelist|adjlist>] "
                        "[-p <number of partitions>] "
                        "[-memsize <memory budget in MB>]";
    google::SetUsageMessage(usage);
    google::ParseCommandLineNonHelpFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = 1; // output log to stderr
    if (FLAGS_help) {
        FLAGS_help = false;
        FLAGS_helpshort = true;
    }
    google::HandleCommandLineHelpFlags();

    Timer timer;
    timer.start();

    std::unique_ptr<PartitionerBase> partitioner = nullptr;
    std::string method = FLAGS_method;
    if (method.substr(0, 3) == "fsm") {
        if (FLAGS_k > 1) {
            partitioner = std::make_unique<FsmPartitioner>(FLAGS_filename);
        } else {
            method = method.substr(4);
        }
    }
    
    if (method == "ne")
        partitioner = std::make_unique<NePartitioner<adj_t>>(FLAGS_filename, false);
    else if (method == "hep")
        partitioner = std::make_unique<HepPartitioner<adj_t>>(FLAGS_filename, false);
    else if (method == "ebv")
        partitioner = std::make_unique<EbvPartitioner>(FLAGS_filename, false);
    else if (method == "dbh")
        partitioner = std::make_unique<DbhPartitioner>(FLAGS_filename, false);
    else if (method == "hdrf")
        partitioner = std::make_unique<HdrfPartitioner>(FLAGS_filename, false);
    else if (method == "fennel")
        partitioner = std::make_unique<FennelPartitioner<adj_t>>(FLAGS_filename, false);
    else if (method == "bpart")
        partitioner = std::make_unique<BPartPartitioner<adj_t>>(FLAGS_filename, false);
    else if (method == "hybridbl")
        partitioner = std::make_unique<HybridBLPartitioner<adj_t>>(FLAGS_filename, false);
    else if (method == "hybrid")
        partitioner = std::make_unique<HybridPartitioner>(FLAGS_filename, false);
    // else if (method == "e2a")
    //     partitioner = std::make_unique<Edgelist2Adjlist>(FLAGS_filename);
    else if (method.substr(0, 3) == "v2e")
        partitioner = std::make_unique<Vertex2EdgePart>(FLAGS_filename, false);
    // else if (method == "test")
    //     partitioner = std::make_unique<Test>(FLAGS_filename);
    // LOG(INFO) << "partition method: " << method;
    partitioner->split();

    timer.stop();
    LOG(INFO) << "total time: " << timer.get_time();
}
