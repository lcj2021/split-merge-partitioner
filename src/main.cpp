#include <string>

#include "util.hpp"
#include "ne_partitioner.hpp"
#include "fsm_partitioner.hpp"
#include "ebv_partitioner.hpp"
#include "dbh_partitioner.hpp"
#include "hdrf_partitioner.hpp"
#include "hep_partitioner.hpp"
#include "edgelist2adjlist.hpp"
#include "vertex2edgepart.hpp"
#include "test.hpp"

DECLARE_bool(help);
DECLARE_bool(helpshort);

DEFINE_int32(p, 32, "number of parititions");
DEFINE_string(filename, "", "the file name of the input graph");
DEFINE_string(filetype, "edgelist",
              "the type of input file (supports 'edgelist' and 'adjlist')");
DEFINE_bool(write, false, "write out partition result");
DEFINE_int32(k, 2, "the exact partitions count / remain partitions");
DEFINE_bool(fastmerge, false, "whether use fast merge");
DEFINE_string(method, "sne",
              "partition method: ne, sne, random, and dbh");

DEFINE_bool(write_low_degree_edgelist, false, "Should the list of edges incident to a low-degree vertex be written out to a file?");
DEFINE_double(hdf, 100, "High-degree factor: hdf * average_degree = high-degree threshold (hdth). Called \tau in the paper. Vertices with than hdth neighbors are treated specially in fast NE");
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

    Partitioner *partitioner = NULL;
    if (FLAGS_method == "ne")
        partitioner = new NePartitioner(FLAGS_filename, false);
    else if (FLAGS_method == "dbh")
        partitioner = new DbhPartitioner(FLAGS_filename, false);
    else if (FLAGS_method == "hdrf")
        partitioner = new HdrfPartitioner(FLAGS_filename, false);
    else if (FLAGS_method == "ebv")
        partitioner = new EbvPartitioner(FLAGS_filename, false);
    else if (FLAGS_method == "hep")
        partitioner = new HepPartitioner(FLAGS_filename, false);
    else if (FLAGS_method.substr(0, 3) == "fsm")
        partitioner = new FsmPartitioner(FLAGS_filename);
    else if (FLAGS_method == "e2a")
        partitioner = new Edgelist2Adjlist(FLAGS_filename);
    else if (FLAGS_method == "v2e")
        partitioner = new Vertex2EdgePart(FLAGS_filename);
    else if (FLAGS_method == "test")
        partitioner = new Test(FLAGS_filename);
    LOG(INFO) << "partition method: " << FLAGS_method;
    partitioner->split();

    timer.stop();
    LOG(INFO) << "total time: " << timer.get_time();
}
