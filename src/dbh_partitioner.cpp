#include <numeric>

#include "dbh_partitioner.hpp"
#include "graph.hpp"
#include "conversions.hpp"

DbhPartitioner::DbhPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename)
{
    if (need_k_split || FLAGS_write == "none") {
        writer = std::make_unique<EdgepartWriterBase<vid_t, bid_t>>(basefilename);
    } else {
        if (FLAGS_write == "onefile") {
            writer = std::make_unique<EdgepartWriterOnefile<vid_t, bid_t>>(basefilename);
        } else if (FLAGS_write == "multifile") {
            writer = std::make_unique<EdgepartWriterMultifile<vid_t, bid_t>>(basefilename);
        }
    }
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
    CHECK_EQ(sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t), filesize);

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    LOG(INFO) << "constructing...";

    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    occupied.assign(num_partitions, 0);
    avg_edge_cnt = (double)num_edges / FLAGS_p;
    // edgelist2bucket.assign(num_edges, kInvalidBid);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
}

void DbhPartitioner::split()
{
    partition_time.start();
    bid_t bucket;

    std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
    eid_t chunk_size = std::min(num_edges, (eid_t)100000);
    eid_t num_remaining_edges = num_edges; // number of edges to be read from file

    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); 

    stream_edges.resize(chunk_size);

    eid_t num_edges_read = 0; // number of edges read from file
    while (num_remaining_edges > 0) { // edges to be read
        fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
        for (eid_t i = 0; i < chunk_size; ++i, ++num_edges_read) {
            const auto& [u, v] = stream_edges[i];
            vid_t w = degrees[u] <= degrees[v] ? u : v;
            bucket = w % num_partitions;
            assign_edge(bucket, u, v, num_edges_read);
            if (num_edges_read % 50000000 == 0) {
                LOG(INFO) << "Processing edges " << num_edges_read;
            }
        }

        num_remaining_edges -= chunk_size;
        if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
            chunk_size = num_remaining_edges;
        }
    }

    partition_time.stop();
    total_time.stop();
    LOG(INFO) << "partition time: " << partition_time.get_time();
    calculate_stats();
}
