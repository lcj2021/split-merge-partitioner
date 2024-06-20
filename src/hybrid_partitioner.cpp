#include "hybrid_partitioner.hpp"
#include "conversions.hpp"

HybridPartitioner::HybridPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd())
    , need_k_split(need_k_split)
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

    average_degree = num_edges * 2.0 / num_vertices;

    LOG(INFO) << "average_degree: " << average_degree;
    assigned_edges = 0;
    occupied.assign(num_partitions, 0);
    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    // edgelist2bucket.assign(num_edges, kInvalidBid);

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();
};

void HybridPartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << num_partitions;

    LOG(INFO) << "partitioning...";
    partition_time.start();

    std::vector<edge_t> stream_edges; // temporary buffer to read edges from file
    eid_t chunk_size = std::min(num_edges, (eid_t)100000);
    eid_t num_remaining_edges = num_edges; // number of edges to be read from file

    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg); 

    stream_edges.resize(chunk_size);

    eid_t eid = 0; // number of edges read from file
    while (num_remaining_edges > 0) { // edges to be read
        fin.read((char *)&stream_edges[0], sizeof(edge_t) * chunk_size);
        for (eid_t i = 0; i < chunk_size; ++i, ++eid) {
            if (eid % 50000000 == 0) {
                LOG(INFO) << "Processing edges " << eid;
            }
            const auto& [uid, vid] = stream_edges[i];
            if (degrees[vid] < degree_threshold) {
                assign_edge(vid % num_partitions, uid, vid, eid);
            } else {
                assign_edge(uid % num_partitions, uid, vid, eid);
            }
        }

        num_remaining_edges -= chunk_size;
        if (num_remaining_edges < chunk_size) { // adapt chunk size for last batch read
            chunk_size = num_remaining_edges;
        }
    }

    partition_time.stop();

    LOG(INFO) << "partition time: " << partition_time.get_time();

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    calculate_stats();
}
