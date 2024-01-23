#include <numeric>

#include "dbh_partitioner.hpp"
#include "conversions.hpp"

DbhPartitioner::DbhPartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), writer(basefilename, !need_k_split && FLAGS_write)
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
    CHECK_EQ(sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t), filesize);

    num_partitions = FLAGS_p;
    if (need_k_split) {
        num_partitions *= FLAGS_k;
    }

    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

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
    // std::shuffle(edges.begin(), edges.end(), rd);
    bid_t bucket;
    for (eid_t eid = 0; eid < edges.size(); ++eid) {
        edge_t e = edges[eid];
        vid_t u = e.first, v = e.second;
        vid_t w = degrees[u] <= degrees[v] ? u : v;
        bucket = w % num_partitions;
        assign_edge(bucket, u, v, eid);
        if (eid % 50000000 == 0) {
            LOG(INFO) << "Processing edges " << eid;
        }
    }

    // eid_t edge_id = 0;
    // while (fin_ptr < fin_end) {
    //     edge_t *e = (edge_t *)fin_ptr;
    //     fin_ptr += sizeof(edge_t);
    //     if (edge_id % 50000000 == 0) {
    //         LOG(INFO) << "Processing edges " << edge_id;
    //     }
    //     ++edge_id;
    //     vid_t u = e->first, v = e->second;
    //     vid_t w = degrees[u] <= degrees[v] ? u : v;
    //     bid_t bucket = w % num_partitions;
    //     assign_edge(bucket, u, v, edge_id);
    // }

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
    for (bid_t b = 0; b < num_partitions; ++b) {
        LOG(INFO) << b << ' ' << is_boundarys[b].popcount() << ' ' << occupied[b];
    }
    calculate_stats();
}
