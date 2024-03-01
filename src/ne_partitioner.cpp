#include "ne_partitioner.hpp"
#include "conversions.hpp"

template <typename TAdj>
NePartitioner<TAdj>::NePartitioner(std::string basefilename, bool need_k_split)
    : basefilename(basefilename), rd(), gen(rd())
    , need_k_split(need_k_split), writer(basefilename, !need_k_split && FLAGS_write)
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

    average_degree = num_edges * 2.0 / num_vertices;
    assigned_edges = 0;
    capacity = static_cast<double>(num_edges) * BALANCE_RATIO / num_partitions + 1;
    occupied.assign(num_partitions, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundarys.assign(num_partitions, dense_bitset(num_vertices));
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    if (need_k_split) {
        edgelist2bucket.assign(num_edges, kInvalidBid);
    }

    Timer read_timer;
    read_timer.start();
    LOG(INFO) << "loading...";
    LOG(INFO) << sizeof(edge_t) * num_edges / 1024.0 / 1024 / 1024 << " G bytes needed for edges " << __FILE__<<":"<<__LINE__;
    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    LOG(INFO) << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    LOG(INFO) << "time used for graph input and construction: " << read_timer.get_time();

    size_t num_out_edges = 0;
    for (vid_t vid = 0; vid < num_vertices; ++vid) {
        num_out_edges += adj_out[vid].size();
    }
    LOG(INFO) << "num_out_edges: " << num_out_edges;
    LOG(INFO) << "num_edges: " << num_edges;
    CHECK_EQ(num_out_edges, num_edges);
};

template <typename TAdj>
void NePartitioner<TAdj>::assign_remaining()
{
    auto &is_boundary = is_boundarys[num_partitions - 1], &is_core = is_cores[num_partitions - 1];
    for (vid_t u = 0; u < num_vertices; ++u) {
        for (auto &i : adj_out[u]) {
            if (edges[i.v].valid()) {
                assign_edge(num_partitions - 1, u, edges[i.v].second, i.v);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }
        }
    }

    for (vid_t i = 0; i < num_vertices; ++i) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            for (bid_t b = 0; b < num_partitions - 1; ++b) {
                if (is_cores[b].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
            }
        }
    }
}

template <typename TAdj>
void NePartitioner<TAdj>::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << (uint32_t)num_partitions;

    min_heap.reserve(num_vertices);

    LOG(INFO) << "partitioning...";
    partition_time.start();
    for (bucket = 0; bucket < num_partitions - 1; ++bucket) {
        std::cerr << (uint32_t)bucket << ", ";
        DLOG(INFO) << "sample size: " << adj_out.num_edges();
        while (occupied[bucket] < capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) {
                if (need_k_split) {
                    if (!get_free_vertex_by_rand(vid)) {
                        DLOG(INFO) << "partition " << bucket
                                << " stop: no free vertices";
                        break;
                    }
                } else {
                    if (!get_free_vertex(vid)) {
                        DLOG(INFO) << "partition " << bucket
                                << " stop: no free vertices";
                        break;
                    }
                }
                
                d = adj_out[vid].size() + adj_in[vid].size();
            } else {
                min_heap.remove(vid);
                /* CHECK_EQ(d, adj_out[vid].size() + adj_in[vid].size()); */
            }

            occupy_vertex(vid, d);
        }
        min_heap.clear();
        for (int direction = 0; direction < 2; ++direction) {
            for (vid_t vid = 0; vid < num_vertices; ++vid) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (vid_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        ++i;
                    } else {
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
        }
    }
    bucket = num_partitions - 1;
    std::cerr << bucket << std::endl;
    assign_remaining();

    partition_time.stop();
    total_time.stop();

    if (!need_k_split) {
        LOG(INFO) << "partitioning time: " << partition_time.get_time();
        LOG(INFO) << "total time: " << total_time.get_time();
    }

    CHECK_EQ(assigned_edges, num_edges);

    calculate_stats();
}
