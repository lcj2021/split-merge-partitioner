#pragma once

#include "util.hpp"
#include "dense_bitset.hpp"
#include "graph.hpp"
#include "partitioner.hpp"
#include "conversions.hpp"

class Edgelist2Adjlist : public Partitioner
{
  private:
    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges;
    int p;

    std::vector<vid_t> degrees;

    std::vector<edge_t> edges;
    graph_t adj_out, adj_in;
    std::vector<dense_bitset> is_boundarys;
    std::ofstream fout;

  public:
    Edgelist2Adjlist(std::string basefilename)
        : basefilename(basefilename), fout(basefilename + ".adjlist")
    {
        CHECK_EQ(FLAGS_filetype, "edgelist");

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
        adj_out.resize(num_vertices);
        adj_in.resize(num_vertices);
        is_boundarys.assign(p, dense_bitset(num_vertices));

        Timer read_timer;
        read_timer.start();
        LOG(INFO) << "loading...";
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
    }

    void split();
};

void 
Edgelist2Adjlist::split()
{
    fout << num_vertices << ' ' << num_edges << " 010 " << 1 << std::endl;
    for (vid_t vid = 0; vid < num_vertices; ++ vid) {
        fout << degrees[vid] << ' ';
        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size(); ++ i) {
                vid_t &uid = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                fout << uid + 1 << ' ';
            }
        }
        fout << std::endl;
    }
}