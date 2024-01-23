#ifndef PARTITIONER_HPP
#define PARTITIONER_HPP

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "dense_bitset.hpp"
#include "hep_graph.hpp"
#include "util.hpp"

class PartitionerBase
{

public:
    Timer total_time;
    Timer partition_time;
    vid_t num_vertices;
    eid_t num_edges;
    bid_t num_partitions;
    virtual void split() = 0;        
    virtual void calculate_stats(bool called_by_fsm) = 0;        
};

class EdgePartitioner: public PartitionerBase
{
public:
    std::vector<eid_t> occupied;
    std::vector<dense_bitset> is_boundarys;
    std::vector<edge_t> edges;
    std::vector<vid_t> degrees;
    std::vector<bid_t> edgelist2bucket;

    void calculate_stats(bool called_by_fsm = false)
    {
        LOG(INFO) << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';
        if (called_by_fsm) {
            num_partitions = FLAGS_p;
        }

        std::vector<vid_t> num_bucket_vertices(num_partitions, 0);
        for (bid_t b = 0; b < num_partitions; ++b) {
            num_bucket_vertices[b] = is_boundarys[b].popcount();
        }
        vid_t max_part_vertice_cnt = *std::max_element(num_bucket_vertices.begin(), num_bucket_vertices.end());
        vid_t all_part_vertice_cnt = accumulate(num_bucket_vertices.begin(), num_bucket_vertices.end(), (vid_t)0);
        eid_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
        eid_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (eid_t)0);

        for (bid_t b = 0; b < num_partitions; ++b) {
            LOG(INFO) << "Bucket_info: " << b
                    << ", vertices: " << num_bucket_vertices[b]
                    << ", edges: " << occupied[b];
        }
        
        double avg_vertice_cnt = static_cast<double>(all_part_vertice_cnt) / num_partitions;
        double avg_edge_cnt = static_cast<double>(all_part_edge_cnt) / num_partitions;

        double std_vertice_deviation = 0.0;
        double std_edge_deviation = 0.0;
        for (bid_t b = 0; b < num_partitions; ++b) {
            std_vertice_deviation += pow(num_bucket_vertices[b] - avg_vertice_cnt, 2);
            std_edge_deviation += pow(occupied[b] - avg_edge_cnt, 2);
        }
        std_vertice_deviation = sqrt(static_cast<double>(std_vertice_deviation) / num_partitions);
        std_edge_deviation = sqrt(static_cast<double>(std_edge_deviation) / num_partitions);

        LOG(INFO) << std::string(20, '#') << "\tVertice    balance\t" << std::string(20, '#');
        LOG(INFO) << "Max vertice count / avg vertice count: "
                << (double)max_part_vertice_cnt / ((double)num_vertices / (num_partitions));
        LOG(INFO) << "Max Vertice count: "
                << max_part_vertice_cnt;
        LOG(INFO) << "Avg Vertice count(No replicate): "
                << num_vertices / num_partitions;
        LOG(INFO) << "Vertice std_vertice_deviation / avg: "
                << std_vertice_deviation / avg_vertice_cnt;

        LOG(INFO) << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#');
        LOG(INFO) << "Max edge count / avg edge count: "
                << (double)max_part_edge_cnt / avg_edge_cnt;
        LOG(INFO) << "Max Edge count: "
                << max_part_edge_cnt;
        LOG(INFO) << "Avg Edge count: "
                << num_edges / num_partitions;
        LOG(INFO) << "Edge std_edge_deviation / avg: "
                << std_edge_deviation / avg_edge_cnt;

        CHECK_EQ(all_part_edge_cnt, num_edges);
        LOG(INFO) << std::string(20, '#') << "\tReplicate    factor\t" << std::string(20, '#');
        LOG(INFO) << "replication factor (final): " << (double)all_part_vertice_cnt / num_vertices;
    }
};

class VertexPartitioner: public PartitionerBase
{
public:
    std::vector<eid_t> occupied;
    std::vector<dense_bitset> is_boundarys;
    std::vector<edge_t> edges;
    std::vector<vid_t> degrees;
    std::vector<bid_t> vertex2bucket;

    void calculate_stats(bool called_by_fsm = false)
    {
        LOG(INFO) << std::string(25, '#') << " Calculating Statistics " << std::string(25, '#') << '\n';

        num_partitions = FLAGS_p;

        std::vector<vid_t> num_bucket_vertices(num_partitions, 0);
        eid_t total_cut_edges = 0;
        for (bid_t b = 0; b < num_partitions; ++b) {
            num_bucket_vertices[b] = is_boundarys[b].popcount();
            total_cut_edges += occupied[b];
        }
        vid_t max_part_vertice_cnt = *std::max_element(num_bucket_vertices.begin(), num_bucket_vertices.end());
        vid_t all_part_vertice_cnt = accumulate(num_bucket_vertices.begin(), num_bucket_vertices.end(), (vid_t)0);
        eid_t max_part_edge_cnt = *std::max_element(occupied.begin(), occupied.end());
        eid_t all_part_edge_cnt = accumulate(occupied.begin(), occupied.end(), (eid_t)0);

        for (bid_t b = 0; b < num_partitions; ++b) 
            LOG(INFO) << "Bucket_info: " << b
                    << ", vertices: " << num_bucket_vertices[b]
                    << ", edges: " << occupied[b];
        
        double avg_vertice_cnt = (double)all_part_vertice_cnt / (num_partitions);
        double avg_edge_cnt = (double)all_part_edge_cnt / (num_partitions);
        double std_vertice_deviation = 0.0;
        double std_edge_deviation = 0.0;
        for (int b = 0; b < num_partitions; ++ b) {
            std_vertice_deviation += pow(num_bucket_vertices[b] - avg_vertice_cnt, 2);
            std_edge_deviation += pow(occupied[b] - avg_edge_cnt, 2);
        }
        std_vertice_deviation = sqrt((double)std_vertice_deviation / num_partitions);
        std_edge_deviation = sqrt((double)std_edge_deviation / num_partitions);
        
        LOG(INFO) << std::string(20, '#') << "\tVertice    balance\t" << std::string(20, '#');
        LOG(INFO) << "Max vertice count / avg vertice count: "
                << (double)max_part_vertice_cnt / ((double)num_vertices / (num_partitions));
        LOG(INFO) << "Max Vertice count: "
                << max_part_vertice_cnt;
        LOG(INFO) << "Avg Vertice count(No replicate): "
                << num_vertices / num_partitions;
        LOG(INFO) << "Vertice std_vertice_deviation / avg: "
                << std_vertice_deviation / avg_vertice_cnt;
        LOG(INFO) << "Vertice jains_fairness: "
                << jains_fairness<vid_t>(num_bucket_vertices);

        LOG(INFO) << std::string(20, '#') << "\tEdge       balance\t" << std::string(20, '#');        
        LOG(INFO) << "Max edge count / avg edge count: "
                << (double)max_part_edge_cnt / avg_edge_cnt;
        LOG(INFO) << "Max Edge count: "
                << max_part_edge_cnt;
        LOG(INFO) << "Avg Edge count: "
                << num_edges / num_partitions;
        LOG(INFO) << "Edge std_edge_deviation / avg: "
                << std_edge_deviation / avg_edge_cnt;
        LOG(INFO) << "Edge jains_fairness: "
                << jains_fairness<eid_t>(occupied);

        CHECK_EQ(all_part_vertice_cnt, num_vertices);

        LOG(INFO) << std::string(20, '#') << "\tEdge  cut  ratio\t" << std::string(20, '#');
        total_cut_edges -= num_edges;
        LOG(INFO) << "Edge cut ratio (final): " << (double)total_cut_edges / num_edges;
    }


        
};


template <typename TAdj>
class AdjListEPartitioner: public EdgePartitioner
{
public:
    mem_graph_t<TAdj> mem_graph;
};

template class AdjListEPartitioner<adj_with_bid_t>;
template class AdjListEPartitioner<adj_t>;

class EdgeListEPartitioner: public EdgePartitioner
{
};



template <typename TAdj>
class AdjListVPartitioner: public VertexPartitioner
{
public:
    mem_graph_t<TAdj> mem_graph;
};

template class AdjListVPartitioner<adj_with_bid_t>;
template class AdjListVPartitioner<adj_t>;

#endif