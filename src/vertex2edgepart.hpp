#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <parallel/algorithm>
#include <boost/unordered_map.hpp>

#include "util.hpp"
#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"

class Vertex2EdgePart : public Partitioner {
private:
	std::string basefilename;

    // vid_t num_vertices;
    // size_t num_edges;
    int k, p;

    int fin;
    off_t filesize;

    std::vector<int16_t> vertex2partition; // maps vertex id to partition id
    vid_t assigned_vertices;
    // std::vector<size_t> occupied;
	// std::vector<dense_bitset> is_boundarys;

    boost::unordered_map<std::tuple<vid_t, vid_t>, size_t> uv2edgeid;
    size_t curr_edge_cnt;

    struct bucket_info_item {
        dense_bitset is_mirror;
        size_t occupied, old_id, replicas;
        bool is_chosen;
        bucket_info_item(vid_t num_edges) {
            is_mirror = dense_bitset(num_edges);
            old_id = occupied = replicas = 0;
            is_chosen = false;
        }
        bool operator < (const bucket_info_item& rhs) const {
            if (is_chosen != rhs.is_chosen) return is_chosen > rhs.is_chosen;
            return old_id < rhs.old_id;
        }
    };
    std::vector<bucket_info_item> bucket_info;

    // Removes \n from the end of line
    void FIXLINE(char *s)
    {
        int len = (int)strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
    }

    void initDataStructures()
    {
    	vertex2partition.resize(num_vertices + 1, -1);
    	is_boundarys.resize(p * k, dense_bitset(num_vertices + 1));
        occupied.assign(p * k, 0);

        bucket_info.assign(p * k, bucket_info_item(num_edges));
        for (int i = 0; i < p * k; ++ i) bucket_info[i].old_id = i;
    }

    size_t rearrange_vertice(std::vector<edge_t> &e, const std::unordered_map<int, int> &valid_bucket)
    {
        size_t curr_assigned_vertices = 0;
        for (vid_t v_id = 1; v_id <= num_vertices; ++ v_id) {
            int16_t &v_bucket = vertex2partition[v_id];
            if (valid_bucket.count(v_bucket)) {
                v_bucket = valid_bucket.at(v_bucket);
                ++ curr_assigned_vertices;
            } else {        // [[unlikely]] No edge should be left
                LOG(ERROR) << "Vertice(id): " << v_id << ", bucket: " << v_bucket << ", should not be left in the final round\n";
            }
        }
        return curr_assigned_vertices;
    }

    std::unordered_map<int, int> mp;
    int curr_bucket_id;
    int get_final_bucket(int bucket_id)
    {
        if (mp.count(bucket_id)) {
            return mp[bucket_id];
        } else {
            return mp[bucket_id] = curr_bucket_id ++;
        }
    }

    void merge();
    void calculate_stats();

    int merge_bucket(int dst, int src, bool &has_intersection);
    std::unordered_map<int, int> merge_by_size();
    std::unordered_map<int, int> merge_by_overlap();

public:
	Vertex2EdgePart(std::string basefilename);
	virtual ~Vertex2EdgePart();

	void split();

	void readVertexPartitioning();

	int findEdgePartition(vid_t u, vid_t v);
};
