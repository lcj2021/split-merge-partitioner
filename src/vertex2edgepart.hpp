#ifndef VERTEX2EDGEPART_HPP
#define VERTEX2EDGEPART_HPP

#include <random>
#include <boost/unordered_map.hpp>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"

class Vertex2EdgePart : public EdgePartitioner
{
private:
	std::string basefilename;
    std::random_device rd;
    std::mt19937 gen;

    bid_t k;

    int fin;
    off_t filesize;

    std::vector<bid_t> vertex2bucket; // maps vertex id to partition id
    vid_t assigned_vertices;

    eid_t curr_edge_cnt;

    edgepart_writer<vid_t, bid_t> writer;

    struct BucketInfo {
        dense_bitset is_mirror;
        size_t occupied, old_id, replicas;
        bool is_chosen{false};
        BucketInfo(eid_t num_edges) {
            is_mirror = dense_bitset(num_edges);
            old_id = occupied = replicas = 0;
        }
        bool operator < (const BucketInfo& rhs) const {
            if (is_chosen != rhs.is_chosen) return is_chosen > rhs.is_chosen;
            return old_id < rhs.old_id;
        }
    };
    std::vector<BucketInfo> bucket_info;

    // Removes \n from the end of line
    void FIXLINE(char *s)
    {
        int len = (int)strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
    }

    void init_datastructures()
    {
    	vertex2bucket.resize(num_vertices + 1, kInvalidBid);
    	is_boundarys.resize(num_partitions * k, dense_bitset(num_vertices + 1));
        occupied.assign(num_partitions * k, 0);

        bucket_info.assign(num_partitions * k, BucketInfo(num_edges));
        for (bid_t b = 0; b < num_partitions * k; ++b) 
            bucket_info[b].old_id = b;
    }

    eid_t rearrange_vertice(std::vector<edge_t> &e, const std::unordered_map<bid_t, bid_t> &valid_bucket)
    {
        eid_t curr_assigned_vertices = 0;
        for (vid_t vid = 1; vid <= num_vertices; ++vid) {
            bid_t &v_bucket = vertex2bucket[vid];
            if (valid_bucket.count(v_bucket)) {
                v_bucket = valid_bucket.at(v_bucket);
                ++curr_assigned_vertices;
            } else {        // [[unlikely]] No edge should be left
                LOG(ERROR) << "Vertice(id): " << vid << ", bucket: " << v_bucket << ", should not be left in the final round\n";
            }
        }
        return curr_assigned_vertices;
    }

    std::unordered_map<bid_t, bid_t> mp;
    bid_t curr_bucket_id;
    bid_t get_final_bucket(bid_t bucket_id)
    {
        if (mp.count(bucket_id)) {
            return mp[bucket_id];
        } else {
            return mp[bucket_id] = curr_bucket_id++;
        }
    }

    void merge();

    size_t merge_bucket(vid_t dst, vid_t src, bool &has_intersection);
    std::unordered_map<bid_t, bid_t> fast_merge();
    std::unordered_map<bid_t, bid_t> precise_merge();

	void read_vertexpart();

	bid_t vertex2edgepartID(vid_t u, vid_t v);

    bool split_in_adjlist()
    {
        // read the metis adjacency list
        curr_edge_cnt = 0;
        FILE *fin = fopen((basefilename + ".adjlist").data(), "r");
        if (fin == NULL) {
            LOG(INFO) << "Could not load:" << basefilename
                        << " error: " << strerror(errno);
            return false;
        }
        LOG(INFO) << "Reading in adjacency list format!";
        LOG(INFO) << basefilename;

        int maxlen = 1000000000;
        char *s = (char *)malloc(maxlen);

        size_t bytesread = 0;

        char delims[] = " \t";
        eid_t linenum = 0;

        while (fgets(s, maxlen, fin) != NULL) {
            ++linenum;

            FIXLINE(s);
            bytesread += strlen(s);

            if (s[0] == '#')            continue; // Comment
            if (s[0] == '%')            continue; // Comment

            if (linenum == 1) {
                char *t = strtok(s, delims);
                if (t == NULL) {
                    LOG(INFO) << "First line must contain num verts and num edges";
                    return false;
                }

                num_vertices = atoi(t);
                t = strtok(NULL, delims);
                if (t != NULL){
                    num_edges = atol(t);
                } else {
                    LOG(INFO) << "First line must contain num verts and num edges";
                    return false;
                }
                LOG(INFO) << "Vertices: " << num_vertices << ", Edges: " << num_edges;

                init_datastructures();
                read_vertexpart();
                continue; 
            }

            vid_t from = linenum - 1; // because first line contained the num of verts and edges
            char *t = strtok(s, delims);
            if (t == NULL)              continue;
            t = strtok(NULL, delims);
            do {
                vid_t to = atoi(t);
                if (from < to) { //ignore one direction, because METIS format contains both directions of an undirected edge
                    edges.emplace_back(edge_t(from, to));
                    bid_t part_u = vertex2bucket[from];
                    bid_t part_v = vertex2bucket[to];

                    bucket_info[part_u].is_mirror.set_bit_unsync(curr_edge_cnt);
                    bucket_info[part_v].is_mirror.set_bit_unsync(curr_edge_cnt);

                    ++curr_edge_cnt;
                }
            } while((t = strtok(NULL, delims)) != NULL);
        }
        free(s);
        fclose(fin);
        return true;
    }

    bool split_in_edgelist()
    {
        curr_edge_cnt = 0;
        std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
        auto filesize = fin.tellg();
        LOG(INFO) << "file size: " << filesize;
        fin.seekg(0, std::ios::beg);

        fin.read((char *)&num_vertices, sizeof(num_vertices));
        fin.read((char *)&num_edges, sizeof(num_edges));

        LOG(INFO) << "num_vertices: " << num_vertices
                << ", num_edges: " << num_edges;
        if (sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t) != (size_t)filesize) {
            LOG(INFO) << "sizeof(vid_t) + sizeof(eid_t) + num_edges * sizeof(edge_t) != filesize";
            return false;
        }

        edges.resize(num_edges);
        init_datastructures();
        read_vertexpart();
        fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

        // std::shuffle(edges.begin(), edges.end(), rd);
        for (eid_t eid = 0; eid < num_edges; ++eid) {
            auto &edge = edges[eid];
            if (eid % 500'000'000 == 0) {
                LOG(INFO) << "edgelist lines read: " << eid;
            }
            vid_t &from = edge.first, &to = edge.second;
            ++from;
            ++to;
            bid_t part_u = vertex2bucket[from];
            bid_t part_v = vertex2bucket[to];

            bucket_info[part_u].is_mirror.set_bit_unsync(curr_edge_cnt);
            bucket_info[part_v].is_mirror.set_bit_unsync(curr_edge_cnt);

            ++curr_edge_cnt;
        }
        fin.close();
        return true;
    }

public:
	Vertex2EdgePart(std::string basefilename, bool need_k_split);

	void split();

};

#endif