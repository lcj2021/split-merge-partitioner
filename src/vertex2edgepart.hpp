#ifndef VERTEX2EDGEPART_HPP
#define VERTEX2EDGEPART_HPP

#include <random>
#include <boost/unordered_map.hpp>

#include "dense_bitset.hpp"
#include "part_writer.hpp"
#include "partitioner.hpp"

class Vertex2EdgePart : public Partitioner {
private:
	std::string basefilename;

    // vid_t num_vertices;
    // size_t num_edges;
    int k, p;
    std::random_device rd;
    std::mt19937 gen;

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

    void init_datastructures()
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
    std::unordered_map<int, int> fast_merge();
    std::unordered_map<int, int> precise_merge();

	void read_vertexpart();

	int vertex2edgepartID(vid_t u, vid_t v);

    bool split_in_adjlist()
    {
        // read the metis adjacency list
        curr_edge_cnt = 0;
        FILE *fin = fopen((basefilename + ".adjlist").c_str(), "r");
        if (fin == NULL) {
            LOG(INFO) << "Could not load:" << basefilename
                        << " error: " << strerror(errno) << std::endl;
            return false;
        }
        LOG(INFO) << "Reading in adjacency list format!" << std::endl;
        LOG(INFO) << basefilename;

        int maxlen = 1000000000;
        char *s = (char *)malloc(maxlen);

        size_t bytesread = 0;

        char delims[] = " \t";
        size_t linenum = 0;

        while (fgets(s, maxlen, fin) != NULL) {
            linenum++;

            FIXLINE(s);
            bytesread += strlen(s);

            if (s[0] == '#')            continue; // Comment
            if (s[0] == '%')            continue; // Comment

            if (linenum == 1) {
                char *t = strtok(s, delims);
                if (t == NULL) {
                    LOG(INFO) << "First line must contain num verts and num edges" << std::endl; // empty line
                    return false;
                }

                num_vertices = atoi(t);
                t = strtok(NULL, delims);
                if (t != NULL){
                    num_edges = atol(t);
                } else {
                    LOG(INFO) << "First line must contain num verts and num edges" << std::endl;
                    return false;
                }
                LOG(INFO) << "Vertices: " << num_vertices << ", Edges: " << num_edges << std::endl;

                init_datastructures();
                read_vertexpart();
                continue; 
            }

            // LOG(INFO) << "***********";
            vid_t from = linenum - 1; // because first line contained the num of verts and edges
            char *t = strtok(s, delims);
            if (t == NULL)              continue;
            t = strtok(NULL, delims);
            do {
                vid_t to = atoi(t);
                if (from < to) { //ignore one direction, because METIS format contains both directions of an undirected edge
                    edges.emplace_back(edge_t(from, to));
                    int part_u = vertex2partition[from];
                    int part_v = vertex2partition[to];

                    bucket_info[part_u].is_mirror.set_bit_unsync(curr_edge_cnt);
                    bucket_info[part_v].is_mirror.set_bit_unsync(curr_edge_cnt);

                    curr_edge_cnt ++;
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
        if (sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t) != (size_t)filesize) {
            LOG(INFO) << "sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t) != filesize";
            return false;
        }

        edges.resize(num_edges);
        init_datastructures();
        read_vertexpart();
        fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

        std::shuffle(edges.begin(), edges.end(), rd);
        for (auto &edge : edges) {
            vid_t &from = edge.first, &to = edge.second;
            from ++, to ++;
            int part_u = vertex2partition[from];
            int part_v = vertex2partition[to];

            bucket_info[part_u].is_mirror.set_bit_unsync(curr_edge_cnt);
            bucket_info[part_v].is_mirror.set_bit_unsync(curr_edge_cnt);

            curr_edge_cnt ++;
        }
        fin.close();
        return true;
    }

public:
	Vertex2EdgePart(std::string basefilename);
	virtual ~Vertex2EdgePart();

	void split();

};

#endif