#ifndef COMMON_HPP
#define COMMON_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define repv(i, n) for (vid_t i = 0; i < n; ++i)

DECLARE_int32(p);
DECLARE_int32(k);
DECLARE_string(filename);
DECLARE_string(filetype);
DECLARE_string(method);
DECLARE_bool(write);

DECLARE_double(hdf);
DECLARE_double(lambda);
DECLARE_bool(write_low_degree_edgelist);
DECLARE_bool(extended_metrics);
DECLARE_bool(random_streaming);
DECLARE_bool(hybrid_NE);

using vid_t = uint32_t;
const vid_t INVALID_VID = -1;
const vid_t offset = (vid_t)1 << 31;
struct edge_with_id_t {
    vid_t first, second;
    vid_t eid;
    edge_with_id_t() : first(0), second(0), eid(0) {}
    edge_with_id_t(vid_t first, vid_t second, vid_t eid) : first(first), second(second), eid(eid) {}
};

struct edge_t {
    vid_t first, second;
    edge_t() : first(0), second(0) {}
    edge_t(vid_t first, vid_t second) : first(first), second(second) {}
    // const bool valid() { return first != INVALID_VID; }
    // void remove() { first = INVALID_VID; }
    bool valid() const { return first < offset; }
    void remove() { first += offset; }
    void recover() {
        if (first >= offset) {
            first -= offset;
        }
    }
    bool operator == (const edge_t & rhs) const {
        return first == rhs.first && second == rhs.second;
    }
};

#endif