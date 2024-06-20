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

DECLARE_int32(p);
DECLARE_int32(k);
DECLARE_string(filename);
DECLARE_string(filetype);
DECLARE_string(method);
DECLARE_string(write);

DECLARE_double(hdf);
DECLARE_double(lambda);
DECLARE_bool(write_low_degree_edgelist);
DECLARE_bool(extended_metrics);
DECLARE_bool(hybrid_NE);

using vid_t = uint32_t;
using eid_t = uint64_t;
using bid_t = uint8_t;
const vid_t kInvalidVid = std::numeric_limits<vid_t>::max();
const bid_t kInvalidBid = std::numeric_limits<bid_t>::max();
const vid_t offset = (vid_t)1 << 31;

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