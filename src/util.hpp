#pragma once

#include <utility>
#include <chrono>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "threadpool11/threadpool11.hpp"

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

typedef uint32_t vid_t;
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

extern threadpool11::Pool pool;

void preada(int f, char *buf, size_t nbytes, size_t off);
void reada(int f, char *buf, size_t nbytes);
void writea(int f, char *buf, size_t nbytes);

inline std::string h2hedgelist_name(const std::string &basefilename)
{
	std::stringstream ss;
	ss << basefilename << ".h2h_edgelist";
	return ss.str();
}

inline std::string lowedgelist_name(const std::string &basefilename)
{
	std::stringstream ss;
	ss << basefilename << ".low_edgelist";
	return ss.str();
}

inline std::string binedgelist_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".binedgelist";
    return ss.str();
}

inline std::string degree_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".degree";
    return ss.str();
}

inline std::string partitioned_name(const std::string &basefilename)
{
    std::string ret = basefilename + ".edgepart.";
    if (FLAGS_method.substr(0, 3) == "smp") {
        std::string split_method = FLAGS_method == "smp" ? "ne" : FLAGS_method.substr(4);
        ret += "smp_" + split_method + "_k_" + std::to_string(FLAGS_k) + ".";
    } else if (FLAGS_method == "hep") {
        ret += "hep_hdf_" + std::to_string((int)FLAGS_hdf) + ".";
    } else {
        ret += FLAGS_method + ".";
    }
    ret += std::to_string(FLAGS_p);
    return ret;
}

inline bool is_exists(const std::string &name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

class Timer
{
  private:
    std::chrono::system_clock::time_point t1, t2;
    double total;

  public:
    Timer() : total(0) {}
    void reset() { total = 0; }
    void start() { t1 = std::chrono::system_clock::now(); }
    void stop()
    {
        t2 = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = t2 - t1;
        total += diff.count();
    }
    double get_time() { return total; }
};
