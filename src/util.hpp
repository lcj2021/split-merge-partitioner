#ifndef UTIL_HPP
#define UTIL_HPP

#include "common.hpp"

inline std::string h2hedgelist_name(const std::string &basefilename)
{
	return basefilename + ".h2h_edgelist";
}

inline std::string lowedgelist_name(const std::string &basefilename)
{
	return basefilename + ".low_edgelist";
}

inline std::string binedgelist_name(const std::string &basefilename)
{
    return basefilename + ".binedgelist";
}

inline std::string degree_name(const std::string &basefilename)
{
    return basefilename + ".degree";
}

inline std::string edge_partitioned_name(const std::string &basefilename)
{
    std::string ret = basefilename + ".edgepart.";
    if (FLAGS_method.substr(0, 3) == "fsm") {
        std::string split_method = FLAGS_method == "fsm" ? "ne" : FLAGS_method.substr(4);
        ret += "fsm_" + split_method + "_k_" + std::to_string(FLAGS_k) + ".";
    } else if (FLAGS_method == "hep") {
        ret += "hep_hdf_" + std::to_string((int)FLAGS_hdf) + ".";
    } else {
        ret += FLAGS_method + ".";
    }
    ret += std::to_string(FLAGS_p);
    return ret;
}

inline std::string vertex_partitioned_name(const std::string &basefilename)
{
    std::string ret = basefilename + ".vertexpart.";
    ret += FLAGS_method + ".";
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

#endif