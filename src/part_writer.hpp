#ifndef PART_WRITER_HPP
#define PART_WRITER_HPP

#include "util.hpp"

template <typename vertex_type, typename proc_type>
struct edgepart_writer {
    char *buffer;
    std::ofstream fout;
    bool write;

    using new_proc_type = std::conditional_t<std::is_same_v<proc_type, uint8_t>, uint16_t, proc_type>;

    edgepart_writer(const std::string &basefilename, bool write)
        : fout(edge_partitioned_name(basefilename), write ? std::ios_base::trunc : std::ios_base::app), write(write)
    {
        size_t s = sizeof(vertex_type) + sizeof(vertex_type) + sizeof(new_proc_type);
        buffer = new char[sizeof(char) + s];
    }

    ~edgepart_writer() { delete[] buffer; }

    void save_edge(vertex_type from, vertex_type to, proc_type proc)
    {
        if (write) {
            fout << from << ' ' << to << ' ' << (new_proc_type)proc << std::endl;
        }
    }
};

template <typename vertex_type, typename proc_type>
struct vertexpart_writer {
    char *buffer;
    std::ofstream fout;
    bool write;

    using new_proc_type = std::conditional_t<std::is_same_v<proc_type, uint8_t>, uint16_t, proc_type>;

    vertexpart_writer(const std::string &basefilename, bool write)
        : fout(vertex_partitioned_name(basefilename), write ? std::ios_base::trunc : std::ios_base::app), write(write)
    {
        size_t s = sizeof(vertex_type) + sizeof(proc_type) + sizeof(new_proc_type);
        buffer = new char[sizeof(char) + s];
    }

    ~vertexpart_writer() { delete[] buffer; }

    void save_vertex(vertex_type vid, proc_type proc)
    {
        if (write) {
            fout << (new_proc_type)proc << std::endl;
        }
    }
};

#endif