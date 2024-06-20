#ifndef PART_WRITER_HPP
#define PART_WRITER_HPP

#include "util.hpp"

template <typename vertex_type, typename proc_type>
class EdgepartWriterBase
{
protected:
    std::vector<std::ofstream> fout_;
    size_t size_pre_line_;

    using new_proc_type = std::conditional_t<std::is_same_v<proc_type, uint8_t>, uint16_t, proc_type>;
    
public:
    EdgepartWriterBase(const std::string &basefilename) 
    {
        size_pre_line_ = sizeof(vertex_type) + sizeof(vertex_type) + sizeof(new_proc_type);
    }

    ~EdgepartWriterBase() {}

    virtual void save_edge(vertex_type from, vertex_type to, proc_type proc)
    {
        return;
    }
};


template <typename vertex_type, typename proc_type>
class EdgepartWriterOnefile : public EdgepartWriterBase<vertex_type, proc_type>
{
    using EdgepartWriterBase<vertex_type, proc_type>::fout_;
    using typename EdgepartWriterBase<vertex_type, proc_type>::new_proc_type;

public:
    EdgepartWriterOnefile(const std::string &basefilename)
        : EdgepartWriterBase<vertex_type, proc_type>(basefilename)
    {
        std::string filename = edge_partitioned_name(basefilename);
        fout_.push_back(std::ofstream(filename, std::ios_base::trunc));
        if (!fout_.back().is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            exit(1);
        }
    }
    
    void save_edge(vertex_type from, vertex_type to, proc_type proc) override
    {
        fout_[0] << from << ' ' << to << ' ' << (new_proc_type)proc << std::endl;
    }
};

template <typename vertex_type, typename proc_type>
class EdgepartWriterMultifile : public EdgepartWriterBase<vertex_type, proc_type>
{
private:
    using EdgepartWriterBase<vertex_type, proc_type>::fout_;
    using typename EdgepartWriterBase<vertex_type, proc_type>::new_proc_type;

public:
    EdgepartWriterMultifile(const std::string &basefilename)
        : EdgepartWriterBase<vertex_type, proc_type>(basefilename)
    {
        for (bid_t b = 0; b < FLAGS_p; ++b) {
            std::string filename = edge_partitioned_name(basefilename) + ".part" + std::to_string(b);
            fout_.push_back(std::ofstream(filename, std::ios_base::trunc));
            if (!fout_.back().is_open()) {
                std::cerr << "Error opening file: " << filename << std::endl;
                exit(1);
            }
        }
    }

    ~EdgepartWriterMultifile() {}

    void save_edge(vertex_type from, vertex_type to, proc_type proc) override
    {
        fout_[(new_proc_type)proc] << from << ' ' << to << std::endl;
    }
};


/// @brief: ONLY SUPPORT ONEFILE MODE
template <typename vertex_type, typename proc_type>
struct VertexpartWriter
{
    std::ofstream fout;
    bool write;

    using new_proc_type = std::conditional_t<std::is_same_v<proc_type, uint8_t>, uint16_t, proc_type>;

    VertexpartWriter(const std::string &basefilename, bool write)
        : fout(write ? std::ofstream(vertex_partitioned_name(basefilename), std::ios_base::trunc) : std::ofstream()), write(write)
    {
    }

    ~VertexpartWriter() {}

    void save_vertex(vertex_type vid, proc_type proc)
    {
        if (write)
        {
            fout << (new_proc_type)proc << std::endl;
        }
    }
};

#endif