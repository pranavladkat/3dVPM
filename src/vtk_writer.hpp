#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <fstream>
#include <sys/stat.h>

#include "surface.hpp"
#include "domain.hpp"

class vtk_writer
{
private:

    /** @brief the file extension will be assigned based on output file format */
    const std::string file_extension;

public:
    vtk_writer();
    ~vtk_writer();

    template <class T>
    void write_surface_data(std::string filename, const std::shared_ptr<Surface> surface, const std::vector<T>& data,std::string name, bool writemesh){

        if(writemesh){
            write_surface_mesh(filename,surface);
        }

        std::ofstream ofile(filename + file_extension, std::ios::app);

        if(data.size() == surface->n_nodes())
            ofile << "POINT_DATA " << surface->nodes.size() << std::endl;

        switch (sizeof(T)) {
        case 4:
            ofile << "SCALARS " << name << " double 1" << std::endl;
            ofile << "LOOKUP_TABLE default" << std::endl;
            break;
        case 8:
            ofile << "SCALARS " << name << " double 1" << std::endl;
            ofile << "LOOKUP_TABLE default" << std::endl;
            break;
        case 24:
            ofile << "VECTORS " << name << " double" << std::endl;
            break;
        default:
            std::cerr << "UNKNOWN DATATYPE!!!" << std::endl;
            break;
        }

        for(int p = 0; p < (int)data.size(); p++)
            ofile << data[p] << std::endl;
        ofile << std::endl;
    }

    void write_surface_mesh(std::string filename,std::shared_ptr<Surface>);

    void write_domain_mesh(std::string filename,std::shared_ptr<Domain>);

    template <class T>
    void write_domain_data(std::string filename, std::shared_ptr<Domain> domain, std::vector<T> data,std::string name, bool writemesh){

        if(writemesh){
            write_domain_mesh(filename,domain);
        }

        std::ofstream ofile(filename + ".vtk", std::ios::app);

        switch (sizeof(T)) {
        case 4:
            ofile << "SCALARS " << name << " double 1" << std::endl;
            ofile << "LOOKUP_TABLE default" << std::endl;
            break;
        case 8:
            ofile << "SCALARS " << name << " double 1" << std::endl;
            ofile << "LOOKUP_TABLE default" << std::endl;
            break;
        case 24:
            ofile << "VECTORS " << name << " double" << std::endl;
            break;
        default:
            std::cerr << "UNKNOWN DATATYPE!!!" << std::endl;
            break;
        }

        for(int p = 0; p < (int)data.size(); p++)
            ofile << data[p] << std::endl;
        ofile << std::endl;
    }

};

#endif // VTK_WRITER_H
