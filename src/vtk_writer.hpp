#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

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
    void write_surface_data(std::string filename, std::shared_ptr<Surface> surface, std::vector<T> data,std::string name, bool writemesh){

        if(writemesh){
            write_surface_mesh(filename,surface);
        }

        std::ofstream ofile(filename + file_extension, std::ios::app);

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

        std::ofstream ofile(filename + ".m", std::ios::app);

        ofile << name << " = [" << std::endl;
        for(int p = 0; p < (int)data.size(); p++)
            ofile << data[p] << std::endl;
        ofile << "];" << std::endl;

        switch (sizeof(T)) {
        case 4:
            std::string scalar_f = "5," + name + "(:));";
            ofile << ofile << "scatter3(x(:,1),x(:,2),x(:,3)," << scalar_f << std::endl;
            ofile << "axis equal tight;" << std::endl;
            break;
        case 8:
            std::string scalar_d = "5," + name + "(:));";
            ofile << ofile << "scatter3(x(:,1),x(:,2),x(:,3)," << scalar_d << std::endl;
            ofile << "axis equal tight;" << std::endl;
            break;
        case 24:
            std::string quiver = name + "(:,1)," + name + "(:,2)," + name + "(:,3));";
            ofile << "quiver(x(:,1),x(:,2),x(:,3)," << quiver << std::endl;
            break;
        default:
            std::cerr << "UNKNOWN DATATYPE!!!" << std::endl;
            break;
        }

    }

};

#endif // VTK_WRITER_H
