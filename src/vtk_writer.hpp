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

/** @brief vtk_writer class to write output in vtk format */

class vtk_writer
{
private:

    /** @brief the file extension will be assigned based on output file format */
    const std::string file_extension;

public:
    vtk_writer();
    ~vtk_writer();

    /** @brief writes surface data to a file
     *
     * @param[in] filename  filename of the output file, extension will be attaced automatically
     * @param[in] surface   surface to write
     * @param[in] data      data to write for a given \b surface
     * @param[in] name      name of the \b data variable
     * @param[in] writemesh true if mesh data needs to be written, else false
     */
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

    /** @brief write mesh data to vtk format */
    void write_surface_mesh(std::string filename,std::shared_ptr<Surface>);

    /** @brief write domain mesh to vtk format */
    void write_domain_mesh(std::string filename,std::shared_ptr<Domain>);

    /** @brief writes surface data to a file
     *
     * @param[in] filename  filename of the output file, extension will be attaced automatically
     * @param[in] domain    domain to write
     * @param[in] data      data to write for a given \b domain
     * @param[in] name      name of the \b data variable
     * @param[in] writemesh true if mesh data needs to be written, else false
     */
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
