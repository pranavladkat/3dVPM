#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

#include "surface.hpp"

class vtk_writer
{
private:

    /** @brief the output file name */
    const std::string filename;

    /** @brief the file extension will be assigned based on output file format */
    const std::string file_extension;

public:
    vtk_writer(std::string);
    ~vtk_writer();

    template <class T>
    void write(std::shared_ptr<Surface> surface, std::vector<T> data,std::string name, bool writemesh){

        if(writemesh){
            write_mesh(surface);
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

    void write_mesh(std::shared_ptr<Surface>);

};

#endif // VTK_WRITER_H
