#ifndef MATLAB_WRITER_H
#define MATLAB_WRITER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <fstream>
#include <sys/stat.h>

#include "surface.hpp"
#include "domain.hpp"

class matlab_writer
{
private:

    /** @brief the file extension will be assigned based on output file format */
    const std::string file_extension;

public:
    matlab_writer();
    ~matlab_writer();

    template <class T>
    void write_surface_data(std::string filename, std::shared_ptr<Surface> surface, std::vector<T>& data,std::string name, bool writemesh){

        assert(surface->n_panels() == data.size());

        std::ofstream ofile(filename + file_extension, std::ios::app);

        if(sizeof(T) == 8){



        }
        else if (sizeof(T) == 24){

            ofile << "x = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << surface->get_collocation_point(i)[0] << "  ";
            ofile << " ];" << endl;

            ofile << "y = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << surface->get_collocation_point(i)[1] << "  ";
            ofile << " ];" << endl;

            ofile << "z = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << surface->get_collocation_point(i)[2] << "  ";
            ofile << " ];" << endl;

            ofile << "u = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << surface->get_collocation_point(i)[0] << "  ";
            ofile << " ];" << endl;

        }

    }

};

#endif // MATLAB_WRITER_H
