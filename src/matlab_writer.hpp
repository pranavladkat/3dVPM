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

using namespace std;

/** @brief writes output in matlab format. DO NOT USE - NOT a standard format  */

class matlab_writer
{
private:

    /** @brief the file extension will be assigned based on output file format */
    const std::string file_extension;

public:
    matlab_writer();
    ~matlab_writer();

    void write_surface_mesh(std::string filename, const std::shared_ptr<Surface> surface){

        std::ofstream ofile(filename + file_extension);

        // write vertices
        ofile << " vert = [" << endl;
        for(int i = 0; i < surface->n_nodes(); i++)
            ofile << surface->nodes[i] << endl;
        ofile << "];" << endl;

        //write faces
        ofile << "faces = [" << endl;
        for(int p = 0; p < surface->n_panels(); p++){
            for(size_t n = 0; n < surface->panels[p].size(); n++)
                ofile << surface->panels[p][n]+1 << "  ";
            ofile << endl;
        }
        ofile << "];" << endl;

        ofile << "patch('Faces',faces,'Vertices',vert);" << endl;
        //ofile << "axis equal tight; \n colorbar;" << endl;

    }


    template <class T>
    void write_surface_data(std::string filename, const std::shared_ptr<Surface> surface, const std::vector<T>& data,std::string name, bool writemesh = false){

        assert(surface->n_panels() == data.size());

        std::ofstream ofile(filename + ".m");

        if(sizeof(T) == 8){

            // write vertices
            ofile << " vert = [" << endl;
            for(size_t i = 0; i < surface->n_nodes(); i++)
                ofile << surface->nodes[i] << endl;
            ofile << "];" << endl;

            //write faces
            ofile << "faces = [" << endl;
            for(size_t p = 0; p < surface->n_panels(); p++){
                for(size_t n = 0; n < surface->panels[p].size(); n++)
                    ofile << surface->panels[p][n]+1 << "  ";
                ofile << endl;
            }
            ofile << "];" << endl;

            ofile << name << " = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << data[i] << endl;
            ofile << "];" << endl;

            ofile << "patch('Faces',faces,'Vertices',vert,'FaceVertexCData',"<< name <<",'FaceColor','flat');" << endl;
            ofile << "axis equal tight; \n colorbar;" << endl;

        }
        else if (sizeof(T) == 24){

            ofile << "x = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << surface->get_collocation_point(i) << endl;
            ofile << " ];" << endl;


            ofile << "u = [" << endl;
            for(size_t i = 0; i < surface->n_panels(); i++)
                ofile << data[i] << endl;
            ofile << " ];" << endl;

            ofile << "quiver3(x(:,1),x(:,2),x(:,3),u(:,1),u(:,2),u(:,3)); axis equal tight;" << endl;

        }

    }


    template <class T>
    void write_domain_data(std::string filename, std::shared_ptr<Domain> domain, std::vector<T> data,std::string name, bool writemesh){

        assert(domain->n_nodes() == data.size());

        std::ofstream ofile(filename + file_extension);

        if (sizeof(T) == 24){

            ofile << "x = [" << endl;
            for(size_t i = 0; i < domain->n_nodes(); i++)
                ofile << domain->nodes[i] << endl;
            ofile << " ];" << endl;


            ofile << "u = [" << endl;
            for(size_t i = 0; i < domain->n_nodes(); i++)
                ofile << data[i] << endl;
            ofile << " ];" << endl;

            ofile << "quiver3(x(:,1),x(:,2),x(:,3),u(:,1),u(:,2),u(:,3)); axis equal tight;" << endl;

        }



    }

};

#endif // MATLAB_WRITER_H
