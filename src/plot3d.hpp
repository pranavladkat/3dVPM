#ifndef PLOT3D_H
#define PLOT3D_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

#include "surface.hpp"
#include "domain.hpp"

/**
   @brief PLOT3D class to read plot3D mesh file and manipulates surface object
*/

class PLOT3D
{
private:

    // filename to read
    std::string surface_filename, domain_filename;

    /** @brief stores the pointer to a surface object */
    std::shared_ptr<Surface> surface;

    /** @brief stores the pointer to a domain object */
    std::shared_ptr<Domain> domain;

    /** @brief each block represents a separate Surface */
    int blocks;

    /** @brief number of nodes in x and y directions. Nodes in z direction will always be 1. */
    int IMAX, JMAX;

    /** @brief decides if normals need to be flipped
      *
      * In case the block orientation is not right, normals may get calculated in opposite direction.
      * Set \flip_normal to true to flip the normal.
      */
    bool flip_normal;

public:
    PLOT3D();
    ~PLOT3D();

    /** @brief set file name to read */
    void set_surface_filename(std::string);

    /** @brief connect to surface object  */
    void set_surface(std::shared_ptr<Surface>);

    /** @brief Flips the normals of the surface */
    void flip_normals(bool);

    /** @brief read the plot3d mesh file */
    void read_surface(std::string name);

    /** @brief buid topology from mesh file */
    void build_topology();

    /** @brief connect to domain object */
    void set_domain(std::shared_ptr<Domain>);

    /** @brief read mesh file for the domain */
    void read_domain(std::string name);


};

#endif // PLOT3D_H
