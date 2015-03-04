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

class PLOT3D
{
private:

    std::string surface_filename, domain_filename;

    std::shared_ptr<Surface> surface;

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

    void set_surface_filename(std::string);

    void set_surface(std::shared_ptr<Surface>);

    /** @brief Flips the normals of the surface */
    void flip_normals(bool);

    void read_surface(std::string name);

    void build_topology();

    void set_domain(std::shared_ptr<Domain>);

    void read_domain(std::string name);


};

#endif // PLOT3D_H
