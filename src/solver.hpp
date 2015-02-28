#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#include "surface.hpp"
#include "vtk_writer.hpp"


class Solver
{
private:

    std::shared_ptr<Surface> surface;

    std::shared_ptr<vtk_writer> log;

public:
    Solver();
    ~Solver();

    void add_surface(const std::shared_ptr<Surface>);

    void add_logger(const std::shared_ptr<vtk_writer>);

};

#endif // SOLVER_H
