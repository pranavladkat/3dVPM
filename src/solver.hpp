#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <iomanip>

#include "surface.hpp"
#include "wake.hpp"
#include "vtk_writer.hpp"

#include "petscksp.h"


class Solver
{
private:

    std::shared_ptr<Surface> surface;

    std::shared_ptr<Wake> wake;

    std::shared_ptr<vtk_writer> log;

    vector3d free_stream_velocity;

    std::vector<double> source_strength;

    double compute_source_strength(const int panel) const;

    std::vector<std::vector<double>> source_influence;

    std::vector<std::vector<double>> doublet_influence;

    Mat Doublet_influence_matrix;

    Vec RHS, solution;

public:
    Solver();
    ~Solver();

    void add_surface(const std::shared_ptr<Surface>);

    void add_wake(const std::shared_ptr<Wake>);

    void add_logger(const std::shared_ptr<vtk_writer>);

    void set_free_stream_velocity(const vector3d&);

    void solve();

};

#endif // SOLVER_H
