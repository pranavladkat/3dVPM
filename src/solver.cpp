#include "solver.hpp"

Solver::Solver()
{

}

Solver::~Solver()
{

}

void Solver :: add_surface(const std::shared_ptr<Surface> surf){
    surface = surf;
}


void Solver :: add_logger(const std::shared_ptr<vtk_writer> writer){
    log = writer;
}
