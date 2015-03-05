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

    std::vector<double> doublet_strength;

    double compute_source_strength(const int panel) const;

    std::vector<std::vector<double>> source_influence;

    std::vector<std::vector<double>> doublet_influence;

    Mat doublet_influence_matrix;

    Vec RHS, solution;

    KSP ksp_doublet;

    void setup_linear_system();

    void initialize_petsc_variables();

    void release_petsc_variables();

    void solve_linear_system();

    int argc; char** args;

    std::vector<vector3d> surface_velocity;

    double compute_surface_velocity(const int panel) const ;

public:
    Solver(int argC,char** argS);
    ~Solver();

    void add_surface(const std::shared_ptr<Surface>);

    void add_wake(const std::shared_ptr<Wake>);

    void add_logger(const std::shared_ptr<vtk_writer>);

    void set_free_stream_velocity(const vector3d&);

    void solve(int iteration = 0);

};


extern void petsc_vec_create(Vec& vec, int size);
extern void petsc_mat_create(Mat& mat, const int rows, const int cols);
extern void WriteMat(Mat& mat,char const *name);

extern "C" void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                         double* b, int* ldb, double* s, double* rcond, int* rank,
                         double* work, int* lwork, int* iwork, int* info );

#endif // SOLVER_H
