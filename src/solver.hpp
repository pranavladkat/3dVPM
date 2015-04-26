#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <iomanip>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "surface.hpp"
#include "wake.hpp"
#include "vtk_writer.hpp"
#include "matlab_writer.hpp"

#include "petscksp.h"


class Solver
{
private:

    std::shared_ptr<Surface> surface;

    std::shared_ptr<Wake> wake;

    std::shared_ptr<vtk_writer> log;

    std::shared_ptr<matlab_writer> mlog;

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

    void solve_linear_system();

    int argc; char** args;

    std::vector<vector3d> surface_velocity;

    vector3d compute_surface_velocity(const int panel) const ;

    std::vector<double> pressure_coefficient;

    std::vector<double> surface_potential;

    std::vector<double> surface_potential_old;

    double compute_pressure_coefficient(const int& panel, const int& iteration, const double& dt) const;

    double compute_surface_potential(const int& panel) const;

    vector3d body_forces, body_force_coefficients;

    vector3d reference_velocity;

    double density;

    std::vector<double> wake_doublet_strength;

    vector3d compute_total_velocity(const vector3d& x) const;

    std::vector<std::vector<double>> wake_doublet_influence;

    vector3d compute_body_forces() const ;

    vector3d compute_body_force_coefficients() const ;

public:
    Solver(int argC,char** argS);
    ~Solver();

    void add_surface(const std::shared_ptr<Surface>);

    void add_wake(const std::shared_ptr<Wake>);

    void add_logger(const std::shared_ptr<vtk_writer>);

    void add_logger(const std::shared_ptr<matlab_writer>);

    void set_free_stream_velocity(const vector3d&);

    void set_reference_velocity(const vector3d&);

    void set_fluid_density(const double value);

    void solve(const double dt, int iteration = 0);

    void convect_wake(const double& dt);

    void compute_domain_velocity(const std::shared_ptr<Domain> domain);

    void finalize_iteration();

    vector3d get_body_forces() const;

    vector3d get_body_force_coefficients() const;

    double get_pressure_coefficient(const int panel) const;

    void write_output(const int& interation) const ;

    void write_matlab_output() const;

};


// global functions
extern void petsc_vec_create(Vec& vec, int size);
extern void petsc_mat_create(Mat& mat, const int rows, const int cols);
extern void WriteMat(Mat& mat,char const *name);

extern "C" void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                         double* b, int* ldb, double* s, double* rcond, int* rank,
                         double* work, int* lwork, int* iwork, int* info );

extern "C" void dgesv_(int *n, int *nrhs,  double *a,  int  *lda,
                       int *ipivot, double *b, int *ldb, int *info);



#endif // SOLVER_H
