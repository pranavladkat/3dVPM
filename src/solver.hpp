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

/**
   Solver class calculates the solution of the given problem.
   Takes input in terms of surface, wake, configurational parameters such as surface velocity, free stream velocity,
   reference velocity, density, etc.
   Calculates the solution in terms of surface velocity, pressure, body forces (and force coefficients).
   Also write output to vtk files.

   @brief Solver class for 3D unsteady panel method
*/

class Solver
{
private:

    /** @brief stores the pointer to a surface object */
    std::shared_ptr<Surface> surface;

    /** @brief stores the pointer to a wake object */
    std::shared_ptr<Wake> wake;

    /** @brief stores the pointer to a vtk_writer object */
    std::shared_ptr<vtk_writer> log;

    /** @brief stores the pointer to a matlab_writer object */
    std::shared_ptr<matlab_writer> mlog;

    /** @brief stores free stream velocity */
    vector3d free_stream_velocity;

    /** @brief stores source strength of each panel in surface object */
    std::vector<double> source_strength;

    /** @brief stores doublet strength of each panel in surface object */
    std::vector<double> doublet_strength;

    /** @brief computes source strength of a panel */
    double compute_source_strength(const int panel) const;

    /** @brief source - influence coefficient matrix */
    std::vector<std::vector<double>> source_influence;

    /** @brief doublet - influence coefficient matrix */
    std::vector<std::vector<double>> doublet_influence;

    /** @brief doublet - influence coefficient matrix, a PETSc object */
    Mat doublet_influence_matrix;

    /** @brief RHS and Solution vector, a PETSc object */
    Vec RHS, solution;

    /** @brief KSP solver context, a PETSc object */
    KSP ksp_doublet;

    /** @brief Sets up a system of equation to be solved using PETSc */
    void setup_linear_system();

    /** @brief allocate memory for PETSc objects */
    void initialize_petsc_variables();

    /** @brief solve the system of equations using PETSc */
    void solve_linear_system();

    /** @brief command line arguments */
    int argc; char** args;

    /** @brief Stores the surface velocity for each panel */
    std::vector<vector3d> surface_velocity;

    /** @brief Compute surface velocity for a panel */
    vector3d compute_surface_velocity(const int panel) const ;

    /** @brief Stores pressure coefficient */
    std::vector<double> pressure_coefficient;

    /** @brief Stores surface potential */
    std::vector<double> surface_potential;

    /** @brief Stores surface potential from previous time step */
    std::vector<double> surface_potential_old;

    /** @brief Stores Calculate pressure coefficient using unsteady bernoulli equation */
    double compute_pressure_coefficient(const int& panel, const int& iteration, const double& dt) const;

    /** @brief compute surface potential of a panel */
    double compute_surface_potential(const int& panel) const;

    /** @brief Stores body forces */
    vector3d body_forces, body_force_coefficients;

    /** @brief Stores reference velocity, used for pressure calculations */
    vector3d reference_velocity;

    /** @brief Stores fluid density */
    double density;

    /** @brief Stores wake panel strengths */
    std::vector<double> wake_doublet_strength;

    /** @brief Compute total velocity at a point */
    vector3d compute_total_velocity(const vector3d& x) const;

    /** @brief doublet - influence coefficient matrix for wake panels */
    std::vector<std::vector<double>> wake_doublet_influence;

    /** @brief computes body forces on surface */
    vector3d compute_body_forces() const ;

    /** @brief computes force coefficients */
    vector3d compute_body_force_coefficients() const ;

public:
    /** @brief constructor - takes command line arguments */
    Solver(int argC,char** argS);
    ~Solver();

    /** @brief attaches surface object with solver */
    void add_surface(const std::shared_ptr<Surface>);

    /** @brief attaches wake object with solver */
    void add_wake(const std::shared_ptr<Wake>);

    /** @brief attaches vtk-logger object with solver */
    void add_logger(const std::shared_ptr<vtk_writer>);

    /** @brief attaches matlab-logger object with solver */
    void add_logger(const std::shared_ptr<matlab_writer>);

    /** @brief Set free stream velocity */
    void set_free_stream_velocity(const vector3d&);

    /** @brief Set reference velocity */
    void set_reference_velocity(const vector3d&);

    /** @brief Set fluid density velocity */
    void set_fluid_density(const double value);

    /** @brief driver function of solver
     *
     * The function is responsible for solution of given problem. The function calculates and influence
     * coefficents, applied boundary conditions, builds system of equations, solves them, post-processes
     * solution, write them in output file. */
    void solve(const double dt, int iteration = 0);

    /** @brief convects the wake with local induced velocity */
    void convect_wake(const double& dt);

    /** @brief compute velocity at each node in the domain */
    void compute_domain_velocity(const std::shared_ptr<Domain> domain);

    /** @brief performs some post-solution operations */
    void finalize_iteration();

    /** @brief Returns body force vectors */
    vector3d get_body_forces() const;

    /** @brief Returns body force coefficents */
    vector3d get_body_force_coefficients() const;

    /** @brief Returns pressure coefficients */
    double get_pressure_coefficient(const int panel) const;

    /** @brief Write output to a file */
    void write_output(const int& interation) const ;

    /** @brief Write output to matlab format - do not use! */
    void write_matlab_output() const;

};


// global functions
/** @brief wrapper class to create PETSc vector */
extern void petsc_vec_create(Vec& vec, int size);
/** @brief wrapper class to create PETSc matrix */
extern void petsc_mat_create(Mat& mat, const int rows, const int cols);
/** @brief write PETSc matrix to a file */
extern void WriteMat(Mat& mat,char const *name);

/** @brief LAPACK solver for least square problem */
extern "C" void dgelsd_( int* m, int* n, int* nrhs, double* a, int* lda,
                         double* b, int* ldb, double* s, double* rcond, int* rank,
                         double* work, int* lwork, int* iwork, int* info );
// solves AX = B , not used in the current work
extern "C" void dgesv_(int *n, int *nrhs,  double *a,  int  *lda,
                       int *ipivot, double *b, int *ldb, int *info);



#endif // SOLVER_H
