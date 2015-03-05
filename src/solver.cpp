#include "solver.hpp"

using namespace std;

Solver::Solver(const int argC,char** argS)
    :argc(argC), args(argS)
{
    free_stream_velocity = 0.0;
}

Solver::~Solver()
{

}

void Solver :: add_surface(const std::shared_ptr<Surface> surf){
    surface = surf;
}


void Solver :: add_wake(const std::shared_ptr<Wake> surf){
    wake = surf;
}


void Solver :: add_logger(const std::shared_ptr<vtk_writer> writer){
    log = writer;
}


void Solver :: set_free_stream_velocity(const vector3d& vel){
    free_stream_velocity = vel;
}

void Solver :: solve(int iteration){

    source_strength.clear();
    source_strength.resize(surface->n_panels());

    // compute source strength
    for(int p = 0; p < surface->n_panels(); p++){
        source_strength[p] = compute_source_strength(p);
        //cout << std::scientific << source_strength[p] << endl;
    }

    // compute source and doublet coefficients and populate matrices
    source_influence.clear();
    source_influence.resize(surface->n_panels(),vector<double>(surface->n_panels()));
    doublet_influence.clear();
    doublet_influence.resize(surface->n_panels(),vector<double>(surface->n_panels()));

    for(int n = 0; n < surface->n_panels(); n++){
        for(int p = 0; p < surface->n_panels(); p++){

            pair<double,double> influence = surface->compute_source_doublet_panel_influence(p,surface->get_collocation_point(n,true));
            if(p == n)
                influence.second = -0.5;

            source_influence[n][p]  = influence.first;
            doublet_influence[n][p] = influence.second;

            //cout << std::scientific << influence.first << "\t" << influence.second << endl;
        }
    }


    // apply Kutta-condition
    int TE_panel_counter = 0;
    for(int wp = 0; wp < wake->n_panels(); wp++){

        if(TE_panel_counter == surface->n_trailing_edge_panels())
            TE_panel_counter = 0;
        int upper_panel = surface->upper_TE_panels[TE_panel_counter];
        int lower_panel = surface->lower_TE_panels[TE_panel_counter];

        for(int sp = 0; sp < surface->n_panels(); sp++){

            // remember to use negative sign when computing doublet coeff of the wake (as normal is in opposite direction)
            double influence = - wake->compute_doublet_panel_influence(wp,surface->get_collocation_point(sp,true));

            doublet_influence[sp][upper_panel] += influence;
            doublet_influence[sp][lower_panel] -= influence;
        }
        TE_panel_counter++;
    }

    // initialize petsc variables
    if(iteration == 0)
        initialize_petsc_variables();

    // setup A and B of the AX = B
    setup_linear_system();

    // solve linear system
    solve_linear_system();

    log->write_surface_data("solver-out",surface,doublet_strength,"mu",true);

    compute_surface_velocity(4);

    release_petsc_variables();
}


double Solver::compute_source_strength(const int panel) const{

        vector3d& node = surface->get_collocation_point(panel,false);
        vector3d vel = free_stream_velocity - surface->get_kinematic_velocity(node);
        return -(vel.dot(surface->get_panel_normal(panel)));
}

void Solver :: initialize_petsc_variables(){

    assert(surface->n_panels() > 0);

    // initialize PETSc
    PetscInitialize(&argc,&args,(char*)0,NULL);

    // create PETSc Vec RHS and solution, and PETSc Mat doublet_influence_matrix
    petsc_vec_create(RHS     ,surface->n_panels());
    petsc_vec_create(solution,surface->n_panels());
    petsc_mat_create(doublet_influence_matrix,surface->n_panels(),surface->n_panels());

    // create KSP solver
    KSPCreate(PETSC_COMM_WORLD,&ksp_doublet);

}

void Solver :: setup_linear_system(){

    // clear previous data, if any.
    VecSet(RHS,0.0);
    VecSet(solution,0.0);
    MatZeroEntries(doublet_influence_matrix);

    int col[surface->n_panels()];
    double val[surface->n_panels()];

    // setup column number for matrix insertion
    for(int i = 0; i < surface->n_panels(); i++)
        col[i] = i;

    /* copy doublet_influence to doublet_influence_matrix */
    for(int i = 0;i < surface->n_panels(); i++){
        for(int j = 0; j < surface->n_panels(); j++){
            val[j] = doublet_influence[i][j];
        }
        MatSetValues(doublet_influence_matrix,1,&i,surface->n_panels(),col,val,INSERT_VALUES);
    }
    MatAssemblyBegin(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);
    MatSetOption(doublet_influence_matrix,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE);

    double *_RHS;
    VecGetArray(RHS,&_RHS);

    for(int i = 0; i < surface->n_panels(); i++){
        for(int j = 0; j < surface->n_panels(); j++)
            _RHS[i] += source_influence[i][j] * source_strength[j];
    }
    VecRestoreArray(RHS,&_RHS);
}

void Solver :: solve_linear_system(){

    int itn;
    PC pc;

    KSPSetOperators(ksp_doublet,doublet_influence_matrix,doublet_influence_matrix);
    KSPGetPC(ksp_doublet,&pc);
    PCSetType(pc,PCLU);
    KSPSetType(ksp_doublet,KSPPREONLY);
    KSPSetTolerances(ksp_doublet,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp_doublet);
    KSPSetUp(ksp_doublet);
    KSPSolve(ksp_doublet,RHS,solution);
    KSPGetIterationNumber(ksp_doublet,&itn);
    PetscPrintf(PETSC_COMM_WORLD,"Iterations taken for KSP: %d\n",itn+1);

    PetscReal *_SOL;

    VecGetArray(solution,&_SOL);

    /* set doublet strength of panels */
    doublet_strength.clear();
    doublet_strength.resize(surface->n_panels());
    for(int p = 0; p < surface->n_panels(); p++)
        doublet_strength[p] = _SOL[p];

    VecRestoreArray(solution,&_SOL);
}

void Solver :: release_petsc_variables(){
    VecDestroy(&RHS);
    VecDestroy(&solution);
    MatDestroy(&doublet_influence_matrix);
    KSPDestroy(&ksp_doublet);
    PetscFinalize();
}


double Solver :: compute_surface_velocity(const int panel) const {

    // computes tangential velocity using Least-Square approach
    // refer : CFD : Principles and Applications, by J. Blazek (pg. 162)

    const vector<int>& neighbour_panels = surface->panel_neighbours[panel];
    assert(neighbour_panels.size() > 0);

    Vec rhs, sol;
    Mat mat;
    KSP ksp;
    PC pc;

    petsc_vec_create(rhs,neighbour_panels.size());
    petsc_vec_create(sol,3);
    petsc_mat_create(mat,neighbour_panels.size(),3);
    KSPCreate(PETSC_COMM_WORLD,&ksp);

    PetscReal *_rhs, *_sol;
    int col[3] = {0,1,2};
    VecGetArray(rhs,&_rhs);
    VecGetArray(sol,&_sol);

    // set RHS
    for(int i = 0; i < (int)neighbour_panels.size(); i++){
        _rhs[i] = doublet_strength[neighbour_panels[i]] - doublet_strength[panel];
        cout << _rhs[i] << endl;
    }

    // set matrix
    for(int i = 0; i < (int)neighbour_panels.size(); i++){
        vector3d node = surface->transform_point_panel(panel,surface->get_collocation_point(neighbour_panels[i],false));
        double val[3];
        val[0] = node[0]; val[1] = node[1]; val[2] = node[2];
        MatSetValues(mat,1,&i,3,col,val,INSERT_VALUES);
    }
    MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);

    //WriteMat(mat,"LSQ");

    KSPSetOperators(ksp,mat,mat);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);
    KSPSetType(ksp,KSPLSQR);
    //KSPSetTolerances(ksp,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    KSPSolve(ksp,rhs,sol);

    for(int i = 0; i < 3; i++)
        cout << _sol[i] << endl;

    VecRestoreArray(rhs,&_rhs);
    VecRestoreArray(sol,&_sol);
    VecDestroy(&rhs);
    VecDestroy(&sol);
    MatDestroy(&mat);
    KSPDestroy(&ksp);

    return 0;
}























// global functions

void petsc_vec_create(Vec& vec, int size){

    VecCreate(PETSC_COMM_WORLD,&vec);
    VecSetSizes(vec,PETSC_DECIDE,size);
    VecSetFromOptions(vec);
    VecSet(vec,0.0);
}


void petsc_mat_create(Mat& mat, const int rows, const int cols){

    MatCreateSeqDense(PETSC_COMM_WORLD,rows,cols,NULL,&mat);
    MatSetFromOptions(mat);
    MatSetUp(mat);
    MatZeroEntries(mat);
}

void WriteMat(Mat& mat,char const *name){

    PetscViewer viewer;
    char filename[64] = "Mat_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);
    PetscObjectSetName((PetscObject)mat,name);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(mat,viewer);
    PetscViewerDestroy(&viewer);
}
