#include "solver.hpp"

using namespace std;

Solver::Solver()
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

            pair<double,double> influence = surface->compute_source_doublet_panel_influence(p,surface->get_collocation_point(n,false));
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
            double influence = - wake->compute_doublet_panel_influence(wp,surface->get_collocation_point(sp,false));

            doublet_influence[sp][upper_panel] += influence;
            doublet_influence[sp][lower_panel] -= influence;
        }
        TE_panel_counter++;
    }

//    for(int p = 0; p < surface->n_panels(); p++){
//        for(int n = 0; n < surface->n_panels(); n++){
//            cout << fixed /*<< setprecision(6)*/ << setw(10) << scientific << doublet_influence[p][n] << "  ";
//        }
//        cout << endl;
//    }

    // initialize petsc variables
    if(iteration == 0)
        initialize_petsc_variables();

    // setup A and B of the AX = B
    setup_linear_system();
    solve_linear_system();

}


double Solver::compute_source_strength(const int panel) const{

        vector3d& node = surface->get_collocation_point(panel,false);
        vector3d vel = free_stream_velocity - surface->get_kinematic_velocity(node);
        return -(vel.dot(surface->get_panel_normal(panel)));
}

void Solver :: initialize_petsc_variables(){

    assert(surface->n_panels() > 0);

    // create PETSc Vec RHS and solution
    VecCreate(PETSC_COMM_WORLD,&RHS);
    VecSetSizes(RHS,PETSC_DECIDE,surface->n_panels());
    VecSetFromOptions(RHS);
    VecSet(RHS,0.0);
    VecDuplicate(RHS,&solution);
    VecSet(solution,0.0);

    // create PETSc Mat doublet_coefficient_matrix
    MatCreateSeqDense(PETSC_COMM_WORLD,surface->n_panels(),surface->n_panels(),NULL,&doublet_influence_matrix);
    MatSetFromOptions(doublet_influence_matrix);
    MatSetUp(doublet_influence_matrix);
    MatZeroEntries(doublet_influence_matrix);

    // create KSP solver
    KSPCreate(PETSC_COMM_WORLD,&ksp);

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

    KSPSetOperators(ksp,doublet_influence_matrix,doublet_influence_matrix);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    KSPSetType(ksp,KSPPREONLY);
    KSPSetTolerances(ksp,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    KSPSolve(ksp,RHS,solution);
    KSPGetIterationNumber(ksp,&itn);
    PetscPrintf(PETSC_COMM_WORLD,"Iterations taken for KSP: %d\n",itn+1);

    PetscReal *_SOL;

    VecGetArray(solution,&_SOL);

    /* set doublet strength of panels */
    doublet_strength.clear();
    doublet_strength.resize(surface->n_panels());
    for(int p = 0; p < surface->n_panels(); p++){
        doublet_strength[p] = _SOL[p];
        cout << doublet_strength[p] << endl;
    }

    VecRestoreArray(solution,&_SOL);

}
