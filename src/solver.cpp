#include "solver.hpp"

using namespace std;

Solver::Solver(const int argC,char** argS)
    :argc(argC), args(argS)
{
    free_stream_velocity = 0;
    reference_velocity = 0;
    density = 0;
}

Solver::~Solver()
{
    VecDestroy(&RHS);
    VecDestroy(&solution);
    MatDestroy(&doublet_influence_matrix);
    KSPDestroy(&ksp_doublet);
    PetscFinalize();
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

void Solver :: set_reference_velocity(const vector3d& vel){
    reference_velocity = vel;
}

void Solver :: set_fluid_density(const double value){
    density = value;
}

void Solver :: solve(const double dt, int iteration){

    cout << "ITERATION : " << iteration + 1 << endl;
    // compute source strength
    cout << "Computing Source Strengths..." ;
    source_strength.clear();
    source_strength.resize(surface->n_panels());
    for(int p = 0; p < surface->n_panels(); p++){
        source_strength[p] = compute_source_strength(p);
        //cout << std::scientific << source_strength[p] << endl;
    }
    cout << "Done." << endl;

    // compute source and doublet coefficients and build Influence matrix
    cout << "Computing Influence Coefficient Matrix...";
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
    cout << "Done." << endl;

    // apply Kutta-condition
    cout << "Applying Kutta-Condition...";

    // compute wake panel range to consider in kutta-condition
    int wake_panel_start, wake_panel_end;
    if(Parameters::unsteady_problem){
        wake_panel_start = wake->n_panels() - surface->n_trailing_edge_panels();
        wake_panel_end = wake->n_panels();
    }else{
        wake_panel_start = 0;
        wake_panel_end = wake->n_panels();
    }
    int TE_panel_counter = 0;

    for(int wp = wake_panel_start; wp < wake_panel_end; wp++){

        if(TE_panel_counter == surface->n_trailing_edge_panels())
            TE_panel_counter = 0;
        int upper_panel = surface->upper_TE_panels[TE_panel_counter];
        int lower_panel = surface->lower_TE_panels[TE_panel_counter];

        for(int sp = 0; sp < surface->n_panels(); sp++){

            // remember to use negative sign when computing doublet coeff of the wake (as normal is in opposite direction)
            double influence = - wake->compute_doublet_panel_influence(wp,surface->get_collocation_point(sp,false));

            doublet_influence[sp][upper_panel] += influence;
            doublet_influence[sp][lower_panel] -= influence;
            //cout << influence << endl;
        }
        TE_panel_counter++;
    }
    cout << "Done." << endl;

    // compute influence coeficient of the wake doublet
    if(wake_doublet_strength.size() > 0){

        cout << "Computing Wake Influence coefficients...";

        wake_doublet_influence.clear();
        wake_doublet_influence.resize(surface->n_panels(),vector<double>(wake_doublet_strength.size()));
        for(int n = 0; n < surface->n_panels(); n++){
            for(int p = 0; p < (int)wake_doublet_strength.size(); p++)
                wake_doublet_influence[n][p] = - wake->compute_doublet_panel_influence(p,surface->get_collocation_point(n,false));
        }
        cout << "Done." << endl;
    }

    // initialize petsc variables
    if(iteration == 0)
        initialize_petsc_variables();

    // setup A and B of the AX = B
    cout << "Solving for unknown doublet strengths..." << endl;
    setup_linear_system();

    // solve linear system & set unknown doublet strengths
    solve_linear_system();

    // compute surface velocity
    cout << "Computing Surface Velocities...";
    if(iteration == 0){
        surface_velocity.clear();
        surface_velocity.resize(surface->n_panels());
    }
    for(int p = 0; p < surface->n_panels(); p++){
        surface_velocity[p] = compute_surface_velocity(p) ;
    }
    cout << "Done." << endl;

    // compute surface potential
    if(iteration == 0){
        surface_potential.clear();
        surface_potential.resize(surface->n_panels());
    }
    for(int p = 0; p < surface->n_panels(); p++){
        surface_potential[p] = compute_surface_potential(p);
        //cout << surface_potential[p] << endl;
    }


    // compute coefficient of pressure
    cout << "Computing Pressure...";
    if(iteration == 0){
        pressure_coefficient.clear();
        pressure_coefficient.resize(surface->n_panels());
    }
    for(int p = 0; p < surface->n_panels(); p++){
        pressure_coefficient[p] = compute_pressure_coefficient(p,iteration,dt) ;
        //cout << pressure_coefficient[p] << endl;
    }
    cout << "Done." << endl;


    // compute wake strength
    if(Parameters::unsteady_problem){
        wake_panel_start = wake->n_panels() - surface->n_trailing_edge_panels();
        wake_panel_end = wake->n_panels();
    }else{
        wake_panel_start = 0;
        wake_panel_end = wake->n_panels();
    }
    TE_panel_counter = 0;

    for(int wp = wake_panel_start; wp < wake_panel_end; wp++){

        assert(doublet_strength.size() > 0);

        if(TE_panel_counter == surface->n_trailing_edge_panels())
            TE_panel_counter = 0;
        int upper_panel = surface->upper_TE_panels[TE_panel_counter];
        int lower_panel = surface->lower_TE_panels[TE_panel_counter];

        double wake_strenth = doublet_strength[upper_panel] - doublet_strength[lower_panel];
        wake_doublet_strength.push_back(wake_strenth);

        TE_panel_counter++;
    }

    // write iteration output
    write_output(iteration);

}


double Solver::compute_source_strength(const int panel) const{

        const vector3d& node = surface->get_collocation_point(panel,false);
        vector3d vel = surface->get_kinematic_velocity(node) - free_stream_velocity;
        return (vel.dot(surface->get_panel_normal(panel)));
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

    // setup RHS
    double *_RHS;
    VecGetArray(RHS,&_RHS);

    // RHS = source_coefficient * source_strength
    for(int i = 0; i < surface->n_panels(); i++){
        for(int j = 0; j < surface->n_panels(); j++)
            _RHS[i] += source_influence[i][j] * source_strength[j];
        //cout << _RHS[i] << endl;
    }

    // RHS -= wake_doublet_coefficient * wake_doublet_strength
    if(wake_doublet_strength.size() > 0){
        for(int i = 0; i < surface->n_panels(); i++){
            for(int j = 0; j < (int)wake_doublet_strength.size(); j++)
                _RHS[i] -= wake_doublet_influence[i][j] * wake_doublet_strength[j] ;
        }
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
    PetscPrintf(PETSC_COMM_WORLD,"Solution converged in %d iterations.\n",itn+1);

    PetscReal *_SOL;

    VecGetArray(solution,&_SOL);

    /* set doublet strength of panels */
    doublet_strength.clear();
    doublet_strength.resize(surface->n_panels());
    for(int p = 0; p < surface->n_panels(); p++){
        doublet_strength[p] = _SOL[p];
        //cout << doublet_strength[p] << endl;
    }

    VecRestoreArray(solution,&_SOL);
}


vector3d Solver :: compute_surface_velocity(const int panel) const {

    vector3d local_velocity = surface->get_kinematic_velocity(surface->get_collocation_point(panel,false)) - free_stream_velocity;

    vector3d local_velocity_transformed = surface->transform_vector_panel(panel,local_velocity);

    /* computes tangential velocity using Least-Square approach
     * refer : CFD : Principles and Applications, by J. Blazek (pg. 162)
     *
     * Least square problem solved with Lapack's divide and conquor routine : dgelsd_
     * Example impliementation of dgelsd_ :
     * https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
     * More info on LLSQ : http://www.netlib.org/lapack/lug/node27.html
     */

    const vector<int>& neighbour_panels = surface->panel_neighbours[panel];
    int neighbour_size = (int)neighbour_panels.size();
    assert(neighbour_size > 0);

    int dim = 2;    // neglect variation in z
    double rhs[neighbour_size];
    double mat[dim*neighbour_size];

    // setup RHS
    for(int i = 0; i < neighbour_size; i++)
        rhs[i] = doublet_strength[neighbour_panels[i]] - doublet_strength[panel];

    // setup matrix (in column major layout)
    for(int i = 0; i < neighbour_size; i++){

        // transform CP of neighbouring node to panel's coordinates
        vector3d neighbour_node = surface->transform_point_panel(panel,surface->get_collocation_point(neighbour_panels[i],false));

        for(int j = 0; j < dim; j++){
            mat[j*neighbour_size+i] = neighbour_node[j];
        }
    }

    /* solve least square problem */
    /* Local variables to dgelsd_ */
    int m = neighbour_size, n = dim, nrhs = 1, lda = m, ldb = max(m,n), info, lwork, rank;
    double rcond = -1.0;
    double wkopt;
    /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
     * where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
     * and smlsiz = 25 */
    int iwork[11*min(m,n)];
    double s[m];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info );
    lwork = (int)wkopt;
    double work[lwork];

    /* Solve the equations A*X = B */
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );

    // check if solution returned success
    assert(info == 0);

    // notice negative sign on rhs terms
    // also notice third component is kept zero
    vector3d total_velocity = vector3d(-rhs[0],-rhs[1],0) - vector3d(local_velocity_transformed[0],local_velocity_transformed[1],0);

    // transform back to global coordinates
    return surface->transform_vector_panel_inverse(panel,total_velocity);
}



double Solver :: compute_pressure_coefficient(const int& panel, const int& iteration, const double &dt) const {

    //compute dphidt
    double dphidt = 0;
    if(iteration > 0 && Parameters::unsteady_problem){
        assert(dt > 0);
        dphidt = (surface_potential[panel] - surface_potential_old[panel]) / dt ;
    }

    assert(reference_velocity.squared_norm() != 0);

    // compute pressure_coefficient
    double Cp = 1.0 - (surface_velocity[panel].squared_norm() + 2.0 * dphidt) / reference_velocity.squared_norm();

    return Cp;
}


double Solver :: compute_surface_potential(const int& panel) const {

    double potential;

    // phi = - doublet_strength
    potential = - doublet_strength[panel];

    // add contribution from free stream velocity and body velocity (phi_infinity = U*x+V*y+W*z)
    vector3d local_velocity = surface->get_kinematic_velocity(surface->get_collocation_point(panel,false)) - free_stream_velocity;
    potential -= local_velocity.dot(surface->get_collocation_point(panel,false));

    return potential;
}


vector3d Solver :: compute_body_forces() const {

    assert(density > 0);

    vector3d Force(0,0,0);

    double dynamic_pressure = 0.5 * density * reference_velocity.squared_norm();

    // compute force
    for(int p = 0; p < surface->n_panels(); p++){
        Force = Force - surface->get_panel_normal(p) * dynamic_pressure * pressure_coefficient[p] * surface->get_panel_area(p);
    }

    return Force;
}


vector3d Solver :: compute_body_force_coefficients() const {

    double dynamic_pressure = 0.5 * density * reference_velocity.squared_norm();

    // compute planform-area
    // S = 0.5 * sum_of ( panel_area * vector_normal_to_free_stream )
    // need to varify for wind-turbine case
    double planform_area = 0;
    for(int p = 0; p < surface->n_panels(); p++){
        planform_area += fabs(surface->get_panel_normal(p)[2] * surface->get_panel_area(p));
    }
    planform_area *= 0.5 ;

    return body_forces / (dynamic_pressure * planform_area);
}


vector3d Solver :: compute_total_velocity(const vector3d& x) const {

    assert(source_strength.size()  > 0);
    assert(doublet_strength.size() > 0);

    vector3d velocity(0,0,0);

    // compute velocity due to surface panels
    for(int sp = 0; sp < surface->n_panels(); sp++){
        velocity += surface->compute_source_panel_unit_velocity(sp,x)  * source_strength[sp];
        velocity += surface->compute_doublet_panel_unit_velocity(sp,x) * doublet_strength[sp];
    }

    // compute velocity due to shedded wake panels
    for(int wp = 0; wp < (int)wake_doublet_strength.size(); wp++){
        velocity -= wake->compute_doublet_panel_unit_velocity(wp,x) * wake_doublet_strength[wp];
    }

    // add free stream velocity
    velocity += free_stream_velocity;

    return velocity;
}



void Solver :: convect_wake(const double& dt){

    // this function only convects nodes which are already present in the wake (except trailing edge nodes)
    // must call shed_wake() function to add new wake panel row

    if(Parameters::static_wake == false){
        assert(wake->n_panels() > 0);

        // nodes_to_convect does not include nodes on trailing edge.
        int nodes_to_convect = wake->n_nodes() - surface->n_trailing_edge_nodes();
        assert(nodes_to_convect > 0);

        // compute velocity at the wake nodes which needs to be convected
        vector<vector3d> wake_velocity(nodes_to_convect);

        for(int wn = 0; wn < nodes_to_convect; wn++)
            wake_velocity[wn] = compute_total_velocity(wake->nodes[wn]);

        // move the nodes with wake velocity
        for(int wn = 0; wn < nodes_to_convect; wn++)
            wake->nodes[wn] += wake_velocity[wn] * dt ;
    }

}



void Solver :: compute_domain_velocity(const std::shared_ptr<Domain> domain){

    assert(domain->n_nodes() > 0);
    vector<vector3d> domain_velocity(domain->n_nodes());

    for(int n = 0; n < domain->n_nodes(); n++)
        //domain_velocity[n] = compute_total_velocity(domain->nodes[n]);
        domain_velocity[n] = surface->get_kinematic_velocity(domain->nodes[n]);

    log->write_domain_data("solver-out-domain",domain,domain_velocity,"V",true);
}




void Solver :: write_output(const int& iteration) const {

    // folder name
    const string foldername = "Output";

    // create directory if does not exist
    struct stat dir_err;
    if(stat(foldername.c_str(),&dir_err) == -1)
        mkdir(foldername.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // set surface and wake names
    string surface_name = foldername + "/" + "surface_" + to_string(iteration);
    string wake_name = foldername + "/" + "wake_" + to_string(iteration);

    // write data
    cout << "writing output...";

    log->write_surface_data(surface_name,surface,surface_velocity,"V",true);
    log->write_surface_data(surface_name,surface,pressure_coefficient,"CP",false);
    log->write_surface_data(surface_name,surface,doublet_strength,"Mu",false);
    log->write_surface_data(surface_name,surface,source_strength,"Sigma",false);
    log->write_surface_mesh(wake_name,wake);
    if(wake_doublet_strength.size() > 0)
        log->write_surface_data(wake_name,wake,wake_doublet_strength,"Mu",false);

    cout << "Done." << endl;

}



void Solver :: finalize_iteration(){

    // if problem is steady
    if(!Parameters::unsteady_problem){

        // clear doublet strength of the wake
        wake_doublet_strength.clear();

    }

    // if problem is unsteady
    else if (Parameters::unsteady_problem) {

        // save surface potential in previous variable
        surface_potential_old = surface_potential;

    }

}














// global functions

// create petsc vec
void petsc_vec_create(Vec& vec, int size){

    VecCreate(PETSC_COMM_WORLD,&vec);
    VecSetSizes(vec,PETSC_DECIDE,size);
    VecSetFromOptions(vec);
    VecSet(vec,0.0);
}

// create petsc mat
void petsc_mat_create(Mat& mat, const int rows, const int cols){

    MatCreateSeqDense(PETSC_COMM_WORLD,rows,cols,NULL,&mat);
    MatSetFromOptions(mat);
    MatSetUp(mat);
    MatZeroEntries(mat);
}

// write petsc mat to a file readable by MATLAB
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
