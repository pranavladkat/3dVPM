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

void Solver :: solve(){

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

}


double Solver::compute_source_strength(const int panel) const{

        vector3d& node = surface->get_collocation_point(panel,false);
        vector3d vel = free_stream_velocity - surface->get_kinematic_velocity(node);
        return -(vel.dot(surface->get_panel_normal(panel)));
}



