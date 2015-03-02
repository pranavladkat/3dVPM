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

    for(int p = 0; p < surface->n_panels(); p++){
        for(int n = 0; n < surface->n_panels(); n++){
            cout << surface->compute_source_panel_influence(p,surface->get_collocation_point(n,false)) << endl;
        }
    }


}


double Solver::compute_source_strength(const int panel) const{

        vector3d& node = surface->get_collocation_point(panel,false);
        vector3d vel = free_stream_velocity - surface->get_kinematic_velocity(node);
        return -(vel.dot(surface->get_panel_normal(panel)));
}
