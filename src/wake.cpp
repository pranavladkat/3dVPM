#include "wake.hpp"

using namespace std;

Wake::Wake()
{
    lifting_surface = NULL;
}

Wake::~Wake()
{

}

void Wake :: add_lifting_surface(const std::shared_ptr<Surface> surf) {
    lifting_surface = surf;
}

void Wake :: initialize(const vector3d& free_stream_velocity, double dt){

    assert(lifting_surface != NULL);
    assert(lifting_surface->n_trailing_edge_nodes() > 0);

    nodes.resize(2*lifting_surface->n_trailing_edge_nodes());

    for(int n = 0; n < lifting_surface->n_trailing_edge_nodes(); n++ ){

        vector3d distance(0,0,0);
        vector3d node(lifting_surface->nodes[lifting_surface->trailing_edge_nodes[n]]);
        vector3d vel = lifting_surface->get_kinematic_velocity(node) + free_stream_velocity;

        if(Parameters::static_wake)
            distance = lifting_surface->get_trailing_edge_bisector(n) * Parameters::static_wake_length;
        else
            distance = vel * dt * Parameters::trailing_edge_wake_shed_factor;

        nodes[n] = node + distance;
    }

    int i = 0;
    for(size_t n = lifting_surface->n_trailing_edge_nodes(); n < nodes.size(); n++ ){
        nodes[n] = lifting_surface->nodes[lifting_surface->trailing_edge_nodes[i]];
        i++;
    }

    // build panels
    build_topology();
    // compute panel components
    compute_panel_components();

}


void Wake :: build_topology(){
    assert(nodes.size() > 0);

    int spanwise_nodes = lifting_surface->n_trailing_edge_nodes();
    int spanwise_panels = lifting_surface->n_trailing_edge_panels();

    int total_nodes = nodes.size();
    int total_panels = n_panels();

    int total_new_panels =  (total_nodes / spanwise_nodes - 1) * spanwise_panels - total_panels;

    int it = 0;
    if(total_panels > 0)
        it++;

    for(int p = 0; p < total_new_panels; p++){

        if(it == spanwise_panels)
            it++;

        vector<int> new_panel;
        new_panel.clear();
        // panel order same as that of surface (anti-clockwise)
        new_panel.push_back((total_panels + spanwise_nodes) + it);
        new_panel.push_back((total_panels) + it);
        new_panel.push_back((total_panels + 1) + it);
        new_panel.push_back((total_panels + 1 + spanwise_nodes) + it);

        // panel order opposite to that of surface (clockwise)
        //new_panel.push_back((total_panels) + it);
        //new_panel.push_back((total_panels + spanwise_nodes) + it);
        //new_panel.push_back((total_panels + 1 + spanwise_nodes) + it);
        //new_panel.push_back((total_panels + 1) + it);

        panels.push_back(new_panel);
        it++;
    }

}


void Wake :: shed_wake(){

}
