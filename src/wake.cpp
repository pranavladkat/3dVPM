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
        vector3d node(lifting_surface->nodes[lifting_surface->trailing_edge_nodes[n]]);
        vector3d vel = lifting_surface->get_kinematic_velocity(node) + free_stream_velocity;
        node[0] = node[0] + vel[0] * dt *  Parameters::trailing_edge_wake_shed_factor;
        node[1] = node[1] + vel[1] * dt *  Parameters::trailing_edge_wake_shed_factor;
        node[2] = node[2] + vel[2] * dt *  Parameters::trailing_edge_wake_shed_factor;

        if(Parameters::wake_shed_along_TE_bisector){
            vector3d vec = lifting_surface->get_trailing_edge_bisector(n);
            node[0] = node[0] * vec[0];
            node[1] = node[1] * vec[1];
            node[2] = node[2] * vec[2];
        }

        nodes[n] = node;
    }

    int i = 0;
    for(size_t n = lifting_surface->n_trailing_edge_nodes(); n < nodes.size(); n++ ){
        nodes[n] = lifting_surface->nodes[lifting_surface->trailing_edge_nodes[i]];
        i++;
    }

}
