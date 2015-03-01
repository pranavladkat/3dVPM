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
        vector3d vec = lifting_surface->get_trailing_edge_bisector(n);
        node[0] = node[0] + free_stream_velocity[0] * dt * vec[0] * Parameters::trailing_edge_wake_shed_factor;
        node[1] = node[1] + free_stream_velocity[1] * dt * vec[1] * Parameters::trailing_edge_wake_shed_factor;
        node[2] = node[2] + free_stream_velocity[2] * dt * vec[2] * Parameters::trailing_edge_wake_shed_factor;

        cout << node << endl;

    }

    //cout << lifting_surface->n_trailing_edge_nodes() << endl;

}
