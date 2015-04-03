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

void Wake :: initialize(const vector3d& free_stream_velocity, const double &dt){

    assert(lifting_surface != NULL);
    assert(lifting_surface->n_trailing_edge_nodes() > 0);

    // compute lifting surface panel components
    lifting_surface->compute_panel_components();

    // allocate memory for nodes
    nodes.resize(2*lifting_surface->n_trailing_edge_nodes());

    // compute node location of the first row of panels (in the free stream)
    for(int n = 0; n < lifting_surface->n_trailing_edge_nodes(); n++ ){

        // get the node at trailing edge
        vector3d node(lifting_surface->nodes[lifting_surface->trailing_edge_nodes[n]]);

        // compute new location of node in the downstream
        vector3d distance(0,0,0);
        // get local velocity at trailing node
        vector3d vel = - lifting_surface->get_kinematic_velocity(node) + free_stream_velocity;

        if(Parameters::static_wake)
            distance = lifting_surface->get_trailing_edge_bisector(n) * Parameters::static_wake_length;
        else
            distance = vel * dt * Parameters::trailing_edge_wake_shed_factor;

        // add computed distance to the node to get new node location
        nodes[n] = node + distance;
    }

    // second row of wake-nodes matches with trailing edge nodes
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

// builds panel neighbors
void Wake :: build_topology(){
    assert(nodes.size() > 0);

    int spanwise_nodes = lifting_surface->n_trailing_edge_nodes();
    int spanwise_panels = lifting_surface->n_trailing_edge_panels();

    int total_nodes = nodes.size();
    int total_panels = n_panels();

    int total_new_panels =  (total_nodes / spanwise_nodes - 1) * spanwise_panels - total_panels;

    for(int p = 0; p < total_new_panels; p++){

        vector<int> new_panel;
        new_panel.clear();

        // compute node numbers of the panels
        int node_a = spanwise_nodes +  n_panels() + n_panels()/spanwise_panels;
        int node_b = n_panels() + n_panels()/spanwise_panels;
        int node_c = node_b + 1;
        int node_d = node_a + 1;

        // add nodes in counter-clockwise
        new_panel.push_back(node_a);
        new_panel.push_back(node_b);
        new_panel.push_back(node_c);
        new_panel.push_back(node_d);

        // push back new panel
        panels.push_back(new_panel);
    }

}


void Wake :: shed_wake(const vector3d &free_stream_velocity, double dt){

    assert(nodes.size() > 0);

    // move nodes on trailing edge with local velocity
    for(int n = nodes.size() - lifting_surface->n_trailing_edge_nodes(); n < (int)nodes.size(); n++ ){

        // get node on trailing edge
        vector3d& node = nodes[n];

        // local wake velocity at trailing edge
        vector3d vel = - lifting_surface->get_kinematic_velocity(node) + free_stream_velocity;

        // compute distance
        vector3d distance(0,0,0);
        if(Parameters::static_wake)
            distance = lifting_surface->get_trailing_edge_bisector(n) * Parameters::static_wake_length;
        else
            distance = vel * dt * Parameters::trailing_edge_wake_shed_factor;

        //update node coordinates
        node += distance;
    }

    // add row of nodes which matches trailing edge
    for(int i = 0; i < lifting_surface->n_trailing_edge_nodes(); i++){
        const vector3d& TE_node = lifting_surface->nodes[lifting_surface->trailing_edge_nodes[i]];
        nodes.push_back(TE_node);
    }

    build_topology();
    compute_panel_components();
}
