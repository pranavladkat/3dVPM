#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"
#include "solver.hpp"
#include "wake.hpp"

using namespace std;

int main()
{

    shared_ptr<Surface> surface(new Surface);

    PLOT3D mesh;

    string filename = "apame.x";

    vector3d free_stream_velocity(1.0,0,0);
    double time_step = 1.0;

    mesh.set_surface(surface);
    mesh.read_mesh(filename);
    mesh.build_topology();

    surface->compute_panel_components();

    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

    shared_ptr<vtk_writer> writer(new vtk_writer());

    Solver solver;
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);

    solver.solve();

    cout << "Hello World!" << endl;
    return 0;
}



/* Code notes:
 * In compute_source_edge_influence(), need to multiply by fabs(z) rather than z.
 * In compute_source_doublet_edge_influence(), need to multiply by fabs(z) rather than z.
 * In Solve(), pass "true" to get_collocation_point() function.
 *
 */
