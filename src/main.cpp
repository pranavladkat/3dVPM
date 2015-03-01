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
    shared_ptr<vtk_writer> writer(new vtk_writer(filename));
    writer->write_mesh(surface);

    Solver solver;
    solver.add_surface(surface);
    solver.add_logger(writer);

    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);
    wake->build_topology();

    shared_ptr<vtk_writer> wake_write(new vtk_writer("wake"));
    wake_write->write_mesh(wake);

    cout << "Hello World!" << endl;
    return 0;
}

