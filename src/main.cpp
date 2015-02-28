#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"
#include "solver.hpp"

using namespace std;

int main()
{

    shared_ptr<Surface> surface(new Surface);

    PLOT3D mesh;

    string filename = "plate_mesh_2.x";

    mesh.set_surface(surface);
    mesh.read_mesh(filename);
    mesh.build_topology();

    surface->compute_panel_components();

    shared_ptr<vtk_writer> writer(new vtk_writer(filename));

    writer->write_mesh(surface);

    Solver solver;

    solver.add_surface(surface);
    solver.add_logger(writer);

    cout << "Hello World!" << endl;
    return 0;
}

