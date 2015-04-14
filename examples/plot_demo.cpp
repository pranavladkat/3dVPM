
#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;


int main(int argc, char** args)
{

    Parameters::unsteady_problem = true;

    // create surface object
    shared_ptr<Surface> surface(new Surface);
    // create domain object
    shared_ptr<Domain> domain(new Domain());

    // read mesh file
    PLOT3D mesh;
    string filename = "1panel.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    PLOT3D domain_mesh;
    filename = "box.x";
    mesh.set_domain(domain);
    mesh.read_domain(filename);

    surface->compute_panel_components();
    shared_ptr<vtk_writer> writer(new vtk_writer());

    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_logger(writer);

    solver.compute_domain_velocity(domain);

    writer->write_surface_mesh("Output/surface-mesh",surface);

    cout << "Exiting program." << endl;

    return 0;
}
