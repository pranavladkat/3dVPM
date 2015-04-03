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

    // read mesh file
    PLOT3D mesh;
    string filename = "blade_edited.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    // set free stream velocity
    vector3d free_stream_velocity(0,7,0);

    //set angular velocity
    vector3d surface_angular_velocity(0,72,0);
    surface->set_angular_velocity(surface_angular_velocity,false);

    double time_step = 0.025;
    double fluid_density = 1.225;

    surface->compute_panel_components();
    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

    shared_ptr<vtk_writer> writer(new vtk_writer());

    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

    vector3d angular_displacement = surface_angular_velocity * (2 * M_PI / 60.0) * time_step;

    for(int i = 0; i < 30; i++){

        solver.solve(time_step,i);
        solver.convect_wake(time_step);
        surface->rotate_surface(angular_displacement,true);
        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();
    }

    cout << "Exiting program." << endl;

    return 0;
}
