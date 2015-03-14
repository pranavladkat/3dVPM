#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"
#include "solver.hpp"
#include "wake.hpp"
#include "domain.hpp"

using namespace std;

int main(int argc, char** args)
{

    Parameters::unsteady_problem = true;
    Parameters::use_vortex_core_model = false;

    // create surface object
    shared_ptr<Surface> surface(new Surface);

    // read mesh file
    PLOT3D mesh;
    string filename = "NACA0012_1.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    // now free stream vel is zero and surface is moving with (-1,0,0)
    vector3d free_stream_velocity(1,0,0);
    //vector3d surface_linear_velocity(0,0,0);
    //vector3d surface_angular_velocity(0,0,72);
    //surface->set_linear_velocity(surface_linear_velocity);
    //surface->set_angular_velocity(surface_angular_velocity,false);

    double time_step = 0.1;
    double fluid_density = 1.225;

    // set blade at AOA
    surface->rotate_surface(vector3d(0,-8,0),false);
    surface->compute_panel_components();

    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

//    shared_ptr<Domain> domain(new Domain());
//    mesh.set_domain(domain);
//    mesh.read_domain("turbine_1.x");

    shared_ptr<vtk_writer> writer(new vtk_writer());

    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

//    vector3d angular_displacement = surface_angular_velocity * (2 * M_PI / 60.0) * time_step;

    for(int i = 0; i < 40; i++){
        solver.solve(time_step,i);
        solver.convect_wake(time_step);
        //surface->translate_surface(surface_linear_velocity * time_step);
        //surface->rotate_surface(angular_displacement,true);
        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();
    }

    cout << "Hello World!" << endl;
    return 0;
}



/* Changes to do in the code:
 * - Test add/subtract kinematic vel in panel coordinate in compute_surface_velociy().
 * - How to compute span_area in wind turbine case?
 */
