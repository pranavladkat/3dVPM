#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"
#include "solver.hpp"
#include "wake.hpp"
#include "domain.hpp"

using namespace std;

int main(int argc, char** args)
{

    Parameters::unsteady_problem = false;
    Parameters::static_wake = true;
    Parameters::static_wake_length = 1;

    shared_ptr<Surface> surface(new Surface);

    PLOT3D mesh;

    string filename = "NACA0012_1.x";

    double speed = 1.0;
    double angle = 10 * M_PI/180.0;
    vector3d free_stream_velocity(speed*cos(angle),0,speed*sin(angle));
    //vector3d free_stream_velocity(1,0,0);

    double time_step = 1;
    double fluid_density = 1.225;

    mesh.set_surface(surface);
    mesh.read_surface(filename);
    mesh.build_topology();

    //surface->rotate_surface(vector3d(0,-5,0),false);

    surface->compute_panel_components();

    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

//    shared_ptr<Domain> domain(new Domain());
//    mesh.set_domain(domain);
//    mesh.read_domain("NACA0012_downstream_domain.x");

    shared_ptr<vtk_writer> writer(new vtk_writer());

    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

    solver.solve(time_step);
    //solver.convect_wake(time_step);


    //solver.compute_domain_velocity(domain);


    cout << "Hello World!" << endl;
    return 0;
}



/* Changes to do in the code:
 * - While applying Kutta-condition, make sure only trailing edge wake panels are considered
 * - Test calc of lower and upper panel during kutta condition for wake panels > spanwise panels
 * - Test add/subtract kinematic vel in panel coordinate in compute_surface_velociy().
 * - How to compute span_area in wind turbine case?
 */
