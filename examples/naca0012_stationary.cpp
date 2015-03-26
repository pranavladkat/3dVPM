#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;

int main(int argc, char** args){

    Parameters::unsteady_problem = true;

    // create surface object
    shared_ptr<Surface> surface(new Surface);

    // read mesh file
    PLOT3D mesh;
    string filename = "NACA0012_1.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    //set free stream velocity
    vector3d free_stream_velocity(1,0,0);

    double time_step = 0.5;
    double fluid_density = 1.225;

    // set blade at AOA
    surface->rotate_surface(vector3d(0,-2.0,0),false);
    surface->compute_panel_components();

    // create wake object
    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

    // create writer
    shared_ptr<vtk_writer> writer(new vtk_writer());
    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

    // solve
    for(int i = 0; i < 50; i++){
        solver.solve(time_step,i);
        solver.convect_wake(time_step);
        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();
    }

    cout << "Body Force Coefficients = " << solver.get_body_force_coefficients() << endl;

    return 0;
}
