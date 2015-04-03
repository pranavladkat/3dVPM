#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;


int main(int argc, char** args){

    Parameters::unsteady_problem = false;
    Parameters::static_wake = true;
    Parameters::static_wake_length = 10.;

    // create surface object
    shared_ptr<Surface> surface(new Surface);

    // read mesh file
    PLOT3D mesh;
    string filename = "apame.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    //set free stream velocity
    vector3d free_stream_velocity(1,0,0);

    double time_step = 0.5;
    double fluid_density = 1.225;

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
    solver.solve(time_step);
    solver.finalize_iteration();
    cout << "Body Force Coefficients = " << solver.get_body_force_coefficients() << endl;



    return 0;
}
