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

    double time_step = 0.025;
    double fluid_density = 1.225;

    /* Refer Low Speed Aerodynamics - J. Katz, A. Plotkin
     * second edition, page 417-418
     */
    double chord_length = 1.;
    //set free stream velocity
    vector3d free_stream_velocity(0.009*chord_length/time_step,0,0);
    //compute amplitude and frequency
    double max_amplitude = 0.019 * chord_length;
    double reduced_frequency = 8.57;
    double omega = reduced_frequency * 2 * free_stream_velocity[0] / chord_length;

    // create wake object
    shared_ptr<Wake> wake(new Wake());
    wake->add_lifting_surface(surface);
    wake->initialize(free_stream_velocity,time_step);

    // create writer
    shared_ptr<vtk_writer> writer(new vtk_writer());

    // create solver
    Solver solver(argc,args);
    solver.add_surface(surface);
    solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

    vector3d surface_velocity(0,0,0);

    double time = 0;

    // solve
    for(int i = 0; i < 300; i++){

        solver.solve(time_step,i);
        solver.convect_wake(time_step);

        surface_velocity = vector3d(0,0,max_amplitude*cos(omega*time));
        surface->set_linear_velocity(surface_velocity);
        surface->translate_surface(surface_velocity*time_step);

        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();

        time += time_step;

    }

    cout << "Exiting program." << endl;

    return 0;
}
