#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;

int main(int argc, char** args){

    Parameters::unsteady_problem = true;
    Parameters::trailing_edge_wake_shed_factor = 1.0;

    // create surface object
    shared_ptr<Surface> surface(new Surface);

    // read mesh file
    PLOT3D mesh;
    string filename = "cylinder.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    double time_step = 0.1;
    double fluid_density = 1.225;

    //set free stream velocity
    vector3d free_stream_velocity(0,0,0);
    vector3d surface_velocity(-1,0,0);
    surface->set_linear_velocity(surface_velocity);

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
    solver.set_fluid_density(fluid_density);


    vector3d acceleration(1.2,0,0);

    for(int i = 0; i < 10; i++){

        cout << "Instataneous Surface Velocity = " << surface_velocity << endl;

        solver.set_reference_velocity(surface_velocity);

        solver.solve(time_step,i);
        surface->translate_surface(surface_velocity*time_step);
        solver.convect_wake(time_step);
        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();

        surface_velocity -= acceleration * time_step;
        surface->set_linear_velocity(surface_velocity);

        cout << "Added mass coefficient = " << solver.get_body_forces()[0]/acceleration[0] << endl;

    }

    cout << "Exiting program." << endl;

    return 0;
}
