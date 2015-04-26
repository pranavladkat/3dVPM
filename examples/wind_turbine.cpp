#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;

// turn off assertions
#define NDEBUG

int main(int argc, char** args)
{

    Parameters::unsteady_problem = true;

    // create surface object
    shared_ptr<Surface> surface(new Surface);

    // read mesh file
    PLOT3D mesh;
    string filename = "blade_5.5_modified.x";
    mesh.set_surface(surface);
    mesh.read_surface(filename);

    // set free stream velocity
    vector3d free_stream_velocity(0,7,0);

    //set angular velocity
    vector3d surface_angular_velocity(0,71.63,0);
    surface->set_angular_velocity(surface_angular_velocity,false);

    double time_step = 0.015;
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

    //write coordinates of the collocation points to data variable
    vector<vector<double>> data (4,vector<double>(surface->n_panels()));
    for(int p = 0; p < surface->n_panels(); p++){
        const vector3d& point = surface->get_collocation_point(p);
        data[0][p] = point[0];
        data[1][p] = point[1];
        data[2][p] = point[2];
    }

    for(int i = 0; i < 90; i++){

        solver.solve(time_step,i);
        solver.convect_wake(time_step);
        surface->rotate_surface(angular_displacement,true);
        wake->shed_wake(free_stream_velocity,time_step);
        solver.finalize_iteration();
    }

    //write Cp values to data
    for(int p = 0; p < surface->n_panels(); p++){
        data[3][p] = solver.get_pressure_coefficient(p);
    }

    // write to file
    ofstream file("Output/pressure_data.csv");
    file << "x,y,z,cp" << endl;
    for(int p = 0; p < surface->n_panels(); p++){
        file << scientific << data[0][p] << "," << data[1][p]
             << "," << data[2][p] << "," << data[3][p] << endl;
    }
    file.close();

    cout << "Exiting program." << endl;

    return 0;
}
