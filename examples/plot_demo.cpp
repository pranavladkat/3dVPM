#include <iostream>

#include "src/plot3d.hpp"
#include "src/vtk_writer.hpp"
#include "src/solver.hpp"
#include "src/wake.hpp"
#include "src/domain.hpp"

using namespace std;

int main(int argc, char** args){

#ifdef _OPENMP
    omp_set_num_threads(3);
#endif

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
    shared_ptr<matlab_writer> writer(new matlab_writer());
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
    for(int i = 0; i < 120; i++){

        solver.solve(time_step,i);
        if(i < 119){
            solver.convect_wake(time_step);

            surface_velocity = vector3d(0,0,max_amplitude*cos(omega*time));
            surface->set_linear_velocity(surface_velocity);
            surface->translate_surface(surface_velocity*time_step);

            wake->shed_wake(free_stream_velocity,time_step);
        }
        solver.finalize_iteration();

        time += time_step;

    }

    solver.write_matlab_output();
    //writer->write_surface_mesh("Output/wake_mesh",wake);

//    ofstream file("Output/pressure_coefficient.csv");
//    file << "x,y,z,cp" << endl;
//    for(int p = 0; p < surface->n_panels(); p++){
//        file << surface->get_collocation_point(p)[0] <<",";
//        file << surface->get_collocation_point(p)[1] <<",";
//        file << surface->get_collocation_point(p)[2] <<",";
//        file << solver.get_pressure_coefficient(p) << endl;
//    }
//    file.close();
    ofstream file("Output/surface_geometry.csv");
    file << "x,y,z" << endl;
    for(int p = 0; p < surface->n_panels(); p++){
        file << surface->get_collocation_point(p)[0] <<",";
        file << surface->get_collocation_point(p)[1] <<",";
        file << surface->get_collocation_point(p)[2] << endl;
    }
    file.close();
    file.open("Output/wake_geometry.csv");
    file << "x,y,z" << endl;
    for(int p = 0; p < wake->n_panels(); p++){
        file << wake->get_collocation_point(p)[0] <<",";
        file << wake->get_collocation_point(p)[1] <<",";
        file << wake->get_collocation_point(p)[2] << endl;
    }
    file.close();

    cout << "Exiting program." << endl;

    return 0;
}
