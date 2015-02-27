#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"

using namespace std;

int main()
{

    shared_ptr<Surface> surface(new Surface);

    PLOT3D mesh;

    mesh.set_surface(surface);
    mesh.read_mesh("plate_mesh_2.x");
    mesh.build_topology();

    surface->compute_panel_components();

    vtk_writer writer("plate");

//    writer.write(surface,surface->panel_normals,"Normals",true);

    cout << "Hello World!" << endl;
    return 0;
}

