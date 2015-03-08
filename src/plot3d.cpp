#include "plot3d.hpp"

using namespace std;

PLOT3D::PLOT3D()
{
    flip_normal = false;
}

PLOT3D::~PLOT3D()
{

}

void PLOT3D :: set_surface_filename(std::string name){
    surface_filename = name;
}


void PLOT3D :: set_surface(std::shared_ptr<Surface> surf){
    surface = surf;
}


void PLOT3D :: flip_normals(bool val){
    flip_normal = val;
}


void PLOT3D :: read_surface(std::string name){

    set_surface_filename(name);

    // open mesh file to read
    ifstream mesh(surface_filename);

    // assert if mesh file is open
    assert(mesh.is_open());

    // read number of blocks
    mesh >> blocks;

    assert(blocks = 1);

    int KMAX;
    // read number of nodes in each direction
    mesh >> IMAX >> JMAX >> KMAX;

    //assert number of nodes in z direction is equal to 1
    assert(KMAX = 1);

    surface->nodes.resize(IMAX*JMAX);

    //read_coordinates
    for(int dim = 0; dim < 3; dim++){
        for(int i = 0; i < IMAX; i++){
            for(int j = 0; j < JMAX; j++){

                // convert 2d index to 1d
                int index = i*JMAX + j;
                mesh >> surface->nodes[index][dim];

            } /* end j loop */
        } /* end i loop */
    } /* end dim loop */

//    for(size_t n = 0; n < surface->nodes.size(); n++)
//        cout << surface->nodes[n] << endl;


    //populate panels with its node numbers
    for(int j = 0; j < JMAX - 1; j++){
        for(int i = 0; i < IMAX - 1; i++){

            // plot3d being structured grid, all panels will have 4 nodes
            vector<int> new_panel(4);

            if(flip_normal){
                new_panel[0] = j*IMAX + i ;
                new_panel[1] = j*IMAX + (i+1) ;
                new_panel[2] = (j+1)*IMAX + (i+1) ;
                new_panel[3] = (j+1)*IMAX + i ;
            }else{
                new_panel[0] = j*IMAX + (i+1) ;
                new_panel[1] = j*IMAX + i ;
                new_panel[2] = (j+1)*IMAX + i ;
                new_panel[3] = (j+1)*IMAX + (i+1) ;
            }

            //cout << new_panel[0] << "\t" << new_panel[1] << "\t" << new_panel[2] << "\t" << new_panel[3] << endl;

            surface->panels.push_back(new_panel);

        } /* end i loop */
    }/* end j loop */

//    for(size_t p = 0; p < surface->panels.size(); p++){
//        for(size_t n = 0; n < surface->panels[p].size(); n++){
//            cout << surface->nodes[surface->panels[p][n]] << endl;
//        }
//        cout << endl;
//    }

}



void PLOT3D :: build_topology(){

    //compute neighbours
    int imax = IMAX - 1;
    int jmax = JMAX - 1;

    for(int j = 0; j < jmax; j++){
        for(int i = 0; i < imax; i++){

            vector<int> new_neighbours;

            if(i == 0){
                new_neighbours.push_back(j * imax + (i+1));
            }else if(i == imax - 1){
                new_neighbours.push_back(j * imax + (i-1));
            }else{
                new_neighbours.push_back(j * imax + (i+1));
                new_neighbours.push_back(j * imax + (i-1));
            }

            if(j == 0){
                new_neighbours.push_back((j+1) * imax + i);
            }else if(j == jmax - 1){
                new_neighbours.push_back((j-1) * imax + i);
            }else{
                new_neighbours.push_back((j+1) * imax + i);
                new_neighbours.push_back((j-1) * imax + i);
            }

//            for(int l = 0; l < new_neighbours.size(); l++){
//                cout << new_neighbours[l] << "\t";
//            } cout << endl;

            surface->panel_neighbours.push_back(new_neighbours);

        } /* end i loop */
    } /* end j loop */


    //compute trailing edge nodes
    for(int j = 0; j < JMAX; j++){
        for(int i = 0; i < IMAX; i++){

            if(i == 0){
                surface->trailing_edge_nodes.push_back(j*IMAX+i);
                //cout << j*IMAX+i << endl;
            }
        } /* end i loop */
    } /* end j loop */

    //compute trailing edge panels
    imax = IMAX - 1;
    jmax = JMAX - 1;

    for(int j = 0; j < jmax; j++){
        for(int i = 0; i < imax; i++){

            if(i == 0){
                if(flip_normal){
                    surface->upper_TE_panels.push_back(j*imax+i);
                }else{
                    surface->lower_TE_panels.push_back(j*imax+i);
                }
            }else if(i == imax - 1){
                if(flip_normal){
                    surface->lower_TE_panels.push_back(j*imax+i);
                }else{
                    surface->upper_TE_panels.push_back(j*imax+i);
                }
            }

        } /* end i loop */
    } /* end j loop */

//    for(int l = 0; l < surface->upper_TE_panels.size(); l++){
//        cout << surface->upper_TE_panels[l] << endl;
//    }

}


void PLOT3D :: set_domain(std::shared_ptr<Domain> dom){
    domain = dom;
}


void PLOT3D :: read_domain(std::string name){

    domain_filename = name;

    // open mesh file to read
    ifstream mesh(domain_filename);

    // assert if mesh file is open
    assert(mesh.is_open());

    // read number of blocks
    mesh >> blocks;

    assert(blocks = 1);

    int KMAX;
    // read number of nodes in each direction
    mesh >> IMAX >> JMAX >> KMAX;

    domain->set_domain_ranges(IMAX,JMAX,KMAX);

    domain->nodes.resize(IMAX*JMAX*KMAX);

    //read_coordinates
    for(int dim = 0; dim < 3; dim++){
        int index = 0;
        for(int k = 0; k < KMAX; k++){
            for(int j = 0; j < JMAX; j++){
                for(int i = 0; i < IMAX; i++){

                    //cout << index << endl;
                    mesh >> domain->nodes[index][dim];
                    index++;

                } /* end i loop */
            } /* end j loop */
        } /* end k loop */
    } /* end dim loop */

//    for(size_t n = 0; n < domain->nodes.size(); n++)
//        cout << domain->nodes[n] << endl;


}
