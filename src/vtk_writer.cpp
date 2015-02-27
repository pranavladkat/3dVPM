#include "vtk_writer.hpp"

using namespace std;

vtk_writer::vtk_writer(string name)
    :file_extension(".vtk"), filename(name)
{

}

vtk_writer::~vtk_writer()
{

}


void vtk_writer :: write_mesh(std::shared_ptr<Surface> surface){

    string name = filename + file_extension;

    ofstream ofile(name);

    //write vtk file info
    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << "OUTPUT by 3D-VPM\n";
    ofile << "ASCII" << endl;
    ofile << "DATASET UNSTRUCTURED_GRID" << endl;

    int total_nodes = surface->nodes.size();
    int total_panels = surface->n_panels();
    int total_panel_nodes = 0;
    for(int i = 0; i < surface->n_panels(); i++)
        total_panel_nodes += surface->panels[i].size() + 1;

    //write node data
    ofile << "POINTS " << total_nodes << " double" << endl;

    for(int i = 0; i < (int)surface->nodes.size(); i++)
        ofile << surface->nodes[i] << endl;
    ofile << endl;

    //write cell data
    ofile << "CELLS " << total_panels << " " << total_panel_nodes << endl;

    for(int p = 0; p < surface->n_panels(); p++){
        ofile << surface->panels[p].size() << " ";
        for(int n = 0; n < (int)surface->panels[p].size(); n++ )
            ofile << surface->panels[p][n] << " ";
        ofile << endl;
    }
    ofile << endl;

    //write cell type
        ofile << "CELL_TYPES " << total_panels << endl;

        for(int p = 0; p < (int)surface->panels.size(); p++){
            int cell_type;
            switch (surface->panels[p].size()) {
            case 3:
                cell_type = 5;
                break;
            case 4:
                cell_type = 9;
                break;
            default:
                cerr << "Panel " << p << ": Unknown cell type." << endl;
                continue;
            }

            ofile << cell_type << endl;
        }
        ofile << endl;

        ofile << "CELL_DATA " << surface->n_panels() << std::endl;
}

