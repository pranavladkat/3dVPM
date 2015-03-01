#include "surface.hpp"

using namespace std;

Surface::Surface()
{
    linear_velocity = 0.0;
    angular_velocity = 0.0;
    surface_origin = 0.0;
    surface_orientation = 0.0;
}

Surface::~Surface()
{

}

int Surface :: n_panels() const{
    return static_cast<int> (panels.size());
}


void Surface :: compute_panel_components(){

    assert(nodes.size() > 0 && panels.size() > 0);

    //compute panel normal and areas
    for(int p = 0; p < n_panels(); p++){

        vector<int> &panel_nodes = panels[p];
        vector3d normal;
        if (panel_nodes.size() == 3) {
            vector3d AB = nodes[panel_nodes[1]] - nodes[panel_nodes[0]];
            vector3d AC = nodes[panel_nodes[2]] - nodes[panel_nodes[0]];
            normal = AB.cross(AC);

        } else { // 4 sides
            vector3d AC = nodes[panel_nodes[2]] - nodes[panel_nodes[0]];
            vector3d BD = nodes[panel_nodes[3]] - nodes[panel_nodes[1]];
            normal = BD.cross(AC);
        }

        // panel area
        panel_areas.push_back(normal.norm()/2.0);
        //cout << panel_areas[p] << endl;

        if(panel_areas[p] < Parameters::inversion_tolerance){
            normal = 0.0;
        }else{
            normal.normalize();
        }
        panel_normals.push_back(normal);
        //cout << normal << endl;
    }



    //compute collocation point
    panel_collocation_points[0].clear();
    panel_collocation_points[1].clear();
    panel_collocation_points[0].resize(n_panels());
    panel_collocation_points[1].resize(n_panels());

    for(int p = 0; p < n_panels(); p++){
        vector3d new_cp(0.0,0.0,0.0);
        for(int n = 0; n < (int)panels[p].size(); n++)
            new_cp = new_cp + nodes[panels[p][n]];

        new_cp = new_cp / (double)panels[p].size();
        panel_collocation_points[0][p] = new_cp;
        //cout << new_cp << endl;

        new_cp = new_cp - panel_normals[p] * Parameters::collocation_point_delta;
        panel_collocation_points[1][p] = new_cp;
        //cout << new_cp << endl;
    }

    // compute transformation
    for(int p = 0; p < n_panels(); p++){

        vector3d longitudinal = nodes[panels[p][0]] - nodes[panels[p][1]];
        if(longitudinal.norm() < Parameters::inversion_tolerance){
            longitudinal = 0.0;
        }else{
            longitudinal.normalize();
        }
        panel_longitudinals.push_back(longitudinal);
        //cout << longitudinal << endl;

        vector3d transverse = panel_normals[p].cross(panel_longitudinals[p]);
        panel_transverse.push_back(transverse);
        //cout << transverse << endl;
    }

    // compute panel local coordinates
    panel_local_coordinates.clear();
    panel_local_coordinates.resize(n_panels());
    for(int p = 0; p < n_panels(); p++){
        for(int n = 0; n < (int)panels[p].size(); n++){
            vector3d transformed_point = transform_point_panel(p,nodes[panels[p][n]]);
            panel_local_coordinates[p].push_back(transformed_point);
            //cout << transformed_point << endl;
        }
        //cout << endl;
    }


    // compute panel farfield distance
    panel_farfield_distance.clear();
    panel_farfield_distance.resize(n_panels());
    for(int p = 0; p < n_panels(); p++){
        vector<int> &panel_nodes = panels[p];
        double d1,d2;
        if (panel_nodes.size() == 3) {
            vector3d AB = nodes[panel_nodes[1]] - nodes[panel_nodes[0]];
            vector3d AC = nodes[panel_nodes[2]] - nodes[panel_nodes[0]];
            d1 = AB.norm();
            d2 = AC.norm();

        } else { // 4 sides
            vector3d AC = nodes[panel_nodes[2]] - nodes[panel_nodes[0]];
            vector3d BD = nodes[panel_nodes[3]] - nodes[panel_nodes[1]];
            d1 = AC.norm();
            d2 = BD.norm();
        }
        panel_farfield_distance[p] = max(d1,d2) * Parameters::farfield_factor;
        //cout << panel_farfield_distance[p] << endl;
    }

}

vector3d& Surface :: get_collocation_point(int panel,bool below_surface) {
    assert(panel < (int)panels.size());
    return panel_collocation_points[below_surface][panel];
}

vector3d Surface :: get_collocation_point(int panel,bool below_surface) const{
    assert(panel < (int)panels.size());
    return panel_collocation_points[below_surface][panel];
}


vector3d Surface :: transform_point_panel(int panel, const vector3d& x){
    vector3d transformed_point, diff;

    vector3d &l = panel_longitudinals[panel];
    vector3d &n = panel_normals[panel];
    vector3d &t = panel_transverse[panel];

    vector3d &cp = get_collocation_point(panel,false);

    diff = x - cp ;

    transformed_point[0] = diff[0]*l[0] + diff[1]*l[1] + diff[2]*l[2];
    transformed_point[1] = diff[0]*t[0] + diff[1]*t[1] + diff[2]*t[2];
    transformed_point[2] = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];

    return transformed_point;

}


void Surface :: translate_surface(const vector3d& dX){

    for(size_t n = 0; n < nodes.size(); n++)
        nodes[n] = nodes[n] + dX;
}


void Surface :: rotate_surface(vector3d dTheta, bool isRadian){


    if(!isRadian){
        dTheta = dTheta * (M_PI/180.0);
    }

    for(size_t n = 0; n < nodes.size(); n++){

        vector3d old_x = nodes[n];

        nodes[n][0] = cos(dTheta[1])*cos(dTheta[2])*old_x[0]
                    + cos(dTheta[1])*sin(dTheta[2])*old_x[1]
                    - sin(dTheta[1])*old_x[2];

        nodes[n][1] = (sin(dTheta[0])*sin(dTheta[1])*cos(dTheta[2]) - cos(dTheta[0])*sin(dTheta[2]))*old_x[0]
                    + (sin(dTheta[0])*sin(dTheta[1])*sin(dTheta[2]) + cos(dTheta[0])*cos(dTheta[2]))*old_x[1]
                    + sin(dTheta[0])*cos(dTheta[1])*old_x[2];

        nodes[n][2] = (cos(dTheta[0])*sin(dTheta[1])*cos(dTheta[2]) + sin(dTheta[0])*sin(dTheta[2]))*old_x[0]
                    + (cos(dTheta[0])*sin(dTheta[1])*sin(dTheta[2]) - sin(dTheta[0])*cos(dTheta[2]))*old_x[1]
                    + cos(dTheta[0])*cos(dTheta[1])*old_x[2];
    }
}


void Surface :: set_linear_velocity(const vector3d& vel){
    linear_velocity = vel;
}

void Surface :: set_angular_velocity(const vector3d& vel){
    angular_velocity = vel;
}


int Surface :: n_trailing_edge_nodes() const {
    assert(trailing_edge_nodes.size() > 0);
    return static_cast<int>(trailing_edge_nodes.size());
}

int Surface :: n_trailing_edge_panels() const {
    assert(upper_TE_panels.size() > 0);
    return static_cast<int>(upper_TE_panels.size());
}


vector3d Surface :: get_trailing_edge_bisector(const int TE_node) const {

    assert(TE_node < n_trailing_edge_nodes());
    assert(panel_longitudinals.size() > 0);

    int panel;

    if(TE_node == n_trailing_edge_nodes() - 1){
        panel = TE_node - 1;
    }else{
        panel = TE_node;
    }

    //cout << upper_TE_panels[panel] << endl;

    //cout << panel_longitudinals[0] << endl;

    vector3d  vec1 = - panel_longitudinals[lower_TE_panels[panel]];
    vector3d  vec2 =   panel_longitudinals[upper_TE_panels[panel]];

    vector3d bisector = vec1 + vec2 ;

    bisector.normalize();

    return bisector;
}


vector3d Surface :: get_kinematic_velocity(const vector3d& x) const {

    vector3d r = x - surface_origin;
    return linear_velocity + angular_velocity.cross(r);
}
