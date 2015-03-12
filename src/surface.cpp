#include "surface.hpp"

using namespace std;

Surface::Surface()
    :fourpi(0.25 * M_1_PI)
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

int Surface :: n_nodes() const{
    return static_cast<int> (nodes.size());
}


void Surface :: compute_panel_components(){

    assert(nodes.size() > 0 && panels.size() > 0);

    //compute panel normal and areas
    panel_areas.clear();
    panel_areas.resize(n_panels());
    panel_normals.clear();
    panel_normals.resize(n_panels());
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
        panel_areas[p] = normal.norm() / 2.0;
        //cout << panel_areas[p] << endl;

        normal.normalize();

        panel_normals[p] = normal;
        //cout << normal << endl;
    }



    //compute collocation point
    panel_collocation_points[0].clear();
    panel_collocation_points[1].clear();
    panel_collocation_points[0].resize(n_panels());
    panel_collocation_points[1].resize(n_panels());

    for(int p = 0; p < n_panels(); p++){
        vector3d new_cp(0,0,0);
        for(int n = 0; n < (int)panels[p].size(); n++)
            new_cp += nodes[panels[p][n]];

        new_cp = new_cp / (double)panels[p].size();
        panel_collocation_points[0][p] = new_cp;
        //cout << new_cp << endl;

        new_cp -= panel_normals[p] * Parameters::collocation_point_delta;
        panel_collocation_points[1][p] = new_cp;
        //cout << new_cp << endl;
    }

    // compute transformation
    panel_longitudinals.clear();
    panel_longitudinals.resize(n_panels());
    panel_transverse.clear();
    panel_transverse.resize(n_panels());
    for(int p = 0; p < n_panels(); p++){

        vector3d longitudinal = nodes[panels[p][0]] - nodes[panels[p][1]];

        longitudinal.normalize();

        panel_longitudinals[p] = longitudinal;
        //cout << longitudinal << endl;

        vector3d transverse = panel_normals[p].cross(panel_longitudinals[p]);
        panel_transverse[p] = transverse;
        //cout << transverse << endl;
    }

    // compute panel local coordinates
    panel_local_coordinates.clear();
    panel_local_coordinates.resize(n_panels());
    for(int p = 0; p < n_panels(); p++){

        panel_local_coordinates[p].clear();
        panel_local_coordinates[p].resize(panels[p].size());

        for(int n = 0; n < (int)panels[p].size(); n++){
            panel_local_coordinates[p][n] = transform_point_panel(p,nodes[panels[p][n]]);
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


vector3d Surface :: transform_point_panel(int panel, const vector3d& x) const{
    vector3d transformed_point, diff;

    const vector3d &l = panel_longitudinals[panel];
    const vector3d &n = panel_normals[panel];
    const vector3d &t = panel_transverse[panel];

    const vector3d &cp = get_collocation_point(panel,false);

    diff = x - cp ;

    transformed_point[0] = diff[0]*l[0] + diff[1]*l[1] + diff[2]*l[2];
    transformed_point[1] = diff[0]*t[0] + diff[1]*t[1] + diff[2]*t[2];
    transformed_point[2] = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];

    return transformed_point;

}


void Surface :: translate_surface(const vector3d& dX){

    for(size_t n = 0; n < nodes.size(); n++)
        nodes[n] += dX;

    // update surface origin
    previous_surface_origin = surface_origin;
    surface_origin += dX;

    // update panel components
    compute_panel_components();
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

    // update surface orientation
    previous_surface_orientation = surface_orientation;
    surface_orientation += dTheta;

    // update panel components
    compute_panel_components();
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

vector3d Surface :: get_panel_normal(const int i) const{
    return panel_normals[i];
}


double Surface :: compute_source_panel_influence(const int panel, const vector3d& node) const{

    vector3d transformed_node = transform_point_panel(panel,node);

    double distance = transformed_node.norm();

    if(distance > panel_farfield_distance[panel])
        return (fourpi * panel_areas[panel] / distance);

    double influence = 0.0;
    for(size_t n = 0; n < panels[panel].size(); n++){

        int next_node = n + 1;
        if(n == panels[panel].size() - 1)
            next_node = 0;

        const vector3d& node_a = panel_local_coordinates[panel][n];
        const vector3d& node_b = panel_local_coordinates[panel][next_node];

        influence += compute_source_edge_influence(node_a, node_b,transformed_node);
    }

    return -influence*fourpi;
}


double Surface :: compute_source_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const{

    double r1 = (x-node_a).norm();
    double r2 = (x-node_b).norm();
    double d12 = (node_b - node_a).norm();

    double influence = 0;

    if(d12 < Parameters::inversion_tolerance || (r1+r2-d12) < Parameters::inversion_tolerance){
        influence = 0;
    }else{
        influence = ((x[0]-node_a[0])*(node_b[1] - node_a[1]) - (x[1]-node_a[1])*(node_b[0]-node_a[0]))
                  / d12 * log((r1+r2+d12) / (r1+r2-d12));
    }

    if(fabs(x[2]) > Parameters::inversion_tolerance){

        double e1 = pow((x[0] - node_a[0]),2) + pow(x[2],2);
        double e2 = pow((x[0] - node_b[0]),2) + pow(x[2],2);
        double h1 = (x[0] -node_a[0])*(x[1] - node_a[1]);
        double h2 = (x[0]- node_b[0])*(x[1] - node_b[1]);
        double m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);

        double F = (m*e1 - h1) / (x[2]*r1) ;
        double G = (m*e2 - h2) / (x[2]*r2) ;

        if(F != G)
            influence -= x[2] * atan2(F-G, 1+F*G);
    }

    return influence;
}



double Surface :: compute_doublet_panel_influence(const int panel, const vector3d& node) const {

    vector3d transformed_node = transform_point_panel(panel,node);

    double distance = transformed_node.norm();

    if(distance > panel_farfield_distance[panel])
        return (fourpi * panel_areas[panel] * transformed_node[2] * pow(distance,-3.0));

    double influence = 0.0;
    for(size_t n = 0; n < panels[panel].size(); n++){

        int next_node = n + 1;
        if(n == panels[panel].size() - 1)
            next_node = 0;

        const vector3d& node_a = panel_local_coordinates[panel][n];
        const vector3d& node_b = panel_local_coordinates[panel][next_node];

        influence += compute_doublet_edge_influence(node_a, node_b,transformed_node);
    }

    return influence*fourpi;
}


double Surface :: compute_doublet_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const{

    double influence = 0;

    double r1 = (x-node_a).norm();
    double r2 = (x-node_b).norm();
    double e1 = pow((x[0] - node_a[0]),2) + pow(x[2],2);
    double e2 = pow((x[0] - node_b[0]),2) + pow(x[2],2);
    double h1 = (x[0] -node_a[0])*(x[1] - node_a[1]);
    double h2 = (x[0]- node_b[0])*(x[1] - node_b[1]);
    double m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);

    double F = (m*e1 - h1) / (x[2]*r1) ;
    double G = (m*e2 - h2) / (x[2]*r2) ;

    if(F != G)
        influence = atan2(F-G, 1+F*G);

    return influence;
}


std::pair<double,double> Surface :: compute_source_doublet_panel_influence(const int panel, const vector3d& node) const {

    //first member is source influence, second is doublet influence
    pair<double,double> influence_coefficient;

    vector3d transformed_node = transform_point_panel(panel,node);

    double distance = transformed_node.norm();

    if(distance > panel_farfield_distance[panel]){
        influence_coefficient.first  = (fourpi * panel_areas[panel] / distance);
        influence_coefficient.second = (fourpi * panel_areas[panel] * transformed_node[2] * pow(distance,-3.0));
        return influence_coefficient;
    }

    for(size_t n = 0; n < panels[panel].size(); n++){

        int next_node = n + 1;
        if(n == panels[panel].size() - 1)
            next_node = 0;

        const vector3d& node_a = panel_local_coordinates[panel][n];
        const vector3d& node_b = panel_local_coordinates[panel][next_node];

        pair<double,double> edge_influence = compute_source_doublet_edge_influence(node_a, node_b,transformed_node);

        influence_coefficient.first  += edge_influence.first;
        influence_coefficient.second += edge_influence.second;
    }

    influence_coefficient.first  *= -fourpi;
    influence_coefficient.second *=  fourpi;

    return influence_coefficient;
}


std::pair<double,double> Surface :: compute_source_doublet_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const {

    pair<double,double> edge_influence(0,0);

    double r1 = (x-node_a).norm();
    double r2 = (x-node_b).norm();
    double d12 = (node_b - node_a).norm();

    if(d12 < Parameters::inversion_tolerance || (r1+r2-d12) < Parameters::inversion_tolerance){
        edge_influence.first = 0;
    }else{
        edge_influence.first = ((x[0]-node_a[0])*(node_b[1] - node_a[1]) - (x[1]-node_a[1])*(node_b[0]-node_a[0]))
                             / d12 * log((r1+r2+d12) / (r1+r2-d12));
    }

    if(fabs(x[2]) > Parameters::inversion_tolerance){

        double e1 = pow((x[0] - node_a[0]),2) + pow(x[2],2);
        double e2 = pow((x[0] - node_b[0]),2) + pow(x[2],2);
        double h1 = (x[0] -node_a[0])*(x[1] - node_a[1]);
        double h2 = (x[0]- node_b[0])*(x[1] - node_b[1]);
        double m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);

        double F = (m*e1 - h1) / (x[2]*r1) ;
        double G = (m*e2 - h2) / (x[2]*r2) ;

        double doublet_coef = 0;

        if(F != G)
            doublet_coef = atan2(F-G, 1+F*G);

        edge_influence.first  -= x[2] * doublet_coef;
        edge_influence.second  = doublet_coef;
    }

    return edge_influence;
}


vector3d Surface :: transform_vector_panel_inverse(int panel, const vector3d& x) const{

    vector3d transformed_vector;

    const vector3d &l = panel_longitudinals[panel];
    const vector3d &n = panel_normals[panel];
    const vector3d &t = panel_transverse[panel];

    transformed_vector[0] = x[0]*l[0] + x[1]*t[0] + x[2]*n[0];
    transformed_vector[1] = x[0]*l[1] + x[1]*t[1] + x[2]*n[1];
    transformed_vector[2] = x[0]*l[2] + x[1]*t[2] + x[2]*n[2];

    return transformed_vector;
}

vector3d Surface :: transform_vector_panel(int panel, const vector3d& x) const{

    vector3d transformed_vector;

    const vector3d &l = panel_longitudinals[panel];
    const vector3d &n = panel_normals[panel];
    const vector3d &t = panel_transverse[panel];

    transformed_vector[0] = x[0]*l[0] + x[1]*l[1] + x[2]*l[2];
    transformed_vector[1] = x[0]*t[0] + x[1]*t[1] + x[2]*t[2];
    transformed_vector[2] = x[0]*n[0] + x[1]*n[1] + x[2]*n[2];

    return transformed_vector;
}


double Surface :: get_panel_area(const int& panel) const{
    return panel_areas[panel];
}


vector3d Surface :: compute_source_panel_unit_velocity(const int& panel, const vector3d& node) const{

    vector3d panel_velocity(0,0,0);
    // transform node in panel coordinates
    vector3d transformed_node = transform_point_panel(panel,node);

    // if z -> +-0 then velocity = (0,0,0) if x lies outside of panel
    if(fabs(transformed_node[2]) < Parameters::inversion_tolerance)
        return panel_velocity;

    double distance = transformed_node.norm();

    // farfield formulation
    if(distance > panel_farfield_distance[panel]){

        panel_velocity[0] = panel_areas[panel] * transformed_node[0] * fourpi * pow(distance,-3.0);
        panel_velocity[1] = panel_areas[panel] * transformed_node[1] * fourpi * pow(distance,-3.0);
        panel_velocity[2] = panel_areas[panel] * transformed_node[2] * fourpi * pow(distance,-3.0);

        // transform velocity back in global coordinates
        panel_velocity = transform_vector_panel_inverse(panel,panel_velocity);
        return panel_velocity;
    }

    // calculate velocity due to panel edges
    for(size_t n = 0; n < panels[panel].size(); n++){

        int next_node = n + 1;
        if(n == panels[panel].size() - 1)
            next_node = 0;

        const vector3d& node_a = panel_local_coordinates[panel][n];
        const vector3d& node_b = panel_local_coordinates[panel][next_node];

        panel_velocity = panel_velocity + compute_source_panel_edge_unit_velocity(node_a, node_b,transformed_node);
    }

    // transform velocity back in global coordinates
    panel_velocity = transform_vector_panel_inverse(panel,panel_velocity);

    return panel_velocity * fourpi;
}


vector3d Surface :: compute_source_panel_edge_unit_velocity(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const{

    vector3d edge_velocity(0,0,0);

    double r1 = (x-node_a).norm();
    double r2 = (x-node_b).norm();
    double d12 = (node_b - node_a).norm();

    double e1 = pow((x[0] - node_a[0]),2) + pow(x[2],2);
    double e2 = pow((x[0] - node_b[0]),2) + pow(x[2],2);
    double h1 = (x[0] -node_a[0])*(x[1] - node_a[1]);
    double h2 = (x[0]- node_b[0])*(x[1] - node_b[1]);
    double m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);

    if(d12 > Parameters::inversion_tolerance && (r1+r2-d12) > Parameters::inversion_tolerance){
        edge_velocity[0] = (node_b[1] - node_a[1]) / d12 * log((r1+r2-d12)/(r1+r2+d12));
        edge_velocity[1] = (node_a[0] - node_b[0]) / d12 * log((r1+r2-d12)/(r1+r2+d12));
    }

    double F = (m*e1 - h1) / (x[2]*r1) ;
    double G = (m*e2 - h2) / (x[2]*r2) ;

    if(F != G)
        edge_velocity[2] = atan2(F-G, 1+F*G);

    return edge_velocity;
}



vector3d Surface :: compute_doublet_panel_unit_velocity(const int& panel, const vector3d& node) const{

    vector3d panel_velocity(0,0,0);

    // transform node in panel coordinates
    //vector3d transformed_node = transform_point_panel(panel,node);

    //double distance = transformed_node.norm();

    // farfield formulation
//    if(distance > panel_farfield_distance[panel]){

//        panel_velocity[0] = 3 * panel_areas[panel] * transformed_node[0] * transformed_node[2] * fourpi * pow(distance,-5.0);
//        panel_velocity[1] = 3 * panel_areas[panel] * transformed_node[1] * transformed_node[2] * fourpi * pow(distance,-5.0);
//        panel_velocity[2] = - panel_areas[panel] * (pow(transformed_node[0],2.0) + pow(transformed_node[1],2.0)
//                          - 2 * pow(transformed_node[2],2.0)) * fourpi * pow(distance,-5.0);
//        // transform to global coordinates
//        panel_velocity = transform_vector_panel_inverse(panel,panel_velocity);
//        return panel_velocity;
//    }

    // calculate velocity due to panel edges using Vortex Line formulation
    // node does not need to be transformed to panel local coordinates
    for(size_t n = 0; n < panels[panel].size(); n++){

        int previous_node = n - 1;
        if(n == 0)
            previous_node = panels[panel].size() - 1;

        const vector3d& node_a = nodes[panels[panel][previous_node]];
        const vector3d& node_b = nodes[panels[panel][n]];

        // send nodes in global coordinates (not in panel local coordinates)
        panel_velocity -= compute_doublet_panel_edge_unit_velocity(node_a, node_b, node);
    }

    return panel_velocity;
}


vector3d Surface :: compute_doublet_panel_edge_unit_velocity(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const{

    // computes velocity using vortex line method

    vector3d edge_velocity(0,0,0);

    double r1r2x, r1r2y, r1r2z, r1r2_sq;
    double r1,r2, r0r1,r0r2, coef;

    r1r2x =   (x[1] - node_a[1]) * (x[2] - node_b[2]) - (x[2] - node_a[2]) * (x[1] - node_b[1]);
    r1r2y = - (x[0] - node_a[0]) * (x[2] - node_b[2]) + (x[2] - node_a[2]) * (x[0] - node_b[0]);
    r1r2z =   (x[0] - node_a[0]) * (x[1] - node_b[1]) - (x[1] - node_a[1]) * (x[0] - node_b[0]);

    r1r2_sq = r1r2x*r1r2x + r1r2y*r1r2y + r1r2z*r1r2z;

    r1 = (x - node_a).norm();
    r2 = (x - node_b).norm();

    if( r1 > Parameters::inversion_tolerance && r2 > Parameters::inversion_tolerance && r1r2_sq > Parameters::inversion_tolerance ){

        r0r1 = (node_b[0] - node_a[0])*(x[0] - node_a[0]) + (node_b[1] - node_a[1])*(x[1] - node_a[1]) + (node_b[2] - node_a[2])*(x[2] - node_a[2]);
        r0r2 = (node_b[0] - node_a[0])*(x[0] - node_b[0]) + (node_b[1] - node_a[1])*(x[1] - node_b[1]) + (node_b[2] - node_a[2])*(x[2] - node_b[2]);

        double Kv = 1.0;
        if(Parameters::use_vortex_core_model){
            double vm = 2.0;    /* decides vortex core model */
            /* h = perpendicular dist of x from line joining x1 and x2,
             * more details: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html */
            double h = sqrt(r1r2_sq) / sqrt(pow((node_b[0]-node_a[0]),2) + pow((node_b[1]-node_a[1]),2) + pow((node_b[2]-node_a[2]),2));
            //rc = sqrt(4*alpha*deltav*visc*t);
            double rc = 0.1;

            /* Kv : parameter to disingularize biot savart law,
             * refer to: Estimating the Angle of Attack from Blade Pressure Measurements on the
             *           NREL Phase VI Rotor Using a Free Wake Vortex Model: Axial Conditions,
             *           Wind Energ. 2006; 9:549–577 */
            Kv = (h*h)/ pow((pow(rc,2*vm) + pow(h,2*vm)),1.0/vm);
        }

        coef = Kv*fourpi/(r1r2_sq)*(r0r1/r1 - r0r2/r2);

        edge_velocity[0] = r1r2x*coef;
        edge_velocity[1] = r1r2y*coef;
        edge_velocity[2] = r1r2z*coef;
    }

    return edge_velocity;
}




