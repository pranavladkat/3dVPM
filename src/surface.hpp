#ifndef SURFACE_H
#define SURFACE_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#include "vector3d.h"
#include "parameters.hpp"

class Surface
{
public:
    Surface();
    ~Surface();

    /** @brief Nodes containing x,y,z space coordinate */
    std::vector<vector3d> nodes;

    /** @brief panel vector containing vector of its nodes */
    std::vector<std::vector<int>> panels;

    /** @brief stores the neighbouring panel ids */
    std::vector<std::vector<int>> panel_neighbours;

    /** @brief nodes on trailing edge */
    std::vector<int> trailing_edge_nodes;

    /** @brief panels connected to upper trailing edge */
    std::vector<int> upper_TE_panels;

    /** @brief panels connected to lower trailing edge */
    std::vector<int> lower_TE_panels;

    //functions

    int n_panels() const;

    void compute_panel_components();

    vector3d& get_collocation_point(int panel,bool below_surface);

    vector3d get_collocation_point(int panel,bool below_surface) const;

    void translate_surface(const vector3d&);

    void rotate_surface(vector3d, bool);

    void set_linear_velocity(const vector3d&);

    void set_angular_velocity(const vector3d&);

    int n_trailing_edge_nodes() const ;

    int n_trailing_edge_panels() const ;

    vector3d get_trailing_edge_bisector(const int) const ;

    vector3d get_kinematic_velocity(const vector3d&) const;

    vector3d get_panel_normal(const int) const;

    double compute_source_panel_influence(const int panel, const vector3d& node) const;


private:

    const double fourpi;

    /** @brief collocation points of the panels */
    std::vector<vector3d> panel_collocation_points[2];

    /** @brief panel normal components */
    std::vector<vector3d> panel_normals;

    std::vector<vector3d> panel_longitudinals;

    std::vector<vector3d> panel_transverse;

    std::vector<double> panel_areas;

    std::vector<std::vector<vector3d>> panel_local_coordinates;

    std::vector<double> panel_farfield_distance;

    vector3d transform_point_panel(int panel, const vector3d& x) const;

    vector3d linear_velocity;

    vector3d angular_velocity;

    vector3d surface_origin, previous_surface_origin;

    vector3d surface_orientation, previous_surface_orientation;

    double compute_source_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;

};

#endif // SURFACE_H
