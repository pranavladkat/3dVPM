#ifndef SURFACE_H
#define SURFACE_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#include "vector3d.h"
#include "parameters.hpp"

/**
   Class for containing geometry entities and functions to perform operations on them

   @brief Represents geometry - takes input from mesh objects.
*/

class Surface
{
public:
    Surface();
    ~Surface();

    /** @brief stores node data containing x,y,z coordinate */
    std::vector<vector3d> nodes;

    /** @brief panel vector containing vector of its nodes
     *
     * represents panel vector. Each panel contains node number associated with it
     * the node numbers are in counter-clockwise order
     */
    std::vector<std::vector<int>> panels;

    /** @brief stores the neighbouring panel ids
     *
     * stores neighbouring panel ids - required in surface velocity calculations
     */
    std::vector<std::vector<int>> panel_neighbours;

    /** @brief nodes on trailing edge */
    std::vector<int> trailing_edge_nodes;

    /** @brief panels connected to upper trailing edge */
    std::vector<int> upper_TE_panels;

    /** @brief panels connected to lower trailing edge */
    std::vector<int> lower_TE_panels;

    /** @brief returns number of panels in the current surface
     *  @param[out] int     total number of panels in the surface */
    int n_panels() const;

    /** @brief returns total number of nodes in the current surface
     *  @param[out] int     total number of nodes in the surface */
    int n_nodes() const;

    /** @brief Computes the panel components
      *
      * computes panel collocation point, panel local coordinates, area,
      * panel transformations and farfield factors
      */
    void compute_panel_components();

    /** @brief returns panel's collocation point by reference
     *
     *  @param[in]  panel       panel number
     *  @param[out] vector3d&   panel collocation point by reference */
    vector3d& get_collocation_point(int panel);

    /** @brief returns panel's collocation point
     *
     *  @param[in]  panel      panel number
     *  @param[out] vector3d   panel collocation point */
    vector3d get_collocation_point(int panel) const;

    /** @brief translates the surface by distance dX
     *
     * @param[in] dX  distance vector */
    void translate_surface(const vector3d& dX);

    /** @brief rotates the surface by an angle
     *
     * @param[in] dTheta    angle by which to rorate the surface
     * @param[in] isRadian  true if \b dTheta is in radians, false otherwise */
    void rotate_surface(vector3d dTheta, bool isRadian);

    /** @brief Set the linear velocity of the surface */
    void set_linear_velocity(const vector3d&);

    /** @brief Set angular velocity of the body
     *
     * @param[in] vel  angular  velocity
     * @param[in] isRadian_sec  true if \b vel is in radians per seconds */
    void set_angular_velocity(vector3d vel, bool isRadian_sec);

    /** @brief returns number of trailing edge nodes */
    int n_trailing_edge_nodes() const ;

    /** @brief returns nuber of trailing edge panels */
    int n_trailing_edge_panels() const ;

    /** @brief returns trailing edge bisecor
     *
     * @param[in]  TE_node  trailing edge node
     * @param[out] vector3d  normalized vector which points in the bisector at \b TE_node */
    vector3d get_trailing_edge_bisector(const int) const ;

    /** @brief returns kinematic velocity at a point on surface
     *
     * @param[in]  x input location
     * @param[out] vector3d  kinematic velocity due to surface's motion */
    vector3d get_kinematic_velocity(const vector3d&) const;

    /** @brief return panel's normal vector */
    vector3d get_panel_normal(const int) const;

    /** @brief computes the influence coefficient due to source panel
     *
     * @param[in] panel  source panel whose influence is sought
     * @param[in] node   point which is being influenced by \b panel
     * @param[in] double         return in the influence due to unit source strength panel */
    double compute_source_panel_influence(const int panel, const vector3d& node) const;

    /** @brief computes the influence coefficient due to doublet panel
     *
     * @param[in] panel  doublet panel whose influence is sought
     * @param[in] node   point which is being influenced by \b panel
     * @param[in] double         return influence due to unit doublet strength panel */
    double compute_doublet_panel_influence(const int panel, const vector3d& node) const;

    /** @brief computes the influence coefficient due to source and doublet panel
     *
     * @param[in] panel  panel whose influence is sought
     * @param[in] node   point which is being influenced by \b panel
     * @param[in] pair<double,double>  return influence due to unit strength source and doublet panel */
    std::pair<double,double> compute_source_doublet_panel_influence(const int panel, const vector3d& node) const;

    /** @brief transforms a point in a panel's local coordinate */
    vector3d transform_point_panel(int panel, const vector3d& x) const;

    /** @brief transforms a point from panel's local coordinate to global coordinate */
    vector3d transform_vector_panel_inverse(int panel, const vector3d& x) const;

    /** @brief transforms a vector in a panel's local coordinate */
    vector3d transform_vector_panel(int panel, const vector3d& x) const;

    /** @brief transforms a vector from panel's local coordinate to global coordinate */
    double get_panel_area(const int& panel) const;

    /** @brief computes the induced velocity due to source panel at a point
     *
     * @param[in] panel  source panel
     * @param[in] node   point at which induced velocity by \b panel is sought
     * @param[in] vector3d    return induced velocity vector */
    vector3d compute_source_panel_unit_velocity(const int& panel, const vector3d& node) const;

    /** @brief computes the induced velocity due to doublet panel at a point
     *
     * @param[in] panel  doublet panel
     * @param[in] node   point at which induced velocity by \b panel is sought
     * @param[in] vector3d    return induced velocity vector */
    vector3d compute_doublet_panel_unit_velocity(const int& panel, const vector3d& node) const;

private:

    /** @brief represents " 1 / (4*pi)" */
    const double fourpi;

    /** @brief collocation points of the panels
     *
     * Stores the collocation point for each panel */
    std::vector<vector3d> panel_collocation_points;

    /** @brief Stores panel normal vector of each panel */
    std::vector<vector3d> panel_normals;

    /** @brief Stores panel longitudinal vector of each panel */
    std::vector<vector3d> panel_longitudinals;

    /** @brief Stores panel transverse vector of each panel */
    std::vector<vector3d> panel_transverse;

    /** @brief Stores area of each panel */
    std::vector<double> panel_areas;

    /** @brief Stores local coordinates of panel's nodes
     *
     * coordinates of nodes forming a panel are transformed to panel's local coordinates */
    std::vector<std::vector<vector3d>> panel_local_coordinates;

    /** @brief Stores farfield distance of each panel */
    std::vector<double> panel_farfield_distance;

    /** @brief Linear velocity of the surface, in m/s */
    vector3d linear_velocity;

    /** @brief Angular velocity of the surface, in rad/s */
    vector3d angular_velocity;

    /** @brief stores surface origin */
    vector3d surface_origin, previous_surface_origin;

    /** @brief saves surface orientation, needed for rotation */
    vector3d surface_orientation, previous_surface_orientation;

    /** @brief computes influence due to source edge at a point
     *
     * @param[in] node_a  one node of the edge
     * @param[in] node_b  other node of the edge
     * @param[in] vector3d    returns influence coefficient due to panel edge */
    double compute_source_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;

    /** @brief computes influence due to doublet edge at a point
     *
     * @param[in] node_a  one node of the edge
     * @param[in] node_b  other node of the edge
     * @param[in] vector3d    returns influence coefficient due to panel edge */
    double compute_doublet_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;

    /** @brief computes influence due to source and doublet edge at a point
     *
     * @param[in] node_a  one node of the edge
     * @param[in] node_b  other node of the edge
     * @param[in] vector3d    returns influence coefficient due to panel edge */
    std::pair<double,double> compute_source_doublet_edge_influence(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;

    /** @brief computes induced velocity due to source edge at a point
     *
     * @param[in] node_a  one node of the edge
     * @param[in] node_b  other node of the edge
     * @param[in] vector3d    return induced velocity vector due to panel edge */
    vector3d compute_source_panel_edge_unit_velocity(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;

    /** @brief computes induced velocity due to doublet edge at a point
     *
     * @param[in] node_a  one node of the edge
     * @param[in] node_b  other node of the edge
     * @param[in] vector3d    return induced velocity vector due to panel edge */
    vector3d compute_doublet_panel_edge_unit_velocity(const vector3d& node_a,const vector3d& node_b,const vector3d& x) const;


};

#endif // SURFACE_H
