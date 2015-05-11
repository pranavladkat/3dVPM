#ifndef WAKE_H
#define WAKE_H

#include "surface.hpp"

/** @brief Wake class to represent wake of a body  */

class Wake : public Surface
{
private:

    /** @brief pointer to a surface  */
    std::shared_ptr<Surface> lifting_surface;

    /** @brief doublet strength os the wake panels */
    std::vector<double> doublet_strength;

public:
    Wake();
    ~Wake();

    /** @brief associate with surface object  */
    void add_lifting_surface(const std::shared_ptr<Surface> surf);

    /** @brief create first row of wake panels */
    void initialize(const vector3d& free_stream_velocity, const double& dt);

    /** @brief build topology of the wake panels, such as connectivity information */
    void build_topology();

    /** @brief adds a row of wake panels */
    void shed_wake(const vector3d& free_stream_velocity, double dt);

};

#endif // WAKE_H
