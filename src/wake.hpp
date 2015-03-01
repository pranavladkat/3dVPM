#ifndef WAKE_H
#define WAKE_H

#include "surface.hpp"

class Wake : public Surface
{
private:

    std::shared_ptr<Surface> lifting_surface;

    std::vector<double> doublet_strength;

public:
    Wake();
    ~Wake();

    void add_lifting_surface(const std::shared_ptr<Surface> surf);

    void initialize(const vector3d& free_stream_velocity, double dt);

    void build_topology();

};

#endif // WAKE_H
