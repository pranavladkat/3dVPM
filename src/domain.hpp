#ifndef DOMAIN_H
#define DOMAIN_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <iomanip>

#include "vector3d.h"

class Domain
{
public:
    Domain();
    ~Domain();

    std::vector<vector3d> nodes;

};

#endif // DOMAIN_H
