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
private:

    int IMAX, JMAX, KMAX;

public:
    Domain();
    ~Domain();

    void set_domain_ranges(const int i,const int j,const int k);

    std::vector<vector3d> nodes;

    int get_IMAX() const ;

    int get_JMAX() const ;

    int get_KMAX() const ;

    int n_nodes() const ;

};

#endif // DOMAIN_H
