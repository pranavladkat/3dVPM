#ifndef DOMAIN_H
#define DOMAIN_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include <iomanip>

#include "vector3d.h"

/** @brief Domain class to compute velocity off-surface  */

class Domain
{
private:

    /** @brief maximum number of nodes in each direction  */
    int IMAX, JMAX, KMAX;

public:
    Domain();
    ~Domain();

    void set_domain_ranges(const int i,const int j,const int k);

    /** @brief point data of the domain mesh file */
    std::vector<vector3d> nodes;

    /** @brief returns maximum nodes in x direction */
    int get_IMAX() const ;

    /** @brief returns maximum nodes in y direction */
    int get_JMAX() const ;

    /** @brief returns maximum nodes in z direction */
    int get_KMAX() const ;

    /** @brief returns total number of nodes */
    int n_nodes() const ;

};

#endif // DOMAIN_H
