#include "domain.hpp"

Domain::Domain()
{

}

Domain::~Domain()
{

}

void Domain :: set_domain_ranges(const int i,const int j,const int k){
    IMAX = i;
    JMAX = j;
    KMAX = k;
}


int Domain :: get_IMAX() const{
    return IMAX;
}

int Domain :: get_JMAX() const{
    return JMAX;
}

int Domain :: get_KMAX() const{
    return KMAX;
}
