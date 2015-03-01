#include "parameters.hpp"

double Parameters :: collocation_point_delta = 1e-12;

double Parameters :: inversion_tolerance = 1e-12;

double Parameters :: farfield_factor = 15.0;

double Parameters :: trailing_edge_wake_shed_factor = 0.25;

bool Parameters :: wake_shed_along_TE_bisector = false;
