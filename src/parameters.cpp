#include "parameters.hpp"

double Parameters :: collocation_point_delta = 1e-12;

double Parameters :: inversion_tolerance = 1e-12;

double Parameters :: farfield_factor = 10.0;

double Parameters :: trailing_edge_wake_shed_factor = 0.25;

bool Parameters :: unsteady_problem = false;

double Parameters :: static_wake_length = 1.0;

bool Parameters :: static_wake = false;

bool Parameters :: use_vortex_core_model = false;
