#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

/** @brief Defines certain parameters used in the panel method */

class Parameters {

public:

    /** @brief tolerance check for division */
    static double inversion_tolerance;

    /** @brief farfield factor */
    static double farfield_factor;

    /** @brief controls the trailing edge wake panel distance */
    static double trailing_edge_wake_shed_factor;

    /** @brief decides whether problem is steady or unsteady */
    static bool unsteady_problem;

    /** @brief decides whether to use static wake or force free wake model */
    static bool static_wake;

    /** @brief sets static wake length */
    static double static_wake_length;

    /** @brief decides whether to use vortex core model */
    static bool use_vortex_core_model;

};



#endif // PARAMETERS

