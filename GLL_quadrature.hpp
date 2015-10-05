#ifndef GLL_QUADRATURE_HPP
#define GLL_QUADRATURE_HPP

#include "mfem.hpp"

/**
 * Initialization of a Gauss-Lobatto-Legendre quadrature rule on a reference
 * segment [0, 1].
 * @param p - order
 * @param segment_GLL - integration rule to be set up
 */
void create_segment_GLL_rule(int p, mfem::IntegrationRule& segment_GLL);

#endif // GLL_QUADRATURE_HPP
