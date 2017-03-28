#include "refinement.h"

/*! Set the coarse Triangulation, i.e. Peclet::triangulation

    In deal.II language, coarse means the coarsest available representation of the geometry. In combination with geometric manifolds, this can indeed be quite coarse, even for curved geometries. For example, see the hemisphere_cylinder_shell in MyGridGenerator.
    
    Unfortunately geometric manifold data is not contained by the Triangulation. This separation allows for flexibility, but also requires extra bookkeeping by the programmer.

    Peclet<1>::create_coarse_grid(), Peclet<2>::create_coarse_grid(), and Peclet<3>::create_coarse_grid() had to be implemented separately because GridTools::rotate() has not be implemented for Triangulation<1>.
    
*/
template<>
void Peclet<1>::create_coarse_grid()
{
    const unsigned int dim = 1;
    // Create grid
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids, this->manifold_descriptors,
        this->params.geometry.grid_name, params.geometry.sizes);
    // Shift and rotate the grid.
    Point<dim> shifted_center;
    for (unsigned int i = 0; i < dim; i++)
    {
        shifted_center[i] = params.geometry.transformations[i];
    }
    GridTools::shift(shifted_center, this->triangulation); 
    this->spherical_manifold_center = shifted_center;
}

/*! Set the coarse Triangulation, i.e. Peclet::triangulation

    See Peclet<1>::create_coarse_grid().

*/
template <>
void Peclet<2>::create_coarse_grid()
{
    const unsigned int dim = 2;
    // Create grid
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids, this->manifold_descriptors,
        params.geometry.grid_name, params.geometry.sizes);
    // Shift and rotate the grid.
    Point<dim> shifted_center;
    GridTools::rotate(params.geometry.transformations[2], this->triangulation);
    for (unsigned int i = 0; i < dim; i++)
    {
        shifted_center[i] = params.geometry.transformations[i];
    }
    GridTools::shift(shifted_center, this->triangulation); 
    this->spherical_manifold_center = shifted_center;
}

/*! Set the coarse Triangulation, i.e. Peclet::triangulation

    See Peclet<1>::create_coarse_grid().

*/
template<>
void Peclet<3>::create_coarse_grid()
{
    Assert(false, ExcNotImplemented()); // Only missing a 3D rotation method
}

/*! Adaptively refine the triangulation based on an error measure.
            
    This is mostly a copy of the routine from deal.II's step-26, which uses the Kelly Error Estimator.
        
*/
template<int dim>
void Peclet<dim>::adaptive_refine()
{
    SolutionTransfer<dim> solution_trans(this->dof_handler);
    Vector<double> previous_solution;
    previous_solution = this->solution;
    Refinement::adaptive_refine_mesh(
        this->triangulation,
        this->dof_handler,
        this->solution,
        solution_trans,
        this->fe,
        this->params.refinement.initial_global_cycles + params.refinement.initial_boundary_cycles,
        this->params.refinement.adaptive.max_level,
        this->params.refinement.adaptive.max_cells,
        this->params.refinement.adaptive.refine_fraction,
        this->params.refinement.adaptive.coarsen_fraction);
    this->setup_system();
    solution_trans.interpolate(previous_solution, this->solution);
    this->constraints.distribute(this->solution);
}
  