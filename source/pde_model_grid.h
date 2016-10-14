#include "pde_refinement.h"

template<>
void Model<1>::create_coarse_grid()
{
    const unsigned int dim = 1;
    // Create grid
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids, this->manifold_descriptors,
        params.geometry.grid_name, params.geometry.sizes);
    // Shift and rotate the grid.
    Point<dim> shifted_center;
    for (unsigned int i = 0; i < dim; i++)
    {
        shifted_center[i] = params.geometry.transformations[i];
    }
    GridTools::shift(shifted_center, this->triangulation); 
    this->spherical_manifold_center = shifted_center;
}

template <>
void Model<2>::create_coarse_grid()
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

template<>
void Model<3>::create_coarse_grid()
{
    Assert(false, ExcNotImplemented()); // Only missing a 3D rotation method
}

template<int dim>
void Model<dim>::adaptive_refine()
{
    SolutionTransfer<dim> solution_trans(dof_handler);
    Vector<double> previous_solution;
    previous_solution = solution;
    Refinement::adaptive_refine_mesh(
        triangulation,
        dof_handler,
        solution,
        solution_trans,
        fe,
        params.refinement.initial_global_cycles + params.refinement.initial_boundary_cycles,
        params.refinement.adaptive.max_level,
        params.refinement.adaptive.max_cells,
        params.refinement.adaptive.refine_fraction,
        params.refinement.adaptive.coarsen_fraction);
    setup_system();
    solution_trans.interpolate(previous_solution, solution);
    constraints.distribute(solution);
}
  