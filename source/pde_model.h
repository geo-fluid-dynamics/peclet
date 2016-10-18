 /*
 * @brief peclet solves the convection-diffusion equation
 *
 * @detail
 *
 *  This solves the dimensionless form of the unsteady convection-diffusion problem,
 *  with space and time dependent convection velocity, but constant diffusivity.
 *  The reference Peclet Number Pe_r = x_r*a_r/nu
 *
 *  Matrix assembly and time stepping are based on deal.II Tutorial 26 by Wolfgang Bangerth, Texas A&M University, 2013
 *
 *  Some of the more notable extensions include:
 *  - Builds convection-diffusion matrix instead of Laplace matrix.
 *  - Supports time-dependent non-zero Dirichlet and Neumann boundary condition.
 *  - Decomposed the classes and methods into multiple files.
 *  - Re-designed parmameter handling
 *  - Generalized boundary condition handling via the parameter input file
 *  - Writes FEFieldFunction to disk, and can read it from disk to initialize a restart
 *  - Added verification via MMS (Method of Manufactured Solutions) with error table based on approach from Tutorial 7
 *
 * @author Alexander Zimmerman <zimmerman@aices.rwth-aachen.de>, RWTH AAchen University, 2016
 */
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/table_handler.h>

#include <iostream>
#include <functional>

#include <assert.h> 
#include <deal.II/grid/manifold_lib.h>

#include "my_functions.h"
#include "extrapolated_field.h"
#include "pde_parameters.h"
#include "my_grid_generator.h"
#include "fe_field_tools.h"
#include "output.h"
#include "my_matrix_creator.h"
#include "my_vector_tools.h"


namespace PDE
{
  using namespace dealii;
  
  template<int dim>
  class Model
  {
  public:
    Model();
    Parameters::StructuredParameters params;
    void read_parameters(std::string file_path);
    void run();

  private:
    void create_coarse_grid();
    void adaptive_refine();
    void setup_system(bool quiet = false);
    void solve_time_step();
    void write_solution();
    
    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> convection_diffusion_matrix;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       old_solution;
    Vector<double>       system_rhs;

    double               time;
    unsigned int         time_step_counter;
    
    Point<dim> spherical_manifold_center;
    
    std::vector<unsigned int> manifold_ids;
    std::vector<std::string> manifold_descriptors;
    
    double reference_peclet_number;
    Function<dim>* convection_velocity_function;
    
    void mms_append_error_table();
    void mms_write_error_table();
    TableHandler mms_error_table;
    std::string mms_error_table_file_name = "mms_error_table.txt";
    
    void append_1D_solution_to_table();
    void write_1D_solution_table(std::string file_name);
    TableHandler solution_table_1D;
    std::string solution_table_1D_file_name = "1D_solution_table.txt";
    
  };

  template<int dim>
  Model<dim>::Model()
    :
    params(Parameters::get_parameters()),
    fe(1),
    dof_handler(triangulation)
  {}
  
  template<int dim>
  void Model<dim>::read_parameters(std::string file_path)
  {
    params = Parameters::get_parameters(file_path);
  }
  
  #include "pde_model_grid.h"
  
  template<int dim>
  void Model<dim>::setup_system(bool quiet)
  {
    dof_handler.distribute_dofs(fe);

    if (!quiet)
    {
        std::cout << std::endl
              << "==========================================="
              << std::endl
              << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;    
    }

    constraints.clear();
    
    DoFTools::make_hanging_node_constraints(
        dof_handler,
        constraints);
        
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    
    DoFTools::make_sparsity_pattern(
        dof_handler,
        dsp,
        constraints,
        /*keep_constrained_dofs = */ true);
        
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    
    convection_diffusion_matrix.reinit(sparsity_pattern);
    
    system_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree+1),
                                      mass_matrix);
                       
    /*
    In the unitless form of the convection-diffusion equation,
    the inverse of the Peclet Number replaces the momentum diffusivity (nu) 
    from the standard formulation.
    */
    ConstantFunction<dim> inverse_reference_peclet_number_function(1./this->reference_peclet_number);
    
    MyMatrixCreator::create_convection_diffusion_matrix(
        dof_handler,
        QGauss<dim>(fe.degree+1),
        convection_diffusion_matrix,
        &inverse_reference_peclet_number_function, 
        &this->convection_velocity_function
        );

    solution.reinit(dof_handler.n_dofs());
    
    old_solution.reinit(dof_handler.n_dofs());
    
    system_rhs.reinit(dof_handler.n_dofs());
    
  }

  template<int dim>
  void Model<dim>::solve_time_step()
  {
    double tolerance = params.solver.tolerance;
    if (params.solver.normalize_tolerance)
    {
        tolerance *= system_rhs.l2_norm();
    }
    SolverControl solver_control(
        params.solver.max_iterations,
        tolerance);
       
    SolverCG<> solver_cg(solver_control);
    SolverBicgstab<> solver_bicgstab(solver_control);

    PreconditionSSOR<> preconditioner;
    
    preconditioner.initialize(system_matrix, 1.0);

    std::string solver_name;
    
    if (params.solver.method == "CG")
    {
        solver_name = "CG";
        solver_cg.solve(
            system_matrix,
            solution,
            system_rhs,
            preconditioner);    
    }
    else if (params.solver.method == "BiCGStab")
    {
        solver_name = "BiCGStab";
        solver_bicgstab.solve(
            system_matrix,
            solution,
            system_rhs,
            preconditioner);
    }

    constraints.distribute(solution);

    std::cout << "     " << solver_control.last_step()
              << " " << solver_name << " iterations." << std::endl;
  }
  
  #include "pde_1D_solution_table.h"
  
  template<int dim>
  void Model<dim>::write_solution()
  {
      
    if (params.output.write_solution_vtk)
    {
        Output::write_solution_to_vtk(
            "solution-"+Utilities::int_to_string(time_step_counter)+".vtk",
            this->dof_handler,
            this->solution);    
    }
    
    if (dim == 1)
    {
        this->append_solution_to_table();
    }
    
  }
  
  template<int dim>
  void Model<dim>::mms_append_error_table()
  {
    assert(params.mms.enabled);
    
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    
    ManufacturedSolution<dim> manufactured_solution;
    manufactured_solution.set_time(time);
    
    VectorTools::integrate_difference()
        dof_handler,
        solution,
        manufactured_solution,
        difference_per_cell,
        QGauss<dim>(3),
        VectorTools::L2_norm);
        
    double L2_norm_error = difference_per_cell.l2_norm();
    
    mms_error_table.add_value("time_step_size", time_step_size);
    mms_error_table.add_value("time", time);
    mms_error_table.add_value("cells", triangulation.n_active_cells());
    mms_error_table.add_value("dofs", dof_handler.n_dofs());
    mms_error_table.add_value("L2_norm_error", L2_norm_error);
  }
  
  template<int dim>
  void Model<dim>::mms_write_error_table()
  {
    const int precision = 14;
    
    this->mms_error_table.set_precision("time", precision);
    this->mms_error_table.set_scientific("time", true);
    
    this->mms_error_table.set_precision("time_step_size", precision);
    this->mms_error_table.set_scientific("time_step_size", true);
    
    this->mms_error_table.set_precision("cells", precision);
    this->mms_error_table.set_scientific("cells", true);
    
    this->mms_error_table.set_precision("dofs", precision);
    this->mms_error_table.set_scientific("dofs", true);
    
    this->mms_error_table.set_precision("L2_norm_error", precision);
    this->mms_error_table.set_scientific("L2_norm_error", true);
    
    std::ofstream out_file(this->mms_error_table_file_name, std::fstream::app);
    assert(out_file.good());
    this->mms_error_table.write_text(out_file);
    out_file.close(); 
  }
  
  template<int dim>
  void Model<dim>::run()
  {
    if (dim == 1)
    {
        std::remove(solution_table_1D_file_name.c_str()); // In 1D, the solution will be appended here at every time step.    
    }        
    
    this->reference_peclet_number = params.pde.reference_peclet_number;
    
    MMS::ConvectionVelocity<dim> mms_convection_velocity_function;
    
    if (params.pde.convection_velocity_function_name == "MMS")
    {
        assert(params.mms.enabled);
        this->convection_velocity_function = &mms_convection_velocity_function;
    }
    else
    {
        Assert(false, ExcNotImplemented);
        // @todo: Implement constant and ramp; shouldn't be much work
    }
    
    create_coarse_grid();
    
    // Attach manifolds
    assert(dim < 3); // @todo: 3D extension: For now the CylindricalManifold is being ommitted.
        // deal.II makes is impractical for a CylindricalManifold to exist in 2D.
    SphericalManifold<dim> spherical_manifold(spherical_manifold_center);
    
    for (unsigned int i = 0; i < manifold_ids.size(); i++)
    {
        if (manifold_descriptors[i] == "spherical")
        {
            triangulation.set_manifold(manifold_ids[i], spherical_manifold);      
        }
    }
    
    //
    
    #include "pde_model_run_initialize_functions.h"
    
    // Initialize refinement
    
    triangulation.refine_global(params.refinement.initial_global_cycles);
    
    Refinement::refine_mesh_near_boundaries(
        triangulation,
        params.refinement.boundaries_to_refine,
        params.refinement.initial_boundary_cycles);
        
    // Initialize system
    setup_system();

    unsigned int pre_refinement_step = 0;
    Vector<double> tmp;
    Vector<double> forcing_terms;
    
    // Iterate
start_time_iteration:

    tmp.reinit (solution.size());

    VectorTools::interpolate(dof_handler,
                             *initial_values_function,
                             old_solution); 
    
    solution = old_solution;
    time_step_counter = 0; // @todo: Expose initial time step as a parameter
    time            = 0; // @todo: Expose initial time as a parameter
    write_solution();
    
    double epsilon = 1e-16;
    while (time < params.time.end_time - epsilon)
    {
        ++time_step_counter;
        time = params.time.step_size*time_step_counter; // Incrementing the time directly would accumulate errors
        
        std::cout << "Time step " << time_step_counter << " at t=" << time << std::endl;

        mass_matrix.vmult(system_rhs, old_solution);

        convection_diffusion_matrix.vmult(tmp, old_solution);
        
        system_rhs.add(-(1. - params.time.semi_implicit_theta) * params.time.step_size, tmp);
        
        // Add source terms
        source_function.set_time(time);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            source_function,
                                            tmp);
        forcing_terms = tmp;
        forcing_terms *= time_step_size * theta;
        
        source_function.set_time(time - time_step_size);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            source_function,
                                            tmp);
        forcing_terms.add(time_step_size * (1 - theta), tmp);
        
        system_rhs += forcing_terms;
        
        // Add natural boundary conditions
        for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
        {
            if ((params.boundary_conditions.implementation_types[boundary] != "natural"))
            {
                continue;
            }
            
            std::set<types::boundary_id> dealii_boundary_id = {boundary}; // @todo: This throws a warning
            
            *boundary_functions[boundary].set_time(time);
            
            MyVectorTools::my_create_boundary_right_hand_side(
                dof_handler,
                QGauss<dim-1>(fe.degree+1),
                *boundary_functions[boundary],
                tmp,
                dealii_boundary_id);
            forcing_terms = tmp;
            forcing_terms *= time_step_size * theta;
                
            *boundary_functions[boundary].set_time(time - time_step_size);
            MyVectorTools::my_create_boundary_right_hand_side(
                dof_handler,
                QGauss<dim-1>(fe.degree+1),
                *boundary_functions[boundary],
                tmp,
                dealii_boundary_id);
             
            forcing_terms.add(time_step_size * (1 - theta), tmp);
            
            system_rhs += forcing_terms;
        }
        
        //
        system_matrix.copy_from(mass_matrix);
        
        system_matrix.add(params.time.semi_implicit_theta * params.time.step_size, convection_diffusion_matrix);

        constraints.condense (system_matrix, system_rhs);

        {
            // Apply strong boundary conditions
            std::map<types::global_dof_index, double> boundary_values;
            for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
            {
                if (params.boundary_conditions.implementation_types[boundary] != "strong") 
                {
                    continue;
                }
                
                *boundary_functions[boundary].set_time(time);
                
                VectorTools::interpolate_boundary_values
                    (
                    dof_handler,
                    boundary,
                    *boundary_functions[boundary],
                    boundary_values
                    );
            }
            MatrixTools::apply_boundary_values(boundary_values,
                                               system_matrix,
                                               solution,
                                               system_rhs);
        }

        solve_time_step();

        write_solution();

        if (params.mms.enabled)
        {
            mms_tabulate_L2_error();
        }
        
        if ((time_step_counter == 1) &&
            (pre_refinement_step < params.refinement.adaptive.initial_cycles))
        {
            adaptive_refine();
            ++pre_refinement_step;

            tmp.reinit (solution.size());

            std::cout << std::endl;

            goto start_time_iteration;
        }
        else if ((time_step_counter > 0) 
                 && (params.refinement.adaptive.interval > 0)  // % 0 (mod 0) is undefined
                 && (time_step_counter % params.refinement.adaptive.interval == 0))
        {
            for (unsigned int cycle = 0;
                 cycle < params.refinement.adaptive.cycles_at_interval; cycle++)
            {
                adaptive_refine();
            }
            tmp.reinit (solution.size());
            
        }
        old_solution = solution;
    }
    
    // Save data for FEFieldFunction so that it can be loaded for initialization
    FEFieldTools::save_field_parts(triangulation, dof_handler, solution);
    
    // Write error table
    if (params.mms.enabled)
    {
        mms_write_error_table();
    }
    
    // Write the solution table containing pointwise values for every timestep.
    if (dim == 1)
    {
        this->write_1D_solution_table(this->solution_table_1D_file_name);
    }
    
    // Cleanup
    this->triangulation.set_manifold(0);
    
  }
}
