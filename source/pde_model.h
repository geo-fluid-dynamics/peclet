 /*
 * Based on deal.II Tutorial 26 by Wolfgang Bangerth, Texas A&M University, 2013
 * Extended for dimice-heat-dealii by Alexander Zimmerman, RWTH AAchen University, 2016
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
    void append_solution_to_table();
    void write_solution_table(std::string file_name);
    
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
    unsigned int         timestep_number;
    
    Point<dim> spherical_manifold_center;
    
    std::vector<unsigned int> manifold_ids;
    std::vector<std::string> manifold_descriptors;
    
    std::vector<double> convection_velocity;
    double diffusivity;
    
    TableHandler solution_table;
    std::string solution_table_file_name = "solution_table.txt";
    
  };

  template<int dim>
  Model<dim>::Model()
    :
    params(Parameters::get_parameters()),
    fe(1),
    dof_handler(triangulation),
    convection_velocity(dim)
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
                                      
    ConstantFunction<dim> diffusivity(this->diffusivity);
    
    ConstantFunction<dim> convection_velocity(this->convection_velocity);
                                      
    MyMatrixCreator::create_convection_diffusion_matrix(
        dof_handler,
        QGauss<dim>(fe.degree+1),
        convection_diffusion_matrix,
        &diffusivity, 
        &convection_velocity
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
  
  #include "pde_solution_table.h"
  
  template<int dim>
  void Model<dim>::write_solution()
  {
      
    if (params.output.write_solution_vtk)
    {
        Output::write_solution_to_vtk(
            "solution-"+Utilities::int_to_string(timestep_number)+".vtk",
            this->dof_handler,
            this->solution);    
    }
    
    if (params.output.write_solution_table)
    {
        this->append_solution_to_table();
    }
    
  }
  
  template<int dim>
  void Model<dim>::run()
  {
    if (dim == 1)
    {
        std::remove(solution_table_file_name.c_str()); // In 1D, the solution will be appended here at every time step.    
    }        
    
    for (unsigned int axis = 0; axis < dim; axis++)
    {
        this->unit_convection_velocity[axis] = params.pde.unit_convection_velocity[axis];
    }
    
    this->peclet_number = params.pde.peclet_number;
    
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
    
    // Make initial values function
    ConstantFunction<dim> constant_function(0.);
    
    Function<dim>* initial_values_function = &constant_function;
    
    MMS::InitialValues mms_initial_values_function(params.mms.iv_perturbation);
    mms_initial_values_function.set_time(0.);
    
    Point<dim> ramp_start_point, ramp_end_point;
    
    double ramp_start_position = 0.,
           ramp_end_position = 0.,
           ramp_start_value = 0.,
           ramp_end_value = 0.;
            
    if (params.initial_values.function_name == "ramp")
    {
        for (unsigned int axis = 0; axis < dim; axis++)
        {
            ramp_start_point[axis] = params.initial_values.function_double_arguments.front();
            params.initial_values.function_double_arguments.pop_front();
        }
        
        for (unsigned int axis = 0; axis < dim; axis++)
        {
            ramp_end_point[axis] = params.initial_values.function_double_arguments.front();
            params.initial_values.function_double_arguments.pop_front();
        }
        
        ramp_start_position = params.initial_values.function_double_arguments.front();
        params.initial_values.function_double_arguments.pop_front();
        
        ramp_end_position = params.initial_values.function_double_arguments.front();
        params.initial_values.function_double_arguments.pop_front();
        
        ramp_start_value = params.initial_values.function_double_arguments.front();
        params.initial_values.function_double_arguments.pop_front();
        
        ramp_end_value = params.initial_values.function_double_arguments.front();
        params.initial_values.function_double_arguments.pop_front();
        
    }
    
    MyFunctions::RampFunctionAlongLine<dim> ramp_function(
            ramp_start_point,
            ramp_end_point,
            ramp_start_position,
            ramp_end_position,
            ramp_start_value,
            ramp_end_value);
            
    
    Triangulation<dim> field_grid;
    DoFHandler<dim> field_dof_handler(field_grid);
    Vector<double> field_solution;
    
    if (params.initial_values.function_name != "interpolate_old_field")
    { // This will write files that need to exist.
        setup_system(true);
        FEFieldTools::save_field_parts(this->triangulation, this->dof_handler, this->solution); 
    }
    
    FEFieldTools::load_field_parts(
        field_grid,
        field_dof_handler,
        field_solution,
        this->fe);
    
    MyFunctions::ExtrapolatedField<dim> field_function(
        field_dof_handler,
        field_solution);
    

    if (params.initial_values.function_name == "interpolate_old_field")
    {
        initial_values_function = & field_function;                      
    }
    else if (params.initial_values.function_name == "constant")
    { 
        constant_function = ConstantFunction<dim,double>(
            params.initial_values.function_double_arguments.front());
        initial_values_function = & constant_function;
                        
    }
    else if (params.initial_values.function_name == "ramp")
    {
        initial_values_function =  & ramp_function;
        
    }
    else if (params.initial_values.function_name == "MMS")
    { 
        assert(params.mms.enabled);
        initial_values_function = & mms_initial_values_function;
                        
    }
    
    // Make RHS functions
    ConstantFunction<dim> zero_source_function(0.);
    
    Function<dim>* source_function = &zero_source_function;
    
    MMS::Source<dim> mms_source_function(params.mms.max_Pe);
    
    if (params.mms.enabled)
    {
        source_function = &mms_source_function;
    }
    
    
    // Make boundary functions
    
    unsigned int boundary_count = params.boundary_conditions.implementation_types.size();
    
    assert(params.boundary_conditions.function_names.size() == boundary_count);

    std::vector<ConstantFunction<dim>> constant_functions;
    
    MMS::ManufacturedSolution dirichlet_boundary_function;
    MMS::NeumannBoundary neumann_boundary_function(this->tria, params.mms.max_Pe);
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
    {
        std::string boundary_type = params.boundary_conditions.implementation_types[boundary];
        std::string function_name = params.boundary_conditions.function_names[boundary];
        
        if ((function_name == "constant"))
        {
            double value = params.boundary_conditions.function_double_arguments.front();
            params.boundary_conditions.function_double_arguments.pop_front();
            if (boundary_type == "natural")
            {
                /*
                @todo
                
                    The natural boundary condition should also be divided by the thermal diffusivity.
                    Note: Changing this will break existing Neumann tests.
                
                */
                
                value *= params.time.time_step;
            }

            constant_functions.push_back(ConstantFunction<dim>(value));
            
        }
    }
        
    // Organize boundary functions to simplify application during the time loop
    
    
    std::vector<Function<dim>*> boundary_functions;
    unsigned int constant_function_index = 0;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)        
    {
        std::string boundary_type = params.boundary_conditions.implementation_types[boundary];
        std::string function_name = params.boundary_conditions.function_names[boundary];

        if ((function_name == "constant"))
        {
            assert(constant_function_index < constant_functions.size());
            boundary_functions.push_back(&constant_functions[constant_function_index]);
            constant_function_index++;
        }
        else if ((function_name == "MMS"))
        {
            assert(params.mms.enabled);
            if (boundary_type == "strong")
            {
                boundary_functions.push_back(&mms_dirichlet);
            }
            else if (boundary_type == "natural")
            {
                boundary_functions.push_back(&mms_neumann);
            }

        }
        
    }
    
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
    timestep_number = 0; // @todo: Expose initial time step as a parameter
    time            = 0; // @todo: Expose initial time as a parameter
    write_solution();
    
    double epsilon = 1e-16;
    while (time < params.time.end_time - epsilon)
    {
        ++timestep_number;
        time = params.time.time_step*timestep_number; // Incrementing the time directly would accumulate errors
        
        std::cout << "Time step " << timestep_number << " at t=" << time << std::endl;

        mass_matrix.vmult(system_rhs, old_solution);

        convection_diffusion_matrix.vmult(tmp, old_solution);
        
        system_rhs.add(-(1. - params.time.semi_implicit_theta) * params.time.time_step, tmp);
        
        // Add source terms
        source_function.set_time(time);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            source_function,
                                            tmp);
        forcing_terms = tmp;
        forcing_terms *= time_step * theta;
        
        source_function.set_time(time - time_step);
        VectorTools::create_right_hand_side(dof_handler,
                                            QGauss<dim>(fe.degree+1),
                                            source_function,
                                            tmp);
        forcing_terms.add(time_step * (1 - theta), tmp);
        
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
            forcing_terms *= time_step * theta;
                
            *boundary_functions[boundary].set_time(time - time_step);
            MyVectorTools::my_create_boundary_right_hand_side(
                dof_handler,
                QGauss<dim-1>(fe.degree+1),
                *boundary_functions[boundary],
                tmp,
                dealii_boundary_id);
             
            forcing_terms.add(time_step * (1 - theta), tmp);
            
            system_rhs += forcing_terms;
        }
        
        //
        system_matrix.copy_from(mass_matrix);
        
        system_matrix.add(params.time.semi_implicit_theta * params.time.time_step, convection_diffusion_matrix);

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

        if ((timestep_number == 1) &&
            (pre_refinement_step < params.refinement.adaptive.initial_cycles))
        {
            adaptive_refine();
            ++pre_refinement_step;

            tmp.reinit (solution.size());

            std::cout << std::endl;

            goto start_time_iteration;
        }
        else if ((timestep_number > 0) 
                 && (params.refinement.adaptive.interval > 0)  // % 0 (mod 0) is undefined
                 && (timestep_number % params.refinement.adaptive.interval == 0))
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
    
    // Write the solution table containing pointwise values for every timestep.
    if (params.output.write_solution_table)
    {
        this->write_solution_table(this->solution_table_file_name);
    }
    
    // Cleanup
    this->triangulation.set_manifold(0);
    
  }
}
