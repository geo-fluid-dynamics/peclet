#ifndef pde_parameters_h
#define pde_parameters_h

#include <vector>
#include <iostream>
#include <fstream>

#include <deal.II/base/parameter_handler.h>

/*
    
    @brief: Encapsulates parameter handling and paramter input file handling.

    @detail:
    
        Originally the ParameterReader from deal.II's step-26 was used;
        but that was too prohibitive. The approach in this file isolates the details
        of the input file handling from the rest of the program,
        allowing for rich data structures and simplifying the user code.

        The goal is to allow one to run the program with different parameters,
        without having to recomplile the program.

        It is valuable to structure the parameter data as done here, to greatly
        simplify the writing and debugging of the code.

        This also makes it simple to insantiate a PDE Model in a user program 
        and then change it's parameters directly without having to use any intermediate text files.
    
    @author: A. Zimmerman <zimmerman@aices.rwth-aachen.de>
    
*/

namespace PDE
{
    // @todo: This should be a Parameters class instead.
    //        In the beginning I didn't realize how involved it would become.
    namespace Parameters
    {   
        using namespace dealii;
    
        struct ConvectionDiffusionEquation
        {
            double reference_peclet_number;
            std::string convection_velocity_function_name;
            std::list<double> convection_velocity_function_double_arguments; 
        };
    
        struct BoundaryConditions
        {
            std::vector<std::string> implementation_types;
            std::vector<std::string> function_names;
            std::list<double> function_double_arguments;
        };
        
        struct InitialValues
        {
            std::string function_name;
            std::list<double> function_double_arguments; 
        };
        
        struct Geometry
        {
            unsigned int dim;
            std::string grid_name;
            std::vector<double> sizes;
            std::vector<double> transformations;
        };
        
        struct AdaptiveRefinement
        {
            unsigned int initial_cycles;
            unsigned int max_level;
            unsigned int max_cells;
            unsigned int interval;
            unsigned int cycles_at_interval;
            double refine_fraction;
            double coarsen_fraction;
        };
            
        struct Refinement
        {
            unsigned int initial_global_cycles;
            unsigned int initial_boundary_cycles;
            std::vector<unsigned int> boundaries_to_refine;
            AdaptiveRefinement adaptive;
        };
        
        struct Time
        {
            double end_time;
            double step_size;
            double semi_implicit_theta;
        };
        
        struct IterativeSolver
        {
            std::string method;
            unsigned int max_iterations;
            double tolerance;
            bool normalize_tolerance;
        };
        
        struct Output
        {
            bool write_solution_vtk;
            bool write_solution_table;
            int time_step_interval;
        };
        
        struct MMS
        {
            bool enabled;
            double initial_values_perturbation;
            std::string name;
            std::list<double> double_arguments;
        };
        
        struct StructuredParameters
        {
            ConvectionDiffusionEquation pde;
            BoundaryConditions boundary_conditions;
            InitialValues initial_values;
            Geometry geometry;
            Refinement refinement;
            Time time;
            IterativeSolver solver;
            Output output;
            MMS mms;
        };
      
        void declare(ParameterHandler &prm)
        {
            prm.enter_subsection("pde");
            {
                prm.declare_entry("reference_peclet_number", "1.",
                    Patterns::Double(0.));
                    
                prm.declare_entry("convection_velocity_function_name", "MMS",
                    Patterns::List(Patterns::Selection("MMS | constant | ramp")));
                    
                prm.declare_entry("convection_velocity_function_double_arguments", "",
                    Patterns::List(Patterns::Double()));
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection ("geometry");
            {
                prm.declare_entry("dim", "2",
                    Patterns::Integer(1, 3));
                    
                prm.declare_entry("grid_name", "hyper_shell",
                     Patterns::Selection("hyper_cube | hyper_shell | hemisphere_cylinder_shell"
                                      " | cylinder | cylinder_with_split_boundaries"
                                      " | hyper_cube_with_cylindrical_hole"),
                     "Select the name of the geometry and grid to generate."
                     "\nhyper_shell"
                     "\n\tInner boundary ID = 0"
                     "\n\tOuter boundary ID = 1"
                     
                     "\nhemisphere_cylinder_shell"
                     
                     "\ncylinder:"
                     "\n\tBoundary ID's"
                     "\n\t\t0: Heat flux"
                     "\n\t\t1: Outflow"
                     "\n\t\t2: Domain sides"
                     "\n\t\t3: Inflow"
                     
                     "\nhyper_cube_with_cylindrical_hole:"
                     "\n\tOuter boundary ID = 0"
                     "\n\tInner spherical boundary ID = 1");
                     
                prm.declare_entry("sizes", "0.5, 1.",
                    Patterns::List(Patterns::Double(0.)),
                    "Set the sizes for the grid's geometry."
                    "\n hyper_shell:"
                        "{inner_radius, outer_radius}"
                    "\n  hemisphere_cylinder_shell: "
                         "{inner_sphere_radius, outer_sphere_radius, "
                         "inner_cylinder_length, outer_cylinder_length}"
                    "\n cylinder: "
                        "{L0, L1, L2}"
                    "\n  hyper_cube_with_cylindrical_hole : {hole_radius, half_of_outer_edge_length}");
                    
                prm.declare_entry("transformations", "0., 0., 0.",
                    Patterns::List(Patterns::Double()),
                    "Set the rigid body transformation vector."
                    "\n  2D : {shift_along_x, shift_along_y, rotate_about_z}"
                    "\n  3D : {shift_along_x, shift_along_y, shift_along_z, "
                              "rotate_about_x, rotate_about_y, rotate_about_z}");
                              
            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("initial_values");
            {
                prm.declare_entry("function_name", "MMS",
                    Patterns::List(Patterns::Selection("constant | MMS | ramp | interpolate_old_field")));
                    
                prm.declare_entry("function_double_arguments", "",
                    Patterns::List(Patterns::Double())); 
                    
            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("boundary_conditions");
            {
                // Each of these lists needs a value for every boundary, in order
                prm.declare_entry("implementation_types", "natural, strong",
                    Patterns::List(Patterns::Selection("natural | strong")),
                    "Type of boundary conditions to apply to each boundary");  
                    
                prm.declare_entry("function_names", "MMS, MMS",
                    Patterns::List(Patterns::Selection("constant | MMS | ramp ")),
                    "Names of functions to apply to each boundary");
                    
                prm.declare_entry("function_double_arguments", "",
                    Patterns::List(Patterns::Double()),
                    "This list of doubles will be popped from front to back as needed."
                    "\nThis puts some work on the user to greatly ease development."
                    "\nHere are some tips:"
                    "\n\t- The function values will only be popped during initialization."
                    "\n\t- Boundaries will be handled in order of their ID's."
                    "\n\t- If a function needs a Point as an argument, then it will pop doubles to make the point in order."); 
    
            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("refinement");
            {
                prm.declare_entry("initial_global_cycles", "2",
                    Patterns::Integer(),
                    "Initially globally refine the grid this many times "
                    "without using any error measure");
                    
                prm.declare_entry("initial_boundary_cycles", "4",
                    Patterns::Integer(),
                    "Initially refine the grid this many times"
                    "near the boundaries that are listed for refinement");
                    
                prm.declare_entry("boundaries_to_refine", "0",
                    Patterns::List(Patterns::Integer()),
                    "Refine cells that contain these boundaries");
                    
                prm.enter_subsection ("adaptive");
                {
                    prm.declare_entry("initial_cycles", "0",
                        Patterns::Integer(),
                        "Refine grid adaptively using an error measure "
                        "this many times before beginning the time stepping.");
                        
                    prm.declare_entry("interval", "0",
                        Patterns::Integer(),
                        "Only refine the grid after every occurence of "
                        "this many time steps.");
                        
                    prm.declare_entry("max_level", "10",
                        Patterns::Integer(),
                        "Max grid refinement level");
                        
                    prm.declare_entry("max_cells", "2000",
                        Patterns::Integer(),
                        "Skip grid refinement if the number of active cells "
                        "already exceeds this");
                        
                    prm.declare_entry("refine_fraction", "0.3",
                        Patterns::Double(),
                        "Fraction of cells to refine");
                        
                    prm.declare_entry("coarsen_fraction", "0.3",
                        Patterns::Double(),
                        "Fraction of cells to coarsen");
                        
                    prm.declare_entry("cycles_at_interval", "5",
                        Patterns::Integer(),
                        "Max grid refinement level");
                        
                }
                prm.leave_subsection();
                
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection ("time");
            {
                prm.declare_entry("end_time", "1.",
                    Patterns::Double(0.),
                    "End the time-dependent simulation once this time is reached.");
                    
                prm.declare_entry("step_size", "0.01",
                    Patterns::Double(1.e-16),
                    "End the time-dependent simulation once this time is reached.");
                    
                prm.declare_entry("semi_implicit_theta", "0.51",
                    Patterns::Double(0., 1.),
                    "This is the theta parameter for the theta-family of "
                    "semi-implicit time integration schemes."
                    " Choose any value between zero and one."
                    " 0 = fully explicit; 0.5 = 'Crank-Nicholson'"
                    " ; 1 = fully implicit");
                    
            }
            prm.leave_subsection();
            
            prm.enter_subsection("solver");
            {
                prm.declare_entry("method", "CG",
                     Patterns::Selection("CG | BiCGStab"));
                     
                prm.declare_entry("max_iterations", "1000",
                    Patterns::Integer(0));
                    
                prm.declare_entry("tolerance", "1e-8",
                    Patterns::Double(0.));
                    
                prm.declare_entry("normalize_tolerance", "true",
                    Patterns::Bool(),
                    "If true, then the residual will be multiplied by the L2-norm of the RHS"
                    " before comparing to the tolerance.");
            }
            prm.leave_subsection();
            
            prm.enter_subsection("output");
            {
                prm.declare_entry("write_solution_vtk", "true", Patterns::Bool());
                prm.declare_entry("write_solution_table", "false", Patterns::Bool(),
                    "This allow for simple export of 1D solutions into a table format"
                    " easily read by MATLAB."
                    "\nThe way this is currently implemented takes a great deal of memory"
                    ", so you should probably only use this in 1D.");
                prm.declare_entry("time_step_interval", "1", Patterns::Integer(0),
                    "Solutions will only be written at every time_step_interval time step."
                    "\nSet to one to output at every time step."
                    "\n Set to zero to output only the final time.");
            }
            prm.leave_subsection();
            
            prm.enter_subsection("mms");
            {
                prm.declare_entry("enabled", "true", Patterns::Bool());
                prm.declare_entry("initial_values_perturbation", "1.001",
                    Patterns::Double());
                prm.declare_entry("name", "ConstantConvection1D",
                    Patterns::Selection("ConstantConvection1D"));
                prm.declare_entry("double_arguments", "-1., -1., 10.",
                    Patterns::List(Patterns::Double()),
                    "ConstantConvection1D:"
                    "\tdouble_arguments: convection_velocity, dirichlet_value, rate_to_steady");
            }
            prm.leave_subsection();
            
        }
    
    
        template<typename ItemType>
        std::vector<ItemType> get_vector(ParameterHandler &prm, std::string parameter_name)
        {
            std::vector<std::string> strings = Utilities::split_string_list(prm.get(parameter_name));
            std::vector<ItemType> items;
            for (auto &string : strings) 
            {
                std::stringstream parser(string);
                ItemType item;
                parser >> item;
                items.push_back(item);
            }
            return items;
        }    
        
        
        StructuredParameters get_parameters(const std::string parameter_file="")
        {
            ParameterHandler prm;
            declare(prm);
            
            if (parameter_file != "")
            {
                prm.read_input(parameter_file);    
            }
            
            // Print a log file of all the ParameterHandler parameters
            std::ofstream parameter_log_file("used_parameters.prm");
            assert(parameter_log_file.good());
            prm.print_parameters(parameter_log_file, ParameterHandler::Text);
            
            // Structure the parameters so that we can stop working with ParameterHandler
            StructuredParameters p;
            
            
            prm.enter_subsection("geometry");
            {
                p.geometry.dim = prm.get_integer("dim");
                p.geometry.grid_name = prm.get("grid_name");
                p.geometry.sizes = get_vector<double>(prm, "sizes");
                p.geometry.transformations = get_vector<double>(prm, "transformations");    
            }
            prm.leave_subsection();
            

            prm.enter_subsection("pde");
            {
                p.pde.reference_peclet_number = 
                    prm.get_double("reference_peclet_number");    
                p.pde.convection_velocity_function_name = 
                    prm.get("convection_velocity_function_name");
                
                std::vector<double> vector = get_vector<double>(
                    prm, "convection_velocity_function_double_arguments");
                
                for (auto v : vector)
                {
                    p.pde.convection_velocity_function_double_arguments.push_back(v);
                }

            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("boundary_conditions");
            {
                p.boundary_conditions.implementation_types = get_vector<std::string>(prm, "implementation_types");
                p.boundary_conditions.function_names = get_vector<std::string>(prm, "function_names");
                
                std::vector<double> vector = get_vector<double>(prm, "function_double_arguments");
                
                for (auto v : vector)
                {
                    p.boundary_conditions.function_double_arguments.push_back(v);
                }
                
            }    
            prm.leave_subsection();
            
            
            prm.enter_subsection("initial_values");
            {               
                p.initial_values.function_name = prm.get("function_name"); 
                
                std::vector<double> vector = get_vector<double>(prm, "function_double_arguments");
                
                for (auto v : vector)
                {
                    p.initial_values.function_double_arguments.push_back(v)    ;
                }
                
            }    
            prm.leave_subsection();
            
            
            prm.enter_subsection("refinement");
            {
                
                p.refinement.initial_global_cycles = prm.get_integer("initial_global_cycles");
                p.refinement.initial_boundary_cycles = prm.get_integer("initial_boundary_cycles");
                p.refinement.boundaries_to_refine = get_vector<unsigned int>(prm, "boundaries_to_refine");
                
                prm.enter_subsection("adaptive");
                {
                    p.refinement.adaptive.initial_cycles = prm.get_integer("initial_cycles");
                    p.refinement.adaptive.max_level = prm.get_integer("max_level");
                    p.refinement.adaptive.max_cells = prm.get_integer("max_cells");
                    p.refinement.adaptive.interval = prm.get_integer("interval");
                    p.refinement.adaptive.cycles_at_interval = prm.get_integer("cycles_at_interval");
                    p.refinement.adaptive.refine_fraction = prm.get_double("refine_fraction");
                    p.refinement.adaptive.coarsen_fraction = prm.get_double("coarsen_fraction");    
                }        
                
                prm.leave_subsection();
                
            }
            prm.leave_subsection();
                
                
            prm.enter_subsection("time");
            {
                p.time.end_time = prm.get_double("end_time");
                p.time.step_size = prm.get_double("step_size");
                p.time.semi_implicit_theta = prm.get_double("semi_implicit_theta");
            }    
            prm.leave_subsection();
            
            
            prm.enter_subsection("solver");
            {
                p.solver.method = prm.get("method");
                p.solver.max_iterations = prm.get_integer("max_iterations");
                p.solver.tolerance = prm.get_double("tolerance");
                p.solver.normalize_tolerance = prm.get_bool("normalize_tolerance");
            }    
            prm.leave_subsection(); 
            
            prm.enter_subsection("output");
            {
                p.output.write_solution_vtk = prm.get_bool("write_solution_vtk");
                p.output.write_solution_table = prm.get_bool("write_solution_table");
                p.output.time_step_interval = prm.get_integer("time_step_interval");
            }
            prm.leave_subsection();
            
            prm.enter_subsection("mms");
            {
                p.mms.enabled = prm.get_bool("enabled");
                p.mms.initial_values_perturbation = 
                    prm.get_double("initial_values_perturbation");
                    
                p.mms.name = prm.get("name");
                
                std::vector<double> vector = get_vector<double>(prm, "double_arguments");
                
                for (auto v : vector)
                {
                    p.mms.double_arguments.push_back(v)    ;
                }

            }
            prm.leave_subsection();
            
            return p;
        }
     
    }
    
}

#endif
