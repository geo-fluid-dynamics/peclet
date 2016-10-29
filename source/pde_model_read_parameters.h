    
    template <int dim>
    void Model<dim>::read_parameters(const std::string parameter_file="")
    {
        ParameterHandler prm;
        Parameters::declare(prm);
        
        // Declare MMS stuff
        /*
        This is separated from declare() because it is a class method.
        The use of ParsedFunction with MMS in this project is new, and it is the first occurrence
        of needing direct access to the Model class when declaring parameters. Maybe there is 
        a better way to implement MMS that avoids all of this; but my other attempts
        have been much more convoluted.
        */
        prm.enter_subsection("mms");
        {
            
            prm.enter_subsection("solution");
            {
                this->mms_solution.declare_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("source");
            {
                this->mms_source.declare_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("neumann");
            {
                this->mms_neumann.declare_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("velocity");
            {
                if (dim == 1)
                {
                    this->mms_velocity = &this->mms_velocity_1D;
                }
                else if (dim == 2)
                {
                    this->mms_velocity = &this->mms_velocity_2D;
                }
                else if (dim == 3)
                {
                    this->mms_velocity = &this->mms_velocity_3D;
                }
                this->mms_velocity->declare_parameters(prm, dim);
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
        
        //
        if (parameter_file != "")
        {
            prm.read_input(parameter_file);    
        }
        
        // Print a log file of all the ParameterHandler parameters
        std::ofstream parameter_log_file("used_parameters.prm");
        assert(parameter_log_file.good());
        prm.print_parameters(parameter_log_file, ParameterHandler::Text);
        
        prm.enter_subsection("geometry");
        {
            this->params.geometry.grid_name = prm.get("grid_name");
            this->params.geometry.sizes = Parameters::get_vector<double>(prm, "sizes");
            this->params.geometry.transformations = 
                Parameters::get_vector<double>(prm, "transformations");    
        }
        prm.leave_subsection();
        

        prm.enter_subsection("pde");
        {
            this->params.pde.reference_peclet_number = 
                prm.get_double("reference_peclet_number");    
                
            this->params.pde.convection_velocity_function_name = 
                prm.get("convection_velocity_function_name");
            
            this->params.pde.convection_velocity_function_double_arguments = 
                Parameters::get_vector<double>(prm, "convection_velocity_function_double_arguments");
                
            this->params.pde.source_function_name = 
                prm.get("source_function_name");
            
            this->params.pde.source_function_double_arguments = 
                Parameters::get_vector<double>(prm, "source_function_double_arguments");
                
        }
        prm.leave_subsection();
        
        
        prm.enter_subsection("boundary_conditions");
        {
            this->params.boundary_conditions.implementation_types = 
                Parameters::get_vector<std::string>(prm, "implementation_types");
            this->params.boundary_conditions.function_names = 
                Parameters::get_vector<std::string>(prm, "function_names");
            
            std::vector<double> vector = 
                Parameters::get_vector<double>(prm, "function_double_arguments");
            
            for (auto v : vector)
            {
                this->params.boundary_conditions.function_double_arguments.push_back(v);
            }
            
        }    
        prm.leave_subsection();
        
        
        prm.enter_subsection("initial_values");
        {               
            this->params.initial_values.function_name = prm.get("function_name"); 
            
            std::vector<double> vector = 
                Parameters::get_vector<double>(prm, "function_double_arguments");
            
            for (auto v : vector)
            {
                this->params.initial_values.function_double_arguments.push_back(v);
            }
            
        }    
        prm.leave_subsection();
        
        
        prm.enter_subsection("refinement");
        {
            
            this->params.refinement.initial_global_cycles = prm.get_integer("initial_global_cycles");
            this->params.refinement.initial_boundary_cycles = prm.get_integer("initial_boundary_cycles");
            this->params.refinement.boundaries_to_refine = 
                Parameters::get_vector<unsigned int>(prm, "boundaries_to_refine");
            
            prm.enter_subsection("adaptive");
            {
                this->params.refinement.adaptive.initial_cycles = prm.get_integer("initial_cycles");
                this->params.refinement.adaptive.max_level = prm.get_integer("max_level");
                this->params.refinement.adaptive.max_cells = prm.get_integer("max_cells");
                this->params.refinement.adaptive.interval = prm.get_integer("interval");
                this->params.refinement.adaptive.cycles_at_interval = prm.get_integer("cycles_at_interval");
                this->params.refinement.adaptive.refine_fraction = prm.get_double("refine_fraction");
                this->params.refinement.adaptive.coarsen_fraction = prm.get_double("coarsen_fraction");    
            }        
            
            prm.leave_subsection();
            
        }
        prm.leave_subsection();
            
            
        prm.enter_subsection("time");
        {
            this->params.time.end_time = prm.get_double("end_time");
            this->params.time.step_size = prm.get_double("step_size");
            this->params.time.global_refinement_levels = 
                prm.get_integer("global_refinement_levels");
            this->params.time.semi_implicit_theta = prm.get_double("semi_implicit_theta");
            this->params.time.stop_when_steady = prm.get_bool("stop_when_steady");
        }    
        prm.leave_subsection();
        
        
        prm.enter_subsection("solver");
        {
            this->params.solver.method = prm.get("method");
            this->params.solver.max_iterations = prm.get_integer("max_iterations");
            this->params.solver.tolerance = prm.get_double("tolerance");
            this->params.solver.normalize_tolerance = prm.get_bool("normalize_tolerance");
        }    
        prm.leave_subsection(); 
        
        prm.enter_subsection("output");
        {
            this->params.output.write_solution_vtk = prm.get_bool("write_solution_vtk");
            this->params.output.write_solution_table = prm.get_bool("write_solution_table");
            this->params.output.time_step_interval = prm.get_integer("time_step_interval");
        }
        prm.leave_subsection();
        
        prm.enter_subsection("mms");
        {
            this->params.mms.enabled = prm.get_bool("enabled");
            
            prm.enter_subsection("solution");
            {
                this->mms_solution.parse_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("source");
            {
                this->mms_source.parse_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("neumann");
            {
                this->mms_neumann.parse_parameters(prm);
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("velocity");
            {
                this->mms_velocity->parse_parameters(prm);
            }
            prm.leave_subsection();
            
            this->params.mms.initial_values_perturbation =
                prm.get_double("initial_values_perturbation");
        }
        prm.leave_subsection();

    }