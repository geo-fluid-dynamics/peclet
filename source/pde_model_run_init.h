        /*
        Originally this file just read input parameters and was named accordingly,
        but it has become necessary to also instantitate some Functions that must survive
        through run time, and so this file is generally named "init". This file needs to be merged
        with pde_model_run_initialize_function.h and overhauled.
        */

        ParameterHandler prm;
        Parameters::declare<dim>(prm);
        
        /*
        @todo
        
        In older versions there was a strict separation between the declaration and the
        reading of parameters. Now there is no strict separation, primarily because often whether
        or not some parameters should be declared depends on the input of other parameters.
        There is still plenty of work to unify this new approach; e.g., as of this writing,
        all of the adaptive grid refinement parameters are being declared even though adaptive
        grid refinement hasn't been used in a very long time (because boundary layer refinement
        has been strictly better).
        
        */
        
        //
        if (parameter_file != "")
        {
            prm.read_input(parameter_file);    
        }
        
        // Print a log file of all the ParameterHandler parameters
        std::ofstream parameter_log_file("used_parameters.prm");
        assert(parameter_log_file.good());
        prm.print_parameters(parameter_log_file, ParameterHandler::Text);
        
        
        Functions::ParsedFunction<dim> parsed_exact_solution_function;
        Functions::ParsedFunction<dim> parsed_initial_values_function;
        Functions::ParsedFunction<dim> parsed_source_function;
        Functions::ParsedFunction<dim> parsed_boundary_function;
        Functions::ParsedFunction<dim> parsed_velocity_function(dim);

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
            this->params.pde.reference_peclet_number = prm.get_double("reference_peclet_number");    
            
            this->params.pde.velocity_function_name = prm.get("velocity_function_name");
            
            prm.enter_subsection("parsed_velocity_function");
            {
                parsed_velocity_function.parse_parameters(prm);    
            }
            prm.leave_subsection();
            
            this->params.pde.velocity_function_double_arguments = 
                Parameters::get_vector<double>(prm, "velocity_function_double_arguments");

            this->params.pde.source_function_name = prm.get("source_function_name");

            prm.enter_subsection("parsed_source_function");
            {
                parsed_source_function.parse_parameters(prm);
            }
            prm.leave_subsection();
            
            this->params.pde.source_function_double_arguments = 
                Parameters::get_vector<double>(prm, "source_function_double_arguments");
            
        }
        prm.leave_subsection();
        
        prm.enter_subsection("verification");
        {
            this->params.verification.enabled = prm.get_bool("enabled");
            
            this->params.verification.exact_solution_function_name = 
                    prm.get("exact_solution_function_name");
                    
            prm.enter_subsection("parsed_exact_solution_function");
            {
                parsed_exact_solution_function.parse_parameters(prm);    
            }
            prm.leave_subsection();
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

            prm.enter_subsection("parsed_function");
            {
                parsed_boundary_function.parse_parameters(prm);
            }
            prm.leave_subsection();
            
        }    
        prm.leave_subsection();
        
        prm.enter_subsection("initial_values");
        {               
            this->params.initial_values.function_name = prm.get("function_name"); 
            
            prm.enter_subsection("parsed_function");
            {
                parsed_initial_values_function.parse_parameters(prm);
            }
            prm.leave_subsection();

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