 
    /*
    This file contains the initialization of many different types of functions that are needed
    throughout the program. Working with deal.II's Function class has been interesting, and I'm 
    sure many of my choices are unorthodox. The most important lesson learned has been that 
    a Function<dim>* can point to any class derived from Function<dim>. The general design pattern
    then is to instantitate all of the functions that might be needed, and then to point to the ones
    actually being used. Since this looks quite messy and could distract from the rest of the program, this design pattern is mostly contained in this file.
    
    Note that I recently discovered ParsedFunction, which obsoletes some of what I had been trying to do in this file, e.g. with manually implemented constant functions and ramp functions. Ultimately the ParsedFunction is not enough. For example this file contains the option for an initial values function that interpolates an old solution loaded from disk. So in most cases one should be able to use a ParsedFunction, but the generality of Function<dim>* (function pointers) allows for a standard way to account for any possible derived class of Function<dim>.
    
    Also this file contains most of what was needed to implement general boundary conditions. I think that the approach here is quite powerful and flexible.
    */
 
    // Velocity function
    
    std::vector<double> constant_velocity(dim);
    if (this->params.pde.velocity_function_name == "constant")
    {
        for (unsigned int i = 0; i < dim; i++)
        {
            constant_velocity[i] =
                this->params.pde.velocity_function_double_arguments[i];
        }    
    }
    
    ConstantFunction<dim> constant_velocity_function(constant_velocity);
    
    if (this->params.pde.velocity_function_name == "parsed")
    {
        this->velocity_function = &parsed_velocity_function;
    }
    else if (this->params.pde.velocity_function_name == "constant")
    {
        this->velocity_function = &constant_velocity_function;
    }
    else
    {
        Assert(false, ExcNotImplemented());
    }
    
    // Diffusivity function
    
    double constant_diffusivity;
    if (this->params.pde.diffusivity_function_name == "constant")
    {
        constant_diffusivity = this->params.pde.diffusivity_function_double_arguments[0];
    }
    
    ConstantFunction<dim> constant_diffusivity_function(constant_diffusivity);
    
    if (this->params.pde.diffusivity_function_name == "parsed")
    {
        this->diffusivity_function = &parsed_diffusivity_function;
    }
    else if (this->params.pde.diffusivity_function_name == "constant")
    {
        this->diffusivity_function = &constant_diffusivity_function;
    }
    else
    {
        Assert(false, ExcNotImplemented());
    }

    // Make initial values function
    ConstantFunction<dim> constant_function(0.);
    
    initial_values_function = &constant_function;

    Point<dim> ramp_start_point, ramp_end_point;
    
    double ramp_start_position = 0.,
           ramp_end_position = 0.,
           ramp_start_value = 0.,
           ramp_end_value = 0.;
            
    if (this->params.initial_values.function_name == "ramp")
    {
        for (unsigned int axis = 0; axis < dim; axis++)
        {
            ramp_start_point[axis] = this->params.initial_values.function_double_arguments.front();
            this->params.initial_values.function_double_arguments.pop_front();
        }
        
        for (unsigned int axis = 0; axis < dim; axis++)
        {
            ramp_end_point[axis] = this->params.initial_values.function_double_arguments.front();
            this->params.initial_values.function_double_arguments.pop_front();
        }
        
        ramp_start_position = this->params.initial_values.function_double_arguments.front();
        this->params.initial_values.function_double_arguments.pop_front();
        
        ramp_end_position = this->params.initial_values.function_double_arguments.front();
        this->params.initial_values.function_double_arguments.pop_front();
        
        ramp_start_value = this->params.initial_values.function_double_arguments.front();
        this->params.initial_values.function_double_arguments.pop_front();
        
        ramp_end_value = this->params.initial_values.function_double_arguments.front();
        this->params.initial_values.function_double_arguments.pop_front();
        
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
    
    if (this->params.initial_values.function_name != "interpolate_old_field")
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
    

    if (this->params.initial_values.function_name == "interpolate_old_field")
    {
        this->initial_values_function = &field_function;                      
    }
    else if (this->params.initial_values.function_name == "constant")
    { 
        constant_function = ConstantFunction<dim,double>(
            this->params.initial_values.function_double_arguments.front());
        this->initial_values_function = &constant_function;
                        
    }
    else if (this->params.initial_values.function_name == "ramp")
    {
        this->initial_values_function =  &ramp_function;
        
    }
    else if (this->params.initial_values.function_name == "parsed")
    { 
        this->initial_values_function = &parsed_initial_values_function;
    }
    
    // Make source functions
    double constant_source_value = 0.;
    
    if (this->params.pde.source_function_name == "constant")
    {
        constant_source_value = params.pde.source_function_double_arguments[0];
    }
    
    ConstantFunction<dim> constant_source_function(constant_source_value);
    
    if (this->params.pde.source_function_name == "parsed")
    {
        this->source_function = &parsed_source_function;
    }
    else if (this->params.pde.source_function_name == "constant")
    {
        this->source_function = &constant_source_function;
    }
    
    // Make boundary functions
    
    unsigned int boundary_count = this->params.boundary_conditions.implementation_types.size();
    
    assert(params.boundary_conditions.function_names.size() == boundary_count);

    std::vector<ConstantFunction<dim>> constant_functions;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];
        
        if (function_name == "constant")
        {
            double value = this->params.boundary_conditions.function_double_arguments.front();
            this->params.boundary_conditions.function_double_arguments.pop_front();
            constant_functions.push_back(ConstantFunction<dim>(value));
        }
    }
        
    // Organize boundary functions to simplify application during the time loop
    
    unsigned int constant_function_index = 0;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)        
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];

        if (function_name == "constant")
        {
            assert(constant_function_index < constant_functions.size());
            this->boundary_functions.push_back(&constant_functions[constant_function_index]);
            constant_function_index++;
        }
        else if (function_name == "parsed")
        {
            this->boundary_functions.push_back(&parsed_boundary_function);
        }
        
    }
    
    // Verification
    if (this->params.verification.enabled)
    {
        if (this->params.verification.exact_solution_function_name == "parsed")
        {
            this->exact_solution_function = &parsed_exact_solution_function;
        }
        else
        {
            Assert(false, ExcNotImplemented());
        }
        
    }