
    // Make initial values function
    ConstantFunction<dim> constant_function(0.);
    
    Function<dim>* initial_values_function = &constant_function;
    
    MMS::InitialValues mms_initial_values_function(params.mms.iv_perturbation);
    mms_initial_values_function.set_time(0.);
    
    MMS::ManufacturedSolution mms_dirichlet_function();
    
    HyperBallBoundary<dim> mms_neumann_boundary(Point<dim>(), params.geometry.sizes[0]);
    
    MMS::NeumannBoundary mms_neumann_function(
        this->triangulation,
        mms_neumann_boundary,
        this->reference_peclet_number);
    
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
    
    MMS::Source<dim> mms_source_function(params.mms.max_peclet_number);
    
    if (params.mms.enabled)
    {
        source_function = &mms_source_function;
    }
    
    
    // Make boundary functions
    
    unsigned int boundary_count = params.boundary_conditions.implementation_types.size();
    
    assert(params.boundary_conditions.function_names.size() == boundary_count);

    std::vector<ConstantFunction<dim>> constant_functions;
    
    MMS::ManufacturedSolution dirichlet_boundary_function;
    MMS::NeumannBoundary neumann_boundary_function(this->tria, params.mms.max_peclet_number);
    
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