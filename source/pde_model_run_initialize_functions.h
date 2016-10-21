    // Manufactured solution and auxiliary conditions
    double mms_convection_velocity = this->params.mms.double_arguments.front();
    this->params.mms.double_arguments.pop_front();
    
    double mms_dirichlet_value = this->params.mms.double_arguments.front();
    this->params.mms.double_arguments.pop_front();
    
    double mms_rate_to_steady = this->params.mms.double_arguments.front();
    this->params.mms.double_arguments.pop_front();
    
    MMS::ConstantConvection1D::ManufacturedSolution<dim> mms_constant_convection_1D(
            this->params.pde.reference_peclet_number,
            mms_convection_velocity,
            mms_dirichlet_value,
            mms_rate_to_steady,
            this->params.mms.initial_values_perturbation);
            
    if (this->params.mms.enabled)
    {
        assert(dim == 1); // @todo: Generalize MMS implementation
        this->mms = &mms_constant_convection_1D;
    }
    
    // Convection velocity function
    
    std::vector<double> constant_convection_velocity(dim);
    if (this->params.pde.convection_velocity_function_name == "constant")
    {
        for (unsigned int i = 0; i < dim; i++)
        {
            constant_convection_velocity[i] =
                this->params.pde.convection_velocity_function_double_arguments.front();
            this->params.pde.convection_velocity_function_double_arguments.pop_front();
        }    
    }
    
    ConstantFunction<dim> 
        constant_convection_velocity_function(constant_convection_velocity);
    
    if (this->params.pde.convection_velocity_function_name == "MMS")
    {
        assert(this->params.mms.enabled);
        assert(dim == 1);
        this->convection_velocity_function = 
            &this->mms->convection_velocity_function;
    }
    else if (this->params.pde.convection_velocity_function_name == "constant")
    {
        this->convection_velocity_function = &constant_convection_velocity_function;
    }
    else
    {
        Assert(false, ExcNotImplemented());
        // @todo: Implement ramp; shouldn't be much work
    }

    // Make initial values function
    ConstantFunction<dim> constant_function(0.);
    
    Function<dim>* initial_values_function = &constant_function;

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
        initial_values_function = &field_function;                      
    }
    else if (this->params.initial_values.function_name == "constant")
    { 
        constant_function = ConstantFunction<dim,double>(
            this->params.initial_values.function_double_arguments.front());
        initial_values_function = &constant_function;
                        
    }
    else if (this->params.initial_values.function_name == "ramp")
    {
        initial_values_function =  &ramp_function;
        
    }
    else if (this->params.initial_values.function_name == "MMS")
    { 
        assert(this->params.mms.enabled);
        initial_values_function = &this->mms->solution_function;
                        
    }
    
    // Make RHS functions
    ConstantFunction<dim> zero_source_function(0.);
    
    Function<dim>* source_function = &zero_source_function;
    
    if (this->params.mms.enabled)
    {
        source_function = &this->mms->source_function;
    }
    
    
    // Make boundary functions
    
    unsigned int boundary_count = this->params.boundary_conditions.implementation_types.size();
    
    assert(params.boundary_conditions.function_names.size() == boundary_count);

    std::vector<ConstantFunction<dim>> constant_functions;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];
        
        if ((function_name == "constant"))
        {
            double value = this->params.boundary_conditions.function_double_arguments.front();
            this->params.boundary_conditions.function_double_arguments.pop_front();
            constant_functions.push_back(ConstantFunction<dim>(value));
        }
    }
        
    // Organize boundary functions to simplify application during the time loop
    
    std::vector<Function<dim>*> boundary_functions;
    unsigned int constant_function_index = 0;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)        
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];

        if ((function_name == "constant"))
        {
            assert(constant_function_index < constant_functions.size());
            boundary_functions.push_back(&constant_functions[constant_function_index]);
            constant_function_index++;
        }
        else if ((function_name == "MMS"))
        {
            assert(this->params.mms.enabled);
            if (boundary_type == "strong")
            {
                boundary_functions.push_back(&this->mms->solution_function);
            }
            else if (boundary_type == "natural")
            {
                boundary_functions.push_back(
                    &this->mms->neumann_boundary_function);
            }

        }
        
    }