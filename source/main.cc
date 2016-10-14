#include "pde_model.h"
int main(int argc, char* argv[])
{
    try
    {   
        /*
        The dimensionality of the model can be specified in the parameter file,
        so we instantiate all possibilities and read the parameters into each
        possible model.
        
        Nothing is assembled until executing model.run(), so the extra models are 
        virtually free.
        */
        
        PDE::Model<1> model_1D;
        PDE::Model<2> model_2D;
        
        if (argc == 2)
        {
            std::string parameter_input_file_path = argv[1];
            model_1D.read_parameters(parameter_input_file_path);
            model_2D.read_parameters(parameter_input_file_path);
        }
        
        /*
        At this line, all three models have the same parameters; but only one of the
        instances is valid.
        */
        
        const unsigned int dim = model_1D.params.geometry.dim;
        switch (dim)
        {
            case 1:
                model_1D.run();
                break;
            case 2:
                model_2D.run();
                break;
        }

    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
        std::cerr << "Exception on processing: " << std::endl << exc.what()
              << std::endl << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
        std::cerr << "Unknown exception!" << std::endl << "Aborting!"
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
        return 1;
    }
    return 0;
}