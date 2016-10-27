#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <numeric>

/*

@brief Manufactured solution for the unsteady convection-diffusion equation.

@detail

    Referencing Sandia Report SAND2000-1444 (June 2000): 
    "Code Verification by the Method of Manufactured Solutions"
*/

namespace MMS
{
    using namespace dealii;
    
    const double EPSILON = 1.e-14;

    template<int dim>
    class InitialValuesFunction : public Function<dim>
    {
    public:
        InitialValuesFunction<dim>()
        :
        Function<dim>()
        {}
        
        FunctionParser<dim>* solution_function;
        
        double perturbation;
            /*
            Allow the user to perturb the initial values.
            SAND2000-1444 says to multiply initial values by a factor 
            "not too close to one" to avoid hiding coding mistakes.
            */
        
        virtual double value(
            const Point<dim> &point,
            const unsigned int component = 0) const;
    };
    
    template<int dim>
    double InitialValuesFunction<dim>::value(
        const Point<dim> &point,
        const unsigned int /* component */) const
    {
        double u = this->solution_function->value(point);
        u *= this->perturbation;
        return u;
    }
    
    
    template <int dim>
    class BaseManufacturedSolution
    {
    public:
        BaseManufacturedSolution() {}
        
        FunctionParser<dim> solution_function;
        FunctionParser<dim> source_function;
        FunctionParser<dim> neumann_boundary_function;
        FunctionParser<dim> convection_velocity_function(3);
        
        InitialValuesFunction<dim> initial_values_function;
        
        void set_time(const double t);
        
        virtual void initialize_functions(const std::vector<double> constants);
    };

    template <int dim>
    void BaseManufacturedSolution<dim>::set_time(const double t)
    {
        this->solution_function.set_time(t);
        this->source_function.set_time(t);
        this->neumann_boundary_function.set_time(t);
    }
    
    
    namespace ConstantConvection1D
    {
        class ManufacturedSolution : public BaseManufacturedSolution<1>
        {
        public:
            ManufacturedSolution() : BaseManufacturedSolution<1>() {};
            
            virtual void initialize_functions(const std::vector<double> c);
        };
        
        void ManufacturedSolution::initialize_functions(
            const std::vector<double> c)
        {
            
            std::string variables = "x,t";
            std::string expression;
            std::map<std::string,double> constants;
            
            constants["Per"] = c[0];
            
            double convection_velocity = c[1];
            constants["a"] = convection_velocity;
            
            constants["g"] = c[2];
            constants["beta"] = c[3];
            
            this->initial_values_function.perturbation = c[4];
    
            // Solution values
            if (abs(convection_velocity) < EPSILON)
            {
                expression = "-g*((exp(-beta*t^2) - 1)*(x - 1) - 1)";
            }
            else
            {
                expression = "-g*(((exp(Per*a*x) - 1)/(exp(Per*a) - 1) - 1)*"
                    "(exp(-beta*t^2) - 1) - 1)";
            }
    
            this->solution_function.initialize(variables, expression, constants, true);
            
            // Initial values function
            this->initial_values_function.solution_function = &this->solution_function;
            
            // Source term
            if (abs(convection_velocity) < EPSILON)
            {
                expression = "2*beta*g*t*exp(-beta*t^2)*(x - 1)";
            }
            else
            {
                expression = "2*beta*g*t*exp(-beta*t^2)*((exp(Per*a*x) - 1)/"
                    "(exp(Per*a) - 1) - 1)";
            }
            
            this->source_function.initialize(variables, expression, constants, true);
            
            // Neumann boundary values
            if (abs(convection_velocity) < EPSILON)
            {
                expression = "(g*(exp(-beta*t^2) - 1))/Per";
            }
            else
            {
                expression = "(a*g*(exp(-beta*t^2) - 1))/(exp(Per*a) - 1)";
            }
            
            this->neumann_boundary_function.initialize(variables, expression, 
                constants, true);
                
            // Convection velocity
            expression = "a; 0; 0";
            this->convection_velocity_function.initialize(variables, expressions, constants);
            
        }

    }
    
    namespace VariableConvection2D
    {
        class ManufacturedSolution : public BaseManufacturedSolution<2>
        {
        public:
            ManufacturedSolution() : BaseManufacturedSolution<2>() {};
            
            virtual void initialize_functions(const std::vector<double> c);
        };
        
        void ManufacturedSolution::initialize_functions(
            const std::vector<double> c)
        {
            
            Assert(false, ExcNotImplemented());
            
        }

    }
    
}
