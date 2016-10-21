#include <deal.II/base/function.h>
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
 
    const double PI = numbers::PI;
    const double EPSILON = 1.e-14;
    
    
    namespace ConstantConvection1D
    {
        
        template <int dim>
        class BaseFunction : public Function<dim>
        {
        public:
            BaseFunction
                (
                double _reference_peclet_number = 1.,
                double _convection_velocity = -1.,
                double _dirichlet_value = -1.,
                double _rate_to_steady = 10.
                )
                : 
                Function<dim>(),
                reference_peclet_number(_reference_peclet_number),
                convection_velocity(_convection_velocity),
                dirichlet_value(_dirichlet_value),
                rate_to_steady(_rate_to_steady)
            {}
        
        protected:
            const double reference_peclet_number;
            const double convection_velocity;
            const double dirichlet_value;
            const double rate_to_steady;
        
        };
        
        template <int dim>
        class SolutionFunction : public BaseFunction<dim>
        {
        public:
            
            SolutionFunction
                (
                const double _reference_peclet_number = 1.,
                const double _convection_velocity = -1.,
                const double _dirichlet_value = -1.,
                const double _rate_to_steady = 10.,
                double _initial_values_perturbation = 1.01
                )
                :
                BaseFunction<dim>
                    (
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady
                    ),
                initial_values_perturbation(_initial_values_perturbation)
            {}
            
            virtual double value(
                const Point<dim> &point,
                const unsigned int component = 0) const;
                
        private:
            const double initial_values_perturbation;
                /*
                Allow the user to perturb the initial values.
                SAND2000-1444 says to multiply initial values by a factor 
                "not too close to one" to avoid hiding coding mistakes.
                */
        };

        template <int dim>
        double SolutionFunction<dim>::value
            (
            const Point<dim> &point,
            const unsigned int /* component */
            ) const
        {
            const double x = point[0];
            const double t = this->get_time();
            const double a = this->convection_velocity;
            const double Pe = this->reference_peclet_number;
            const double g = this->dirichlet_value;
            const double beta = this->rate_to_steady;
            
            double u = g*(1. + (exp(a*Pe) - exp(a*Pe*x))/
                (1. - exp(a*Pe))*(1. - exp(-beta*t*t)));
            
            if (t < EPSILON)
            {
                u *= this->initial_values_perturbation;    
            }
            
            return u;
        }
        
        template <int dim>
        class SourceFunction : public BaseFunction<dim>
        {
        public:
            SourceFunction
                (
                const double _reference_peclet_number = 1.,
                const double _convection_velocity = -1.,
                const double _dirichlet_value = -1.,
                const double _rate_to_steady = 10.
                )
                :
                BaseFunction<dim>
                    (
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady
                    )
            {}
        
            virtual double value(
                const Point<dim> &point,
                const unsigned int component = 0) const;
        };
        
        template <int dim>
        double SourceFunction<dim>::value
            (
            const Point<dim> &point,
            const unsigned int /* component */
            ) const
        {
            const double x = point[0];
            const double t = this->get_time();
            const double a = this->convection_velocity;
            const double Pe = this->reference_peclet_number;
            const double g = this->dirichlet_value;
            const double beta = this->rate_to_steady;
            
            double s = g/(exp(a*Pe) - 1.)*
                (
                  a*exp(a*Pe*x)*(exp(-beta*t*t) - 1.)*(1. - a*Pe) 
                  - 2.*beta*t*exp(-beta*t*t)*(exp(a*Pe) + exp(a*Pe*x))
                );
            
            return s;
        }
        
        template <int dim>
        class NeumannBoundaryFunction : public BaseFunction<dim>
        {
        public:
            NeumannBoundaryFunction
                (
                const double _reference_peclet_number = 1.,
                const double _convection_velocity = -1.,
                const double _dirichlet_value = -1.,
                const double _rate_to_steady = 10.
                )
                :
                BaseFunction<dim>
                    (
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady
                    )
            {}
            
            virtual double value(
                const Point<dim> &point,
                const unsigned int component = 0) const;
        };
        
        template <int dim>
        double NeumannBoundaryFunction<dim>::value
            (
            const Point<dim> &point,
            const unsigned int /* component */
            ) const
        {
            const double x = point[0];
            assert(x < EPSILON);
            
            const double t = this->get_time();
            const double a = this->convection_velocity;
            const double Pe = this->reference_peclet_number;
            const double g = this->dirichlet_value;
            const double beta = this->rate_to_steady;
            
            double h = a*Pe*g*(1. - exp(-beta*t*t))/(1. - exp(Pe*a));
            
            return h;
        }
        
        template <int dim>
        class ManufacturedSolution
        {
        public:
            ManufacturedSolution
                (
                const double _reference_peclet_number = 1.,
                const double _convection_velocity = -1.,
                const double _dirichlet_value = -1.,
                const double _rate_to_steady = 10.,
                const double _initial_values_perturbation = 1.01
                )
            :  
            
                solution_function(SolutionFunction<dim>(
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady,
                    _initial_values_perturbation)),
                    
                source_function(SourceFunction<dim>(
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady)),
                    
                neumann_boundary_function(NeumannBoundaryFunction<dim>(
                    _reference_peclet_number,
                    _convection_velocity,
                    _dirichlet_value,
                    _rate_to_steady)),
                    
                convection_velocity_function(ConstantFunction<dim>(_convection_velocity))
            {}  
            
            SolutionFunction<dim> solution_function;
            SourceFunction<dim> source_function;
            NeumannBoundaryFunction<dim> neumann_boundary_function;
            
            ConstantFunction<dim> convection_velocity_function;
            
            void set_time(double t);
            
        };
    
        template <int dim>
        void ManufacturedSolution<dim>::set_time(double t)
        {
            this->solution_function.set_time(t);
            this->source_function.set_time(t);
            this->neumann_boundary_function.set_time(t);
        }
        
    }
    
}
