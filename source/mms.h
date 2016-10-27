#ifndef mms_h
#define mms_h

/*
@brief Methods helpful for MMS (Method of Manufactured Solutions) verification

@detail
    Referencing Sandia Report SAND2000-1444 (June 2000): 
    "Code Verification by the Method of Manufactured Solutions"
    
@author A. Zimmerman <zimmerman@aices.rwth-aachen.de>
*/

namespace MMS
{
    template <int dim>
    class InitialValuesFunction : public Function<dim>
    {
    public:
        InitialValuesFunction(Function<dim>* _solution_function) : Function<dim>(), 
            solution_function(_solution_function),
            perturbation(1. + 1.e-9)
        {}
        Function<dim>* solution_function;
        double perturbation;
            /*
            Allow the user to perturb the initial values.
            SAND2000-1444 says to multiply initial values by a factor 
            "not too close to one" to avoid hiding coding mistakes.
            */
        
        virtual double value(
            const Point<dim>  &point,
            const unsigned int component = 0) const;
    };
    
    template<int dim>
    double InitialValuesFunction<dim>::value(
        const Point<dim> &point,
        const unsigned int component) const
    {
        double val = this->solution_function->value(point, component);
        val *= this->perturbation;
        return val;
    }
    
}

#endif
