#include <deal.II/base/function.h>
#include <numeric>

namespace MyFunctions
{
    using namespace dealii;
    template<int dim>
    class RampFunctionAlongLine : public Function<dim>
    {
    public:
        RampFunctionAlongLine
            (
            const Point<dim> _line_start,
            const Point<dim> _line_end,
            const double _ramp_start_position,
            const double _ramp_end_position,
            const double _start_value,
            const double _end_value
            )
          : Function<dim>(),
            line_start(_line_start),
            line_end(_line_end),
            ramp_start_position(_ramp_start_position),
            ramp_end_position(_ramp_end_position),
            start_value(_start_value),
            end_value(_end_value),
            direction((_line_end - _line_start)/(_line_end - _line_start).norm()),
            length((_line_end - _line_start).norm()
            )
        {
            assert(_ramp_start_position >= 0.);
            assert(_ramp_end_position <= 1.);
            assert(_ramp_end_position >= ramp_start_position);
        }
        virtual double value(const Point<dim>  &point,
                             const unsigned int component = 0) const;
    private:
        const Point<dim> line_start;
        const Point<dim> line_end;
        const double ramp_start_position;
        const double ramp_end_position;
        const double start_value;
        const double end_value;
        const Point<dim> direction;
        const double length;
    };

    template<int dim>
    double RampFunctionAlongLine<dim>::value
        (
        const Point<dim> &point,
        const unsigned int // component
        ) const
    {
        //Assert(component == 0, ExcInternalError());
        double val;
        double epsilon = 1e-8;
        // Get distance along line
        // Note: This method projects the point onto the line,
        //       and does nothing to guarantee that the point is actually "near" the line.
        Point<dim> position_vector;
        double inner_product = 0;
        for (int i = 0; i < dim; i++)
        {
            position_vector[i] = point[i] - this->line_start[i];
            inner_product += position_vector[i]*this->direction[i];
        }
        double parametric_position = inner_product/this->length;
        assert(parametric_position >= 0 - epsilon);
        assert(parametric_position <= 1 + epsilon);
        // Return an end value if not on the ramp
        if (parametric_position <= ramp_start_position)
        {
            val = start_value;
        }
        else if (parametric_position >= ramp_end_position)
        {
            val = end_value;
        }
        // Interpolate if on ramp
        else
        {
            val = start_value + (parametric_position - ramp_start_position)*
                  (end_value - start_value)/(ramp_end_position - ramp_start_position);
        }
        return val;
    }
    
}
