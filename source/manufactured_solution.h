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
    
    
    template<int dim>
    class ManufacturedSolution : public Function<dim>
    {
    public:
        ManufacturedSolution() : Function<dim>() {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    };

    template<int dim>
    double ManufacturedSolution<dim>::value
        (
        const Point<dim> &x,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        
        double product = 1;
        for (unsigned int i = 0; i < dim; i++)
        {
            double s = sin((i + 1.)*pi*x[i]);
            product *= s*s;
        }
        
        double u = exp(-10.*t*t)*(1. + product)
        return u;
    }
    
    
    template<int dim>
    class PecletNumber : public Function<dim>
    {
    public:
        PecletNumber(double _max_value) : Function<dim>(), max_value(_max_value) {}
        virtual double value(const Point<dim> &x,
                             const unsigned int component = 0) const;
    private:
        double max_value;
    };

    template<int dim>
    double PecletNumbern<dim>::value
        (
        const Point<dim> /* &x */,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        
        double Pe = 0.5*this->max_value*(1. + sin(-pi/2. + pi*t))
        return Pe;
    }
    
    
    template<int dim>
    class AdvectionVelocity : public Function<dim>
    {
    public:
        AdvectionVelocity() : Function<dim>() {}
        virtual double value(const Point<dim> &x,
                             const unsigned int component = 0) const;
    };

    template<int dim>
    double AdvectionVelocity<dim>::value
        (
        const Point<dim> /* &x */,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        
        Tensor<dim, 1> a;
        
        a[0] = cos(2.*pi*t);
        if (dim > 1)
        {
            a[1] = sin(2.*pi*t);
        }
        if (dim > 2)
        {
            a[2] = 0.;
        }
        
        return a;
    }
    
    
    template<int dim>
    class InitialValues : public Function<dim>
    {
    public:
        InitialValues(double _perturbation) : Function<dim>(), perturbation(_perturbation) {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    private:
        ManufacturedSolution<dim> u;
        double perturbation; // SAND2000-1444 says to multiply initial values by a factor "not too close to one" to avoid hiding coding mistakes.
    };

    template<int dim>
    double InitialValues<dim>::value
        (
        const Point<dim> &x,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        this->u.set_time(t);

        return perturbation*this->u(x);
    }
    
    
    template<int dim>
    class Source : public Function<dim>
    {
    public:
        Source(double _max_Pe) : Function<dim>(), max_Pe(_max_Pe) {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    private:
        double max_Pe;
    };

    template<int dim>
    double Source<dim>::value
        (
        const Point<dim> &point,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        double x, y, z;
        x = point[0];
        if (dim > 1)
        {
            y = point[1];
        }
        if (dim > 2)
        {
            z = point[2];
        }
        
        double e = exp(-10.*t*t);
        double pi2e = pi*pi*e, pit = pi*t, pix = pi*x, twopiy = 2*pi*y, threepiz = 3*pi*z;
        double sinpix2 = sin(pix)*sin(pix), sin2piy2 = sin(twopiy)*sin(twopiy), 
            sin3piz2 = sin(threepiz)*sin(threepiz);
        
        double s = 
            (2.*(2.*pi2e*cos(pix)*cos(pix)*sintwopiy2*sinthreepiz2 
            + 8.*pi2e*cos(twopiy)*cos(twopiy)*sinpix2*sin3piz2 
            + 18.*pi2e*cos(threepiz)*cos(threepiz)*sinpix2*sintwopiy2 
            - 28.*pi2e*sinpix2*sintwopiy2*sin3piz2))/(this->max_Pe*(sin(pi/2. - pit) - 1.))
            - 20.*t*e*(sinpix2*sintwopiy2*sin3piz2 + 1.) 
            + 2.*pi*e*cos(2.*pit)*cos(pix)*sin(pix)*sintwopiy2*sinthreepiz2 
            + 4.*pi*e*cos(twopiy)*sin(2.*pit)*sinpix2*sin(twopiy)*sinthreepiz2;
        
        return s;
    }
    
    
    template<int dim>
    class NeumannBoundary : public Function<dim>
    {
    public:
        NeumannBoundary(
            const Triangulation<dim> & _tria,
            const Boundary<dim> _Gamma,
            const double _max_Pe) 
        : Function<dim>(),
            tria(&_tria, "MeltFilmHeatFluxFunction"),
            Gamma(_Gamma),
            max_Pe(_max_Pe)
        {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    private:
        SmartPointer<const Triangulation<dim>,NeumannBoundary<dim>> tria;
        Boundary<dim> Gamma;
        double max_Pe;
    };

    template<int dim>
    double NeumannBoundary<dim>::value
        (
        const Point<dim> &point,
        const unsigned int /* component */
        ) const
    {
               
        /*
        This is a terribly inefficient approach for obtaining the normal vector,
        but I do not know any other non-invasive way to do it. By non-invasive I mean that this function
        should have no knowledge of the quadrature points or fe_face_values, since we want to be able
        to simply hand this function as an argument to high level methods that build system matrices or 
        right-hand-sides.
        */
        auto cell = GridTools::find_active_cell_around_point(*tria, point);
        assert(cell->at_boundary());
        
        Tensor<dim, 1> n;
        
        unsigned int f;
        for (f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {

            if (cell->face(f)->at_boundary())
            {
                this->Gamma.normal_vector(cell->face(f), point);
                break;
            }
            
        }
        
        n = n/n.norm(); // We want the unit normal
        
        
        double t = this->get_time();
        double x, y, z;
        x = point[0];
        if (dim > 1)
        {
            y = point[1];
        }
        if (dim > 2)
        {
            z = point[2];
        }
        
        
        double e = exp(-10.*t*t);
        double pie = pi*e, pit = pi*t, pix = pi*x, twopiy = 2.*pi*y, threepiz = 3.*pi*z;
        double sinpix = sin(pi*x), sin2piy = sin(2.*pi*y), sin3piz = sin(3.*pi*z);
        double sinpix2 = sin(pix)*sin(pix), sin2piy2 = sin(twopiy)*sin(twopiy), 
            sin3piz2 = sin(threepiz)*sin(threepiz);
        
        double h = 
            -(2.*(2.*n[0]*pie)*cos(pi*x)*sinpix*sin2piy*sin2piy*sin3piz*sin3piz
            + 4.*n[1]*pie*cos(2.*pi*y)*sinpix*sinpix*sin2piy*sin3piz*sin3piz
            + 6.*n[2]*pie*cos(3.*pi*z)*sinpix*sinpix*sin2piy*sin2piy*sin3piz))
            /(this->max_Pe*(sin(pi/2. - pi*t) - 1.));
        
        return h;
    }
    
}
