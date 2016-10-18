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
    class ConvectionVelocity : public Function<dim>
    {
    public:
        ConvectionVelocity() : Function<dim>() {}
        virtual double value(const Point<dim> &x,
                             const unsigned int component = 0) const;
    };

    template<int dim>
    double ConvectionVelocity<dim>::value
        (
        const Point<dim> &x,
        const unsigned int /* component */
        ) const
    {
        
        double t = this->get_time();
        
        Tensor<dim, 1> a;
        
        a[0] = x[0];
        if (dim > 1)
        {
            a[1] = x[1]/2.;
        }
        if (dim > 2)
        {
            a[2] = x[2]/3.;
        }
        
        return t*a;
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
        Source(double _max_Pe) : Function<dim>(), Pe_r(_Pe_r) {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    private:
        double Pe_r;
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
        double pie = pie, pi2e = pi*pie, pit = pi*t, pix = pi*x,
            twopiy = 2.*pi*y, threepiz = 3.*pi*z;
        double sinpix2 = sin(pix)*sin(pix), sin2piy2 = sin(twopiy)*sin(twopiy), 
            sin3piz2 = sin(threepiz)*sin(threepiz);
        double cos2piy = cos(2.*pi*y), cos3piz = cos(3.*pi*z);
        
        double s = 
            2.*t*x*pie*cos(pi*x)*sin(pi*x)*sin2piy2*sin3piz2 
            - 20.*t*e*(sinpix2*sin2piy2*sin3piz2 + 1.)
             - (2.*pi2e*cos(pi*x)*cos(pi*x)*sin2piy2*sin3piz2 
            + 8.*pi2e*cos2piy*cos2piy*sinpix2*sin3piz2
             + 18.*pi2e*cos3piz*cos3piz*sinpix2*sin2piy2 
            - 28.*pi2e*sinpix2*sin2piy2*sin3piz2)/this->Pe_r 
            + 2.*t*y*pie*cos2piy*sinpix2*sin(twopiy)*sin3piz2 
            + 2.*t*z*pie*cos3piz*sinpix2*sin2piy2*sin(threepiz);
        
        return s;
    }
    
    
    template<int dim>
    class NeumannBoundary : public Function<dim>
    {
    public:
        NeumannBoundary(
            const Triangulation<dim> & _tria,
            const Boundary<dim> _Gamma,
            const double _Pe_r) 
        : Function<dim>(),
            tria(&_tria, "NeumannBoundary"),
            Gamma(_Gamma),
            Pe_r(_Pe_r)
        {}
        virtual double value(const Point<dim>  &x,
                             const unsigned int component = 0) const;
    private:
        SmartPointer<const Triangulation<dim>,NeumannBoundary<dim>> tria;
        Boundary<dim> Gamma;
        double Pe_r;
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
        double pie = pie, pit = pi*t, pix = pi*x, twopiy = 2.*pi*y, threepiz = 3.*pi*z;
        double sinpix = sinpix, sin2piy = sin(2.*pi*y), sin3piz = sin(3.*pi*z);
        double sinpix2 = sin(pix)*sin(pix), sin2piy2 = sin(twopiy)*sin(twopiy), 
            sin3piz2 = sin(threepiz)*sin(threepiz);
        
        double h = 
            (2.*n0*pie*cos(pi*x)*sin(pi*x)*sin2piy2*sin3piz2
            + 4.*n1*pie*cos(twopiy)*sinpix2*sin(twopiy)*sin3piz2
            + 6.*n2*pie*cos(threepiz)*sinpix2*sin2piy2*sin(threepiz))/this->Pe_r;
        
        return h;
    }
    
}
