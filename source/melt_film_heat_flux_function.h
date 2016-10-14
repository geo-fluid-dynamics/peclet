#ifndef _melt_film_function_h_
#define _melt_film_function_h_

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary.h>

#include "pde_parameters.h"

namespace CloseContactMelting
{
    using namespace dealii;
    
    template<int dim>
    class MeltFilmHeatFluxFunction : public Function<dim>
    {
    public:
        MeltFilmHeatFluxFunction(
            const Triangulation<dim> & _tria,
            PDE::Parameters::StructuredParameters & _params)
        : Function<dim>(),
            tria(&_tria, "MeltFilmHeatFluxFunction"),
            params(&_params)
        {}
        
        virtual double value(
            const Point<dim> &point,
            const unsigned int component = 0) const;
            
    private:
        SmartPointer<const Triangulation<dim>,MeltFilmHeatFluxFunction<dim>> tria;
        PDE::Parameters::StructuredParameters * params;
        Boundary<dim> boundary;
        
    };
    
    template<int dim>
    double MeltFilmHeatFluxFunction<dim>::value(
        const Point<dim> &point,
        const unsigned int /*component*/) const
    {
       
        auto cell = GridTools::find_active_cell_around_point(*tria, point);
        assert(cell->at_boundary());
        
        Tensor<1,dim> normal;
        
        bool on_melt_boundary = false;
        unsigned int f;
        for (f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {

            if (!cell->face(f)->at_boundary())
            {
                continue;
            }
            
            if (std::find(params->boundary_conditions.melt_film.boundary_ids.begin(),
                          params->boundary_conditions.melt_film.boundary_ids.end(),
                          cell->face(f)->boundary_id())
                    < params->boundary_conditions.melt_film.boundary_ids.end())
            {
                on_melt_boundary = true;
                break;
            }
            
        }
        
        assert(on_melt_boundary);
        
        Tensor<1,dim> unit_normal;
        
        if (dim > 1)
        {
            normal = boundary.normal_vector(cell->face(f), point);            
            unit_normal = normal/normal.norm();    
        }
        else
        {
            unit_normal[0] = 1.;
        }
        
        double q_minus = -params->liquid.heat_conductivity *
                (params->solid.melt_temperature - 
                 params->boundary_conditions.melt_film.wall_temperature) /
                params->boundary_conditions.melt_film.thickness;
                
        assert(q_minus >= 0.);
        
        Tensor<1,dim> pci_velocity;
        Tensor<1,dim> melting_heat_values;
        Tensor<1,dim> qplus;
        
        for (unsigned int i = 0; i < dim; i++)
        {   
    
            pci_velocity[i] = - params->pde.convection_velocity[i]; // The velocity for the Stefan condition has the opposite sign
    
            melting_heat_values[i] = params->solid.density * 
                params->solid.latent_heat_of_melting * pci_velocity[i];    
                    
            qplus[i] =  q_minus - melting_heat_values[i];
            
        }
            
        double normal_flux = qplus*unit_normal;
        
        /*
        Convert the heat flux into a Neumann boundary conditions with units [K / s]
        */
        double neumann_bc = normal_flux / 
            (params->solid.density * params->solid.specific_heat_capacity);
            
        double unsteady_neumann_bc = params->time.time_step * neumann_bc;
            
        return unsteady_neumann_bc;
    }
     
}

#endif
