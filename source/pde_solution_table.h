#include <deal.II/base/table_handler.h>

template<int dim>
void Model<dim>::append_solution_to_table()
{
    /*
    I haven't yet figured out how to simply get the solution at each vertex.
    The following method is going to write duplicates.
    */
    
    for (auto cell = this->dof_handler.begin_active(); cell != this->dof_handler.end(); ++cell)
    {
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
            this->solution_table.add_value("t", time);
            
            for (unsigned int axis = 0; axis < dim; axis++)
            {
                this->solution_table.add_value("x"+std::to_string(axis),
                    cell->vertex(v)[axis]);
            }
            
            for (unsigned int dof = 0; dof < fe.dofs_per_vertex; dof++)
            {
                this->solution_table.add_value("u"+std::to_string(dof),
                    this->solution(cell->vertex_dof_index(v, 0)));
            }
            
        }
    }

}

template<int dim>
void Model<dim>::write_solution_table(const std::string file_name)
{
    const int precision = 14;
    this->solution_table.set_precision("t", precision);
    this->solution_table.set_scientific("t", true);
    for (unsigned int axis = 0; axis < dim; axis++)
    {
        this->solution_table.set_precision("x"+std::to_string(axis), precision);
        this->solution_table.set_scientific ("x"+std::to_string(axis), true);
    }
    for (unsigned int dof = 0; dof < this->fe.dofs_per_vertex; dof++)
    {
        this->solution_table.set_precision("u"+std::to_string(dof), precision);
        this->solution_table.set_scientific ("u"+std::to_string(dof), true);
    }

    std::ofstream out_file(file_name, std::fstream::app);
    assert(out_file.good());
    this->solution_table.write_text(out_file);
    out_file.close();        
}
