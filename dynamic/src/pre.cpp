/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: zelu Xu, 2018
 */

// @sect3{Include files}



#include "../include/includeall.h"
#include "../include/pre.h"
#include "../include/user_input.h"
namespace incompressible
{
    using namespace dealii;
    void createmesh(Triangulation<2> &triangulation,MappingQ<2> mapping,int n_refines)
    {
        // create a 2d rectangle
        /*
         GridGenerator::hyper_cube(triangulation);
         triangulation.refine_global(RefineLevel);
         std::ofstream out("grid-1.eps");
         GridOut       grid_out;
         grid_out.write_eps(triangulation, out);
         std::cout << "Grid written to grid-1.eps" << std::endl;
         */
        
        /*
        std::vector<unsigned int> subdivisions (2, 2);
        double x1;double x2;double y1;double y2;
        
        defineMesh(subdivisions, x1,x2,y1,y2);
        
        const Point<2> bottom_left= (Point<2>(x1,y1));
        const Point<2> top_right= (Point<2>(x2,y2));
        
        GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                   subdivisions,
                                                   bottom_left,
                                                   top_right);
        std::ofstream out("grid-1.eps");
        GridOut       grid_out;
        grid_out.write_eps(triangulation, out);
        std::cout << "Grid written to grid-1.eps" << std::endl;
        
        
        // assign boundary ids
        Triangulation<2>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        double xcord; double ycord;
        for (; cell!=endc; ++cell)
            for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
                if (cell->face(face)->at_boundary())
                {
                    
                    xcord = cell->face(face)->center()[0];
                    ycord = cell->face(face)->center()[1];
                    if (ycord==y1 ) //y=0 pressure boundary
                    {
                        cell->face(face)->set_boundary_id(1);
                    }
                    
                    //else if (ycord==y2) //y=0
                    //{cell->face(face)->set_boundary_id(BC_ID[2]);}
                    //else if (xcord==x1 )//x=0
                    //{cell->face(face)->set_boundary_id(BC_ID[2]);}
                    //else if (xcord==x2 ) //x=1
                    //{cell->face(face)->set_boundary_id(BC_ID[3]);}
                    
                    else //velocity boundary
                    {cell->face(face)->set_boundary_id(0);}
                }
        
        */
        GridIn<2> grid_in;
        grid_in.attach_triangulation(triangulation);
        {
            //std::string   filename = "cylinderquadcc.inp";
            std::string   filename = "./mesh/gmsh/squarecylinder.msh";
            std::ifstream file(filename.c_str());
            Assert(file, ExcFileNotOpen(filename.c_str()));
            //grid_in.read_ucd(file);
            grid_in.read_msh(file);
        }
        //std::cout << "Number of refines = " << n_refines << std::endl;
        
        //int n_refines = 3;
        triangulation.refine_global(n_refines);
        std::vector<types::boundary_id>           boundary_ids;
        boundary_ids = triangulation.get_boundary_ids();
        
        std::ofstream out("grid-flowovercylinder.eps");
        GridOut       grid_out;
        GridOutFlags::Gnuplot gnuplot_flags(false, 60);
        grid_out.set_flags(gnuplot_flags);
        
        grid_out.write_eps(triangulation, out, &mapping);
        std::cout << "Grid written to grid-flowoversquare.eps" << std::endl;
        
        
    }
    
    
    double Dirichlet_BC::value(const Point<2> & p, const unsigned int component) const
    {
        //std::cout<<"BC_id "<<boundaryID<<" point "<<p[0]<<" "<<p[1]<<" component "<<component<<std::endl;
        
        Assert(component < this->n_components,
               ExcIndexRange(component, 0, this->n_components));
        
        //if (component == 0 && std::abs(p[1]-1.0) < 1e-10)
        //    return 1.0;
        
        if (component == 0 )
            return Dirichlet_u(p,boundaryID,t);
        else if (component == 1 )
            return Dirichlet_v(p,boundaryID,t);
        else if (component == 2 )
            return Dirichlet_p(p,boundaryID,t);
        else
            std::cout<<"which component?? "<<std::endl;
        
        return 0;
    }
    
    void Dirichlet_BC::vector_value(const Point<2> &p,
                                      Vector<double> &  values) const
    {
        for (unsigned int c = 0; c < this->n_components; ++c)
            values(c) = Dirichlet_BC::value(p, c);
    }
    

    /*
    double BoundaryValues::value(const Point<2> &p,
                                      const unsigned int component) const
    {
        Assert (component < this->n_components,
                ExcIndexRange (component, 0, this->n_components));
        if (component == 0 && std::abs(p[1]-1.0) < 1e-10)
            return 1.0;
        
        return 0;
    }
    
    void BoundaryValues::vector_value ( const Point<2> &p,
                                            Vector<double>   &values ) const
    {
        for (unsigned int c = 0; c < this->n_components; ++c)
            values(c) = BoundaryValues::value (p, c);
    }
    */
    void getflux(const Point<2> & p1, const Point<2> & p2, char flux_type,
                 bool &ifconsist, double &valuex, double &valuey )
    {
        // Check index of the current boundaryID
        //std::vector<int>::iterator it = std::find(BC_ID.begin(), BC_ID.end(), boundaryID);
        
        if (flux_type == 'd')
        {
            diff_flux(p1,p2,ifconsist,valuex,valuey);
        }
        
        else if (flux_type == 'p')
        {
            pressure_flux(p1,p2,ifconsist,valuex,valuey);
        }
        else
            std::cout << " what flux? " << '\n';
        
    }
}


