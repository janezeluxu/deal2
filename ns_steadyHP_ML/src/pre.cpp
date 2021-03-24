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
#include <fstream>
#include <iostream>
#include <sstream>
namespace incompressible
{
    using namespace dealii;
    void createmesh(Triangulation<2> &triangulation)
    {

        GridIn<2> gridin;
        gridin.attach_triangulation(triangulation);
        
        //const std::__1::string meshfilename = "./mesh/lidCavity001.msh";
        std::ifstream f(meshfile);
        //std::ifstream f("./mesh/lidCavity001.msh");
        gridin.read_msh(f);
        //gridin.read_ucd(f);

        std::vector<types::boundary_id>           boundary_ids;
        boundary_ids = triangulation.get_boundary_ids();
        
        std::ofstream out("grid-lidcavity.eps");
        GridOut       grid_out;
        grid_out.write_eps(triangulation, out);
        std::cout << "Grid written to grid-lidcavity.eps" << std::endl;
        
    }
    
    
    double Dirichlet_BC::value(const Point<2> & p, const unsigned int component) const
    {
        //std::cout<<"BC_id "<<boundaryID<<" point "<<p[0]<<" "<<p[1]<<" component "<<component<<std::endl;
        
        Assert(component < this->n_components,
               ExcIndexRange(component, 0, this->n_components));
        
        //if (component == 0 && std::abs(p[1]-1.0) < 1e-10)
        //    return 1.0;
        
        if (component == 0 )
            return Dirichlet_u(p,boundaryID);
        else if (component == 1 )
            return Dirichlet_v(p,boundaryID);
        else if (component == 2 )
            return Dirichlet_p(p,boundaryID);
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


