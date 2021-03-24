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

// As usual, we start by including some well-known files:

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

//#include "../include/readin.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <string.h>
using namespace std;
namespace incompressible
{
    using namespace dealii;

    double Re = 40;
    double gamma = 1;
    double dt = 1;
    double alphaf = 1;
    double alpham = 1;
    int maxNewtonIter = 20;
    double tolerance = 1e-16;
    
    int degree = 1;
    int meshRefineLevel = 3;
    unsigned int nstp = 1;
    
    void defineMesh(std::vector<unsigned int> &subdivisions,double &x1,double &x2,double &y1,double &y2)
    {
        subdivisions[0] = 40;
        subdivisions[1] = 40;
        
        x1=0.0;x2=1.0;y1=0.0;y2=1.0;
    }
    
    
    double Dirichlet_u(const dealii::Point<2> & p,int BC_ID)
    {
        double u_out;
        double u_inlet;
        
        double pi = M_PI;
        double miu = 1/Re;
        double x = p[0];
        double y = p[1];
        double lam = 0.5*Re-sqrt((0.25*Re*Re+4*pi*pi));
        u_inlet =  (1-exp(lam*x)*cos(2*pi*y));
        
        u_out =  (1-exp(lam*x)*cos(2*pi*y));
        
        
        if (BC_ID==3) //left
            return u_inlet;
        else if (BC_ID==4) //right
            return u_out;
        return 0;
    }
    
    double Dirichlet_v(const dealii::Point<2> & p,int BC_ID)
    {
        double v_inlet;
        double v_out;
        double x = p[0];
        double y = p[1];
        double pi = M_PI;
        double lam = 0.5*Re-sqrt((0.25*Re*Re+4*pi*pi));
        v_inlet =  ((lam/(2*pi))*exp(lam*x)*sin(2*pi*y));
        v_out =  (lam/(2*pi))*exp(lam*x)*sin(2*pi*y);
        if (BC_ID==1) //top
            return 0;
        else if (BC_ID==2) //buttom
            return 0;
        else if (BC_ID==3) //left
            return v_inlet;
        else if (BC_ID==4) //right
            return v_out;
        
        return 0;
    }
    
    double Dirichlet_p(const dealii::Point<2> & p,int BC_ID)
    {
        double pi = M_PI;
        double lam = 0.5*Re-sqrt((0.25*Re*Re+4*pi*pi));
        double p_out;
        double x = p[0];
        p_out = (0.5*(1-exp(2*lam*x)));
        if (BC_ID==4) //right
            return p_out;
        
        return 0;
    }
    
    void diff_flux(const Point<2> & p1, const Point<2> & p2,
                 bool &ifconsist, double &valuex, double &valuey)
    {
        ifconsist = false;
        valuex = 0;
        valuey = 0;
    }
    
    void pressure_flux(const Point<2> & p1, const Point<2> & p2,
                   bool &ifconsist, double &valuex, double &valuey)
    {
        ifconsist = false;
        valuex = 0;
        valuey = 0;
    }
    
}


