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
#include "../include/user_input.h"
using namespace std;
namespace incompressible
{
    using namespace dealii;

    double Re = 100;
    double dt = 0.125;
    
    double ri = 0.25;
    double alpham = 0.5*(3- ri)/(1+ri);
    double alphaf = 1/(1+ri);
    double gamma = 0.5+alpham-alphaf;
    
    //double gamma = 1;
    //double alphaf = 1;
    //double alpham = 1;
    
    int maxNewtonIter = 10;
    double tolerance = 1e-6;
    
    int degree = 1;
    int mapdegree = 1;
    int meshRefineLevel = 2;
    unsigned int nstp = 820;
    unsigned int num_start = 0;
    
    //dynamic info
    int surfaceTag = 1;
    //int dynamicBC = 50;
    
    bool dynamicFlag = false;
    int dynamicBC = 50;
    double c1_start = 0.3;
    int stage1_outter = 1;
    int stage1_inner = 1;
    int stage2 = 1;

    void defineMesh(std::vector<unsigned int> &subdivisions,double &x1,double &x2,double &y1,double &y2)
    {
        subdivisions[0] = 40;
        subdivisions[1] = 40;
        
        x1=0.0;x2=1.0;y1=0.0;y2=1.0;
    }
    
    
    double Dirichlet_u(const dealii::Point<2> & p,int BC_ID,double t)
    {
        if (BC_ID==6) //left
            return 1;//+0.001*p[1];
        else if (BC_ID==50) // surface
        {
            return 0;
        }
        return 0;
    }
    
    double Dirichlet_v(const dealii::Point<2> & p,int BC_ID,double t)
    {
        if (BC_ID==6) //left
            return 0;
        else if (BC_ID==24) // top
            return 0;
        else if (BC_ID==23) //buttom
            return 0;
        else if (BC_ID==50) //surface
            return 0;
        return 0;
    }
    
    double Dirichlet_p(const dealii::Point<2> & p,int BC_ID,double t)
    {
        if (BC_ID==14) //right
            return 0;
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


