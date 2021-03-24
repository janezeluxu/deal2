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
#include "../include/includeall.h"
#include "../include/NavierStoke.h"
//#include "../include/generalizeAlpha.h"
//#include "../include/linear.h"
//#include "../include/newton.h"
#include "../include/pre.h"
#include "../include/user_input.h"
#include "../include/post.h"

int main()
{
    
  try
    {
        using namespace dealii;
        using namespace incompressible;
        
        //Triangulation<2> triangulation;
        //DoFHandler<2> dof_handler(triangulation);
        //extern int degree;
        
        NavierStokes flow(degree,mapdegree);
        flow.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
