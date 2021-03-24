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

#include "../include/generalizeAlpha.h"
#include "../include/newton.h"
#include "../include/user_input.h"
#include <fstream>
#include <iostream>
#include <sstream>
namespace incompressible
{
    using namespace dealii;
    
    generalizeAlpha::generalizeAlpha(double viscosity)
    :
    viscosity(viscosity)
    {}
    
    void generalizeAlpha::initialize(hp::DoFHandler<2> &dof_handler,ConstraintMatrix &nonzero_constraints, std::vector<types::global_dof_index> &dofs_per_block)
    {
        system_matrix.clear();
        //pressure_mass_matrix.clear();
        {
            BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
            DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
            sparsity_pattern.copy_from(dsp);
        }
        system_matrix.reinit(sparsity_pattern);
        
        newton_update.reinit(dofs_per_block);
        system_rhs.reinit(dofs_per_block);
        
    }
    
    
    void generalizeAlpha::corrector(int maxNewtonIter,double tolerance, hp::DoFHandler<2> &dof_handler,hp::FECollection<2> fe,ConstraintMatrix &zero_constraints, ConstraintMatrix &nonzero_constraints,BlockVector<double> &n_solution)
    {
        
        // newtons loop
        
        double current_res=0;
        //double last_res;

        NewtonSolve newtonSystem(viscosity,
                                 dof_handler,fe,
                                 sparsity_pattern,
                                 zero_constraints,
                                 n_solution,
                                 system_matrix,system_rhs);
        
        nonzero_constraints.distribute(n_solution);
        
        for(int i=0; i<maxNewtonIter; ++i)
        {
            /*
            Vector<double> du_n = n_solution.block(0);
            Vector<double> dp_n = n_solution.block(1);
            std::cout << "n_solution init " <<std::endl;
            for(std::size_t i=0; i<du_n.size(); ++i)
            {
                std::cout<<du_n[i]<<std::endl;
            }
            
            for(std::size_t i=0; i<dp_n.size(); ++i)
            {
                std::cout<<dp_n[i]<<std::endl;
            }
            */
            std::cout << "----start iteration Number ----" <<i<< std::endl;
            newtonSystem.getNewtonUpdate(n_solution,newton_update);
            
            /*
             du_n = newton_update.block(0);
             dp_n = newton_update.block(1);
            std::cout << "newton_update " <<std::endl;
            for(std::size_t i=0; i<du_n.size(); ++i)
            {
                std::cout<<du_n[i]<<std::endl;
            }
            
            for(std::size_t i=0; i<dp_n.size(); ++i)
            {
                std::cout<<dp_n[i]<<std::endl;
            }
            */
            std::cout << "-------now update flow-------" << std::endl;
            newtonSystem.update_flow(newton_update,n_solution);
            
            newtonSystem.getResidual(n_solution,current_res);
            
            
            if (current_res < tolerance)
                break;
        }
        
        
    }
    
    
}


