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
namespace incompressible
{
    using namespace dealii;
    
    generalizeAlpha::generalizeAlpha(int degree, double viscosity, double gamma, double dt,double alphaf, double alpham)
    :
    degree(degree),
    viscosity(viscosity),
    gamma(gamma),
    dt(dt),
    alphaf(alphaf),
    alpham(alpham)
    {}
    
    void generalizeAlpha::initialize(DoFHandler<2> &dof_handler,ConstraintMatrix &nonzero_constraints, std::vector<types::global_dof_index> &dofs_per_block)
    {
        system_matrix.clear();
        //pressure_mass_matrix.clear();
        {
            BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
            DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
            sparsity_pattern.copy_from(dsp);
        }
        system_matrix.reinit(sparsity_pattern);
        
        //n_solution.reinit(dofs_per_block);
        //n_solution_time_derivative.reinit(dofs_per_block);
        //np1_solution.reinit(dofs_per_block);
        //np1_solution_time_derivative.reinit(dofs_per_block);
        
        npaf_solution.reinit(dofs_per_block);
        npam_solution_time_derivative.reinit(dofs_per_block);
        
        newton_update.reinit(dofs_per_block);
        system_rhs.reinit(dofs_per_block);
        
    }
    
    void generalizeAlpha::getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij)
    {   // get element-level metric tensor
        double Jp11,Jp12,Jp21,Jp22;
        
        Jp11 = JacobiInverse[0][0];
        Jp12 = JacobiInverse[0][1];
        Jp21 = JacobiInverse[1][0];
        Jp22 = JacobiInverse[1][1];
        
        gij[0] =(Jp11)*(Jp11)+(Jp21)*(Jp21);
        gij[1] =(Jp11)*(Jp12)+(Jp21)*(Jp22);
        gij[2] =(Jp12)*(Jp11)+(Jp22)*(Jp21);
        gij[3] =(Jp12)*(Jp12)+(Jp22)*(Jp22);
        //std::cout << " j "<<gij[0] <<" "<<gij[1] <<" "<< gij[2] <<" "<<gij[3] << '\n';
    }
    
    void generalizeAlpha:: pre_compute(MappingQ<2> mapping, FESystem<2> fe,
                                       DoFHandler<2> &dof_handler,
                        int EleNum,double *** &diffusive_k,
                        double*** &div_SF_u,double**** &grad_SF_u,
                        double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,
                                       double*** &gijG)
    {
        
        //int cellnum = 0;
        QGauss<2> quadrature_formula(degree*2);
        
        FEValues<2> fe_values(mapping, fe, quadrature_formula,
                              update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_inverse_jacobians);
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points    = quadrature_formula.size();
        std::cout << "n_q_points  " << n_q_points<<std::endl;
        
        std::vector<double>         div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
        std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
        std::vector<double>         phi_p(dofs_per_cell);
        std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
        
        FullMatrix<double> lhstest(dofs_per_cell, dofs_per_cell);
        
        double* gij = new double[4];
        int cellnum = 0;
        double cc = (alphaf*gamma*dt);
        for (; cell != endc; ++cell)
        {
            /*
            int vertexNum = GeometryInfo<2>::vertices_per_cell;
            double* coordx = new double[vertexNum];
            double* coordy = new double[vertexNum];
            
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                Point<2> &v = cell->vertex(i);
                coordx[i] =  v[0];
                coordy[i] = v[1];
            }
            
            std::cout << "cellnum  " << cellnum<<std::endl;
            
            std::cout << "coordinatepre  " << coordx[0]<<" "<< coordx[1]<<" "<< coordx[2]<<" "<< coordx[3]<<std::endl;
            std::cout << "coordinatepre  " << coordy[0]<<" "<< coordy[1]<<" "<< coordy[2]<<" "<< coordy[3]<<std::endl;
            */
            lhstest = 0;
            fe_values.reinit(cell);
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                getGij(fe_values.inverse_jacobian(q),gij);
                //std::cout << "Jppre  " << fe_values.inverse_jacobian(q)[0][0]<<std::endl;
                //std::cout << "gij  " << gij[0] <<" "<< gij[1]<<" "<<gij[2]<<" "<<gij[3]<<std::endl;
                gijG[cellnum][q][0]=gij[0];
                gijG[cellnum][q][1]=gij[1];
                gijG[cellnum][q][2]=gij[2];
                gijG[cellnum][q][3]=gij[3];
                
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                    
                    div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                    grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                    phi_u[k]      =  fe_values[velocities].value(k, q);
                    phi_p[k]      =  fe_values[pressure].value(k, q);
                    grad_phi_p[k] =  fe_values[pressure].gradient(k, q);
                    
                    div_SF_u[cellnum][q][k] = fe_values[velocities].divergence (k, q);
                    grad_SF_u[cellnum][q][k][0] =grad_phi_u[k][0][0];
                    grad_SF_u[cellnum][q][k][1] =grad_phi_u[k][0][1];
                    grad_SF_u[cellnum][q][k][2] =grad_phi_u[k][1][0];
                    grad_SF_u[cellnum][q][k][3] =grad_phi_u[k][1][1];
                    
                    SF_u[cellnum][q][k][0] = phi_u[k][0];
                    SF_u[cellnum][q][k][1] = phi_u[k][1];
                    
                    SF_p[cellnum][q][k] = fe_values[pressure].value(k, q);
                    grad_SF_p[cellnum][q][k][0] = grad_phi_p[k][0];
                    grad_SF_p[cellnum][q][k][1] = grad_phi_p[k][1];
                    //div_phi_u[k];
                    //std::cout<<"grad_SF_u "<<div_SF_u[cellnum][q][k]<<std::endl;
                }
                
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                        //std::cout << "grad_phi_u  " << grad_phi_u[j] <<std::endl;
                        
                         lhstest(i,j) +=
                         alpham*phi_u[j]*phi_u[i]*fe_values.JxW(q)
                         +(viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])
                           - div_phi_u[i]*phi_p[j]
                           + phi_p[i]*div_phi_u[j]
                           )
                         *cc*fe_values.JxW(q);
                        
                        //lhstest(i,j) += viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])*cc*fe_values.JxW(q);
                        
                    }
                }
                
            }
            
            //std::cout<<"cellnum "<<cellnum<<std::endl;
            for ( unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                //local_rhs(k) = rhs[k];
                for ( unsigned int l = 0; l < dofs_per_cell; ++l)
                {
                    //local_matrix(k,l) = lhs[k][l]+diffusive_k[cellnum-1][k][l];
                    diffusive_k[cellnum][k][l] = lhstest(k,l);
                    //std::cout <<" "<< diffusive_k[cellnum][k][l]<<" ";
                    //std::cout <<" diff"<< diffusive_k[cellnum-1][k][l]<<" ";
                }
                //std::cout<<std::endl;
                
            }
            cellnum++;
        }
        
        
    }
    
    void generalizeAlpha::predictor(BlockVector<double> &n_solution,
                                    BlockVector<double> &n_solution_time_derivative,
                                    BlockVector<double> &np1_solution,
                                    BlockVector<double> &np1_solution_time_derivative)
    {
        
        // guess on ut_np1
        //printf ("sizeu: %lu  \n", n_solution_time_derivative.block(0).size());
        //printf ("sizep: %lu  \n", n_solution_time_derivative.block(1).size());
        
        Vector<double> ut_n (n_solution_time_derivative.block(0).size());
        Vector<double> pt_n (n_solution_time_derivative.block(1).size());
        Vector<double> ut_np1 (np1_solution_time_derivative.block(0).size());
        Vector<double> pt_np1 (np1_solution_time_derivative.block(1).size());
        
        Vector<double> u_n (n_solution.block(0).size());
        Vector<double> p_n (n_solution.block(1).size());
        Vector<double> u_np1 (np1_solution.block(0).size());
        Vector<double> p_np1 (np1_solution.block(1).size());
        
        
        ut_n =n_solution_time_derivative.block(0);
        pt_n =n_solution_time_derivative.block(1);
        ut_np1 =np1_solution_time_derivative.block(0);
        pt_np1 =np1_solution_time_derivative.block(1);
        u_n =n_solution.block(0);
        p_n =n_solution.block(1);
        u_np1 =np1_solution.block(0);
        p_np1 =np1_solution.block(1);
        
        // use ODE1 to guess u_np1, p_np1
        for(std::size_t i=0; i<ut_n.size(); ++i)
        {
            ut_np1[i] = ((gamma-1)/gamma)*ut_n[i];
            //std::cout<<ut_np1[i]<<std::endl;
            u_np1[i] = u_n[i]+dt*ut_n[i]+gamma*dt*(ut_np1[i]-ut_n[i]);
            //std::cout<<u_np1[i]<<std::endl;
        }
        np1_solution_time_derivative.block(0) = ut_np1;
        np1_solution.block(0) =u_np1;
        
        for(std::size_t i=0; i<pt_n.size(); ++i)
        {
            pt_np1[i] = (0)*pt_n[i];
            //std::cout<<pt_np1[i]<<std::endl;
            p_np1[i] = p_n[i]+dt*pt_n[i]+dt*(pt_np1[i]-pt_n[i]);
            //std::cout<<p_np1[i]<<std::endl;
        }
        
        np1_solution_time_derivative.block(1) = pt_np1;
        np1_solution.block(1) =p_np1;
        
    }
    
    void generalizeAlpha::corrector(MappingQ<2> mapping,
                                    int maxNewtonIter,double tolerance,
                                    DoFHandler<2> &dof_handler,FESystem<2> fe,
                                    ConstraintMatrix &zero_constraints,
                                    ConstraintMatrix &nonzero_constraints,
                                    BlockVector<double> &n_solution,
                                    BlockVector<double> &n_solution_time_derivative,
                                    BlockVector<double> &np1_solution,
                                    BlockVector<double> &np1_solution_time_derivative,
                                    double *** &diffusive_k,
                                    double*** &div_SF_u,double**** &grad_SF_u,
                                    double**** &SF_u, double*** &SF_p,
                                    double**** &grad_SF_p,double*** &gijG)
    {
        /*
         newtons loop with a very easy line search algorithm
         difference between first iteration and not first is for wich constraint to use
         */
        
        double current_res=0;
        //double last_res;
        std::cout << "viscosity  in generalizeA " << viscosity<<degree<<gamma<<dt<<alphaf<<alpham<<std::endl;
        NewtonSolve newtonSystem(degree,viscosity, gamma,dt,alphaf,alpham,
                                 dof_handler,fe,mapping,sparsity_pattern,
                                 zero_constraints,
                                 system_matrix,system_rhs);
        
        for(int i=0; i<maxNewtonIter; ++i)
        {
            //bool isfirst = FALSE;
            std::cout << "----start iteration Number ----" <<i<< std::endl;
        
            newtonSystem.iterPC(n_solution,n_solution_time_derivative,
                                np1_solution,np1_solution_time_derivative,
                                npaf_solution,npam_solution_time_derivative);
            
            
            std::cout << "----get newton update ----" << std::endl;
            newtonSystem.getNewtonUpdate(npaf_solution,npam_solution_time_derivative,newton_update,diffusive_k,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG);
            
            std::cout << "-------now update flow-------" << std::endl;
            newtonSystem.update_flow(newton_update,np1_solution,np1_solution_time_derivative);
            /*
            Vector<double> du_n = np1_solution.block(0);
            Vector<double> dp_n = np1_solution.block(1);
            std::cout << "np1_solution " <<std::endl;
            for(std::size_t i=0; i<du_n.size(); ++i)
            {
                std::cout<<du_n[i]<<std::endl;
            }
            
            for(std::size_t i=0; i<dp_n.size(); ++i)
            {
                std::cout<<dp_n[i]<<std::endl;
            }
            */
            std::cout << "-------now get residual-------" << std::endl;
            newtonSystem.getResidual(current_res);
                
            //std::cout << "residue is "<<current_res << std::endl;
            
            /*
             Vector<double> du_n = np1_solution.block(0);
             Vector<double> dp_n = np1_solution.block(1);
             std::cout << "updated_flow " <<std::endl;
             for(std::size_t i=0; i<du_n.size(); ++i)
             {
             std::cout<<du_n[i]<<std::endl;
             }
             
            for(std::size_t i=0; i<dp_n.size(); ++i)
            {
                std::cout<<dp_n[i]<<std::endl;
            }
            */
            if (current_res < tolerance)
                break;
        }
        
        
    }
    
    void generalizeAlpha::updator(BlockVector<double> &n_solution,
                                  BlockVector<double> &n_solution_time_derivative,
                                  BlockVector<double> &np1_solution,
                                  BlockVector<double> &np1_solution_time_derivative)
    {
        
        n_solution_time_derivative.block(0) = np1_solution_time_derivative.block(0);
        n_solution.block(0) =np1_solution.block(0);
        
        n_solution_time_derivative.block(1) = np1_solution_time_derivative.block(1);
        n_solution.block(1) =np1_solution.block(1);
        
        /* update flow 
        u_n = u_np1;
        ut_n = ut_np1;
        p_n = p_np1;
        pt_n = pt_np1;
         */
        
    }
    
}


