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

#include "../include/newton.h"
#include "../include/pre.h"
#include "../include/user_input.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
namespace incompressible
{
    using namespace dealii;
    NewtonSolve::NewtonSolve(double viscosity,
                             hp::DoFHandler<2> &dof_handler,
                             hp::FECollection<2> &fe,
                             BlockSparsityPattern &sparsity_pattern,
                             AffineConstraints<double> &zero_constraints,
                             BlockVector<double> &n_solution,
                             BlockSparseMatrix<double> &system_matrix,
                             BlockVector<double> &system_rhs)
    :
    viscosity(viscosity),
    dof_handler(dof_handler),
    fe(fe),
    sparsity_pattern(sparsity_pattern),
    zero_constraints(zero_constraints),
    n_solution(n_solution),
    system_matrix(system_matrix),
    system_rhs(system_rhs)
    {}
    
    
    
    void NewtonSolve::setup_system()
    {
        for (unsigned int deg=min_degree; deg<=max_degree; ++deg)
        {
            quadrature_collection.push_back (QGauss<2>(deg*2));
        }
        
        hp::FEValues<2> hp_fe_values(fe,
                                quadrature_collection,
                                update_values | update_quadrature_points |
                                update_JxW_values | update_gradients |
                                update_inverse_jacobians);
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        FullMatrix<double> local_matrix;
        Vector<double>     local_rhs;
        Vector<double> rhstest;
        FullMatrix<double> lhstest;
        
        std::vector<types::global_dof_index> local_dof_indices;

        hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        hp::DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

        int cellnum = 0;
        
        double* gij = new double[4];
        double rc,tauM,tauC,tauBar;
        std::cout << "start cell iteration " <<std::endl;
        
        system_matrix = 0;
        system_rhs = 0;
        
        for (; cell != endc; ++cell)
        {
            cellnum++;
            //std::cout<<"cellnum "<<cellnum<<std::endl;
            const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
            local_matrix.reinit (dofs_per_cell, dofs_per_cell);
            local_matrix = 0;
            
            local_rhs.reinit (dofs_per_cell);
            local_rhs = 0;
            hp_fe_values.reinit(cell);
            
            std::vector<double>         div_phi_u(dofs_per_cell);
            std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
            std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
            std::vector<double>         phi_p(dofs_per_cell);
            std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
            
            const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
            unsigned int n_q_points  = fe_values.n_quadrature_points;
            rhstest.reinit (dofs_per_cell);
            rhstest = 0;
            lhstest.reinit (dofs_per_cell, dofs_per_cell);
            lhstest = 0;
            
            std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
            std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
            std::vector<double>       present_pressure_values(n_q_points);
            std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
            
            fe_values[velocities].get_function_values(n_solution,present_velocity_values);
            fe_values[velocities].get_function_gradients(n_solution, present_velocity_gradients);
            fe_values[pressure].get_function_values(n_solution,
                present_pressure_values);
            fe_values[pressure].get_function_gradients(n_solution,
                present_pressure_gradients);
            
            std::vector<Tensor<1, 2> > rm(n_q_points);
            std::vector<Tensor<1, 2> > velocitybar(n_q_points);
            
            std::vector<Tensor<1, 2> >        lhs1(dofs_per_cell);
            std::vector<Tensor<1, 2> >        lhs2(dofs_per_cell);
            std::vector<Tensor<1, 2> >        lhs3(dofs_per_cell);
            std::vector<double >              lhs4(dofs_per_cell);
            std::vector<Tensor<1, 2> >        lhs5(dofs_per_cell);
            std::vector<Tensor<1, 2> >        rhs1(n_q_points);
            
            int celldeg = cell->active_fe_index()+1;
            //std::cout<<"celldeg "<<celldeg<<std::endl;
            //celldeg = 1;
            // loop through intergral points
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                
                const DerivativeForm<1, 2, 2> &JacobiInverse = fe_values.inverse_jacobian(q);
                getGij(JacobiInverse,gij);
                
                rm[q] = present_velocity_gradients[q]*present_velocity_values[q]+present_pressure_gradients[q];//-sourceTerm[q];
                
                rc = present_velocity_gradients[q][0][0]+present_velocity_gradients[q][1][1];
                
                getTaum(celldeg,gij,present_velocity_values[q][0],present_velocity_values[q][1],rm[q][0],rm[q][1],tauM,tauC,tauBar);
                //tauC = (1)/(3*(degree*degree)*tauM*(gij[0]+gij[3]));
                
                velocitybar[q] = present_velocity_values[q]-tauM*rm[q];
                double Jx = fe_values.JxW(q);
                double present_velocity_divergence =trace(present_velocity_gradients[q]);
                
                /*
                if (cellnum == 3819){
                    std::cout<<"gij "<<gij[0]<<std::endl;
                    std::cout<<"tauM "<<tauM<<std::endl;
                    std::cout<<"tauBar "<<tauBar<<std::endl;
                }
                 */
                rhs1[q]=-present_velocity_gradients[q]*velocitybar[q];
                
                //std::cout<<"tauM "<<tauM<<std::endl;
                //std::cout<<"rhs1[q] "<<rhs1[q]<<std::endl;
                // shape functions
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                    div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                    grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                    phi_u[k]      =  fe_values[velocities].value(k, q);
                    phi_p[k]      =  fe_values[pressure].value(k, q);
                    grad_phi_p[k] =  fe_values[pressure].gradient(k, q);
                    
                    //std::cout<<"grad_phi_p "<<grad_phi_p[k]<<std::endl;
                    
                    lhs1[k] = grad_phi_u[k]*velocitybar[q]+
                            present_velocity_gradients[q]*phi_u[k];
                    lhs2[k] = grad_phi_u[k]*present_velocity_values[q];
                    lhs3[k] = tauM*lhs2[k];
                    lhs4[k] = tauC*div_phi_u[k];
                    lhs5[k] = tauM*grad_phi_p[k];
                }
                // i ,j loop
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    
                    local_rhs(i) += (rhs1[q]*phi_u[i]
                        - viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                        + (present_pressure_values[q]- rc*tauC)*div_phi_u[i]
                        - present_velocity_divergence*phi_p[i]
                        - tauM*grad_phi_u[i]*present_velocity_values[q]*rm[q]
                        - tauM*rm[q]*grad_phi_p[i])* Jx
                        -grad_phi_u[i]*rm[q]*tauBar*
                        present_velocity_gradients[q]*rm[q]*Jx;
                     
                    //rhstest(i) +=-grad_phi_u[i]*rm[q]*tauBar*
                     //present_velocity_gradients[q]*rm[q]*Jx;
                    
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            local_matrix(i,j)+= (lhs1[j]*phi_u[i]
                                        +lhs3[j]*lhs2[i]
                                        +lhs4[j]*div_phi_u[i]
                                        +lhs5[j]*grad_phi_p[i]
                                        + viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])
                                        - div_phi_u[i]*phi_p[j]
                                        + phi_p[i]*div_phi_u[j])*Jx
                                        +grad_phi_u[i]*rm[q]*tauBar
                                        *grad_phi_u[j]*rm[q]*Jx;
                        //lhstest(i,j) += grad_phi_u[i]*rm[q]*tauBar
                          //  *grad_phi_u[j]*rm[q]*Jx;
                        //lhs3[k] = tauM*lhs2[k];* fe_values.JxW(q);
                        
                        }
                }
            }
            /*
            if (cellnum == 3819){
            for ( unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                std::cout << "rhstest  " << rhstest(k) <<std::endl;
                
            }
            
            //std::cout << "local_matrix  ";
            
            
            for ( unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                for ( unsigned int l = 0; l < dofs_per_cell; ++l)
                {
                    std::cout <<" "<< lhstest(k,l)<<" ";
                }
                std::cout<<std::endl;
                
            }
            
            }
            */
            //std::cout << "cellnum bc " << cellnum <<std::endl;
         
            local_dof_indices.resize (dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);
            
            //nonzero_constraints.clear();
            const AffineConstraints<double> &constraints_used = zero_constraints;
            
            constraints_used.distribute_local_to_global(local_matrix,
                                                        local_rhs,
                                                        local_dof_indices,
                                                        system_matrix,
                                                        system_rhs);
            
            
            
        }
        /*
        std::cout << "rhs=[ " <<std::endl;
        for(std::size_t i=0; i<system_rhs.size(); ++i)
        {
            
            std::cout<<system_rhs[i]<<std::endl;
        }
         */
        //add_flux();
       
    }
    
    void NewtonSolve::getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij)
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
    
    void NewtonSolve::getTaum(int degree, double* gij,double uele,double vele, double rum, double rvm,double &tauM, double &tauC, double &tauBar)
    {

        double C2 = (3*degree*degree)*(3*degree*degree);
        double gij0 = 4*gij[0];
        double gij1 = 4*gij[1];
        double gij2 = 4*gij[2];
        double gij3 = 4*gij[3];
        double tau1sqinv = uele*(uele*gij0+vele*gij2)+vele*(uele*gij1+vele*gij3);
        double A11 =gij0*gij0+gij2*gij2;
        double A22 =gij1*gij1+gij3*gij3;
        double tau2sqinv = C2*viscosity*viscosity*(A11+A22);
        tauM = 1/sqrt(tau1sqinv+tau2sqinv);
        tauC = (1)/((degree*degree)*tauM*(gij0+gij3));
        double tauRsqinv = rum*(rum*gij0+rvm*gij2)+rvm*(rum*gij1+rvm*gij3);
        if (tauRsqinv >0)
            tauBar = tauM/sqrt(tauRsqinv);
        else
            tauBar = tauRsqinv;
        
        //std::cout << "tau1sqinv  " << tau1sqinv<<std::endl;
        //std::cout << "A11  " << A11<<std::endl;
        //std::cout << "A22  " << A22<<std::endl;
    }
    
    void NewtonSolve::add_flux()
    {
        
        bool ifconsist_diff; double value_diffx; double value_diffy;
        bool ifconsist_p; double value_px; double value_py;
        
        for (unsigned int deg=min_degree; deg<=max_degree; ++deg)
            face_quadrature_collection.push_back (QGauss<1>(deg*2));
        
        hp::FEFaceValues<2> hp_fe_face_values(fe,
                                     face_quadrature_collection,
                                     update_values | update_quadrature_points |
                                          update_normal_vectors |
                                     update_JxW_values | update_gradients );
        
        
        hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        hp::DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        std::vector<types::global_dof_index> local_dof_indices;
        
        
        for (; cell != endc; ++cell)
        {
            const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();
            unsigned int n_face_q_points  = fe_face_values.n_quadrature_points;
            
            std::vector<Tensor<1, 2>> face_velocity_values(n_face_q_points);
            std::vector<Tensor<2, 2>> face_velocity_gradients(n_face_q_points);
            std::vector<double>       face_pressure_values(n_face_q_points);
            std::vector<Tensor<1, 2>> face_pressure_gradients(n_face_q_points);
            
            const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
            FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
            Vector<double>     local_rhs(dofs_per_cell);

            std::vector<double>         div_phi_u(dofs_per_cell);
            std::vector<Tensor<1, 2>> phi_u(dofs_per_cell);
            std::vector<Tensor<2, 2>> grad_phi_u(dofs_per_cell);
            std::vector<double>         phi_p(dofs_per_cell);
            std::vector<Tensor<1, 2>> grad_phi_p(dofs_per_cell);
            
            
            std::vector<Tensor<1, 2>>        diff_flux(dofs_per_cell);
            std::vector<Tensor<1, 2>>        pressure_flux(dofs_per_cell);
        
            
            for (unsigned int face_number = 0;
                 face_number < GeometryInfo<2>::faces_per_cell;
                 ++face_number)
            {
                //add diffusive flux
                
                if (cell->face(face_number)->at_boundary())
                {
                    hp_fe_face_values.reinit(cell, face_number);
                    const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();
                    
                    fe_face_values[velocities].get_function_values(n_solution,face_velocity_values);
                    fe_face_values[velocities].get_function_gradients(n_solution, face_velocity_gradients);
                    fe_face_values[pressure].get_function_values(n_solution,
                                                                 face_pressure_values);
                    fe_face_values[pressure].get_function_gradients(n_solution,
                                                    face_pressure_gradients);
                    
                    Point<2> &v1 = cell->face(face_number)->vertex(0);
                    Point<2> &v2 = cell->face(face_number)->vertex(1);
                    
                    //std::cout << "coordinate  " << v1[0]<<" "<< v1[1]<<" "<<std::endl;
                    //std::cout << "coordinate  " << v2[0]<<" "<< v2[1]<<" "<<std::endl;
                    
                    getflux(v1,v2,'d',ifconsist_diff,value_diffx,value_diffy);
                    getflux(v1,v2,'p',ifconsist_p,value_px,value_py);
                    
                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         ++q_point)
                    {
                        //std::cout << "face_velocity_values  " << face_velocity_values[q_point]<<std::endl;
                        //std::cout << "face_velocity_gradients  " << face_velocity_gradients[q_point]<<std::endl;
                        Tensor<1, 2> normal_vector =fe_face_values.normal_vector(q_point);
                        
                        for (unsigned int k=0; k<dofs_per_cell; ++k)
                        {
                            div_phi_u[k]  =  fe_face_values[velocities].divergence (k, q_point);
                            grad_phi_u[k] =  fe_face_values[velocities].gradient(k, q_point);
                            phi_u[k]      =  fe_face_values[velocities].value(k, q_point);
                            phi_p[k]      =  fe_face_values[pressure].value(k, q_point);
                            grad_phi_p[k] =  fe_face_values[pressure].gradient(k, q_point);
                            
                            //std::cout<<"grad_p "<<grad_phi_p[k]<<std::endl;
                        }
                        
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            
                            if (ifconsist_diff == TRUE)
                            { diff_flux[i][0]=viscosity*(face_velocity_gradients[q_point][0][0]+face_velocity_gradients[q_point][0][0])*normal_vector[0]+viscosity*(face_velocity_gradients[q_point][0][1]+face_velocity_gradients[q_point][1][0])*normal_vector[1];
                                
                                diff_flux[i][1]=viscosity*(face_velocity_gradients[q_point][1][0]+face_velocity_gradients[q_point][0][1])*normal_vector[0]+viscosity*(face_velocity_gradients[q_point][1][1]+face_velocity_gradients[q_point][1][1])*normal_vector[1];
                                
                            }
                            else
                            {
                                diff_flux[i][0] = value_diffx;
                                diff_flux[i][1] = value_diffy;
                            }
                            
                            if (ifconsist_p == TRUE)
                            {  pressure_flux[i][0]=face_pressure_values[q_point]*normal_vector[0];
                                pressure_flux[i][1]=face_pressure_values[q_point]*normal_vector[1];
                                //std::cout << "face_pressure_values  " << face_pressure_values[q_point] <<std::endl;
                            }
                            else
                            {
                                pressure_flux[i][0] = value_px;
                                pressure_flux[i][1] = value_py;
                            }
                            
                        }
                        //getdiffFlux(ifconsist,value);
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            
                            local_rhs(i) +=
                            (diff_flux[i] *                          // g(x_q)
                             phi_u[i] * // phi_i(x_q)
                             fe_face_values.JxW(q_point));            // dx
                            
                            local_rhs(i) -=
                            (pressure_flux[i] *                          // g(x_q)
                             phi_u[i] * // phi_i(x_q)
                             fe_face_values.JxW(q_point));            // dx
                            //std::cout << "rhstest  " << rhstest(i) <<std::endl;
                            
                        }
                    }
                    
                }
                
                
            }
            
            local_dof_indices.resize (dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);
            const AffineConstraints<double> &constraints_used = zero_constraints;
            constraints_used.distribute_local_to_global(local_rhs,
                                                        local_dof_indices,
                                                        system_rhs);
            
        }
        
    }
    void NewtonSolve::linear_solve(BlockVector<double> &newton_update)
    {
        //call linear solver
        
        //system_matrix,system_rhs,newton_update
        
        const AffineConstraints<double> &constraints_used = zero_constraints;
        
        /*
        SolverControl  solver_control(system_matrix.m(),1e-4 * system_rhs.l2_norm(),true);
        SolverFGMRES<BlockVector<double>> gmres(solver_control);
        
        SparseILU<double> pmass_preconditioner;
        pmass_preconditioner.initialize(pressure_mass_matrix,SparseILU<double>::AdditionalData());
        
        const BlockSchurPreconditioner<SparseILU<double>> preconditioner(gamma,viscosity,system_matrix,pressure_mass_matrix,pmass_preconditioner);
        
        gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
        */

        
        newton_update=system_rhs;
        SparseDirectUMFPACK directsolver;
        directsolver.initialize(system_matrix);
        directsolver.solve(newton_update);
        
        //std::cout << " ****FGMRES steps: " << solver_control.last_step()
        //<< std::endl;
        
        constraints_used.distribute(newton_update);
        
    }
    
    void NewtonSolve::update_flow(BlockVector<double> &newton_update,BlockVector<double> &n_solution)
    {
        /*
        // update ut, pt
        u_n = u_n+du;
        p_n = p_n+dp;
        */
                
        Vector<double> u_n (n_solution.block(0).size());
        Vector<double> p_n (n_solution.block(1).size());
        Vector<double> du (newton_update.block(0).size());
        Vector<double> dp (newton_update.block(1).size());
        
        u_n = n_solution.block(0);
        du = newton_update.block(0);
        
        p_n = n_solution.block(1);
        dp = newton_update.block(1);
        
        for(std::size_t i=0; i<u_n.size(); ++i)
        {
            u_n[i] = u_n[i]+du[i];
        }
        
        for(std::size_t i=0; i<p_n.size(); ++i)
        {
            p_n[i] = p_n[i]+dp[i];
        }
        
        n_solution.block(0) = u_n;
        n_solution.block(1) = p_n;
    }

    void NewtonSolve::assemble_rhs()
    {
        setup_system();
    }
    
    
    void NewtonSolve::getNewtonUpdate( BlockVector<double> &n_solution, BlockVector<double> &newton_update)
    {
        
        setup_system();
        linear_solve(newton_update);
    }
    
    void NewtonSolve::getResidual(BlockVector<double> &n_solution, double &current_res)
    {
        //assemble_rhs();
        
        current_res = system_rhs.l2_norm();
        
        std::cout << "******************************" << std::endl;
        std::cout << " The residual of this guess is " << current_res << std::endl;
        
    }
    
    void NewtonSolve::check_convergence( double &res)
    {
     //assemble_rhs();
     
     res = system_rhs.l2_norm();
     
     std::cout << "******************************" << std::endl;
     std::cout << " The residual of this guess is " << res << std::endl;
     std::cout << " Initialization complete!  " << std::endl;
     
    }
    
}


