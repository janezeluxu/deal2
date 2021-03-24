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
#include "../include/newton.h"
#include "../include/linear.h"
#include "../include/pre.h"
namespace incompressible
{
    using namespace dealii;
    NewtonSolve::NewtonSolve(int degree,double viscosity, double gamma, double dt,
                             double alphaf, double alpham,
                             DoFHandler<2> &dof_handler,FESystem<2> &fe,
                             MappingQ<2> &mapping,
                             BlockSparsityPattern &sparsity_pattern,
                             ConstraintMatrix &zero_constraints,
                             //SparseMatrix<double> &pressure_mass_matrix, //temp
                             BlockSparseMatrix<double> &system_matrix,
                             BlockVector<double> &system_rhs)
    :
    degree(degree),
    viscosity(viscosity),
    gamma(gamma),
    dt(dt),
    alphaf(alphaf),
    alpham(alpham),
    dof_handler(dof_handler),
    fe(fe),
    mapping(mapping),
    sparsity_pattern(sparsity_pattern),
    zero_constraints(zero_constraints),
    //pressure_mass_matrix(pressure_mass_matrix),
    system_matrix(system_matrix),
    system_rhs(system_rhs)
    {}
    
    
    void NewtonSolve::iterPC(BlockVector<double> &n_solution,BlockVector<double> &n_solution_time_derivative, BlockVector<double> &np1_solution,BlockVector<double> &np1_solution_time_derivative,BlockVector<double> &npaf_solution,BlockVector<double> &npam_solution_time_derivative)
    {
        //printf ("sizeuiterPC: %lu  \n", n_solution.block(0).size());
        //printf ("sizepiterPC: %lu  \n", n_solution.block(1).size());
        
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
        
        
        Vector<double> u_npaf (n_solution.block(0).size());
        Vector<double> ut_npam (n_solution_time_derivative.block(0).size());
        
        //printf ("sizeuiterPCu_npaf: %lu  \n", u_npaf.size());
        //printf ("sizepiterPCut_npam: %lu  \n", ut_npam.size());
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
        for(std::size_t i=0; i<ut_n.size(); ++i)
        {
            u_npaf[i] = u_n[i]+alphaf*(u_np1[i]-u_n[i]);
            //u_npaf[i] = 1;
            //std::cout<<u_npaf[i]<<std::endl;
            ut_npam[i] = ut_n[i]+alpham*(ut_np1[i]-ut_n[i]);
            //std::cout<<ut_npam[i]<<std::endl;
        }
        
        npaf_solution.block(0) = u_npaf;
        npaf_solution.block(1) = p_np1;
        
        npam_solution_time_derivative.block(0) = ut_npam;
        npam_solution_time_derivative.block(1) = pt_np1;
    }
    
    void NewtonSolve::setup_system(BlockVector<double> &npaf_solution,
                                   BlockVector<double> &npam_solution_time_derivative,
                                   double *** diffusive_k,
                                   double*** &div_SF_u,double**** &grad_SF_u,
                                   double**** &SF_u, double*** &SF_p,
                                   double**** &grad_SF_p,double*** &gijG)
    {
        
        //int dim = 2;
        system_matrix = 0;
        system_rhs = 0;
        QGauss<2> quadrature_formula(degree*2);
        
        
        FEValues<2> fe_values(mapping, fe,quadrature_formula,
                                update_values | update_quadrature_points |
                                update_JxW_values | update_gradients |
                                update_inverse_jacobians);
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points    = quadrature_formula.size();
        
        
        FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     local_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        
        std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
        std::vector<Tensor<1, 2> > present_vt_values(n_q_points);
        
        std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
        std::vector<double>       present_pressure_values(n_q_points);
        std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
        
        std::vector<double>         div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
        std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
        std::vector<double>         phi_p(dofs_per_cell);
        std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
        
        std::vector<double>         SF(dofs_per_cell);
        std::vector<Tensor<1, 2> > SFgrad(dofs_per_cell);
        std::vector<Tensor<1, 2> > rm(n_q_points);
        //std::vector<double>       rc(n_q_points);
        std::vector<Tensor<1, 2> > sourceTerm(n_q_points);
        std::vector<Tensor<1, 2> > velocitybar(n_q_points);
        
        std::vector<Tensor<1, 2> >        diff_flux(dofs_per_cell);
        std::vector<Tensor<1, 2> >        pressure_flux(dofs_per_cell);
        //double viscosity = 0.01;
        //int totalSF = dofs_per_cell/3;
        //Vector<double>     Rum(totalSF);
        //Vector<double>     Rvm(totalSF);
        //Vector<double>     Rc(totalSF);
        Vector<double>     rhstest(dofs_per_cell);
        //std::vector<Tensor<1, 2>>   Rm(dofs_per_cell);
        FullMatrix<double> lhstest(dofs_per_cell, dofs_per_cell);
        std::vector<Tensor<1, 2> >        lhs1(dofs_per_cell);
        std::vector<Tensor<1, 2> >        lhs2(dofs_per_cell);
        std::vector<Tensor<1, 2> >        lhs3(dofs_per_cell);
        std::vector<double >        lhs4(dofs_per_cell);
        std::vector<Tensor<1, 2> >        lhs5(dofs_per_cell);
        std::vector<Tensor<1, 2> >        rhs1(n_q_points);
        
        
        /*
        Vector<double> du_n = npaf_solution.block(0);
        Vector<double> dp_n = npaf_solution.block(1);
        std::cout << "npaf_solution " <<std::endl;
        for(std::size_t i=0; i<du_n.size(); ++i)
        {
            std::cout<<du_n[i]<<std::endl;
        }
        
        for(std::size_t i=0; i<dp_n.size(); ++i)
        {
            std::cout<<dp_n[i]<<std::endl;
        }
        */
        
        //typename DoFHandler<2>::active_cell_iterator cell =
        //dof_handler.begin_active(),
        //endc = dof_handler.end();
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        int cellnum = 0;
        
        double* gij = new double[4];
        double rc,tauM,tauC;
        std::cout << "start cell iteration " <<std::endl;
        //std::cout << "start cell iteration " <<std::endl;
        double cc = (alphaf*gamma*dt);
        
        
        double** lhs = new double*[dofs_per_cell];
        for(unsigned int i=0; i<dofs_per_cell; i++)
        {
            lhs[i] = new double[dofs_per_cell];
        }
        //double* rhs = new double[dofs_per_cell];
        
        for (; cell != endc; ++cell)
        {
            /*
            //std::cout << "cc " <<cc<<std::endl;
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
            
            std::cout << "coordinate  " << coordx[0]<<" "<< coordx[1]<<" "<< coordx[2]<<" "<< coordx[3]<<std::endl;
            std::cout << "coordinate  " << coordy[0]<<" "<< coordy[1]<<" "<< coordy[2]<<" "<< coordy[3]<<std::endl;
            */
            //std::cout << "cellnum  " << cellnum<<std::endl;
            cellnum++;
            
            //ele();
            fe_values.reinit(cell);
            local_matrix = 0;
            local_rhs    = 0;
            rhstest = 0;
            lhstest = 0;
            
            for (unsigned int i = 0; i<dofs_per_cell; i++)
            {
                //rhs[i] = 0.0;
                for (unsigned int j = 0; j<dofs_per_cell; j++)
                {
                    lhs[i][j]=0.0;
                }
                
            }
            
            fe_values[velocities].get_function_values(npaf_solution,present_velocity_values);
            fe_values[velocities].get_function_gradients(npaf_solution, present_velocity_gradients);
            fe_values[pressure].get_function_values(npaf_solution,
                present_pressure_values);
            fe_values[pressure].get_function_gradients(npaf_solution,
                present_pressure_gradients);
            
            fe_values[velocities].get_function_values(npam_solution_time_derivative,present_vt_values);
            /*
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                sourceTerm[q][0] = 0;
                sourceTerm[q][1] = 0;
            }
            */
            //std::cout << "n_q_points  " << n_q_points<<std::endl;
            
            // loop through intergral points
            for (unsigned int q=0; q<n_q_points; ++q)
            {

                
                //getGij(fe_values.inverse_jacobian(q),gij);
                //std::cout << "gij  " << gij[0] <<" "<< gij[1]<<" "<<gij[2]<<" "<<gij[3]<<std::endl;
                
                gij[0] = gijG[cellnum-1][q][0];
                gij[1] = gijG[cellnum-1][q][1];
                gij[2] = gijG[cellnum-1][q][2];
                gij[3] = gijG[cellnum-1][q][3];
                //std::cout << "gijpre  " << gij[0] <<" "<< gij[1]<<" "<<gij[2]<<" "<<gij[3]<<std::endl;
                
                rm[q] = present_vt_values[q] +
                present_velocity_gradients[q]*present_velocity_values[q]+present_pressure_gradients[q];//-sourceTerm[q];
                
                rc = present_velocity_gradients[q][0][0]+present_velocity_gradients[q][1][1];
                
                 getTaum(gij,dt,present_velocity_values[q][0],present_velocity_values[q][1],tauM,tauC);
                
                velocitybar[q] = present_velocity_values[q]-tauM*rm[q];

				double cJx = cc*fe_values.JxW(q);
				double Jx = fe_values.JxW(q);
				
				rhs1[q]=-present_vt_values[q]
                     -present_velocity_gradients[q]*velocitybar[q];
                     
                double present_velocity_divergence =  trace(present_velocity_gradients[q]);
                // shape functions
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                    
                    div_phi_u[k]  = div_SF_u[cellnum-1][q][k];
                    grad_phi_u[k][0][0] = grad_SF_u[cellnum-1][q][k][0];
                    grad_phi_u[k][0][1] = grad_SF_u[cellnum-1][q][k][1];
                    grad_phi_u[k][1][0] = grad_SF_u[cellnum-1][q][k][2];
                    grad_phi_u[k][1][1] = grad_SF_u[cellnum-1][q][k][3];
                    
                    phi_u[k][0] = SF_u[cellnum-1][q][k][0];
                    phi_u[k][1] = SF_u[cellnum-1][q][k][1];
                    phi_p[k]  = SF_p[cellnum-1][q][k];
                    grad_phi_p[k][0] =grad_SF_p[cellnum-1][q][k][0];
                    grad_phi_p[k][1] =grad_SF_p[cellnum-1][q][k][1];
                    
                    
                    lhs1[k] = grad_phi_u[k]*velocitybar[q]+
								present_velocity_gradients[q]*phi_u[k];
					lhs2[k] = grad_phi_u[k]*present_velocity_values[q];
					lhs3[k] = tauM*lhs2[k];
					lhs4[k] = tauC*div_phi_u[k];
					lhs5[k] = tauM*grad_phi_p[k];
                    
                    /*
                    std::cout << "phi_u_pre  " << phi_u[k] <<std::endl;
                    div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                    grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                    phi_u[k]      =  fe_values[velocities].value(k, q);
                    phi_p[k]      =  fe_values[pressure].value(k, q);
                    grad_phi_p[k] =  fe_values[pressure].gradient(k, q);
                    std::cout << "phi_u  " << phi_u[k] <<std::endl;
                    */
                }
                // i ,j loop
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    //double present_velocity_divergence =  trace(present_velocity_gradients[q]);
                    
                    //local_rhs(i)
                    //rhs[i]
                    /*
                    local_rhs(i)+=
                    (-present_vt_values[q]*phi_u[i]
                     -present_velocity_gradients[q]*velocitybar[q]*phi_u[i]
                        - viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                        + present_pressure_values[q]*div_phi_u[i]
                        - present_velocity_divergence*phi_p[i]
                        - tauM*grad_phi_u[i]*present_velocity_values[q]*rm[q]
                        - rc*tauC*div_phi_u[i]
                        - tauM*rm[q]*grad_phi_p[i])
                        * fe_values.JxW(q);
                    */  
                    local_rhs(i)+=
						(rhs1[q]*phi_u[i]
                        - viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                        + (present_pressure_values[q]- rc*tauC)*div_phi_u[i]
                        - present_velocity_divergence*phi_p[i]
                        -( lhs3[i]+lhs5[i])*rm[q]
                        )
                        * Jx;
                    
                    //rhstest(i) += (- present_velocity_divergence*phi_p[i]-tauM*rm[q]*grad_phi_p[i])*fe_values.JxW(q);
                    
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
							/*
                            //lhs[i][j]
                            //local_matrix(i,j)
                            lhs[i][j]+=
                            //alpham*phi_u[i]*phi_u[j]*fe_values.JxW(q)
                            + (grad_phi_u[j]*velocitybar[q]*phi_u[i]
                            + present_velocity_gradients[q]*phi_u[j]*phi_u[i]
                            //+ viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])
                            //- div_phi_u[i]*phi_p[j]
                            //+ phi_p[i]*div_phi_u[j]
                            + tauM*(grad_phi_u[j]*present_velocity_values[q])*
                            (grad_phi_u[i]*present_velocity_values[q])
                            + tauC*div_phi_u[j]*div_phi_u[i]
                            + tauM*grad_phi_p[i]*grad_phi_p[j]
                            )
                            * cc*fe_values.JxW(q);
                            */
                            
                            
                            lhs[i][j]+= (lhs1[j]*phi_u[i]
												+lhs3[j]*lhs2[i]
												+lhs4[j]*div_phi_u[i]
												+lhs5[j]*grad_phi_p[i])
												*cJx;
							
                            //lhstest(i,j) += scalar_product(grad_phi_u[j], grad_phi_u[i]);
                            //std::cout << "lhs  " << lhs[i][j] <<std::endl;
                        }
                }
            }
            
            
            /*
             for ( unsigned int k = 0; k < dofs_per_cell; ++k)
             {
             std::cout << "local_rhs  " << local_rhs(k) <<std::endl;
             
             }
             
             std::cout << "local_matrix  ";
             */
            
            //std::cout<<"cellnum "<<cellnum<<std::endl;
             for ( unsigned int k = 0; k < dofs_per_cell; ++k)
             {
                 //local_rhs(k) = rhs[k];
                 for ( unsigned int l = 0; l < dofs_per_cell; ++l)
                 {
                     local_matrix(k,l) = lhs[k][l]+diffusive_k[cellnum-1][k][l];
                     
                     //std::cout <<" "<< lhstest(k,l)<<" ";
                     //std::cout <<" "<< local_matrix(k,l)<<" ";
                 }
                 //std::cout<<std::endl;
             
             }
            
            cell->get_dof_indices(local_dof_indices);
            
            //nonzero_constraints.clear();
            const ConstraintMatrix &constraints_used = zero_constraints;
            
            constraints_used.distribute_local_to_global(local_matrix,
                                                        local_rhs,
                                                        local_dof_indices,
                                                        system_matrix,
                                                        system_rhs);
        }
        
            //std::cout << "cellnum bc " << cellnum <<std::endl;
            
            /* Neumann boundary condition
            First we have to find out whether the intersection of the
            faces of this cell with the newmann boundary part is nonzero.
             */
        
        //add_flux(npaf_solution,npam_solution_time_derivative);
        
        
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
    
    void NewtonSolve::getTaum(double* gij,double dt, double uele,double vele,
                                double &tauM, double &tauC)
    {
        //double uele
        //double vele
        double C1 = 1;
        double C2 = 9;//(3*degree*degree)*(3*degree*degree);
        
        double gij0 = 4*gij[0]*degree*degree;
        double gij1 = 4*gij[1]*degree*degree;
        double gij2 = 4*gij[2]*degree*degree;
        double gij3 = 4*gij[3]*degree*degree;
        
        double tau1sqinv = uele*(uele*gij0+vele*gij2)+vele*(uele*gij1+vele*gij3);
        double A11 =gij0*gij0+gij2*gij2;
        double A22 =gij1*gij1+gij3*gij3;
        
        //double tau1sqinv = uele*(uele*gij[0]+vele*gij[2])+vele*(uele*gij[1]+vele*gij[3]);
        //double A11 =gij[0]*gij[0]+gij[2]*gij[2];
        //double A22 =gij[1]*gij[1]+gij[3]*gij[3];
        
        double tau2sqinv = C2*viscosity*viscosity*(A11+A22);
        double taut = (2*C1/dt)*(2*C1/dt);
        tauM = 1/sqrt(taut+tau1sqinv+tau2sqinv);
        double tauMtilta = 1/sqrt(tau1sqinv+tau2sqinv);
        tauC = (1)/((3)*tauMtilta*(gij0+gij3));
        //double tauM =0;
        //std::cout << "degree  " << degree<<std::endl;
        //std::cout << "A11  " << A11<<" A22 "<<A22<<std::endl;
        //std::cout << "viscosity  " << viscosity<< "C2  " << C2<<std::endl;
        //std::cout << "tau2sqinv  " << tau2sqinv<<std::endl;
    }
    
    
    void NewtonSolve::add_flux(BlockVector<double> &npaf_solution)
    {
        bool ifconsist_diff; double value_diffx; double value_diffy;
        bool ifconsist_p; double value_px; double value_py;
        QGauss<1> face_quadrature_formula(degree*2);
        FEFaceValues<2> fe_face_values (fe, face_quadrature_formula,
                                        update_values | update_quadrature_points  |
                                        update_normal_vectors | update_JxW_values |update_gradients );
        const unsigned int n_face_q_points = face_quadrature_formula.size();
        std::vector<Tensor<1, 2> > face_velocity_values(n_face_q_points);
        std::vector<Tensor<2, 2> > face_velocity_gradients(n_face_q_points);
        std::vector<double>       face_pressure_values(n_face_q_points);
        std::vector<Tensor<1, 2> > face_pressure_gradients(n_face_q_points);
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     local_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        
        std::vector<double>         div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
        std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
        std::vector<double>         phi_p(dofs_per_cell);
        std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
    
        
        std::vector<Tensor<1, 2> >        diff_flux(dofs_per_cell);
        std::vector<Tensor<1, 2> >        pressure_flux(dofs_per_cell);

        
        
        for (; cell != endc; ++cell)
        {
            
            for (unsigned int face_number = 0;
                 face_number < GeometryInfo<2>::faces_per_cell;
                 ++face_number)
            {
                //add diffusive flux
                
                if (cell->face(face_number)->at_boundary())
                {
                    fe_face_values.reinit(cell, face_number);
                    fe_face_values[velocities].get_function_values(npaf_solution,face_velocity_values);
                    fe_face_values[velocities].get_function_gradients(npaf_solution, face_velocity_gradients);
                    fe_face_values[pressure].get_function_values(npaf_solution,
                                                                 face_pressure_values);
                    fe_face_values[pressure].get_function_gradients(npaf_solution,
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
                            
                            if (ifconsist_diff == true)
                            { diff_flux[i][0]=viscosity*(face_velocity_gradients[q_point][0][0]+face_velocity_gradients[q_point][0][0])*normal_vector[0]+viscosity*(face_velocity_gradients[q_point][0][1]+face_velocity_gradients[q_point][1][0])*normal_vector[1];
                                
                                diff_flux[i][1]=viscosity*(face_velocity_gradients[q_point][1][0]+face_velocity_gradients[q_point][0][1])*normal_vector[0]+viscosity*(face_velocity_gradients[q_point][1][1]+face_velocity_gradients[q_point][1][1])*normal_vector[1];
                                
                            }
                            else
                            {
                                diff_flux[i][0] = value_diffx;
                                diff_flux[i][1] = value_diffy;
                            }
                            
                            if (ifconsist_p == true)
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
            const ConstraintMatrix &constraints_used = zero_constraints;
            constraints_used.distribute_local_to_global(local_rhs,
                                                        local_dof_indices,
                                                        system_rhs);
            
        }
    }
    void NewtonSolve::linear_solve(BlockVector<double> &newton_update)
    {
        //call linear solver
        
        //system_matrix,system_rhs,newton_update
        
        const ConstraintMatrix &constraints_used = zero_constraints;
        
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
    
    void NewtonSolve::update_flow(BlockVector<double> &newton_update,BlockVector<double> &np1_solution,BlockVector<double> &np1_solution_time_derivative)
    {
        /*
        // update ut, pt
        ut_np1 = ut_np1+dut;
        pt_np1 = pt_np1+alphaf*gamma*dpt;
        
        // use ODE1 to update u_np1, p_np1
        u_np1 = u_np1+gamma*dt*dut;
        p_np1 = p_np1+alphaf*gamma*dt*dpt;
        */
        
        //npaf_solution = newton_update;
        //nonzero_constraints.distribute(npaf_solution);
        
        Vector<double> u_np1 (np1_solution.block(0).size());
        Vector<double> p_np1 (np1_solution.block(1).size());
        Vector<double> ut_np1 (np1_solution_time_derivative.block(0).size());
        Vector<double> pt_np1 (np1_solution_time_derivative.block(1).size());
        Vector<double> dut (np1_solution.block(0).size());
        Vector<double> dpt (np1_solution.block(1).size());
        
        u_np1 = np1_solution.block(0);
        ut_np1 = np1_solution_time_derivative.block(0);
        dut = newton_update.block(0);
        
        p_np1 = np1_solution.block(1);
        pt_np1 = np1_solution_time_derivative.block(1);
        dpt = newton_update.block(1);
        
        for(std::size_t i=0; i<ut_np1.size(); ++i)
        {
            ut_np1[i] = ut_np1[i]+dut[i];
            //std::cout<<u_npaf[i]<<std::endl;
            u_np1[i] = u_np1[i]+gamma*dt*dut[i];
            //std::cout<<ut_npam[i]<<std::endl;
        }
        
        for(std::size_t i=0; i<pt_np1.size(); ++i)
        {
            pt_np1[i] = pt_np1[i]+alphaf*gamma*dpt[i];
            //std::cout<<u_npaf[i]<<std::endl;
            p_np1[i] = p_np1[i]+alphaf*gamma*dt*dpt[i];
            //std::cout<<ut_npam[i]<<std::endl;
        }
        
        np1_solution.block(0) = u_np1;
        np1_solution.block(1) = p_np1;
        np1_solution_time_derivative.block(0) = ut_np1;
        np1_solution_time_derivative.block(1) = pt_np1;
    }

    
    void NewtonSolve::assemble_rhs(BlockVector<double> &npaf_solution,
                                   BlockVector<double> &npam_solution_time_derivative,
                                   double *** diffusive_k,
                                   double*** &div_SF_u,double**** &grad_SF_u,
                                   double**** &SF_u, double*** &SF_p,
                                   double**** &grad_SF_p,double*** &gijG)
    {
        setup_system(npaf_solution,npam_solution_time_derivative,diffusive_k,
                     div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG);
    }
    
    
    void NewtonSolve::getNewtonUpdate(BlockVector<double> &npaf_solution,
                                    BlockVector<double> &npam_solution_time_derivative,
                                    BlockVector<double> &newton_update,
                                    double *** diffusive_k,
                                    double*** &div_SF_u,double**** &grad_SF_u,
                                    double**** &SF_u, double*** &SF_p,
                                      double**** &grad_SF_p,double*** &gijG)
    {

        setup_system(npaf_solution,npam_solution_time_derivative,diffusive_k
                     ,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG);
        linear_solve(newton_update);
        
        /*
        Vector<double> du_n = newton_update.block(0);
        Vector<double> dp_n = newton_update.block(1);
        std::cout << "newton_update in newton" <<std::endl;
        for(std::size_t i=0; i<du_n.size(); ++i)
        {
            std::cout<<du_n[i]<<std::endl;
        }
        
        for(std::size_t i=0; i<dp_n.size(); ++i)
        {
            std::cout<<dp_n[i]<<std::endl;
        }
        */
    }
    
    void NewtonSolve::getResidual(double &current_res)
    {
        //iterPC(n_solution,n_solution_time_derivative,
         //      np1_solution,np1_solution_time_derivative);
        
        //assemble_rhs();
        //setup_system();
        
        
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


