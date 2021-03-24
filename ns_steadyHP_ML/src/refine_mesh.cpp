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

#include "../include/refine_mesh.h"
#include "../include/user_input.h"

#include <fstream>
#include <iostream>
#include <sstream>
namespace incompressible
{
    using namespace dealii;

    refine_mesh::refine_mesh(Triangulation<2> &triangulation,
                             hp::DoFHandler<2> &dof_handler,
                             BlockVector<double> &present_solution,
                             hp::FECollection<2> &fe)
    :
    triangulation(triangulation),
    dof_handler(dof_handler),
    present_solution(present_solution),
    fe(fe)
    {}
    
    void refine_mesh:: perfrom_refinement(hp::QCollection<1> &
                        face_quadrature_collection, int cycle,
                        SolutionTransfer<2, BlockVector<double>,
                        hp::DoFHandler<2>  > &solution_transfer,
                        Vector<float> &estimated_error_per_cell,
                        Vector<float> &true_error_per_cell)
    {
        //Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
        error_indicator(face_quadrature_collection,estimated_error_per_cell);
        
        Vector<float> smoothness_indicators(triangulation.n_active_cells());
        //estimate_smoothness(smoothness_indicators);
        
        savemesh(cycle,estimated_error_per_cell,true_error_per_cell,
                 smoothness_indicators);
        
        
        //GridRefinement::refine_and_coarsen_fixed_number(triangulation,
          //                                         estimated_error_per_cell,
            //                                        0.3,
              //                                      0.03);
        
        // change element degree based on smoothness
        /*
        float max_smoothness = *std::min_element(smoothness_indicators.begin(),
                                                 smoothness_indicators.end()),
        min_smoothness = *std::max_element(smoothness_indicators.begin(),
                                           smoothness_indicators.end());
        for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->refine_flag_set())
            {
                max_smoothness =
                std::max(max_smoothness,
                         smoothness_indicators(cell->active_cell_index()));
                min_smoothness =
                std::min(min_smoothness,
                         smoothness_indicators(cell->active_cell_index()));
            }
        const float threshold_smoothness = (max_smoothness + min_smoothness) / 2;
        std::cout<<"threshold_smoothness "<<threshold_smoothness<<std::endl;
        
        for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->refine_flag_set() &&
                (smoothness_indicators(cell->active_cell_index()) >
                 threshold_smoothness) &&
                (cell->active_fe_index() + 1 < fe.size()))
            {
                cell->clear_refine_flag();
                cell->set_active_fe_index(cell->active_fe_index() + 1);
            }
        */
        //triangulation.prepare_coarsening_and_refinement();
        //solution_transfer.prepare_for_pure_refinement();
        
        
        std::string inputfilename ="./input/refine_element/Shape1-1cycle";
        const std::string inputfile =
        inputfilename.c_str() + Utilities::int_to_string(cycle, 1) + ".txt";
        std::cout << "inputfilename "<<inputfile<<std::endl;
        
        std::ifstream MyReadFile(inputfile);
        std::cout << "read input refine_ele" <<std::endl;
        std::vector<int>         Refine_ele_list;
        // Use a while loop together with the getline() function to read the file line by line
        std::string myText;
        while (getline (MyReadFile, myText)) {
            int ele = std::stod(myText);
            //std::cout << "refine_ele "<<ele<<std::endl;
            Refine_ele_list.push_back(ele);
        }
        
        
        
        int cellnum = 0;
        std::vector<int>::iterator it;
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (Refine_ele_list[cellnum]==1)
                cell->set_refine_flag();
            cellnum = cellnum+1;
        }
        
        triangulation.prepare_coarsening_and_refinement();
        solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
        std::vector<bool> flag;
        triangulation.save_refine_flags(flag);
        
        //std::vector<bool> flag1;
        //triangulation.save_coarsen_flags(flag1);
        
        const std::string filename =
        "./output/mesh/refineFlag" + Utilities::int_to_string(cycle, 1) + ".txt";
        FILE * refineFile;
        refineFile = fopen(filename.c_str(), "w");
        
        int refine_ele = 0;
        for (int i = 0; i<flag.size();i++)
        {
            if (flag[i]>0)
            {
                refine_ele++;
                fprintf(refineFile,"%d \n",1);
            }
            else
            {
                //refine_ele++;
                fprintf(refineFile,"%d \n",0);
            }
            
            //std::cout << "refine_flag "<<flag[i]<<std::endl;
        }
        fclose(refineFile);
        
        
        const std::string errorfilename =
        "./output/mesh/estimated_error_per_cell" + Utilities::int_to_string(cycle, 1) + ".txt";
        FILE * errorFile;
        errorFile = fopen(errorfilename.c_str(), "w");
        for (int i = 0; i<estimated_error_per_cell.size();i++)
        {
            fprintf(errorFile,"%f \n",estimated_error_per_cell[i]);
        }
        fclose(errorFile);
        
        triangulation.execute_coarsening_and_refinement ();
        int total_ele = flag.size()/2-refine_ele/2+2*refine_ele;
        std::cout<<"adapt mesh finished "<<refine_ele<< " "<<flag.size()<<" "<<total_ele<<std::endl;
        
        std::cout << "   Number of active cells after adapt " << triangulation.n_active_cells()
        << std::endl;
        
    }
    
    void refine_mesh::error_indicator(hp::QCollection<1> &face_quadrature_collection,Vector<float> &estimated_error_per_cell)
    {
        
        if (KellyError)
          error_Kelly(face_quadrature_collection,estimated_error_per_cell);
        else{
            double Error_total=0;
            error_vms(estimated_error_per_cell,Error_total);
            std::cout<<"total_estimate_error "<<Error_total<<std::endl;
        }
        
    }
    
    void refine_mesh::error_Kelly(hp::QCollection<1> &face_quadrature_collection,
                                  Vector<float> &estimated_error_per_cell)
    {
        FEValuesExtractors::Vector velocity(0);
        KellyErrorEstimator<2>::estimate (dof_handler,
                                          face_quadrature_collection,
                                          typename FunctionMap<2>::type(),
                                          present_solution,
                                          estimated_error_per_cell,
                                          fe.component_mask(velocity));
    }
    
    void refine_mesh::error_vms(Vector<float> &estimated_error_per_cell,double &Error_total)
    {
        for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
        {
            quadrature_collection.push_back (QGauss<2>(degree*2));
        }
        
        hp::FEValues<2> hp_fe_values(fe,
                                     quadrature_collection,
                                     update_values | update_quadrature_points |
                                     update_JxW_values | update_gradients |
                                     update_inverse_jacobians);
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        hp::DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        int cellnum = 0;
        
        double* gij = new double[4];
        double rc,tauM,tauC;
        double c = 2/pow(3,0.5);
        //std::cout << "start cell iteration " <<c<<std::endl;
        
        Error_total=0;
        for (; cell != endc; ++cell)
        {
            cellnum++;
            const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
            hp_fe_values.reinit(cell);
            
            std::vector<double>         div_phi_u(dofs_per_cell);
            std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
            std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
            std::vector<double>         phi_p(dofs_per_cell);
            std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
            
            const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
            unsigned int n_q_points  = fe_values.n_quadrature_points;
            
            std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
            std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
            std::vector<double>       present_pressure_values(n_q_points);
            std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
            
            fe_values[velocities].get_function_values(present_solution,present_velocity_values);
            fe_values[velocities].get_function_gradients(present_solution, present_velocity_gradients);
            fe_values[pressure].get_function_values(present_solution,
                                                    present_pressure_values);
            fe_values[pressure].get_function_gradients(present_solution,
                                                       present_pressure_gradients);
            std::vector<Tensor<1, 2> > rm(n_q_points);
            
            int celldeg = cell->active_fe_index()+1;
            double L2_error = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                const DerivativeForm<1, 2, 2> &JacobiInverse = fe_values.inverse_jacobian(q);
                
                getGij(JacobiInverse,gij);
                
                rm[q] = present_velocity_gradients[q]*present_velocity_values[q]+present_pressure_gradients[q];//-sourceTerm[q];
                
                rc = present_velocity_gradients[q][0][0]+present_velocity_gradients[q][1][1];
                
                getTaum(celldeg,gij,present_velocity_values[q][0],present_velocity_values[q][1],tauM,tauC);
                
                double magrm = rm[q][0]*rm[q][0]+rm[q][1]*rm[q][1];
                L2_error = L2_error + c*tauM*c*tauM*magrm*fe_values.JxW(q);
                //std::cout << "rm " <<magrm<<" "<<tauM<<" "<<fe_values.JxW(q)<<std::endl;
                Error_total = Error_total+c*tauM*c*tauM*magrm*fe_values.JxW(q);
            }
            //L2_error = pow(L2_error,0.5);
            estimated_error_per_cell[cellnum-1] = pow(L2_error,0.5);
            //std::cout<<" L2_est_u "<<estimated_error_per_cell[cellnum-1]<<std::endl;
        }
        Error_total = pow(Error_total,0.5);
    }
    
    void refine_mesh::getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij)
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
    
    void refine_mesh::getTaum(int degree,double* gij,double uele,double vele,double &tauM, double &tauC)
    {
        double viscosity = 1/Re;
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
        tauC = (1)/(3*(degree*degree)*tauM*(gij0+gij3));
        //std::cout << "degree  " << degree<<std::endl;
        //std::cout << "A22  " << A22<<std::endl;
    }

    
    void refine_mesh::estimate_smoothness(Vector<float> &smoothness_indicators)
    {
        const unsigned int N = max_degree;
        QGauss<1>      base_quadrature(2);
        QIterated<2> quadrature(base_quadrature, N);
        for (unsigned int i = 0; i < fe.size(); i++)
            fourier_q_collection.push_back(quadrature);
        
        fourier = std::make_shared<FESeries::Fourier<2>>(N,fe,fourier_q_collection);
        resize(fourier_coefficients, N);
        
        Vector<double> local_dof_values;
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            
            local_dof_values.reinit(cell->get_fe().dofs_per_cell);
            cell->get_dof_values(present_solution, local_dof_values);
            
            fourier->calculate(local_dof_values,
                               cell->active_fe_index(),
                               fourier_coefficients);
            
            std::pair<std::vector<unsigned int>, std::vector<double>> res =
            FESeries::process_coefficients<2>(fourier_coefficients,
                                std::bind(&refine_mesh::predicate,this,std::placeholders::_1),
                                VectorTools::Linfty_norm);
            
            Assert(res.first.size() == res.second.size(), ExcInternalError());
            
            if (ln_k.size() == 0)
            {
                ln_k.resize(res.first.size(), 0);
                for (unsigned int f = 0; f < ln_k.size(); f++)
                    ln_k[f] = std::log(2.0 * numbers::PI * std::sqrt(1. * res.first[f]));
            }
            
            for (double &residual_element : res.second)
                residual_element = std::log(residual_element);
            std::pair<double, double> fit = FESeries::linear_regression(ln_k, res.second);
            smoothness_indicators(cell->active_cell_index()) = -fit.first - 1. * dim / 2;
            
        }
    }
    void refine_mesh::savemesh(int cycle,Vector<float> estimated_error_per_cell,
                               Vector<float> true_error_per_cell,
                               Vector<float> smoothness_indicators)
    {
        
        //int dim = 2;
         Vector<float> error_efficiency_per_cell (triangulation.n_active_cells());
         Vector<float> fe_degrees(triangulation.n_active_cells());
        
         int cellnum = 0;
         for (const auto &cell : dof_handler.active_cell_iterators())
         {
             cellnum++;
             fe_degrees(cell->active_cell_index()) =fe[cell->active_fe_index()].degree;
             error_efficiency_per_cell[cellnum-1] =estimated_error_per_cell[cellnum-1]/true_error_per_cell[cellnum-1];
             //std::cout<<" estimated_error_per_cell "<<estimated_error_per_cell[cellnum-1]<<std::endl;
             //std::cout<<" true_error_per_cell "<<true_error_per_cell[cellnum-1]<<std::endl;
             //std::cout<<" error_efficiency_per_cell "<<error_efficiency_per_cell[cellnum-1]<<std::endl;
         }
        
        double estimated_error = estimated_error_per_cell.l2_norm();
        
        std::cout << "******************************" << std::endl;
        std::cout << " The estimated_error is " << estimated_error << std::endl;
        
        double true_error = true_error_per_cell.l2_norm();
        
        std::cout << "******************************" << std::endl;
        std::cout << " The true_error is " << true_error << std::endl;
        
        std::cout << "******************************" << std::endl;
        std::cout << " The error_efficiency is " << estimated_error/true_error << std::endl;
        
         DataOut<2, hp::DoFHandler<2>> data_out;
         data_out.attach_dof_handler(dof_handler);
         data_out.add_data_vector(present_solution, "solution");
         data_out.add_data_vector(estimated_error_per_cell, "estimate_error");
         //data_out.add_data_vector(true_error_per_cell, "true_error");
         //data_out.add_data_vector(error_efficiency_per_cell, "error_efficiency");
         //data_out.add_data_vector(smoothness_indicators, "smoothness");
         //data_out.add_data_vector(fe_degrees, "fe_degree");
         data_out.build_patches();
        
         const std::string filename =
         vtkfilename.str().c_str() + Utilities::int_to_string(cycle, 2) + ".vtk";
         std::ofstream output(filename);
         data_out.write_vtk(output);
    }
    
    template <typename T>
    void refine_mesh::resize(Table<2, T> &coeff, const unsigned int N)
    {
        TableIndices<2> size;
        for (unsigned int d = 0; d < 2; d++)
            size[d] = N;
        coeff.reinit(size);
    }
    
    std::pair<bool, unsigned int>
    refine_mesh::predicate(const TableIndices<2> &ind)
    {
        unsigned int dim = 2;
        unsigned int v = 0;
        for (unsigned int i = 0; i < dim; i++)
            v += ind[i] * ind[i];
        if (v > 0 && v < max_degree * max_degree)
            return std::make_pair(true, v);
        else
            return std::make_pair(false, v);
    }
    
}


