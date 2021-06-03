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

#include "../include/dynamic.h"
//#include "../include/newton.h"
#include "../include/Eigen/Dense"
#include "../include/user_input.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
namespace incompressible
{
    using namespace dealii;
    using namespace Eigen;
    using namespace std;
    dynamic::dynamic(double viscosity,
                             DoFHandler<2> &dof_handler,
                             FESystem<2> &fe,
                             Triangulation<2> &triangulation)
    :
    viscosity(viscosity),
    dof_handler(dof_handler),
    fe(fe),
    triangulation(triangulation)
    {}
    
    void dynamic::getInitialtau(int Grid_size,double vmax,double *cellArea, std::vector<std::vector<DoFHandler<2>::
    active_cell_iterator> > &v_to_e_list,BlockVector<double> &n_solution,
                                double *taumELE, double **taum)
    {
        //taumELE = ones(Grid_size,1)*0.002;
        //std::cout << "Grid_size "<<Grid_size<<std::endl;
        
        for (int i = 0; i<Grid_size;i++)
        {
            //Area = cellArea[i];
            taumELE[i] =0.33*cellArea[i];
            for (int j = 0; j<4;j++)
                taum[i][j] = 0.33*cellArea[i];
        }
        
        /*
        QGauss<2> quadrature_formula(degree*2);
        
        DoFHandler<2>::active_cell_iterator cell;
        auto vertices = triangulation.get_vertices();
        double *c1node = new double [vertices.size()];
        double *c2Dist = new double [vertices.size()];
        
        FEValues<2> fe_values(fe,quadrature_formula,
                    update_values | update_quadrature_points |
                    update_JxW_values | update_gradients |
                        update_inverse_jacobians);
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        Vector<float> fe_degrees(triangulation.n_active_cells());
        
        unsigned int n_q_points  = fe_values.n_quadrature_points;
        std::vector<int> pH;
        std::vector<double> PatchArea;
        std::vector<std::vector<double> > uEle;
        std::vector<std::vector<double> > vEle;
        
        for (unsigned long int node = 0; node<vertices.size(); node++)
        {
            c1node[node] = 1.0/vmax;
            int nElement = v_to_e_list[node].size();
            pH.resize(nElement);
            PatchArea.resize(nElement);
            //std::cout <<" node "<<node<< " nElement "<<nElement<<std::endl;
            uEle.resize(nElement);
            vEle.resize(nElement);
            for (int j = 0; j<nElement; j++)
            {
                cell = v_to_e_list[node][j];
                fe_values.reinit(cell);
                int cell_index = cell->active_cell_index();
                //std::cout <<" cell_index "<<cell_index<<std::endl;
                // get least square H grid constructions
                pH[j] = fe[cell->active_fe_index()].degree;
                PatchArea[j] = cellArea[cell_index];
                //std::cout <<" pH "<<pH[j]<<" PatchArea "<<PatchArea[j]<<std::endl;
                
                //fe_values.reinit(cell);
                //unsigned int n_q_points  = fe_values.n_quadrature_points;
                std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
                fe_values[velocities].get_function_values(n_solution,present_velocity_values);
                //std::cout <<" n_q_points "<<n_q_points<<std::endl;
                for (unsigned int q=0; q<n_q_points; ++q)
                {
                    //std::cout <<" present_velocity_values "<<present_velocity_values[q]<<std::endl;
                    // get u,v,p at quadrature points
                    uEle[j].push_back(present_velocity_values[q][0]);
                    vEle[j].push_back(present_velocity_values[q][1]);
                    //std::cout <<" uEle "<<uEle[j][q]<<std::endl;
                }
         } // end patch element loop
        int patchOrder = *max_element(pH.begin(), pH.end());
        int c2 = determinc2(uEle,vEle,patchOrder,PatchArea,nElement);
        c2Dist[node] = c2;
        //std::cout <<node<< " c1node[node] "<<c1node[node] <<c2Dist[node]<<std::endl;
        }
        
        cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        int cellnum = 0;
        for (; cell != endc; ++cell)
         {
            int cell_index = cell->active_cell_index();
            double h = sqrt(cellArea[cell_index]);
            
            //taumELE[cellnum] = 0.02;
            //std::cout << "taumELE "<<taumELE[i]<<std::endl;
            for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
            {
                //taum[i][v] = 0.02;
                int v_index = cell->vertex_index(v);
                taum[cellnum][v]= c1node[v_index]*pow(h,c2Dist[v_index]);
                //std::cout << "c1node[v_index] "<<c1node[v_index]<<c2Dist[v_index]<<std::endl;
            }
             
            double *i1 = std::max_element(taum[cellnum], taum[cellnum] + 4);
            //std::cout << " taum max "<<*i1<<std::endl;
            taumELE[cellnum] = *i1;
             
            cellnum++;
        }
         */
    }
    
    void dynamic::getTave(int Grid_size, BlockVector<double> &n_solution, double *Tave)
    {
        for (int i = 0; i<Grid_size; i++)
        {
            Tave[i] = 0.1;
        }
        /*
        QGauss<2> quadrature_formula(degree*2);
        
        FEValues<2> fe_values(fe,quadrature_formula,
                              update_values | update_quadrature_points |
                              update_JxW_values);
        const FEValuesExtractors::Vector velocities(0);
        
        // compute element-wise tau, loop through element
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        int cellnum = 0;
        double* gij = new double[4];
        double uquad; double vquad;
        double tau1sqinv = 0;
        double C1 = 1;
        double taut = (2*C1/dt)*(2*C1/dt);
        for (; cell != endc; ++cell)
        {
            cellnum++;
            fe_values.reinit(cell);
            const unsigned int n_q_points    = quadrature_formula.size();
            std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
            
            fe_values[velocities].get_function_values(n_solution,present_velocity_values);
            double TT = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                uquad = present_velocity_values[q][0];
                vquad = present_velocity_values[q][1];
                gij[0] = gijG[cellnum-1][q][0];
                gij[1] = gijG[cellnum-1][q][1];
                gij[2] = gijG[cellnum-1][q][2];
                gij[3] = gijG[cellnum-1][q][3];
                
                tau1sqinv = uquad*(uquad*gij[0]+vquad*gij[2])+vquad*(uquad*gij[1]+vquad*gij[3]);
                TT = TT+10/sqrt((taut+tau1sqinv));
            }
            Tave[cellnum-1] = TT/n_q_points;
            
        } // end loop for taumele_g
        */
        
    }
    
    void dynamic::relxation_ele(int Grid_size, double **taum_g,double *Tave,double **taum, double *taumEle)
    {
        
        for (int ele = 0; ele<Grid_size; ele++)
        {
            for (int i = 0; i<4; i++)
            {
                taum[ele][i] = (taum[ele][i]-taum_g[ele][i])*exp(-dt/Tave[ele])+taum_g[ele][i];
                //std::cout << " taum update "<<i<<" "<<taum_g[ele][i]<<" "<<taum[ele][i]<<std::endl;
            }
            double *i1 = std::max_element(taum[ele], taum[ele] + 4);
            //std::cout << " taum max "<<*i1<<std::endl;
            taumEle[ele] = *i1;
        }
    
    
    }
    
    void dynamic::getBCinfo(std::vector<int> &IBC, double *bcTag)
    {
        // loop through element
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        //QGauss<2> quadrature_formula(degree*2);
        
        //FEValues<2> fe_values(fe,quadrature_formula,
          //                    update_values | update_quadrature_points |
            //                  update_JxW_values);
        int cellnum = 0;
        
        // get IBC and BCtag
        cellnum = 0;
        for (; cell != endc; ++cell)
        {
            cellnum++;
            for (unsigned int face_number = 0;
                 face_number < GeometryInfo<2>::faces_per_cell;
                 ++face_number)
            {
                if (cell->face(face_number)->at_boundary())
                {
                    int mID = cell->face(face_number)->boundary_id();
                    int vID0 = cell->face(face_number)->vertex_index(0);
                    int vID1 = cell->face(face_number)->vertex_index(1);
                    //Point<2> &v1 = cell->face(face_number)->vertex(0);
                    //Point<2> &v2 = cell->face(face_number)->vertex(1);
                    IBC.push_back(vID0);
                    IBC.push_back(vID1);
                    bcTag[vID0] = mID;
                    bcTag[vID1] = mID;
                    //std::cout<<"mID "<<vID0<<" "<<vID1<<" "<<mID<<std::endl;
                }
            }
        }// end cell loop
        
        //IBC.erase( unique( IBC.begin(), IBC.end() ), IBC.end() );
        /*
        for (int i = 0; i<IBC.size(); i++)
          std::cout<<"IBC "<<IBC[i]<<std::endl;
        
        auto vertices = triangulation.get_vertices();
        for (int i = 0; i<vertices.size(); i++)
            std::cout<<"i "<<bcTag[i]<<std::endl;
        */
    }
    
    void dynamic::additionalDynamicinfo(std::vector<std::vector<DoFHandler<2>::
                                        active_cell_iterator> > &v_to_e_list,
                                        double *cellArea,
                                        double **lhsG, double **rhsG,
                                        double **xyuHG, double**xyduHxG, double **xyduHyG)
    {
            
            DoFHandler<2>::active_cell_iterator cell;
            auto vertices = triangulation.get_vertices();
        
            QGauss<2> quadrature_formula(degree*2);
            FEValues<2> fe_values(fe,quadrature_formula,
                                  update_values | update_quadrature_points |
                                  update_JxW_values | update_gradients |
                                  update_inverse_jacobians);
            
            Vector<float> fe_degrees(triangulation.n_active_cells());
        
            unsigned int n_q_points  = fe_values.n_quadrature_points;
            //std::cout<<" n_q_points "<<n_q_points<<std::endl;
        
            int sizeA = (degree+1)*(degree+2)/2;
            int sizeElenum = 10;
            for(unsigned long i = 0; i <vertices.size(); i++)
            {
                lhsG[i] = new double[sizeA*sizeA];
                rhsG[i] = new double[sizeA*n_q_points*sizeElenum];
                xyuHG[i] = new double[sizeElenum*n_q_points*sizeA];
                xyduHxG[i] = new double[sizeElenum*n_q_points*sizeA];
                xyduHyG[i] = new double[sizeElenum*n_q_points*sizeA];
            }
        
            double *xyuH = new double[sizeA];
            double *xyduHx = new double[sizeA];
            double *xyduHy = new double[sizeA];
        
            for (unsigned long int i = 0; i<vertices.size(); i++)
            {
                //std::cout << "----------v -------------  " << i<<std::endl;
                
                // loop though patch elements get least square integral as H grid construction
                int nElement = v_to_e_list[i].size();
                
                std::vector<std::vector<double> > xEle;
                std::vector<std::vector<double> > yEle;
                xEle.resize(nElement);
                yEle.resize(nElement);
                
                std::vector<std::vector<double> > JxEle;
                JxEle.resize(nElement);
                
                std::vector<int> pH;
                std::vector<double> PatchArea;
                pH.resize(nElement);
                PatchArea.resize(nElement);
                
                for (int j = 0; j<nElement; j++)
                {
                    cell = v_to_e_list[i][j];
                    int cell_index = cell->active_cell_index();
                    //std::cout << "v_to_e_indices  "<<cell_index<<std::endl;
                    
                    //const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
                    //std::cout<<" dofs_per_cell "<<dofs_per_cell<<std::endl;
                    fe_values.reinit(cell);
                    
                    //const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
                    //unsigned int n_q_points  = fe_values.n_quadrature_points;
                    //std::cout<<" n_q_points "<<n_q_points<<std::endl;
                    
                    // get least square H grid constructions
                    pH[j] = fe[cell->active_fe_index()].degree;
                    //long double area = 0;
                    //std::cout<<" celldeg "<<" "<<fe_degrees[cell_index]<<std::endl;
                    for (unsigned int q=0; q<n_q_points; ++q)
                    {
                        //std::cout<<"present_velocity_values "<<uEle[j][q]<<std::endl;
                        // get x y at quadrature points
                        const Point<2> &vPoint=fe_values.quadrature_point(q);
                        xEle[j].push_back(vPoint[0]);
                        yEle[j].push_back(vPoint[1]);
                        //std::cout << "vPoint  " << vPoint[0]<<" "<< vPoint[1]<<std::endl;
                        //double Jx = fe_values.JxW(q);
                        JxEle[j].push_back(fe_values.JxW(q));
                        //area += fe_values.JxW (q);
                    }
                    
                    PatchArea[j] = cellArea[cell_index];
                    //std::cout << "cellArea  " << cellArea[cell_index]<<" "<< area<<std::endl;
                    
                } // end patch element loop
                int patchOrder = *max_element(pH.begin(), pH.end());

                std::vector<double> xmaxEle;
                std::vector<double> xminEle;
                std::vector<double> ymaxEle;
                std::vector<double> yminEle;
                xmaxEle.resize(nElement);
                xminEle.resize(nElement);
                ymaxEle.resize(nElement);
                yminEle.resize(nElement);
                
                for (int j = 0; j<nElement; j++)
                {
                    xmaxEle[j] = *max_element(xEle[j].begin(), xEle[j].end());
                    xminEle[j] = *min_element(xEle[j].begin(), xEle[j].end());
                    ymaxEle[j] = *max_element(yEle[j].begin(), yEle[j].end());
                    yminEle[j] = *min_element(yEle[j].begin(), yEle[j].end());
                    
                }
                double xmax = *max_element(xmaxEle.begin(), xmaxEle.end());
                double xmin = *min_element(xminEle.begin(), xminEle.end());
                double ymax = *max_element(ymaxEle.begin(), ymaxEle.end());
                double ymin = *min_element(yminEle.begin(), yminEle.end());
                
                //Point<2> &vcord_u = vertices[vertex_value_u_map[i][1]];
                
                //FullMatrix<double> LSlhs;
                int sizeA = (patchOrder+1)*(patchOrder+2)/2;
                double *lhs1D = new double[sizeA*sizeA];
                double *rhs1D = new double[sizeA*n_q_points*nElement];
                getLSintlhs(xmax,xmin,ymax,ymin,nElement,xEle,yEle,patchOrder,JxEle,lhs1D,rhs1D);
                
                for ( int lv = 0; lv<sizeA*sizeA;lv++)
                    lhsG[i][lv] = lhs1D[lv];
                for (unsigned int lv = 0; lv<sizeA*n_q_points*nElement;lv++)
                    rhsG[i][lv] = rhs1D[lv];
                
                int sizeElenum = v_to_e_list[i].size();
                int count = 0;
                for (int j = 0; j<sizeElenum; j++)
                {
                    cell = v_to_e_list[i][j];
                    fe_values.reinit(cell);
                    
                    //const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
                    //unsigned int n_q_points  = fe_values.n_quadrature_points;
                
                    for (unsigned int q=0; q<n_q_points; ++q)
                    {
                        //std::cout <<"q"<<q<<std::endl;
                        // get x y at quadrature points
                        const Point<2> &vPoint = fe_values.quadrature_point(q);
                        
                        //double *xyuH = new double[sizeA];
                        //double *xyduHx = new double[sizeA];
                        //double *xyduHy = new double[sizeA];
                        
                        calcuHxycoeff(vPoint,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy);
                        
                        for (int np=0; np<sizeA;np++)
                        {
                            xyuHG[i][count] = xyuH[np];
                            xyduHxG[i][count] = xyduHx[np];
                            xyduHyG[i][count] = xyduHy[np];
                            //std::cout <<"i "<<i<<" np "<<np<< " xyuHG  " << xyuHG[i][count]<<std::endl;
                            count = count+1;
                        }
                    
                    }// end quadrature
                    
                } //end patch iteration
                
            } // end vertex iteration
        
    }
    
    void dynamic::buildDynamicMeshinfo(std::map<int,int> &vertices_to_dofs,
                                       std::vector<std::vector<DoFHandler<2>::active_cell_iterator> > &v_to_e_list,
                                       std::vector<std::vector<int> > &v_to_e_indices,
                                       AffineConstraints<double> &nonzero_constraints,
                                       double **vertex_value_u_map,double **vertex_value_v_map,
                                       double **vertex_value_p_map,double *cellArea,
                                       std::vector<int> &IBC, double *bcTag)
    {
        //int *v_flag;
        auto vertices = triangulation.get_vertices();
        v_flag = new int [vertices.size()];

        // loop through element
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        QGauss<2> quadrature_formula(degree*2);
        
        FEValues<2> fe_values(fe,quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_JxW_values);
        int cellnum = 0;
        
        getBCinfo(IBC, bcTag);
        // loop through element
        for (; cell != endc; ++cell)
        {
            cellnum++;
            // what about hanging nodes????
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                //Point<2> &v = cell->vertex(i);
                //cell->vertex(i).set_manifold(10);
                //cell->vertex(i)->get_manifold();
                //std::cout << "coordinate  " << v[0]<<" "<< v[1]<<std::endl;
                vertices_to_dofs[cell->vertex_index(i)] = cell->vertex_dof_index(i,0,0);
                int v_index = cell->vertex_index(i);
                v_flag[v_index] = 0;
                //std::cout << "v_index  "<<v_index<<" bcTag "<<bcTag[v_index]<<std::endl;
                if (bcTag[v_index] == dynamicBC )
                {
                    v_flag[v_index] = 1;
                    //std::cout << "v_index  "<<v_index<<" v_flag "<<v_flag[v_index]<<std::endl;
                }
                //std::cout << "v_index  "<<v_index<<" dof "<<cell->vertex_dof_index(i,0,0)<<std::endl;
                v_to_e_list[v_index].push_back(cell);
                v_to_e_indices[v_index].push_back(cell->active_cell_index());
                //std::cout << "v_to_e_indices  " << cell->index()<<std::endl;
            }
            
            fe_values.reinit(cell);
            
            //const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
            unsigned int n_q_points  = fe_values.n_quadrature_points;
            double area = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                // get cell area
                area += fe_values.JxW (q);
            }
            cellArea[cellnum-1] = area;
            //std::cout << "cellArea  "<<cellArea[cellnum-1]<<std::endl;
        } //end cell iteration
        
        
        double **DOF_value;
        DOF_value = new double *[dof_handler.n_dofs()];
        for(unsigned long i = 0; i <dof_handler.n_dofs(); i++)
        {
            DOF_value[i] = new double[2];
        }
        
        for (unsigned long int i = 0; i<dof_handler.n_dofs(); i++)
        {
            DOF_value[i][0] = 0;
            DOF_value[i][1] = 0;
            //std::cout << "DOF_value " << i<<" "<<DOF_value[i][0]<<" "<<DOF_value[i][1]<<std::endl;
        }
        
        std::cout << "nonzero_constraints lines in dynamic"<<std::endl;
        for (const auto &line : nonzero_constraints.get_lines())
        {
            //add boundary DOFs
            DOF_value[line.index][0] = 1.0;
            DOF_value[line.index][1] = line.inhomogeneity;
            //std::cout << "DoF in dynamic" << line.index<<std::endl;
            // remove hanging nodes, they are not at boundaries
            for (const auto &entry : line.entries)
            {
                DOF_value[line.index][0] = 0;
                //std::cout << "DoF " << line.index <<" "<<entry.first<<std::endl;
            }
            
        }
        /*
        for (unsigned long int i = 0; i<dof_handler.n_dofs(); i++)
        {
            std::cout << "DOF_value " << i<<" "<<DOF_value[i][0]<<" "<<DOF_value[i][1]<<std::endl;
        }
        */
        //std::cout << "vertices.size() " <<vertices.size()<<std::endl;
        
        int verticesSize = vertices.size();
        // get bcnode fit
        for (int i = 0; i<verticesSize; i++)
        {
            //std::cout << "v  " << i<<std::endl;
            //vertex_value_u_map[i][0] = 0;
            Point<2> &v = vertices[i];
            double xc = v[0];
            double yc = v[1];
            //std::cout << "xc  " << xc<<" yc "<<yc<<std::endl;
            
            DoFHandler<2>::active_cell_iterator cell;
            int nElement = v_to_e_list[i].size();
            // get surronding element
            //int count = 0;
            std::vector<double> bc_test_value_u;
            std::vector<int> bc_test_vertex_u;
            
            // output flag, fit_vertex_index and fit_vertex_value
            bool vertex_flag_u; int max_index_u;
            bool vertex_flag_v; int max_index_v;
            bool vertex_flag_p; int max_index_p;
            getBCnode(nElement,v_to_e_list,DOF_value,i,vertices_to_dofs,verticesSize,
                       xc, yc,vertex_flag_u,max_index_u,vertex_flag_v,max_index_v,
                      vertex_flag_p,max_index_p);
            
            vertex_value_u_map[i][0] = vertex_flag_u;
            vertex_value_u_map[i][1] = max_index_u; // vertex
            vertex_value_u_map[i][2] = 0; //need to convert to DOF
            //std::cout << "max_index_u  " << max_index_u<<std::endl;
            if (vertex_flag_u>0)
            {
                int max_index_u_dof = vertices_to_dofs[max_index_u];
                vertex_value_u_map[i][2] = DOF_value[max_index_u_dof][1];
            }
            
            vertex_value_v_map[i][0] = vertex_flag_v;
            vertex_value_v_map[i][1] = max_index_v;
            vertex_value_v_map[i][2] = 0;
            if (vertex_flag_v>0)
            {
                int max_index_u_dof = vertices_to_dofs[max_index_v];
                int max_index_v_dof = max_index_u_dof+1;
                vertex_value_v_map[i][2] = DOF_value[max_index_v_dof][1];
            }

            vertex_value_p_map[i][0] = vertex_flag_p;
            vertex_value_p_map[i][1] = max_index_p;
            vertex_value_p_map[i][2] = 0;
            if (vertex_flag_p>0)
            {
                int max_index_u_dof = vertices_to_dofs[max_index_p];
                int max_index_p_dof = max_index_u_dof/2+2*verticesSize;
                vertex_value_p_map[i][2] = DOF_value[max_index_p_dof][1];
            }
            
            //if (vertex_value_u_map[i][0]>0){
              //  std::cout << "xc  " << xc<<" yc "<<yc<<" i "<<i<<std::endl;
                //std::cout << "vertex_value_u_map_value  " << vertex_value_u_map[i][0]<<" "<<vertex_value_u_map[i][1]<<" "<<vertex_value_u_map[i][2]<<std::endl;
           // }
            //std::cout << "vertex_value_v_map_value  " << vertex_value_v_map[i][0]<<" "<<vertex_value_v_map[i][1]<<" "<<vertex_value_v_map[i][2]<<std::endl;
            //std::cout << "vertex_value_p_map_value  " << vertex_value_p_map[i][0]<<" "<<vertex_value_p_map[i][1]<<" "<<vertex_value_p_map[i][2]<<std::endl;
            
        }//end vertex loop
        
    }
    
    void dynamic::getBCnode(int nElement,
                            std::vector<std::vector<DoFHandler<2>::
                            active_cell_iterator> > &v_to_e_list,
                            double **DOF_value,int i,
                            std::map<int,int> &vertices_to_dofs,
                            int verticesSize,
                            double xc, double yc,
                            bool &vertex_flag_u,int  &max_index_u,
                            bool &vertex_flag_v,int  &max_index_v,
                            bool &vertex_flag_p,int  &max_index_p)
    {
        //double tmp = 0;
        DoFHandler<2>::active_cell_iterator cell;
        
        std::vector<int> bc_test_flag_u;
        std::vector<int> bc_test_vertex_u;
        std::vector<double> bc_test_value_u;
        
        std::vector<int> bc_test_flag_v;
        std::vector<int> bc_test_vertex_v;
        std::vector<double> bc_test_value_v;
        
        std::vector<int> bc_test_flag_p;
        std::vector<int> bc_test_vertex_p;
        std::vector<double> bc_test_value_p;
        
        int ui_vertex_dof = vertices_to_dofs[i];
        int vi_vertex_dof = ui_vertex_dof+1;
        int pi_vertex_dof = ui_vertex_dof/2+2*verticesSize;
        
        for (int j = 0; j<nElement; j++)
        {
            cell = v_to_e_list[i][j];
            // get all vertex of surrounding cell
            for (unsigned int k=0; k<GeometryInfo<2>::vertices_per_cell; ++k)
            {
                int v_index = cell->vertex_index(k);
                int u_vertex_dof = vertices_to_dofs[v_index];
                int v_vertex_dof = u_vertex_dof+1;
                int p_vertex_dof = u_vertex_dof/2+2*verticesSize;
                
                Point<2> &v = cell->vertex(k);
                int u_flag = DOF_value[u_vertex_dof][0];
                int v_flag = DOF_value[v_vertex_dof][0];
                int p_flag = DOF_value[p_vertex_dof][0];
                
                // get bc_test_value_u and bc_test_vertex_u for u, v and p
                getbc_test_value(u_flag,v_index,i,ui_vertex_dof,DOF_value,v,xc,yc,
                                 u_vertex_dof,v_vertex_dof,
                                 bc_test_flag_u,bc_test_vertex_u,bc_test_value_u);
                
                getbc_test_value(v_flag,v_index,i,vi_vertex_dof,DOF_value,v,xc,yc,
                                 u_vertex_dof,v_vertex_dof,
                                 bc_test_flag_v,bc_test_vertex_v,bc_test_value_v);
                
                getbc_test_value(p_flag,v_index,i,pi_vertex_dof,DOF_value,v,xc,yc,
                                 u_vertex_dof,v_vertex_dof,
                                 bc_test_flag_p,bc_test_vertex_p,bc_test_value_p);
                
            }
        }//finish element loop
        
        vertex_flag_u = !bc_test_flag_u.empty();
        vertex_flag_v = !bc_test_flag_v.empty();
        vertex_flag_p = !bc_test_flag_p.empty();
        //vertex_flag_u = vertex_flag_p;
        //vertex_flag_v = vertex_flag_p;
        
        //std::cout << "flag  " << vertex_flag_u<<std::endl;
        
        // get max_index_u, max_index_v, max_index_p
        //int max_index_u = -1; int max_index_v = -1; int max_index_p = -1;
         max_index_u = getMaxIndex(vertex_flag_u,bc_test_vertex_u,bc_test_value_u);
         max_index_v = getMaxIndex(vertex_flag_v,bc_test_vertex_v,bc_test_value_v);
         max_index_p = getMaxIndex(vertex_flag_p,bc_test_vertex_p,bc_test_value_p);
        //std::cout << "index  " << max_index_u<<max_index_v<<max_index_p<<std::endl;
        
    }
    
    void dynamic::getbc_test_value(bool flag, int v_index, int i, int ui_vertex_dof,
                                   double **DOF_value,Point<2> &v,
                                   double xc,double yc,
                                   int u_vertex_dof,
                                   int v_vertex_dof,
                                   std::vector<int> &bc_test_flag,
                                   std::vector<int> &bc_test_vertex,
                                   std::vector<double> &bc_test_value)
    {
        //int tmp;
        int node_flag = DOF_value[ui_vertex_dof][0];
        int dynamic_bc = v_flag[v_index];
        //std::cout << " v_flag " << v_flag[v_index]<<" vx "<<v[0]<<std::endl;
        
        if (flag && v_index!= i && !node_flag && dynamic_bc )
        {
            // it is adjacent to a bc_node but not a bc_node
            int flag_u = 1;
            //double bc_value_u = DOF_value[u_vertex_dof][1];
            //double bc_value_v = n_solution[v_vertex_dof];
            
            double xs = v[0];
            double ys = v[1];
            double Dx = xs-xc;
            double Dy = ys-yc;
            double Ds = sqrt(Dx*Dx+Dy*Dy);
            /*
            double velocity_norm = sqrt(bc_value_u*bc_value_u+bc_value_v*bc_value_v);
            if (Ds*velocity_norm!=0)
            {
                tmp = (Dx*bc_value_u+Dy*bc_value_v)/(Ds*velocity_norm);
            }
            else
            {
                tmp = 0;
            }
            */
            //tmp = Ds;
            //count = count+1;
            bc_test_flag.push_back(flag_u);
            bc_test_value.push_back(Ds);
            bc_test_vertex.push_back(v_index);
            //std::cout << "flag_u  " << flag_u<<std::endl;
            //std::cout << "xs  " << xs<<" ys "<<ys<<" xc "<<xc<<" yc "<<yc<<std::endl;
            //std::cout << "Dx  " << Dx<<" Dy "<<Dy<<" Ds "<<Ds<<std::endl;
            //std::cout << "v_index  " << v_index<<" i "<<i<<std::endl;
        }
    }
    
    int dynamic::getMaxIndex(int vertex_flag,
                              std::vector<int> bc_test_vertex,
                              std::vector<double> bc_test_value)
    {
        //std::cout << "vertex_flag  " << vertex_flag<<std::endl;
        if (vertex_flag)
        {
            int min_index = bc_test_vertex[0];
            double min_value = bc_test_value[0];
        
            for (unsigned long tmpcount = 0; tmpcount<bc_test_value.size();tmpcount++)
            {
                
                if (bc_test_value[tmpcount]<min_value)
                {
                    min_value = bc_test_value[tmpcount];
                    min_index = bc_test_vertex[tmpcount];
                }
            }
            //std::cout << "max_index  " << max_index<<" max_value "<<max_value<<std::endl;
            return min_index;
            
        }
        return -1;
    }
    void dynamic::getc1c2(BlockVector<double> &n_solution,
                          BlockVector<double> &n_solution_time_derivative,
                          std::vector<std::vector<DoFHandler<2>::
                        active_cell_iterator> > &v_to_e_list,
                        double **vertex_value_u_map,
                        double **vertex_value_v_map,
                        double **vertex_value_p_map,
                        double *cellArea, std::vector<int> IBC,
                        double **lhsG, double **rhsG,
                        double **xyuHG, double**xyduHxG, double **xyduHyG,
                        double **taumele_g)
    {
        /*
         Vector<double> du_n = n_solution.block(0);
         Vector<double> dp_n = n_solution.block(1);
         std::cout << "n_solution " <<std::endl;
         for(std::size_t i=0; i<du_n.size(); ++i)
         {
         std::cout<<du_n[i]<<std::endl;
         }
         
         for(std::size_t i=0; i<dp_n.size(); ++i)
         {
         std::cout<<dp_n[i]<<std::endl;
         }
        */
        
        
        std::cout << "Im here " <<std::endl;
        QGauss<2> quadrature_formula(degree*2);
        
        DoFHandler<2>::active_cell_iterator cell;
        auto vertices = triangulation.get_vertices();
        double *c1Node = new double [vertices.size()];
        double *c2Dist = new double [vertices.size()];
        double *LL = new double [vertices.size()];
        double *MM = new double [vertices.size()];
        
        FEValues<2> fe_values(fe,quadrature_formula,
                    update_values | update_quadrature_points |
                    update_JxW_values | update_gradients |
                        update_inverse_jacobians);
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        Vector<float> fe_degrees(triangulation.n_active_cells());
        
        unsigned int n_q_points  = fe_values.n_quadrature_points;
        
        std::cout << "----------v -------------  " <<std::endl;
        for (unsigned long int i = 0; i<vertices.size(); i++)
        {
            std::cout << "----------v -------------  " << i<<std::endl;
            
            // loop though patch elements get least square integral as H grid construction
            int nElement = v_to_e_list[i].size();
            
            std::vector<std::vector<double> > uEle;
            std::vector<std::vector<double> > vEle;
            std::vector<std::vector<double> > pEle;
            uEle.resize(nElement);
            vEle.resize(nElement);
            pEle.resize(nElement);
            
            std::vector<std::vector<double> > utEle;
            std::vector<std::vector<double> > vtEle;
            utEle.resize(nElement);
            vtEle.resize(nElement);
            
            std::vector<std::vector<double> > xEle;
            std::vector<std::vector<double> > yEle;
            xEle.resize(nElement);
            yEle.resize(nElement);
            
            std::vector<std::vector<double> > JxEle;
            JxEle.resize(nElement);
            
            std::vector<int> pH;
            std::vector<double> PatchArea;
            pH.resize(nElement);
            PatchArea.resize(nElement);
            
            std::cout << "resize finished " << i<<std::endl;
            for (int j = 0; j<nElement; j++)
             {
                 std::cout << "element loop" << j<<std::endl;
                 cell = v_to_e_list[i][j];
                 int cell_index = cell->active_cell_index();
                 //std::cout << "v_to_e_indices  "<<cell_index<<std::endl;
                 
                 //const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
                 //std::cout<<" dofs_per_cell "<<dofs_per_cell<<std::endl;
                 fe_values.reinit(cell);
                 
                 //const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
                 unsigned int n_q_points  = fe_values.n_quadrature_points;
                 //std::cout<<" n_q_points "<<n_q_points<<std::endl;
                 
                 std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
                 //std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
                 std::vector<double>       present_pressure_values(n_q_points);
                 //std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
                 std::vector<Tensor<1, 2> > present_vt_values(n_q_points);
                 
                 fe_values[velocities].get_function_values(n_solution,present_velocity_values);
                 //fe_values[velocities].get_function_gradients(n_solution, present_velocity_gradients);
                 fe_values[pressure].get_function_values(n_solution,present_pressure_values);
                 //fe_values[pressure].get_function_gradients(n_solution,present_pressure_gradients);
                 fe_values[velocities].get_function_values(n_solution_time_derivative,present_vt_values);
                 // get least square H grid constructions
                 pH[j] = fe[cell->active_fe_index()].degree;
                 //long double area = 0;
                 /*
                 if (i == 2)
                 {
                 std::cout << " j "<<j<<" v_to_e_indices "<<cell_index<<std::endl;
                 //std::cout<<" present_velocity_values "<<" "<<present_velocity_values<<std::endl;
                 //std::cout<<" present_pressure_values "<<" "<<present_pressure_values<<std::endl;
                 }
                  */
                 for (unsigned int q=0; q<n_q_points; ++q)
                 {
                     // get u,v,p at quadrature points
                     //std::cout << "quadrature loop " << q<<n_q_points<<std::endl;
                     uEle[j].push_back(present_velocity_values[q][0]);
                     vEle[j].push_back(present_velocity_values[q][1]);
                     
                     utEle[j].push_back(present_vt_values[q][0]);
                     vtEle[j].push_back(present_vt_values[q][1]);
                     
                     pEle[j].push_back(present_pressure_values[q]);
                     //std::cout<<"present_velocity_values "<<uEle[j][q]<<std::endl;
                     // get x y at quadrature points
                     const Point<2> &vPoint=fe_values.quadrature_point(q);
                     xEle[j].push_back(vPoint[0]);
                     yEle[j].push_back(vPoint[1]);
                     //std::cout << "vPoint  " << vPoint[0]<<" "<< vPoint[1]<<std::endl;
                     //double Jx = fe_values.JxW(q);
                     JxEle[j].push_back(fe_values.JxW(q));
                     //area += fe_values.JxW (q);
                     
                     /*
                     if (i == 2)
                     {
                         std::cout << " q "<<q<<" xEle "<<vPoint[0]<<" yEle" <<vPoint[1]<<"JxEle "<<fe_values.JxW(q)<<std::endl;
                         
                     std::cout<<" uEle "<<present_velocity_values[q]<<" utEle "<<present_vt_values[q] <<std::endl;
                
                    std::cout<<" pEle " <<present_pressure_values[q]<<std::endl;

                     }
                     */
                 }
                 std::cout << "end quadrature loop " <<n_q_points<<std::endl;
                 PatchArea[j] = cellArea[cell_index];
                 //std::cout << "cellArea  " << cellArea[cell_index]<<" "<< area<<std::endl;
                 
             } // end patch element loop
            int patchOrder = *max_element(pH.begin(), pH.end());
            int c2 = determinc2(uEle,vEle,patchOrder,PatchArea,nElement);
            c2Dist[i] = c2;
            /*
            for (int j = 0; j<nElement; j++)
            {
                std::cout<<"pH "<<pH[j]<<std::endl;
                std::cout<<"PatchArea "<<PatchArea[j]<<std::endl;
                for (unsigned int q=0; q<uEle[j].size(); ++q)
                {
                    std::cout<<"present_velocity_values "<<uEle[j][q]<<std::endl;
                    std::cout << "vPoint  " << xEle[j][q]<<" "<< yEle[j][q]<<std::endl;
                    std::cout<<"JxEle "<<JxEle[j][q]<<std::endl;
                }
            }
             */
            //int flag = vertex_value_u_map[i][0];
            //double value = vertex_value_u_map[i][1];
            //std::cout<<"flag "<<flag<<std::endl;
            
            Vector<double> LSu;
            Vector<double> LSv;
            Vector<double> LSp;
            Vector<double> LSut;
            Vector<double> LSvt;
            //std::cout << "vertex_value_u_map_value  " << vertex_value_u_map[i][0]<<" "<<vertex_value_u_map[i][1]<<" "<<vertex_value_u_map[i][2]<<std::endl;
            //std::cout << "vertex_value_v_map_value  " << vertex_value_v_map[i][0]<<" "<<vertex_value_v_map[i][1]<<" "<<vertex_value_v_map[i][2]<<std::endl;
            //std::cout << "vertex_value_p_map_value  " << vertex_value_p_map[i][0]<<" "<<vertex_value_p_map[i][1]<<" "<<vertex_value_p_map[i][2]<<std::endl;
            
            std::vector<double> xmaxEle;
            std::vector<double> xminEle;
            std::vector<double> ymaxEle;
            std::vector<double> yminEle;
            xmaxEle.resize(nElement);
            xminEle.resize(nElement);
            ymaxEle.resize(nElement);
            yminEle.resize(nElement);
            
            std::cout << "here 1 " <<n_q_points<<std::endl;
            for (int j = 0; j<nElement; j++)
            {
                xmaxEle[j] = *max_element(xEle[j].begin(), xEle[j].end());
                xminEle[j] = *min_element(xEle[j].begin(), xEle[j].end());
                ymaxEle[j] = *max_element(yEle[j].begin(), yEle[j].end());
                yminEle[j] = *min_element(yEle[j].begin(), yEle[j].end());
                
            }
            double xmax = *max_element(xmaxEle.begin(), xmaxEle.end());
            double xmin = *min_element(xminEle.begin(), xminEle.end());
            double ymax = *max_element(ymaxEle.begin(), ymaxEle.end());
            double ymin = *min_element(yminEle.begin(), yminEle.end());
            std::cout << "here 1.1 " <<n_q_points<<std::endl;
            int sizeA = (patchOrder+1)*(patchOrder+2)/2;
            FullMatrix<double> lhs;
            lhs.reinit(sizeA, sizeA);
            lhs = 0;
            int count = 0;
            for (int k1=0; k1<sizeA; k1++)
            {
                for (int k2=0; k2<sizeA; k2++)
                {
                    lhs[k1][k2] = lhsG[i][count];
                    //std::cout << "lhs  "<<lhsG[i][count]<<std::endl;
                    count++;
                }
                //std::cout << "lhs  "<<lhsG[i][k]<<std::endl;
            }
            std::cout << "here 1.2 " <<n_q_points<<std::endl;
            count = 0;
            double **rhsxy = new double *[sizeA];
            for(int i1 = 0; i1 < sizeA; ++i1)
            {
                rhsxy[i1] = new double[n_q_points*nElement];
            }
            for (int k1 = 0; k1<sizeA;k1++)
            {
                for (unsigned int k2 = 0; k2<n_q_points*nElement;k2++)
                {
                    rhsxy[k1][k2] = rhsG[i][count];
                    //std::cout << "rhs  "<<rhsxy[k1][k2]<<std::endl;
                    count++;
                }
                //std::cout << "rhs  "<<rhsG[i][lv]<<std::endl;
            }
            std::cout << "here 1.3 " <<n_q_points<<std::endl;
            Point<2> &vcord_u = vertices[vertex_value_u_map[i][1]];
            std::cout << "here 1.4 " <<n_q_points<<std::endl;
            getLSint(xmax,xmin,ymax,ymin,nElement,uEle,xEle,yEle,vertex_value_u_map[i][0],vcord_u,vertex_value_u_map[i][2],patchOrder,JxEle,lhs,rhsxy,LSu);
            std::cout << "here 1.5 " <<n_q_points<<std::endl;
            Point<2> &vcord_v = vertices[vertex_value_v_map[i][1]];
            getLSint(xmax,xmin,ymax,ymin,nElement,vEle,xEle,yEle,vertex_value_v_map[i][0],vcord_v,vertex_value_u_map[i][2],patchOrder,JxEle,lhs,rhsxy,LSv);
            Point<2> &vcord_p = vertices[vertex_value_p_map[i][1]];
            getLSint(xmax,xmin,ymax,ymin,nElement,pEle,xEle,yEle,vertex_value_p_map[i][0],vcord_p,vertex_value_u_map[i][2],patchOrder,JxEle,lhs,rhsxy,LSp);
            
            getLSint(xmax,xmin,ymax,ymin,nElement,utEle,xEle,yEle,vertex_value_u_map[i][0],vcord_u,vertex_value_u_map[i][2],patchOrder,JxEle,lhs,rhsxy,LSut);
            getLSint(xmax,xmin,ymax,ymin,nElement,vtEle,xEle,yEle,vertex_value_v_map[i][0],vcord_v,vertex_value_u_map[i][2],patchOrder,JxEle,lhs,rhsxy,LSvt);
            //std::cout << "here 1.4 " <<n_q_points<<std::endl;
            /*
            if (i == 2)
            {
                std::cout<<" xmax "<<xmax<<" xmin "<<xmin<<" ymax "<<ymax<<" ymin "<<ymin<<std::endl;
            std::cout <<" bc flag "<<vertex_value_u_map[i][0]<<" fitvID"<<vertex_value_u_map[i][1]<<std::endl;
            std::cout<<" vcord_u "<<vcord_u[0]<<" "<<vcord_u[1]<<" bc value "<<vertex_value_u_map[i][2]<<std::endl;
            //std::cout<<" present_pressure_values "<<" "<<present_pressure_values<<std::endl;
            
            std::cout << "LSu  "<<LSu[0]<<" "<<LSu[1]<<" "<<LSu[2]<<std::endl;
            std::cout << "LSv  "<<LSv[0]<<" "<<LSv[1]<<" "<<LSv[2]<<std::endl;
            std::cout << "LSp  "<<LSp[0]<<" "<<LSp[1]<<" "<<LSp[2]<<std::endl;
            std::cout << "LSut  "<<LSut[0]<<" "<<LSut[1]<<" "<<LSut[2]<<std::endl;
            std::cout << "LSvt  "<<LSvt[0]<<" "<<LSvt[1]<<" "<<LSvt[2]<<std::endl;
            }
            */
            // loop though patch elements again get LH_conv, MH_c1 and MH_c2
            double Lh_conv = 0; double Mh_c1 = 0; double Mh_c2 = 0;
            double LH_conv = 0; double MH_c1 = 0; double MH_c2 = 0;
            
            double WHx[4];
            double pArea = 0.0;
            pArea = accumulate(PatchArea.begin(), PatchArea.end(), pArea);
            double HH = sqrt(pArea);
            
            count = 0;
            //std::cout << "pArea  "<<pArea<<" "<<HH<<std::endl;
            double *xyuH = new double[sizeA];
            double *xyduHx = new double[sizeA];
            double *xyduHy = new double[sizeA];
            
            std::cout << "here 2 " <<n_q_points<<std::endl;
            for (unsigned long int j = 0; j<v_to_e_list[i].size(); j++)
            {
                cell = v_to_e_list[i][j];
                //int cell_index = cell->active_cell_index();
                //std::cout << "v_to_e_indices  "<<cell_index<<std::endl;
                
                //const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
                //std::cout<<" dofs_per_cell "<<dofs_per_cell<<std::endl;
                fe_values.reinit(cell);
                
                //const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
                unsigned int n_q_points  = fe_values.n_quadrature_points;
                //std::cout<<" n_q_points "<<n_q_points<<std::endl;
                
                std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
                std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
                //std::vector<double>       present_pressure_values(n_q_points);
                std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
                std::vector<Tensor<1, 2> > present_vt_values(n_q_points);
                
                //fe_values[velocities].get_function_values(n_solution,present_velocity_values);
                fe_values[velocities].get_function_gradients(n_solution, present_velocity_gradients);
                //fe_values[pressure].get_function_values(n_solution,present_pressure_values);
                fe_values[pressure].get_function_gradients(n_solution,
                                                           present_pressure_gradients);
                fe_values[velocities].get_function_values(n_solution_time_derivative,present_vt_values);
                
                std::vector<Tensor<1, 2> > rm(n_q_points);
                double hh = sqrt(PatchArea[j]);
                std::cout << "hh  " << hh<<std::endl;
                for (unsigned int q=0; q<n_q_points; ++q)
                {
                    present_velocity_values[q][0] = uEle[j][q];
                    present_velocity_values[q][1] = vEle[j][q];
                    
                    present_vt_values[q][0] = utEle[j][q];
                    present_vt_values[q][1] = vtEle[j][q];
                    // get u,v,p at quadrature points
                    //std::cout<<"present_velocity_values "<< " "<<uEle[j][q]<< " "<<vEle[j][q]<<" "<<pEle[j][q]<<std::endl;
                    // get x y at quadrature points
                    const Point<2> &vPoint = fe_values.quadrature_point(q);
                    //std::cout << "vPoint  " << vPoint[0]<<" "<< vPoint[1]<<std::endl;
                    
                    double Jx = fe_values.JxW(q);
                    
                    for (int np=0; np<sizeA;np++)
                    {
                        xyuH[np] = xyuHG[i][count] ;
                        xyduHx[np] = xyduHxG[i][count] ;
                        xyduHy[np] = xyduHyG[i][count] ;
                        //std::cout <<"i "<<i<<" np "<<np<< " xyuHG  " << xyuHG[i][count]<<std::endl;
                        count = count+1;
                    }
                    
                    double uH = 0; double duH[2] = {};
                    calcuH(vPoint,LSu,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy,uH,duH);
                    double vH = 0; double dvH[2] = {};
                    calcuH(vPoint,LSv,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy,vH,dvH);
                    double pressureH = 0; double dpH[2] = {};
                    calcuH(vPoint,LSp,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy,pressureH,dpH);
                   
                    double uHt = 0; double duHt[2] = {};
                    calcuH(vPoint,LSut,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy,uHt,duHt);
                    double vHt = 0; double dvHt[2] = {};
                    calcuH(vPoint,LSvt,patchOrder,xmin,xmax,ymin,ymax,xyuH,xyduHx,xyduHy,vHt,dvHt);
                    
                    //std::cout << "duH  " << duH[0]<<" "<< duH[1]<<std::endl;
                    //std::cout << "dvH  " << dvH[0]<<" "<< dvH[1]<<std::endl;
                    
                    WHx[0] = duH[0];
                    WHx[1] = duH[1];
                    WHx[2] = dvH[0];
                    WHx[3] = dvH[1];
                    
                    
                    double rumH = uHt+uH*duH[0]+vH*duH[1]+dpH[0];
                    double rvmH = vHt+uH*dvH[0]+vH*dvH[1]+dpH[1];
                    
                    //std::cout << "rumH  " << rumH<<" rvmH "<< rvmH<<std::endl;
                    
                    LH_conv += (WHx[0]*uH*uH+WHx[1]*uH*vH+WHx[2]*vH*uH+WHx[3]*vH*vH)*Jx
                                +(uH*uHt+uH*uHt+vH*vHt+vH*vHt)*Jx;
                    //std::cout << "LH_conv  " << LH_conv<<std::endl;
                    
                    MH_c1 += pow(HH,c2)*(WHx[0]*uH*rumH+WHx[1]*uH*rumH+WHx[2]*vH*rvmH+WHx[3]*vH*rvmH)*Jx;
                    MH_c2 += pow(HH,c2)*(WHx[0]*uH*rumH+WHx[1]*vH*rumH+WHx[2]*uH*rvmH+WHx[3]*vH*rvmH)*Jx;
                    //std::cout << "MH_c1  " << MH_c1<<std::endl;
                    //std::cout << "MH_c2  " << MH_c2<<std::endl;
                    
                    rm[q] = present_vt_values[q]+present_velocity_gradients[q]*present_velocity_values[q]+present_pressure_gradients[q];//-sourceTerm[q];
                    
                    
                    double uele = uEle[j][q]; double vele = vEle[j][q];
                    double utele = utEle[j][q]; double vtele = vtEle[j][q];
                    double ux = present_velocity_gradients[q][0][0];
                    double uy = present_velocity_gradients[q][0][1];
                    double vx = present_velocity_gradients[q][1][0];
                    double vy = present_velocity_gradients[q][1][1];
                    Lh_conv +=(uH*uele*ux+uH*vele*uy+vH*uele*vx+vH*vele*vy)*Jx+(uH*utele+uH*utele+vH*vtele+vH*vtele)*Jx;
                    
                    //Lh_conv += (WHx[0]*uele*uele+WHx[1]*uele*vele+WHx[2]*vele*uele+WHx[3]*vele*vele)*Jx+(uH*utele+uH*utele+vH*vtele+vH*vtele)*Jx;
                    //std::cout << "Lh_conv  " << Lh_conv<<std::endl;
                    
                    double rum = rm[q][0];
                    double rvm = rm[q][1];
                    
                    //uH*ux*rum+uH*uy*rum+vH*vx*rvm+vH*vy*rvm
                    Mh_c1 += pow(hh,c2)*(uH*ux*rum+uH*uy*rum+vH*vx*rvm+vH*vy*rvm)*Jx;
                    Mh_c2 += pow(hh,c2)*(duH[0]*uele*rum+duH[1]*vele*rum+dvH[0]*uele*rvm+dvH[1]*vele*rvm)*Jx;
                    //std::cout << "Mh_c1  " << Mh_c1<<std::endl;
                    //std::cout << "Mh_c2  " << Mh_c2<<std::endl;
                    
                }// end quadrature
                std::cout << "here 3 " <<n_q_points<<std::endl;
                
            } //end patch iteration
        
             LL[i] = -(Lh_conv-LH_conv);
             MM[i] = -(Mh_c1-MH_c1)-(Mh_c2-MH_c2);
            
            /*
            if (i == 2)
            {
                std::cout<<" Lh_conv "<<Lh_conv<<" LH_conv "<<LH_conv<<std::endl;
                std::cout<<" Mh_c1 "<<Mh_c1<<" MH_c1 "<<MH_c1<<std::endl;
                std::cout<<" Mh_c2 "<<Mh_c2<<" MH_c2 "<<MH_c2<<std::endl;
                std::cout<<" LL "<<LL[i]<<" MM "<<MM[i]<<std::endl;
            }
            */
            } // end vertex iteration
        
            //std::cout << "LL  " << LL<<" Lh_conv "<<Lh_conv<<" LH_conv "<<LH_conv<<std::endl;
            //std::cout << "MM  " << MM<<std::endl;
            std::cout << "end vertex iteration  " << MM<<std::endl;
            // add volume average
            double tmp; double aref = 1; int pH;
            for (unsigned long int node = 0; node<vertices.size(); node++)
            {   int nElement = v_to_e_list[node].size();
                //std::cout << "node  " << node <<std::endl;
                
                double Lg = 0; double Mg = 0;
                for (int j = 0; j<nElement; j++)
                {
                    cell = v_to_e_list[node][j];
                    fe_values.reinit(cell);
                    pH = fe[cell->active_fe_index()].degree;
                    
                    // get vertex for patch cell
                    for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
                    {   int v_index = cell->vertex_index(i);
                        Lg = Lg+std::abs(LL[v_index]);
                        Mg = Mg+std::abs(MM[v_index]);
                        //std::cout << "v_index  " << v_index<<std::endl;
                    }
                
                } // end patch element loop
                
                if (Mg<1e-10)
                {
                    tmp = 0;
                }
                else
                {
                    tmp = std::abs(Lg/Mg);
                }
                
                if (c2Dist[node]>1)
                {
                    if (tmp > 10/viscosity/pH/pH)
                        tmp = 10/viscosity/pH/pH;
                    
                         if (tmp < 0.01/viscosity/pH/pH)
                            tmp = 0.01/viscosity/pH/pH;
                }
                else
                {
                    if (tmp > 10/aref/pH)
                            tmp = 10/aref/pH;
                    
                    if (tmp < 0.01/aref/pH)
                            tmp = 0.01/aref/pH;
                }
                
                c1Node[node] = tmp;

                /*
                if (node == 2)
                {
                    std::cout<<" c1Node[node] "<<c1Node[node]<<std::endl;
                }
                */
            } // end node loop
            
            // add boundary condition
            for (unsigned long int i = 0; i<IBC.size(); i++)
            {
                c1Node[IBC[i]] = 0;
                //std::cout<<"IBC "<<IBC[i]<<std::endl;
            }
            
            // compute element-wise tau, loop through element
            DoFHandler<2>::active_cell_iterator celliter = dof_handler.begin_active();
            DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
            int cellnum = 0;
            for (; celliter != endc; ++celliter)
            {   double hh = sqrt(cellArea[cellnum]);
                //std::cout<<"taumele_g cellnum "<<cellnum<<std::endl;
                // get vertex for patch cell
                for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
                {   int v_index = celliter->vertex_index(i);
                    taumele_g[cellnum][i] = c1Node[v_index]*pow(hh,c2Dist[v_index]);
                    //std::cout<<taumele_g[cellnum][i]<<std::endl;
                }
                cellnum++;
            } // end loop for taumele_g

    } // end function
    
    
    int dynamic::determinc2(std::vector<std::vector<double> > uEle,
                             std::vector<std::vector<double> > vEle,
                             int patchOrder,
                             std::vector<double> PatchArea,
                             int nElement)
    {
        int count = 0;
        std::vector<double> velocity_scalar;
        double Parea = 0;
        //int patchOrder = *max_element(pH.begin(), pH.end());
        for (int j = 0; j<nElement; j++)
        {
            //std::cout<<"pH "<<pH[j]<<std::endl;
            //std::cout<<"PatchArea "<<PatchArea[j]<<std::endl;
            Parea = Parea+PatchArea[j];
            for (unsigned int q=0; q<uEle[j].size(); ++q)
            {
                velocity_scalar.push_back(sqrt(uEle[j][q]*uEle[j][q]+vEle[j][q]*vEle[j][q]));
                count++;
                //std::cout<<"velocity_scalar "<<sqrt(uEle[j][q]*uEle[j][q]+vEle[j][q]*vEle[j][q])<<std::endl;
            }
        }
        
        //std::cout<<"pArea "<<pArea<<" "<<Parea<<std::endl;
        
        double max_velocity = *max_element(velocity_scalar.begin(), velocity_scalar.end());
        
        double constant = 4/(sqrt(3));
        double PeH = sqrt(constant*Parea)*max_velocity/(2*(viscosity)*patchOrder);
        
        //std::cout<<"PeH "<<PeH<<std::endl;
        //std::cout<<"pH "<<patchOrder<<std::endl;
        //std::cout<<"max_velocity "<<max_velocity<<std::endl;
        int c2;
        if (PeH > 3)
             c2 = 1;
        else
             c2 = 2;
        
        //std::cout<<"c2 "<<c2<<std::endl;
        return c2;
    }
        
        void dynamic::getLSintlhs(double xmax,double xmin,double ymax,double ymin,
                               int nElement,
                               std::vector<std::vector<double> > xEle,
                               std::vector<std::vector<double> > yEle,
                               int p,
                               std::vector<std::vector<double> > JxEle,
                               double *lhs1D, double *rhs1D)
        {
            double **x;
            double **y;
            
            x = new double *[nElement];
            y = new double *[nElement];
            int q_points = xEle[0].size();
            for( int i = 0; i <nElement; i++)
            {
                x[i] = new double[q_points];
                y[i] = new double[q_points];
            }
            
            mapxy(nElement,xEle,yEle, xmax,xmin,ymax,ymin,x,y);
            
            int sizeA = (p+1)*(p+2)/2;
            
            FullMatrix<double> lhs;
            FullMatrix<double> rhs;
            
            lhs.reinit(sizeA, sizeA);
            lhs = 0;
            
            rhs.reinit(sizeA,q_points*nElement);
            rhs = 0;
            
            int r=0;
            int c;
            double sumxyn;
            for (int k1 = 0; k1<=p; k1++) //k1 = 0:p
            {
                for (int m1 = 0; m1<=k1; m1++) //m1 = 0:k1
                {
                    r=r+1;
                    c=0;
                    for (int k2 = 0; k2<=p; k2++) //k2 = 0:p
                    {
                        for (int m2 = 0; m2<=k2; m2++) //m2 = 0:k2
                        {
                            c=c+1;
                            sumxyn = 0;
                            for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                            {
                                for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                                {
                                    sumxyn = sumxyn+pow(x[ele][l],(k1-m1))*pow(y[ele][l],(m1))*pow(x[ele][l],(k2-m2))*pow(y[ele][l],(m2))*JxEle[ele][l];
                                }
                                
                            } // l
                            lhs(r-1,c-1) = sumxyn;
                        }// m2
                    }// k2
                }//m1
            }//k1
            
            
            int i;
            if (p==1)
            {
                i = 0;
                for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                    {
                    rhs[0][i] = JxEle[ele][l];
                    rhs[1][i] = x[ele][l]*JxEle[ele][l];
                    rhs[2][i] = y[ele][l]*JxEle[ele][l];
                    i++;
                    }

                }//q_points
            }
            else if (p==2)
            {
                i = 0;
                for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                    {
                        rhs[0][i] = JxEle[ele][l];
                        rhs[1][i] = x[ele][l]*JxEle[ele][l];
                        rhs[2][i] = y[ele][l]*JxEle[ele][l];
                        rhs[3][i] = x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[4][i] = x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[5][i] = y[ele][l]*y[ele][l]*JxEle[ele][l];
                        i++;
                    }
                    
                }//q_points
            }
            else if (p==3)
            {
                i = 0;
                for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                    {
                        rhs[0][i] = JxEle[ele][l];
                        rhs[1][i] = x[ele][l]*JxEle[ele][l];
                        rhs[2][i] = y[ele][l]*JxEle[ele][l];
                        rhs[3][i] = x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[4][i] = x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[5][i] = y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[6][i] = x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[7][i] = x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[8][i] = x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[9][i] = y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        i++;
                    }
                    
                }//q_points
            }
            else if (p==4)
            {
                i = 0;
                for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                    {
                        rhs[0][i] = JxEle[ele][l];
                        rhs[1][i] = x[ele][l]*JxEle[ele][l];
                        rhs[2][i] = y[ele][l]*JxEle[ele][l];
                        rhs[3][i] = x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[4][i] = x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[5][i] = y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[6][i] = x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[7][i] = x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[8][i] = x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[9][i] = y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        
                        rhs[10][i] = x[ele][l]*x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[11][i] = x[ele][l]*x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[12][i] = x[ele][l]*x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[13][i] = x[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[14][i] = y[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        i++;
                    }
                    
                }//q_points
            }//if else
            else if (p==5)
            {
                i = 0;
                for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                    {
                        rhs[0][i] = JxEle[ele][l];
                        rhs[1][i] = x[ele][l]*JxEle[ele][l];
                        rhs[2][i] = y[ele][l]*JxEle[ele][l];
                        rhs[3][i] = x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[4][i] = x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[5][i] = y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[6][i] = x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[7][i] = x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[8][i] = x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[9][i] = y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        
                        rhs[10][i] = x[ele][l]*x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[11][i] = x[ele][l]*x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[12][i] = x[ele][l]*x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[13][i] = x[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[14][i] = y[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        
                        rhs[15][i] = x[ele][l]*x[ele][l]*x[ele][l]*x[ele][l]*x[ele][l]*JxEle[ele][l];
                        rhs[16][i] = x[ele][l]*x[ele][l]*x[ele][l]*x[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[17][i] = x[ele][l]*x[ele][l]*x[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[18][i] = x[ele][l]*x[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[19][i] = x[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        rhs[20][i] = y[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*y[ele][l]*JxEle[ele][l];
                        i++;
                    }
                    
                }//q_points
            }
                else
                {
                    std::cout<<"not support p>5 "<<std::endl;
                }
            
            int k = 0;
            for (int i = 0; i<sizeA;i++)
            {
                for(int j = 0; j<sizeA; j++)
                {
                    lhs1D[k] = lhs[i][j];
                    //std::cout<<"lhs "<<lhs[i][j]<<std::endl;
                    k++;
                }
            }
            
            k = 0;
            for (int i = 0; i<sizeA;i++)
            {
                for(int j = 0; j<q_points*nElement; j++)
                {
                    rhs1D[k] = rhs[i][j];
                    //std::cout<<"rhs "<<rhs[i][j]<<std::endl;
                    k++;
                }
            }
        }

    void dynamic::calcuHxycoeff(Point<2> vpoint,
                         int p,double xmin,double xmax,double ymin,double ymax,
                         double *xyuH,double *xyduHx,double *xyduHy)
    {
        double x = vpoint[0];
        double y = vpoint[1];
        
        double xi; double eta;
        // mapping to a box
        mapxy(x,y,xmax,xmin,ymax,ymin,xi,eta);
        
        // uH = a0+a1x+a2y+a3x^2+a4xy+a5y^2+...
        if (p==1)
        {
            xyuH[0] = 1;
            xyduHx[0] = 0;
            xyduHy[0] = 0;
                
            xyuH[1] = xi;
            xyduHx[1] = 1;
            xyduHy[1] = 0;
            
            xyuH[2] = eta;
            xyduHx[2] = 0;
            xyduHy[2] = 1;
            
        }
        else if (p==2)
        {
            xyuH[0] = 1;
            xyduHx[0] = 0;
            xyduHy[0] = 0;
            
            xyuH[1] = xi;
            xyduHx[1] = 1;
            xyduHy[1] = 0;
            
            xyuH[2] = eta;
            xyduHx[2] = 0;
            xyduHy[2] = 1;
            
            xyuH[3] = xi*xi;
            xyduHx[3] = 2*xi;
            xyduHy[3] = 0;
            
            xyuH[4] = xi*eta;
            xyduHx[4] = eta;
            xyduHy[4] = xi;
            
            xyuH[5] = eta*eta;
            xyduHx[5] = 0;
            xyduHy[5] = 2*eta;
            
        }
        else if (p==3)
        {
            xyuH[0] = 1;
            xyduHx[0] = 0;
            xyduHy[0] = 0;
            
            xyuH[1] = xi;
            xyduHx[1] = 1;
            xyduHy[1] = 0;
            
            xyuH[2] = eta;
            xyduHx[2] = 0;
            xyduHy[2] = 1;
            
            xyuH[3] = xi*xi;
            xyduHx[3] = 2*xi;
            xyduHy[3] = 0;
            
            xyuH[4] = xi*eta;
            xyduHx[4] = eta;
            xyduHy[4] = xi;
            
            xyuH[5] = eta*eta;
            xyduHx[5] = 0;
            xyduHy[5] = 2*eta;
            
            //x.^3,x.^2.*y,x.*y.^2,y.^3
            //3*x.^2,2*x.*y,y.^2,0
            //0,x.^2,x.*2*y,3*y.^2
            xyuH[6] = xi*xi*xi;
            xyduHx[6] = 3*xi*xi;
            xyduHy[6] = 0;
            
            xyuH[7] = xi*xi*eta;
            xyduHx[7] = 2*xi*eta;
            xyduHy[7] = xi*xi;
            
            xyuH[8] = xi*eta*eta;
            xyduHx[8] = eta*eta;
            xyduHy[8] = 2*xi*eta;
            
            xyuH[9] = eta*eta*eta;
            xyduHx[9] = 0;
            xyduHy[9] = 3*eta*eta;
        }
        else if (p==4)
        {
            xyuH[0] = 1;
            xyduHx[0] = 0;
            xyduHy[0] = 0;
            
            xyuH[1] = xi;
            xyduHx[1] = 1;
            xyduHy[1] = 0;
            
            xyuH[2] = eta;
            xyduHx[2] = 0;
            xyduHy[2] = 1;
            
            xyuH[3] = xi*xi;
            xyduHx[3] = 2*xi;
            xyduHy[3] = 0;
            
            xyuH[4] = xi*eta;
            xyduHx[4] = eta;
            xyduHy[4] = xi;
            
            xyuH[5] = eta*eta;
            xyduHx[5] = 0;
            xyduHy[5] = 2*eta;
            
            xyuH[6] = xi*xi*xi;
            xyduHx[6] = 3*xi*xi;
            xyduHy[6] = 0;
            
            xyuH[7] = xi*xi*eta;
            xyduHx[7] = 2*xi*eta;
            xyduHy[7] = xi*xi;
            
            xyuH[8] = xi*eta*eta;
            xyduHx[8] = eta*eta;
            xyduHy[8] = 2*xi*eta;
            
            xyuH[9] = eta*eta*eta;
            xyduHx[9] = 0;
            xyduHy[9] = 3*eta*eta;
            
            //x.^4,x.^3.*y,x.^2*y.^2,x.*y.^3,y.^4
            //4*x.^3,3*x.^2*y,2*x*y.^2,y.^3,0
            //0,x.^3,x.^2*2*y,x.*3*y.^2,4*y.^3
            xyuH[10] = xi*xi*xi*xi;
            xyduHx[10] = 4*xi*xi*xi;
            xyduHy[10] = 0;
            
            xyuH[11] = xi*xi*xi*eta;
            xyduHx[11] = 3*xi*xi*eta;
            xyduHy[11] = xi*xi*xi;
            
            xyuH[12] = xi*xi*eta*eta;
            xyduHx[12] = 2*xi*eta*eta;
            xyduHy[12] = 2*xi*xi*eta;
            
            xyuH[13] = xi*eta*eta*eta;
            xyduHx[13] = eta*eta*eta;
            xyduHy[13] = 3*eta*eta*xi;
            
            xyuH[14] = eta*eta*eta*eta;
            xyduHx[14] = 0;
            xyduHy[14] = 4*eta*eta*eta;
            
        }//if else
        else if (p==5)
        {
            xyuH[0] = 1;
            xyduHx[0] = 0;
            xyduHy[0] = 0;
            
            xyuH[1] = xi;
            xyduHx[1] = 1;
            xyduHy[1] = 0;
            
            xyuH[2] = eta;
            xyduHx[2] = 0;
            xyduHy[2] = 1;
            
            xyuH[3] = xi*xi;
            xyduHx[3] = 2*xi;
            xyduHy[3] = 0;
            
            xyuH[4] = xi*eta;
            xyduHx[4] = eta;
            xyduHy[4] = xi;
            
            xyuH[5] = eta*eta;
            xyduHx[5] = 0;
            xyduHy[5] = 2*eta;
            
            xyuH[6] = xi*xi*xi;
            xyduHx[6] = 3*xi*xi;
            xyduHy[6] = 0;
            
            xyuH[7] = xi*xi*eta;
            xyduHx[7] = 2*xi*eta;
            xyduHy[7] = xi*xi;
            
            xyuH[8] = xi*eta*eta;
            xyduHx[8] = eta*eta;
            xyduHy[8] = 2*xi*eta;
            
            xyuH[9] = eta*eta*eta;
            xyduHx[9] = 0;
            xyduHy[9] = 3*eta*eta;
            
            xyuH[10] = xi*xi*xi*xi;
            xyduHx[10] = 4*xi*xi*xi;
            xyduHy[10] = 0;
            
            xyuH[11] = xi*xi*xi*eta;
            xyduHx[11] = 3*xi*xi*eta;
            xyduHy[11] = xi*xi*xi;
            
            xyuH[12] = xi*xi*eta*eta;
            xyduHx[12] = 2*xi*eta*eta;
            xyduHy[12] = 2*xi*xi*eta;
            
            xyuH[13] = xi*eta*eta*eta;
            xyduHx[13] = eta*eta*eta;
            xyduHy[13] = 3*eta*eta*xi;
            
            xyuH[14] = eta*eta*eta*eta;
            xyduHx[14] = 0;
            xyduHy[14] = 4*eta*eta*eta;
            
            xyuH[15] = xi*xi*xi*xi*xi;
            xyduHx[15] = 5*xi*xi*xi*xi;
            xyduHy[15] = 0;
            
            xyuH[16] = xi*xi*xi*xi*eta;
            xyduHx[16] = 4*xi*xi*xi*eta;
            xyduHy[16] = xi*xi*xi*xi;
            
            xyuH[17] = xi*xi*xi*eta*eta;
            xyduHx[17] = 3*xi*xi*eta*eta;
            xyduHy[17] = 2*xi*xi*xi*eta;
            
            xyuH[18] = xi*xi*eta*eta*eta;
            xyduHx[18] = 2*xi*eta*eta*eta;
            xyduHy[18] = 3*eta*eta*xi*xi;
            
            xyuH[19] = xi*eta*eta*eta*eta;
            xyduHx[19] = eta*eta*eta*eta;
            xyduHy[19] = 4*eta*eta*eta*xi;
            
            xyuH[20] = eta*eta*eta*eta*eta;
            xyduHx[20] = 0;
            xyduHy[20] = 5*eta*eta*eta*eta;
            
        }//if else
        
        else
        {
            std::cout<<"not support p>5 "<<std::endl;
        }
        //std::cout<<"uH "<<uH<<std::endl;
        //std::cout<<"duH "<<duH[0]<<" "<<duH[1]<<std::endl;
        
    }
    
    
    void dynamic::getLSint(double xmax,double xmin,double ymax,double ymin,
                           int nElement,std::vector<std::vector<double> > uEle,
                  std::vector<std::vector<double> > xEle,
                  std::vector<std::vector<double> > yEle,
                  bool flag,Point<2> & bcvertex,double bcvalue,
                  int p,
                  std::vector<std::vector<double> > JxEle,
                  FullMatrix<double> lhs, double **rhsxy,
                  Vector<double> &LS)
    {
        double **x;
        double **y;
        
        x = new double *[nElement];
        y = new double *[nElement];
        int q_points = xEle[0].size();
        for( int i = 0; i <nElement; i++)
        {
            x[i] = new double[q_points];
            y[i] = new double[q_points];
        }
        
        mapxy(nElement,xEle,yEle, xmax,xmin,ymax,ymin,x,y);
        
        int sizeA = (p+1)*(p+2)/2;
        
        //FullMatrix<double> lhs;
        //Vector<double>     rhs;
        
        //lhs.reinit(sizeA, sizeA);
        //lhs = 0;
        
        //rhs.reinit(sizeA);
        //rhs = 0;
        
        /*
        int r=0;
        int c;
        double sumxyn;
        for (int k1 = 0; k1<=p; k1++) //k1 = 0:p
        {
            for (int m1 = 0; m1<=k1; m1++) //m1 = 0:k1
            {
                r=r+1;
                c=0;
                for (int k2 = 0; k2<=p; k2++) //k2 = 0:p
                {
                    for (int m2 = 0; m2<=k2; m2++) //m2 = 0:k2
                    {
                        c=c+1;
                        sumxyn = 0;
                        for (int l = 0; l<q_points; l++) //l = 1:size(x,1)
                        {
                            for (int ele = 0; ele<nElement; ele++) //ele = 1:nEle
                            {
                                sumxyn = sumxyn+pow(x[ele][l],(k1-m1))*pow(y[ele][l],(m1))*pow(x[ele][l],(k2-m2))*pow(y[ele][l],(m2))*JxEle[ele][l];
                            }
        
                        }
                        lhs(r-1,c-1) = sumxyn;
                    }
                }
            }
        }
        */
        /*
        for (int k1=0; k1<sizeA; k1++)
        {
            for (int k2=0; k2<sizeA; k2++)
            {
                std::cout << "lhs  "<<lhs[k1][k2]<<std::endl;
                //std::cout << "lhstest  "<<lhstest[k1][k2]<<std::endl;
            }
            //std::cout << "lhs  "<<lhsG[i][k]<<std::endl;
        }
         */
        
        /*
        int c=0;
        double sumeleu;
        for (int k1 = 0; k1<=p; k1++)//k1 = 0:p
        {
            for (int m1 = 0; m1<=k1; m1++)//m1 = 0:k1
            {
                c=c+1;
                sumeleu = 0;
                for (int l = 0; l<q_points; l++)//l = 1:size(x,1)
                {
                    for (int ele = 0; ele<nElement; ele++)//ele = 1:nEle
                        {
                        sumeleu = sumeleu+uEle[ele][l]*pow(x[ele][l],(k1-m1))*pow(y[ele][l],(m1))*JxEle[ele][l];
                        }
                }
                rhs(c-1) = sumeleu;
            }
        }
        */
        Vector<double>     rhs;
        rhs.reinit(sizeA);
        rhs = 0;
        double sumeleu;
        std::cout << "here 1.4.1 " <<std::endl;
        for(int i1 = 0; i1<sizeA; i1++)
        {
            sumeleu = 0;
            int count = 0;
            for (int l = 0; l<q_points; l++)//l = 1:size(x,1)
            {
                for (int ele = 0; ele<nElement; ele++)//ele = 1:nEle
                {
                    sumeleu = sumeleu+uEle[ele][l]*rhsxy[i1][count];
                    count++;
                }
            }
            rhs[i1] = sumeleu;
        }
        std::cout << "here 1.4.2 " <<std::endl;
        
        for(int i1 = 0; i1<sizeA; i1++)
        {
            std::cout<<"rhs "<<rhs[i1]<<std::endl;
        }
        
        //std::cout<<"lhssize "<<lhs[0].size()<<std::endl;
        //for(int i1 = 0; i1<sizeA; i1++)
        //{
          //  for(int i2 = 0; i2<lhs[i1].size(); i1++)
            //{
              //  std::cout<<"lhs "<<lhs[i1][i2]<<std::endl;
            //}
        //}
        
        LS.reinit(sizeA);
        LS = 0;
        std::cout << "here 1.4.3 " <<sizeA<<std::endl;
        matrixSolve(sizeA,lhs,rhs,LS);
        std::cout << "here 1.4.3 " <<std::endl;
        if (flag)
        {
            // exact fit bcNode
            double atest = 0;
            int c=0;
            double xbcnode = bcvertex[0];
            double ybcnode = bcvertex[1];

            
            double xbc;
            double ybc;
            // need mapping
            mapxy(xbcnode,ybcnode,xmax,xmin,ymax,ymin,xbc,ybc);
            
            for (int k1 = 0; k1<=p; k1++) //k1 = 0:p
            {
                for (int m1 = 0; m1<=k1; m1++) //m1 = 0:k1
                {
                    c=c+1;
                    atest = atest-LS(c-1)*pow(xbc,(k1-m1))*pow(ybc,(m1));
                }
            }
        
            double a = atest+bcvalue+LS(0);
            LS(0) = a;
            
        }
        std::cout << "here 1.4.4 " <<std::endl;
        
    }
    
  
    void dynamic::calcuH(Point<2> vpoint, Vector<double> LS,
                         int p,double xmin,double xmax,double ymin,double ymax,
                         double *xyuH,double *xyduHx,double *xyduHy,
                         double &uH, double *duH)
    {
        double x = vpoint[0];
        double y = vpoint[1];
        
        double xi; double eta;
        // mapping to a box
        mapxy(x,y,xmax,xmin,ymax,ymin,xi,eta);
        // uH = a0+a1x+x2y+a3x^2+a4xy+a5y^2+...
        
        //double duHxi = 0;
        //double duHeta = 0;
        //int c=0;
        //uH = 0;
        
        //double addduHx; double addduHy;
        
        //std::cout<<"x "<<x<<" y "<<y<<std::endl;
        
        //std::cout<<"xmin "<<xmin<<"xmax "<<xmax<<std::endl;
        //std::cout<<"ymin "<<ymin<<"ymax "<<ymax<<std::endl;
        
        //std::cout<<"xi "<<xi<<" eta "<<eta<<std::endl;
        
        uH = 0; double duHxi = 0; double duHeta = 0;
        for (unsigned long i = 0; i<LS.size(); i++) //m1 = 0:k1
        {
            uH = uH+LS(i)*xyuH[i];
            duHxi = duHxi+LS(i)*xyduHx[i];
            duHeta = duHeta+LS(i)*xyduHy[i];
            //std::cout<<"LS "<<LS(i)<<std::endl;
        }
        
        /*
        for (int k1 = 0; k1<=p; k1++)  //k1 = 0:p
        {
            for (int m1 = 0; m1<=k1; m1++)  //m1 = 0:k1
            {
                c=c+1;
                uH = uH+LS(c-1)*pow(xi,(k1-m1))*pow(eta,(m1));
                //std::cout<<"uH "<<LS(c-1)<<" "<<pow(xi,(k1-m1))<<" "<<pow(eta,(m1))<<" "<<uH<<std::endl;
                if (k1-m1 ==0 )
                    addduHx=0;
                else
                    addduHx = (k1-m1)*LS(c-1)*pow(xi,(k1-m1-1))*pow(eta,(m1));
        
                duHxi = duHxi+addduHx;
         
                if (m1 == 0)
                    addduHy=0;
                else
                    addduHy = (m1)*LS(c-1)*pow(xi,(k1-m1))*pow(eta,(m1-1));
        
                duHeta = duHeta+addduHy;
            } //end m1
        }//end k1
        */
        //std::cout<<" uH "<<uH<<std::endl;
        //std::cout<<"duHxitest "<<duHxitest<<" duHxi "<<duHxi<<std::endl;
        //std::cout<<"duHetatest "<<duHetatest<<" duHeta "<<duHeta<<std::endl;
        
        double duHx = duHxi/(xmax-xmin);
        double duHy = duHeta/(ymax-ymin);
        duH[0] = duHx;
        duH[1] = duHy;
        
        //std::cout<<"uH "<<uH<<std::endl;
        //std::cout<<"duH "<<duH[0]<<" "<<duH[1]<<std::endl;
        
    }
    void dynamic::mapxy(int nElement,std::vector<std::vector<double> > xEle,
                        std::vector<std::vector<double> > yEle,
                        double xmax,double xmin,double ymax,double ymin,
                        double **x,
                        double **y)
    {
        for (int j = 0; j<nElement; j++)
        {
            //std::cout<<"j "<<j<<" xEleSize "<<xEle[j].size()<<std::endl;
            for (unsigned int q=0; q<xEle[j].size(); ++q)
            {
                //mapx[j].push_back((xEle[j][q]-xmin)/(xmax-xmin));
                //mapy[j].push_back((yEle[j][q]-ymin)/(ymax-ymin));
                
                x[j][q] = (xEle[j][q]-xmin)/(xmax-xmin);
                y[j][q] = (yEle[j][q]-ymin)/(ymax-ymin);
                //std::cout << "xmapymap  " << x[j][q]<<" "<< y[j][q]<<std::endl;
                //std::cout << "xEle  " << xEle[j][q]<<" "<< yEle[j][q]<<std::endl;
                //std::cout << "xmin  " << xmin<<" "<< xmax<<std::endl;
                //std::cout << "ymin  " << ymin<<" "<< ymax<<std::endl;
            }
        }
        
    }
    
    void dynamic::mapxy(double xEle,
                        double yEle,
                        double xmax,double xmin,double ymax,double ymin,
                        double &x,
                        double &y)
    {
        x = (xEle-xmin)/(xmax-xmin);
        y = (yEle-ymin)/(ymax-ymin);
        
    }
    
    void dynamic::matrixSolve(int mSize,FullMatrix<double> lhs,Vector<double> rhs,
                              Vector<double> &sol)
    {
        // use Eigen to solve dense matrix
        std::cout << "here 1.4.3.1 " <<mSize<<std::endl;
        Eigen::MatrixXf leftHS(mSize,mSize);
        Eigen::VectorXf rightHS(mSize);
        std::cout << "here 1.4.3.1 " <<std::endl;
        for (int i = 0; i<mSize; i++)
        {
            rightHS(i) = rhs(i);
            for (int j = 0; j<mSize; j++)
            {
                leftHS(i,j) = lhs(i,j);
            }
        }
        std::cout << "here 1.4.3.2 " <<std::endl;
        Eigen::VectorXf x = leftHS.colPivHouseholderQr().solve(rightHS);

        std::cout << "here 1.4.3.3 " <<std::endl;
        for (int i = 0; i<mSize; i++)
        {
            sol[i] = x(i);
        }

    }
    
    void dynamic::matrixSolve(int mSize,Vector<double> lhs,Vector<double> rhs,
                              double &sol)
    {
        
        // use Eigen to solve least square dense matrix
        MatrixXf A(mSize,1);
        VectorXf b(mSize);
        for (int i = 0; i<mSize; i++)
        {
            A(i) = lhs(i);
            b(i) = rhs(i);
        }
        
        MatrixXf x = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
        sol = x(0);
        
        //MatrixXf A = MatrixXf::Random(4, 1);
        //cout << "Here is the matrix A:\n" << A << endl;
        //VectorXf b = VectorXf::Random(4);
        //cout << "Here is the right hand side b:\n" << b << endl;
        //cout << "The least-squares solution is:\n"
        //<< x << endl;
        
        
    }
    
    double dynamic::getnorm(int sizeM,Vector<double> MM)
    {
        Eigen::VectorXf v(sizeM);
        for (int i = 0; i<sizeM; i++)
        {
            v[i] = MM(i);
        }
        return v.norm();
    }

}//end class


