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

#include "../include/post.h"
#include "../include/user_input.h"

#include <fstream>
#include <iostream>
#include <sstream>
namespace incompressible
{
    using namespace dealii;

    void Post::restart_save()
    {
    }
    
    void Post::savetovtk(hp::DoFHandler<2> &dof_handler,BlockVector<double> &present_solution)
    {
        
        int dim = 2;
        std::vector<std::string> solution_names(dim, "velocity");
        solution_names.push_back("pressure");
        
        std::vector<DataComponentInterpretation::DataComponentInterpretation>data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
        
        data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
        
        DataOut<2, hp::DoFHandler<2>> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(present_solution,
                                 solution_names,
                                 DataOut<2, hp::DoFHandler<2>>::type_dof_data,
                                 data_component_interpretation);
        data_out.build_patches();
        std::ofstream output("-solution-.vtk");
        data_out.write_vtk(output);
        
    }
    
    void Post::savesolution(BlockVector<double> present_solution, int istp)
    {
        using namespace std;
        
        const std::string filename =
        velocity_sol.str().c_str() + Utilities::int_to_string(istp, 1) + "-.txt";
        //std::ofstream output(filename);
        
        //const std::string filename =
        //vtkfilename.str().c_str() + Utilities::int_to_string(istp, 2) + ".vtk";
        
        //filename<< Utilities::int_to_string (istp, 1);
        //<< ".txt";
        
        //std::cout << filename.str() <<std::endl;
        FILE * pFile;
        pFile = fopen(filename.c_str(), "w");
        //std::cout << "save u_solution2 " <<std::endl;
        
        Vector<double> du_n = present_solution.block(0);
        
        //std::cout << "save u_solution3 " <<std::endl;
        for(size_t i=0; i<du_n.size(); ++i)
        {
            //cout<<du_n[i]<<endl;
            fprintf(pFile,"%6.16e ",du_n[i]);
            //cout<<du_n[i]<<endl;
        }
        fclose(pFile);
        //std::cout << "save u_solution " <<std::endl;
        
        //ostringstream filename1;
        //filename1 << Utilities::int_to_string (istp, 1);
        //<< ".txt";
        //std::cout << filename1.str() <<std::endl;
        
        const std::string filename1 =
        pressure_sol.str().c_str() + Utilities::int_to_string(istp, 1) + "-.txt";
        
        FILE * pFile1;
        pFile1 = fopen(filename1.c_str(), "w");
        
        Vector<double> p_n = present_solution.block(1);
        //std::cout << "present_solution " <<std::endl;
        
        for(size_t i=0; i<p_n.size(); ++i)
        {
            fprintf(pFile1,"%6.16e ",p_n[i]);
        }
        
        fclose(pFile1);
        //std::cout << "save p_solution " <<std::endl;
        
    }

    void Post::savemesh(int cycle,int vertices_size,int edge_size,hp::DoFHandler<2> &dof_handler,hp::FECollection<2> &fe, AffineConstraints<double> &nonzero_constraints)
    {
        using namespace std;
        
        //std::ostringstream filenamemesh("./mesh/", ios_base::app);
        //filenamemesh<< Utilities::int_to_string (cycle, 1)
        //<< "mesh.txt";
        const std::string filenamemesh =
        outputmeshfilename.str().c_str() + Utilities::int_to_string(cycle, 1) + "mesh.txt";
        
        std::cout << filenamemesh <<std::endl;
        FILE * pFile;
        pFile = fopen(filenamemesh.c_str(), "w");
        
        
        int cellnum = 0;
        //QGauss<2> quadrature_formula(degree*2);
        hp::QCollection<2>     quadrature_collection;
        for (unsigned int deg=min_degree; deg<=max_degree; ++deg)
        {
            quadrature_collection.push_back (QGauss<2>(deg*2));
        }
        
        hp::FEValues<2> hp_fe_values(fe,
        quadrature_collection,
        update_values | update_quadrature_points |
        update_JxW_values | update_gradients |
        update_inverse_jacobians);
        
        //const unsigned int dofs_per_cell = fe.dofs_per_cell;
        //std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        
        std::vector<types::global_dof_index> local_dof_indices;
        
        hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        hp::DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        std::vector<std::vector<int> > v_to_e_indices;
        v_to_e_indices.clear();
        v_to_e_indices.resize(vertices_size);
        
        std::vector<std::vector<double> > v_coordinate;
        v_coordinate.clear();
        v_coordinate.resize(vertices_size);
        
        //for (const auto &cell : dof_handler.active_cell_iterators())
        for (; cell != endc; ++cell)
        {
            cellnum++;
            //fe_values.reinit(cell);
            //cell->get_dof_indices(local_dof_indices);
            const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
            hp_fe_values.reinit(cell);
            
            //const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
            const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
            //std::cout<<"cellnum "<<cellnum-1<<" "<<dofs_per_cell<<std::endl;
            unsigned int n_q_points  = fe_values.n_quadrature_points;
            local_dof_indices.resize (dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);
            
            fprintf(pFile,"%s %d %d ","cellnum ",cellnum-1, dofs_per_cell);
            fprintf(pFile,"\n ");
            
            for(std::size_t i=0; i<dofs_per_cell; ++i)
            {
                
                //std::cout<<local_dof_indices[i]<<std::endl;
                fprintf(pFile,"%d ",local_dof_indices[i]);
                fprintf(pFile,"\n ");
            }
            
            
            
            int vertexNum = GeometryInfo<2>::vertices_per_cell;
            double* coordx = new double[vertexNum];
            double* coordy = new double[vertexNum];
            int* vertex_index = new int[vertexNum];
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                Point<2> &v = cell->vertex(i);
                coordx[i] =  v[0];
                coordy[i] = v[1];
                
                //vertices_to_dofs[cell->vertex_index(i)] = cell->vertex_dof_index(i,0,0);
                int v_index = cell->vertex_index(i);
                vertex_index[i] = v_index;
                //v_to_e_list[v_index].push_back(cell);
                v_to_e_indices[v_index].push_back(cell->active_cell_index());
                v_coordinate[v_index].push_back(v[0]);
                v_coordinate[v_index].push_back(v[1]);
                //std::cout << "v_index  " << v_index<<std::endl;
            }
            
            //std::cout<<"cellnum "<<cellnum-1<<std::endl;
            //std::cout << coordx[0]<<" "<< coordx[1]<<" "<< coordx[2]<<" "<< coordx[3]<<std::endl;
            //std::cout << coordy[0]<<" "<< coordy[1]<<" "<< coordy[2]<<" "<< coordy[3]<<std::endl;
            
            fprintf(pFile,"%f %f %f %f",coordx[0],coordx[1],coordx[2],coordx[3]);
            fprintf(pFile,"\n ");
            fprintf(pFile,"%f %f %f %f",coordy[0],coordy[1],coordy[2],coordy[3]);
            //fprintf(pFile,"\n ");
            //fprintf(pFile,"%d %d %d %d",vertex_index[0],vertex_index[1],vertex_index[2],vertex_index[3]);
            
            fprintf(pFile,"\n");
            
        } //end cell loop
        
        //fprintf(pFile,"%s %d","Grid_size ",cellnum-1);
        //fprintf(pFile,"\n");
        //fprintf(pFile,"%s %d","vertices_size ",vertices_size);
        //fprintf(pFile,"\n");
        //fprintf(pFile,"%s %d","total_DOF ",dof_handler.n_dofs());
        fclose(pFile);
        
        //std::ostringstream filenamev_to_e_indices("./mesh/", ios_base::app);
        //filenamev_to_e_indices<< Utilities::int_to_string (cycle, 1)
        //<< "v_to_e_indices.txt";
        const std::string filenamev_to_e_indices =
        outputmeshfilename.str().c_str() + Utilities::int_to_string(cycle, 1) + "v_to_e_indices.txt";
        
        //std::cout << filenamev_to_e_indices.str() <<std::endl;
        FILE * pFilev_to_e_indices;
        pFilev_to_e_indices = fopen(filenamev_to_e_indices.c_str(), "w");
        
        for (int i = 0; i<vertices_size; i++)
        {
            fprintf(pFilev_to_e_indices,"%d ",i);
            
            fprintf(pFilev_to_e_indices,"%f %f ",v_coordinate[i][0],v_coordinate[i][1]);
            
            for (unsigned long int j = 0; j<v_to_e_indices[i].size(); j++)
            {
                //std::cout << "v_to_e_indices  " << v_to_e_indices[i][j]<<std::endl;
                fprintf(pFilev_to_e_indices,"%d ",v_to_e_indices[i][j]);
            }
            
            fprintf(pFilev_to_e_indices,"\n ");
        }
        fclose(pFilev_to_e_indices);
        
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
        
        const std::string filenameIBC =
        outputmeshfilename.str().c_str() + Utilities::int_to_string(cycle, 1) + "IBC.txt";
        
        //std::cout << "nonzero_constraints lines in post"<<std::endl;
        std::filebuf fb;
        fb.open (filenameIBC,std::ios::out);
        std::ostream os(&fb);
        nonzero_constraints.print(os);
        fb.close();
        
        for (const auto &line : nonzero_constraints.get_lines())
        {
            //add boundary DOFs
            DOF_value[line.index][0] = 1.0;
            DOF_value[line.index][1] = line.inhomogeneity;
            //std::cout << "DoF in post" << line.index<<std::endl;
            // remove hanging nodes, they are not at boundaries
            for (const auto &entry : line.entries)
            {
                DOF_value[line.index][0] = 0;
                //std::cout << "DoF " << line.index <<" "<<entry.first<<std::endl;
            }
            
        }
        
        //for (unsigned long int i = 0; i<dof_handler.n_dofs(); i++)
        //{
            //std::cout << "DOF_value in post" << i<<" "<<DOF_value[i][0]<<" "<<DOF_value[i][1]<<std::endl;
        //}
        
        //std::ostringstream filenamebc("./mesh/", ios_base::app);
        //filenamebc<< Utilities::int_to_string (cycle, 1)
        //<< "bcEle.txt";
        
        const std::string filenameedge =
        outputmeshfilename.str().c_str() + Utilities::int_to_string(cycle, 1) + "EdgeInfo.txt";
        
        std::cout << filenameedge<<std::endl;
        FILE * EdgeFile;
        EdgeFile = fopen(filenameedge.c_str(), "w");
        
        const std::string filenamebc =
        outputmeshfilename.str().c_str() + Utilities::int_to_string(cycle, 1) + "bc.txt";
        
        std::cout << filenamebc<<std::endl;
        FILE * bcFile;
        bcFile = fopen(filenamebc.c_str(), "w");
        
        std::cout << "start saving mesh edge info"<<std::endl;
        std::vector<std::vector<int> > e_to_node;
        e_to_node.clear();
        e_to_node.resize(edge_size);
        
        cell = dof_handler.begin_active();
        endc = dof_handler.end();
        double xcord; double ycord;
        
        fprintf(EdgeFile,"%s %d","Grid_size ",cellnum);
        fprintf(EdgeFile,"\n");
        fprintf(EdgeFile,"%s %d","vertices_size ",vertices_size);
        fprintf(EdgeFile,"\n");
        fprintf(EdgeFile,"%s %d","Edge_size ",edge_size);
        fprintf(EdgeFile,"\n");
        fprintf(EdgeFile,"%s %d","total_DOF ",dof_handler.n_dofs());
        fprintf(EdgeFile,"\n");
        fprintf(EdgeFile,"Face_Vertex_Edge ");
        fprintf(EdgeFile,"\n");
        for (; cell!=endc; ++cell){
            int element_id = cell->active_cell_index();
            //std::cout << "element_id " << element_id<<std::endl;
            fprintf(EdgeFile,"%d ",element_id);
            fprintf(EdgeFile,"\n ");
            
            int vertexNum = GeometryInfo<2>::vertices_per_cell;
            int* vertex_index = new int[vertexNum];
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                int v_index = cell->vertex_index(i);
                vertex_index[i] = v_index;
            }
            fprintf(pFile,"%d %d %d %d",vertex_index[0],vertex_index[1],vertex_index[2],vertex_index[3]);
            fprintf(EdgeFile,"\n ");
            for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
            {
                int edge_id = cell->line_index(face_number);
                //std::cout << "edge_id " <<edge_id<<std::endl;
                fprintf(EdgeFile,"%d ",edge_id);
                int edge_vertex0 = cell->face(face_number)->vertex_index(0);
                int edge_vertex1 = cell->face(face_number)->vertex_index(1);
                    //std::cout <<edge_vertex0<<" "<<edge_vertex1<<std::endl;
                e_to_node[edge_id].push_back(edge_vertex0);
                e_to_node[edge_id].push_back(edge_vertex1);
                if (cell->face(face_number)->at_boundary())
                {   // get cellID
                    int cell_index = cell->active_cell_index();
                    int bcID = cell->face(face_number)->boundary_id();
                    fprintf(bcFile,"%d %d %d %d %d %d",cell_index,bcID,face_number,edge_id,edge_vertex0,edge_vertex1);
                    fprintf(bcFile,"\n");
                }
            }
            fprintf(EdgeFile,"\n");
        }
        
        /*
        fprintf(EdgeFile,"Edge_Vertex ");
        fprintf(EdgeFile,"\n ");
        for (int i = 0; i<edge_size; i++)
        {
            //fprintf(EdgeFile,"%d ",i);
            
            for (unsigned long int j = 0; j<2; j++)
            {
                //std::cout << "v_to_e_indices  " << v_to_e_indices[i][j]<<std::endl;
                fprintf(EdgeFile,"%d ",e_to_node[i][j]);
            }
            
            fprintf(EdgeFile,"\n");
        }
        */
        fclose(bcFile);
        std::cout << "finish saving bc mesh"<<std::endl;
        fclose(EdgeFile);
        std::cout << "finish saving mesh"<<std::endl;
        
    }// end function

    void Post::getTrueError(hp::FECollection<2> &fe,
                            Triangulation<2> &triangulation,
                            hp::DoFHandler<2> &dof_handler,
                            BlockVector<double> &n_solution,
                            Vector<float> &true_error_per_cell)
    {
        //create reference mesh
        Triangulation<2> triangulationRef;
        DoFHandler<2> dof_handlerRef(triangulationRef);
        //int degreeRef = 2;
        //int n_Ref = 80;
        
        BlockVector<double> ref_solution;
        createRefmesh(triangulationRef,dof_handlerRef,degreeRef,n_Ref,ref_solution);
        
        //read in reference solution
        readinRefSolution(ref_solution);
        
        //comput difference between current-solution and reference solution
        //double Error_u=0;
        //double Error_v = 0;
        double Error_p = 0;

        hp::QCollection<2>     quadrature_collection;
        
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
        hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        hp::DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        int cellnum = 0;

        std::vector<types::global_dof_index> local_dof_indices;
        for (; cell != endc; ++cell)
        {
            double Error_u = 0;
            double Error_v = 0;
            cellnum++;
            /*
            std::cout << "coordinate  ";
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                Point<2> &v = cell->vertex(i);
                std::cout << v[0]<<" "<< v[1]<<" ";
            }
            std::cout<<std::endl;
            */
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
            
            fe_values[velocities].get_function_values(n_solution,present_velocity_values);
            fe_values[velocities].get_function_gradients(n_solution, present_velocity_gradients);
            fe_values[pressure].get_function_values(n_solution,
                                                    present_pressure_values);
            fe_values[pressure].get_function_gradients(n_solution,
                                                       present_pressure_gradients);

            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                
                
                const Point<2> &vPoint=fe_values.quadrature_point(q);
                
                
                //std::cout << "vPoint  " << vPoint[0]<<" "<< vPoint[1]<<std::endl;
                
                //std::cout << "solution_value  " << present_velocity_values[q][0]<<" "<< present_velocity_values[q][1]<<" "<<present_pressure_values[q]<<std::endl;
                
                // map the point to reference mesh
                /*
                MappingQ1<2> qmap;
                std::pair<DoFHandler<2>::active_cell_iterator, Point<2> > cell_iterator;
                
                cell_iterator = GridTools::find_active_cell_around_point(qmap, dof_handlerRef, vPoint);
                
                DoFHandler<2>::active_cell_iterator ref_cell = cell_iterator.first;
                
                std::cout << "ref_cell coordinate  ";
                for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
                {
                    Point<2> &v = ref_cell->vertex(i);
                    std::cout << v[0]<<" "<< v[1]<<" ";
                }
                std::cout<<std::endl;
                */
                
                Vector<double> ref_velocity_value(3);
                VectorTools::point_value(dof_handlerRef,
                             ref_solution,
                             vPoint,
                             ref_velocity_value);
                
                
                
                //std::cout << "referenc_value  " << ref_velocity_value[0]<<" "<< ref_velocity_value[1]<<std::endl;
                
                //std::cout << "present_velocity_values  " << present_velocity_values[q][0]<<" "<< present_velocity_values[q][1]<<std::endl;
                double u = ref_velocity_value[0];
                double v = ref_velocity_value[1];
                double p = ref_velocity_value[2];
                
                Error_u = Error_u+(present_velocity_values[q][0]-ref_velocity_value[0])*(present_velocity_values[q][0]-ref_velocity_value[0])*fe_values.JxW(q);
                
                Error_v = Error_v+(present_velocity_values[q][1]-ref_velocity_value[1])*(present_velocity_values[q][1]-ref_velocity_value[1])*fe_values.JxW(q);
                
                //std::cout<<" Error_quad "<<(present_velocity_values[q][1]-ref_velocity_value[1])*(present_velocity_values[q][1]-ref_velocity_value[1])*fe_values.JxW(q)<<std::endl;
                
                //std::cout << "differ_value  " << present_velocity_values[q][0]-ref_velocity_value[0]<<" "<< present_velocity_values[q][1]-ref_velocity_value[1]<<" "<<ref_velocity_value[2]-p<<std::endl;
                
                //std::cout << "Error_u  " <<Error_u<<std::endl;
                Error_p = Error_p+(present_pressure_values[q]-p)*(present_pressure_values[q]-p)*fe_values.JxW(q);
                
                /*
                Error_ux = Error_ux+(present_velocity_gradients[q][0][0]-ux)*(present_velocity_gradients[q][0][0]-ux)*fe_values.JxW(q);
                Error_vx = Error_vx+(present_velocity_gradients[q][1][0]-vx)*(present_velocity_gradients[q][1][0]-vx)*fe_values.JxW(q);
                Error_px = Error_px+(present_pressure_gradients[q][0]-px)*(present_pressure_gradients[q][0]-px)*fe_values.JxW(q);
                
                Error_uy = Error_uy+(present_velocity_gradients[q][0][1]-uy)*(present_velocity_gradients[q][0][1]-uy)*fe_values.JxW(q);
                Error_vy = Error_vy+(present_velocity_gradients[q][1][1]-vy)*(present_velocity_gradients[q][1][1]-vy)*fe_values.JxW(q);
                Error_py = Error_py+(present_pressure_gradients[q][1]-py)*(present_pressure_gradients[q][1]-py)*fe_values.JxW(q);
                
                */
            }
            double L2_ref_u = pow((Error_u+Error_v),0.5);
            double L2_ref_p = pow((Error_p),0.5);
            //std::cout<<" L2_true_u "<<L2_ref_u<<std::endl;
            true_error_per_cell[cellnum-1] = L2_ref_u;
            //std::cout<<" Error_v "<<Error_v<<std::endl;
        }
        
        //double L2_ref_u = pow((Error_u+Error_v),0.5);
        //double L2_ref_p = pow((Error_p),0.5);
        
        //H1_velocity  = pow((Error_ux+Error_uy+Error_vx+Error_vy),0.5);
        //H1_pressure  = pow((Error_px+Error_py),0.5);
        //std::cout<<" total_true_error "<<L2_ref_u<<" "<<L2_ref_p<<std::endl;
    }
    
    void Post::createRefmesh(Triangulation<2> &triangulation,DoFHandler<2> &dof_handler,int degreeRef,int n_Ref,BlockVector<double> &n_solution)
    {
        std::vector<unsigned int> subdivisions (2, 2);
        double x1;double x2;double y1;double y2;
        
        //defineMesh(subdivisions, x1,x2,y1,y2);
        subdivisions[0] = n_Ref;
        subdivisions[1] = n_Ref;
        
        x1=0.0;x2=1.0;y1=0.0;y2=1.0;
        const Point<2> bottom_left= (Point<2>(x1,y1));
        const Point<2> top_right= (Point<2>(x2,y2));
        
        GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                   subdivisions,
                                                   bottom_left,
                                                   top_right);
        std::ofstream out("grid-1.eps");
        GridOut       grid_out;
        grid_out.write_eps(triangulation, out);
        std::cout << "Grid written to grid-1.eps" << std::endl;
        
        FE_Q<2> u(degreeRef);
        FE_Q<2> p(degreeRef);
        FESystem<2> fe(u,2, p,1);
        std::vector<types::global_dof_index> dofs_per_block;
        
        dof_handler.distribute_dofs(fe);
        std::vector<unsigned int> block_component(dim + 1, 0);
        block_component[dim] = 1;
        DoFRenumbering::component_wise(dof_handler, block_component);
        dofs_per_block.resize(2);
        DoFTools::count_dofs_per_block(dof_handler, dofs_per_block, block_component);
        unsigned int dof_u = dofs_per_block[0];
        unsigned int dof_p = dofs_per_block[1];
        
        std::cout << "   Number of active cells in reference solution: " << triangulation.n_active_cells()
        << std::endl
        << "   Number of degrees of freedom in reference solution: " << dof_handler.n_dofs()
        << " (" << dof_u << '+' << dof_p << ')' << std::endl;
        
        n_solution.reinit(dofs_per_block);

    }
    
    void Post::readinRefSolution(BlockVector<double> &n_solution)
    {
        using namespace std;
        //ostringstream filenameRefu;
        //filenameRefu << "./solution/"
        //<< "test-un"
        //<< ".txt";
        
        FILE * pFile;
        pFile = fopen(filenameRefu.str().c_str(), "r+");
        
        Vector<double> u_n = n_solution.block(0);
        
        //cout<<" u_n.size"<<u_n.size()<<endl;
        
        double solution_u;
        for(size_t i=0; i<u_n.size(); ++i)
        {
            int ret = fscanf(pFile, "%le", &solution_u);
            u_n[i] = solution_u;
            //cout<<" u_n[i] "<<u_n[i]<<endl;
            if(ret == EOF) {
                break;
            }
        }
        fclose(pFile);
        
        //ostringstream filenameRefp;
        //filenameRefp << "./solution/"
        //<< "test-pn"
        //<< ".txt";
        
        FILE * pFile1;
        pFile1 = fopen(filenameRefp.str().c_str(), "r+");
        
        Vector<double> p_n = n_solution.block(1);
        
        //cout<<" p_n size "<<p_n.size()<<endl;
        
        double solution_p;
        for(size_t i=0; i<p_n.size(); ++i)
        {
            int ret = fscanf(pFile1, "%le", &solution_p);
            p_n[i] = solution_p;
            //cout<<" p_n[i] "<<p_n[i]<<endl;
            if(ret == EOF) {
                break;
            }
        }
        fclose(pFile1);

        n_solution.block(0) = u_n;
        n_solution.block(1) = p_n;
    }
    /*
    void Post::ErrorEstimate(FESystem<2> fe, DoFHandler<2> &dof_handler,
                             int degree,
                             BlockVector<double> n_solution,
                             double &L2_velocity,double &H1_velocity,
                             double &L2_presure,double &H1_pressure)
    {
        
        QGauss<2> quadrature_formula(degree + 2);
        
        FEValues<2> fe_values(fe, quadrature_formula,
                              update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_inverse_jacobians);
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        //const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points    = quadrature_formula.size();
        
        std::vector<Tensor<1, 2>> present_velocity_values(n_q_points);
        std::vector<Tensor<2, 2>> present_velocity_gradients(n_q_points);
        std::vector<double>       present_pressure_values(n_q_points);
        std::vector<Tensor<1, 2>> present_pressure_gradients(n_q_points);
        
        double u=0.0; double v=0.0; double p=0.0;
        double ux=0.0; double vx=0.0; double px=0.0;
        double uy=0.0; double vy=0.0; double py=0.0;
        double Error_u =0.0; double Error_v =0.0; double Error_p =0.0;
        double Error_ux =0.0; double Error_vx =0.0; double Error_px =0.0;
        double Error_uy =0.0; double Error_vy =0.0; double Error_py =0.0;
        
        for (; cell != endc; ++cell)
        {
            fe_values.reinit(cell);
            
            fe_values[velocities].get_function_values(n_solution,present_velocity_values);
            fe_values[velocities].get_function_gradients(n_solution, present_velocity_gradients);
            fe_values[pressure].get_function_values(n_solution,
                                                    present_pressure_values);
            fe_values[pressure].get_function_gradients(n_solution,
                                                       present_pressure_gradients);
            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                //getCoordinate(coordx,coordy,q,x,y);
                const Point<2> &vPoint=fe_values.quadrature_point(q);
                
                getExact(vPoint[0],vPoint[1],u,v,p,ux,vx,px,uy,vy,py);
                
                //std::cout << "coordinate  " << vPoint[0]<<" "<< vPoint[1]<<std::endl;
                
                Error_u = Error_u+(present_velocity_values[q][0]-u)*(present_velocity_values[q][0]-u)*fe_values.JxW(q);
                Error_v = Error_v+(present_velocity_values[q][1]-v)*(present_velocity_values[q][1]-v)*fe_values.JxW(q);
                Error_p = Error_p+(present_pressure_values[q]-p)*(present_pressure_values[q]-p)*fe_values.JxW(q);
                
                Error_ux = Error_ux+(present_velocity_gradients[q][0][0]-ux)*(present_velocity_gradients[q][0][0]-ux)*fe_values.JxW(q);
                Error_vx = Error_vx+(present_velocity_gradients[q][1][0]-vx)*(present_velocity_gradients[q][1][0]-vx)*fe_values.JxW(q);
                Error_px = Error_px+(present_pressure_gradients[q][0]-px)*(present_pressure_gradients[q][0]-px)*fe_values.JxW(q);
                
                Error_uy = Error_uy+(present_velocity_gradients[q][0][1]-uy)*(present_velocity_gradients[q][0][1]-uy)*fe_values.JxW(q);
                Error_vy = Error_vy+(present_velocity_gradients[q][1][1]-vy)*(present_velocity_gradients[q][1][1]-vy)*fe_values.JxW(q);
                Error_py = Error_py+(present_pressure_gradients[q][1]-py)*(present_pressure_gradients[q][1]-py)*fe_values.JxW(q);
                
                
            }
        }
        
        L2_velocity = pow((Error_u+Error_v),0.5);
        L2_presure = pow((Error_p),0.5);
        
        H1_velocity  = pow((Error_ux+Error_uy+Error_vx+Error_vy),0.5);
        H1_pressure  = pow((Error_px+Error_py),0.5);
        
    }
    
    void Post::getExact(double x,double y,
                        double &u,double &v,double &p,
                        double &ux,double &vx,double &px,double &uy,double &vy,
                        double &py)
    {
        
        double pi = M_PI;
        double lam = 0.5*Re-sqrt((0.25*Re*Re+4*pi*pi));
        
        u = 1-exp(lam*x)*cos(2*pi*y);
        v = (lam/(2*pi))*exp(lam*x)*sin(2*pi*y);
        p = 0.5*(1-exp(2*lam*x));
        ux = -lam*exp(lam*x)*cos(2*pi*y);
        uy = 2*pi*exp(lam*x)*sin(2*pi*y);
        vx = (lam*lam*exp(lam*x)*sin(2*pi*y))/(2*pi);
        vy = lam*exp(lam*x)*cos(2*pi*y);
        px = -lam*exp(2*lam*x);
        py = 0;
    }
     */
}


