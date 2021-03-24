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
#include "../include/post.h"
#include <math.h>
namespace incompressible
{
    using namespace dealii;

    Post::Post(int degree)
    :
    degree(degree)
    {}
    
    void Post::restart_save()
    {
    }
    
    void Post::savetovtk(DoFHandler<2> &dof_handler,BlockVector<double> present_solution,unsigned int istp)
    {
        
        int dim = 2;
        std::vector<std::string> solution_names(dim, "velocity");
        solution_names.push_back("pressure");
        
        std::vector<DataComponentInterpretation::DataComponentInterpretation>data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
        
        data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
        
        DataOut<2> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(present_solution,
                                 solution_names,
                                 DataOut<2>::type_dof_data,
                                 data_component_interpretation);
        data_out.build_patches();
        
        std::ostringstream filename;
        filename << "./paraview/"
        << 100
        << "-solution-"
        << Utilities::int_to_string (istp, 4)
        << ".vtk";
        //std::ofstream output("-solution-.vtk");
        std::ofstream output (filename.str().c_str());
        data_out.write_vtk(output);
        
    }
    
    void Post::savesolution(BlockVector<double> present_solution,BlockVector<double> n_solution_time_derivative,unsigned int istp)
    {
        using namespace std;
        ostringstream filename;
        filename << "./solution/"
        << 100
        << "-un-"
        << Utilities::int_to_string (istp, 1)
        << ".txt";
        
        //std::cout << "save u_solution1 " <<std::endl;
        FILE * pFile;
        pFile = fopen(filename.str().c_str(), "w");
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
        
        ostringstream filename1;
        filename1 << "./solution/"
        << 100
        << "-pn-"
        << Utilities::int_to_string (istp, 1)
        << ".txt";
        
        FILE * pFile1;
        pFile1 = fopen(filename1.str().c_str(), "w");
        
        Vector<double> p_n = present_solution.block(1);
        //std::cout << "present_solution " <<std::endl;
        
         for(size_t i=0; i<p_n.size(); ++i)
         {
         fprintf(pFile1,"%6.16e ",p_n[i]);
         }
         
        fclose(pFile1);
        //std::cout << "save p_solution " <<std::endl;
        
        ostringstream filename2;
        filename2 << "./solution/"
        << 100
        << "-ut_n-"
        << Utilities::int_to_string (istp, 1)
        << ".txt";
        
        FILE * pFile2;
        pFile2 = fopen(filename2.str().c_str(), "w");
        
        Vector<double> ut_n = n_solution_time_derivative.block(0);
        //std::cout << "present_solution " <<std::endl;
        for(size_t i=0; i<ut_n.size(); ++i)
        {
            fprintf(pFile2,"%6.16e ",ut_n[i]);
            //cout<<du_n[i]<<endl;
        }
        fclose(pFile2);
        //std::cout << "save ut_solution " <<std::endl;
        
        ostringstream filename3;
        filename3 << "./solution/"
        << 100
        << "-pt_n-"
        << Utilities::int_to_string (istp, 1)
        << ".txt";
        
        FILE * pFile3;
        pFile3 = fopen(filename3.str().c_str(), "w");
        
        Vector<double> pt_n = n_solution_time_derivative.block(1);
        
         for(size_t i=0; i<pt_n.size(); ++i)
         {
         fprintf(pFile3,"%6.16e ",pt_n[i]);
         }
        
        fclose(pFile3);
        //std::cout << "save pt_solution " <<std::endl;
        
    }
    
    void Post::savemesh(int cycle,DoFHandler<2> &dof_handler,FESystem<2> &fe,MappingQ<2> mapping, AffineConstraints<double> &nonzero_constraints)
    {
        
        //int dim = 2;
        //Vector<float> error_efficiency_per_cell (triangulation.n_active_cells());
        //Vector<float> fe_degrees(triangulation.n_active_cells());
        
        using namespace std;
        
        std::ostringstream filenamemesh("./mesh/", ios_base::app);
        filenamemesh<< Utilities::int_to_string (cycle, 1)
        << "mesh.txt";
        
        std::cout << filenamemesh.str() <<std::endl;
        FILE * pFile;
        pFile = fopen(filenamemesh.str().c_str(), "w");
        
        
        int cellnum = 0;
        QGauss<2> quadrature_formula(degree*2);
        
        
        FEValues<2> fe_values(mapping, fe,quadrature_formula,
                              update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_inverse_jacobians);
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        
        
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            cellnum++;
            fe_values.reinit(cell);
            cell->get_dof_indices(local_dof_indices);
            
            //std::cout<<"cellnum "<<cellnum-1<<" "<<dofs_per_cell<<std::endl;
            
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
            
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                Point<2> &v = cell->vertex(i);
                coordx[i] =  v[0];
                coordy[i] = v[1];
            }
            
            //std::cout<<"cellnum "<<cellnum-1<<std::endl;
            //std::cout << coordx[0]<<" "<< coordx[1]<<" "<< coordx[2]<<" "<< coordx[3]<<std::endl;
            //std::cout << coordy[0]<<" "<< coordy[1]<<" "<< coordy[2]<<" "<< coordy[3]<<std::endl;
            
            fprintf(pFile,"%f %f %f %f",coordx[0],coordx[1],coordx[2],coordx[3]);
            fprintf(pFile,"\n ");
            fprintf(pFile,"%f %f %f %f",coordy[0],coordy[1],coordy[2],coordy[3]);
            fprintf(pFile,"\n");
            
        } //end cell loop
        
        fclose(pFile);
        
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
        
        
        std::cout << "nonzero_constraints lines in post"<<std::endl;
        std::filebuf fb;
        fb.open ("./mesh/pbc.txt",std::ios::out);
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
        
        /*
        for (unsigned long int i = 0; i<dof_handler.n_dofs(); i++)
        {
            //std::cout << "DOF_value in post" << i<<" "<<DOF_value[i][0]<<" "<<DOF_value[i][1]<<std::endl;
        }
        */
        const std::string filenamebc ="./mesh/bcEle.txt";
        
        std::cout << filenamebc<<std::endl;
        FILE * bcFile;
        bcFile = fopen(filenamebc.c_str(), "w");
        
        //QGauss<2> quadrature_formula(degree*2);
        
        //FEValues<2> fe_values(mapping, fe, quadrature_formula,
          //                    update_values | update_quadrature_points |
            //                  update_JxW_values | update_gradients |
              //                update_inverse_jacobians);
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        for (; cell!=endc; ++cell){
            int element_id = cell->active_cell_index();
            
            int vertexNum = GeometryInfo<2>::vertices_per_cell;
            int* vertex_index = new int[vertexNum];
            for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            {
                int v_index = cell->vertex_index(i);
                vertex_index[i] = v_index;
            }
            for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
            {
                int edge_id = cell->line_index(face_number);
                int edge_vertex0 = cell->face(face_number)->vertex_index(0);
                int edge_vertex1 = cell->face(face_number)->vertex_index(1);
                if (cell->face(face_number)->at_boundary())
                {   // get cellID
                    int cell_index = cell->active_cell_index();
                    int bcID = cell->face(face_number)->boundary_id();
                    fprintf(bcFile,"%d %d %d ",cell_index,bcID,face_number);
                    fprintf(bcFile,"\n");
                }
            }
        }

        fclose(bcFile);
        
    }// end function
    
    void Post::ErrorEstimate(MappingQ<2> mapping, FESystem<2> fe, DoFHandler<2> &dof_handler,
                             int EleNum,double*** div_SF_u,double**** grad_SF_u,
                             double**** SF_u, double*** SF_p,double**** grad_SF_p,BlockVector<double> n_solution,
                             double &L2_velocity,double &H1_velocity,
                             double &L2_presure,double &H1_pressure)
    {
        
        QGauss<2> quadrature_formula(degree*2);
        
        FEValues<2> fe_values(mapping, fe, quadrature_formula,
                              update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_inverse_jacobians);
        
        DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
        //DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
        
        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(2);
        
        //const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points    = quadrature_formula.size();
        
        std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
        std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
        std::vector<double>       present_pressure_values(n_q_points);
        std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
        
        double u=0.0; double v=0.0; double p=0.0;
        double ux=0.0; double vx=0.0; double px=0.0;
        double uy=0.0; double vy=0.0; double py=0.0;
        double Error_u =0.0; double Error_v =0.0; double Error_p =0.0;
        double Error_ux =0.0; double Error_vx =0.0; double Error_px =0.0;
        double Error_uy =0.0; double Error_vy =0.0; double Error_py =0.0;
        
        for (int cellnum=0; cellnum<EleNum; ++cellnum)
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
        L2_presure =pow((Error_p),0.5);
        
        H1_velocity  = pow((Error_ux+Error_uy+Error_vx+Error_vy),0.5);
        H1_pressure  = pow((Error_px+Error_py),0.5);
        
    }
    
    void Post::getExact(double x,double y,double &u,double &v,double &p,
                        double &ux,double &vx,double &px,double &uy,double &vy,
                        double &py)
    {
        double U = 0.2;
        double L = 0.2;
        double miut = 0.1;
        
        double sumsine = 0.0;
        double diffsumsine = 0.0;
        for ( int n=1; n<100; ++n)
        {
            sumsine = sumsine+
            sin(n*M_PI/L)*y*exp(-(n*M_PI)*(n*M_PI)*miut)*(2/(n*M_PI))*pow(-1,n+1);
            
            diffsumsine = diffsumsine+
            sin(n*M_PI/L)*exp(-(n*M_PI)*(n*M_PI)*miut)*(2/(n*M_PI))*pow(-1,n+1);
        }
        u = (y/L-sumsine)*U;
        v = 0;
        p = 0;
        ux = 0;
        vx = 0;
        px = 0;
        uy = (1/L-diffsumsine)*U;
        vy = 0;
        py = 0;
    }
}


