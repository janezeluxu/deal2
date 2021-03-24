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
#include "../include/pre.h"
#include "../include/dynamic.h"
#include "../include/user_input.h"
#include "../include/post.h"
namespace incompressible
{
    using namespace dealii;
    NavierStokes::NavierStokes(const unsigned int degree,const unsigned int mapdegree)
    :
    //triangulation(triangule),
    fe(FE_Q<2>(degree), 2, FE_Q<2>(degree), 1),
    triangulation(Triangulation<2>::maximum_smoothing),
    mapping(MappingQ<2> (mapdegree)),
    dof_handler(triangulation),
    degree(degree),
    mapdegree(mapdegree)
    //dof_handler(triangulation)
    {}
    
    void NavierStokes::setup_meshDOF()
    {
        extern int meshRefineLevel ;
        createmesh(triangulation,mapping,meshRefineLevel);
        setDof();
        //initialize(dof_handler);
    }
    
    void NavierStokes::setDof()
    {
        /*
         set up degree of freedoms and Dirichlet BCs using boundary ID
         */
        int dim = 2;
        //system_matrix.clear();
        //pressure_mass_matrix.clear();
        dof_handler.distribute_dofs(fe);
        std::vector<unsigned int> block_component(dim + 1, 0);
        block_component[dim] = 1;
        DoFRenumbering::component_wise(dof_handler, block_component);
        dofs_per_block.resize(2);
        DoFTools::count_dofs_per_block(dof_handler, dofs_per_block, block_component);
        unsigned int dof_u = dofs_per_block[0];
        unsigned int dof_p = dofs_per_block[1];
        
        std::cout << "   Number of active cells: " << triangulation.n_active_cells()
        << std::endl
        << "   Number of degrees of freedom: " << dof_handler.n_dofs()
        << " (" << dof_u << '+' << dof_p << ')' << std::endl;
        
        system_matrix.clear();

        BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
        DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
        sparsity_pattern.copy_from(dsp);
        system_matrix.reinit(sparsity_pattern);
        
        n_solution.reinit(dofs_per_block);
        n_solution_time_derivative.reinit(dofs_per_block);
        np1_solution.reinit(dofs_per_block);
        np1_solution_time_derivative.reinit(dofs_per_block);
        
        npaf_solution.reinit(dofs_per_block);
        npam_solution_time_derivative.reinit(dofs_per_block);
        
        newton_update.reinit(dofs_per_block);
        system_rhs.reinit(dofs_per_block);
    }
    
void NavierStokes::timeloop()
{
    viscosity = 1/Re;
    extern unsigned int nstp ;
    extern unsigned int num_start ;
    std::cout << "viscosity in NavierStoke " << viscosity<<std::endl;
    
    int ele_num = triangulation.n_active_cells();
    double*** diffusive_k = new double**[ele_num];
    double*** div_SF_u = new double**[ele_num];
    double**** grad_SF_u = new double***[ele_num];
    double**** SF_u = new double***[ele_num];
    double*** SF_p = new double**[ele_num];
    double**** grad_SF_p = new double***[ele_num];
    double*** gijG = new double**[ele_num];
    
    initialize(diffusive_k,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG);
    pre_compute(diffusive_k,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG);
    
    startingTimeloop(num_start);
    Post processResult(degree);
    int cycle = 1;
    applyBC(0);
    processResult.savetovtk(dof_handler,n_solution,0000);
    //std::cout << "nonzero_constraints lines "<<std::endl;
    //nonzero_constraints.print(std::cout);
    processResult.savemesh(cycle,dof_handler,fe,mapping, nonzero_constraints);
    
    // -------- prepare dynamic mesh info---------- //
    double *cellArea;
    int Grid_size = triangulation.n_active_cells();
    cellArea = new double [Grid_size];
    double *Tave = new double [Grid_size];
    std::cout << "dynamic system "<<std::endl;
    dynamic tauSystem(viscosity,
                      dof_handler,fe,
                      triangulation);
    
    std::map<int,int> vertices_to_dofs;
    std::vector<std::vector<DoFHandler<2>::active_cell_iterator> > v_to_e_list;
    std::vector<std::vector<int> > v_to_e_indices;
    
    v_to_e_list.clear();
    v_to_e_indices.clear();
    auto vertices = triangulation.get_vertices();
    v_to_e_list.resize(vertices.size());
    v_to_e_indices.resize(vertices.size());
    //double vertex_value_u_map[vertices.size()][2];
    double **vertex_value_u_map;
    double **vertex_value_v_map;
    double **vertex_value_p_map;
    vertex_value_u_map = new double *[vertices.size()];
    vertex_value_v_map = new double *[vertices.size()];
    vertex_value_p_map = new double *[vertices.size()];
    for(unsigned long i = 0; i <vertices.size(); i++)
    {
        vertex_value_u_map[i] = new double[3];
        vertex_value_v_map[i] = new double[3];
        vertex_value_p_map[i] = new double[3];
    }
    
    std::vector<int> IBC;
    double *bcTag = new double[vertices.size()];
    for(unsigned long i = 0; i <vertices.size(); i++)
    {
        bcTag[i] = surfaceTag;
    }
    
    tauSystem.buildDynamicMeshinfo(vertices_to_dofs,v_to_e_list,v_to_e_indices,nonzero_constraints,
        vertex_value_u_map,vertex_value_v_map,vertex_value_p_map,cellArea,IBC,bcTag);
    
    double **lhsG = new double *[vertices.size()];
    double **rhsG = new double *[vertices.size()];
    double **xyuHG = new double *[vertices.size()];
    double **xyduHxG = new double *[vertices.size()];
    double **xyduHyG = new double *[vertices.size()];

    tauSystem.additionalDynamicinfo(v_to_e_list, cellArea,
                                    lhsG, rhsG,
                                    xyuHG, xyduHxG, xyduHyG);
    
    double **taum_g = new double *[Grid_size];
    double ** taum = new double*[Grid_size];
    for(int i = 0; i < Grid_size; ++i)
    {
        taum[i] = new double[4];
        taum_g[i] = new double[4];
    }
        
    double *taumEle = new double [Grid_size];
    std::cout << "initialize taum "<<std::endl;
    double vmax = 1.0;
    tauSystem.getInitialtau(Grid_size, vmax,cellArea,v_to_e_list,n_solution,taumEle, taum);
    tauSystem.getTave(Grid_size,gijG,n_solution,Tave);
    /*
    for (int i = 0; i<Grid_size; i++)
    {
        for (int v = 0; v<4; v++)
        {
            std::cout << "initialize taum "<<taum[i][v]<<std::endl;
        }
    }
     */
    // -------- prepare dynamic mesh info end ---------- //
    for (unsigned int istp = num_start; istp < nstp+num_start; ++istp)
    {
        std::cout << "-------start timestep-------" <<istp<< std::endl;
        double t = istp*dt;
        applyBC(t);

        //processResult.savetovtk(dof_handler,n_solution,1000);
        predictor();
        //processResult.savetovtk(dof_handler,np1_solution,2000);
        //Newton's iteration
        corrector(diffusive_k,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG,taumEle);

        //processResult.savetovtk(dof_handler,n_solution,3000);
        
        updator();
        //processResult.savetovtk(dof_handler,n_solution,4000);

        
        //----------- dynamic update tau --------------//
        if (dynamicFlag)
        {
            std::cout << "-------start dynamic-------" << std::endl;
            if (istp+1 % 10 == 0)
            {
                std::cout << "-------get Tave-------" << std::endl;
                tauSystem.getTave(Grid_size,gijG,n_solution,Tave);
            }
            
            tauSystem.getc1c2(n_solution,n_solution_time_derivative,
                              v_to_e_list,vertex_value_u_map,
                              vertex_value_v_map,vertex_value_p_map,
                              cellArea,IBC,lhsG, rhsG,
                              xyuHG, xyduHxG, xyduHyG,taum_g);
            
            // update taum element based on relaxation
            tauSystem.relxation_ele(Grid_size,taum_g,Tave,taum,taumEle);
        }
        //----------- dynamic update end --------------//
        
        
        if (istp % 10 == 0)
        {
        //Post processResult(degree);
        std::cout << "-------save vtk-------" << std::endl;
        processResult.savetovtk(dof_handler,n_solution,istp);
        std::cout << "-------save solution-------" << std::endl;
        processResult.savesolution(n_solution,n_solution_time_derivative,istp+1);
        }
        
    }
}

void NavierStokes::predictor()
{

    //ut_np1[i] = ((gamma-1)/gamma)*ut_n[i];
    np1_solution_time_derivative = n_solution_time_derivative;
    np1_solution_time_derivative.block(0) *= (gamma-1)/gamma;
    np1_solution_time_derivative.block(1) =0;
    
    BlockVector<double>      tmp;
    tmp = n_solution_time_derivative;
    np1_solution = np1_solution_time_derivative;
    
    //u_np1[i] = u_n[i]+dt*ut_n[i]+gamma*dt*(ut_np1[i]-ut_n[i]);
    np1_solution.block(0) -= n_solution_time_derivative.block(0);
    np1_solution.block(0) *= gamma*dt;
    tmp.block(0) *= dt;
    np1_solution.block(0) += tmp.block(0);
    np1_solution.block(0) += n_solution.block(0);
    
    //p_np1[i] = p_n[i]+dt*pt_n[i]+dt*(pt_np1[i]-pt_n[i]);
    //np1_solution.block(1) = np1_solution_time_derivative.block(1);
    np1_solution.block(1) -= n_solution_time_derivative.block(1);
    np1_solution.block(1) *= dt;
    tmp.block(1) *=dt;
    np1_solution.block(1) += tmp.block(1);
    np1_solution.block(1) += n_solution.block(1);
    
}

void NavierStokes::corrector(double *** &diffusive_k,
                                double*** &div_SF_u,double**** &grad_SF_u,
                                double**** &SF_u, double*** &SF_p,
                                double**** &grad_SF_p,double*** &gijG,
                                double *taumEle)
{
    /*
     newtons loop with a very easy line search algorithm
     difference between first iteration and not first is for wich constraint to use
     */
    
    double current_res=0;
    viscosity = 1/Re;
    //double last_res;
    //std::cout << "viscosity  in generalizeA " << viscosity<<std::endl;
    for(int i=0; i<maxNewtonIter; ++i)
    {
        //bool isfirst = FALSE;
        std::cout << "----start iteration Number ----" <<i<< std::endl;

        iterPC();
        //std::cout << "----get newton update ----" << std::endl;
        
        setup_system(diffusive_k,div_SF_u,grad_SF_u,SF_u,SF_p,grad_SF_p,gijG,taumEle);
        linear_solve(newton_update);
        
        //Post processResult(degree);
        //processResult.savetovtk(dof_handler,newton_update,100001);
        
        //std::cout << "-------now update flow-------" << std::endl;
        update_flow();
        
        //std::cout << "-------now get residual-------" << std::endl;
        
        getResidual(current_res);
        if (current_res < tolerance)
            break;
    }
    
    
}

void NavierStokes::iterPC()
{
    
    npaf_solution = np1_solution;
    npam_solution_time_derivative= np1_solution_time_derivative;
    
    npaf_solution.block(0) -= n_solution.block(0);
    npaf_solution.block(0) *= alphaf;
    npaf_solution.block(0) += n_solution.block(0);
    npaf_solution.block(1) = np1_solution.block(1);
    
    npam_solution_time_derivative.block(0) -= n_solution_time_derivative.block(0);
    npam_solution_time_derivative.block(0) *= alpham;
    npam_solution_time_derivative.block(0) += n_solution_time_derivative.block(0);
    npam_solution_time_derivative.block(1) = np1_solution_time_derivative.block(1);

}

void NavierStokes::setup_system(double *** diffusive_k,
                               double*** &div_SF_u,double**** &grad_SF_u,
                               double**** &SF_u, double*** &SF_p,
                               double**** &grad_SF_p,double*** &gijG,
                               double *taumEle)
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

    Vector<double>     rhstest(dofs_per_cell);
    //std::vector<Tensor<1, 2>>   Rm(dofs_per_cell);
    FullMatrix<double> lhstest(dofs_per_cell, dofs_per_cell);
    std::vector<Tensor<1, 2> >        lhs1(dofs_per_cell);
    std::vector<Tensor<1, 2> >        lhs2(dofs_per_cell);
    std::vector<Tensor<1, 2> >        lhs3(dofs_per_cell);
    std::vector<double >        lhs4(dofs_per_cell);
    std::vector<Tensor<1, 2> >        lhs5(dofs_per_cell);
    std::vector<Tensor<1, 2> >        rhs1(n_q_points);
    
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
    DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    int cellnum = 0;
    
    double* gij = new double[4];
    double rc,tauM,tauC;
    //std::cout << "start cell iteration " <<std::endl;
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
        
        double tau_dyn = taumEle[cellnum-1];
        
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
        
        
        //std::cout << "cc " <<cc<<std::endl;
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
        */
        /*
        if (cellnum ==214)
        {
            std::cout << "cellnum  "<<cellnum<<std::endl;
        std::cout << "cellnum  " << cellnum<<" n_q_points "<<n_q_points<<std::endl;
        
        std::cout << "coordinate  " << coordx[0]<<" "<< coordx[1]<<" "<< coordx[2]<<" "<< coordx[3]<<std::endl;
        std::cout << "coordinate  " << coordy[0]<<" "<< coordy[1]<<" "<< coordy[2]<<" "<< coordy[3]<<std::endl;
        }
        */
        // loop through intergral points
        for (unsigned int q=0; q<n_q_points; ++q)
        {
            //gij[0] = gijG[cellnum-1][q][0];
            //gij[1] = gijG[cellnum-1][q][1];
            //gij[2] = gijG[cellnum-1][q][2];
            //gij[3] = gijG[cellnum-1][q][3];
            
            getGij(fe_values.inverse_jacobian(q),gij);
            rm[q] = present_vt_values[q] +
            present_velocity_gradients[q]*present_velocity_values[q]+present_pressure_gradients[q];//-sourceTerm[q];
            
            rc = present_velocity_gradients[q][0][0]+present_velocity_gradients[q][1][1];
            
             getTaum(tau_dyn,gij,dt,present_velocity_values[q][0],present_velocity_values[q][1],tauM,tauC);
            
            velocitybar[q] = present_velocity_values[q]-tauM*rm[q];

            double cJx = cc*fe_values.JxW(q);
            double Jx = fe_values.JxW(q);
            
            rhs1[q]=-present_vt_values[q]
                 -present_velocity_gradients[q]*velocitybar[q];
             
            /*
            if (cellnum <=1)
            {
                std::cout << "cellnum  "<<cellnum<<std::endl;
            //std::cout << "present_vt_values  " << present_vt_values[q]<<" present_velocity_values "<< present_velocity_values[q]<<" present_velocity_gradients "<< present_velocity_gradients[q]<<"present_pressure_values  "<< present_pressure_values[q]<<std::endl;
            std::cout << "gij  " << 4*gij[0]<<" "<< 4*gij[1]<<" "<< 4*gij[2]<<" "<< 4*gij[3]<<std::endl;
            std::cout << "rm  " << rm[q]<<" rc "<< rc<<" tauM "<< tauM<<" tauC "<< tauC<<std::endl;
            std::cout << "velocitybar  " << velocitybar[q]<<" Jx "<< Jx<<std::endl;
            }
            */
            double present_velocity_divergence =  trace(present_velocity_gradients[q]);
            // shape functions
            for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
                
                //div_phi_u[k]  = div_SF_u[cellnum-1][q][k];
                //grad_phi_u[k][0][0] = grad_SF_u[cellnum-1][q][k][0];
                //grad_phi_u[k][0][1] = grad_SF_u[cellnum-1][q][k][1];
                //grad_phi_u[k][1][0] = grad_SF_u[cellnum-1][q][k][2];
                //grad_phi_u[k][1][1] = grad_SF_u[cellnum-1][q][k][3];
                
                //phi_u[k][0] = SF_u[cellnum-1][q][k][0];
                //phi_u[k][1] = SF_u[cellnum-1][q][k][1];
                //phi_p[k]  = SF_p[cellnum-1][q][k];
                //grad_phi_p[k][0] =grad_SF_p[cellnum-1][q][k][0];
                //grad_phi_p[k][1] =grad_SF_p[cellnum-1][q][k][1];
                
                div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                phi_u[k]      =  fe_values[velocities].value(k, q);
                phi_p[k]      =  fe_values[pressure].value(k, q);
                grad_phi_p[k] =  fe_values[pressure].gradient(k, q);
                
                
                lhs1[k] = grad_phi_u[k]*velocitybar[q];
                //+present_velocity_gradients[q]*phi_u[k];
                lhs2[k] = grad_phi_u[k]*present_velocity_values[q];
                lhs3[k] = tauM*lhs2[k];
                lhs4[k] = tauC*div_phi_u[k];
                lhs5[k] = tauM*grad_phi_p[k];
                
                /*
                if (cellnum ==214)
                {
                    std::cout << "cellnum  "<<cellnum<<std::endl;
                std::cout << "div_phi_u  " << div_phi_u[k]<<" grad_phi_u "<< grad_phi_u[k]<<" phi_u "<< phi_u[k]<<"phi_p  "<< phi_p[k]<<std::endl;
                std::cout << "grad_phi_p  " << grad_phi_p[k]<<std::endl;
                std::cout << "lhs1  " << lhs1[k]<<" lhs2 "<< lhs2[k]<<" lhs3 "<< lhs3[k]<<" lhs4 "<< lhs4[k]<<std::endl;
                std::cout << "lhs5  " << lhs5[k]<<std::endl;
                }
                 */
            }
            // i ,j loop
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                local_rhs(i)+=
                    (
                     rhs1[q]*phi_u[i]
                    -viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                    + present_pressure_values[q]*div_phi_u[i]
                     - rc*tauC*div_phi_u[i]
                    - present_velocity_divergence*phi_p[i]
                    - lhs3[i]*rm[q]
                    -lhs5[i]*rm[q]
                    )
                    * Jx;
 
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {

                        lhs[i][j]+= (lhs1[j]*phi_u[i]
                                            +lhs3[j]*lhs2[i]
                                            +lhs4[j]*div_phi_u[i]
                                            +lhs5[j]*grad_phi_p[i])
                                            *cJx;
                        /*
                        local_matrix(i,j) +=
                        alpham*phi_u[j]*phi_u[i]*Jx // time
                        +viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])*cJx
                        - (div_phi_u[i]*phi_p[j] - phi_p[i]*div_phi_u[j])*cJx
                        +(lhs1[j]*phi_u[i] //convection
                        +lhs3[j]*lhs2[i] // stabilization
                        +lhs4[j]*div_phi_u[i]
                        +lhs5[j]*grad_phi_p[i]
                          )
                        *cJx;
                        */
                    }
            }
        }
        
        for ( unsigned int k = 0; k < dofs_per_cell; ++k)
            for ( unsigned int l = 0; l < dofs_per_cell; ++l)
                local_matrix(k,l) = lhs[k][l]+diffusive_k[cellnum-1][k][l];
        
        /*
        if (cellnum <=1)
        {
        std::cout << "local_rhs  "<<std::endl;;
         for ( unsigned int k = 0; k < dofs_per_cell; ++k)
         {
         std::cout  << local_rhs(k) <<std::endl;
         
         }
         
         std::cout << "local_matrix  "<<cellnum<<std::endl;
         
        //std::cout <<"system_matrix "<< system_matrix(0,0)<<" ";
        //std::cout<<"cellnum "<<cellnum<<std::endl;
        
         for ( unsigned int k = 0; k < dofs_per_cell; ++k)
         {
             //local_rhs(k) = rhs[k];
             for ( unsigned int l = 0; l < dofs_per_cell; ++l)
             {
                 //std::cout <<" "<< lhstest(k,l)<<" ";
                 std::cout <<" "<< local_matrix(k,l)<<" ";
             }
             std::cout<<std::endl;
         }        
        }
        */
        cell->get_dof_indices(local_dof_indices);
        
        //nonzero_constraints.clear();
        const AffineConstraints<double> &constraints_used = zero_constraints;
        
        constraints_used.distribute_local_to_global(local_matrix,
                                                    local_rhs,
                                                    local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);
        //if (cellnum <=4)
          //  std::cout <<"system_matrix "<< system_matrix(0,0)<<" ";
    }
    /*
    std::cout <<"system_matrix final "<< system_matrix(0,0)<<" ";
    std::cout << "system_rhs  ";
    
    for ( unsigned int k = 0; k < 11760; ++k)
    {
        if (system_rhs(k)>0){
            std::cout  << system_rhs(k) <<std::endl;
        }
    }
    */
        //std::cout << "cellnum bc " << cellnum <<std::endl;
        
        /* Neumann boundary condition
        First we have to find out whether the intersection of the
        faces of this cell with the newmann boundary part is nonzero.
         */
    
    //add_flux(npaf_solution,npam_solution_time_derivative);
    
    
}
void NavierStokes::getTaum(double tau_dyn, double* gij,double dt, double uele,double vele, double &tauM, double &tauC)
{
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
    tauC = (1)/(3*tauMtilta*(gij0+gij3));
    
    if (dynamicFlag)
        tauMtilta = tau_dyn;
        tauM = 1/sqrt(taut+(1/tauMtilta)*(1/tauMtilta));
    if (tauM>0)
        tauC = (1)/(3*tauMtilta*(gij0+gij3));
    else
        tauC=0;
    
    //double tauM =0;
    //std::cout << "degree  " << degree<<std::endl;
    //std::cout << "A11  " << A11<<" A22 "<<A22<<std::endl;
    //std::cout << "viscosity  " << viscosity<< " tau1sqinv  " << tau1sqinv<<std::endl;
    //std::cout << "tau2sqinv  " << tau2sqinv<<std::endl;
    //std::cout << "tauMtilta  " << tauMtilta<<std::endl;
    //std::cout << "gij0  " << gij0<<" gij3"<<gij3<<std::endl;
}
void NavierStokes::linear_solve(BlockVector<double> &newton_update)
{
    //call linear solver
    
    //system_matrix,system_rhs,newton_update
    
    //const ConstraintMatrix &constraints_used = zero_constraints;

    newton_update=system_rhs;
    SparseDirectUMFPACK directsolver;
    directsolver.initialize(system_matrix);
    directsolver.solve(newton_update);
    
    zero_constraints.distribute(newton_update);
    
}
void NavierStokes::update_flow()
{
    /*
    // update ut, pt
    ut_np1 = ut_np1+dut;
    pt_np1 = pt_np1+alphaf*gamma*dpt;
    
    // use ODE1 to update u_np1, p_np1
    u_np1 = u_np1+gamma*dt*dut;
    p_np1 = p_np1+alphaf*gamma*dt*dpt;
    */
    BlockVector<double>       tmp;
    tmp = newton_update;
    
    np1_solution_time_derivative.block(0) += newton_update.block(0);
    
    tmp.block(1) *= alphaf*gamma;
    np1_solution_time_derivative.block(1) += tmp.block(1);
        
    tmp.block(0) *= gamma*dt;
    np1_solution.block(0) += tmp.block(0);
    
    tmp.block(1) *= dt;
    np1_solution.block(1) += tmp.block(1);
}

void NavierStokes::getResidual(double &current_res)
{

    current_res = system_rhs.l2_norm();
    
    std::cout << "******************************" << std::endl;
    std::cout << " The residual of this guess is " << current_res << std::endl;
    
}

void NavierStokes::updator()
{
    
    n_solution_time_derivative.block(0) = np1_solution_time_derivative.block(0);
    n_solution.block(0) =np1_solution.block(0);
    
    n_solution_time_derivative.block(1) = np1_solution_time_derivative.block(1);
    n_solution.block(1) =np1_solution.block(1);
    
}
void NavierStokes::getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij)
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
void NavierStokes:: pre_compute(double *** &diffusive_k,
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
    //std::cout << "n_q_points  " << n_q_points<<std::endl;
    
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
        lhstest = 0;
        fe_values.reinit(cell);
        for (unsigned int q=0; q<n_q_points; ++q)
        {
            getGij(fe_values.inverse_jacobian(q),gij);
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
                diffusive_k[cellnum][k][l] = lhstest(k,l);
                //std::cout <<" diff"<< diffusive_k[cellnum-1][k][l]<<" ";
            }
            //std::cout<<std::endl;
            
        }
        cellnum++;
    }
    
    
}

    void NavierStokes::initialize(double*** &diffusive_k,double*** &div_SF_u,
                                  double**** &grad_SF_u,double**** &SF_u,
                                  double*** &SF_p,double**** &grad_SF_p,
                                  double*** &gijG)
    {
        
        //n_solution.reinit(dofs_per_block);
        //n_solution_time_derivative.reinit(dofs_per_block);
        //np1_solution.reinit(dofs_per_block);
        //np1_solution_time_derivative.reinit(dofs_per_block);
        
        QGauss<2> quadrature_formula(degree*2);
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points    = quadrature_formula.size();
        std::cout << "dofs_per_cell " << dofs_per_cell <<std::endl;
        
        int ele_num = triangulation.n_active_cells();
        for( int i=0; i<ele_num; i++)
        {
            diffusive_k[i] = new double*[dofs_per_cell];
            for (unsigned int j = 0; j<dofs_per_cell; j++)
            {
                diffusive_k[i][j] = new double[dofs_per_cell];
            }
        }
        
        for (int i = 0; i<ele_num; i++)
        {
            for (unsigned int j = 0; j<dofs_per_cell; j++)
            {
                for (unsigned int k = 0; k<dofs_per_cell; k++)
                {
                    diffusive_k[i][j][k]=0;
                }
            }
            
        }
        
        int EleNum = ele_num;
        for( int i=0; i<EleNum; i++)
        {
            div_SF_u[i] = new double*[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                div_SF_u[i][j] = new double[dofs_per_cell];
            }
        }
        
        for( int i=0; i<EleNum; i++)
        {
            grad_SF_u[i] = new double**[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                grad_SF_u[i][j] = new double*[dofs_per_cell];
                for (unsigned int k = 0; k<dofs_per_cell; k++)
                {
                    grad_SF_u[i][j][k] = new double[4];
                }
            }
        }
        
        for( int i=0; i<EleNum; i++)
        {
            SF_u[i] = new double**[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                SF_u[i][j] = new double*[dofs_per_cell];
                for (unsigned int k = 0; k<dofs_per_cell; k++)
                {
                    SF_u[i][j][k] = new double[2];
                }
            }
        }
        
        for( int i=0; i<EleNum; i++)
        {
            SF_p[i] = new double*[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                SF_p[i][j] = new double[dofs_per_cell];
            }
        }
        
        for( int i=0; i<EleNum; i++)
        {
            grad_SF_p[i] = new double**[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                grad_SF_p[i][j] = new double*[dofs_per_cell];
                for (unsigned int k = 0; k<dofs_per_cell; k++)
                {
                    grad_SF_p[i][j][k] = new double[2];
                }
            }
        }
        
        
        for( int i=0; i<EleNum; i++)
        {
            gijG[i] = new double*[n_q_points];
            for (unsigned int j = 0; j<n_q_points; j++)
            {
                gijG[i][j] = new double[4];
            }
        }

        
    }
    
    void NavierStokes::applyBC(double t)
    {
        int dim = 2;
        
        FEValuesExtractors::Vector velocities(0);
        FEValuesExtractors::Scalar pressure (dim);
        FEValuesExtractors::Scalar uy (1);
        nonzero_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
        
        ComponentMask mask2 = fe.component_mask(velocities);
        mask2.set(0,false);
        mask2.set(1,true);
        //std::cout <<"component_mask " <<mask2<<std::endl;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 6,
                                                 Dirichlet_BC(6,t),
                                                 nonzero_constraints,
                                                 fe.component_mask(velocities));
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 50,
                                                 Dirichlet_BC(50,t),
                                                 nonzero_constraints,
                                                 fe.component_mask(velocities));
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 24,
                                                 Dirichlet_BC(24,t),
                                                 nonzero_constraints,
                                                 mask2);
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 23,
                                                 Dirichlet_BC(23,t),
                                                 nonzero_constraints,
                                                 mask2);
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 14,
                                                 Dirichlet_BC(14,t),
                                                 nonzero_constraints,
                                                 fe.component_mask(pressure));
        //fe.component_mask(pressure));
        nonzero_constraints.close();
        
        
        zero_constraints.clear();
        
        DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 6,
                                                 ZeroFunction<2>(dim+1),
                                                 zero_constraints,
                                                 fe.component_mask(velocities));
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 50,
                                                 ZeroFunction<2>(dim+1),
                                                 zero_constraints,
                                                 fe.component_mask(velocities));
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 24,
                                                 ZeroFunction<2>(dim+1),
                                                 zero_constraints,
                                                 mask2);
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 23,
                                                 ZeroFunction<2>(dim+1),
                                                 zero_constraints,
                                                 mask2);
        
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 14,
                                                 ZeroFunction<2>(dim+1),
                                                 zero_constraints,
                                                 fe.component_mask(pressure));
        //fe.component_mask(velocities));
        zero_constraints.close();
        
        nonzero_constraints.distribute(n_solution);
        zero_constraints.distribute(n_solution_time_derivative);
        //std::cout << "nonzero_constraints lines in applyBC"<<std::endl;
        //nonzero_constraints.print(std::cout);
    }

    
    void NavierStokes::startingTimeloop(int num_start)
    
    {
        using namespace std;
        
        if (num_start!=0)              //readin starting points
        {
            ostringstream filename;
            filename << "./solution/"
            << 100
            << "-un-"
            << Utilities::int_to_string (num_start, 1)
            << ".txt";
            
            FILE * pFile;
            pFile = fopen(filename.str().c_str(), "r+");
            
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
            
            ostringstream filename1;
            filename1 << "./solution/"
            << 100
            << "-pn-"
            << Utilities::int_to_string (num_start, 1)
            << ".txt";
            
            FILE * pFile1;
            pFile1 = fopen(filename1.str().c_str(), "r+");
            
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
            
            ostringstream filename2;
            filename2 << "./solution/"
            << 100
            << "-ut_n-"
            << Utilities::int_to_string (num_start, 1)
            << ".txt";
            
            FILE * pFile2;
            pFile2 = fopen(filename2.str().c_str(), "r+");
            
            Vector<double> ut_n = n_solution_time_derivative.block(0);
            
            //cout<<" ut_n"<<ut_n.size()<<endl;
            
            double solution_ut;
            for(size_t i=0; i<ut_n.size(); ++i)
            {
                int ret = fscanf(pFile2, "%le", &solution_ut);
                ut_n[i] = solution_ut;
                //cout<<" ut_n[i] "<<ut_n[i]<<endl;
                if(ret == EOF) {
                    break;
                }
            }
            fclose(pFile2);
            
            n_solution.block(0) = u_n;
            n_solution.block(1) = p_n;
            n_solution_time_derivative.block(0) = ut_n;
        }
            
    }

    
    void NavierStokes::run()
    {
        setup_meshDOF();
        timeloop();
        //postprocess();
    }
     
}



