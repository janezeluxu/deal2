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
#include "../include/user_input.h"
#include "../include/generalizeAlpha.h"
#include "../include/post.h"
#include "../include/refine_mesh.h"
namespace incompressible
{
    using namespace dealii;
    
    NavierStokes::NavierStokes()
    :
    triangulation(Triangulation<2>::maximum_smoothing),
    dof_handler(triangulation)
    {
        for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
        {
            FE_Q<2> u(degree);
            FE_Q<2> p(degree);
            FESystem<2> fe(u,2, p,1);
            
            fe_collection.push_back (fe);
        }
    }

    
    void NavierStokes::setup_meshDOF()
    {
        //createmesh(triangulation);
        setDof();
        //initialize(dof_handler);
    }
    
    void NavierStokes::setDof()
    {
         //set up degree of freedoms and Dirichlet BCs using boundary ID
        
        dof_handler.distribute_dofs(fe_collection);
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
        
        for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
            face_quadrature_collection.push_back (QGauss<1>(degree*2));
        
    }
    
    void NavierStokes::applyBC()
    {
        FEValuesExtractors::Vector velocities(0);
        FEValuesExtractors::Scalar pressure (dim);
        FEValuesExtractors::Scalar uy (1);
        nonzero_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
        
        for (int i = 0; i < nBC_com; i++)
        {
            //std::cout << component_mask[i] << "\n";
            if (component_mask[i] == 0) // u
            {
                ComponentMask mask01 = fe_collection.component_mask(velocities);
                mask01.set(1,false);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                Dirichlet_BC(component_geotag[i]),
                nonzero_constraints,
                mask01);
            }
            else if (component_mask[i] == 1)// v
            {
                ComponentMask mask01 = fe_collection.component_mask(velocities);
                mask01.set(0,false);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                Dirichlet_BC(component_geotag[i]),
                nonzero_constraints,
                mask01);
            }
            else if (component_mask[i] == 2)// p
            {
                //ComponentMask mask01 = fe_collection.component_mask(pressure);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                Dirichlet_BC(component_geotag[i]),
                nonzero_constraints,
                fe_collection.component_mask(pressure));
            }
            else if (component_mask[i] == 3) //uv
            {
                //ComponentMask mask01 = fe_collection.component_mask(velocities);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                Dirichlet_BC(component_geotag[i]),
                nonzero_constraints,
                fe_collection.component_mask(velocities));
            }
            else if (component_mask[i] == 4) //uvp
            {
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                Dirichlet_BC(component_geotag[i]),
                nonzero_constraints);
            }
        }
        nonzero_constraints.close();
        
        
        zero_constraints.clear();
        
        DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
        
        for (int i = 0; i < nBC_com; i++)
        {
            //std::cout << component_mask[i] << "\n";
            if (component_mask[i] == 0) // u
            {
                ComponentMask mask01 = fe_collection.component_mask(velocities);
                mask01.set(1,false);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                ZeroFunction<2>(dim+1),
                zero_constraints,
                mask01);
            }
            else if (component_mask[i] == 1)// v
            {
                ComponentMask mask01 = fe_collection.component_mask(velocities);
                mask01.set(0,false);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                ZeroFunction<2>(dim+1),
                zero_constraints,
                mask01);
            }
            else if (component_mask[i] == 2)// p
            {
                //ComponentMask mask01 = fe_collection.component_mask(pressure);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                ZeroFunction<2>(dim+1),
                zero_constraints,
                fe_collection.component_mask(pressure));
            }
            else if (component_mask[i] == 3) //uv
            {
                //ComponentMask mask01 = fe_collection.component_mask(velocities);
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                ZeroFunction<2>(dim+1),
                zero_constraints,
                fe_collection.component_mask(velocities));
            }
            else if (component_mask[i] == 4) //uvp
            {
                VectorTools::interpolate_boundary_values(dof_handler,
                component_geotag[i],
                ZeroFunction<2>(dim+1),
                zero_constraints);
            }
        }
        zero_constraints.close();
    }
    void NavierStokes::initialize()
    {
        
        n_solution.reinit(dofs_per_block);
    }
    
    void NavierStokes::loop()
    {
        double viscosity = 1/Re;
        //extern unsigned int nstp ;
        std::cout << "viscosity in NavierStoke " << viscosity<<std::endl;
        
        
        setup_meshDOF();
        initialize();
        applyBC();
        
        generalizeAlpha multicorrector(viscosity);
        
        for(int cycle=1; cycle<=max_cycle; ++cycle)
        {
            std::cout << "cycle " << cycle<<std::endl;
            multicorrector.initialize(dof_handler,nonzero_constraints,dofs_per_block);
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
            //Newton's iteration
            multicorrector.corrector(maxNewtonIter,tolerance,dof_handler,fe_collection,zero_constraints,nonzero_constraints,n_solution);
            
            /*
             Vector<double> du_n = n_solution.block(0);
             Vector<double> dp_n = n_solution.block(1);
             std::cout << "n_solution check " <<std::endl;
             for(std::size_t i=0; i<du_n.size(); ++i)
             {
             std::cout<<du_n[i]<<std::endl;
             }
             
             for(std::size_t i=0; i<dp_n.size(); ++i)
             {
             std::cout<<dp_n[i]<<std::endl;
             }
            */
            
            Post processResult;
            processResult.savesolution(n_solution,cycle);
            processResult.savetovtk(dof_handler,n_solution);
            auto vertices = triangulation.get_vertices();
            int nEdge = triangulation.n_lines();
            processResult.savemesh(cycle,vertices.size(),nEdge,dof_handler,fe_collection, nonzero_constraints);
            
            std::cout<<"save mesh finished"<<std::endl;
            Vector<float> true_error_per_cell (triangulation.n_active_cells());
            if (calcError)
            {
                processResult.getTrueError(fe_collection,triangulation,dof_handler,n_solution,true_error_per_cell);
            }
            //hp
            std::cout<<"adapt mesh"<<std::endl;
            mesh_adapt(cycle,true_error_per_cell);
            std::cout<<"adapt mesh finished"<<std::endl;
            
            
        }
        
        
        
    }

    void NavierStokes::mesh_adapt(int cycle,Vector<float> &true_error_per_cell)
    {
        refine_mesh hp_refinement(triangulation,dof_handler,n_solution,fe_collection);
        SolutionTransfer<2, BlockVector<double>,hp::DoFHandler<2>  > solution_transfer(dof_handler);
        Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
        hp_refinement.perfrom_refinement(face_quadrature_collection, cycle,solution_transfer,estimated_error_per_cell,true_error_per_cell);
        
        setup_meshDOF();
        BlockVector<double> tmp (dofs_per_block);
        //n_solution.reinit(dof_handler.n_dofs());
        solution_transfer.interpolate(n_solution, tmp);
        n_solution = tmp;
        nonzero_constraints.distribute(tmp);
        
        //initialize();
        applyBC();
        //n_solution = tmp;
    }
    
    void NavierStokes::run()
    {
        createmesh(triangulation);
        loop();
        //postprocess();
    }
     
}



