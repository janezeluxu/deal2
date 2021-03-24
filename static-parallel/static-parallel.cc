/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Zelu Xu, RPI, 2020
 */


// @sect3{Include files}
//
// Most of the include files we need for this program have already been
// discussed in previous programs. In particular, all of the following should
// already be familiar friends:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

// uncomment the following #define if you have PETSc and Trilinos installed
// and you prefer using Trilinos in this example:
// #define FORCE_USE_OF_TRILINOS

// This will either import PETSc or TrilinosWrappers into the namespace
// LA. Note that we are defining the macro USE_PETSC_LA so that we can detect
// if we are using PETSc (see solve() for an example where this is necessary)
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/conditional_ostream.h>

#include <fstream>
#include <iostream>
#include <sstream>

namespace incompressible
{
    /*
     This is a paralleled version of the same SUPG based incompressible Navier Stoke solver. and we solve the same flow over cylinder problem.
     The linear solver is an iterative solver, we first form the auxiliary probelem (see equation 31-33 in doc). The preconditioner for the auxiliary form, I only used the diagonal terms for K^-1
     */
  using namespace dealii;
  using namespace std;
  template <int dim>
  class NavierStokes
  {
  public:
    NavierStokes();

    void run();

  private:
    void createmesh(int n_refines);
    void restart_time(int num_start);
    void setup_system();
    void timeloop(int num_start, int nstp);
    void predictor();
    void corrector();
    void updator();
    void iterPC();
    void assemble_system();
    void update_flow();
    void iterative_solve();
    //void direct_solve();
    void getResidual(double &current_res);
    void applyBC(double t);
    void getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij);
    void getTaum(double* gij,double dt, double uele,double vele,
                 double &tauM, double &tauC);
      void getdata(std::string datatype, int num_start, double* datanum);
    //void solve();
    //void output_results(const unsigned int cycle) const;
    void savetovtk(int istp) const;
    void save_solution(int istp);
    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    FESystem<dim>      fe;
    DoFHandler<dim> dof_handler;
    std::vector<types::global_dof_index> dofs_per_block;
      
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;
      
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;
      
    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockVector       system_rhs;
      
    
    LA::MPI::BlockVector       n_solution_time_derivative;
    LA::MPI::BlockVector       n_solution;
    LA::MPI::BlockVector       np1_solution_time_derivative;
    LA::MPI::BlockVector       np1_solution;
    LA::MPI::BlockVector       npaf_solution;
    LA::MPI::BlockVector       npam_solution_time_derivative;
     
    LA::MPI::BlockVector       locally_relevant_npaf;
    LA::MPI::BlockVector       locally_relevant_npam;
    LA::MPI::BlockVector       locally_relevant_solution_output;
    LA::MPI::BlockVector       newton_update;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
      
    // parameters
    double Re = 100;
    double viscosity = 1/Re;
    double dt = 0.03125;
    double ri = 0.25;
    double alpham = 0.5*(3- ri)/(1+ri);
    double alphaf = 1/(1+ri);
    double gamma = 0.5+alpham-alphaf;
    int maxNewtonIter = 10;
    double tolerance = 1e-6;
  };

  template <int dim>
  NavierStokes<dim>::NavierStokes()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(FE_Q<dim>(1), 2, FE_Q<dim>(1), 1)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {}

  template <int dim>
  void NavierStokes<dim>::createmesh(int n_refines)
  {
      GridIn<2> grid_in;
      grid_in.attach_triangulation(triangulation);
      {
          std::string   filename = "cylinderquadcc.inp";
          std::ifstream file(filename.c_str());
          Assert(file, ExcFileNotOpen(filename.c_str()));
          grid_in.read_ucd(file);
      }
      //int n_refines = 3;
      triangulation.refine_global(n_refines);
      std::vector<types::boundary_id>           boundary_ids;
      boundary_ids = triangulation.get_boundary_ids();
      
      std::ofstream out("grid-flowovercylinder.eps");
      GridOut       grid_out;
      grid_out.write_eps(triangulation, out);
      
      //GridOutFlags::Gnuplot gnuplot_flags(false, 60);
      //grid_out.set_flags(gnuplot_flags);
      //grid_out.write_eps(triangulation, out, &mapping);
      std::cout << "Grid written to grid-flowoversquare.eps" << std::endl;
  }

  template <int dim>
  void NavierStokes<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup");

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
      
    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, dof_u);
    owned_partitioning[1] =
    dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);
      
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, dof_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(dof_u, dof_u + dof_p);
      
    // with ghost cell, need in function evaluation
    locally_relevant_npaf.reinit(owned_partitioning,
                                     relevant_partitioning,
                                     mpi_communicator);
    locally_relevant_npam.reinit(locally_relevant_npaf);
      
    locally_relevant_solution_output.reinit(locally_relevant_npaf);
    // no ghost cell
    n_solution_time_derivative.reinit(owned_partitioning,
    mpi_communicator);
    np1_solution_time_derivative.reinit(n_solution_time_derivative);
    n_solution.reinit(n_solution_time_derivative);
    np1_solution.reinit(n_solution_time_derivative);
    newton_update.reinit(n_solution_time_derivative);
      
    system_rhs.reinit(owned_partitioning, mpi_communicator);

    applyBC(0);
      
    BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
    SparsityTools::distribute_sparsity_pattern(
    dsp,
    dof_handler.locally_owned_dofs_per_processor(),
    mpi_communicator,
    locally_relevant_dofs);
      
    system_matrix.reinit(owned_partitioning,
    owned_partitioning,
    dsp,
    mpi_communicator);
  }

template <int dim>
void NavierStokes<dim>::timeloop(int num_start,int nstp)
{   //double viscosity = 1/Re;
    //std::cout << "viscosity in NavierStoke " << viscosity<<std::endl;
    //int Grid_size = triangulation.n_active_cells();
    restart_time(num_start);
    
    for (int istp = num_start; istp < nstp+num_start; ++istp)
    {
        std::cout << "-------start timestep-------" <<istp<< std::endl;
        double t = istp*dt;
        applyBC(t);
    
       //std::cout << "-------predictor-------" <<istp<< std::endl;
       predictor();

       //Newton's iteration
       corrector();
       //std::cout << "-------updator-------" <<istp<< std::endl;
       updator();

        
        if (istp % 1 == 0)
        {
        //Post processResult(degree);
        locally_relevant_solution_output = n_solution;
        //std::cout << "-------save vtk-------" << std::endl;
        savetovtk(istp);
        //std::cout << "-------save solution-------" << std::endl;
        save_solution(istp+1);
        }
        
    }
}

template <int dim>
void NavierStokes<dim>::predictor()
{
    //ut_np1[i] = ((gamma-1)/gamma)*ut_n[i];
    np1_solution_time_derivative = n_solution_time_derivative;
    np1_solution_time_derivative.block(0) *= (gamma-1)/gamma;
    np1_solution_time_derivative.block(1) =0;
    
    LA::MPI::BlockVector       tmp;
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

template <int dim>
void NavierStokes<dim>::corrector()
{
    
    double current_res=0;
    
    for(int i=0; i<maxNewtonIter; ++i)
    {
        //bool isfirst = FALSE;
        //std::cout << "----start iteration Number ----" <<i<< std::endl;
    
        iterPC();

        //std::cout << "----get newton update ----" << std::endl;
        assemble_system();
        //direct_solve();
        iterative_solve();

        //std::cout << "-------now update flow-------" << std::endl;
        update_flow();

        //std::cout << "-------now get residual-------" << std::endl;
        getResidual(current_res);

        if (current_res < tolerance)
            break;
    }

}

template <int dim>
void NavierStokes<dim>::updator()
{
    /* update flow
       u_n = u_np1;
       ut_n = ut_np1;
       p_n = p_np1;
       pt_n = pt_np1;
    */
    n_solution_time_derivative.block(0) = np1_solution_time_derivative.block(0);
    n_solution.block(0) =np1_solution.block(0);
    
    n_solution_time_derivative.block(1) = np1_solution_time_derivative.block(1);
    n_solution.block(1) =np1_solution.block(1);
    
}

template <int dim>
void NavierStokes<dim>::iterPC()
{
    //u_npaf[i] = u_n[i]+alphaf*(u_np1[i]-u_n[i]);
    //ut_npam[i] = ut_n[i]+alpham*(ut_np1[i]-ut_n[i]);
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
    
    
    locally_relevant_npam = npam_solution_time_derivative;
    locally_relevant_npaf = npaf_solution;
}

template <int dim>
void NavierStokes<dim>::assemble_system()
{
    system_matrix = 0;
    system_rhs = 0;
  TimerOutput::Scope t(computing_timer, "assembly");

  const QGauss<dim> quadrature_formula(fe.degree*2);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values |
                          update_inverse_jacobians);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(2);
    
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    
  std::vector<Tensor<1, 2> > present_velocity_values(n_q_points);
  std::vector<Tensor<1, 2> > present_vt_values(n_q_points);

  std::vector<Tensor<2, 2> > present_velocity_gradients(n_q_points);
  std::vector<double>       present_pressure_values(n_q_points);
  std::vector<Tensor<1, 2> > present_pressure_gradients(n_q_points);
    
  std::vector<Tensor<1, 2> > rm(n_q_points);
  std::vector<Tensor<1, 2> > velocitybar(n_q_points);
    
  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, 2> > phi_u(dofs_per_cell);
  std::vector<Tensor<2, 2> > grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);
  std::vector<Tensor<1, 2> > grad_phi_p(dofs_per_cell);
    
  std::vector<Tensor<1, 2> >        lhs1(dofs_per_cell);
  std::vector<Tensor<1, 2> >        lhs2(dofs_per_cell);
  std::vector<Tensor<1, 2> >        lhs3(dofs_per_cell);
  std::vector<double >        lhs4(dofs_per_cell);
  std::vector<Tensor<1, 2> >        lhs5(dofs_per_cell);
  std::vector<Tensor<1, 2> >        rhs1(n_q_points);
    
  int cellnum = 0;
  double* gij = new double[4];
  double rc = 0;
  double tauM = 0;
  double tauC = 0;

    
  double cc = (alphaf*gamma*dt);
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        int vertexNum = GeometryInfo<2>::vertices_per_cell;
        double* coordx = new double[vertexNum];
        double* coordy = new double[vertexNum];
        
        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
        {
            Point<2> &v = cell->vertex(i);
            coordx[i] =  v[0];
            coordy[i] = v[1];
        }

        cellnum++;
        cell_matrix = 0.;
        cell_rhs    = 0.;

        fe_values.reinit(cell);

        fe_values[velocities].get_function_values(locally_relevant_npaf,present_velocity_values);
         
        fe_values[velocities].get_function_gradients(locally_relevant_npaf, present_velocity_gradients);
          
        fe_values[pressure].get_function_values(locally_relevant_npaf,present_pressure_values);
        fe_values[pressure].get_function_gradients(locally_relevant_npaf,present_pressure_gradients);
        
        fe_values[velocities].get_function_values(locally_relevant_npam,present_vt_values);
          
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
 
            getGij(fe_values.inverse_jacobian(q),gij);
              
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
                
                div_phi_u[k]  =  fe_values[velocities].divergence (k, q);
                grad_phi_u[k] =  fe_values[velocities].gradient(k, q);
                phi_u[k]      =  fe_values[velocities].value(k, q);
                phi_p[k]      =  fe_values[pressure].value(k, q);
                grad_phi_p[k] =  fe_values[pressure].gradient(k, q);
                
                lhs1[k] = grad_phi_u[k]*velocitybar[q]+
                            present_velocity_gradients[q]*phi_u[k];
                lhs2[k] = grad_phi_u[k]*present_velocity_values[q];
                lhs3[k] = tauM*lhs2[k];
                lhs4[k] = tauC*div_phi_u[k];
                lhs5[k] = tauM*grad_phi_p[k];

                
            }
              
            for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                     
                    cell_rhs(i)+=
                        (rhs1[q]*phi_u[i]
                        - viscosity*scalar_product(present_velocity_gradients[q],grad_phi_u[i])
                        + (present_pressure_values[q]- rc*tauC)*div_phi_u[i]
                        - present_velocity_divergence*phi_p[i]
                        -( lhs3[i]+lhs5[i])*rm[q]
                        )
                        * Jx;
                    
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            cell_matrix(i, j) += alpham*phi_u[j]*phi_u[i]*Jx
                            +(viscosity*scalar_product(grad_phi_u[j], grad_phi_u[i])
                              - div_phi_u[i]*phi_p[j]
                              + phi_p[i]*div_phi_u[j])*cJx
                            +(lhs1[j]*phi_u[i]+lhs3[j]*lhs2[i]+lhs4[j]*div_phi_u[i]
                                                +lhs5[j]*grad_phi_p[i])
                                                *cJx;
                        }
                }// end i j
            }// end n_q_points


        cell->get_dof_indices(local_dof_indices);
        zero_constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
      }

  // Notice that the assembling above is just a local operation. So, to
  // form the "global" linear system, a synchronization between all
  // processors is needed. This could be done by invoking the function
  // compress(). See @ref GlossCompress "Compressing distributed objects"
  // for more information on what is compress() designed to do.
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

}

template <int dim>
void NavierStokes<dim>::update_flow()
{
    /*
    // update ut, pt
    ut_np1 = ut_np1+dut;
    pt_np1 = pt_np1+alphaf*gamma*dpt;
    
    // use ODE1 to update u_np1, p_np1
    u_np1 = u_np1+gamma*dt*dut;
    p_np1 = p_np1+alphaf*gamma*dt*dpt;
    */
    
    LA::MPI::BlockVector       tmp;
    tmp = newton_update;
    
    np1_solution_time_derivative.block(0) += newton_update.block(0);
    
    tmp.block(1) *= alphaf*gamma;
    np1_solution_time_derivative.block(1) += tmp.block(1);
        
    tmp.block(0) *= gamma*dt;
    np1_solution.block(0) += tmp.block(0);
    
    tmp.block(1) *= dt;
    np1_solution.block(1) += tmp.block(1);
}


template <int dim>
void NavierStokes<dim>::getResidual(double &current_res)
   {
       current_res = system_rhs.l2_norm();
       
       std::cout << "******************************" << std::endl;
       std::cout << " The residual of this guess is " << current_res << std::endl;
       
   }


template <int dim>
void NavierStokes<dim>::getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij)
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

template <int dim>
void NavierStokes<dim>::getTaum(double* gij,double dt, double uele,double vele,
                            double &tauM, double &tauC)
{

    double C1 = 1;
    double C2 = 9;//(3*degree*degree)*(3*degree*degree);

    int deg = fe.degree;
    double gij0 = 4*gij[0]*deg*deg;
    double gij1 = 4*gij[1]*deg*deg;
    double gij2 = 4*gij[2]*deg*deg;
    double gij3 = 4*gij[3]*deg*deg;
    
    double tau1sqinv = uele*(uele*gij0+vele*gij2)+vele*(uele*gij1+vele*gij3);
    double A11 =gij0*gij0+gij2*gij2;
    double A22 =gij1*gij1+gij3*gij3;
    
    double tau2sqinv = C2*viscosity*viscosity*(A11+A22);
    double taut = (2*C1/dt)*(2*C1/dt);
    tauM = 1/sqrt(taut+tau1sqinv+tau2sqinv);
    double tauMtilta = 1/sqrt(tau1sqinv+tau2sqinv);
    tauC = (1)/((3)*tauMtilta*(gij0+gij3));
}

template <int dim>
void NavierStokes<dim>::savetovtk(int istp) const
{
  std::vector<std::string> solution_names(dim, "velocity");
  solution_names.push_back("pressure");
  

  std::vector<DataComponentInterpretation::DataComponentInterpretation>data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(locally_relevant_solution_output,
  solution_names,
  DataOut<2>::type_dof_data,
  data_component_interpretation);
    
    //std::cout << "----stop3 ----" << std::endl;
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches();

  // The next step is to write this data to disk. We choose file names of
  // the form <code>solution-XX.PPPP.vtu</code> where <code>XX</code>
  // indicates the refinement cycle, <code>PPPP</code> refers to the
  // processor number (enough for up to 10,000 processors, though we hope
  // that nobody ever tries to generate this much data -- you would likely
  // overflow all file system quotas), and <code>.vtu</code> indicates the
  // XML-based Visualization Toolkit for Unstructured grids (VTU) file
  // format.
  const std::string filename =
    ("./solution-" + Utilities::int_to_string(istp, 2) + "." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));
  std::ofstream output(filename + ".vtu");
  data_out.write_vtu(output);
    
  // The last step is to write a "master record" that lists for the
  // visualization program the names of the various files that combined
  // represents the graphical data for the entire domain. The
  // DataOutBase::write_pvtu_record does this, and it needs a list of
  // filenames that we create first. Note that only one processor needs to
  // generate this file; we arbitrarily choose processor zero to take over
  // this job.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
        filenames.push_back("./solution-" + Utilities::int_to_string(istp, 2) +
                            "." + Utilities::int_to_string(i, 4) + ".vtu");

      std::ofstream master_output(
        "./solution-" + Utilities::int_to_string(istp, 2) + ".pvtu");
      data_out.write_pvtu_record(master_output, filenames);
    }
     
}

template <int dim>
void NavierStokes<dim>::save_solution(int istp)
{
    using namespace std;
    
    // print n_solution values to file
    const std::string filename =
    ("./solution/" + Utilities::int_to_string(istp, 2) + "un." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 2)+".txt");
    
    std::filebuf fb;
    fb.open (filename,std::ios::out);
    std::ostream os(&fb);
    
    LA::MPI::Vector du_n = n_solution.block(0);
    du_n.print(os,10);
    
    const std::string filenamep =
    ("./solution/" + Utilities::int_to_string(istp, 2) + "pn." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 2)+".txt");
    
    std::filebuf fbp;
    fbp.open (filenamep,std::ios::out);
    std::ostream osp(&fbp);
    
    LA::MPI::Vector dp_n = n_solution.block(1);
    dp_n.print(osp,10);
    
    // print n_solution_time_derivative values to file
    const std::string filenameut =
    ("./solution/" + Utilities::int_to_string(istp, 2) + "utn." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 2)+".txt");
    
    std::filebuf fbut;
    fbut.open (filenameut,std::ios::out);
    std::ostream osut(&fbut);
    
    LA::MPI::Vector dut_n = n_solution_time_derivative.block(0);
    dut_n.print(osut,10);
    
    const std::string filenamept =
    ("./solution/" + Utilities::int_to_string(istp, 2) + "ptn." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 2)+".txt");
    
    std::filebuf fbpt;
    fbpt.open (filenamept,std::ios::out);
    std::ostream ospt(&fbpt);
    
    LA::MPI::Vector dpt_n = n_solution_time_derivative.block(1);
    dpt_n.print(ospt,10);
}

template <int dim>
void NavierStokes<dim>::getdata(std::string datatype, int num_start, double* datanum)
{
    string filename;
    for (unsigned int i = 0;
         i < Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
    {
      filename = ("./solution/" + Utilities::int_to_string(num_start, 2) +
                          datatype + Utilities::int_to_string(i, 2) + ".txt");
      string line;
      ifstream myfile (filename);
      if (myfile.is_open())
      {
        getline (myfile,line); // file is only one line
        myfile.close();
      }
        std::string delimiter = "]";
        std::size_t pos = line.find(delimiter);
        std::string token = line.substr(0, pos);
        std::string s = line.substr(pos+2,line.length()-1);
        
        delimiter = " ";
        pos = token.find(delimiter);
        std::string token1 = token.substr(0, pos);
        std::string token2 = token.substr(pos);
        
        delimiter = "-";
        pos = token2.find(delimiter);
        std::string Range1 = token2.substr(0, pos);
        std::string Range2 = token2.substr(pos+1);
        
        std::pair< int, int > dataRange;
        dataRange.first = stoi(Range1);
        dataRange.second = stoi(Range2);
        
        string to;
        pos = 0;
        delimiter = " ";

        int count = dataRange.first;
        std::cout << " count " << count << std::endl;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            datanum[count] = stod(token);
            s.erase(0, pos + delimiter.length());
            count++;
        } // end while
    } // end for
    
}
template <int dim>
void NavierStokes<dim>::restart_time(int num_start)
  
  {
      //using namespace std;
      
      std::cout << "   Number of total DOF: " << n_solution.block(0).size()<<" "<<n_solution.block(1).size()
      << std::endl;
      
      if (num_start!=0)              //readin starting points
      {
          int nuSize = n_solution.block(0).size();
          double* data_un = new double[nuSize];
          double* data_utn = new double[nuSize];
          
          int npSize = n_solution.block(1).size();
          double* data_pn = new double[npSize];
        
            getdata("un.", num_start, data_un);
            getdata("pn.", num_start, data_pn);
            getdata("utn.", num_start, data_utn);
              
            for (int i = 0; i < 1125; ++i)
            {
                std::cout << data_pn[i]<< std::endl;
            }
          
          // distribute data to different processors
          LA::MPI::Vector u_n = n_solution.block(0);
          std::pair< int, int > Mrange =  u_n.local_range();
          for (int i = Mrange.first; i < Mrange.second; ++i)
                u_n[i] = data_un[i];
          u_n.compress(VectorOperation::insert);
          n_solution.block(0) = u_n;
          
          LA::MPI::Vector ut_n = n_solution_time_derivative.block(0);
          Mrange =  ut_n.local_range();
          for (int i = Mrange.first; i < Mrange.second; ++i)
                ut_n[i] = data_utn[i];
          ut_n.compress(VectorOperation::insert);
          n_solution_time_derivative.block(0) = ut_n;
          
          LA::MPI::Vector p_n = n_solution.block(1);
          Mrange =  p_n.local_range();
          for (int i = Mrange.first; i < Mrange.second; ++i)
                p_n[i] = data_pn[i];
          p_n.compress(VectorOperation::insert);
          n_solution.block(1) = p_n;
        
      }
          
  }
double Dirichlet_u(const dealii::Point<2> & p,int BC_ID,double t)
{
    if (BC_ID==6) //left
        return 1;
    else if (BC_ID==50) // surface
    {
        return 0;
    }
    return 0;
}

double Dirichlet_v(const dealii::Point<2> & p,int BC_ID,double t)
{
    if (BC_ID==6) //left
        return 0;
    else if (BC_ID==24) // top
        return 0;
    else if (BC_ID==23) //buttom
        return 0;
    else if (BC_ID==50) //surface
        return 0;
    return 0;
}

double Dirichlet_p(const dealii::Point<2> & p,int BC_ID,double t)
{
    if (BC_ID==14) //right
        return 0;
    return 0;
}


class Dirichlet_BC : public Function<2>
{
public:
    Dirichlet_BC(int bc_id,double time) : Function<2>(3)
    {
        boundaryID=bc_id;
        t = time;
    }
    
    int boundaryID;
    double t;
    
    virtual double value( const Point<2> & p,
                         const unsigned int component) const;
    virtual void vector_value( const Point<2> &p,
                              Vector<double> &  values) const;
};

double Dirichlet_BC::value(const Point<2> & p, const unsigned int component) const
{

    Assert(component < this->n_components,
           ExcIndexRange(component, 0, this->n_components));

    if (component == 0 )
        return Dirichlet_u(p,boundaryID,t);
    else if (component == 1 )
        return Dirichlet_v(p,boundaryID,t);
    else if (component == 2 )
        return Dirichlet_p(p,boundaryID,t);
    else
        std::cout<<"which component?? "<<std::endl;
    
    return 0;
}

void Dirichlet_BC::vector_value(const Point<2> &p,
                                  Vector<double> &  values) const
{
    for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = Dirichlet_BC::value(p, c);
}

 template <int dim>
 void NavierStokes<dim>::applyBC(double t)
 {
     //int dim = 2;
     
     FEValuesExtractors::Vector velocities(0);
     FEValuesExtractors::Scalar pressure (dim);
     FEValuesExtractors::Scalar uy (1);
     nonzero_constraints.clear();
     //nonzero_constraints.reinit(locally_relevant_dofs);
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
     //zero_constraints.reinit (locally_relevant_dofs);
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
     //std::cout << "applyBC finish" <<std::endl;
 }

template <class MatrixType>
class InverseMatrix : public Subscriptor
{
public:
    InverseMatrix(const MatrixType &m,
    const IndexSet           &locally_owned,
    const MPI_Comm           &mpi_communicator);
    //template <typename VectorType>
    void vmult(LA::MPI::Vector      &dst,
               const LA::MPI::Vector &src) const;
private:
    const SmartPointer<const MatrixType> matrix;
    const MPI_Comm *mpi_communicator;
    mutable LA::MPI::Vector tmp;
};

template <class MatrixType>
InverseMatrix<MatrixType>::InverseMatrix
(const MatrixType &m,
const IndexSet           &locally_owned,
const MPI_Comm           &mpi_communicator)
:
matrix (&m),
mpi_communicator (&mpi_communicator),
tmp(locally_owned, mpi_communicator)
{}

template <class MatrixType>
//template <typename VectorType>
void InverseMatrix<MatrixType>::vmult (LA::MPI::Vector        &dst,
                                       const LA::MPI::Vector  &src) const
{
    SolverControl solver_control (std::max<unsigned int>(src.size(), 200),
                                  1e-8*src.l2_norm());
    SolverFGMRES<LA::MPI::Vector> gmres(solver_control);
    //tmp = 0.;
    gmres.solve (*matrix, dst, src, PreconditionIdentity());
    //dst = tmp;
    
}

class SchurComplement : public Subscriptor
{
public:
    SchurComplement (const LA::MPI::BlockSparseMatrix            &A,
                     const InverseMatrix<LA::MPI::SparseMatrix > &Minv,
                     const IndexSet   &owned_vel0,
                     const IndexSet   &owned_vel1,
                     const MPI_Comm &mpi_communicator);
    //template <typename VectorType>
    void vmult (LA::MPI::Vector      &dst,
                const LA::MPI::Vector &src) const;
private:
    const SmartPointer<const LA::MPI::BlockSparseMatrix  > system_matrix;
    const SmartPointer<const InverseMatrix<LA::MPI::SparseMatrix > > m_inverse;
    //template <typename VectorType>
    mutable LA::MPI::Vector tmp1, tmp2, tmpm,tmpc;
};

//template <typename VectorType>
SchurComplement
::SchurComplement (const LA::MPI::BlockSparseMatrix           &A,
                   const InverseMatrix<LA::MPI::SparseMatrix > &Minv,
                   const IndexSet                            &owned_vel0,
                   const IndexSet                            &owned_vel1,
                   const MPI_Comm                            &mpi_communicator)
:
system_matrix (&A),
m_inverse (&Minv),
tmp1 (owned_vel0,mpi_communicator),
tmp2 (tmp1),
tmpm (owned_vel1,mpi_communicator),
tmpc (tmpm)

{}

//template <typename VectorType>
void SchurComplement::vmult (LA::MPI::Vector      &dst,
                             const LA::MPI::Vector &src) const
{
    system_matrix->block(0,1).vmult (tmp1, src);// multiply with the top right block: B
    m_inverse->vmult (tmp2, tmp1); // multiply with M^-1
    system_matrix->block(1,0).vmult (tmpm, tmp2); // multiply with the bottom left block: -B^T
    system_matrix->block(1,1).vmult (tmpc, src); // multiply with the bottom right block: C
    dst = tmpc;
    dst -= tmpm;
}
class ApproximateSchurComplement : public Subscriptor
{
public:
    ApproximateSchurComplement (const LA::MPI::BlockSparseMatrix &A,
                                const IndexSet   &owned_vel0,
                                const IndexSet   &owned_vel1,
                                const MPI_Comm &mpi_communicator);
    //template <typename VectorType>
    void vmult (LA::MPI::Vector       &dst,
                const LA::MPI::Vector &src) const;
private:
    const SmartPointer<const LA::MPI::BlockSparseMatrix> system_matrix;
    mutable LA::MPI::Vector tmp1, tmp2,tmpm,tmpc;
};

ApproximateSchurComplement::ApproximateSchurComplement
(const LA::MPI::BlockSparseMatrix &A,const IndexSet   &owned_vel0,
const IndexSet   &owned_vel1,
const MPI_Comm &mpi_communicator) :
system_matrix (&A),
tmp1 (owned_vel0,mpi_communicator),
tmp2 (tmp1),
tmpm (owned_vel1,mpi_communicator),
tmpc (tmpm)

{}

//template <typename VectorType>
void ApproximateSchurComplement::vmult
(LA::MPI::Vector       &dst,
 const LA::MPI::Vector  &src) const
{
    system_matrix->block(0,1).vmult (tmp1, src); // multiply with the top right block: G
    
    //multiply with block diag(K)^-1
    const auto &M = system_matrix->block(0, 0);
    std::pair< int, int > Mrange =  M.local_range();
    for (int i = Mrange.first; i < Mrange.second; ++i)
    tmp2(i) = tmp1(i) / M.diag_element(i);
    
    //multiply with bottom left block G^T
    system_matrix->block(1,0).vmult (tmpm, tmp2);

    system_matrix->block(1,1).vmult (tmpc, src); // multiply with the bottom right block: C

    dst = tmpc;
    dst -= tmpm;
     
    
}

template <int dim>
void NavierStokes<dim>::iterative_solve ()
{
    TimerOutput::Scope              t(computing_timer, "solve_iterative");
    InverseMatrix<LA::MPI::SparseMatrix > inverse_mass (system_matrix.block(0,0),owned_partitioning[0],mpi_communicator);
    LA::MPI::Vector tmp (owned_partitioning[0],mpi_communicator);
    {
        
        LA::MPI::Vector schur_rhs (owned_partitioning[1],mpi_communicator);
        inverse_mass.vmult (tmp, system_rhs.block(0));
        
        system_matrix.block(1,0).vmult (schur_rhs, tmp);
        
        schur_rhs -= system_rhs.block(1);
        
        schur_rhs *= -1;
        
        SolverControl  solver_control(system_matrix.m(),1e-4 * schur_rhs.l2_norm(),true);
        SolverFGMRES<LA::MPI::Vector> gmres(solver_control);
        
        
        ApproximateSchurComplement approximate_schur (system_matrix,owned_partitioning[0],owned_partitioning[1],
        mpi_communicator);
        
        InverseMatrix<ApproximateSchurComplement> preconditioner (approximate_schur,owned_partitioning[1],mpi_communicator);
        
        //PreconditionIdentity preconditioner1;
        
        SchurComplement schur_complement (system_matrix, inverse_mass,owned_partitioning[0],owned_partitioning[1],
        mpi_communicator);
        //std::cout << " gmres solve Schur complement"<< std::endl;
        gmres.solve(schur_complement, newton_update.block(1), schur_rhs, preconditioner);
        std::cout << solver_control.last_step()
        << " gmres Schur complement iterations to obtain convergence."
        << std::endl;
        
    }
    
    
    system_matrix.block(0,1).vmult (tmp, newton_update.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    inverse_mass.vmult (newton_update.block(0), tmp);
}

  template <int dim>
  void NavierStokes<dim>::run()
  {
    pcout << "Running with "
#ifdef USE_PETSC_LA
          << "PETSc"
#else
          << "Trilinos"
#endif
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

        int meshRefineLevel = 0;
        createmesh(meshRefineLevel);
      
        setup_system();

        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

        timeloop(0,1);

      /*
        int cycle = 0;
        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          {
            TimerOutput::Scope t(computing_timer, "output");
            //output_results(cycle);
          }
       */
        computing_timer.print_summary();
        computing_timer.reset();

        pcout << std::endl;
  }
} // namespace incompressible

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace incompressible;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      NavierStokes<2> flow;
      flow.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
