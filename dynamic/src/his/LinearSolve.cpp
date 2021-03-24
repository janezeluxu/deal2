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
#include "../include/linear.h"
namespace incompressible
{
    using namespace dealii;
    
    /*
    template <class PreconditionerMp>
    BlockSchurPreconditioner<PreconditionerMp>::BlockSchurPreconditioner(double                           gamma,double viscosity,const BlockSparseMatrix<double> &S,const SparseMatrix<double> & P,const PreconditionerMp & Mppreconditioner)
    : gamma(gamma),
    viscosity(viscosity),
    stokes_matrix(S),
    pressure_mass_matrix(P),
    mp_preconditioner(Mppreconditioner)
    {
        A_inverse.initialize(stokes_matrix.block(0, 0));
    }
    
    
    template <class PreconditionerMp>
    void BlockSchurPreconditioner<PreconditionerMp>::vmult(
                                                           BlockVector<double> &      dst,
                                                           const BlockVector<double> &src) const
    {
        Vector<double> utmp(src.block(0));
        {
            SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
            SolverCG<>    cg(solver_control);
            dst.block(1) = 0.0;
            cg.solve(pressure_mass_matrix,
                     dst.block(1),
                     src.block(1),
                     mp_preconditioner);
            dst.block(1) *= -(viscosity + gamma);
        }
        {
            stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
            utmp *= -1.0;
            utmp += src.block(0);
        }
        A_inverse.vmult(dst.block(0), utmp);
    }
     */
    
}


