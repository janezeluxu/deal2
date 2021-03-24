#ifndef linear_H_   /* Include guard */
#define linear_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class LinearSolve
{
public:
    LinearSolve();
    void solve();
    void preconditioner();
};
 
    template <class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
    public:
        BlockSchurPreconditioner(double                           gamma,
                                 double                           viscosity,
                                 const BlockSparseMatrix<double> &S,
                                 const SparseMatrix<double> &     P,
                                 const PreconditionerMp &         Mppreconditioner);
        void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;
    private:
        const double                     gamma;
        const double                     viscosity;
        const BlockSparseMatrix<double> &stokes_matrix;
        const SparseMatrix<double> &     pressure_mass_matrix;
        const PreconditionerMp &         mp_preconditioner;
        SparseDirectUMFPACK              A_inverse;
    };

    
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
    
}
#endif // FOO_H_
