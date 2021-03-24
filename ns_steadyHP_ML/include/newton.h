#ifndef Newton_H_   /* Include guard */
#define Newton_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class NewtonSolve
{
public:
    NewtonSolve(double viscosity,
                hp::DoFHandler<2> &dof_handler,
                hp::FECollection<2> &fe,
                BlockSparsityPattern &sparsity_pattern,
                AffineConstraints<double> &zero_constraints,
                BlockVector<double> &n_solution,
                BlockSparseMatrix<double> &system_matrix,
                BlockVector<double> &system_rhs);
    
    
    void getGij(const DerivativeForm<1, 2, 2> &Jacobi,double* &gij);
    void setup_system();
    void getTaum(int degree, double* gij,double uele,double vele,double rum, double rvm,double &tauM, double &tauC,double &tauBar);
    
    //void getflux(int boundaryID, char flux_type,bool &ifconsist, double &valuex, double &valuey);
    
    void assemble_rhs();
    
    void linear_solve( BlockVector<double> &newton_update);
    
    void update_flow(BlockVector<double> &newton_update,BlockVector<double> &n_solution);
    
    void check_convergence( double &res);
    void add_flux();
    void getNewtonUpdate(BlockVector<double> &n_solution,BlockVector<double> &newton_update);
    
    void getResidual(BlockVector<double> &n_solution,double &current_res);
    
private:
    double           viscosity;
    hp::DoFHandler<2> &dof_handler;
    hp::FECollection<2> &fe;
    hp::QCollection<2>     quadrature_collection;
    hp::QCollection<1>     face_quadrature_collection;

    BlockSparsityPattern &sparsity_pattern;
    AffineConstraints<double> &zero_constraints;
    
    BlockVector<double> &n_solution;
    BlockSparseMatrix<double> &system_matrix;
    BlockVector<double> &system_rhs;
    
    //BlockSparseMatrix<double> &system_matrix;
    //BlockVector<double> &newton_update;
    //BlockVector<double> &system_rhs;
};
    
}
#endif // FOO_H_
