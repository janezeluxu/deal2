#ifndef Newton_H_   /* Include guard */
#define Newton_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class NewtonSolve
{
public:
    NewtonSolve(int degree, double viscosity, double gamma, double dt,double alphaf,
                double alpham,
                DoFHandler<2> &dof_handler,
                FESystem<2> &fe,
                MappingQ<2> &mapping,
                BlockSparsityPattern &sparsity_pattern,
                ConstraintMatrix &zero_constraints,
                //SparseMatrix<double> &pressure_mass_matrix, //temp
                BlockSparseMatrix<double> &system_matrix,
                BlockVector<double> &system_rhs);
    
    void iterPC(BlockVector<double> &n_solution,BlockVector<double> &n_solution_time_derivative, BlockVector<double> &np1_solution,BlockVector<double> &np1_solution_time_derivative,BlockVector<double> &npaf_solution,BlockVector<double> &npam_solution_time_derivative);
    
    void getGij(const DerivativeForm<1, 2, 2> &Jacobi,double* &gij);
    void setup_system(BlockVector<double> &npaf_solution,
                      BlockVector<double> &npam_solution_time_derivative,
                      double *** diffusive_k,double*** &div_SF_u,double**** &grad_SF_u,
                      double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,double*** &gijG);
    void getTaum(double* gij,double t, double uele,double vele,double &tauM, double &tauC);
    
    //void getflux(int boundaryID, char flux_type,bool &ifconsist, double &valuex, double &valuey);
    
    void assemble_rhs(BlockVector<double> &npaf_solution,
                      BlockVector<double> &npam_solution_time_derivative,
                      double *** diffusive_k,double*** &div_SF_u,double**** &grad_SF_u,
                      double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,double*** &gijG);
    void add_flux(BlockVector<double> &npaf_solution);
    void linear_solve( BlockVector<double> &newton_update);
    
    void update_flow(BlockVector<double> &newton_update,BlockVector<double> &np1_solution,BlockVector<double> &np1_solution_time_derivative);
    
    void check_convergence( double &res);
    
    void getNewtonUpdate(BlockVector<double> &npaf_solution,
                         BlockVector<double> &npam_solution_time_derivative,
                         BlockVector<double> &newton_update,
                         double *** diffusive_k,
                         double*** &div_SF_u,double**** &grad_SF_u,
                         double**** &SF_u, double*** &SF_p,
                         double**** &grad_SF_p,double*** &gijG);
    
    void getResidual(double &current_res);
    
private:
    int              degree;
    double           viscosity;
    double           gamma;
    double           dt;
    double           alphaf;
    double           alpham;
    DoFHandler<2> &dof_handler;
    FESystem<2> &fe;
    MappingQ<2> &mapping;
    BlockSparsityPattern &sparsity_pattern;
    ConstraintMatrix &zero_constraints;
    
    //SparseMatrix<double> &pressure_mass_matrix; //temp
    BlockSparseMatrix<double> &system_matrix;
    BlockVector<double> &system_rhs;
    
    //BlockSparseMatrix<double> &system_matrix;
    //BlockVector<double> &newton_update;
    //BlockVector<double> &system_rhs;
};
    
}
#endif // FOO_H_
