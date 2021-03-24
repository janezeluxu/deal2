#ifndef generalizeAlpha_H_   /* Include guard */
#define generalizeAlpha_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class generalizeAlpha
{
public:
    generalizeAlpha(int degree, double viscosity, double gamma, double dt,double alphaf, double alpham);
    void initialize(DoFHandler<2> &dof_handler,ConstraintMatrix &nonzero_constraints,
                    std::vector<types::global_dof_index> &dofs_per_block);
    void  pre_compute(MappingQ<2> mapping, FESystem<2> fe, DoFHandler<2> &dof_handler,int EleNum, double ***
                      &diffusive_k,double*** &div_SF_u,double**** &grad_SF_u,
                      double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,
                      double*** &gijG);
    void predictor(BlockVector<double> &n_solution,
                   BlockVector<double> &n_solution_time_derivative,
                   BlockVector<double> &np1_solution,
                   BlockVector<double> &np1_solution_time_derivative);
    
    void getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij);
    
    void corrector(MappingQ<2> mapping, int maxNewtonIter, double tolerance, DoFHandler<2> &dof_handler,FESystem<2> fe,ConstraintMatrix &zero_constraints, ConstraintMatrix &nonzero_constraints,BlockVector<double> &n_solution,
                   BlockVector<double> &n_solution_time_derivative,
                   BlockVector<double> &np1_solution,
                   BlockVector<double> &np1_solution_time_derivative,
                   double *** &diffusive_k,
                   double*** &div_SF_u,double**** &grad_SF_u,
                   double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,
                   double*** &gijG);
    void updator(BlockVector<double> &n_solution,
                 BlockVector<double> &n_solution_time_derivative,
                 BlockVector<double> &np1_solution,
                 BlockVector<double> &np1_solution_time_derivative);
    
    std::vector<types::global_dof_index> dofs_per_block;
    
private:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    //SparseMatrix<double>      pressure_mass_matrix;
    
    BlockVector<double> newton_update;
    BlockVector<double> system_rhs;
    
    BlockVector<double> npaf_solution;
    BlockVector<double> npam_solution_time_derivative;
    
    
    int              degree;
    double           viscosity;
    double           gamma;
    double           dt;
    double           alphaf;
    double           alpham;
};
}
#endif // FOO_H_
