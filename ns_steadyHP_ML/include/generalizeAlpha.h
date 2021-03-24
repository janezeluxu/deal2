#ifndef generalizeAlpha_H_   /* Include guard */
#define generalizeAlpha_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class generalizeAlpha
{
public:
    generalizeAlpha( double viscosity);
    void initialize(hp::DoFHandler<2> &dof_handler,
                AffineConstraints<double> &nonzero_constraints,
                std::vector<types::global_dof_index> &dofs_per_block);
    void corrector(int maxNewtonIter, double tolerance, hp::DoFHandler<2> &dof_handler,hp::FECollection<2> fe,AffineConstraints<double> &zero_constraints, AffineConstraints<double> &nonzero_constraints,BlockVector<double> &n_solution);
    
    std::vector<types::global_dof_index> dofs_per_block;
    
private:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    //SparseMatrix<double>      pressure_mass_matrix;
    
    BlockVector<double> newton_update;
    BlockVector<double> system_rhs;
    
    double           viscosity;
};
}
#endif // FOO_H_
