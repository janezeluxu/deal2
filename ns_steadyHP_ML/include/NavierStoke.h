#ifndef NavierStoke_H_   /* Include guard */
#define NavierStoke_H_

//#include "../include/meshdof.h"
//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    
class NavierStokes
{
public:
    NavierStokes();
    void run();
    std::vector<types::global_dof_index> dofs_per_block;
    //FESystem<2>      fe;
    //Triangulation<2> triangulation;
    //DoFHandler<2>    dof_handler;
    //FESystem<2>         fe;
    hp::FECollection<2>    fe_collection;
    Triangulation<2> triangulation;
    hp::QCollection<1>     face_quadrature_collection;
    hp::DoFHandler<2>      dof_handler;
    //hp::QCollection<2>     quadrature_collection;
    //hp::QCollection<1> face_quadrature_collection;
    
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;
    //BlockSparsityPattern      sparsity_pattern;
    //BlockSparseMatrix<double> system_matrix;
    //SparseMatrix<double>      pressure_mass_matrix;
    
    BlockVector<double> n_solution;
    
    //BlockVector<double> newton_update;
    //BlockVector<double> system_rhs;
    BlockVector<double> evaluation_point;
    
//private:
    void setup_meshDOF();
    void loop();
    void mesh_adapt(int cycle,Vector<float> &true_error_per_cell);
    void postprocess();
    
    void applyBC();
    void setDof();
    void initialize();
    void BCDof();
    void output_results () const;

    int                          RefineLevel;
    
};
}
#endif // FOO_H_
