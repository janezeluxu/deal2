#ifndef NavierStoke_H_   /* Include guard */
#define NavierStoke_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    
class NavierStokes
{
public:
    NavierStokes(const unsigned int degree,const unsigned int mapdegree);
    void run();
    
    std::vector<types::global_dof_index> dofs_per_block;
    FESystem<2>      fe;
    Triangulation<2> triangulation;
     MappingQ<2> mapping;
    DoFHandler<2>    dof_handler;
    
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;
    //BlockSparsityPattern      sparsity_pattern;
    //BlockSparseMatrix<double> system_matrix;
    //SparseMatrix<double>      pressure_mass_matrix;
    
    BlockVector<double> n_solution;
    BlockVector<double> n_solution_time_derivative;
    BlockVector<double> np1_solution;
    BlockVector<double> np1_solution_time_derivative;
    
    //BlockVector<double> newton_update;
    //BlockVector<double> system_rhs;
    BlockVector<double> evaluation_point;
    
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    //SparseMatrix<double>      pressure_mass_matrix;
    
    BlockVector<double> newton_update;
    BlockVector<double> system_rhs;
    
    BlockVector<double> npaf_solution;
    BlockVector<double> npam_solution_time_derivative;
    
//private:
    void preprocess();
    void setup_meshDOF();
    void timeloop();
    void postprocess();
    
    void setDof();
    void initialize(double*** &diffusive_k,double*** &div_SF_u,
                    double**** &grad_SF_u,double**** &SF_u,
                    double*** &SF_p,double**** &grad_SF_p,
                    double*** &gijG);
    void pre_compute(double *** &diffusive_k,
    double*** &div_SF_u,double**** &grad_SF_u,
    double**** &SF_u, double*** &SF_p,double**** &grad_SF_p,
                     double*** &gijG);
    void getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij);
    void startingTimeloop(int num_start);
    void applyBC(double t);
    void predictor();
    void corrector(double *** &diffusive_k,
    double*** &div_SF_u,double**** &grad_SF_u,
    double**** &SF_u, double*** &SF_p,
    double**** &grad_SF_p,double*** &gijG,double *taumEle);
    void updator();
    void output_results () const;
    void iterPC();
    void setup_system(double *** diffusive_k,
    double*** &div_SF_u,double**** &grad_SF_u,
    double**** &SF_u, double*** &SF_p,
    double**** &grad_SF_p,double*** &gijG,double *taumEle);
    void getTaum(double tau_dyn, double* gij,double dt,
                 double uele,double vele,
                 double &tauM, double &tauC);
    void linear_solve(BlockVector<double> &newton_update);
    void update_flow();
    void getResidual(double &current_res);
    
    const unsigned int           degree;
    const unsigned int           mapdegree;
    int                          RefineLevel;
    
    double           viscosity;
    
};
}
#endif // FOO_H_
