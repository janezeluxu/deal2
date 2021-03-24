#ifndef refine_mesh_H_   /* Include guard */
#define refine_mesh_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    
class refine_mesh
{
public:
    refine_mesh(Triangulation<2> &triangulation,
                hp::DoFHandler<2> &dof_handler,
                BlockVector<double> &present_solution,
                hp::FECollection<2> &fe);
    
    void perfrom_refinement(hp::QCollection<1> &face_quadrature_collection,int cycle,
        SolutionTransfer<2,
        BlockVector<double>,hp::DoFHandler<2>  > &solution_transfer,
        Vector<float> &estimated_error_per_cell,
        Vector<float> &true_error_per_cell);
    void error_indicator(hp::QCollection<1> &face_quadrature_collection,
                         Vector<float> &estimated_error_per_cell);
    void error_vms(Vector<float> &estimated_error_per_cell,double &Error_total);
    void error_Kelly(hp::QCollection<1> &face_quadrature_collection,
                Vector<float> &estimated_error_per_cell);
    
    void getGij(const DerivativeForm<1, 2, 2> &JacobiInverse, double* &gij);
    void getTaum(int degree, double* gij,double uele,double vele,double &tauM, double &tauC);
    void estimate_smoothness(Vector<float> &smoothness_indicators);

    void savemesh(int cycle,Vector<float> estimated_error_per_cell,
                  Vector<float> true_error_per_cell,
                  Vector<float> smoothness_indicators);
    
    template <typename T>
    void resize(Table<2, T> &coeff, const unsigned int N);
    
    std::pair<bool, unsigned int>
    predicate(const TableIndices<2> &ind);
    
private:
    Triangulation<2> &triangulation;
    hp::DoFHandler<2> &dof_handler;
    BlockVector<double> &present_solution;
    hp::FECollection<2> &fe;
    hp::QCollection<2>     quadrature_collection;
    //Vector<float> &estimated_error_per_cell;
    //SolutionTransfer<2, BlockVector<double>,hp::DoFHandler<2>  > &solution_transfer;
    
    hp::QCollection<2>                    fourier_q_collection;
    std::shared_ptr<FESeries::Fourier<2>> fourier;
    std::vector<double>                   ln_k;
    Table<2, std::complex<double>>        fourier_coefficients;
    //Vector<float> &smoothness_indicators;
};
}
#endif // FOO_H_
