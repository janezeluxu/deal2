#ifndef dynamic_H_   /* Include guard */
#define dynamic_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

class dynamic
{
public:
    dynamic(double viscosity,
            DoFHandler<2> &dof_handler,
            FESystem<2> &fe,
            Triangulation<2> &triangulation);
    void getInitialtau(int Grid_size,double vmax,double *cellArea, std::vector<std::vector<DoFHandler<2>::
    active_cell_iterator> > &v_to_e_list, BlockVector<double> &n_solution, double *taumELE, double **taum);
    void getTave(int Grid_size,double*** &gijG,
                 BlockVector<double> &n_solution, double *Tave);
    void relxation_ele(int Grid_size, double **taum_g,double *Tave,double **taum, double *taumEle);
    void getBCinfo(std::vector<int> &IBC, double *bcTag);
    void additionalDynamicinfo(std::vector<std::vector<DoFHandler<2>::
                            active_cell_iterator> > &v_to_e_list,
                               double *cellArea,
                               double **lhsG, double **rhsG,
                               double **xyuHG, double**xyduHxG, double **xyduHyG);
                    //double **lhsG,double**xpowerG, double **ypowerG);
    void getLSintlhs(double xmax,double xmin,double ymax,double ymin,
                    int nElement,std::vector<std::vector<double> > xEle,
                    std::vector<std::vector<double> > yEle, int p,
                    std::vector<std::vector<double> > JxEle,
                    double *lhs1D, double *rhs1D);
    void calcuHxycoeff(Point<2> vpoint,
                int p,double xmin,double xmax,double ymin,double ymax,
                       double *xyuH,double *xyduHx,double *xyduHy);
    void buildDynamicMeshinfo(std::map<int,int> &vertices_to_dofs,std::vector<std::vector<DoFHandler<2>::active_cell_iterator> > &v_to_e_list,std::vector<std::vector<int> > &v_to_e_indices, AffineConstraints<double> &nonzero_constraints,
                              double **vertex_value_u_map,double **vertex_value_v_map,double **vertex_value_p_map,double *cellArea,std::vector<int> &IBC, double *bcTag);
    void getBCnode(int nElement,
                std::vector<std::vector<DoFHandler<2>::
                active_cell_iterator> > &v_to_e_list,
                double **DOF_value,int i,
                std::map<int,int> &vertices_to_dofs,
                   int verticesSiz,double xc, double yc,
                   bool &vertex_flag_u,int  &max_index_u,
                   bool &vertex_flag_v,int  &max_index_v,
                   bool &vertex_flag_p,int  &max_index_p);
    int getMaxIndex(int vertex_flag,
                    std::vector<int> bc_test_vertex,
                    std::vector<double> bc_test_value);
    void getbc_test_value(bool u_flag,int v_index,int i,int ui_vertex_dof,
                        double **DOF_value,Point<2> &v,
                          double xc,double yc,int u_vertex_dof,
                          int v_vertex_dof,
                          std::vector<int> &flag,
                        std::vector<int> &bc_test_vertex_u,
                        std::vector<double> &bc_test_value_u);
    void getc1c2(BlockVector<double> &n_solution,
                 BlockVector<double> &n_solution_time_derivative,
                 std::vector<std::vector<DoFHandler<2>::
               active_cell_iterator> > &v_to_e_list,
               double **vertex_value_u_map,
               double **vertex_value_v_map,
               double **vertex_value_p_map,
               double *cellArea,std::vector<int> IBC,
                 double **lhsG, double **rhsG,
                 double **xyuHG, double**xyduHxG, double **xyduHyG,
               double **taumele_g);
    int determinc2(std::vector<std::vector<double> > uEle,
                    std::vector<std::vector<double> > vEle,
                    int patchOrder,
                    std::vector<double> PatchArea,
                    int nElement);
    void getLSint(double xmax,double xmin,double ymax,double ymin,
                  int nElement,std::vector<std::vector<double> > uEle,
                  std::vector<std::vector<double> > xEle,
                  std::vector<std::vector<double> > yEle,
                  bool flag,Point<2> & bcvertex,double bcvalue,
                  int pH,
                  std::vector<std::vector<double> > JxEle,
                  FullMatrix<double> lhs,double **rhsxy,
                  Vector<double> &LSu);
    void calcuH(Point<2> x, Vector<double> LS,int patchOrder,
                double xmin,double xmax,double ymin,double ymax,
                double *xyuH,double *xyduHx,double *xyduHy,
                double &uH, double *duH);
    void mapxy(int nElement,std::vector<std::vector<double> > xEle,
                        std::vector<std::vector<double> > yEle,
                        double xmax,double xmin,double ymax,double ymin,
                        double **mapx,
                        double **mapy);
    void mapxy(double xEle,double yEle,
               double xmax,double xmin,double ymax,double ymin,
               double &x,double &y);
    void matrixSolve(int mSize,FullMatrix<double> lhs,
                     Vector<double> rhs,Vector<double> &sol);
    void matrixSolve(int mSize,Vector<double> lhs,Vector<double> rhs,double &sol);
    double getnorm(int sizeM,Vector<double> MM);
    std::vector<types::global_dof_index> dofs_per_block;
    
private:
    double           viscosity;
    DoFHandler<2> &dof_handler;
    //FECollection<2> &fe;
    FESystem<2> &fe;
    Triangulation<2> &triangulation;
    //QCollection<2>     quadrature_collection;
    //QGauss<2> quadrature_formula;
    //AffineConstraints<double> &nonzero_constraints;
    int *v_flag;
    //BlockVector<double> &n_solution;
    //BlockVector<double> &n_solution_time_derivative;
};
}
#endif // FOO_H_
