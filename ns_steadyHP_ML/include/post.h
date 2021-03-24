#ifndef post_H_   /* Include guard */
#define post_H_

//#include "../include/meshdof.h"
//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    
class Post
{
public:
    static void error_estimator();
    void restart_save();
    static void aerofynamic_force();
    void savetovtk(hp::DoFHandler<2> &dof_handler,BlockVector<double> &present_solution);
    void savesolution(BlockVector<double> present_solution, int istp);
    void savemesh(int cycle,int vertices_size,int edgeSize,hp::DoFHandler<2> &dof_handler,hp::FECollection<2> &fe, AffineConstraints<double> &nonzero_constraints);
    void ErrorEstimate(FESystem<2> fe, DoFHandler<2> &dof_handler,
                       int degree,
                       BlockVector<double> n_solution,
                       double &L2_velocity,double &H1_velocity,
                       double &L2_presure,double &H1_pressure);
    void getExact(double x,double y,
             double &u,double &v,double &p,
             double &ux,double &vx,double &px,double &uy,double &vy,
                  double &py);
    
    void createRefmesh(Triangulation<2> &triangulation,DoFHandler<2> &dof_handler,int degreeRef,int n_Ref,BlockVector<double> &n_solution);
    void readinRefSolution(BlockVector<double> &n_solution);
    void getTrueError(hp::FECollection<2> &fe,
                      Triangulation<2> &triangulation,
                      hp::DoFHandler<2> &dof_handler,
                      BlockVector<double> &n_solution,
                      Vector<float> &true_error_per_cell);
};
}
#endif // FOO_H_
