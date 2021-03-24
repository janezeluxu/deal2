#ifndef post_H_   /* Include guard */
#define post_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    
class Post
{
public:
    Post(int degree);
    static void error_estimator();
    void restart_save();
    static void aerofynamic_force();
    void savetovtk(DoFHandler<2> &dof_handler,BlockVector<double> present_solution, unsigned int istp);
    void savesolution(BlockVector<double> present_solution,BlockVector<double> n_solution_time_derivative,unsigned int istp);
    void savemesh(int cycle,DoFHandler<2> &dof_handler,FESystem<2> &fe,MappingQ<2> mapping, AffineConstraints<double> &nonzero_constraints);
    void ErrorEstimate(MappingQ<2> mapping,
                       FESystem<2> fe, DoFHandler<2> &dof_handler,
                  int EleNum,double*** div_SF_u,double**** grad_SF_u,
                       double**** SF_u, double*** SF_p,double**** grad_SF_p, BlockVector<double> present_solution,
                       double &L2_velocity,double &H1_velocity,
                       double &L2_presure,double &H1_pressure);
    
    void getExact(double x,double y,double &u,double &v,double &p,
                  double &ux,double &vx,double &px,double &uy,double &vy,
                  double &py);
private:
    int              degree;
};
}
#endif // FOO_H_
