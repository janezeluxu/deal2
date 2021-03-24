#ifndef user_input_H_   /* Include guard */
#define user_input_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    extern int degree;
    extern int mapdegree;
    extern std::vector<int> BC_ID;
    extern std::vector<char> diff_type;
    extern std::vector<char> pressure_type;
    extern std::vector<double> diff_value_x;
    extern std::vector<double> diff_value_y;
    extern std::vector<double> pressure_value_x;
    extern std::vector<double> pressure_value_y;
    extern double Re ;
    extern double gamma ;
    extern double dt;
    extern double alphaf;
    extern double alpham ;
    extern int maxNewtonIter ;
    extern double tolerance ;
    
    // dynamic setups
    extern int surfaceTag;
    extern bool dynamicFlag;
    extern int dynamicBC;
    
    extern int stage1_outter;
    extern int stage1_inner;
    extern int stage2;
    extern double c1_start;

    void defineMesh(std::vector<unsigned int> &subdivisions,double &x1,double &x2,double &y1,double &y2);
    double Dirichlet_u(const dealii::Point<2> & p,int BC_ID,double t);
    double Dirichlet_v(const dealii::Point<2> & p,int BC_ID,double t);
    double Dirichlet_p(const dealii::Point<2> & p,int BC_ID,double t);
    void diff_flux(const Point<2> & p1, const Point<2> & p2,
                   bool &ifconsist, double &valuex, double &valuey);
    void pressure_flux(const Point<2> & p1, const Point<2> & p2,
                       bool &ifconsist, double &valuex, double &valuey);

}
#endif // FOO_H_
