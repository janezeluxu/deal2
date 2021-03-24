#ifndef user_input_H_   /* Include guard */
#define user_input_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;
    extern int unsigned max_degree;
    extern int max_cycle;
    extern std::vector<int> BC_ID;
    extern std::vector<char> diff_type;
    extern std::vector<char> pressure_type;
    extern std::vector<double> diff_value_x;
    extern std::vector<double> diff_value_y;
    extern std::vector<double> pressure_value_x;
    extern std::vector<int> component_mask;
    extern std::vector<int> component_geotag;
    extern int nBC_com;

    extern std::vector<int> BC_ID;
    extern double Re ;
    extern unsigned int dim;
    //extern double gamma ;
    //extern double dt;
    //extern double alphaf;
    //extern double alpham ;
    extern int maxNewtonIter ;
    extern double tolerance ;
    //extern int nx ;
    extern int unsigned min_degree;
    extern std::ostringstream velocity_sol;
    extern std::ostringstream pressure_sol;
    extern std::ostringstream outputmeshfilename;

    extern std::ostringstream filenameRefu;
    extern std::ostringstream filenameRefp;
    extern int degreeRef;
    extern int n_Ref;
    extern bool calcError;
    extern std::ostringstream vtkfilename;
    extern std::string meshfile;
    extern bool KellyError;
    
    void read_input(std::string inputfile);
    double Dirichlet_u(const dealii::Point<2> & p,int BC_ID);
    double Dirichlet_v(const dealii::Point<2> & p,int BC_ID);
    double Dirichlet_p(const dealii::Point<2> & p,int BC_ID);
    void diff_flux(const Point<2> & p1, const Point<2> & p2,
                   bool &ifconsist, double &valuex, double &valuey);
    void pressure_flux(const Point<2> & p1, const Point<2> & p2,
                       bool &ifconsist, double &valuex, double &valuey);

}
#endif // FOO_H_
