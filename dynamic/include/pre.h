#ifndef pre_H_   /* Include guard */
#define pre_H_

//using namespace incompressible;
namespace incompressible
{
using namespace dealii;

    void createmesh(Triangulation<2> &triangulation,MappingQ<2> mapping,int n_refines);
    
    void getflux(const Point<2> & p1, const Point<2> & p2, char flux_type,bool &ifconsist, double &valuex, double &valuey);
    
    class Dirichlet_BC : public Function<2>
    {
    public:
        Dirichlet_BC(int bc_id,double time) : Function<2>(3)
        {
            boundaryID=bc_id;
            t = time;
        }
        
        int boundaryID;
        double t;
        
        virtual double value( const Point<2> & p,
                             const unsigned int component) const;
        virtual void vector_value( const Point<2> &p,
                                  Vector<double> &  values) const;
    };
    
    /*
    class BoundaryValues : public Function<2>
    {
    public:
        BoundaryValues() : Function<2>(3) {}
        virtual double value(const Point<2> &p,
                             const unsigned int component) const;
        
        virtual void   vector_value(const Point <2>    &p,
                                    Vector<double> &values) const;
    };
     */
}
#endif // FOO_H_
