#ifndef ENO_ADVECTION_H
#define ENO_ADVECTION_H
#include <vector>
#include <grid2d.h>


class ENO_Advection
{
public: 
    ENO_Advection();
    ENO_Advection(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_,std::vector<double> &tn_,std::vector<double> &true_,std::vector<double> &Velx_, std::vector<double> &Vely_);
    std::vector<double> one_Step_Upwind(double dt);
    void one_Step_ENO(double dt);
    void save_vtk(std::string file_name,std::string file_name2);
    void Central_Advection(double dt);
    void get_norms();
    //std::vector<double> initial_cond(std::vector<double> X, std::vector<double> Y);

private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    std::vector<double> tn,True,Velx, Vely;
};

#endif // ENO_ADVECTION_H
