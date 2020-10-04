#include <iostream>
#include <cmath>
#include <cf_2.h>
#include <vector>
#include <grid2d.h>
#include <eno_advection.h>
#include <setup.h>

using namespace std;

//Setup Velocity Stuff
class velocity_X :CF_2
{

public: double operator()(double x, double y) const{
        return -0.01;
    }
};


class velocity_Y :CF_2
{
public: double operator()(double x, double y) const{
        return 0.01;
    }
};





//_____________________
int main()
{

    int N = 51;
    int M = 51;
    double xmin = -1;
    double xmax = 1;
    double ymin = -1;
    double ymax = 1;
    double coeff = 0.25;
    double TF = 2*3.1415; //2*3.1415 the usual


    // lets use the fancy constructor for setup and velocities
    Setup setup(N, M, xmin, xmax, ymin, ymax);
    velocity_X vx;
    velocity_Y vy;

    vector<double> VelX;
    vector<double> VelY;

    //Setup mesh in X
    vector<double> MeshX = setup.meshgridX();

    //Setup mesh in Y
    vector<double> MeshY = setup.meshgridY();

    //Generate initial condition
    vector<double> t0 = setup.initial_cond(MeshX, MeshY);

    //Generate velocity in x
     for (int i=0; i < N*M;i++){
         VelY.push_back(vy.operator()(MeshX[i], MeshY[i]));
         VelX.push_back(vx.operator()(MeshX[i], MeshY[i]));
     }

     //Generate true solution
     vector<double> True = setup.true_solution(MeshX, MeshY, VelX, VelY, TF);


     cout<<"Setup Complete"<<endl;

     // lets use the fancy constructor for solver
     ENO_Advection Eno(N, M, xmin, xmax, ymin, ymax, t0, True, VelX, VelY);
     // lets use the fancy constructor for solver
     //ENO_Advection EnoCenter(N, M, xmin, xmax, ymin, ymax, t0, True, VelX, VelY);

     // calculate dt
     double dt = setup.calc_dt(coeff);
     double j = 0;

     //Timesteping
     while (j <= TF){
         Eno.one_Step_ENO(dt);
        // EnoCenter.Central_Advection(dt);
         j=j+dt;
         cout << "Timestep at t= " <<j<< " complete" << endl;
     }

    cout<<"Solution at Final Time Reached"<<endl;

    //Check norm
    cout<<"Norms for ENO scheme:"<<endl;
    Eno.get_norms();

    //cout<<"Norms for CENTRAL scheme:"<<endl;
    //EnoCenter.get_norms();


    // Save Solution
    //Eno.save_vtk("True_ENO_Const","Est_ENO_Const");
    //EnoCenter.save_vtk("True_CENTER","Est_CENTER");

    return 0;
};


