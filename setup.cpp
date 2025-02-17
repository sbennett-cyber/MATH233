#include "setup.h"
#include <iostream>
#include <cmath>
#include <grid2d.h>
#include <vector>

using namespace std;
Setup::Setup()
{

}

Setup::Setup(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{

    N=N_; //number of internal points
    M=M_; //number of internal points
    xmin=xmin_;
    xmax=xmax_;
    ymin=ymin_;
    ymax=ymax_;
    //N_in=N-2;
    //M_in=M-2;

}

// Meshgrid in X //
vector<double> Setup::meshgridX(){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    vector<double> X;
    double valX = xmin+dx;
    //cout<<valX<<endl;

    for (int j=0;j<M;j++){
        for (int i=0;i<N;i++){
            X.push_back(valX);
            //cout<<valX<<endl;
            valX=valX+dx;
        }
        valX = xmin+dx;
    }
    return X;
}
//______________//
// Meshgrid in Y
vector<double> Setup::meshgridY(){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    vector<double> Y;
    double valY = ymin+dx;

    for (int j=0;j<M;j++){
        for (int i=0;i<N;i++){
            Y.push_back(valY);
            //cout<<valY<<endl;
        }
    valY=valY+dx;
    }
    return Y;
}
//__________________//
// Initial Condition
vector<double> Setup::initial_cond(vector<double> X, vector<double> Y){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        if ((sqrt((X[j]-0.5)*(X[j]-0.5)+Y[j]*Y[j]))-0.2<=0){
            t0.push_back(1);
        }
        else{
            t0.push_back(0);
        }

    }
     return t0;
}

//______________//
// True Solution
vector<double> Setup::true_solution(vector<double> X, vector<double> Y, vector<double> VelX, vector<double> VelY, double dt){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        if ((sqrt((X[j]-0.5)*(X[j]-0.5)+Y[j]*Y[j]))-0.2<=0){
            t0.push_back(1);
        }
        else{
            t0.push_back(0);
        }

    }
     return t0;
}

//________________________//
// Calculate and return dt
double Setup::calc_dt(double coef){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    cout << "The value for dx = " << dx << " and dt =  " << coef*dx <<endl;
    return coef*dx;
}

