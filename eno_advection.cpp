#include "eno_advection.h"
#include <iostream>
#include <cmath>
#include <grid2d.h>
#include <vector>
#include <omp.h>

using namespace std;

ENO_Advection::ENO_Advection()
{

}

ENO_Advection::ENO_Advection(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_, vector<double> &tn_,vector<double> &true_,vector<double> &Velx_, vector<double> &Vely_)
{

    N=N_; //number of internal points
    M=M_; //number of internal points
    xmin=xmin_;
    xmax=xmax_;
    ymin=ymin_;
    ymax=ymax_;
    Velx=Velx_;
    Vely=Vely_;
    tn=tn_;
    True=true_;
    //N_in=N-2;
    //M_in=M-2;

}
//______________________________________________________________//
// Solve the solution over one timestep using the Upwind scheme //
vector<double> ENO_Advection::one_Step_Upwind(double dt){
    //constructor for grid functions
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double inX;
    double inY;
    vector<double> t_n1(N*M);
    //loop through space
    for (int j = 0; j < N*M; j++){
         //check the sign velocity in the x direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Velx[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inX = grid.dx_forward(tn,j);
        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inX = grid.dx_backward(tn,j);
        }
        // if velocity is zero, then not sure what to

        //check the sign velocity in the y direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Vely[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inY = grid.dy_forward(tn,j);
        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inY = grid.dy_backward(tn,j);
        }
        // if velocity is zero, then not sure what to

        //put it all together    solution_t+(velocity_x*r*discretization_x+velocity_y*r*discretization_y)
        t_n1[j]=tn[j]-(dt*Velx[j]*inX + dt*Vely[j]*inY);

    }
    tn= t_n1;
    return t_n1;
}


//______________________________________________________________//
// Solve the solution over one timestep using the ENO scheme //
void ENO_Advection::one_Step_ENO(double dt){
    //constructor for grid functions
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx =grid.get_dx();
    omp_set_num_threads(4);
    vector<double> t_n1(N*M);
    //loop through space
#pragma omp parallel for
    for (int j = 0; j < N*M; j++){
        cout<<"Thread number " << omp_get_thread_num() <<endl;
        double inX;
        double inY;
        double Center1;
        double Center2;
        double CenterX;
        double CenterY;
         //check the sign velocity in the x direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Velx[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inX = grid.dx_forward(tn,j);

            //Calculate second derivative at j and j-1
            Center1 = grid.dx_centered(tn,j);
            Center2 = grid.dy_centered(tn,j-1);

            //Use minmid to choose the correct second derivative and multiply by coefficent dx/2
            CenterX = grid.minmod(Center1, Center2)*(dx/2);


        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inX = grid.dx_backward(tn,j);

            //Calculate second derivative at j and j+1
            Center1 = grid.dx_centered(tn,j);
            Center2 = grid.dy_centered(tn,j+1);

            //Use minmid to choose the correct second derivative and multiply by coefficent dx/2
            CenterX = grid.minmod(Center1, Center2)*(-dx/2);
        }
        // if velocity is zero, then not sure what to

        //check the sign velocity in the y direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Vely[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inY = grid.dy_forward(tn,j);

            //Calculate second derivative at j and j-M
            Center1 = grid.dx_centered(tn,j);
            Center2 = grid.dy_centered(tn,j-M);

            //Use minmid to choose the correct second derivative and multiply by coefficent dx/2
            CenterY = grid.minmod(Center1, Center2)*(dx/2);
        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inY = grid.dy_backward(tn,j);

            //Calculate second derivative at j and j+M
            Center1 = grid.dx_centered(tn,j);
            Center2 = grid.dy_centered(tn,j+M);

            //Use minmid to choose the correct second derivative and multiply by coefficent dx/2
            CenterY = grid.minmod(Center1, Center2)*(-dx/2);
        }
        // if velocity is zero, then not sure what to

        //put it all together    solution_t+(velocity_x*r*discretization_x+velocity_y*r*discretization_y)
        t_n1[j]=tn[j]-(dt*Velx[j]*(inX+CenterX) + dt*Vely[j]*(inY+CenterY));

    }

    tn= t_n1;
}

// Solve the solution over one timestep using the ENO scheme //
void ENO_Advection::Central_Advection(double dt){
    //constructor for grid functions
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx =grid.get_dx();
    omp_set_num_threads(4);
    vector<double> t_n1(N*M);
    //loop through space
#pragma omp parallel for
    for (int j = 0; j < N*M; j++){
        cout<<"Thread number " << omp_get_thread_num() <<endl;
        double inX;
        double inY;
        double Center1;
        double Center2;
        double CenterX;
        double CenterY;
         //check the sign velocity in the x direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Velx[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inX = grid.dx_forward(tn,j);

            //Calculate second derivative at j and j-1
            CenterX = grid.dx_centered(tn,j);

        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inX = grid.dx_backward(tn,j);

            //Calculate second derivative at j and j+1
            CenterX = grid.dx_centered(tn,j);
        }
        // if velocity is zero, then not sure what to

        //check the sign velocity in the y direction -- WHAT IF THE VELOCITY IS EXACTLY ZERO???
        if (Vely[j]>=0){
            //if velocity is positive, the front is moving to the right and need forward method
            inY = grid.dy_forward(tn,j);

            //Calculate second derivative at j and j-M
            CenterY = grid.dx_centered(tn,j);
        }
        else{
            //if velocity is negative, the front is moving to the right and need backward method
            inY = grid.dy_backward(tn,j);

            //Calculate second derivative at j and j+M
            CenterY = grid.dx_centered(tn,j);
        }
        // if velocity is zero, then not sure what to

        //put it all together    solution_t+(velocity_x*r*discretization_x+velocity_y*r*discretization_y)
        t_n1[j]=tn[j]-(dt*Velx[j]*(inX+CenterX) + dt*Vely[j]*(inY+CenterY));

    }

    tn= t_n1;
}

void ENO_Advection::get_norms(){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    grid.two_Norm(True, tn);
    grid.Inf_Norm(True, tn);
}

void ENO_Advection::save_vtk(std::string file_name,std::string file_name2){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    //save as .vtk file
    grid.initialize_VTK_file(file_name);
    grid.print_VTK_Format(True,"true",file_name);
    grid.initialize_VTK_file(file_name2);
    grid.print_VTK_Format(tn,"tn",file_name2);
}

//______________________________

