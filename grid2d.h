#ifndef GRID2D_H
#define GRID2D_H

#include <vector>



class Grid2D
{
private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    double dx,dy;
public:
    Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
 //   Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double get_dx();
    int i_from_n(int n);
    int j_from_n(int n);
    int n_from_ij(int i, int j);
    double dx_forward(std::vector<double> function, int n);
    double dx_backward(std::vector<double> function, int n);
    double dy_forward(std::vector<double> function, int n);
    double dy_backward(std::vector<double> function, int n);
    void print();
    void initialize_VTK_file(std::string file_name);
    double x_from_n(int n);
    double y_from_n(int n);
    double minmod(double D1, double D2);
    double dx_centered(std::vector<double> function, int n);
    double dy_centered(std::vector<double> function, int n);
    void two_Norm(std::vector<double> True, std::vector<double> Est);
    void Inf_Norm(std::vector<double> True, std::vector<double> Est);
    void print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name );

};

#endif // GRID2D_H
