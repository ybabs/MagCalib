#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "Eigen/Dense"
#include <string>
using Eigen::EigenBase;   

const double PI = std::atan(1.0)*4;

template <typename Derived>
std::string get_shape(const EigenBase<Derived>& x)
{
    std::ostringstream oss;
    oss  << "(" << x.rows() << ", " << x.cols() << ")";
    return oss.str();
}

class Mag
{
public:
    explicit Mag() = default;
    ~Mag();
    void processFile();
    void fitEllipsoid(Eigen::VectorXd& magX, Eigen::VectorXd&  magY, Eigen::VectorXd&  magZ );    void computeError();

   private:
   std::ifstream inputFile;
   Eigen::MatrixXd Q;
   Eigen::VectorXd n;
   Eigen::VectorXd mag_x;
   Eigen::VectorXd mag_y;
   Eigen::VectorXd mag_z;
   double d;


};
