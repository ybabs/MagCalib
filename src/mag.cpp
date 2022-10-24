#include <mag.hpp>
#include <algorithm>
#include <cmath>
// #include <matplot/matplot.h>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <unsupported/Eigen/MatrixFunctions>


using namespace Spectra;

// using namespace matplot;

Mag::~Mag()
{

}

void Mag::processFile()
{
   inputFile.open("../magnetometer.csv");

   if(inputFile)
   {
            std::vector<double> magXVec;
        std::vector<double> magYVec;
        std::vector<double> magZVec;

        std::string line = "";
        while(getline(inputFile, line))
        {
            float accX;
            float accY;
            float accZ;
            float gyroX;
            float gyroY;
            float gyroZ;
            float magX;
            float magY;
            float magZ;

            std::string tempString = "";

            std::stringstream inputString(line);

            getline(inputString, tempString, ',');
            accX = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            accY = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            accZ = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            gyroX = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            gyroY = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            gyroZ= atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            magX = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            magY = atof(tempString.c_str());
            tempString.clear();
            getline(inputString, tempString, ',');
            magZ = atof(tempString.c_str());

            magXVec.push_back(magX);
            magYVec.push_back(magY);
            magZVec.push_back(magZ);

            line = "";
        }

        std::transform(magXVec.begin(), magXVec.end(), magXVec.begin(), 
                        [](float val) -> float {return val * 0.08;});

            std::transform(magYVec.begin(), magYVec.end(), magYVec.begin(), 
                        [](float val) -> float {return val * 0.08;});

            std::transform(magZVec.begin(), magZVec.end(), magZVec.begin(), 
                        [](float val) -> float {return val * 0.08;});


            mag_x = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(magXVec.data(), magXVec.size());
            mag_y = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(magYVec.data(), magYVec.size());
            mag_z = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(magZVec.data(), magZVec.size());
                    

        std::cout << "MAG X Size: " << magXVec.size() << std::endl;
        //    plt::plot(mag_x);
        //     plt::show();

        //    Eigen::VectorXd y = mag_x;
        //    plt::plot(y);
        //    plt::show();
            fitEllipsoid(mag_x,mag_y, mag_z);
   }

   else
   {
      std::cout<< "Unable to Open file" << std::endl;
   }

   
   
}


// m rows x n cols (10 x 10)
void Mag::fitEllipsoid(Eigen::VectorXd& magX, Eigen::VectorXd&  magY, Eigen::VectorXd&  magZ )
{
    Eigen::VectorXd a1 = magX.array().square();
    Eigen::VectorXd a2 = magY.array().square();
    Eigen::VectorXd a3 = magZ.array().square();

    //std::cout << "Shape a1" <<  get_shape(a1) << std::endl;

    //std::cout << "Shape a2" <<  get_shape(a2) << std::endl;

    Eigen::VectorXd a4 = 2 * (magY.array() * magZ.array());

    //std::cout << "Shape a4" <<  get_shape(a4) << std::endl;
    Eigen::VectorXd a5 = 2 * (magX.array() * magZ.array());
    Eigen::VectorXd a6 = 2 * (magX .array()* magY.array());

    Eigen::VectorXd a7 = 2 * magX;
    Eigen::VectorXd a8 = 2 * magY;
    Eigen::VectorXd a9 = 2 * magZ;
    Eigen::VectorXd a10 = Eigen::VectorXd::Ones(magX.size()).transpose();

    // std::cout << "Shape a4" <<  get_shape(a10) << std::endl;

    //Eigen::MatrixXd D(magX.size(), 10);
    Eigen::MatrixXd D(10, magX.size());

    D << a1.transpose(), a2.transpose(), a3.transpose(), a4.transpose(), a5.transpose(), a6.transpose(), a7.transpose(), a8.transpose(), a9.transpose(), a10.transpose(); 

    // std::cout << "Shape D " <<  get_shape(D) << std::endl;

    //  std::cout << "Shape D Transpose " <<  get_shape(D.transpose()) << std::endl;

    Eigen::MatrixXd C {
                   {-1, 1, 1, 0, 0, 0},
                   {1, -1, 1, 0, 0, 0},
                   {1, 1, -1, 0, 0, 0},
                   {0, 0, 0, -4, 0, 0},
                   {0, 0, 0, 0, -4, 0},
                   {0, 0, 0, 0, 0, -4}
    };


    Eigen::MatrixXd S = D * D.transpose();
    Eigen::MatrixXd S11 =  S.topLeftCorner(6,6);
    Eigen::MatrixXd S12 =  S.topRightCorner(6,4);
    Eigen::MatrixXd S21 =  S.bottomLeftCorner(4,6);
    Eigen::MatrixXd S22 = S.bottomRightCorner(4,4); 

    Eigen::MatrixXd temp = C.inverse() * (S11 - (S12 *(S22.inverse() * S21)));

    Spectra::DenseGenMatProd<double> op(temp);

    Spectra::GenEigsSolver<DenseGenMatProd<double>> eigs(op, 1, 6);

    eigs.init();
    int nconv = eigs.compute(SortRule::LargestReal);

    Eigen::VectorXcd evalues; 
    //auto evectors; 
    if(eigs.info() == CompInfo::Successful)
    {
        evalues = eigs.eigenvalues();
    }

    Eigen::VectorXd u1 = eigs.eigenvectors().real();
    Eigen::VectorXd u2 = ((-1* (S22.inverse() * S21)) * u1 );

    //std::cout << "U2 found:\n" << u2  << std::endl;

    Eigen::VectorXd u(u1.size() + u2.size());
    u << u1, u2;
    u = u.transpose();

    Q.resize(3, 3);
    n.resize(3);

    Q << u(0), u(5), u(4),
         u(5), u(1), u(3),
         u(4), u(3), u(2);

    n << u(6), u(7), u(8);

    d = u(9);

    computeError();

}

void Mag::computeError()
{

    auto Qinv = Q.inverse();

    Eigen::VectorXd b = -1 * (Qinv * n);

    // std::cout << "Qinv: \n" << Qinv << std::endl;

     std::cout << "b : = \n" << get_shape(b) << std::endl;

    Eigen::MatrixXd Ainv = (1 / (std::sqrt(n.transpose() * (Qinv * n) - d)) * Q.sqrt()).real();

    //std::cout<< "tmp: " << Ainv << std::endl;

    Eigen::VectorXd calibX(mag_x.size());
    Eigen::VectorXd calibY(mag_y.size());
    Eigen::VectorXd calibZ(mag_z.size());

    float  totalError = 0;
    for (int i = 0; i < mag_y.size(); i++)
    {
        Eigen::VectorXd h(3);
        h << mag_x(i), mag_y(i) , mag_z(i);
        //std::cout << mag_x(i) << " " <<  mag_y(i) << " " << mag_z(i) << std::endl;
        //std::cout << "h: " << get_shape(h) << std::endl;
        Eigen::VectorXd h_b = h - b;
        // std::cout <<  "hb Shape:" << get_shape(h_b) << std::endl;
        auto hHat = Ainv * h_b;
        //std::cout <<  "hhat Shape:" << get_shape(hHat) << std::endl;
        calibX(i) = hHat(0);
        calibY(i) = hHat(1);
        calibZ(i) = hHat(2);
        auto mag = hHat.transpose().dot(hHat);
        //std::cout << "Mag" << mag << std::endl;
        float error = (mag -1) * (mag - 1);
        //std::cout << "Error " << error << std::endl;
        totalError+=error;
        //std::cout << "Total Error " << totalError << std::endl;


    }

 
   // copy vectors back
   std::vector<float> x(&calibX[0], calibX.data()+calibX.cols()*calibX.rows());
   std::vector<float> y(&calibY[0], calibY.data()+calibY.cols()*calibY.rows());
   std::vector<float> z(&calibZ[0], calibZ.data()+calibZ.cols()*calibZ.rows());

    // auto l = plot3(x, y, z);
    // show();

    //plt::scatter(calibX, calibY, calibZ, {{"color", "red"}, {"label", "a circle!"}});

    std::cout<< "Error: " <<totalError << std::endl;
    

}
