// ./testfunction 0.1 0.001 1 1 1 1 1                  [this argurements are not needed just there because was adapted from another main.cc]
// -106.794678 6.27930438 7.05394284


// first make dlib
// mkdir build
// cd build
// cmake ..
// cmake --build . --config Release
// 


// mpicxx -std=c++11 -O3 -I./dlib-19.17/ testfunction.cpp ./dlib-19.17/dlib/all/source.cpp -lpthread -lX11  -o testfunction

#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include <iostream>

#include <complex.h>
#include <math.h>
#include <cstdlib>

#define NOUTPUT 4   // the outputs: V0, tau4, tau5, rho1

using namespace std;
using namespace dlib;

typedef matrix<double,0,1> column_vector;

// ----------------------------------------------------------------------------------------
double testfunction(double x0, double x1);   


int main(int argc, char **argv){

  //  double AbsW0, AbsA1, a1, theta0, phi1, gstring;
  double theta[6];
  int index; 
  if(argc < 8){ 
    printf(" This program needs 6 input parameters and an index_int:\n" );
    exit(1); 
  } else{ 
    theta[0] = atof(argv[1]);      //AbsW0  
    theta[1] = atof(argv[2]);	   //AbsA1  
    theta[2] = atof(argv[3]);	   //a1     
    theta[3] = atof(argv[4]);	   //theta0 
    theta[4] = atof(argv[5]);	   //phi1   
    theta[5] = atof(argv[6]);	   //gstring
    index = atol(argv[7]);
  }

  auto complex_holder_table = [theta](double x0, double x1){
    return testfunction(x0, x1);              // For Steve's dummy function
  };

  double pi = 4* atan(1.0);
  auto result = find_min_global(complex_holder_table,   
				{5.0,5.0}, 
				{10.0,10.0}, 
				max_function_calls(300));

      column_vector take_results = result.x; 
      cout.precision(9);
      cout << result.y << " " << result.x(0,0) << " " << result.x(0,1) << "\n" << endl << endl;

    return 0; 
}


double testfunction(double x0, double x1){   
  return -(36.0*sin(2*x1)*cos(2*x0) + 12*(x0 + x1) - x1*x1 - x0*x0); 
}

