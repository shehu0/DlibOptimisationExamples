// first make dlib
// mkdir build
// cd build
// cmake ..
// cmake --build . --config Release
// 


// mpicxx -std=c++11 -O3 -I./dlib-19.17/ lv3terms.cpp ./dlib-19.17/dlib/all/source.cpp -lpthread -lX11 -o a_4par

// with foru parameters
// ./lv3terms 0.305175782249969465E-03 0.228686632513999966E+01 0.419310240354582575E+01 0.122912548184394829E+00 1
// -2.24617385e-121 23.5041921 8.47336553e+25


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
double p11169E_lv(double tau4, double tau5, const double *theta);

int main(int argc, char **argv){

  //  double AbsW0, AbsA1, a1, theta0, phi1, gstring;
  double theta[4];
  int index; 
  if(argc < 6){ 
    printf(" This program needs 6 input parameters and an index_int:\n" );
    exit(1); 
  } else{ 
    theta[0] = atof(argv[1]);      //AbsW0  
    theta[1] = atof(argv[2]);	   //AbsA1  
    theta[2] = atof(argv[3]);	   //a1     
    theta[3] = atof(argv[4]);	   //gstring
    index = atol(argv[5]);
  }

  auto complex_holder_table = [theta](double x0, double x1){ 
    return p11169E_lv(x0, x1, theta);
  };

  double pi = 3.141592654;//  double pi = 4* atan(1.0);
  auto result = find_min_global(complex_holder_table, 
				{0.0,0.0}, // lower bounds
				{40.0,100000000000000000000000000.0}, 
				max_function_calls(300));
  
  cout.precision(9);
  cout << result.y << " " << result.x(0,0) << " " << result.x(0,1) << "\n" << endl;

    return 0; 
}


double p11169E_lv(double tau1, double tau2, const double *theta){
  double pi = 3.141592654;//  double pi = 4* atan(1.0);
  double Kcs = 0.1;  // could be another parameter.. fit to COBE?
  double k111 = 9.0, k222 = 9.0;

  // the parameters
  double AbsW0, AbsA1, a1, theta0, phi1, gstring;

  // variable to find at the minimum of potential
  //  double tau4, tau5, rho1; 

  //  double logZero = -DBL_MAX;
  //  double loglike = 0.0;

  // transfer the parameters
  AbsW0  = theta[0];
  AbsA1  = theta[1];
  a1     = theta[2];
  //  theta0 = theta[3];
  theta0 = 0.0;
  //  phi1   = theta[4];
  phi1   = 0.0;
  gstring= theta[3];

  double chiiP11169 = -540.0; // ??? what is the value for p11169 ?
  double xi = -1.2*chiiP11169/(16.0 * pi*pi*pi); // fixed 
  double xiHat = xi/pow(gstring,3.0/2.0);

  double Vol = sqrt(2.0)*(pow(tau2,3.0/2.0) - pow(tau1,3.0/2.0))/9.0;
  double expKhalerPot = exp( Kcs - 2.0*log(Vol+xiHat/2.0) );  // use log to base e 
  double xi2 = xiHat*xiHat;
  double V2 = Vol*Vol;

  double valpp = expKhalerPot * AbsW0 * AbsW0 *
    3.0*xiHat*(xi2+7.0*xiHat*Vol+V2)/
    ( (Vol - xiHat)*(xiHat + 2*Vol)*(xiHat + 2*Vol) );
  

  double vnp1 = 2.0 * expKhalerPot * AbsW0 * AbsA1 * exp(-a1*tau1) * (-1.0) * // cos(a1*rho1 + theta0 - phi1) *
    ( 2.0*a1*tau1*(4.0*xi2 + xiHat*Vol + 4.0* V2)/(2.0*(Vol - xiHat)*(xiHat + 2*Vol)) + 
      3.0*xiHat*(xi2 + 7.0*xiHat*Vol + V2)/((Vol - xiHat)*(xiHat + 2*Vol)*(xiHat + 2*Vol)) );
  
  double randpm = -1.0;
  //  if(rand() < 0.5) randpm = 1.0;
  //  else randpm = -1.0;
  double vnp2 = expKhalerPot * AbsA1 * AbsA1 * exp(-2.0*a1*tau1) *
    (
     ( -4.0*a1*a1*randpm*sqrt(2*tau1*k111)*(Vol + xiHat/2.0) ) + 
     ( 2.0*a1*tau1*2.0*a1*tau1*(4.0*Vol - xiHat)/(4.0*(Vol - xiHat)) ) + 
     2.0*a1*tau1*(4.0*xi2 + xiHat*Vol + 4.0* V2)/((Vol - xiHat)*(xiHat + 2*Vol)) + 
     3.0*xiHat*(xi2 + 7.0*xiHat*Vol + V2)/((Vol - xiHat)*(xiHat + 2*Vol)*(xiHat + 2*Vol))
     );

  double V = valpp + vnp1 + vnp2;

  return V;
}
