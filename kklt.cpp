// first make dlib
// mkdir build
// cd build
// cmake ..
// cmake --build . --config Release
// 


// mpicxx -std=c++11 -O3 -I./dlib-19.17/ kklt.cpp ./dlib-19.17/dlib/all/source.cpp -lpthread -lX11 -o a_kklt 


#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include <iostream>

#include <complex.h>
#include <math.h>
#include <cstdlib>

#define NOUTPUT 3   // the outputs: V0, tau4

using namespace std;
using namespace dlib;

typedef matrix<double,0,1> column_vector;

// ----------------------------------------------------------------------------------------
double kklt(double tau1);
double kkltup(double tau4, double );
double kkltup1par(double tau1);

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


    // (3) KKLT with uplift term; one parameters
  auto complex_holder_table = [theta](double x0){
    //  auto complex_holder_table = [theta](double x0){ // x0=tau4
    //    return kklt(x0);
    return kkltup1par(x0);  
// ./a_kklt 1 1 1 1 1 1 1 
// 2.5835963e-07 20.4419235
  };

  double pi = 4* atan(1.0);
  auto result = find_min_global(complex_holder_table, 
				{15}, // lower bounds
				{50}, // upper bounds 
				max_function_calls(500));


      column_vector take_results = result.x; 
      cout.precision(9);
      cout << result.y << " " << result.x(0,0) << endl << endl;


    // (2) KKLT with uplift term; two parameters
  auto complex_holder_table2 = [theta](double x0, double x1){
    return kkltup(x0, x1);
// ./a_kklt 1 1 1 1 1 1 1 
// 2.5835963e-07 20.4419234 -1.48109985e-07
// rho1 should be fixed at 0, so here -1.48109985e-07 is accepted as zero
  };


  auto result2 = find_min_global(complex_holder_table2, 
  				{15, -3.14 }, // lower bounds
  				{50.0, 3.14}, // upper bounds 
  				max_function_calls(500));
      cout.precision(9);
      cout << result2.y << " " << result2.x(0,0) << " " << result2.x(0,1) << endl << endl;


      //     (1) STANDARD KKLT
  auto complex_holder_table3 = [theta](double x0){ // x0=tau4
    return kklt(x0);
  };


  auto result3 = find_min_global(complex_holder_table3, 
  				{0,  }, // lower bounds
  				{200.0}, // upper bounds 
  				max_function_calls(500));
      cout.precision(9);
      cout << result3.y << " " << result3.x(0,0) << endl << endl;




    return 0; 
}

double kklt(double tau1){
  return (0.15000000000000002*exp(0.1 - 0.2*tau1)*(3 - 0.0003*exp(0.1*tau1) + 0.1*tau1))/pow(tau1,2);

}

double kkltup1par(double tau1){

  return (0.11051709180756478*((22.44075251576936*(29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/
			       ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) - 
			       (2.35*((16.254053429258022*(29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/
				      ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) + 
				      (0.15*tau1*(117.41966794718861 + 1.1422184454783528*pow(tau1,1.5) + 0.17777777777777778*pow(tau1,3)))/
				      ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5)))))/
			       pow(2.718281828459045,0.15*tau1) + (-0.28460498941515416*sqrt(tau1)*(2.7090089048763373 + 0.21081851067789195*pow(tau1,1.5)) + 
								     (0.0225*pow(tau1,2)*(-5.4180178097526746 + 0.8432740427115678*pow(tau1,1.5)))/
								     (-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5)) + 
								     (16.254053429258022*(29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/
								     ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) + 
								     (0.3*tau1*(117.41966794718861 + 1.1422184454783528*pow(tau1,1.5) + 0.17777777777777778*pow(tau1,3)))/
								     ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5))))/
			       pow(2.718281828459045,0.3*tau1)))/pow(2.7090089048763373 + 0.21081851067789195*pow(tau1,1.5),2);

  }

double kkltup(double tau1, double rho1){
  double V =  (0.11051709180756478*((22.44075251576936*(29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/ ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) + (-0.28460498941515416*sqrt(tau1)*(2.7090089048763373 + 0.21081851067789195*pow(tau1,1.5)) +           (0.0225*pow(tau1,2)*(-5.4180178097526746 + 0.8432740427115678*pow(tau1,1.5)))/           (-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5)) +           (16.254053429258022*(29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/           ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) +           (0.3*tau1*(117.41966794718861 + 1.1422184454783528*pow(tau1,1.5) + 0.17777777777777778*pow(tau1,3)))/           ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5))))/        pow(2.718281828459045,0.3*tau1) - (2.35*((16.254053429258022*               (29.354916986797154 + 7.995529118348469*pow(tau1,1.5) + 0.044444444444444446*pow(tau1,3)))/             ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*pow(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5),2)) +             (0.15*tau1*(117.41966794718861 + 1.1422184454783528*pow(tau1,1.5) + 0.17777777777777778*pow(tau1,3)))/             ((-5.4180178097526746 + 0.21081851067789195*pow(tau1,1.5))*(5.4180178097526746 + 0.4216370213557839*pow(tau1,1.5))))*          cos(0.15*rho1))/pow(2.718281828459045,0.15*tau1)))/pow(2.7090089048763373 + 0.21081851067789195*pow(tau1,1.5),2); 

  return V;
}
