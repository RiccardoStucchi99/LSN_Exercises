#include "random.h"
#include <iostream>
#include <float.h>
#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <iterator>
#include "math.h"


using namespace std;

double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}

//_______________________________________________________
// COMMENT:
// The main problem related to the codification of the buffon's algorythm arises when randomly generating the needle direction. In fact, the easiest way to go, namely uniform(0,pi) is inconsistent
// from a logical standpoint: we would be computing pi by already using pi in the algorythm. We hence need a pi independent angle generator. The adopted method consists in randomlu generating
// (uniformly) a point in a 1*1 square. When the distance from the origin is smaller than 1, the angle is given as acos(x/sqrt(x2+y2)) (equivalently we can directly use the cosine). However, there is
// a second possible solution, even if it is less rigorous. It simply consists in uniformly generating a number from the largest possible interval (ideally infinite). This second method exploits the
// periodicity of angles and the fact that even if the maximum interval width is not an integer multiple of pi (indeed it is not), the remainder of width/pi would be negligible with respect to the
// integer part of width/pi
//________________________________________________________



int main (int argc, char *argv[]){

 //_____________________________________________________Random generator__________________________________________________
  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  ofstream buffon("Buffon.dat");
  
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  //___________________________________________________________________________________________________________________________
  

  int n_per_cycle=10000;
  int cycles=1000;
  double L=0.08;
  double D=0.1;

  double prog_sum=0;
  double prog_sum2=0;
 
  for(int i=0; i<cycles; i++)
    {
      int n_inter=0;
  
      for(int j=0; j<n_per_cycle; j++)
	{
	  double x_c=myRandom.Rannyu();
	  double y_c=myRandom.Rannyu();
	  //double teta=DBL_MAX*myRandom.Rannyu();
	  double teta=myRandom.Unif_0_PI();

	  double l=fabs(L*sin(teta));

	  if( fmod(y_c, D) < l/2 ||  fmod(y_c, D) > D - l/2 )
	    n_inter++;
	}

      double pi=2*L*double(n_per_cycle)/(double(n_inter)*D);
      double pi2= pi*pi;
      
      prog_sum += pi;
      prog_sum2 += pi2;


      double prog_av = prog_sum / (i+1);
      double prog_av2 = prog_sum2 / (i+1);
      double prog_error = error(prog_av,  prog_av2, i);

      buffon << prog_av << " " << prog_error << endl; 

     }
  cout << "pi code succesfully executed!" << endl;
 
  return 0;

}
  


  
 
