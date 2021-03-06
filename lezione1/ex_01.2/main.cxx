#include "random.h"
#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <iterator>
#include <iomanip>
#include "math.h"


using namespace std;

double error(double av, double av2, double n)
{ 
  return sqrt( (av2 - av*av)/ n );
}


int main (int argc, char *argv[]){


  //_____________________________________________________Random generator__________________________________________________
  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  //________________________________________________________________
  int n_distr[]={1,2,10,100};
  ofstream Dice("ex_01.2_dice.dat");
  ofstream Exp("ex_01.2_exp.dat");
  ofstream Cauchy("ex_01.2_cauchy.dat");
  
  int n_throws = 10000;

 
  for(int j=0; j<n_throws; j++)
    {
      for(int k=0; k<4; k++)
	{
	  double sum_d=0;
	  double sum_e=0;
	  double sum_c=0;
	       
	  for(int h=0; h<n_distr[k]; h++)
	    {
	      double x=myRandom.Dice();
	      double y=myRandom.Exp(1);
	      double z=myRandom.Cauchy(0,1);
	      sum_d+=x;
	      sum_e+=y;
	      // cout<<y<<endl;
	      sum_c+=z;
	    }
	  
	  Dice << sum_d/n_distr[k] << " ";
	  Exp << sum_e/n_distr[k] <<  " ";
	  Cauchy << sum_c/n_distr[k] <<  " ";
	}

      Dice << endl;
      Exp <<  endl;
      Cauchy <<  endl;
    }

  cout<<"distributions printed and code executed"<<endl;
   return 0;

}
