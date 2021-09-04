#include "math.h"
#include "random.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>


using namespace std;

double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}

double compute_displ(vector<double> p)
{
  double sum=0;

  for(int j=0; j<3; j++)
      sum+= pow(p[j],2);

  return sum;
}
      

vector<double> block_av_sigma(vector<double> displ)
{
  int N = displ.size();
  int n_blocks = 100;
  int L = N/n_blocks;

  double prog_sum=0;
  double prog_sum2=0;

  for(int i=0; i<n_blocks; i++)
    {
      double sum=0;
      
      for(int j=0; j<L; j++)
	{
	  sum+=displ[L*i + j];
	}

      double block_value = sqrt(sum/double(L));                                          //x variable,          in this case: x= sqrt( < |position|**2 > )
      double block_value2 = block_value*block_value;                                     //x2 variable

      prog_sum+=block_value;
      prog_sum2+=block_value2;
    }

  double average = prog_sum/double(n_blocks);                                            //<x> = sum(x_i)/n_blocks
  double average2 = prog_sum2/double(n_blocks);                                          //<x2> = sum( x2_i ) / n_blocks
  double sigma = error(average, average2, n_blocks-1);                                   //unc on the mean =sqrt[ (<x2> -<x>2)/ n-1 ]
  
  vector<double> results={average, sigma};

  return results;
}

      
	
      




int main (int argc, char *argv[]){


  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());
  ofstream rw_discr("discrete_rw.dat");


  int n_steps=100;
  int n_rw=10000;

  vector<vector<double>> positions( n_rw, vector<double> (3, 0.));
  vector<double> displacements2( n_rw, 0.);
  rw_discr << 0. << " " << 0. << endl;                                    //initial displacement=0
  double passo = 4.;

  for(int i=0; i<n_steps; i++)
    {
      for(int j=0; j<n_rw; j++)
	{
	  vector<double> step(3, 0.);
	  
	  int dim=floor(myRandom.Rannyu(0,3));
	  int verso=myRandom.sign();
	  step[dim]=passo*double(verso);

	  positions[j][dim]+=step[dim];
	  displacements2[j]=compute_displ(positions[j]);
	}

    
      vector<double> res = block_av_sigma(displacements2);
      rw_discr << res[0] << " " << res[1] << endl;
    }
  

  cout << "discrete random walk printed with success" << endl;


  Random myRandom2;
  myRandom2.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  vector<vector<double>> positions_cont( n_rw, vector<double> (3, 0.));
  vector<double> displacements2_cont( n_rw, 0.);
  ofstream rw_cont("continuous_rw.dat");
  rw_cont << 0. << " " << 0. << endl;

  for(int i=0; i<n_steps; i++)
    {
      for(int j=0; j<n_rw; j++)
	{
	  double theta=acos(1 - 2*myRandom2.Rannyu());
	  double phi=myRandom.Rannyu(0., 2*acos(-1));
	  vector<double> step={passo*sin(theta)*cos(phi), passo*sin(theta)*sin(phi), passo*cos(theta)};

	  for(int k=0; k<3; k++) positions_cont[j][k]+=step[k];
	  displacements2_cont[j]=compute_displ(positions_cont[j]);
	}

      vector<double> res = block_av_sigma(displacements2_cont);
      rw_cont << res[0] << " " << res[1] << endl;
    }
  

  cout << "continuous random walk printed with success" << endl;

 
  return 0;

}
