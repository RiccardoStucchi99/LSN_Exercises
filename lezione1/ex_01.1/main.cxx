#include "random.h"
#include <iostream>
#include <fstream>
#include "TH1F.h"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class FunzioneBase{ public:  virtual double Eval(double x) const =0;};
class x_function:public FunzioneBase{ public:  double Eval(double x) const{ return x;};};
class sigma_function:public FunzioneBase{ public: double Eval(double x) const{ return pow(x-0.5,2);};};
double error(double av, double av2, double n);
void prog_av_unc(Random myRandom, int block_size, int n_blocks, FunzioneBase *f, const char* out_file);

/*_________________________________________
COMMENT TO THE CODE:
the main task of the exercise, which is the data blocking method, is performed by means of a function "prov_av_unc" that takes as arguments a random generator and a class that basically tells how to 
use the random values generated.
For the ex_01.1.3 the root libraries are exploited. Instead of creating arrays of bins, the histo class is used and the data are automatically binned.  
  ___________________________________________*/

int main (){

  string ex_01="ex_01.1.1.dat";
  string ex_02="ex_01.1.2.dat";
  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());
  FunzioneBase *r=new x_function();
  FunzioneBase *sigma=new sigma_function();

  int n_throws=100000;
  int block_size=1000;
  int L=int(n_throws/block_size);            //n_blocks

  ofstream test("test.dat");
  for(int i=0;i<n_throws;i++)
    test<<myRandom.Rannyu()<<endl;
  test.close();
  
  //punto1
  prog_av_unc(myRandom, block_size, L, r, ex_01.c_str());
  cout<<"progressive average and uncertainty of <r> printed"<<endl;

  //punto2
  prog_av_unc(myRandom, block_size, L, sigma, ex_02.c_str());
  cout<<"progressive average and uncertainty of <(r-0.5)^2> printed"<<endl;
  
 
  //_______________________________________________________________________________________
  Random myRandom3;
  myRandom3.MySetRandom(infile_primes.c_str(), infile_seed.c_str());
  ofstream punto3("ex_01.1.3.dat");

  int n_bins=100;
  int n_throws_per_cycle=10000;
  int n_cycles= 100;
  double exp_ev= n_throws_per_cycle/n_bins;
  //int L=int(n_throws/block_size);

  vector<double> chi_squared;
  
  for(int j=0; j<n_cycles; j++)
    {
      TH1F histo("histo", "histo", 100, 0, 1);

      for(int i=0; i<n_throws_per_cycle; i++)
	{
	  double x=myRandom.Rannyu();
	  histo.Fill(x);
	}
      

      
      double chi=0;

      for(int i=0; i<n_bins; i++)
	{
	  chi+= pow(histo.GetBinContent(i+1)-exp_ev, 2) / double(exp_ev);
	}
  
      chi_squared.push_back(chi);
      punto3 << chi << endl;
   }

  cout<<"chi^2 values printed"<< endl;

  return 0;

}





void prog_av_unc(Random myRandom, int block_size, int n_blocks, FunzioneBase *f, const char* out_file)
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
 
  for(int i=0; i<n_blocks; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	sum+=f->Eval(myRandom.Rannyu());

      double block_av = sum/double(block_size);
      double block_av2 = block_av*block_av;
	
	
      prog_sum += block_av;
      prog_sum2 += block_av2;

      double prog_av = prog_sum / (i+1);
      double prog_av2 = prog_sum2 / (i+1);
      double prog_error = error(prog_av,  prog_av2, i);
      file << prog_av << " " << prog_error << endl;
    }

  file.close();
}



double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}
