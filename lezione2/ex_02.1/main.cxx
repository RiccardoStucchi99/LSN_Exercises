#include "math.h"
#include "random.h"
#include <iostream>
#include <string>
#include <fstream> 

using namespace std;

double const pi=acos(-1);
class FunzioneBase{ public:  virtual double Eval(double x) const =0;};
class integranda1:public FunzioneBase{ public:  double Eval(double x) const{ return pi*0.5*cos(0.5*pi*x);};};
class integranda2:public FunzioneBase{ public: double Eval(double x) const{ return (pi*0.5*cos(0.5*pi*x))/(-2*x+2);};};
class unif_sampling_function:public FunzioneBase{ public: double Eval(double x) const{ return x;};};
  class sampling_function_1:public FunzioneBase{ public: double Eval(double x) const{ return 1-sqrt(1-x);};};
double error(double av, double av2, double n);
void prog_av_unc(Random &myRandom, int block_size, int n_blocks, FunzioneBase *f, FunzioneBase *is_func, const char* out_file);



int main (){

  string ex_01="ex_02.1.1.dat";
  string ex_02="ex_02.1.2.dat";
  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  int n_throws=100000;
  int block_size=1000;
  int L=int(n_throws/block_size);            //n_blocks


  //punto1
  FunzioneBase *f1=new integranda1();
  FunzioneBase *p1=new unif_sampling_function();
   
  prog_av_unc(myRandom, block_size, L, f1, p1, ex_01.c_str());
  cout<<"progressive average and uncertainty of Integral computed with uniform distribution"<<endl;

  //punto2
  FunzioneBase *f2=new integranda2();
  FunzioneBase *p2=new sampling_function_1();
   
  prog_av_unc(myRandom, block_size, L, f2, p2, ex_02.c_str());
  cout<<"progressive average and uncertainty of Integral computed with importance sampling"<<endl;

  return 0;
}




  
void prog_av_unc(Random &myRandom, int block_size, int n_blocks, FunzioneBase *f, FunzioneBase *is_func,  const char* out_file)
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
 
  for(int i=0; i<n_blocks; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	{
	  //sum+=f->Eval(myRandom.Rannyu());
	  double x=is_func->Eval(myRandom.Rannyu());
	  sum+=f->Eval(x);
	}
	  
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
