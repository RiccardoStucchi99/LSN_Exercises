#include "random.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <iterator>
#include <math.h>
#include <algorithm>
#include "classes.h"


using namespace std;


double error(double av, double av2, double n);
void print_radius_prog_av( Metropolis& M, int n_steps, int block_size, const char* out_file);
void print_radius_prog_av_Gauss( Metropolis& M, int n_steps, int block_size, const char* out_file);

//___________________________________________________
// CODE DESCRIPTION
// the main code is divided into two parts: test and averages computation. After initializing the variables (note that Metropolis is implemented through a class, its details such as members and methods
// can be found in the classes.* files) a test run of test_steps is made to check the main properties of the Metropolis. By changing the metropolis step method called (gauss or not), two different
// transition probabilities are tried.
// the second part consist of only two functions that compute progressive averages and uncertaintes. Also in this case, if we change the functions name by adding Gauss, a different step transition
// probability is used
//_______________________________________________


int main(){

 
  

  // ofstream file2("ex_05.1.psi_210.dat");
  // ofstream file3("ex_05.1.positions.dat");
  
  int n_steps=1000000;          //*10 gauss
  int n_blocks=100;
  int block_size=n_steps/n_blocks;


  string file_100 = "ex_05.1.psi_100.dat";
  FunzioneBase *psi_100 = new psi_100_squared;
  double w=1.;
  vector<double> x0(3, w);
  double unif_step_100=1.;                                     //1 for gauss and uniform
  Metropolis M_100(x0, unif_step_100, psi_100);
  double sum_acc=0.;


  string file_210 = "ex_05.1.psi_210.dat";
  FunzioneBase *psi_210 = new psi_210_squared;
  double q=1.5;
  vector<double> p0(3, q);
  double unif_step_210=2;                                   // 2 for uniform
  Metropolis M_210(p0, unif_step_210, psi_210);
  double sum_acc_2=0.;



  ofstream block_file("ex_05.1.block_test.dat");

  int test_steps=1000000;
  for(int j=0; j<test_steps; j++)
    {
      vector<double> x_100=M_100.Get_x();
      vector<double> x_210=M_210.Get_x();
      
      M_100.Metropolis_step();                            //gauss_step
      M_210.Metropolis_step();                            //gauss_step
      
      copy(x_100.begin(), x_100.end(), ostream_iterator<double>(block_file, " "));
      block_file << M_100.Get_radius() << " ";

      copy(x_210.begin(), x_210.end(), ostream_iterator<double>(block_file, " "));
      block_file << M_210.Get_radius() << endl;


      sum_acc+=M_100.Get_A_prob();
      sum_acc_2+=M_210.Get_A_prob();
    } 


  cout << "average acceptance probability of psi100:   " << sum_acc/double(test_steps) << endl;
  cout << "average acceptance probability of psi210:   " << sum_acc_2/double(test_steps) << endl;
  
  //_____________________________________________________________________________________________



  print_radius_prog_av(M_100, n_steps, block_size, file_100.c_str());            //gauss
  print_radius_prog_av(M_210, n_steps, block_size, file_210.c_str());
  cout << "progressive averages and uncertaintes printed"<<endl;

  
  return 0;
}






//_______________________funzioni
double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}


void print_radius_prog_av_Gauss( Metropolis& M, int n_steps, int block_size, const char* out_file)
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
  int L=int(n_steps/block_size);
  
  for(int i=0; i<L; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	{
	  M.Metropolis_step_Gauss();
	  sum+=M.Get_radius();
	} 
	 
      double block_av = sum/double(block_size);
      double block_av2 = block_av*block_av;
      	
      prog_sum += block_av;
      prog_sum2 += block_av2;
            
      double prog_av = prog_sum / (i+1);
      double prog_av2 = prog_sum2/ (i+1);
      double prog_error = error(prog_av,  prog_av2, i);
      
      file << prog_av << " " << prog_error << endl;
      //cout << i+1 << " blocks completed" << endl;
    }

  file.close();
}


void print_radius_prog_av( Metropolis& M, int n_steps, int block_size, const char* out_file)
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
  int L=int(n_steps/block_size);
  
  for(int i=0; i<L; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	{
	  M.Metropolis_step();
	  sum+=M.Get_radius();
	} 
	 
      double block_av = sum/double(block_size);
      double block_av2 = block_av*block_av;
      	
      prog_sum += block_av;
      prog_sum2 += block_av2;
            
      double prog_av = prog_sum / (i+1);
      double prog_av2 = prog_sum2/ (i+1);
      double prog_error = error(prog_av,  prog_av2, i);
      
      file << prog_av << " " << prog_error << endl;
      //cout << i+1 << " blocks completed" << endl;
    }

  file.close();
}
