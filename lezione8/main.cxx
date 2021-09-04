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
#include <iomanip>


using namespace std;


double error(double av, double av2, double n);
double estimate_energy(double mu, double sigma, int n_steps);
void print_energy_prog_av( Metropolis M, int n_blocks, int block_size, const char* out_file, FunzioneBase *f);
void test_Metropolis(void);
void energies_mu_sigma(void);
void ground_state_prog_ave(int opt);
void plot_psi2();

double best_mu=0.8113;
double best_sigma=0.6168;
//0.82855241 0.64951924

int main(){
  //test metropolis, check acceptance average
  //test_Metropolis();

  //ex_01
  //energies_mu_sigma();

  //ex_02.1
  ground_state_prog_ave(0);

  //ex_02.2
  //plot_psi2();

  //using fit parameters
  //best_mu=0.82855241; best_sigma=0.64951924; 
  //ground_state_prog_ave(1);
  
  return 0;
}






//_______________________funzioni
void test_Metropolis(void)
{
  //check that metropolis works and set parameters
  int eq_steps=1000000;
  string file = "test.dat";
  double _mu=0.8;
  double _sigma=0.6;
  FunzioneBase *psi2_trial = new psi_trial_squared(_mu, _sigma);
  double x0=1.;
  double unif_step=2.5;
  Metropolis Metro(x0, unif_step, psi2_trial);

  cout << "Metropolis test, parameters: mu=" << _mu << ",          sigma=" << _sigma << ",          starting point x0=" << x0 << ",             unif.step=" << unif_step << endl;  
  ofstream test_file(file);
  double sum=0.;
  for(int i=0; i<eq_steps; i++)
    {
      Metro.Metropolis_step();
      test_file<<Metro.Get_x()<<endl;
      sum+=Metro.Get_A_prob();
    }

  cout<<"average acceptance "<<sum/eq_steps<<endl<<endl;;
}


//_________________________________________________________________________________________________________

void energies_mu_sigma(void)
{
  ofstream en_m_s;
  en_m_s.open("energy_mu_sigma.dat");
  int wd=20;
  
  int n_points=100;
  int n_est=100000;
  double sigma0=0.6; double sigma1=0.63; double sigma_delta=sigma1-sigma0; double d_s=sigma_delta/n_points;
  double mu0=0.79; double mu1=0.82; double mu_delta=mu1-mu0;  double d_m=mu_delta/n_points;
  cout << "estimating ground state energies at different mu and sigma: " << n_points << "^2 points,          mu_interval=["<<mu0<<","<<mu1<<"],       sigma interval=["<<sigma0<<","<<sigma1<<"],                estimations per point=" << n_est <<endl;
  
  
  double sigma; double mu;
  for(int i=0; i<n_points; i++)
    {
      sigma=sigma0 + i*d_s;

      for(int j=0; j<n_points; j++)
	{
	  mu=mu0 + j*d_m;
	  en_m_s << setw(wd) << mu << setw(wd) << sigma << setw(wd) << estimate_energy(mu, sigma, n_est) << endl;
	}
      if(((i+1)*n_points)%1000==0) cout << (i+1)*n_points << " points done"<<endl;
    }
  en_m_s.close();
  cout << "ground state energies <H>(mu,sigma) printed" << endl << endl;
}



double estimate_energy(double mu, double sigma, int n_steps)
{
  FunzioneBase *psi2_trial = new psi_trial_squared(mu, sigma);
  FunzioneBase *f = new H_psi_f_psi(mu, sigma);
  Metropolis Metro(1., 2.5, psi2_trial);     //x0,unifstep,f_prob

  double sum=0.;
  for(int i=0; i<n_steps; i++)
    {
      Metro.Metropolis_step();
      double x=Metro.Get_x();
      sum+=f->Eval(x);
    }
  return sum/n_steps;
}


//________________________________________________________________________________________



void ground_state_prog_ave(int opt)
{
  string file_en="energy_best_ave.dat";
  if(opt==1) file_en="energy_ave_fit.dat";
  FunzioneBase *psi2_best = new psi_trial_squared(best_mu, best_sigma);
  FunzioneBase *f_best = new H_psi_f_psi(best_mu, best_sigma);
  Metropolis Metro_best(1., 2.5, psi2_best);
  int n_blocks=1000;   int block_size=10000;

  cout << "computing average and uncertainty of the best ground state energy estimation:       best_mu=" << best_mu << ",          best_sigma=" << best_sigma << ",          n_blocks=" << n_blocks << ",             block_size=" << block_size << endl;  

  if(opt==0) cout << "using parameters minimizing energy estimation" << endl;
  if(opt==1) cout << "using optimized parameters from fit" << endl;
  print_energy_prog_av( Metro_best, n_blocks, block_size, file_en.c_str(), f_best);    //( Metropolis& M, int n_blocks, int block_size, const char* out_file, FunzioneBase *f)
}


//______________________________________________________________________

void plot_psi2(void)
{
  int eq_steps=10000000;
  string file = "psi2.dat";
  FunzioneBase *psi2 = new psi_trial_squared(best_mu, best_sigma);
  double x0=1.;
  double unif_step=2.;
  Metropolis Metro(x0, unif_step, psi2);

  ofstream psi2_file(file);
  for(int i=0; i<eq_steps; i++)
    {
      Metro.Metropolis_step();
      psi2_file<<Metro.Get_x()<<endl;
    }
}




void print_energy_prog_av( Metropolis M, int n_blocks, int block_size, const char* out_file, FunzioneBase *f)  //class as argument function w/ or w/out &
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
  //int L=int(n_steps/block_size);
  
  for(int i=0; i<n_blocks; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	{
	  M.Metropolis_step();
	  double x=M.Get_x();
	  sum+=f->Eval(x);
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

  cout << "progressive average of ground state best estimation printed" << endl;
  file.close();
}


double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}
