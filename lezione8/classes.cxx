#include "classes.h"
#include "random.h"
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;


Metropolis :: Metropolis(double x0, double unif_passo, const FunzioneBase* f)
{
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  
  m_Random.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  m_x0=x0;
  m_x =x0;
  m_unif_step=unif_passo;
  m_f=f;
}


void Metropolis :: Metropolis_step()
{
  double new_x = m_x + m_Random.Rannyu(- m_unif_step, m_unif_step);

  double q_prob = m_f->Eval(new_x)/m_f->Eval(m_x);
  double A_prob = min(1., q_prob);
  m_A_prob=A_prob;

  double r = m_Random.Rannyu();
      
  if(r<A_prob){ m_x=new_x; }
}


void Metropolis :: Equilibrate(int n_steps)
{
  for(int i=0; i<n_steps; i++)
    {
      Metropolis_step();
    }
}




//_____________________________________________________
/*
double psi_trial :: Eval(double x) const
{
  double a_exp= 0.5*pow((x-mu)/sigma, 2);
  double b_exp= 0.5*pow((x+mu)/sigma, 2);
 
  return exp(-a_exp)+exp(-b_exp);
}
*/

double psi_trial_squared :: Eval(double x) const
{
  double minus_x=pow((x-m_mu)/m_sigma,2);
  double plus_x=pow((x+m_mu)/m_sigma,2);
  double exp_minus=exp(-0.5*minus_x);
  double exp_plus=exp(-0.5*plus_x);
 
  return pow(exp_minus + exp_plus ,2);
}




double H_psi_f_psi :: Eval(double x) const
{
  double minus_x=pow((x-m_mu)/m_sigma,2);
  double plus_x=pow((x+m_mu)/m_sigma,2);
  double exp_minus=exp(-0.5*minus_x);
  double exp_plus=exp(-0.5*plus_x);
  double coeff_minus=  (pow(m_mu-x,2) - pow(m_sigma,2))/pow(m_sigma,4);
  double coeff_plus=  (pow(m_mu+x,2) - pow(m_sigma,2))/pow(m_sigma,4);
  
  
  double Tpsi =-0.5 * ( coeff_minus * exp_minus + coeff_plus * exp_plus ) ;

  double psi=exp_minus+exp_plus;
  double V_x=pow(x,4) -2.5*pow(x,2);
  double Vpsi = V_x*psi;
 
  return (Tpsi + Vpsi)/psi;
}


