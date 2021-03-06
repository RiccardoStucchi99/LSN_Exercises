#include "classes.h"
#include "random.h"
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;


Metropolis :: Metropolis(const vector<double> &x0, double unif_passo, const FunzioneBase* f)
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
  vector<double> new_x = {m_x[0] + m_Random.Rannyu(- m_unif_step, m_unif_step),    m_x[1] +  m_Random.Rannyu(- m_unif_step, m_unif_step),        m_x[2] +  m_Random.Rannyu(- m_unif_step, m_unif_step)};

  double q_prob = m_f->Eval(new_x)/m_f->Eval(m_x);
  double A_prob = min(1., q_prob);
  m_A_prob=A_prob;

  double r = m_Random.Rannyu();
      
  if(r<A_prob){ m_x=new_x; }
}


void Metropolis :: Metropolis_step_Gauss()
{
  vector<double> new_x = {m_x[0] + m_Random.Gauss(0,m_unif_step),    m_x[1] +  m_Random.Gauss(0, m_unif_step),        m_x[2] +  m_Random.Gauss(0, m_unif_step)};

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


double Metropolis :: Get_radius()
{
  return sqrt(inner_product( m_x.begin(), m_x.end(), m_x.begin(), 0.));
}



//_____________________________________________________

double psi_100_squared :: Eval(const vector<double> &x) const
{
  double r=sqrt(inner_product( x.begin(), x.end(), x.begin(), 0.));

  double N = 1/acos(-1);

  return N*exp(-2*r);
}



double psi_210_squared :: Eval(const vector<double> &x) const
{
  double r=sqrt(inner_product( x.begin(), x.end(), x.begin(), 0.));

  double N = 1/(32.*acos(-1));

  return N*exp(-r)*pow(x[2],2);
}
