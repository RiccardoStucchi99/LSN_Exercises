#ifndef __classes_h__
#define __classes_h__

#include "random.h"
#include <math.h>
#include <vector>

using namespace std;

class FunzioneBase
{
 public:
  virtual double Eval(double x) const =0;
};


class psi_trial_squared:public FunzioneBase
{
 public:
  psi_trial_squared(double mu, double sigma) { m_mu=mu; m_sigma=sigma;};
  double Eval(double x) const;
  
 private:
  double m_mu;
  double m_sigma;
};


class H_psi_f_psi:public FunzioneBase
{
 public:
  H_psi_f_psi(double mu, double sigma) { m_mu=mu; m_sigma=sigma;};
  double Eval(double x) const;
  
 private:
  double m_mu;
  double m_sigma;
};






class Metropolis
{
 public:
  Metropolis(double x0, double unif_step, const FunzioneBase* f);
  ~Metropolis(){};

  //void Metropolis_reset(const vector<double> &x0, double unif_step) { m_x0=x0; m_x=x0; m_unif_step=unif_step;};
  void Metropolis_step();
  void Metropolis_step_Gauss();
  void Equilibrate(int n_steps);
  double Get_A_prob() { return m_A_prob; };
  double Get_x() { return m_x; };

 private:
  Random m_Random;
  const FunzioneBase *m_f;
  double m_unif_step;
  double m_x0;                               
  double m_x;
  double m_A_prob;
};



#endif
