#ifndef __classes_h__
#define __classes_h__

#include "random.h"
#include <math.h>
#include <vector>

using namespace std;

class FunzioneBase
{
 public:
  virtual double Eval(const vector<double> &x) const =0;
};


class psi_100_squared:public FunzioneBase
{
 public:
  double Eval(const vector<double> &x) const;
};

class psi_210_squared:public FunzioneBase
{
 public:
  double Eval(const vector<double> &x) const;
};



class Metropolis
{
 public:
  Metropolis(const vector<double> &x0, double unif_step, const FunzioneBase* f);
  ~Metropolis(){};

  void Metropolis_reset(const vector<double> &x0, double unif_step) { m_x0=x0; m_x=x0; m_unif_step=unif_step;};
  void Metropolis_step();
  void Metropolis_step_Gauss();
  void Equilibrate(int n_steps);
  double Get_A_prob() { return m_A_prob; };
  double Get_radius();
  vector<double> Get_x() { return m_x; };

 private:
  Random m_Random;
  const FunzioneBase *m_f;
  double m_unif_step;
  vector<double> m_x0;                               
  vector<double> m_x;
  double m_A_prob;
};



#endif
