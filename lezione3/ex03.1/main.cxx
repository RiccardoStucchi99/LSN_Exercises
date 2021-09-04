#include "math.h"
#include "random.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>


using namespace std;

double error(double av, double av2, double n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ n );
}


class Option{
protected:
  double m_S_0;
  double m_del_time;
  double m_K;
  double m_rate;
  double m_sigma;
  Random m_rnd;
public:
  Option(double S_0, double del_time, double K, double rate, double sigma, Random rnd){m_S_0=S_0, m_del_time=del_time, m_K=K, m_rate=rate, m_sigma=sigma,  m_rnd=rnd;}
  ~Option();
  virtual double Price_eval() =0;
};
  
class Call: public Option{
public:
  Call(double S_0, double del_time, double K, double rate, double sigma, Random m_rnd) : Option(S_0, del_time, K, rate, sigma,  m_rnd){};
  double Price_eval();
};


class Put: public Option{
public:
  Put(double S_0, double del_time, double K, double rate, double sigma, Random m_rnd) : Option(S_0, del_time, K, rate, sigma,  m_rnd){};
  double Price_eval();
};


class Discr_Call: public Option{
public:
  Discr_Call(double S_0, double del_time, double K, double rate, double sigma, Random m_rnd, double n_t) : Option(S_0, del_time, K, rate, sigma,  m_rnd){m_n_t=n_t; m_dt=m_del_time / double(n_t);};
  double Price_eval();

private:
  double m_n_t;
  double m_dt;
};


class Discr_Put: public Option{
public:
  Discr_Put(double S_0, double del_time, double K, double rate, double sigma, Random m_rnd, double n_t) : Option(S_0, del_time, K, rate, sigma,  m_rnd){m_n_t=n_t; m_dt=m_del_time / double(n_t);};
  double Price_eval();

private:
  double m_n_t;
  double m_dt;
};


void prog_av_unc(int block_size, int n_blocks, Option *opt, const char* out_file);








  
int main (){

  Random myRandom;
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  myRandom.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  string call_file="call.dat";
  string d_call_file="discr_call.dat";
  string put_file="put.dat";
  string d_put_file="discr_put.dat";
  
  double S_0 = 100.;
  double del_time = 1.;
  double K = 100.;
  double rate = 0.1;
  double sigma = 0.25;
  double n_t=100;
  
  int n_throws=100000;
  int block_size=1000;
  int n_blocks= int(n_throws/block_size);

   
  Option *my_call = new Call(S_0, del_time, K, rate, sigma, myRandom);
  prog_av_unc(block_size, n_blocks, my_call, call_file.c_str());
  cout << "call price printed"  << endl;
  
  Option *my_put = new Put(S_0, del_time, K, rate, sigma, myRandom);
  prog_av_unc(block_size, n_blocks, my_put, put_file.c_str());
  cout << "put price printed"  << endl;

  Option *my_call_discr = new Discr_Call(S_0, del_time, K, rate, sigma, myRandom, n_t);
  prog_av_unc(block_size, n_blocks, my_call_discr, d_call_file.c_str());
  cout << "discrete call price printed"  << endl;

  Option *my_put_discr = new Discr_Put(S_0, del_time, K, rate, sigma, myRandom, n_t);
  prog_av_unc(block_size, n_blocks, my_put_discr, d_put_file.c_str());
  cout << "discrete put price printed"  << endl;
    
  return 0;

}




double Call::Price_eval()
{
  double z = m_rnd.Gauss(0,1);
  double S_T = m_S_0 * exp( (m_rate-m_sigma*m_sigma/2.)*m_del_time  + m_sigma*z*sqrt(m_del_time) );
  double C = exp(-m_rate*m_del_time)*max(0., S_T - m_K);

  return C;
}



double Discr_Call::Price_eval()
{
  double S_i=m_S_0;
   
  for(int k=0; k<m_n_t; k++)
    {
      double z = m_rnd.Gauss(0,1);
      S_i *= exp( (m_rate-m_sigma*m_sigma/2.)*m_dt  + m_sigma*z*sqrt(m_dt) );
    }
      
  double C = exp(-m_rate*m_del_time)*max(0., S_i - m_K);

  return C;
}




double Put::Price_eval()
{
  double z = m_rnd.Gauss(0,1);
  double S_T = m_S_0 * exp( (m_rate-m_sigma*m_sigma/2.)*m_del_time  + m_sigma*z*sqrt(m_del_time) );
  double P = exp(-m_rate*m_del_time)*max(0., -S_T + m_K);

  return P;
}



double Discr_Put::Price_eval()
{
  double S_i=m_S_0;
  
  for(int k=0; k<m_n_t; k++)
    {
      double z = m_rnd.Gauss(0,1);
      S_i *= exp( (m_rate-m_sigma*m_sigma/2.)*m_dt  + m_sigma*z*sqrt(m_dt) );
    }
      
  double P = exp(-m_rate*m_del_time)*max(0., -S_i + m_K);

  return P;
}



void prog_av_unc(int block_size, int n_blocks, Option *opt, const char* out_file)
{
  ofstream file(out_file);
  double prog_sum=0;
  double prog_sum2=0;
 
  for(int i=0; i<n_blocks; i++)
    {
      double sum=0;
      
      for(int j=0; j<block_size; j++)
	sum+=opt->Price_eval();

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
