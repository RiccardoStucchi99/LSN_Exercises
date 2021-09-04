/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <vector>
#include <string>
using namespace std;

string folder;

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;
//string s1="epot"; string s2="ekin"; string s3="temp"; string s4="etot";
vector<string> file_title = {"epot", "ekin", "temp", "etot"};

// averages
vector<double> obs_sums(4, 0.);
vector<double> prog_sum(4, 0.);
vector<double> prog_sum2(4, 0.);
double acc,att;
int block_index=0;
int block_size=100;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut,target_temp;

// simulation
int nstep, iprint, seed; int nstep_in;
double delta;

//functions
void Input(void);
void Restart_Input(void);
void Move(void);
void PrintConf(const char*);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Print_Averages(void);
void Update_Averages_sums(void);
double error(double, double, int);
double Compute_gofr_coeff(double r);




//gofr measurement
const int nbins=100;
double bin_size;

vector<double> gofr_sums(nbins, 0.);
vector<double> stima_gofr(nbins, 0.);
vector<double> gofr_prog_sum(nbins, 0.);
vector<double> gofr_prog_sum2(nbins, 0.);
vector<double> g(nbins, 0.);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
