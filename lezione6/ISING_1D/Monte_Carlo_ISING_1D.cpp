/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Monte_Carlo_ISING_1D.h"
#include <vector>
 
using namespace std;


/* COMMENTS ABOUT THE CODE
To test the code several input commands are available. However, in order to complete the required exercise we need only two of them plus some changes to the input file. Here are the required operations for ex06.1:
1) execute  ./clean.sh to remove any previous simulations data;
2) set the input file with a combination (algorythm - h) chosen among the followings: (1-0), (1-0.02), (0-0), (0-0.02), Where 0 stands for Gibbs, 1 for Metropolis. The other input data do not need to be modified;
3) check how many equilibration steps are needed: ./exe 0 n__equil_steps   (25 are enough). A file with the instantaneous values of internal energy will be created at Metropolis/h=0. THe temperatuere in given by the input (0.8 seems fine).
4) start simulations at different temperatures. (./exe 3 n__equil_steps) The code is written so that there is no need to modify the input temperature. Moreover the equilibration is automatically done, provided a proper number of eq. steps are provided (as said before 25 equil_steps are enough) and there is no need to restart the system.
5)repeat steps 2 to 5 with a different combination (algorythm - h)
*/


int main(int argc, char** argv)
{
  if(argc != 3)
    {
      cerr << " for " << argv[0] << " first argument: <0> only equil, <1> single run with restart, <2> single run, <3> temperature run; second argument:< > equilibration steps (<-> no equilibration)" << endl;
      return 1;
    }

  prog_option=atof(argv[1]);
  
  if(prog_option==0){ Input(0,0.); Equilibration(atof(argv[2]));};                               //only equilibration
  
  if(prog_option==1){ Input(1,0.); Run_Simulation();};                                           //restart from an already equilibrated config

  if(prog_option==2){ Input(0,0.); Equilibration(atof(argv[2])); Run_Simulation();};             //full run, equilibration+simulation

  if(prog_option==3)
    {
      cout<<"Temperature run at h="<<h<<"in progress"<<endl;
      double T0=0.5; double T1=2.0; int n_d=15; double d_T=(T1-T0)/double(n_d);
      for(int i=0; i<=n_d; i++)
	{
	  double T=T0+i*d_T;
	  temp_str="T"+to_string(i);
	  Input(0,T);
	  Equilibration(atof(argv[2]));
	  Run_Simulation();
	}
      cout<<"Temperature run at h="<<h<<" completed"<<endl;
    }

  return 0;
}





void Input(int restart_opt, double T)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  if(T!=0.) temp=T; 
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
  if(h==0.) h_str="h=0";
  else h_str="h=0.02";
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  if(metro==1) alg_str="Metropolis";
  else alg_str="Gibbs";
  
  ReadInput >> nblk;

  ReadInput >> nstep;

  //folder instruction
  folder=alg_str+"/"+h_str+"/";
  

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  ix = 2; //Magnetization
  im = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables


  if(restart_opt==0){
    for (int i=0; i<nspin; ++i)
      {
	if(rnd.Rannyu() >= 0.5) s[i] = 1;
	else s[i] = -1;
      }
  }
  else if(restart_opt==1){
    ifstream ReadConfig(folder+"config.final");
    if(!ReadConfig.good()) {cerr << "error in opening config.final"<<endl; exit(-1);};
    for (int i=0; i<nspin; ++i)
      {
	ReadConfig >> s[i];
      }
  }
  else cerr<<"input option error"<<endl;
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}



void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;
 
  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    sm=s[o];
    energy_old=Boltzmann(sm, o);
    energy_new=Boltzmann(-sm, o);
    double delta_energy=-energy_old+energy_new;                               
    attempted++;
    
    if(metro==1) //Metropolis
    {
      double a=min(1., exp(-beta*delta_energy));
      double r=rnd.Rannyu();
      if(r<a) {s[o]=-sm; accepted++;} //cout<<"flipped"<<endl;}
	// else {cout<<"non flipped" << endl;}
    }
    else //Gibbs sampling
    {
      accepted++;                                     //Gibbs sampling has a=1
      double p=1./(1.+exp(beta*delta_energy));                           //flipping probability        
      double r=rnd.Rannyu();
      if(r<p) s[o]=-sm;
    }
  }
}



void Equilibration(int eq_step)
{
  ofstream inst_ene;  ofstream inst_mag;
  if(prog_option==0){inst_ene.open(folder+"inst_ene.dat");  inst_mag.open(folder+"inst_mag.dat");};
  
  for(int i=1; i <= eq_step; ++i) //Simulation
  {
    Move(metro);
    Measure();
    if(prog_option==0){	inst_ene << walker[iu] << endl;  inst_mag << walker[im] << endl;   }
  }
  
  if(prog_option==0) ConfFinal(); //Write final configuration

  cout << "Equilibration completed, equilibration steps="<<eq_step<<endl;
}


void Run_Simulation(void)
{
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
	{
	  Move(metro);
	  Measure();
	  Accumulate(); //Update block averages
	}
      Averages(iblk);   //Print results for current block
    }
   //ConfFinal(); //Write final configuration
}





double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
// INCLUDE YOUR CODE HERE
    m += s[i];
  }
  walker[iu] = u;                          
// INCLUDE YOUR CODE HERE
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{
   for(int i=0; i<4; i++)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    if(metro==1) cout << "Acceptance rate " << accepted/attempted << endl << endl;

    /*
     string folder="h=0";
    Ene.open(temp_string+"output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per particle
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();
    */
    

    
    stima_u = blk_av[iu]/blk_norm;                                            //energy
    stima_c = beta*beta * (blk_av[ic]/blk_norm - stima_u*stima_u);            //heat capacity
    stima_m = blk_av[im]/blk_norm;                                            //magnetization
    stima_x = beta * blk_av[ix]/blk_norm;                                   //susceptibility
    vector<string> obs_string = {"ene", "heat", "chi", "mag"};
    vector<double> obs_stima = {stima_u, stima_c, stima_x, stima_m};
    vector<double> obs_err = {err_u, err_c, err_x, err_m};
    int idx_0=0; int idx_f=3; 
    if(h!=0) {idx_0=3; idx_f=4;};

    
    ofstream obs_file;
    obs_file.open(folder+temp_str+".dat",ios::app);
    if(!obs_file.good()) {cerr << "error in opening ob_file"<<endl;exit(-1);};

    obs_file << setw(wd) << iblk;
    for(int i=idx_0; i<idx_f; i++)
      {	
	glob_av[i]  += obs_stima[i];
	glob_av2[i] += obs_stima[i]*obs_stima[i];
	obs_err[i]=Error(glob_av[i],glob_av2[i],iblk);
	obs_file <<  setw(wd) << obs_stima[i] << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << obs_err[i];
      }
    obs_file<<endl;
    obs_file.close();
   
    cout<< "data printed in " <<folder+temp_str+".dat"<<endl;
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(folder+temp_str+"config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
