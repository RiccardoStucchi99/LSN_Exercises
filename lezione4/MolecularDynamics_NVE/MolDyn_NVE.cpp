/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

//g++ -o prog1 MolDyn_NVE.cpp
int main(int argc, char** argv)
{
  if(argc != 5)
    {
      cerr << " for " << argv[0] << " first argument: <0> first run, <1> restart, second argument: <n steps> equilibration steps,  <1> compute averages, third argument: s,l of g for state of matter, fourth argument: <frames> to print frames, <0> no frames" << endl;
      return 1;
    }

  nstep_in=atof(argv[2]);
  if(argv[3][0]=='s') folder="solid";
  if(argv[3][0]=='l') folder="liquid";
  if(argv[3][0]=='g') folder="gas";
  cout<<"simulation of " << folder << " phase" << endl;
  string str(argv[4]);

  int nconf = 1;
  
  if(atof(argv[1])==0)
    Input();             //Inizialization
  else if(atof(argv[1])==1)
    Restart_Input();
  else cout << "input error" << endl;
 
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm

     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;

     if(istep%(block_size*10)==0 && atof(argv[2])==1)
       {
	 Print_Averages();
	 cout << "block index: " << block_index << endl;
	 block_index++;
       }
	 
     if(istep%10 == 0)
       {
	 Measure();     //Properties measurement
	 
	 if(str=="frames"){ ConfXYZ(nconf);    nconf += 1;};//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
	 
	 if(atof(argv[2])==1){ Update_Averages_sums();}
       }

     if(istep==nstep-1){string conf=folder+"/old.final"; PrintConf(conf.c_str());};
  }

  string conf=folder+"/old.0";
  PrintConf(conf.c_str());         //Write final configuration to restart
  
  return 0;
}












void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open(folder+"/input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  if(nstep_in!=1) nstep=nstep_in;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open(folder+"/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }


   //measurement of g(r)
   bin_size = (box/2.0)/(double)nbins;
   //r_array
   ofstream r_arr;
   r_arr.open(folder+"/r_array.dat");
   for(int i=0;i<nbins;i++) r_arr<<(i+0.5)*bin_size<<endl;
   r_arr.close();
   
   return;
}







void Restart_Input(void){ //Prepare all stuff for the simulation after restart
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl << endl;
  cout << "Simulation after restart" << endl;

  ReadInput.open(folder+"/input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> target_temp;

  if(nstep_in!=1) nstep=nstep_in;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read r(t-dt) configuration from previous simulation
  cout << "Read r(t-dt) configuration from file old.final " << endl << endl;
  ReadConf.open(folder+"/old.final");
  if(!ReadConf.good()){ cerr << "error in opening old.final" << endl; exit(-1); } 
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();

//Read r(t) configuration from previous simulation
  cout << "Read r(t) configuration from file old.0 " << endl << endl;
  ReadConf.open(folder+"/old.0");
  if(!ReadConf.good()){ cerr << "error in opening old.0" << endl; exit(-1); } 
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//compute r(t+dt) with one step of the Verlet algorythm,         !!! r(t)->x,y,z_old, r(t+dt)->x,y,z !!!
  Move(); 
  
  
//Prepare initial velocities
   cout << "Prepare velocities v(t+dt/2) and compute T(t+dt/2)" << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = (x[i] - xold[i])/delta;
     vy[i] = (y[i] - yold[i])/delta;
     vz[i] = (z[i] - zold[i])/delta;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
 
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * target_temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }



  //measurement of g(r)
  bin_size = (box/2.0)/(double)nbins;
   
   return;
}

  


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open(folder+"/output_epot.dat",ios::app);
  Ekin.open(folder+"/output_ekin.dat",ios::app);
  Temp.open(folder+"/output_temp.dat",ios::app);
  Etot.open(folder+"/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  for(int i=0; i<nbins; i++) g[i]=0.;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     if(dr<box/2.0)
       {
	 bin=int(dr/bin_size);                   //index to update
	 g[bin]+=2.;         
       }
     
     
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

   for(int i=0; i<nbins; i++) stima_gofr[i]=g[i]*Compute_gofr_coeff(bin_size*i);

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void PrintConf(const char *conf_file){ //Write final configuration
  ofstream WriteConf;

  cout << "Print configuration to file " << conf_file << " " << endl << endl;
  WriteConf.open(conf_file);

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open(folder+"/frames/config_" + to_string(nconf) + ".xyz");
  
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


void Update_Averages_sums(void)
{
  vector<double> obs = {stima_pot, stima_kin, stima_temp, stima_etot};
    
  for(int i=0; i<obs.size(); i++)
    {
      obs_sums[i]+=obs[i];
    }

  for(int i=0; i<nbins; i++)
    {
      gofr_sums[i]+=stima_gofr[i];
    }
}



void Print_Averages(void)
{
 for(int i=0; i<obs_sums.size(); i++)
   { 
    ofstream file;
    file.open(folder+"/ave_"+file_title[i]+".dat",ios::app);
    
    double block_av = obs_sums[i]/double(block_size);
    double block_av2 = block_av*block_av;
	
    prog_sum[i] += block_av;
    prog_sum2[i] += block_av2;

    double prog_av = prog_sum[i] / (block_index+1);
    double prog_av2 = prog_sum2[i] / (block_index+1);
    double prog_error = error(prog_av,  prog_av2, block_index);
    file << prog_av << " " << prog_error << endl;

    file.close();
   }



 //g(r)
 int wd=12;
 ofstream g_file;
 g_file.open(folder+"/gave.dat",ios::app);
 g_file << setw(wd) << block_index;   

 for(int i=0; i<nbins; i++)
   { 
    double block_av = gofr_sums[i]/double(block_size);
    double block_av2 = block_av*block_av;
	
    gofr_prog_sum[i] += block_av;
    gofr_prog_sum2[i] += block_av2;

    double prog_av = gofr_prog_sum[i] / (block_index+1);
    double prog_av2 = gofr_prog_sum2[i] / (block_index+1);
    double prog_error = error(prog_av,  prog_av2, block_index);
    g_file << setw(wd) << prog_av << setw(wd) << prog_error;
   }

 g_file<<endl;
 g_file.close();


 

 //reset sums
 vector<double> zeros(4, 0.);
 obs_sums=zeros;
 vector<double> zeros_g(nbins, 0.);
 gofr_sums=zeros_g;
}


double error(double av, double av2, int n)
{
  if(n==0)
    return 0;
  
  return sqrt( (av2 - av*av)/ double(n) );
}


double Compute_gofr_coeff(double r)
{
  double deltaV = (pow(r+bin_size,3)-pow(r,3))*4.*acos(-1)/3.;
  return 1./(deltaV*rho*npart);
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
