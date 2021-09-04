#include "classes.h"
#include <iostream>
#include <fstream>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <numeric>
#include <cmath>
#include <iterator>
#include <sstream>
#include <string>

using namespace std;


void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, Random &rnd)
{
  int n = (last-first);

  for (int i=n-1; i>0; --i) {
    swap(first[i],first[int(rnd.Rannyu()*(i+1))]);
  }
} 


Map :: Map(int n_cit, int opt)
{
  n_cities=n_cit;

  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  m_rnd.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  //circle option
  if(opt==0)
    {
      m_folder="circle/";
      ofstream file(m_folder+"cities.dat");
  
      for(int i=0; i<n_cities; i++)
	{
	  double teta= 2*acos(-1)*m_rnd.Rannyu();
	  double x=cos(teta); double y=sin(teta);
	  m_x_cities.push_back(x);    m_y_cities.push_back(y);
	  file<<x<<' '<<y<<endl;
	}

      file.close();
      cout<<"cities generated on an unitary circle"<<endl;
    }


  //random in square option
  else if(opt==1)
    {
      m_folder="Rnd_Square/";
      ofstream file(m_folder+"cities.dat");
  
      for(int i=0; i<n_cities; i++)
	{
	  double x=m_rnd.Rannyu(-1,1); double y=m_rnd.Rannyu(-1,1);
	  m_x_cities.push_back(x);    m_y_cities.push_back(y);
	  file<<x<<' '<<y<<endl;
	}

      file.close();
      cout<<"cities generated in a square"<<endl;
    }

  else{ cout << "map option error" << endl; exit(-1);};
      
}


void Map :: Set_Folder(int opt)
{ 
  if(opt==0){m_folder="circle/";}
  else if(opt==1){m_folder="Rnd_Square/";}
  else{ cout << "map option error" << endl; exit(-1);}
}


Chromosome :: Chromosome(Map map)
{
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  c_rnd.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  c_map=map;
  c_n_genes=map.Get_n_cities();
  
  vector<int> vect(c_n_genes);
  iota(vect.begin(),vect.end(),0);
  rnd_shuffle(vect.begin()+1, vect.end(), c_rnd);
  
  c_genes=vect;
  c_cost_func=Compute_Cost_Function();
  cout << "chromosome created" << endl;
}




double Chromosome :: Compute_Cost_Function(void)
{
  double sum=0.;
  int dim= c_n_genes;
  
  for(int i=0;i<dim;i++)
    {
      int n_city1=c_genes[i];
      int n_city2=c_genes[(i+1)%dim];
      
      double delta_x_2=pow( c_map.Get_x(n_city1) - c_map.Get_x(n_city2), 2.);
      double delta_y_2=pow( c_map.Get_y(n_city1) - c_map.Get_y(n_city2), 2.);      

      sum+=sqrt(delta_x_2+delta_y_2);
    }

  return sum;
}





double Chromosome :: Compute_Cost_Function(vector<int> new_genes)
{
  double sum=0.;
  if(c_n_genes!=int(new_genes.size())){ cout << "new genes size error" << endl; exit(-1);};
  int dim= c_n_genes;
  
  for(int i=0;i<dim;i++)
    {
      int n_city1=new_genes[i];
      int n_city2=new_genes[(i+1)%dim];
      
      double delta_x_2=pow( c_map.Get_x(n_city1) - c_map.Get_x(n_city2), 2.);
      double delta_y_2=pow( c_map.Get_y(n_city1) - c_map.Get_y(n_city2), 2.);      

      sum+=sqrt(delta_x_2+delta_y_2);
    }

  return sum;
}




void Chromosome :: Print_Chromosome(void)
{
  for (std::vector<int>::iterator it=c_genes.begin(); it!=c_genes.end(); ++it)
	cout << *it << ' ';
  cout << Get_Cost_Function() << '\n';
}



void Chromosome :: Check_Bonds(void)
{
  vector<int> to_check=Get_Genes();
  vector<int>::iterator it;

  for(int value=1; value<c_n_genes; value++)
    {
      it = find(to_check.begin(), to_check.end(), value);
      if (it != to_check.end())
	continue;
      else
	{
	  cerr << "bonds violation" << endl;
	  exit(-1);
	}
    }
  cout << "check completed with success" << endl;
}




vector<int> Chromosome :: Mutation()
{
  double r=c_rnd.Rannyu();
  vector<int> genes=c_genes;
  
  //mutation algorythm 1
  if(r<0.25)
    {
      int idx1=int(c_rnd.Rannyu(1,c_n_genes));
      int idx2=int(c_rnd.Rannyu(1,c_n_genes));
      // cout<<idx1<<' '<<idx2<<endl;
      // cout<<"mutation1"<<endl;
      
      swap(genes[idx1],genes[idx2]);
    }
  
  
  //mutation algorythm 2
  if(r>0.25 && r<0.5)
    {
       
      int n=int(c_rnd.Rannyu(1,c_n_genes));
      //cout << n << endl;
      //cout<<"mutation2"<<endl;
      for(int i=1;i<c_n_genes;i++)
	{
	  if(genes[i]+n<c_n_genes)
	    {
	      genes[i]+=n;
	    }
	  else
	    {
	      genes[i]=(genes[i]+n)%(c_n_genes) +1;
	    }
	}
    }
  
  if(r>0.5 && r<0.75)
    {
      int m=int(c_rnd.Rannyu(2,c_n_genes/2));
      int idx1=int(c_rnd.Rannyu(1,c_n_genes-2*m));
      int idx2=int(c_rnd.Rannyu(idx1+m,c_n_genes-m));
      //cout << m << ' ' << idx1 << ' ' << idx2 << endl;
      //cout<<"mutation3"<<endl;
      
      swap_ranges(genes.begin()+idx1,genes.begin()+idx1+m,genes.begin()+idx2);
    }

  
  if(r>0.75)
    {
      int m=int(c_rnd.Rannyu(2,c_n_genes));
      int idx1=int(c_rnd.Rannyu(1,c_n_genes-m));
      //cout<<"mutation4"<<endl;
      //cout << m << ' ' << idx1 << endl;
      
      reverse(genes.begin()+idx1, genes.begin()+m+idx1);
    }  

  return genes;
}


void Chromosome :: SA_step(double temp)
{
  vector<int> new_genes=Mutation();
 
  double new_cost=Compute_Cost_Function(new_genes);
  double old_cost=c_cost_func;
  double beta=1./temp;
  double prob_exp=exp(-beta*(new_cost-old_cost));
  double prob=min(1.,prob_exp);

  double r=c_rnd.Rannyu();
  if(r<prob){ c_genes=new_genes; c_cost_func=new_cost;};
  
  ofstream cost_file;
  cost_file.open(c_map.Get_Folder()+"cost_file.dat",ios::app);
  cost_file << c_cost_func << endl;
  cost_file.close();
}


void Chromosome :: SA_evolution(int n_steps, double temp)
{
  for(int i=0;i<n_steps;i++)
    {
      SA_step(temp);
    }
}

void Chromosome :: Print_Conf_Final(void)
{
  ofstream file(c_map.Get_Folder()+"cromo_final.dat");

  for(int i=0;i<c_map.Get_n_cities();i++)
    {
      file << c_map.Get_x(i) << ' ';
    }

  file<<endl;
  for(int i=0;i<c_map.Get_n_cities();i++)
    {
      file << c_map.Get_y(i) << ' ';
    }

  file<<endl;
  for(int i=0;i<c_n_genes;i++)
    {
      file << c_genes[i] << ' ';
    }

  cout << "final chromosome printed in " << c_map.Get_Folder() << "cromo_final.dat" << endl << endl;

}



void Chromosome :: Restart(int opt)
{
  c_map.Set_Folder(opt);
  ifstream file(c_map.Get_Folder()+"cromo_final.dat");
  if(!file.good()) {cerr<<"error in opening final config file"<<endl; exit(-1);};

  vector<double> x; vector<double> y; vector<int> genes;
 
  string tx;
  getline(file, tx); 
  istringstream inx(tx);
  copy(istream_iterator<double>(inx), 
       istream_iterator<double>(), 
       back_inserter(x));


  string ty;
  getline(file, ty); 
  istringstream iny(ty);
  copy(istream_iterator<double>(iny), 
       istream_iterator<double>(), 
       back_inserter(y));

  string tg;
  getline(file, tg); 
  istringstream ing(tg);
  copy(istream_iterator<double>(ing), 
       istream_iterator<double>(), 
       back_inserter(genes));

  
  
  c_map.Set_X_Y(x,y);
  
  c_n_genes=genes.size();
  c_genes=genes;
  c_cost_func=Compute_Cost_Function();

  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  c_rnd.MySetRandom(infile_primes.c_str(), infile_seed.c_str());

  cout << "restart completed with success" << endl;
}




