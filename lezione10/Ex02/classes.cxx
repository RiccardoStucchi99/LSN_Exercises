#include "classes.h"
#include <iostream>
#include <fstream>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <numeric>
#include <cmath>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;


void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, Random &rnd)
{
  int n = (last-first);

  for (int i=n-1; i>0; --i) {
    swap(first[i],first[int(rnd.Rannyu()*(i+1))]);
  }
} 


void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, int seed)
{
  srand ( seed );
  int n = (last-first);

  for (int i=n-1; i>0; --i) {
    swap(first[i],first[rand()%(i+1)]);
  }
}


Map :: Map(int n_cit, int opt, int rank)
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
      m_folder="RiS_rank"+to_string(rank)+"/";
      ofstream file(m_folder+"cities.dat");
  
      for(int i=0; i<n_cities; i++)
	{
	  double x=m_rnd.Rannyu(-1,1); double y=m_rnd.Rannyu(-1,1);
	  m_x_cities.push_back(x);    m_y_cities.push_back(y);
	  file<<x<<' '<<y<<endl;
	}

      file.close();
      cout<<"Rank= "<<rank<< ": cities generated in a square"<<endl;
    }

  else{ cout << "map option error" << endl; exit(-1);};
      
}



double Chromosome :: Compute_Cost_Function(void)
{
  double sum=0.;
  int dim= c_genes.size();
  
  for(int i=0;i<dim;i++)
    {
      int n_city1=c_genes[i];
      int n_city2=c_genes[(i+1)%dim];
      
      //sum_x+=pow( c_map->Get_x(n_city1) - c_map->Get_x(n_city2), 2.);  //L2
      //sum_y+=pow( c_map->Get_y(n_city1) - c_map->Get_y(n_city2), 2.);      

      double delta_x_2=pow( c_map->Get_x(n_city1) - c_map->Get_x(n_city2), 2.);
      double delta_y_2=pow( c_map->Get_y(n_city1) - c_map->Get_y(n_city2), 2.);      

      sum+=sqrt(delta_x_2+delta_y_2);
    }

  return sum;
}



void Chromosome :: Check_Bonds(void)
{
  vector<int> genome_to_check=Get_Genes();
  vector<int>::iterator it;

  for(int value=1; value<int(genome_to_check.size()); value++)
    {
      it = find(genome_to_check.begin(), genome_to_check.end(), value);
      if (it != genome_to_check.end())
	continue;
      else
	{
	  cerr << "bonds violation" << endl;
	  exit(-1);
	}
    }
  //cout << "check completed with success" << endl;
}




Population :: Population(int n_ind, int n_genes, Map mappa, int random_seed)
{
  string infile_primes = "Primes";
  string infile_seed = "seed.in";
  p_rnd.MySetRandom(infile_primes.c_str(), infile_seed.c_str(), random_seed);

  p_n_ind=n_ind;
  p_n_genes=n_genes;
  p_map=mappa;
   
  for(int i=0;i<p_n_ind;i++)
    {
      vector<int> vect(p_n_genes);
      iota(vect.begin(),vect.end(),0);

      rnd_shuffle(vect.begin()+1, vect.end(), p_rnd);
      //random_shuffle(vect.begin()+1, vect.end());
      Chromosome cromo(vect,&p_map); 
      
      p_pop.push_back(cromo);
    }
}


void Population :: Check_Bonds(void)
{
  for(int i=0;i<p_n_ind;i++)
    {
      vector<int> genome_to_check=Get_Genome(i).Get_Genes();
      vector<int>::iterator it;

      for(int value=1; value<p_n_genes; value++)
	{
	  it = find(genome_to_check.begin(), genome_to_check.end(), value);
	  if (it != genome_to_check.end())
	    continue;
	  else
	    {
	      cerr << "bonds violation" << endl;
	      exit(-1);
	    }
	}
    }
  //cout << "check completed with success" << endl;
}


void Population :: print_pop(const char *outfile)
{
  ofstream file(outfile);

  for(int i=0;i<p_n_ind;i++)
    {
      Chromosome cromo = Get_Genome(i);
      vector<int> genome=cromo.Get_Genes();
      for (std::vector<int>::iterator it=genome.begin(); it!=genome.end(); ++it)
	file << *it << ' ';
      file << cromo.Get_Cost_Function() << '\n';
    }

  file.close();
  cout << "population printed in " << outfile << endl;
}



void Population :: Print_Best_Chromo()
{
  ofstream file(p_map.Get_Folder()+"/best_chromo.dat",ios::app);

  Chromosome cromo = Get_Genome(0);
  vector<int> genome=cromo.Get_Genes();
  for (std::vector<int>::iterator it=genome.begin(); it!=genome.end(); ++it)
    file << *it << ' ';
  file << '\n';
    

  file.close();
}


void Population :: Print_Best_Chromo_Terminal()
{
  Chromosome cromo = Get_Genome(0);
  vector<int> genome=cromo.Get_Genes();
  for (std::vector<int>::iterator it=genome.begin(); it!=genome.end(); ++it)
    cout << *it << ' ';
  cout << '\n';
}



void Population :: Fitness_Sorting(void)
{
  sort(p_pop.begin(), p_pop.end(), [] (const Chromosome &c1, const Chromosome &c2) -> bool
				   {
				     return c1.Get_Cost_Function() < c2.Get_Cost_Function();
				   });

  //cout << "population sorted" << endl;
}



//_____________________EVOLUTION OPERATORS

vector<Chromosome> Population :: Selection(void)
{
  double r=p_rnd.Rannyu();
  double s=p_rnd.Rannyu();
  p_exp=2;
  int j=int(p_n_ind*pow(r,p_exp)); 
  int k=int(p_n_ind*pow(s,p_exp));
  vector<Chromosome> selected{p_pop[j], p_pop[k]};

  //check if the exponential is right
  /*
  ofstream file;
  file.open("selection_index.dat",ios::app);
  file << j << endl;  file << k << endl;
  file.close();
  */
  return selected;
}


vector<Chromosome> Population :: Crossover(vector<Chromosome> selected)
{
  double r=p_rnd.Rannyu();

  vector<int> par1=selected[0].Get_Genes();
  vector<int> par2=selected[1].Get_Genes();
  vector<int> son1;
  vector<int> son2;

  
  if(r>p_co_prob)
    {
      //cout<<"nocross"<<endl;
      son1=par1;
      son2=par2;
    }
  else
    {
      //cout<<"cross"<<endl;
      int cut_idx=int(p_rnd.Rannyu(1,p_n_genes));
      // cout<<cut_idx<<endl;
      son1.insert(son1.begin(), par1.begin(), par1.begin()+cut_idx);
      son2.insert(son2.begin(), par2.begin(), par2.begin()+cut_idx);

      //add tail to son 1
      vector<int>::iterator it1;

      for(std::vector<int>::iterator it=par2.begin(); it!=par2.end(); ++it)
	{
	  int new_value=*it;
	  it1 = find(son1.begin(), son1.end(), new_value);
	  if (it1 == son1.end())
	    son1.push_back(new_value);
	}


      //add tail to son 2
      vector<int>::iterator it2;

      for(std::vector<int>::iterator it=par1.begin(); it!=par1.end(); ++it)
	{
	  int new_value=*it;
	  it1 = find(son2.begin(), son2.end(), new_value);
	  if (it1 == son2.end())
	    son2.push_back(new_value);
	}
      
    }

  Chromosome offspr1(son1,&p_map);
  Chromosome offspr2(son2,&p_map);
  vector<Chromosome> crossed_cromos{offspr1, offspr2};

  return crossed_cromos;
}

 
Chromosome Population :: Mutation(Chromosome cromo)
{
  double r=p_rnd.Rannyu();
  vector<int> genes=cromo.Get_Genes();
  vector<int> new_genes;
  
  //mutation algorythm 1
  if(r<0.1)
    {
     
      int idx1=int(p_rnd.Rannyu(1,genes.size()));
      int idx2=int(p_rnd.Rannyu(1,genes.size()));
      // cout<<idx1<<' '<<idx2<<endl;
      // cout<<"mutation1"<<endl;
      
      swap(genes[idx1],genes[idx2]);
    }
  
  
  //mutation algorythm 2
  if(r>0.1 && r<0.2)
    {
       
      int n=int(p_rnd.Rannyu(1, genes.size()));
      //cout << n << endl;
      //cout<<"mutation2"<<endl;
      for(int i=1;i<int(genes.size());i++)
	{
	  if(genes[i]+n<int(genes.size()))
	    {
	      genes[i]+=n;
	    }
	  else
	    {
	      genes[i]=(genes[i]+n)%(genes.size()) +1;
	    }
	}
    }
  
  if(r>0.2 && r<0.3)
    {
      int m=int(p_rnd.Rannyu(2,genes.size()/2));
      int idx1=int(p_rnd.Rannyu(1,genes.size()-2*m));
      int idx2=int(p_rnd.Rannyu(idx1+m,genes.size()-m));
      //cout << m << ' ' << idx1 << ' ' << idx2 << endl;
      //cout<<"mutation3"<<endl;
      
      swap_ranges(genes.begin()+idx1,genes.begin()+idx1+m,genes.begin()+idx2);
    }

  
  if(r>0.3 && r<0.4)
    {
      int m=int(p_rnd.Rannyu(2, genes.size()));
      int idx1=int(p_rnd.Rannyu(1,genes.size()-m));
      //cout<<"mutation4"<<endl;
      //cout << m << ' ' << idx1 << endl;
      
      reverse(genes.begin()+idx1, genes.begin()+m+idx1);
    }  
  
  //if(r>0.4) cout << "no mutation" << endl;
  new_genes=genes;
  Chromosome new_cromo(new_genes, &p_map);
  return new_cromo;
}




void Population :: Evolution_Step(void)
{
  vector<Chromosome> new_cromos;

  while(int(new_cromos.size())<p_n_ind)
    {
      vector<Chromosome> parents=Selection();
      vector<Chromosome> crossed_parents=Crossover(parents);
      Chromosome son1=Mutation(crossed_parents[0]);
      Chromosome son2=Mutation(crossed_parents[1]);

      new_cromos.push_back(son1); new_cromos.push_back(son2);
    }

  p_pop=new_cromos;
}



void Population :: Print_Cost_Function()
{
  string folder=p_map.Get_Folder();

  ofstream best_path;
  best_path.open(folder+"best_path.dat",ios::app);
  
  ofstream best_half;
  best_half.open(folder+"best_half.dat",ios::app);

  double sum=0.;
  for(int i=0; i<p_n_ind/2; i++)
    {
      sum+=p_pop[i].Get_Cost_Function();
    }

  best_path << p_pop[0].Get_Cost_Function() << endl;
  best_half << sum/double(p_n_ind/2) << endl;

  best_path.close();
  best_half.close();
}

  
void Population :: Evolution(int n_steps)
{
  for(int i=0;i<n_steps;i++)
    {
      Evolution_Step();
      Fitness_Sorting();
      Print_Cost_Function();
      //Print_Best_Chromo();
      Check_Bonds();
    }
}


vector<vector<int>> Population :: Best_n_Chromo(int n)
{
  vector<vector<int>> best_n;

  for(int i=0; i<n; i++)
    {
      Chromosome c_i=Get_Genome(i);
      vector<int> genes= c_i.Get_Genes();
      best_n.push_back(genes);
    }

  return best_n;
}


void Population :: Upload_Bests(vector<vector<int>> new_paths, int n_bests)
{
  for(int i=0; i<n_bests; i++)
    {
      Chromosome new_chromo(new_paths[i], &p_map);
      p_pop[i]=new_chromo;
    }
  Check_Bonds();
  Fitness_Sorting();
}

/*
  

void Population :: Print_Cost_Function()
{
  ofstream file("cost_function.dat");
  for(int i=0; i<m_ind; i++)
    {
      double cf=Compute_Cost_Function(i);
      cost_func[i]=cf;
      file << cf << endl;
    }
  file.close();
}







double Population :: Compute_Cost_Function(int ind)
{
  if(ind>=m_ind){cerr<<"index error in computing cost function"<<endl; exit(-1);};

  double sum_x=0.;
  double sum_y=0.;
  
  for(int i=0;i<m_crom-1;i++)
    {
      int n_city1=m_pop[ind][i];
      int n_city2=m_pop[ind][i+1];
      
      sum_x+=pow( m_map.Get_x(n_city1) - m_map.Get_x(n_city2), 2.);
      sum_y+=pow( m_map.Get_y(n_city1) - m_map.Get_y(n_city2), 2.);      
    }

  return sum_x+sum_y;
}
*/

