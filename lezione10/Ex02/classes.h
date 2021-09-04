#ifndef __classes_h__
#define __classes_h__

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <numeric>
#include "random.h"
#include <string>

using namespace std;


void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, Random &rnd);
void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, int seed);

class Map
{
public:
  Map(){};
  Map(int n_cit, int opt, int rank);
  ~Map(){};
  double Get_x(int city_idx){return m_x_cities[city_idx];};
  double Get_y(int city_idx){return m_y_cities[city_idx];};
  string Get_Folder(void){ return m_folder;};

private:
  int n_cities;
  vector<double> m_x_cities;
  vector<double> m_y_cities;
  Random m_rnd;
  string m_folder;
};


class Chromosome
{
public:
  Chromosome(){};
  Chromosome(vector<int> _genes, Map *_map){c_genes=_genes; c_map=_map; c_cost_func=Compute_Cost_Function();};
  ~Chromosome(){};
  double Compute_Cost_Function(void);
  double Get_Cost_Function(void) const {return c_cost_func;};
  vector<int> Get_Genes(void){return c_genes;};
  void Check_Bonds(void); 


private:
  vector<int> c_genes;
  double c_cost_func;
  Map *c_map;
};


class Population
{
public:
  Population(int n_ind, int n_genes, Map _map,int seed);  //da modificaere
  ~Population(){};
  Chromosome Get_Genome(int i){return p_pop[i];};
  void print_pop(const char *file);
  void Compute_Cost_Function(void);
  double Compute_Cost_Function(int i);
  void Fitness_Sorting(void);
  void Check_Bonds(void);
  void Print_Best_Chromo(void);
  void Print_Best_Chromo_Terminal(void);
  vector<vector<int>> Best_n_Chromo(int);
  
  //evolution
  vector<Chromosome> Selection(void);
  vector<Chromosome> Crossover(vector<Chromosome>);
  Chromosome Mutation(Chromosome);
  void Evolution_Step(void);
  void Evolution(int n_steps);
  
  void Print_Cost_Function(void);
  void Upload_Bests(vector<vector<int>>,int);
  

private:
  int p_n_ind;
  int p_n_genes;
  Random p_rnd;
  Map p_map;
  vector<Chromosome> p_pop;
  vector<double> p_cost_func;
  double p_exp;    //selection exponential
  double p_co_prob=0.6;  //crossover probability
};

#endif
