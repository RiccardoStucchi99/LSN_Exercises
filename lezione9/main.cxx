#include <iostream>     // std::cout
#include <fstream>
#include "classes.h"
#include <string>

using namespace std;


int main(int argc, char** argv)
{
  if(argc != 3) { cout << "first argument: <0> circle map, <1> random 2D map, second argument: evolution steps " << endl; return -1; }

  int opt=atoi(argv[1]);
  
  
  int n_ind=1000;
  int n_cities=32;

  Map myMap(n_cities, opt);
  Population myPop(n_ind,n_cities,myMap);

  //  myPop.Compute_Cost_Function();
  myPop.Check_Bonds();
  myPop.Fitness_Sorting();
  string folder=myMap.Get_Folder();
  //string outfile=folder+"test.dat";
  //myPop.print_pop(outfile.c_str());

  myPop.Evolution(atoi(argv[2]));
  string outfile_fin=folder+"final.dat";
  myPop.print_pop(outfile_fin.c_str());

  
  
  return 0;
}




//selection check
/*

  for(int i=0; i<2000;i++)
    {
      vector<Chromosome> vett=myPop.Selection();
    }
*/


//test crossover
 /*
  for(int h=0; h<10; h++)
    {
      vector<Chromosome> vett0=myPop.Selection();
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<5; j++)
	    {cout << vett0[i].Get_Genes()[j] << ' ';}
	  cout << endl;
	}
      


      vector<Chromosome> vett1=myPop.Crossover(vett0);
      for(int i=0; i<2; i++)
	{
	  for(int j=0; j<5; j++)
	    {cout << vett1[i].Get_Genes()[j] << ' ';}
	  cout << endl;
	}
      cout << endl;
      
    }
  */

//mutation test
  /*
  for(int h=0; h<10; h++)
    {
      vector<Chromosome> vett0=myPop.Selection();
      for(int j=0; j<n_cities; j++)
	{cout << vett0[0].Get_Genes()[j] << ' ';}
      cout << endl;
    
      Chromosome mut_cromo=myPop.Mutation(vett0[0]);
     
      for(int j=0; j<n_cities; j++)
	{cout << mut_cromo.Get_Genes()[j] << ' ';}
      cout << endl;
      
      cout << endl;
    }
  */
