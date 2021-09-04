#include <iostream>     // std::cout
#include <fstream>
#include "classes.h"
#include <string>
#include "mpi.h"
#include <stdlib.h>
#include <algorithm>
#include "random.h"

using namespace std;

//input istruction: mpiexec main.exe -np 4 a.out

int main(int argc, char** argv)
{

  //___________MPI_INITIALIZATION
  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  
  
  //____________PREPARING_MAP_AND_POPULATION__________________________________________  
  int n_ind=2000;
  int n_cities=32;
  int n_bests=5;
  int ev_steps=200;
  int delta_steps=20;
  srand(rank);
  
  Map myMap(n_cities, 1, rank);  //second argument: map_OPtion=1=random_in_square
  Population myPop(n_ind,n_cities,myMap,rank);
  myPop.Check_Bonds();
  myPop.Fitness_Sorting();
  cout << "rank "<<rank<< ": map and population ready for evolution" << endl;

  for(int o=0;o<int(ev_steps/delta_steps);o++)
    {
      myPop.Evolution(delta_steps);
      // cout << "prima dello scambio, rank = " <<rank<< ": ";

      vector<vector<int>> best_ind=myPop.Best_n_Chromo(n_bests);
      int a[n_bests][n_cities], buf[n_bests][n_cities];
      
      for(int k=0; k<n_bests; k++)
	{
	  for (int i=0; i<n_cities; i++){
	    //if(k==0) cout << "\t" << best_ind[k][i];
	    a[k][i]=best_ind[k][i];
	  }
	}
  
  
      MPI_Status stats[2];
      MPI_Request reqs[2];
      int tag1=1;
      
      vector<int> sequence(4);
      iota(sequence.begin(), sequence.end(), 0);
      rnd_shuffle(sequence.begin(), sequence.end(),1);
      if (rank==0)
	{
	  for (int i=0; i<4; i++)
	    cout << "\t" << sequence[i];
	}
  
      vector<int> sender(4); 
      vector<int> recip(4); 
      
      for(int i=0; i<4; i++)
	{
	  vector<int>::iterator itr = find(sequence.begin(), sequence.end(), i);
	  int idx= distance(sequence.begin(), itr);
	  sender[i]=sequence[(idx+3)%4];
	  recip[i]=sequence[(idx+1)%4];
	}
      

      MPI_Irecv(&buf, n_cities*n_bests, MPI_INT, sender[rank], tag1, MPI_COMM_WORLD, &reqs[0]);
      MPI_Isend(&a, n_cities*n_bests, MPI_INT, recip[rank], tag1, MPI_COMM_WORLD, &reqs[1]);
      
      MPI_Waitall(2, reqs, stats);
      

      vector<vector<int>> new_bests;
      for(int i=0; i<n_bests; i++)
	{
	  vector<int> i_best(n_cities);
	  for(int j=0; j<n_cities; j++)
	    {
	      i_best[j]=buf[i][j];
	    }
	  new_bests.push_back(i_best);
	}
  
	  
      myPop.Upload_Bests(new_bests,n_bests);
      cout << "Rank= "<<rank<< ": "<<delta_steps<< " steps evolution executed and "<<n_bests<<" best path exchanged"<<endl; 
    }


  
  MPI_Finalize();
  
  return 0;
}


