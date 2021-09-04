#include "mpi.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {
  int numtasks, rank, next, prev, buf[2], tag1=1,tag2=2;
  MPI_Request reqs[2]; // required variable for non-blocking calls
  MPI_Status stats[2]; // required variable for Waitall routine
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // determine left and right neighbors 
  prev = rank-1;
  next = rank+1;
  if (rank == 0) prev = numtasks - 1;
  if (rank == (numtasks - 1)) next = 0;

  //int message[2][2]={{rank,rank*10},{rank*20,rank*30}};
  vector<int> sequence = {2,0,1,3};
  vector<int> sender(4); 
  vector<int> recip(4); 
    
    
  for(int i=0; i<4; i++)
    {
      vector<int>::iterator itr = find(sequence.begin(), sequence.end(), i);
      int idx= distance(sequence.begin(), itr);
      sender[i]=sequence[(idx+3)%4];
      recip[i]=sequence[(idx+1)%4];
    }
  
  
  
  // post non-blocking receives and sends for neighbors
  MPI_Irecv(&buf[0], 1, MPI_INT, sender[rank], tag1, MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(&rank, 1, MPI_INT, recip[rank], tag1, MPI_COMM_WORLD, &reqs[1]);
  
  MPI_Waitall(2, reqs, stats);
  
  
  cout << "rank = " <<rank<<", buff = " << buf[0] << endl;
  
  MPI_Finalize();

}


/*

int main(int argc, char *argv[]) {
  int numtasks, rank, next, prev, tag1=1;
  MPI_Request reqs[2]; // required variable for non-blocking calls
  MPI_Status stats[2]; // required variable for Waitall routine
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // determine left and right neighbors
  int message[2][2]={{rank,rank*10},{rank*20,rank*30}};
  int gift[2][2];

  
  prev = rank-1;
  next = rank+1;
  if (rank == 0) prev = numtasks - 1;
  if (rank == (numtasks - 1)) next = 0;
  // post non-blocking receives and sends for neighbors
  MPI_Irecv(&gift[0], 4, MPI_INT, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
  //MPI_Irecv(&gift[1], 4, MPI_INT, next, tag2, MPI_COMM_WORLD, &reqs[1]);
  //MPI_Isend(&rank, 1, MPI_INT, prev, tag2, MPI_COMM_WORLD, &reqs[2]);
  MPI_Isend(&message, 1, MPI_INT, next, tag1, MPI_COMM_WORLD, &reqs[1]);
  // do some work while sends/receives progress in background
  // wait for all non-blocking operations to complete
  MPI_Waitall(4, reqs, stats);
  // continue - do more work

  
  cout << "rank = " <<rank<<", buff = " << gift[0][0] << " " << gift[1][1] << endl;
  
  MPI_Finalize();

}

 */
