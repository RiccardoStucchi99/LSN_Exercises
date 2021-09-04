#include <iostream>     // std::cout
#include <fstream>
#include "classes.h"

using namespace std;


int main(int argc, char** argv)
{
  if(argc != 5)
    {
      cerr << " for " << argv[0] << " first argument: <0> first run, <1> restart, second: temperaure, third: nsteps, fourth: <0> circle, <1> square" << endl;
      return 1;
    }

  int map_opt=atoi(argv[4]);
  int opt=atoi(argv[1]);
  double temp=atof(argv[2]);
  int n_steps=atoi(argv[3]);
  if(opt==0)
    {
      int n_cities=32;

      Map myMap(n_cities, map_opt);
      Chromosome myCromo(myMap);
      cout<<"first run ready to go!"<<endl;

      cout<<endl;
      cout << "initial path and cost: ";
      myCromo.Print_Chromosome();
      myCromo.Check_Bonds();
      cout << endl;
      
      myCromo.SA_evolution(n_steps,temp);
      cout << "simulation completed with (n_steps,temp) = (" << n_steps << "," << temp << ")" << endl << endl;

      cout << "final path and cost: ";
      myCromo.Print_Chromosome();
      myCromo.Print_Conf_Final();
    }
  
  if(opt==1)
    {
      Chromosome myCromo;
      myCromo.Restart(map_opt);
   
      

       cout<<endl;
      cout << "initial path and cost: ";
      myCromo.Print_Chromosome();
      myCromo.Check_Bonds();
      cout << endl;

      myCromo.SA_evolution(n_steps,temp);
      cout << "simulation completed with (n_steps,temp) = (" << n_steps << "," << temp << ")" << endl  << endl;
      
      cout << "final path and cost: ";
      myCromo.Print_Chromosome();
      myCromo.Print_Conf_Final();
    }
  
  return 0;
}





