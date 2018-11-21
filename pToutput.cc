#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "crossSections.h"

using namespace Pythia8;
using namespace std;

int main() {

  // Loop over all events read in by file
  int event, flick;
  vector<Vec4> p, pb, n, nb;
  Particle_def plist, pblist, nlist, nblist;
  string path, file, ppTOut, pbpTOut, npTOut, nbpTOut, oP, oPb, oN, oNb, ifile, dat;
  path = "/media/thomhaa/StorageDrive/Documents/results/";
  ifile = "results/data";
  oP = "protpTOut", oPb = "protbpTOut", oN = "neutpTOut", oNb = "neutbpTOut";
  dat = ".dat";
    
  int MPI = 64;
  for (int mpi = 0; mpi < MPI; mpi++)
    {
      stringstream mpinum; mpinum << mpi;
      file = path + ifile + mpinum.str() + dat;
      ppTOut = path + oP + mpinum.str() + dat;
      pbpTOut = path + oPb + mpinum.str() + dat;
      npTOut = path + oN + mpinum.str() + dat;
      nbpTOut = path + oNb + mpinum.str() + dat;
      ifstream infile;
      infile >> setprecision(64) >> fixed;
      string line;
      infile.open(file.c_str());
      if (infile)
	{
	  while (getline(infile, line) )
	    {
	      if (!line.empty())
		{
		  if (line[1] == 'p')
		    {
		      // Fill p vector
		      if (line[2] == 'b') { flick = 1; }
		      else { flick = 0; }
		    }
		  else if (line[1] == 'n')
		    {
		      if (line[2] == 'b') { flick = 3; }
		      else { flick = 2; }
		    }
		  // Execute event operation and Start new event
		  else if (line[0] == '#')
		    { 
		      // Checks and balances
		      line.erase(0, 1);
		      if (event%1000 == 0 ) { cout << "Event: " << event << endl; }
		      istringstream iss(line);
		      iss >> event;
		      
		      // Do event related things here with filled vetors
		      
		      plist.set(p ), pblist.set(pb), nlist.set(n), nblist.set(nb);
		      
		      
		      // Write pT to file:
		      plist.print_pT( (char*) ppTOut.c_str(), event);
		      pblist.print_pT( (char*) pbpTOut.c_str(), event);
		      nlist.print_pT( (char*) npTOut.c_str(), event);
		      nblist.print_pT( (char*) nbpTOut.c_str(), event);
		      
		      
		      
		      // Fresh event from here on
		      p.clear(), pb.clear(), n.clear(), nb.clear();
		    }
		  else
		    {
		      istringstream iss(line);
		      Vec4 temp;
		      long double px, py, pz, E;
		      iss >> px >> py >> pz >> E;
		      temp.p(px,py,pz,E);
		      switch (flick)
			{
			case 0:
			  p.push_back(temp);
			  break;
			case 1:
			  pb.push_back(temp);
			  break;
			case 2:
			  n.push_back(temp);
			  break;
			case 3:
			  nb.push_back(temp);
			  break;
			}
		    }
		}
	    }
	}
      else cout << "Unable to open file." << endl;
    }

  
  return 0;
}
