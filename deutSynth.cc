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

  // Particle masses
  //const double mn = 0.93956536;   // Neutron
  //const double mp = 0.93827203;   // Proton
  //const double md = 1.875612793;  // Deuteron
  
  // Loop over all events read in by file
  int event, flick; flick = 999;
  vector<Vec4> p, pb, n, nb;
  Particle_def plist, pblist, nlist, nblist;
  string path, file, ppTOut, pbpTOut, npTOut, nbpTOut, deCopTOut, deXspTOut, debCopTOut, debXspTOut;
  string oP, oPb, oN, oNb, odeCo, odeXs, odebCo, odebXs, ifile, dat;
  path = "/media/thomhaa/StorageDrive/Documents/results/";
  ifile = "results/data";
  oP = "protpTOutwW", oPb = "protbpTOutwW", oN = "neutpTOutwW", oNb = "neutbpTOutwW";
  odeCo = "deutCopTOut", odeXs = "deutXspTOut", odebCo = "deutbarCopTOut", odebXs = "deutbarXspTOut";
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
      deCopTOut = path + odeCo + mpinum.str() + dat;
      deXspTOut = path + odeXs + mpinum.str() + dat;
      debCopTOut = path + odebCo + mpinum.str() + dat;
      debXspTOut = path + odebXs + mpinum.str() + dat;
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
		      istringstream iss(line);
		      iss >> event;
		      if (event%10000 == 0 ) { cout << "Event: " << mpi*1000000 + event << endl; }
		      
		      // Do event related things here with filled vetors
		      
		      //plist.set(p);
		      pblist.set(pb);
		      //nlist.set(n);
		      nblist.set(nb);
		      //plist.assign_weights_from_file((char*) "../../plotting/protonWeights.dat");
		      //nlist.assign_weights_from_file((char*) "../../plotting/protonWeights.dat");
		      pblist.assign_weights_from_file((char*) "../../plotting/protonbarWeights.dat");
		      nblist.assign_weights_from_file((char*) "../../plotting/protonbarWeights.dat");
		      if (nblist.see_weight().size() > 0) {cout << "Weight: " << nblist.see_weight()[0] << endl;}
		      
		      // Make deuterons
		      Particle_def deCo, debCo, deXs, debXs;
		      //deCo.coalescence(plist, nlist);
		      //deXs.x_sec(plist, nlist, 100);
		      debCo.coalescence(pblist, nblist);
		      debXs.x_sec(pblist, nblist, 100);
		      
		      // Write pT to file:
		      //deCo.print_pT( (char*) deCopTOut.c_str(), event);
		      //deXs.print_pT( (char*) deXspTOut.c_str(), event);
		      debCo.print_pT( (char*) debCopTOut.c_str(), event);
		      debXs.print_pT( (char*) debXspTOut.c_str(), event);
		      //plist.print_pT( (char*) ppTOut.c_str(), event);
		      //pblist.print_pT( (char*) pbpTOut.c_str(), event);
		      //nlist.print_pT( (char*) npTOut.c_str(), event);
		      //nblist.print_pT( (char*) nbpTOut.c_str(), event);
		      
		      
		      
		      // Clean the event
		      p.clear(), pb.clear(), n.clear(), nb.clear();
		    }
		  // Isolate the event
		  else
		    {
		      istringstream iss(line);
		      Vec4 temp;
		      long double px, py, pz, E;
		      iss >> px >> py >> pz >> E;
		      temp.p(px,py,pz,E);

		      // cases: 0 = p, 1 = pd, 2 = n, 3 = nd
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
