#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "crossSections.h"

using namespace Pythia8;
using namespace std;

int main()
{

  //Permanent declarations
  vector<Particle_def> protons;

  //Create weights for protons and Neutrons
  Particle_def permProt, permNeut;
  cout << "Making Weights" << endl;
  for (int i = 0; i < 64; i++)
    {
      // i needs to be used in filename, so needs to be converted to string
      stringstream num; num << i;
      
      // Fetch data from files, by first creating the location and filename, then fetching
      string dataProt = "../../results17Juni/results/data" + num.str() + "_prot_Vec4.dat";
      string dataNeut = "../../results17Juni/results/data" + num.str() + "_neut_Vec4.dat";

      permProt.set_read((char*) dataProt.c_str(), true);
      cout << "Number of protons loaded: " << permProt.Extract_Id().size() << "- file number: " << i << endl;

      permNeut.set_read((char*) dataNeut.c_str(), true);
      cout << "Number of neutrons loaded: " << permNeut.Extract_Id().size() << "- file number: " << i << endl;
      
    }
  vector < vector <double> > protWeightHist, neutWeightHist;
  protWeightHist = permProt.hist_weighting((char*) "../../plotting/ALICE_pbar_data_7.dat", 64*1000000); 
  neutWeightHist = permNeut.hist_weighting((char*) "../../plotting/ALICE_pbar_data_7.dat", 64*1000000); 
  cout << "Weights created" << endl;

  permProt.print_pT((char*) "../../plotting/protPT.dat" );
  permNeut.print_pT((char*) "../../plotting/neutPT.dat" );
  
  ofstream outputFileProt;
  
  outputFileProt.open((char*) "../../plotting/weightsForProt.dat");
  outputFileProt << setprecision(32) << fixed;
  for (int i = 0; i < abs((protWeightHist[0]).size()); i++)
    {
      outputFileProt << (protWeightHist[1])[i+1] << " " << (protWeightHist[0])[i] << endl;
    }
  outputFileProt.close();
  cout << "written to file ../../plotting/weightsForProt.dat" << endl;
  
  ofstream outputFileNeut;
  outputFileNeut.open((char*) "../../plotting/weightsForNeut.dat");
  outputFileNeut << setprecision(32) << fixed;
  for (int i = 0; i < abs((neutWeightHist[0]).size()); i++)
    {
      outputFileNeut << (neutWeightHist[1])[i+1] << " " << (neutWeightHist[0])[i] << endl;
    }
  outputFileNeut.close();
  cout << "written to file ../../plotting/weightsForNeut.dat" << endl;

  return 0;
}
