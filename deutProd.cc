//////////////////////////////////////////////////////////////////
// Strategy:                                                    //
// 1) make a function to read in data from file                 //
// 2) remembering which event it comes from                     //
// 3) make deutrons plainly                                     //
// 4) make a fit with the available data and distribute weights //
// 4) make deutrons with assigned weights                       //
//////////////////////////////////////////////////////////////////

// Headers and namespace
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "crossSections.h"

using namespace Pythia8;
using namespace std;

//int main(int argc, char* argv[] ) {
int main() {

  ///////////////////////////////////////////////////////////
  // Need to read in protons and neutrons produced by Abel //
  ///////////////////////////////////////////////////////////

  // Permanent declarations
  //vector<Particle_def> protons;

  // Create weights for protons and Neutrons
  // Particle_def permProt, permNeut;
  // cout << "Making Weights" << endl;
  // for (int i = 0; i < 64; i++)
  //   {
  //     // i needs to be used in filename, so needs to be converted to string
  //     stringstream num; num << i;
      
  //     // Fetch data from files, by first creating the location and filename, then fetching
  //     string dataProt = "../../results17Juni/results/data" + num.str() + "_prot_Vec4.dat";
  //     string dataNeut = "../../results17Juni/results/data" + num.str() + "_neut_Vec4.dat";

  //     permProt.set_read((char*) dataProt.c_str(), true);
  //     cout << "Number of protons loaded: " << permProt.Extract_Id().size() << endl;

  //     permNeut.set_read((char*) dataNeut.c_str(), true);
  //     cout << "Number of neutrons loaded: " << permNeut.Extract_Id().size() << endl;
      
  //   }
  // vector < vector <double> > protWeightHist, neutWeightHist;
  // protWeightHist = permProt.hist_weighting((char*) "../../plotting/ALICE_pbar_data_7.dat", 64*1000000);
  // neutWeightHist = permNeut.hist_weighting((char*) "../../plotting/ALICE_pbar_data_7.dat", 64*1000000);
  // cout << "Weights created" << endl;

  
  // A loop running over instances previously created by MPI
  unsigned int MPIs = 64;
  for (unsigned int i = 0; i < MPIs; i++)
    {
      // i needs to be used in filename, so needs to be converted to string
      stringstream num; num << i;

      // temporary declarations which needs to be remembered in the entire outside loop
      Particle_def tempProt, tempNeut;
      
      // Fetch data from files, by first creating the location and filename, then fetching
      //string dataProt = "../../results/results/data" + num.str() + "_protb_Vec4.dat";
      //string dataNeut = "../../results/results/data" + num.str() + "_neutb_Vec4.dat";
      string dataProt = "/media/thomhaa/StorageDrive/Documents/results/results/data" + num.str() + "_prot_Vec4.dat";
      string dataNeut = "/media/thomhaa/StorageDrive/Documents/results/results/data" + num.str() + "_neut_Vec4.dat";

      cout << "Reading file: " << dataProt << "\r" << endl;
      tempProt.set_read((char*) dataProt.c_str());
      tempProt.assign_weights_from_file((char*) "../../plotting/protonWeights.dat");
      
      cout << "Reading file: " << dataNeut << "\r" << endl;
      tempNeut.set_read((char*) dataNeut.c_str());
      tempNeut.assign_weights_from_file((char*) "../../plotting/neutronWeights.dat");
      
      // Extract the event number as a separate vector
      vector<int> protId = tempProt.Extract_Id();
      vector<int> neutId = tempNeut.Extract_Id();
      
      vector<Vec4> protVec4 = tempProt.vec4();
      vector<Vec4> neutVec4 = tempNeut.vec4();

      vector<long double> protWeights = tempProt.see_weight();
      vector<long double> neutWeights = tempNeut.see_weight();
      
      // Use an eventSorter to make sure the deutronProduction methods only use one event at a time
      
      string locationCo, locationCross, locationProt, locationNeut, path;
      path = "/media/thomhaa/StorageDrive/Documents/resultsDeuts7/";
      locationCo = path + "deutronCoalesce_100_" + num.str() + ".dat";
      locationCross = path + "deutronCrossSection_100_" + num.str() + ".dat";
      locationProt = path + "protWeighted" + num.str() + ".dat";
      locationNeut = path + "neutWeighted" + num.str() + ".dat";

      
      // A loop running over all events and sorting out the same event to separate lists, performed once per event
      cout << "Starting eventSorter: " << endl;

      // A loop running over all events and creating deuterons on a per event basis
      int start;
      start = 0;
      
      int protIter = protId[start];
      int neutIter = neutId[start];
      vector<Vec4> protEvent, neutEvent;
      vector< vector<Vec4> >  protEvents, neutEvents;
      vector<int> protIdList, neutIdList;
      vector<long double> protEventWeight, neutEventWeight;
      vector< vector<long double> > protEventsWeight, neutEventsWeight;
      // Run through all the events and sort on a per event basis
      for (int n = 0; n < abs(protId.size()+1); n++)
	{
	  if (protId[n] == protIter) { protEvent.push_back(protVec4[n]); protEventWeight.push_back(protWeights[n]); }
	  if (n == abs(protId.size()))
	    {
	      protEvents.push_back(protEvent);
	      protEventsWeight.push_back(protEventWeight);
	      protIdList.push_back(protIter);
	    }
	  else if (protId[n] > protIter)
	    {
	      protEvents.push_back(protEvent);
	      protEvent.clear();
	      protEvent.push_back(protVec4[n]);
	      
	      protEventsWeight.push_back(protEventWeight);
	      protEventWeight.clear();
	      protEventWeight.push_back(protWeights[n]);
	      
	      protIdList.push_back(protIter);
	      protIter = protId[n];
	    }
	  else if (protId[n] < protIter) { continue; }
	}
      for (int n = 0; n < abs(neutId.size()+1); n++)
	{
	  if (neutId[n] == neutIter) { neutEvent.push_back(neutVec4[n]); neutEventWeight.push_back(neutWeights[n]); }
	  if (n == abs(neutId.size()))
	    {
	      neutEvents.push_back(neutEvent);
	      neutEventsWeight.push_back(neutEventWeight);
	      neutIdList.push_back(neutIter);
	    }
	  else if (neutId[n] > neutIter)
	    {
	      neutEvents.push_back(neutEvent);
	      neutEvent.clear();
	      neutEvent.push_back(neutVec4[n]);
	      neutEventsWeight.push_back(neutEventWeight);
	      neutEventWeight.clear();
	      neutEventWeight.push_back(neutWeights[n]);
	      neutIdList.push_back(neutIter);
	      neutIter = neutId[n];
	    }
	  else if (neutId[n] < neutIter) { continue; }
	}
      
      vector<Vec4> prots, neuts;

      cout << "Generating deuterons" << endl;
      ////////////////////////////////////////////////////////////////////////////////////////////
      //                                                                                        //
      // protEvents and neutEvents now contain particles sorted after eventorder good particles //
      //                                                                                        //
      ////////////////////////////////////////////////////////////////////////////////////////////
      

      int protIntId, neutIntId;
      int m = 0;
      vector<long double> pT_bins;
      double bin_array [] = { 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.7, 2.9};
      pT_bins.insert( pT_bins.begin(), bin_array, bin_array+21);
      int append = 0;
      for (int n = 0; n < abs(protIdList.size()); n++)
       	{
	  protIntId = protIdList[n];
	  neutIntId = neutIdList[m];
	  while (protIntId >= neutIntId)
 	    {
 	      if ((protIntId == neutIdList[m]) & (protIntId != neutIdList[m-1]))
 		{
		  Particle_def protsP, neutsP;
		  Particle_def deutCo, deutXsec;
 		  // Perform event calculation here
 		  protsP.set(protEvents[n], protEventsWeight[n], protIntId);
		  neutsP.set(neutEvents[m], neutEventsWeight[m], neutIdList[m]);
 		  deutCo.coalescence(protsP, neutsP); deutXsec.x_sec(protsP, neutsP, 100);
		  //deutXsec.assign_weights(pT_bins, 100);
 		  deutCo.print_pT((char*) locationCo.c_str(), append, abs(protEvents.size()));
		  deutXsec.print_pT((char*) locationCross.c_str(), append, abs(protEvents.size()));
		  protsP.print_pT((char*) locationProt.c_str(), append, abs(protEvents.size()));
		  neutsP.print_pT((char*) locationNeut.c_str(), append, abs(neutEvents.size()));
 		  prots.clear(); neuts.clear();
 		}
 	      m++;
	      if (m < abs(neutIdList.size())) { neutIntId = neutIdList[m]; }
	      else { break; }
 	    }
	  append = 1;
 	  if ( n%10000 == 0 ) { cout << "Progress: " << round( n/abs(protIdList.size())*100 ) << "%         \r" << endl; }
 	}
      cout << locationCo << " have been produced ---- " << locationCross << " have been produced \r" << endl;
    }
  cout << endl;
  cout << "Program completed successfully" << endl;
  return 0;
}
