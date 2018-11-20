#ifndef CROSSSECTIONS_H //Making sure functions are declared only once
#define CROSSSECTIONS_H

#include <cmath>
#include "Pythia8/Pythia.h"
#include <limits>

// Declare functions
using std::sqrt;
using std::pow;
using std::exp;
using namespace Pythia8;

// Particle masses (GeV)
const double mn = 0.93956536;   // Neutron
const double mp = 0.93827203;   // Proton
const double md = 1.875612793;  // Deuteron
const double mpic = 0.13957018; // Charged pion
const double mpi0 = 0.1349766;  // Neutral pion

// constants
const double sigma_0_inv = 2.58 / 1000000; // units is microbarn &try 2.58 for db, 2.63 for d at 7TeV
const double pi = 3.14159265359;

// Fit function used for pion processes
double fitFunc(double x, double a, double b, double c, double d, double e)
{
    return a*pow(x,b)/(pow((c-exp(d*x)),2)+e);
}

// N N -> d pi helper function
double xs_pp_dpip_q(double q)
{
    double eta = q/mpic;
    double a[5] = {0.17, 1.34, 1.77, 0.38, 0.096};    
    return 1e3 * fitFunc(eta, a[0],a[1],a[2],a[3],a[4]);
}

/////////////////////////////////////////////////////////////////
//  Cross section parameterizations                            //
/////////////////////////////////////////////////////////////////

// p n -> d gamma
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pn_dgamma(double k)
{
    double a[12] = {2.3034605532591175,  -93.663463313902028, 2565.3904680353621, 
                    -25594.100560137995, 143513.10872427333,  -503572.89020794741, 
                    1149248.0196165806,  -1723683.9119787284, 1679348.7891145353, 
                    -1019888.5470232342, 349840.35161061864,  -51662.760038375141};
    double b[2]  = {-5.1885266705385051, 2.9195632726211609};
    if(k<1.28)
    {
        double result = 0;        
        for(int i=0;i<12;i++)
        {
            result += a[i] * pow(k,i-1);
        }
        return result;        
    }
    else
    {
        return exp(-b[0]*k -b[1]*k*k);
    }
}
// p n -> d pi0
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pn_dpi0(double k)
{
    double E_CoM = sqrt(mp*mp+0.25*k*k) + sqrt(mn*mn+0.25*k*k);
    double s = E_CoM*E_CoM;
    if(E_CoM < md+mpi0) 
        return 0;
    double q = sqrt(0.25*pow(s+mpi0*mpi0-md*md,2)/s - mpi0*mpi0);
    return 0.5*xs_pp_dpip_q(q); //**
}
// p n -> d pi+ pi-
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pn_dpippim(double k)
{
    double E_CoM = sqrt(mp*mp+0.25*k*k) + sqrt(mn*mn+0.25*k*k);
    if(E_CoM < md+2*mpic) 
        return 0;    
    double a[10] = {6.46455516e+06, 1.05136338e+01, 1.97924778e+03, 5.36301369e+00, 6.04534114e+05, 2.54935423e+15, 1.65669163e+01, 2.32961298e+07, 1.11937373e+01, 2.86815089e+16};
    return fitFunc(k, a[0],a[1],a[2],a[3],a[4]) + fitFunc(k, a[5],a[6],a[7],a[8],a[9]);    
}
// p n -> d pi0 pi0
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pn_dpi0pi0(double k)
{
    double E_CoM = sqrt(mp*mp+0.25*k*k) + sqrt(mn*mn+0.25*k*k);
    if(E_CoM < md+2*mpi0) 
        return 0;     
    double a[5] = {2.85519622e+06, 1.31114126e+01, 2.96145497e+03, 5.57220777e+00, 1.46051932e+06};
    return fitFunc(k, a[0],a[1],a[2],a[3],a[4]);
}
// p p -> d pi+
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pp_dpip(double k)
{
    double E_CoM = 2*sqrt(mp*mp+0.25*k*k);
    double s = E_CoM*E_CoM;
    if(E_CoM < md+mpic) 
        return 0;
    double q = sqrt(0.25*pow(s+mpic*mpic-md*md,2)/s - mpic*mpic);
    return xs_pp_dpip_q(q); //**
}
// p p -> d pi+ pi0
// Returns cross section in microbarn. k must be in units of GeV.
double xs_pp_dpippi0(double k)
{
    double E_CoM = 2*sqrt(mp*mp+0.25*k*k);
    if(E_CoM < md+mpic+mpi0) 
        return 0;     
    double a[5] = {5.09870846e+15, 1.65581228e+01, 2.33337076e+07, 1.13304315e+01, 2.86815089e+16};
    return fitFunc(k, a[0],a[1],a[2],a[3],a[4]);    
}
// n n -> d pi-
// Returns cross section in microbarn. k must be in units of GeV.
double xs_nn_dpim(double k)
{
    double E_CoM = 2*sqrt(mn*mn+0.25*k*k);
    double s = E_CoM*E_CoM;
    if(E_CoM < md+mpic) 
        return 0;  
    double q = sqrt(0.25*pow(s+mpic*mpic-md*md,2)/s - mpic*mpic);
    return xs_pp_dpip_q(q); //**
}
// n n -> d pi- pi0
// Returns cross section in microbarn. k must be in units of GeV.
double xs_nn_dpimpi0(double k)
{
    double E_CoM = 2*sqrt(mn*mn+0.25*k*k);
    if(E_CoM < md+mpic+mpi0) 
        return 0;      
    double a[5] = {5.09870846e+15, 1.65581228e+01, 2.33337076e+07, 1.13304315e+01, 2.86815089e+16};
    return fitFunc(k, a[0],a[1],a[2],a[3],a[4]);   
}

/////////////////////////////////////////////////////////////////////
/////	My Own 				   //////////////////////////
/////////////////////////////////////////////////////////////////////

//Variables

long double p_abs(Vec4 a, Vec4 b, double m1, double m2)
{
  long double s = (a + b).m2Calc();
  long double p_abs = sqrt((s - pow((m1 + m2),2))*(s - pow((m1 - m2),2))) / (2*sqrt(s));
  return p_abs;
}

// Need lab fram vector a and b from the original particles
// But also the original particle mass
vector< vector<long double> > Gamma(Vec4 a, Vec4 b, long double ma, long double mb)
{
	vector< vector<long double> > G;
	vector<long double> g;
	Vec4 p = a + b;
	long double pCM = p_abs(a, b, ma, mb);
	long double E = sqrt(pCM*pCM + ma*ma) + sqrt(pCM*pCM + mb*mb);
	long double gamma = (p.e())/E;
	
	long double beta = sqrt(1-1/(gamma*gamma));
	long double nx = p.px()/(E*gamma*beta);
	long double ny = p.py()/(E*gamma*beta);
	long double nz = p.pz()/(E*gamma*beta);

	long double g0[] = {gamma, gamma*beta*nx, gamma*beta*ny, gamma*beta*nz};
	long double g1[] = {g0[1], 1+(gamma-1)*nx*nx, (gamma-1)*nx*ny, (gamma-1)*nx*nz};
	long double g2[] = {g0[2], g1[2], 1+(gamma-1)*ny*ny, (gamma-1)*ny*nz};
	long double g3[] = {g0[3], g1[3], g2[3], 1+(gamma-1)*nz*nz};
	
	for (int i = 0; i<4; i++) {
		g.push_back (g0[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g1[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g2[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g3[i]);
	}
	G.push_back (g); g.clear();

	return G;
}


///////////////////////////////////////
// This is the original gamma matrix //
///////////////////////////////////////

/*
vector< vector<long double> > Gamma(Vec4 a, Vec4 b)
{
	vector< vector<long double> > G;
	vector<long double> g;
	Vec4 p = a + b;
	long double E = p.mCalc();
	long double gamma = (p.e())/E;
	long double beta = sqrt(1-1/(gamma*gamma));
	long double nx = p.px()/(E*gamma*beta);
	long double ny = p.py()/(E*gamma*beta);
	long double nz = p.pz()/(E*gamma*beta);

	long double g0[] = {gamma, gamma*beta*nx, gamma*beta*ny, gamma*beta*nz};
	long double g1[] = {g0[1], 1+(gamma-1)*nx*nx, (gamma-1)*nx*ny, (gamma-1)*nx*nz};
	long double g2[] = {g0[2], g1[2], 1+(gamma-1)*ny*ny, (gamma-1)*ny*nz};
	long double g3[] = {g0[3], g1[3], g2[3], 1+(gamma-1)*nz*nz};
	
	for (int i = 0; i<4; i++) {
		g.push_back (g0[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g1[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g2[i]);
	}
	G.push_back (g); g.clear();

	for (int i = 0; i<4; i++) {
		g.push_back (g3[i]);
	}
	G.push_back (g); g.clear();

	return G;
}
*/

//Declare classes

class Particle_list
{
		vector<Vec4> p, pb, n, nb;
	public:

		// Sorts an eventlist into p, pb, n, nb
		void sort (Event event)
		{
			for (int i = 0; i < event.size(); ++i)
			{
				if (event[i].isFinal()) 
				{
					switch (event[i].id())
					{
						case 2212:
							p.push_back (event[i].p());
							break;
						case -2212:
							pb.push_back (event[i].p());
							break;
						case 2112:
							n.push_back (event[i].p());
							break;
						case -2112:
							nb.push_back (event[i].p());
							break;
					}
				}
			}
		}
		
		vector<Vec4> pret(){return p;}
		vector<Vec4> pbret(){return pb;}
		vector<Vec4> nret(){return n;}
		vector<Vec4> nbret(){return nb;}
};

// Must be a Vec4 for now...
class Particle_def
{
		vector<Vec4> part;
		//vector<Vec4> part_hidden;
		vector<int> partId;
		vector<long double> weight;
	public:
		// functions to extract elements in the class
		vector<Vec4> vec4() { return part; }
		vector<int> Extract_Id() { return partId; }
		vector<long double> see_weight() {return weight; }
		
		// Use a vector filled with Vec4 to fill the Particle class
		void set(vector<Vec4> ret, int add = 0) // insert iEvent here if you want to include all finalstate events
		{
		  if (add == 0) { part = ret; }
		  else
		    {
		      for (int i = 0; i < abs((ret).size()); i++)
			{
			  part.push_back((ret)[i]);
			}
		    }
		}
				
		void set_wW(vector<Vec4> ret, vector<long double> weighted, int add = 0) // insert iEvent here if you want to include all finalstate events
		{
		  if (add == 0) { part = ret; weight = weighted; }
		  else
		    {
		      for (int i = 0; i < abs((ret).size()); i++)
			{
			  part.push_back((ret)[i]);
			  weight.push_back(weighted[i]);
			}
		    }
		}

		void set_read(char* file)
		{
		  ifstream infile;
		  infile >> setprecision(64) >> fixed;
		  string line;
		  infile.open(file);
		  if (infile)
		    {
		      while ( getline(infile, line) )
			{
			  if (!line.empty())
			    {
			      Vec4 temp;
			      int Id;
			      long double x, y, z, t;
			      istringstream iss(line);
			      iss >> Id >> x >> y >> z >> t;
			      partId.push_back(Id);
			      temp.p(x,y,z,t);
			      part.push_back(temp);
			    }
			}
		    }
		  else cout << "Unable to open file." << endl;
		}

		// read in and set a particle from data
		// assign each particle a weight per Event, N_tot is number of times a each baryon pair is used to make a deutron per Event
		void assign_weights_from_file(char* file)
		{
		  ifstream infile;
		  infile >> setprecision(64) >> fixed;
		  string line;
		  vector<long double> pT_bins;
		  vector<long double> bin_weight;
		  infile.open(file);
		  while ( getline(infile, line) )
		    {
		      if (!line.empty())
			{
			  long double x, y;
			  istringstream iss(line);
			  iss >> x >> y;
			  pT_bins.push_back(x); bin_weight.push_back(y);
			}
		    }
		  pT_bins.insert(pT_bins.begin(), pT_bins[0]-(pT_bins[1]-pT_bins[0]));
		  
		  // test this weight code
		  vector<long double> tWeight((part).size(), 0);
		  for (int i = 0; i < abs((part).size()); i++)
		    {
		      for (int j = 0; j < abs(pT_bins.size()-1); j++)
			{
			  if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) )
			    {
			      tWeight[i] = bin_weight[j];
			    }
			}
		      if ( (part[i].pT() < pT_bins[0]) ) { tWeight[i] = 1; }
		      if ( (part[i].pT() >= pT_bins.back() ) ) { tWeight[i] = 1; }
		    }
		  weight = tWeight;
		}

		void weight_push_back(long double w)
		{
		  weight.push_back(w);
		}

		// print pT for this particle
		void print_pT (char* location, int iEvent = 0, int nEvent = 1)
		{
			ofstream outputFile;
			if (iEvent == 0) outputFile.open(location, ios::out );
			else outputFile.open(location, ios::out | ios::app );
			outputFile << setprecision(64) << fixed;
			for (int i = 0; i < abs((part).size()); i++)
			{
			  if (abs((part)[i].rap()) < 0.5)
				{
				  if ( !weight.empty() )
				    {
				      outputFile << iEvent << " "
						 << part[i].pT() << " " << weight[i] << endl;
				    }
				  else
				    {
				      outputFile << iEvent << " "
						 << part[i].pT() << endl;
				    }
				}
			}
			if (iEvent == nEvent-1) outputFile.close();
		}
		
		// print Vec4 for this particle
		void print_Vec4 (char* location, int iEvent=0, int nEvent=1)
		{
			ofstream outputFile;
			outputFile << setprecision(64) << fixed;
			if (iEvent == 0) outputFile.open(location);
			else outputFile.open(location, ios::out | ios::app);
			for (int i = 0; i < abs((part).size()); i++)
			{
			  outputFile
			    << (part)[i] << endl;
			}
			if (iEvent == nEvent-1) outputFile.close();

		}

		void cout_Vec4 ()
		{
			for (int i = 0; i < abs((part).size()); i++)
			{
				cout << (part)[i] << endl;
			}
		}

		// Coalescence model

		void coalescence(Particle_def partAi, Particle_def partBi, bool iD = false)
		{
			vector<Vec4> partA = partAi.vec4(); vector<Vec4> partB = partBi.vec4();
			vector<int> partAiId = partAi.Extract_Id();
			vector<int> partBiId = partBi.Extract_Id();
			if (iD == true)
			  {
			    int p = 0; int n = 0; int nInit = 0;
			    while ( p < abs( partAiId.size() ) )
			      {
				if (n == abs(partBiId.size()) ) {break;}
				while ( n < abs(partBiId.size()) )
				  {
				    if (partAiId[p] > partBiId[n]) {n++; nInit = n;}
				    else if (partAiId[p] < partBiId[n]) {p++; n = nInit;}
				    else if (partAiId[p] == partBiId[n])
				      {
					double k = 2 * p_abs((partA)[p], (partB)[n], (partA)[p].mCalc(), (partB)[n].mCalc());
					if (k < 0.2)
					  {
					    if ( ( !(partAi.see_weight()).empty() ) & ( !(partBi.see_weight()).empty() ) )
					      {
						weight.push_back(partAi.see_weight()[p]*partBi.see_weight()[n]);
						part.push_back ((partA)[p] + (partB)[n]);
					      }
					    else
					      {
						part.push_back ((partA)[p] + (partB)[n]);
						cout << "No weights" << endl;
					      }
					  }
					n++;
				      }
				  }
			      }
			  }
			else
			  {
			    for (int p = 0; p < abs((partA).size()); ++p) 
			      {
				for (int n = 0; n < abs((partB).size()); ++n)
				  {
				    double k = 2 * p_abs((partA)[p], (partB)[n], (partA)[p].mCalc(), (partB)[n].mCalc());
				    if (k < 0.2)
				      {
					if ( ( !(partAi.see_weight()).empty() ) & ( !(partBi.see_weight()).empty() ) )
					  {
					    weight.push_back(partAi.see_weight()[p]*partBi.see_weight()[n]);
					    part.push_back ((partA)[p] + (partB)[n]);
					  }
					else
					  {
					    part.push_back ((partA)[p] + (partB)[n]);
					  }
				      }
				  }
			      }
			  }
		}


		
		// Original

		// Introducing weights. To activate the weighting function, need to add an integer N_MC as the last input
		void x_sec(Particle_def partAi, Particle_def partBi, int N_MC = 0)
		{
		  if ( ( !(partAi.vec4()).empty() ) & ( !(partBi.vec4()).empty() ) )
		    {
		      weight.clear();
		      vector<Vec4> partA = partAi.vec4(); vector<Vec4> partB = partBi.vec4();
                      vector<long double> partA_weight, partB_weight, temp_weight, m1, m2;

		      partA_weight = partAi.see_weight();
		      partB_weight = partBi.see_weight();
		      // New take on the problem... Need to make a list which includes event ID and particle info and crosssection.
		      // 3 lists maybe to be nested into a single one?
		      vector<Vec4> particleStats1; // Save the particle stats so it can be associated with a particle ID.
		      vector<Vec4> particleStats2; // Save the particle stats so it can be associated with a particle ID.
		      vector<long double> pdA;

		      // pn -> dg
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = 0; n < abs((partB).size()); ++n)
			    {
			      long double k = 2 * p_abs((partB)[n], (partA)[p], (partB)[n].mCalc(), (partA)[p].mCalc());
			      long double pd = p_abs((partB)[n], (partA)[p], 0, md);
			      long double prob = xs_pn_dgamma(k) * sigma_0_inv;
			      for (int N = 0; N < N_MC; N++)
				{
				  long double r = ((long double) rand()/RAND_MAX);
				  if (prob > r) // If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partA_weight[p] * partB_weight[n]);
                                      //m1.push_back(mp), m2.push_back(mn);
				    }
				}
			    }
			}

		      // pn -> dpi0
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = 0; n < abs((partB).size()); ++n)
			    {
			      double k = 2 * p_abs((partB)[n], (partA)[p], (partB)[n].mCalc(), (partA)[p].mCalc());
			      double pd = p_abs((partB)[n], (partA)[p], mpi0, md);
			      double prob = xs_pn_dpi0(k) * sigma_0_inv; //**
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r)  // If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partA_weight[p] * partB_weight[n]);
				      //m1.push_back(mp), m2.push_back(mn);
				    }
				}
			    }
			}
		      
		      
		      // pn -> dpi+pi-
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = 0; n < abs((partB).size()); ++n)
			    {
			      double k = 2 * p_abs((partB)[n], (partA)[p], (partB)[n].mCalc(), (partA)[p].mCalc());
			      double s = ((partB)[n] + (partA)[p]).m2Calc();
			      double x = ((double) rand()/RAND_MAX);
			      double m_min = sqrt(s) - md; double m_max = 2*mpic;
			      //double mpipi = x*(m_max-m_min) + m_min;
			      double mpipi = sqrt(x*(m_max*m_max-m_min*m_min) + m_min*m_min);
			      double pd = p_abs((partB)[n], (partA)[p], mpipi, md); 
			      double prob = xs_pn_dpippim(k) * sigma_0_inv;
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) // If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partA_weight[p] * partB_weight[n]);
				      //m1.push_back(mp), m2.push_back(mn);
				    }
				}
			    }
			}
		      
		      
		      
		      // pn -> dpi0pi0
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = 0; n < abs((partB).size()); ++n)
			    {
			      double k = 2 * p_abs((partB)[n], (partA)[p], (partB)[n].mCalc(), (partA)[p].mCalc()); 
			      double s = ((partB)[n] + (partA)[p]).m2Calc();
			      double x = ((double) rand()/RAND_MAX);
			      double m_min = sqrt(s) - md; double m_max = 2*mpi0;
			      //double mpipi = x*(m_max-m_min) + m_min;
			      double mpipi = sqrt(x*(m_max*m_max-m_min*m_min) + m_min*m_min);
			      double pd = p_abs((partB)[n], (partA)[p], mpipi, md);
			      double prob = xs_pn_dpi0pi0(k) * sigma_0_inv;

			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) // If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partA_weight[p] * partB_weight[n]);
				      //m1.push_back(mp), m2.push_back(mn);
				    }
				}
			    }
			}

		      // pp -> d pi+
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = p+1; n < abs((partA).size()); ++n)
			    {
			      double k = 2 * p_abs((partA)[n], (partA)[p], (partA)[n].mCalc(), (partA)[p].mCalc());
			      double pd = p_abs((partA)[n], (partA)[p], md, mpic);
			      double prob = xs_pp_dpip(k) * sigma_0_inv; //**
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) // If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partA)[n]); 
				      temp_weight.push_back(partA_weight[p] * partA_weight[n]);
				      //m1.push_back(mp), m2.push_back(mp);
				    }
				}
			    }
			}
		      
		      // pp -> dpi+pi0
		      for (int p = 0; p < abs((partA).size()); ++p)
			{
			  for (int n = p+1; n < abs((partA).size()); ++n)
			    {
			      double k = 2 * p_abs((partA)[n], (partA)[p], (partA)[n].mCalc(), (partA)[p].mCalc());
			      double s = ((partA)[n] + (partA)[p]).m2Calc();
			      double x = ((double) rand()/RAND_MAX);
			      double m_min = sqrt(s) - md; double m_max = mpic + mpi0;
			      //double mpipi = x*(m_max-m_min) + m_min;
			      double mpipi = sqrt(x*(m_max*m_max-m_min*m_min) + m_min*m_min);
			      double pd = p_abs((partA)[n], (partA)[p], mpipi, md);
			      double prob = xs_pp_dpippi0(k) * sigma_0_inv;
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) 			// If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partA)[p]); 
				      pdA.push_back(pd);// particle A info
                                      particleStats2.push_back((partA)[n]); 
				      temp_weight.push_back(partA_weight[p] * partA_weight[n]);
				      //m1.push_back(mp), m2.push_back(mp);
				    }
				}
			    }
			}
		      
		      // nn -> dpi-
		      for (int p = 0; p < abs((partB).size()); ++p)
			{
			  for (int n = p+1; n < abs((partB).size()); ++n)
			    {
			      double k = 2 * p_abs((partB)[n], (partB)[p], (partB)[n].mCalc(), (partB)[p].mCalc());
			      double pd = p_abs((partB)[n], (partB)[p], md, mpic);
			      double prob = xs_nn_dpim(k) * sigma_0_inv; //**
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) 			// If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partB)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partB_weight[p] * partB_weight[n]);
				      //m1.push_back(mn), m2.push_back(mn);
				    }
				}
			    }
			}
		      

		      // nn -> dpi-pi0
		      // SPECIAL CASE: Process: nn -> d + pi + pi- (3 particle decay, no massless)
		      for (int p = 0; p < abs((partB).size()); ++p)
			{
			  for (int n = p+1; n < abs((partB).size()); ++n)
			    {	
			      double k = 2 * p_abs((partB)[n], (partB)[p], (partB)[n].mCalc(), (partB)[p].mCalc());
			      double s = ((partB)[n] + (partB)[p]).m2Calc();
			      double x = ((double) rand()/RAND_MAX);
			      double m_min = sqrt(s) - md; double m_max = mpic + mpi0;
			      //double mpipi = x*(m_max-m_min) + m_min;
			      double mpipi = sqrt(x*(m_max*m_max-m_min*m_min) + m_min*m_min);
			      double pd = p_abs((partB)[n], (partB)[p], mpipi, md);
			      double prob = xs_nn_dpimpi0(k) * sigma_0_inv;
			      for (int N = 0; N < N_MC; N++)
				{
				  double r = ((double) rand()/RAND_MAX); 
				  if (prob > r) 			// If the process happens, save all the different relevant information to two lists A and B
				    {
				      particleStats1.push_back((partB)[p]); 
				      pdA.push_back(pd);// particle A info
				      particleStats2.push_back((partB)[n]); 
				      temp_weight.push_back(partB_weight[p] * partB_weight[n]);
				      //m1.push_back(mn), m2.push_back(mn);
				    }
				}
			    }
			}
		      
		      // Run through the list of deuterons and trim away excess momentum after 
		      
		      if (abs(particleStats1.size()) >= 1 && abs(particleStats2.size()) >= 1)
			{
			  for (int l = 0; l < abs(particleStats1.size()); ++l)
			    {	
			      vector<double> pd;
			      // Trim away the escaping particle. This is a sketch to be implemented in a better way.
			      // in CMS and for a two particle final state, the magnitude momentum of each particle must be
			      double p_d = pdA[l];
			      // However, for threeparticle end states, things change.
			      // Choose two random numbers, one for theta, one for phi.
			      double theta = ((double) rand())/RAND_MAX;
			      double zs = theta*2-1;
			      double phi = ((double ) rand())/RAND_MAX;
			      double ts = phi*2*pi;
			      double xs = sqrt(1-zs*zs)*cos(ts);
			      double ys = sqrt(1-zs*zs)*sin(ts);
			      // Find the four vector
			      double E = sqrt(p_d*p_d + md*md);
			      double pdx = p_d*xs;
			      double pdy = p_d*ys;
			      double pdz = p_d*zs;
			      
			      pd.push_back (E); pd.push_back (pdx); pd.push_back (pdy); pd.push_back (pdz);
			      
			      // Multiply Gamma and thus boost to the lab system
			      //vector<vector <long double> > Ga = Gamma(particleStats1[l],particleStats2[l], m1[l], m2[l]);
			      vector<vector <long double> > Ga = Gamma(particleStats1[l],particleStats2[l], particleStats1[l].mCalc(), particleStats2[l].mCalc()); //m1[l], m2[l]);
			      double x = 0; double y = 0; double z = 0; double t = 0;
			      for (int j = 0; j < 4; j++) {
				t += Ga[j][0]*pd[j];
				x += Ga[j][1]*pd[j];
				y += Ga[j][2]*pd[j];
				z += Ga[j][3]*pd[j];
			      }
			      Vec4 momenta_d(x, y, z, t);
			      part.push_back (momenta_d);
			      weight.push_back(temp_weight[l]);
			    }
			}
		    }
		}
};



/*
// Test the different functions:
void p_abs_test(Vec4 a, Vec4 b)
{
	long double p = p_abs(a, b, a.mCalc(), b.mCalc());
	cout << (a+b).m2Calc() << " after: " << pow(sqrt(p*p + a.m2Calc()) + sqrt(p*p + b.m2Calc()),2) << endl;
}


void Gamma_test(Vec4 a, Vec4 b)
{
	double theta = ((double) rand()/RAND_MAX)*2*pi;
	double phi = ((double) rand()/RAND_MAX)*pi;

	vector<vector <double> > G = Gamma(a, b);
	double p_a = p_abs(a, b, a.mCalc(), b.mCalc());
	double Ea = sqrt(p_a*p_a + a.m2Calc());
	double pax = p_a*cos(theta)*sin(phi);
	double pay = p_a*sin(theta)*sin(phi);
	double paz = p_a*cos(phi);

	vector<double> pa;
	pa.push_back(Ea); pa.push_back(pax); pa.push_back(pay); pa.push_back(paz); 

	double Eb = sqrt(p_a*p_a + b.m2Calc());
	vector<double> pb;
	pb.push_back(Eb); pb.push_back(-pax); pb.push_back(-pay); pb.push_back(-paz); 

	//cout << "p_abs: " << p_a*p_a << " p_abs after: " << pax*pax + pay*pay + paz*paz << endl;

	//cout << (a+b).m2Calc() << " -> " << pow(pa[0]+pb[0],2) - pow(pa[1]+pb[1],2) - pow(pa[2]+pb[2],2) - pow(pa[3]+pb[3],2) << endl;

	double xa = 0; double ya = 0; double za = 0; double ta = 0;
	double xb = 0; double yb = 0; double zb = 0; double tb = 0;
	for (int j = 0; j < 4; j++) 
	{
		ta += G[j][0]*pa[j];
		xa += G[j][1]*pa[j];
		ya += G[j][2]*pa[j];
		za += G[j][3]*pa[j];

		tb += G[j][0]*pb[j];
		xb += G[j][1]*pb[j];
		yb += G[j][2]*pb[j];
		zb += G[j][3]*pb[j];
	}
	Vec4 pa_test(xa, ya, za, ta);
	Vec4 pb_test(xb, yb, zb, tb);
	
	cout << "before: " << (a+b).m2Calc() << " after: " << (pa_test + pb_test).m2Calc() << endl;
	cout << "a: " << a << " b: " << b << endl;
	cout << "pa: " << pa_test << " pb: " << pb_test << endl;
	cout << "a+b: " << a+b << " pa+pb: " << pa_test + pb_test << endl;
}

*/

		// Make a nested list containing histogram weights with pT location. Need an external file to weight against
		/* vector< vector<double> > hist_weighting(char* file, int resolution) */
		/* { */
		/*   weight.clear(); */
		/*   ifstream infile; */
		/*   string line; */
		/*   vector<double> pT_bins, data_hist, model_hist, bin_weight; */
		/*   double weightNum, tempWeight; */
		/*   infile.open(file); */
		/*   while ( getline(infile, line) ) */
		/*     { */
		/*       if (!line.empty()) */
		/* 	{ */
		/* 	  double x, y; */
		/* 	  istringstream iss(line); */
		/* 	  iss >> x >> y; */
		/* 	  pT_bins.push_back(x); */
		/* 	  data_hist.push_back(y); */
		/* 	} */
		/*     } */
		/*   pT_bins.insert(pT_bins.begin(), pT_bins[0]-(pT_bins[1]-pT_bins[0])); */
		/*   cout << "pT bin size: " << pT_bins.size() << endl; */
		/*   for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/*     { */
		/*       double N_bj = 0; */
		/*       for (int i = 0; i < abs((part).size()); i++) */
		/* 	{ */
		/* 	  if ( (abs((part)[i].rap()) < 0.5) & (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 	    { */
		/* 	      N_bj += 1; */
		/* 	    } */
		/* 	} */
		/*       tempWeight = N_bj/(2*pi*pT_bins[j+1]*(pT_bins[j+1]-pT_bins[j])*resolution); */
		/*       cout << "tempWeight: " << tempWeight << " bin: " << j << endl; */
		/*       model_hist.push_back(tempWeight); */
		/*     } */
		/*   for (int i = 0; i < abs(model_hist.size()); i++) */
		/*     { */
		/*       weightNum = data_hist[i]/(model_hist[i]); //0.861 */
		/*       cout << "Weights: " << weightNum << endl; */
		/*       bin_weight.push_back(weightNum); */
		/*     } */
		/*   vector < vector <double> > outPut; */
		/*   outPut.push_back(bin_weight); outPut.push_back(pT_bins); */
		/*   return outPut; */
		/* } */

		/* void weights_from_variable(vector < vector < double > > weight_pT_bins) */
		/* { */
		/*   vector < double > bin_weight; */
		/*   vector < double > pT_bins; */
		/*   bin_weight = weight_pT_bins[0]; */
		/*   pT_bins = weight_pT_bins[1]; */
		/*   weight.clear(); */
		/*   for (int i = 0; i < abs((part).size()); i++) */
		/*     { */
		/*       for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	{ */
		/* 	  if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 	    { */
		/* 	      weight.push_back(bin_weight[j]); */
		/* 	    } */
		/* 	} */
		/*       if ( (part[i].pT() < pT_bins[0]) ) { weight.push_back(1); } */
		/*       if ( (part[i].pT() >= pT_bins.back() ) ) {weight.push_back(1); } */
		/*     } */
		/* } */

		// read in and set a particle from data
		// assign weights conforming to an external datafile. The file which is read is the data we need to conform to.		
		/* void conform_weights(char* file, int resolution) */
		/* { */
		/*   weight.clear(); */
		/*   ifstream infile; */
		/*   string line; */
		/*   vector<double> pT_bins; */
		/*   vector<double> data_hist; */
		/*   vector<double> model_hist; */
		/*   vector<double> bin_weight; */
		/*   double weightNum; */
		/*   infile.open(file); */
		/*   while ( getline(infile, line) ) */
		/*     { */
		/*       if (!line.empty()) */
		/* 	{ */
		/* 	  double x, y; */
		/* 	  istringstream iss(line); */
		/* 	  iss >> x >> y; */
		/* 	  pT_bins.push_back(x); */
		/* 	  data_hist.push_back(y); */
		/* 	} */
		/*     } */
		/*   pT_bins.insert(pT_bins.begin(), pT_bins[0]-(pT_bins[1]-pT_bins[0])); */
		/*   for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/*     { */
		/*       double N_bj = 0; */
		/*       for (int i = 0; i < abs((part).size()); i++) */
		/* 	{ */
		/* 	  if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 	    { */
		/* 	      N_bj += 1; */
		/* 	    } */
		/* 	} */
		/*       model_hist.push_back(N_bj/(2*pi*pT_bins[j+1]*(pT_bins[j+1]-pT_bins[j])*resolution)); */
		/*     } */
		/*   for (int i = 0; i < abs(model_hist.size()); i++) */
		/*     { */
		/*       weightNum = data_hist[i]/(model_hist[i]/0.861); */
		/*       bin_weight.push_back(weightNum); */
		/*     } */
		/*   for (int i = 0; i < abs((part).size()); i++) */
		/*     { */
		/*       for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	{ */
		/* 	  if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 	    { */
		/* 	      weight.push_back(bin_weight[j]); //Possible multiply with 0.861 */
		/* 	    } */
		/* 	} */
		/*       if ( (part[i].pT() < pT_bins[0]) ) { weight.push_back(1); } */
		/*       if ( (part[i].pT() >= pT_bins.back() ) ) {weight.push_back(1); } */
		/*     } */
		/* } */
		
		/* void assign_weights(vector<double> pT_bins, int N_tot = 1) */
		/* { */
		/*   vector<double> bin_weight; */
		/*   bin_weight.clear(); */
		/*   vector<int> temp_hist; */
		/*   temp_hist.clear(); */
		/*   if (weight.empty()) */
		/*     { */
		/*       for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	{ */
		/* 	  double N_bj = 0; */
		/* 	  for (int i = 0; i < abs((part).size()); i++) */
		/* 	    { */
		/* 	      if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 		{ */
		/* 		  N_bj += 1; */
		/* 		} */
		/* 	    } */
		/* 	  bin_weight.push_back(N_bj/N_tot); */
		/* 	} */
		/*       for (int i = 0; i < abs(part.size()); i++) */
		/* 	{ */
		/* 	  for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	    { */
		/* 	      if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 		{ */
		/* 		  //			      cout << bin_weight[j] << endl; */
		/* 		  weight.push_back(bin_weight[j]); */
		/* 		  //			      cout << part[i].pT() << " " << weight.back() << endl; */
		/* 		  break; */
		/* 		} */
		/* 	      else if ( (part[i].pT() < pT_bins[0]) ) */
		/* 		{ */
		/* 		  weight.push_back(1); */
		/* 		  //			      cout << part[i].pT() << " " << weight.back() << endl; */
		/* 		  break; */
		/* 		} */
		/* 	      else if (  (part[i].pT() >= pT_bins.back()) ) */
		/* 		{ */
		/* 		  weight.push_back(1); */
		/* 		  break; */
		/* 		} */
		/* 	      else if ( isnan( part[i].pT()) ) */
		/* 		{ */
		/* 		  weight.push_back(1); // 0? */
		/* 		  break; */
		/* 		} */
		/* 	    } */
		/* 	} */
		/*     } */
		/*   else if (!weight.empty()) */
		/*     { */
		/*       for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	{ */
		/* 	  double N_bj = 0; */
		/* 	  for (int i = 0; i < abs((part).size()); i++) */
		/* 	    { */
		/* 	      if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 		{ */
		/* 		  N_bj += 1; */
		/* 		} */
		/* 	    } */
		/* 	  bin_weight.push_back(N_bj/N_tot); */
		/* 	} */
		/*       for (int i = 0; i < abs(part.size()); i++) */
		/* 	{ */
		/* 	  for (int j = 0; j < abs(pT_bins.size()-1); j++) */
		/* 	    { */
		/* 	      if ( (part[i].pT() >= pT_bins[j]) & (part[i].pT() < pT_bins[j+1]) ) */
		/* 		{ */
		/* 		  //			      cout << bin_weight[j] << endl; */
		/* 		  weight[i] = weight[i]*bin_weight[j]; //weight[i]*bin_weight[j]; */
		/* 		  //			      cout << part[i].pT() << " " << weight.back() << endl; */
		/* 		  break; */
		/* 		} */
		/* 	      else if ( isnan( part[i].pT()) ) */
		/* 		{ */
		/* 		  weight[i] = 1; //0? */
		/* 		  cout << "WARNING: ISNAN!" << endl; */
		/* 		  break; */
		/* 		} */
		/* 	    } */
		/* 	} */
		/*     } */
		/*   if (abs(part.size()) != abs(weight.size())) */
		/*     { */
		/*       //		      cout << part.size() << " " << weight.size() << endl; */
		/*       for (int k = 0; k < abs(part.size()); k++) */
		/* 	{ */
		/* 	  cout << part[k].pT() << endl; */
		/* 	} */
		/*     } */
		/* } */


#endif
