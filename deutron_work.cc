// Headers and namespace
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Pythia8/Pythia.h"
#include <ctime>
#include <vector>
#include "crossSections.h"
using namespace Pythia8;
using namespace std;

// Define and Declare functions and methods

// p_abs gives the absolute value of the d and x particles 3momentum in CMS. Must send in two Vec4 containing the relevant particles. 
// m1 is first particle mass, m2 is second particle mass and f_states are the number of final states. For three-body decay, m1 = m12 and m2 = m3.
// However, the sequence of m1 and m2 should not matter.
double p_abs(Vec4 a, Vec4 b, double m1, double m2)
{

	double s = (a + b).m2Calc();
	double p_abs = sqrt((s - pow((m1 + m2),2))*(s - pow((m1 - m2),2))) / (2*sqrt(s));

//	if (isnan(p_abs)) cout << "You got NaN: " << (s - pow((m1 + m2),2))*(s - pow((m1 - m2),2))/(2*sqrt(s)) << " (m1+m2)^2: " << pow(m1+m2,2) << " s: " << s << endl; 

	return p_abs;	
}

// constructs the Gamma matrix boosting from CMS to lab, a and b is two Vec4 vectors in the lab system.
vector<vector<double> > Gamma(Vec4 a, Vec4 b)
{
	vector< vector<double> > G;
	vector<double> g;
	Vec4 p = a + b;
	double E = p.mCalc();
	double gamma = (p.e())/E;
	double beta = sqrt(1-1/gamma);
	double nx = p.px()/(E*gamma*beta);
	double ny = p.py()/(E*gamma*beta);
	double nz = p.pz()/(E*gamma*beta);

	double g0[] = {gamma, gamma*beta*nx, gamma*beta*ny, gamma*beta*nz};
	double g1[] = {g0[1], 1+(gamma-1)*nx*nx, (gamma-1)*nx*ny, (gamma-1)*nx*nz};
	double g2[] = {g0[2], g1[2], 1+(gamma-1)*ny*ny, (gamma-1)*ny*nz};
	double g3[] = {g0[3], g1[3], g2[3], 1+(gamma-1)*nz*nz};
	
	for (int i = 0; i<4; i++) {
		g.push_back (g0[i]);
	}
	G.push_back (g);

	for (int i = 0; i<4; i++) {
		g.push_back (g1[i]);
	}
	G.push_back (g);

	for (int i = 0; i<4; i++) {
		g.push_back (g2[i]);
	}
	G.push_back (g);

	for (int i = 0; i<4; i++) {
		g.push_back (g3[i]);
	}
	G.push_back (g);

	return G;
}

// Declare Global constants
double sigma_0_inv = 2.58 / 1000000; // units is microbarn
double pi = 3.14159265359;

//int main(int argc, char* argv[]) {
int main() {
	
	// Generate!
	Pythia pythia;
	Event& event = pythia.event;

	// Commandfile! :-D	
	pythia.readFile("project4.cmnd");
	// Extract number of generated events from command file
	int nEvent = pythia.mode("Main:numberOfEvents");

	// Prepare a file to write to

	// AntiProton
	ofstream outputFileAntiProton;
	outputFileAntiProton.open("outAntiProton.dat");

	// AntiDeuterium
	ofstream outputFile;
	outputFile.open("outDeuteron.dat");
	ofstream outputFile_cross;
	outputFile_cross.open("outDeuteron_cross.dat");

	// Initialize Pythia!
	pythia.init();

	// --------BEGIN EVENT LOOP--------------//
	// --------------------------------------//
	for (int iEvent = 0; iEvent < nEvent; ++iEvent)
	{
		if (!pythia.next()) continue;

		// Production of antideuterium
		// Check if there are any antiprotons and antineutrons and try to make antideuterons
		vector<Particle> protons; // A list of all the protons in their particle state
		vector<Particle> neutrons; // A list of all the neutrons in their particle state
		vector<Vec4> deuterons_coalesc; // List to be filled by neutrons using the coalescence model, its Vec4 so must use Vec4 commands
		vector<Vec4> deuterons_cross;  // List to be filled by deutrons using cross section fit, its Vec4 so must use Vec4 commands
		vector<Vec4> deuterons_list; // List to be filled with deuterons after trimming away superflous particles

		// antiProton spectra
		int ipb = 0;
		double yb_sq = 0;

		cout << "1" << endl;
		for (int i = 0; i < event.size(); ++i)
		{
			if (event[i].isFinal())
			{
				if (pythia.event[i].id() == -2212)
				{ 
					ipb = i;
					yb_sq = pythia.event[ipb].y()*pythia.event[ipb].y();
					if (yb_sq < 0.5*0.5) outputFileAntiProton << event[ipb].pT() << endl;
				}
			}
		}
		cout << "Not 1" << endl;


		cout << "2" << endl;
		// Create lists of protons and neutrons		
		for (int i = 0; i < event.size(); ++i)
		{
			if (event[i].isFinal())
			{
				switch (event[i].id())
				{
					case -2212:
						protons.push_back (event[i]);
						break;
					case -2112:
						neutrons.push_back (event[i]);
						break;
				}
			}
		}
		cout << "Not 2" << endl;
		
		// Create Deuterons using the coalescence model.
		// NOTE: Switched to pAbs2() from m2Calc() since I should look at p3
		// NOTE2: For best accuracy, I should use the p_abs function, which finds the absolute value of p in CMS, using pAbs2() is WRONG!!

		cout << "3" << endl;
		for (int p = 0; p < abs(protons.size()); ++p) 
		{
			for (int n = 0; n < abs(neutrons.size()); ++n)
			{
				double k = 2 * p_abs(neutrons[n].p(), protons[p].p(), neutrons[n].mCalc(), protons[p].mCalc());
				if (k < 0.2)
				{
					deuterons_coalesc.push_back (protons[p].p() + neutrons[n].p());
				}
			}
		}
		cout << "Not 3" << endl;

		

		// New take on the problem... Need to make a list which includes event ID and particle info and crosssection.
		// 3 lists maybe to be nested into a single one?
		vector<int> eventIdA; // Which event is this? This number tells specifically which particle is used in the process.
		vector<int> eventIdB; // Which event is this? This number tells specifically which particle is used in the process.
		vector<Particle> particleStatsA; // Save the particle stats so it can be associated with a particle ID.
		vector<Particle> particleStatsB; // Save the particle stats so it can be associated with a particle ID.
		vector<long double> crossSectionA; // Need to save the probability for the event to happen.
		vector<long double> crossSectionB; // Need to save the probability for the event to happen.
		vector<bool> booleanA;		// To keep track of deletions
		vector<bool> booleanB;		// To keep track of deletions
		vector<double> pdA;
		vector<double> pdB;


		// This is used many times, so should turn this into a method
		// k is |p_1 - p_2|_COM initial particles(!!)


		cout << "4" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = 0; n < abs(neutrons.size()); ++n)
			{
				double k = 2 * p_abs(neutrons[n].p(), protons[p].p(), neutrons[n].p().mCalc(), protons[p].p().mCalc());
				double pd = p_abs(neutrons[n].p(), protons[p].p(), md, 0);
				long double prob = xs_pn_dgamma(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r) 			// If the process happens, save all the different relevant information to two lists A and B
				{	
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 4" << endl;

		cout << "5" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = 0; n < abs(neutrons.size()); ++n)
			{
				double k = 2 * p_abs(neutrons[n].p(), protons[p].p(), neutrons[n].mCalc(), protons[p].mCalc());
				double pd = p_abs(neutrons[n].p(), protons[p].p(), md, mpi0);
				long double prob = xs_pn_dpi0(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 5" << endl;

		cout << "6" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = 0; n < abs(neutrons.size()); ++n)
			{
				long double k = 2 * p_abs(protons[p].p(), neutrons[n].p(), protons[p].mCalc(), neutrons[n].mCalc());
				long double s = (neutrons[n].p() + protons[p].p()).m2Calc();
				long double x = ((long double) rand()/RAND_MAX);
				long double m_max = sqrt(s) - md; long double m_min = 2*mpic;
				long double mpipi = x*(m_max*m_max-m_min*m_min) + m_min*m_min;
				long double pd = p_abs(neutrons[n].p(), protons[p].p(), md, mpipi); 
				long double prob = xs_pn_dpippim(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 6" << endl;

		cout << "7" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = 0; n < abs(neutrons.size()); ++n)
			{
				double k = 2 * p_abs(protons[p].p(), neutrons[n].p(), protons[p].mCalc(), neutrons[n].mCalc()); 
				long double s = (neutrons[n].p() + protons[p].p()).m2Calc();
				long double x = ((long double) rand()/RAND_MAX);
				long double m_max = sqrt(s) - md; long double m_min = 2*mpi0;
				long double mpipi = x*(m_max*m_max-m_min*m_min) + m_min*m_min;
				long double pd = p_abs(neutrons[n].p(), protons[p].p(), md, mpipi);
				long double prob = xs_pn_dpi0pi0(k) * sigma_0_inv;
				long double r = ((double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 7" << endl;

		cout << "8" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = p+1; n < abs(protons.size()); ++n)
			{
				double k = 2 * p_abs(protons[n].p(), protons[p].p(), protons[n].mCalc(), protons[p].mCalc());
				double pd = p_abs(protons[n].p(), protons[p].p(), md, mpic);
				long double prob = xs_pp_dpip(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(protons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 8" << endl;


		cout << "9" << endl;
		for (int p = 0; p < abs(protons.size()); ++p)
		{
			for (int n = p+1; n < abs(protons.size()); ++n)
			{
				double k = 2 * p_abs(protons[n].p(), protons[p].p(), protons[n].mCalc(), protons[p].mCalc());
				long double s = (protons[n].p() + protons[p].p()).m2Calc();
				long double x = ((long double) rand()/RAND_MAX);
				long double m_max = sqrt(s) - md; long double m_min = mpic + mpi0;
				long double mpipi = x*(m_max*m_max-m_min*m_min) + m_min*m_min;
				long double pd = p_abs(protons[n].p(), protons[p].p(), md, mpipi);
				long double prob = xs_pp_dpippi0(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(protons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(protons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 9" << endl;

		cout << "10" << endl;
		for (int p = 0; p < abs(neutrons.size()); ++p)
		{
			for (int n = p+1; n < abs(neutrons.size()); ++n)
			{
				double k = 2 * p_abs(neutrons[n].p(), neutrons[p].p(), neutrons[n].mCalc(), neutrons[p].mCalc());
				double pd = p_abs(neutrons[n].p(), neutrons[p].p(), md, mpic);
				long double prob = xs_nn_dpim(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(neutrons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}
			}
		}
		cout << "Not 10" << endl;

		cout << "11" << endl;
		for (int p = 0; p < abs(neutrons.size()); ++p)
		{
			for (int n = p+1; n < abs(neutrons.size()); ++n)
			{	
				double k = 2 * p_abs(neutrons[n].p(), neutrons[p].p(), neutrons[n].mCalc(), neutrons[p].mCalc());
				long double s = (neutrons[n].p() + neutrons[p].p()).m2Calc();
				long double x = ((long double) rand()/RAND_MAX);
				long double m_max = sqrt(s) - md; long double m_min = mpic + mpi0;
				long double mpipi = x*(m_max*m_max-m_min*m_min) + m_min*m_min;
				long double pd = p_abs(neutrons[n].p(), neutrons[p].p(), md, mpipi);
				long double prob = xs_nn_dpimpi0(k) * sigma_0_inv;
				long double r = ((long double) rand()/RAND_MAX); 
				if (prob > r)
				{
					eventIdA.push_back(p); particleStatsA.push_back(neutrons[p]); 
					crossSectionA.push_back(prob); booleanA.push_back(true); pdA.push_back(pd);// particle A info
					eventIdB.push_back(n); particleStatsB.push_back(neutrons[n]); 
					crossSectionB.push_back(prob); booleanB.push_back(true); pdB.push_back(pd);// particle B info
				}

			}

		}
		cout << "Not 11" << endl;

		// Need to sort through the listed particles and check to see if any1 has been used twice and save the events to a list.
		vector<int> twoTimersA;
		vector<int> twoTimersB;

		cout << "12" << endl;
		for (int a = 0; a < abs(eventIdA.size()); ++a)
		{	
			if (booleanA[a] == true)
			{
				twoTimersA.push_back(a);
				int a2 = a+1;
				if (a2 < abs(eventIdA.size()))
				{
					for (a2 = a+1; a2 < abs(eventIdA.size()); ++a2)
					{
						if (booleanA[a2] == true && eventIdA[a] == eventIdA[a2]) twoTimersA.push_back(a2);
					}
				}
				if (twoTimersA.size() == 1) twoTimersA.pop_back();
				// I want to select one process out of the available, so must make a list of the probabilities
				double choice = ((double) rand()/RAND_MAX);
				vector<long double> probabilities;
				probabilities.push_back(0);
				long double cross = 0;
				for (int t1 = 0; t1 < abs(twoTimersA.size()); ++t1)
				{
					cross += crossSectionA[twoTimersA[t1]];
					probabilities.push_back(cross);
				}
				for (int t2 = 0; t2 < abs(probabilities.size())-1; ++t2)
				{
					if (choice >= probabilities[t2]/cross && choice < probabilities[t2+1]/cross){}
					else booleanA[twoTimersA[t2]] = false;
				}
			}

		}
		cout << "Not 12" << endl;

		cout << "13" << endl;
		for (int a = 0; a < abs(eventIdB.size()); ++a)
		{	
			if (booleanB[a] == true)
			{
				twoTimersB.push_back(a);
				int a2 = a+1;
				if (a2 < abs(eventIdB.size()))
				{
					for (a2 = a+1; a2 < abs(eventIdB.size()); ++a2)
					{
						if (booleanB[a2] == true && eventIdB[a] == eventIdB[a2]) twoTimersB.push_back(a2);
					}
				}
				if (twoTimersB.size() == 1) twoTimersB.pop_back();
				// I want to select one process out of the available, so must make a list of the probabilities
				double choice = ((double) rand()/RAND_MAX);
				vector<long double> probabilities;
				probabilities.push_back(0);
				long double cross = 0;
				for (int t1 = 0; t1 < abs(twoTimersB.size()); ++t1)
				{
					cross += crossSectionB[twoTimersB[t1]];
					probabilities.push_back(cross);
				}
				for (int t2 = 0; t2 < abs(probabilities.size())-1; ++t2)
				{
					if (choice >= probabilities[t2]/cross && choice < probabilities[t2+1]/cross){}
					else booleanB[twoTimersB[t2]] = false;
				}
			}
		}
		cout << "Not 13" << endl;

		cout << "14" << endl;
		
		int q = 0;
		while (q < abs(booleanA.size()))
		{
			if (booleanA[q] == 0)
			{
				booleanA[q] = booleanA.back(); 			eventIdA[q] = eventIdA.back(); 
				particleStatsA[q] = particleStatsA.back();	crossSectionA[q] = crossSectionA.back(); 
				pdA[q] = pdA.back();

				booleanA.pop_back(); 				eventIdA.pop_back(); 
				particleStatsA.pop_back(); 			crossSectionA.pop_back();
				pdA.pop_back();
			}
			else q++;
		}

		cout << "Not 14" << endl;

		cout << "15" << endl;
		int qq = 0;
		while (qq < abs(booleanB.size()))
		{
			if (booleanB[qq] == 0)
			{
				booleanB[qq] = booleanB.back(); 			eventIdB[qq] = eventIdB.back(); 
				particleStatsB[qq] = particleStatsB.back();	crossSectionB[qq] = crossSectionB.back(); 
				pdB[qq] = pdB.back();

				booleanB.pop_back(); 				eventIdB.pop_back(); 
				particleStatsB.pop_back(); 			crossSectionB.pop_back();
				pdB.pop_back();
			}
			else qq++;
		}
		cout << "Not 15" << endl;

		// Run through the list of deuterons and trim away excess momentum after 
		cout << "16" << endl;
	if (abs(particleStatsA.size()) >= 1 && abs(particleStatsB.size()) >= 1)
	{
		for (int l = 0; l < abs(particleStatsA.size()); ++l)
		{	
			vector<double> pd;
			vector<double> deuteron_p;
			// Trim away the escaping particle. This is a sketch to be implemented in a better way.
			// in CMS and for a two particle final state, the magnitude momentum of each particle must be
			double p_d = pdA[l];
			// However, for threeparticle end states, things change.
			// Choose two random numbers, one for theta, one for phi.
			double theta = ((double) rand()/RAND_MAX)*2*pi;
			double phi = ((double) rand()/RAND_MAX)*pi;
			// Find the four vector
			double E = sqrt(p_d*p_d + md*md);
			double pdx = p_d*cos(theta)*sin(phi);
			double pdy = p_d*sin(theta)*sin(phi);
			double pdz = p_d*cos(phi);
			
			pd.push_back (E); pd.push_back (pdx); pd.push_back (pdy); pd.push_back (pdz);

			// Multiply Gamma and thus boost to the lab system
			vector<vector <double> > Ga = Gamma(particleStatsA[l].p(),particleStatsB[l].p());
			double x = 0; double y = 0; double z = 0; double t = 0;
			for (int j = 0; j < 4; j++) {
				t += Ga[j][0]*pd[j];
				x += Ga[j][1]*pd[j];
				y += Ga[j][2]*pd[j];
				z += Ga[j][3]*pd[j];
			}
			Vec4 momenta_d(x, y, z, t);
			deuterons_list.push_back (momenta_d);
		}
	}
		cout << "Not 16" << endl;
		 
		// Finally sort according to rapidity and write to file

		cout << "17" << endl;
		for (int d = 0; d < abs(deuterons_coalesc.size()); ++d)
		{
			if (abs(deuterons_coalesc[d].rap()) < 0.5)
			{
				outputFile << abs(deuterons_coalesc[d].pT()) << endl;
			}
		}
		cout << "Not 17" << endl;

		// Finally sort according to rapidity and write to file
		cout << "18" << endl;
		for (int d = 0; d < abs(deuterons_list.size()); ++d)
		{
			if (abs(deuterons_list[d].rap()) < 0.5)
			{
				outputFile_cross << abs(deuterons_list[d].pT()) << endl;
			}
		}
		cout << "Not 18" << endl;

	}
	outputFileAntiProton.close();
	outputFile.close();
	outputFile_cross.close();
	return 0;
	}

	
