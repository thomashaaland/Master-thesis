import numpy as np
import matplotlib.pylab as plt

pi = np.pi

# Define program to read file

# string is the filename, numberOfHeaders is number of headerlines the file starts with, 
# 	variable is which variable you want numbered
def extractData(string,numberOfHeaders,variable):
	f = open(string,'r')
	i = 0
	gold = []
	while i < numberOfHeaders-1:
		f.readline()
		i += 1
	for line in f:
		line = line.strip()
		coloumn = line.split()
		gold += [float(coloumn[variable])]
	return gold

# Normalised to the initial points
def binCount(data, x):
	#Set up the bins
	dx = np.zeros(len(x)-1)
	new_x = np.zeros(len(x)-1)
	for n in range(0,len(x)-1):
		dx[n] = x[n+1]-x[n]
		new_x[n] = (x[n]+x[n+1])/2
	y = np.zeros(len(x)-1)
	for j in range(0,len(x)-1):
		for i in data:
			bot = x[j]
			top = x[j+1]
			if i < top and i >= bot:
				y[j] += 1
		y[j] = y[j]/(2*pi*new_x[j]*dx[j]*4000000)   # size data should be 1M
	return new_x, y
"""	
# ALICE proton
x = extractData("ALICE_p_data_7.dat",0,0)
y = extractData("ALICE_p_data_7.dat",0,1)
error = extractData("ALICE_p_data_7.dat",0,2)

plt.errorbar(x,np.array(y),error)

# ALICE antiproton
x = extractData("ALICE_pbar_data_7.dat",0,0)
y = extractData("ALICE_pbar_data_7.dat",0,1)
error = extractData("ALICE_pbar_data_7.dat",0,2)

plt.semilogy(x,y,'s')

# ALICE deutron
x = extractData("ALICE_d_data_7.dat",0,0)
y = extractData("ALICE_d_data_7.dat",0,1)
error = extractData("ALICE_d_data_7.dat",0,2)

plt.errorbar(x,np.array(y),error)

#Plot proton
pTP = extractData("outProjectP4P.dat",0,0)
#x = extractData("ALICE_p_data_7.dat",0,0)
x = np.linspace(0.25,3.5,101)
y = np.array(binCount(pTP,x)[1])
x = np.array(binCount(pTP,x)[0])

plt.semilogy(x,y)
plt.ylabel('1/N * d^2N/pTdpTdy')
plt.xlabel('p [GeV]')

#Plot antiproton
pTP = extractData("outProject4Pb.dat",0,0)
x = np.linspace(0.25,3.5,31)
y = np.array(binCount(pTP,x)[1])
x = np.array(binCount(pTP,x)[0])

plt.semilogy(x,y*0.831)
plt.semilogy(x,y*0.861)
plt.semilogy(x,y*0.891)

plt.ylabel('1/N * d^2N/pT dpT dy')
plt.xlabel('p [GeV]')

plt.legend(['Experiment','f=0.831','f=0.861','f=0.891'])
#plt.legend(['ALICE deuterium','ALICE anti deuterium','Antideuterium'])
plt.title('Spectrum, |y|<0.5, Pythia 8.186, 100 bins')
plt.axis([0,4,10**(-5),10**(-1)])
plt.show()
"""

#ALICE antideutron
x = extractData("ALICE_dbar_data_7.dat",0,0)
y = extractData("ALICE_dbar_data_7.dat",0,1)
error = extractData("ALICE_dbar_data_7.dat",0,2)

plt.errorbar(x,y,error)

# Plot Antideutron

pT = extractData("outDeuteron.dat",0,0)
x = np.linspace(0.1,4,31)
y = np.array(binCount(pT,x)[1])
x = np.array(binCount(pT,x)[0])
plt.semilogy(x,y*0.861,'o')

pT = extractData("outDeuteron_cross.dat",0,0)
x = np.linspace(0.1,4,31)
y = np.array(binCount(pT,x)[1])
x = np.array(binCount(pT,x)[0])
plt.semilogy(x,y*0.861,'s')

plt.legend(['Experiment', 'Coalescence', 'Crosssection fit'])
plt.title('Spectrum, |y|<0.5, Pythia 8.186, 30 bins')
plt.ylabel('1/N * d^2N/pT dpT dy')
plt.xlabel('p [GeV]')
plt.show()
