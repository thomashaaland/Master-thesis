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

# Assembly
def assemble(string):
	pT1 = extractData("./pythia_1/pythia8186/examples/"+string,0,0)
	pT2 = extractData("./pythia_2/pythia8186/examples/"+string,0,0)
	pT3 = extractData("./pythia_3/pythia8186/examples/"+string,0,0)
	pT4 = extractData("./pythia_4/pythia8186/examples/"+string,0,0)	
	return pT1 + pT2 + pT3 + pT4	


#ALICE antideutron
x = extractData("ALICE_dbar_data_7.dat",0,0)
y = extractData("ALICE_dbar_data_7.dat",0,1)
error = extractData("ALICE_dbar_data_7.dat",0,2)

plt.errorbar(x,y,error)

# Plot Antideutron

pT_1 = assemble("outDeuteron.dat")

x1 = np.linspace(0.1,4,31)
y1 = np.array(binCount(pT_1,x1)[1])
x1 = np.array(binCount(pT_1,x1)[0])
plt.semilogy(x1,y1*0.861,'o')



pT_2 = assemble("outDeuteron_cross.dat")

x2 = np.linspace(0.1,4,31)
y2 = np.array(binCount(pT_2,x2)[1])
x2 = np.array(binCount(pT_2,x2)[0])
plt.semilogy(x2,y2*0.861,'s')

plt.legend(['Experiment', 'Coalescence', 'Crosssection fit'])
plt.title('Spectrum, |y|<0.5, Pythia 8.186, 30 bins')
plt.ylabel('1/N * d^2N/pT dpT dy')
plt.xlabel('p [GeV]')
plt.show()
