import pickle
import numpy as np
import astropy.io.fits as fits
from input_parameters import pathto
#set limits on the LX band
UpLimit = 2.0 #keV
LowLimit = 0.5 #keV

#constants
hc = 12.398#keV A
eVkB = 11604.505 #1eV/kB = 11604.505 K
erg2keV = 6.242e8 

class run:
    def __init__(self, vol, ni, ne, temp, Fe, Ox, H, He):
        self.vol = vol
        self.ni = ni
	self.ne = ne
	self.temp = temp
	self.Fe = Fe
	self.Ox = Ox
	self.H = H
	self.He = He

def open_files(name):
	#arrays of quantities for gas particles within (0,240)kpc around a main galaxy

    filenames = ["Vol.p","ni.p", 'ne.p', "Temp.p","Fe.p", "Ox.p","H.p","He.p"]
    buff = {}
    for ii in range(8):
        with open(pathto+name+filenames[ii], 'rb') as f:
            if ii==4: #Fe
                buff[ii] = pickle.load(f)/0.01267

            elif ii==5: #Ox
                buff[ii] = pickle.load(f)/0.009618

            else:
	        buff[ii] = pickle.load(f)

    return run(buff[0],buff[1], buff[2], buff[3], buff[4], buff[5], buff[6], buff[7])

def lines(sim, line):
	""" pass class sim file and line fits file, returns epsilon"""
	#1. Read-off temp. cuts
	cut = line[1].data['kT']*eVkB*1e3 #result in K
	eps = np.zeros((4, len(sim.vol)))#initialize arrays for the readout of emissivities per element
	for i in range(len(cut)):
		if i==0:
			index = np.where(sim.temp < cut[i])[0]
		elif i==(len(cut)-1):
			index = np.where(sim.temp>cut[i])[0]
		else:
			index = np.where((sim.temp > cut[i-1])& (sim.temp < cut[i]))[0]
	#The position of the element of an array "cut" defines which subtable of line the data belongs to (except line[0]: empty, line[1]: initial parameter block)


	#2. Consider the lines with energies in the band (LowLimit, UpLimit)	
		line_energy = hc/line[i+2].data['Lambda'] # result in keV
		acceptID = np.where((line_energy > LowLimit) & (line_energy < UpLimit))[0]
	#3. Fish out the right elements
		atomic = [8, 26, 1, 2] # O,Fe,H,He

		for j,a in enumerate(atomic):
			idA = np.where(line[i+2].data['Element']== a )[0]
			mask = np.in1d(idA, acceptID) # Test whether each element of "idA" indices array is also present in a second array "acceptID"
	
			eps[j][index] = erg2keV**(-1)*(line_energy[idA[mask]]*line[i+2].data['Epsilon'][idA[mask]]).sum() #units: phot erg cm^3 s^-1
	return eps

def continua(sim, continuum):

	""" pass class sim file and line fits file, returns epsilon"""
	#1. Read-off temp. cuts
	N = len(continuum)
	eps = np.zeros((4, len(sim.vol)))#initialize arrays for the readout of emissivities per element
	peps = np.zeros((4, len(sim.vol))) #pseudoeps
	#1. Read-off temp. cuts
	cut = continuum[1].data['kT']*eVkB*1e3 #result in K

	for i in range(len(cut)):
		if i==0:
			index = np.where(sim.temp < cut[i])[0]
		elif i==(len(cut)-1):
			index = np.where(sim.temp>cut[i])[0]
		else:
			index = np.where((sim.temp > cut[i-1])& (sim.temp < cut[i]))[0]

	#The position of the element of an array "cut" defines which subtable of line the data belongs to (except line[0]: empty, line[1]: initial parameter block)
	#2. Consider the lines with energies in the band (LowLimit, UpLimit)	
		Econt = continuum[i+2].data['E_cont']
		Cont = continuum[i+2].data['Continuum']
		acceptID = np.where((Econt>LowLimit)&(Econt<UpLimit))[0]
		Epseudo = continuum[i+2].data['E_Pseudo']
		Pseudo = continuum[i+2].data['Pseudo']
		pacceptID = np.where((Epseudo>LowLimit)&(Epseudo<UpLimit))[0] #keV

	#3. Fish out the right elements
		atomic = [8, 26, 1, 2] # O,Fe,H,He

		for j,a in enumerate(atomic):
			idA = np.where(continuum[i+2].data['Z']== a )[0]
			mask = np.in1d(idA, acceptID) # Test whether each element of "idA" indices array is also present in a second array "acceptID"
			pmask = np.in1d(idA, pacceptID)

#integration using the trapezoidal rule Lambda = 0.5*delta_Continuum*delta_Energy*mean_Energy_betweenPoints

			E_mean = np.zeros(len(idA[mask]))
			Ps_mean = np.zeros(len(idA[pmask]))
			suma = 0 #lambda continuum 
			Psuma = 0 #lambda pseudocont
			for k in range(len(idA[mask])-1):
				E_mean[k] = 0.5*(Econt[idA[mask]][k]+Econt[idA[mask]][k+1])
				Ps_mean[k] = 0.5*(Epseudo[idA[pmask]][k]+Epseudo[idA[pmask]][k+1])

				suma = suma + 0.5*E_mean[k]*(Cont[idA[mask]][k+1]+Cont[idA[mask]][k])*(Econt[idA[mask]][k+1]-Econt[idA[mask]][k])

				Psuma = Psuma + 0.5*E_mean[k]*(Pseudo[idA[mask]][k+1]+Pseudo[idA[mask]][k])*(Epseudo[idA[mask]][k+1]-Epseudo[idA[mask]][k])

	
			eps[j][index] = erg2keV**(-1)*suma #units: phot erg cm^3 s^-1
			peps[j][index] = erg2keV**(-1)*Psuma #units: phot erg cm^3 s^-1

	return eps+peps

