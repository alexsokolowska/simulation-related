
""" loading snapshots, halos, centering; calculates rvir and rho_crit"""

import pynbody
import math, numpy as np

def open_and_center(f):
	"""open_and_center(f) loads and centers sim on the first halo; f: path to the output"""
	s = pynbody.load(f)
	h = s.halos()
	h[1].physical_units()
	pynbody.analysis.halo.center(h[1], mode='ssc')
	return s
	
def open_and_center2(f):
	"""open_and_center(f) loads and centers sim on the first halo; f: path to the output"""
	s = pynbody.load(f)
	h = s.halos()
	h[1].physical_units()
	pynbody.analysis.halo.center(h[1], mode='ssc')
	return s,h


def rvir(s, delta_c=None):
	"""rho_vir(s) calculates the critical density; s: simulation snapshot; apply when centered on a given halo; optional keywords: delta_c (float)"""
	#loading cosmological parameters
	aexp        =  s.properties['a']
	H0          =  s.properties['h']*100 #km/s/Mpc
	omega_m     =  s.properties['omegaM0']
	omega_l     =  s.properties['omegaL0']
#calculations of rho_vir
	z = (1/aexp)-1
	omega_m_z=(omega_m*(1+z)**3)/(1-omega_m+omega_m*(1+z)**3)
	if delta_c is None:
		delta_c=49+96*omega_m_z+(200*omega_m_z)/(1+5*omega_m_z)
	print 'Control: delta_c \n ', delta_c
	H=H0*np.sqrt(omega_m*(1+z)**3+omega_l)
	H = pynbody.units.Unit("km s^-1 Mpc^-1")*H
	rho_crit=3*H**2/(8*math.pi*pynbody.units.G)
	rho_vir=rho_crit*delta_c
	rho_vir_Msol=rho_vir.in_units("Msol kpc^-3")
	p = pynbody.analysis.profile.Profile(s,type='log',min=2.,max=500,ndim=3)
	rho_tot=p['mass_enc'].in_units("Msol")[::-1][:-1]/((4./3*math.pi)*(p['rbins'].in_units("kpc")[::-1][:-1]**3))
	return np.interp(rho_vir_Msol, rho_tot, p['rbins'].in_units("kpc")[::-1][:-1])

def rho_crit(s, delta_c=None):
	#loading cosmological parameters
	aexp        =  s.properties['a']
	H0          =  s.properties['h']*100 #km/s/Mpc
	omega_m     =  s.properties['omegaM0']
	omega_l     =  s.properties['omegaL0']
#calculations of rho_vir
	z = (1/aexp)-1
	omega_m_z=(omega_m*(1+z)**3)/(1-omega_m+omega_m*(1+z)**3)
	if delta_c is None:
		delta_c=49+96*omega_m_z+(200*omega_m_z)/(1+5*omega_m_z)
	print 'Control: delta_c \n ', delta_c
	H=H0*np.sqrt(omega_m*(1+z)**3+omega_l)
	H = pynbody.units.Unit("km s^-1 Mpc^-1")*H
	rho_crit=3*H**2/(8*math.pi*pynbody.units.G)
	return rho_crit
	
def new_arrays(s):
        s.g['N_HI'] = s.g['HI']*s.g['mass'].in_units("m_p")
        s.g['N_HeI'] = s.g['HeI']*s.g['mass'].in_units("m_p")/4.
        s.g['N_ion_HII'] = s.g['mass'].in_units("m_p")*s.g['HII']               
        s.g['N_ion_HeII'] = s.g['mass'].in_units("m_p")*s.g['HeII']/4. 
        s.g['N_ion_HeIII'] = s.g['mass'].in_units("m_p")*s.g['HeIII']/4.
        s.g['N_ions'] = s.g['N_ion_HII'] + s.g['N_ion_HeII'] + s.g['N_ion_HeIII']
        s.g['N_electrons'] = s.g['N_ion_HII'] + s.g['N_ion_HeII'] + 2*s.g['N_ion_HeIII'] # in units m_p
        s.g['N_cool'] = s.g['N_ions'] + s.g['N_electrons']
        s.g['all_particles']= s.g['N_HI'] +s.g['N_HeI'] +s.g['N_ions']+s.g['N_electrons']
        s.g['ia'] = s.g['N_HI'] +s.g['N_HeI'] +s.g['N_ions']

def number_of_particles(s):
        """creates new arrays for a snapshot s"""
        s.g['N_ion_HII'] = s.g['mass'].in_units("m_p")*s.g['HII']
        s.g['N_ion_HeII'] = s.g['mass'].in_units("m_p")*s.g['HeII']/4. 
        s.g['N_ion_HeIII'] = s.g['mass'].in_units("m_p")*s.g['HeIII']/4.
        s.g['N_electrons'] = s.g['N_ion_HII'] + s.g['N_ion_HeII'] + 2*s.g['N_ion_HeIII'] # in units m_p
        s.g['N_particles'] = s.g['N_ion_HII'] + s.g['N_ion_HeII'] + s.g['N_ion_HeIII']+s.g['N_electrons']
        return s.g['N_particles']	
