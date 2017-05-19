import loading as ld
import pickle, gc, pynbody
import numpy as np
from input_parameters import pathto, fs, name

keywords = ['vol', 'hydrogen', 'hetot', 'temp', 'OxMassFrac', 'FeMassFrac', 'nee', 'ni']
filenames = ["Vol.p","H.p","He.p", "Temp.p","Ox.p","Fe.p","ne.p","ni.p"]
for i in range(len(fs)):
#for i in range(1,2):
	s = ld.open_and_center(fs[i])
	ld.new_arrays(s)
	rvir = ld.rvir(s)
	sphere = pynbody.filt.Sphere(rvir[i])
	s.g['vol'] = (s.g['mass']/s.g['rho']).in_units("cm^3")
	s.g['ni'] = s.g['N_ions']/s.g['vol'].in_units("cm^3")
	s.g['nee'] =  s.g['N_electrons']/s.g['vol'].in_units("cm^3")

	for j in range(8):
	
	    pickle.dump(np.array(s.g[sphere][keywords[j]]), open(pathto+name[i]+filenames[j], "wb"))

        close(pathto+name[i]+filenames[j])
        
    print "Extracted files from your snapshot."
