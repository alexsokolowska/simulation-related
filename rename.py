""" 
rename.py
=========

This program renames simulation files. If you output every Xth run of your simulation, i.e. the snapshot numbers follow a sequence with a step X (start, start+X, start+2X, ...), with this program you can change the ordering to be an array of integers with a step X=1. 

"""

import os, glob, numpy as np
from sys import argv
try:
	script, sim_name, start, end, new_start, step = argv
except ValueError:
	print "\n \n USAGE: script sim_name_without_dot start_file_number end_file_number new_start_file_nuber step \n \n "
	raise
	
start = int(start)
end = int(end)
new_start = int(new_start)
step = int(step)

N = (end - start)/step + 1 

print "Number of outputs %d" %N

vec = np.arange(new_start, new_start + N)

for i in range(N):
	if (start+N*i<10):
		string = "0000"+str(start+i*step)
	elif (start+N*i<100):
		string = "000"+str(start+i*step)
	else:
		string = "00"+str(start+i*step)

	path = sim_name + '.' + string + '*'
	LoR = glob.glob(path) #list of runs with the same old number
	print LoR
	
	if vec[i]<10:
		string2 = "0000"+str(vec[i])
	elif vec[i]<100:
		string2 = "000"+str(vec[i])
	else:
		string2 = "00"+str(vec[i])
	for j in range(len(LoR)):
		old_name = LoR[j]
		new_name = sim_name + '.' + string2+LoR[j][len(sim_name)+6:]
		os.rename(old_name, new_name)
		print "Moved ", old_name, "into ", new_name

	
