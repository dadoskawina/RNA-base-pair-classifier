""" Zliczanie danych testowych dotyczacych par RNA, zrodlo FR3D """

""" Output:
 in each row columns:
 [0:2] pairing type
 [3]   first base name
 [4]   second base name
 [6:]  counted number of basepairs

 Note that according to FR3D file structure, each basepair is counted twice.
"""

__author__ = "Jacek Smietanski"
__email__ = "jacek.smietanski@ii.uj.edu.pl" 

import os.path
path = "..\\data\\fr3d\\"
counter = {}  # dictionary, with counts each base in all FR3D files
for f in os.listdir(path): 
   filename = os.path.join(path, f)
   if os.path.isfile(filename):
       try:
           file = open(filename,"r")
       except IOError, error:
           print "Open file %s error (%s)" % (filename,error)
       else:
           for i,line in enumerate(file.readlines()):
                first = line[6]
                second = line[17]
                typ = line[31:34]
                key = typ + first + second
                try:
                   counter[key]+=1
                except KeyError:
                   counter[key]=1
           file.close()
    
""" Summary """
for key in sorted(counter.keys()):
    print "%s %d" % (key,counter[key])
    

