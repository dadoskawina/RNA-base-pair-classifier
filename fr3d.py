# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Parse FR3D files """

import os.path

class FR3D:
    """ Parse FR3D files """

    def __init__(self):
        """ Initialize
        :param pairs - list of pairs in analyzed structure """
        self.pairs = []

    def getPair(self,line):
        """ Read one row from FR3D file """
        first = line[6]
        if line[12]=='(':  # conditions for manage some exceptions
            residue1 = line[7:12].strip()
            chain1 = line[13]
            second = line[18]
            if line[24]=='(':
                typ = line[33:36]
                residue2 = line[19:24].strip()
                chain2 = line[25]
            else:
                typ = line[32:35]
                residue2 = line[19:23].strip()
                chain2 = line[24]
        elif line[23]=='(':
            second = line[17]
            typ = line[32:35]
            residue1 = line[7:11].strip()
            chain1 = line[12]
            residue2 = line[18:23].strip()
            chain2 = line[24]
        else:
            second = line[17]
            typ = line[31:34]
            residue1 = line[7:11].strip()
            chain1 = line[12]
            residue2 = line[18:22].strip()
            chain2 = line[23]
        pair = first + second 

        try:
            resSeq1 = int(residue1)
            iCode1= " "
        except:
            resSeq1 = int(residue1[0:-1]) 
            iCode1 = residue1[-1]

        try:
            resSeq2 = int(residue2)
            iCode2= " "
        except:
            resSeq2 = int(residue2[0:-1]) 
            iCode2 = residue2[-1]

        return({'chain1' : chain1, 
                'chain2' : chain2,
                'resseq1' : resSeq1,
                'resseq2' : resSeq2,
                'icode1' : iCode1,
                'icode2' : iCode2,
                'first' : first,
                'second' : second,
                'basepair' : pair,
                'type' : typ})


    def parseFR3D (self,path):
        """ Analyze the whole FR3D file
        
        :param path: input fileNAME """

        self.pairs = []
        if os.path.isfile(path):
            try:
                f = open(path,"r")
            except IOError:
                print("Open file %s error (%s)" % path)
            else:
                for line in f:
                    self.pairs.append(self.getPair(line))
                f.close()
                
    def get_basepairs(self,path):
        self.parseFR3D(path)
        return self.pairs
         

if __name__ == '__main__':
    # test
    path = "..//..//data//fr3d//157D.FR3D"
    p = FR3D().get_basepairs(path)
    for i in range(len(p)):
        print(p[i])
