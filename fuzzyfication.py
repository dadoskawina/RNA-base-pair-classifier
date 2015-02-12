# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Fuzzyfication functions """

import math
import os.path
import random
from Bio.PDB.PDBParser import PDBParser

class Fuzzyfication:

    def move(self,coord):
        """ Returns fuzzy atom coordinates.
        Ccoordinates given in input are moved in random direction 
        and random distance (max 0,5A)
    
        :param coord - tuple with oryginal coordinates (coord[0],coord[1],coord[2])
    
        :return - new coordinates (x,y,z) """

        # Movement vector in spherical coordinates
        r =  random.uniform(-0.5, 0.5)
        fi = random.uniform(0, 2*math.pi)
        psi = random.uniform(-math.pi*0.5, math.pi*0.5)
        # Movement vector: transformation to cartesian coordinates
        x1 = r*math.cos(fi)*math.cos(psi)
        y1 = r*math.cos(fi)*math.sin(psi)
        z1 = r*math.sin(fi)
        # Sum of base coords and movement
        x = float(coord[0])+x1
        y = float(coord[1])+y1
        z = float(coord[2])+z1
        # Return new coordinates
        return (x,y,z)

    def move_all_atoms(self,atoms1,atoms2):
        """ Returns fuzzy all atoms coordinates
        
        :param atoms1 - a tuple of tuples of (x,y,z) first residue atom coordinates
        :param atoms2 - a tuple of tuples of (x,y,z) second residue atom coordinates
        
        :return - fuzzy coordinates
        
        SimBase3 model takes the following atoms:
        for A and G residues: N9, C2, C6
        for C and U residues: N1, C2, C4

        This method assumes that input atoms are in the above sequence """
        
        moved_atoms1 = []
        for i in range(3):
            moved_atoms1.append(self.move(atoms1[i]))
        moved_atoms2 = []
        for j in range(3):
            moved_atoms2.append(self.move(atoms2[j]))
        return (moved_atoms1,moved_atoms2)


    def atom_distance(self,atom1,atom2):
        """ Calculate a distance between two given atoms
        
        :param atom1 - a tuple of (x,y,z) first atom coordinates
        :param atom2 - a tuple of (x,y,z) second atom coordinates
        
        :return - distance between atom1 and atom2"""
        from math import sqrt
        
        dist = 0
        for i in range(3):
            dist += (atom1[i]-atom2[i])*(atom1[i]-atom2[i])
        return sqrt(dist)

    def SimBase3_distances(self,atoms1,atoms2):
        """ Calculate distances between each pair of given atoms in SimRNA3 model
        
        :param atoms1 - a tuple of tuples of (x,y,z) first residue atom coordinates
        :param atoms2 - a tuple of tuples of (x,y,z) second residue atom coordinates
        
        :return - a tuple with distances between all atoms pairs.
        
        SimBase3 model takes the following atoms:
        for A and G residues: N9, C2, C6
        for C and U residues: N1, C2, C4

        This method assumes that input atoms are in the above sequence
        
        Return tuple: (N9/N1-N9/N1, N9/N1-C2, N9/N1-C6/C4, C2-N9/N1, ...) """
        
        dists = []
        for i in range(3):
            for j in range(3):
                dists.append(self.atom_distance(atoms1[i],atoms2[j]))
        return dists


    def get_data_for_bp (self,nucleotides,path_PDB,mode="train"):
        """ Reads SimBase3 data from PDB files for one nucleotides pair  
        
        :param nucleotides -  eg. AA, AC, AG, AU, CA, CC, ...
        :param path_PDB - path to directory with PDB data
        :param mode - 'train' for training; 'validate' for validation
       
        :return list [[pairtype, coords-atom1, ..., coords-atom9],...] """
        
        result = []
        path = os.path.join(path_PDB,nucleotides)
        for pairtype in ['cWW', 'tWW', 'cWH', 'cHW', 'tWH', 'tHW',
                         'cWS', 'cSW', 'tWS', 'tSW', 'cHH', 'tHH',
                         'cHS', 'cSH', 'tHS', 'tSH', 'cSS', 'tSS']:
            PDB_filename = os.path.join(path,pairtype+".pdb")
        
            parser=PDBParser(QUIET=True)
            try:
                structure=parser.get_structure(pairtype,PDB_filename)   # PDB structure
            except(IOError):
                print("Open file %s error" % (PDB_filename))
            else: # structure is correctly loaded
                #print(nucleotides+" "+pairtype)
                for s_model in structure:
                    # residues info from PDB
                    s1_chain = s_model["A"]
                    s2_chain = s_model["B"]
                    try:
                        if mode=='train':
                            s1_residue = s1_chain[0]
                            s2_residue = s2_chain[1]
                        else:
                            rt = [r for r in s1_chain]
                            s1_residue = rt[1]
                            rt = [r for r in s2_chain]
                            s2_residue = rt[1]
                    except(IndexError):
                        pass
                    else:
                        # first residue 
                        if s1_residue.resname.strip() in ("U", "C"):
                            try:
                                atoms1 = (s1_residue["N1"].get_coord(),
                                          s1_residue["C2"].get_coord(),
                                          s1_residue["C4"].get_coord()
                                         )
                            except(KeyError):
                                break
                        else: # A or G
                            try:
                                atoms1 = (s1_residue["N9"].get_coord(),
                                          s1_residue["C2"].get_coord(),
                                          s1_residue["C6"].get_coord()
                                          )
                            except(KeyError):
                                break
                            
                        # second residue
                        if s2_residue.resname.strip() in ("U", "C"):
                            try:
                                atoms2 = (s2_residue["N1"].get_coord(),
                                          s2_residue["C2"].get_coord(),
                                          s2_residue["C4"].get_coord()
                                         )
                            except(KeyError):
                                break
                        else:  # A or G
                            try:
                                atoms2 = (s2_residue["N9"].get_coord(),
                                          s2_residue["C2"].get_coord(),
                                          s2_residue["C6"].get_coord()
                                         )
                            except(KeyError):
                                break
    
                        # To result
                        result.append([pairtype,atoms1,atoms2])
                            
        return result


    
        
if __name__ == '__main__':
    path_PDB = "..//..//data//ClaRNA//training_pdb"
    for nucleotides in ['AA','AC','AG','AU','CA','CC','CG','CU',
                        'GA','GC','GG','GU','UA','UC','UG','UU',]:
        print(nucleotides)
        print(Fuzzyfication().get_data_for_bp (nucleotides,path_PDB))
                