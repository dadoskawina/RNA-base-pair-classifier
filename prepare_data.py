# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Reads coordinates from PDB files and creates files with distances """

import os.path
from Bio.PDB.PDBParser import PDBParser
from src.data.fr3d import FR3D 
class PrepareData:
    """ Reads coordinates from PDB files and creates files with distances """

    
    def square_atom_distance(self,atom1,atom2):
        """ Calculate a square of distance between two given atoms
        
        :param atom1 - a tuple of (x,y,z) first atom coordinates
        :param atom2 - a tuple of (x,y,z) second atom coordinates
        
        :return - square of distance between atom1 and atom2"""
        
        sqdist = 0
        for i in range(3):
            sqdist += (atom1[i]-atom2[i])*(atom1[i]-atom2[i])
        return sqdist

    
    def SimRNA3_square_distances(self,atoms1,atoms2):
        """ Calculate square of distances between each pair of given atoms
        in SimRNA3 model
        
        :param atoms1 - a tuple of tuples of (x,y,z) first residue atom coordinates
        :param atoms2 - a tuple of tuples of (x,y,z) second residue atom coordinates
        
        :return - a tuple with square of distances between all atoms pairs.
        
        SimRNA3 model takes the following atoms:
        for A and G residues: P, C4’, N9, C2, C6
        for C and U residues: P, C4’, N1, C2, C4

        This method assumes that input atoms are in the above sequence
        
        Return tuple: (P-P, P-C4', P-N9/N1, P-C2, P-C6/C4, C4'-P, ...) """
        
        sqdists = []
        for i in range(5):
            for j in range(5):
                sqdists.append(self.square_atom_distance(atoms1[i],atoms2[j]))
        return sqdists


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

    
    def SimRNA3_distances(self,atoms1,atoms2):
        """ Calculate distances between each pair of given atoms in SimRNA3 model
        
        :param atoms1 - a tuple of tuples of (x,y,z) first residue atom coordinates
        :param atoms2 - a tuple of tuples of (x,y,z) second residue atom coordinates
        
        :return - a tuple with distances between all atoms pairs.
        
        SimRNA3 model takes the following atoms:
        for A and G residues: P, C4’, N9, C2, C6
        for C and U residues: P, C4’, N1, C2, C4

        This method assumes that input atoms are in the above sequence
        
        Return tuple: (P-P, P-C4', P-N9/N1, P-C2, P-C6/C4, C4'-P, ...) """
        
        dists = []
        for i in range(5):
            for j in range(5):
                dists.append(self.atom_distance(atoms1[i],atoms2[j]))
        return dists


    
    def create_training_file(self,nucleotides,pairtype,structlist,path_FR3D,path_PDB,path_res,mode="pdb.sim3"):
        """ Creates a distance file for one given nucleodides- and pair type 
        
        :param nucleotides -  eg. AA, AC, AG, AU, CA, CC, ...
        :param pairtype - eg. cWW, cWS, ...
        :param structlist - list of PDBid-s of structures to analyze
        :param path_FR3D - path to directory with FR3D data
        :param path_PDB - path to directory with PDB data
        :param path_res - path to directory where results should be saved
        :param mode - extension of pdb files (default: pdb.sim3)  """
        
        # TODO - co z plikami FR3D, gdzie jest więcej niż jeden model?
        
        result_filename =os.path.join(path_res,nucleotides+'.'+pairtype+'.data')
        result_file = open(result_filename,"w")
        
        for PDBid in structlist:
            FR3D_filename = os.path.join(path_FR3D,PDBid+'.fr3d')
            basepairs = FR3D().get_basepairs(FR3D_filename)
            
            PDB_filename = os.path.join(path_PDB,PDBid+'.'+mode)
            parser=PDBParser(QUIET=True)
            try:
                structure=parser.get_structure(PDBid,PDB_filename)   # PDB structure
            except(IOError):
                print("Open file %s error" % (PDB_filename))
                return
            else: # structure is correctly loaded
                print(PDBid)
                for basepair in basepairs:
                    if basepair['basepair']==nucleotides and basepair['type']==pairtype:
                        s_model = structure[0]        #always first model
                        # first residue info from PDB
                        s1_chain = s_model[basepair['chain1']]
                        s1_residue = s1_chain[basepair['resseq1']]
                        # second residue info from PDB
                        s2_chain = s_model[basepair['chain2']]
                        s2_residue = s2_chain[basepair['resseq2']]

                        # first residue 
                        if s1_residue.resname.strip() in ("U", "C"):
                            try:
                                atoms1 = (s1_residue["P"].get_coord(),
                                          s1_residue["C4'"].get_coord(),
                                          s1_residue["N1"].get_coord(),
                                          s1_residue["C2"].get_coord(),
                                          s1_residue["C4"].get_coord()
                                         )
                            except(KeyError):
                                break
                        else: # A or G
                            try:
                                atoms1 = (s1_residue["P"].get_coord(),
                                          s1_residue["C4'"].get_coord(),
                                          s1_residue["N9"].get_coord(),
                                          s1_residue["C2"].get_coord(),
                                          s1_residue["C6"].get_coord()
                                         )
                            except(KeyError):
                                break
                        
                        # second residue
                        if s2_residue.resname.strip() in ("U", "C"):
                            try:
                                atoms2 = (s2_residue["P"].get_coord(),
                                          s2_residue["C4'"].get_coord(),
                                          s2_residue["N1"].get_coord(),
                                          s2_residue["C2"].get_coord(),
                                          s2_residue["C4"].get_coord()
                                         )
                            except(KeyError):
                                break
                        else: # A or G
                            try:
                                atoms2 = (s2_residue["P"].get_coord(),
                                          s2_residue["C4'"].get_coord(),
                                          s2_residue["N9"].get_coord(),
                                          s2_residue["C2"].get_coord(),
                                          s2_residue["C6"].get_coord()
                                         )
                            except(KeyError):
                                break

                        # Distance calculation
                        distances = self.SimRNA3_distances(atoms1,atoms2)
                        line = [PDBid,s1_chain.id,s1_residue.id[1],
                                s2_chain.id,s2_residue.id[1] ] + distances
                        for element in line:
                            result_file.write(str(element)+"\t")
                        result_file.write("\n")
                        
        result_file.close()



    def create_all_training_files(self,structlist,path_FR3D,path_PDB,path_res,mode="pdb.sim3"):
        """ Creates a distance file for all possible nucleodides- and pair types. 
        
        :param structlist - list of PDBid-s of structures to analyze
        :param path_FR3D - path to directory with FR3D data
        :param path_PDB - path to directory with PDB data
        :param path_res - path to directory where results should be saved
        :param mode - extension of pdb files (default: pdb.sim3)  """

        for nucleotides in ['AA','AC','AG','AU','CA','CC','CG','CU',
                            'GA','GC','GG','GU','UA','UC','UG','UU',]:
            for pairtype in ['cWW', 'cWw', 'cwW', 'tWW', 'tWw', 'twW',
                             'cWH', 'cHW', 'tWH', 'tHW',
                             'cWS', 'cSW', 'tWS', 'tSW',
                             'cHH', 'cHh', 'chH', 'tHH', 'tHh', 'thH',
                             'cHS', 'cSH', 'tHS', 'tSH',
                             'cSS', 'cSs', 'CsS', 'tSS', 'tSs', 'tsS', 'bif']:
                self.create_training_file(nucleotides,pairtype,structlist,path_FR3D,path_PDB,path_res,mode)

       

    def create_list_assymetric(self,path_FR3D,path_PDB,path_biounit):
        """ Creates a list of PDBids for assymetric structures
        
        :param path_FR3D - path to directory with FR3D data
        :param path_PDB - path to directory with PDB assymetric data
        :param path_biounit - path to directory with PDB biounit data
        
        :return - list of PDBids 
        
        For each FR3D non-empty data (FR3D identifies at least one pair) search
        for structures, where there is no biounit data (it means that the
        assymetric unit is also a biological assemble). """

        from os import listdir
        # ids of structures for wchich at least one pair is identified
        not_empty_ids = [ f[:4] for f in listdir(path_FR3D) 
                         if os.path.isfile(os.path.join(path_FR3D,f)) 
                         and os.path.getsize(os.path.join(path_FR3D,f))>0 ]
        # ids of structures for wchich the bio unit exists
        bio_unit_ids = [ f[:4].upper() for f in listdir(path_biounit) 
                         if os.path.isfile(os.path.join(path_biounit,f)) ]
        # create result
        return [ name for name in not_empty_ids if name not in bio_unit_ids ]



    def create_list_1biounit(self,path_FR3D,path_PDB,path_biounit):
        """ Creates a list of PDBids for structures where exactly one biounit file is given
        
        :param path_FR3D - path to directory with FR3D data
        :param path_PDB - path to directory with PDB assymetric data
        :param path_biounit - path to directory with PDB biounit data
        
        :return - list of PDBids 
        
        For each FR3D non-empty data (FR3D identifies at least one pair) seach
        for structures, where there is exactly one biounit file (it means that the
        assymetric unit differs from a biological assemble and ). """

        from os import listdir
        # ids of structures for wchich at least one pair is identified
        not_empty_ids = [ f[:4] for f in listdir(path_FR3D) 
                         if os.path.isfile(os.path.join(path_FR3D,f)) 
                         and os.path.getsize(os.path.join(path_FR3D,f))>0 ]
        # ids of structures for wchich at least one (may be more) bio unit exists
        bio_unit_ids = [ f[:4].upper() for f in listdir(path_biounit) 
                         if os.path.isfile(os.path.join(path_biounit,f)) ]
        # ids of structures for wchich at least two bio units exist
        bio_unit_more_ids = [ f[:4].upper() for f in listdir(path_biounit) 
                         if (os.path.isfile(os.path.join(path_biounit,f)) and f[8]=="2") ]
        # create result
        return [ name for name in not_empty_ids if (name in bio_unit_ids) and (name not in bio_unit_more_ids) ]



    def create_list_morebiounits(self,path_FR3D,path_PDB,path_biounit):
        """ Creates a list of PDBids for structures where more than one biounit file exists
        
        :param path_FR3D - path to directory with FR3D data
        :param path_PDB - path to directory with PDB assymetric data
        :param path_biounit - path to directory with PDB biounit data
        
        :return - list of PDBids 
        
        For each FR3D non-empty data (FR3D identifies at least one pair) search
        for structures, where there are at least two biounit files """

        from os import listdir
        # ids of structures for wchich at least one pair is identified
        not_empty_ids = [ f[:4] for f in listdir(path_FR3D) 
                         if os.path.isfile(os.path.join(path_FR3D,f)) 
                         and os.path.getsize(os.path.join(path_FR3D,f))>0 ]
        # ids of structures for wchich at least two bio units exist
        bio_unit_more_ids = [ f[:4].upper() for f in listdir(path_biounit) 
                         if (os.path.isfile(os.path.join(path_biounit,f)) and f[8]=="2") ]
        # create result
        return [ name for name in not_empty_ids if name in bio_unit_more_ids ]


    def select_FR3D_PDB_to_analyze(self,struct_list,source_paths,result_paths):
        """ Combines PDB and Biounit data to select good data source.
        
        :param struct_list - list of structures considered to analyze (eg. nonredundant only)
        :param source_paths - path to source folders
        :param result_paths - path to result folders
        ... """
        from shutil import copyfile
        
        source_path_FR3D = os.path.join(source_paths,"FR3D")
        source_path_PDB = os.path.join(source_paths,"PDB")
        source_path_biounit = os.path.join(source_paths,"biounits")

        result_path_FR3D = os.path.join(result_paths,"FR3D")
        result_path_PDB = os.path.join(result_paths,"PDB")
        
        assym = self.create_list_assymetric(source_path_FR3D,source_path_PDB,source_path_biounit)
        biounit1 = self.create_list_1biounit(source_path_FR3D,source_path_PDB,source_path_biounit)
        biounit2 = self.create_list_morebiounits(source_path_FR3D,source_path_PDB,source_path_biounit)
        
        for PDBid in struct_list:
            source_FR3D = os.path.join(source_path_FR3D,PDBid+".fr3d")
            dest_FR3D = os.path.join(result_path_FR3D,PDBid+".fr3d")
            dest_PDB = os.path.join(result_path_PDB,PDBid+".pdb")
            if PDBid in assym:
                source_PDB = os.path.join(source_path_PDB,PDBid+".pdb")
                try:
                    copyfile(source_FR3D,dest_FR3D)
                except(FileNotFoundError):
                    pass
                try:
                    copyfile(source_PDB,dest_PDB)
                except(FileNotFoundError):
                    print ("Not found PDB "+PDBid)
            elif PDBid in biounit1:
                source_PDB = os.path.join(source_path_biounit,PDBid.lower()+".pdb1.unit.pdb")
                try:
                    copyfile(source_FR3D,dest_FR3D)
                except(FileNotFoundError):
                    pass
                try:
                    copyfile(source_PDB,dest_PDB)
                except(FileNotFoundError):
                    print ("Not found PDB "+PDBid)
            else: #PDBid in biounit2:
                #TODO
                #print("Pomijamy "+PDBid) #strukturę pomijamy: wyswietlamy co to jest, moze recznie dodamy
                pass
       
      
        
if __name__ == '__main__':
    """
    structfile = "..//..//data//nonredundant-lists//Nonredundant_4A_2014-12-05_list.pdb"
    struct_list = []
    f = open(structfile,"r")
    for line in f:
        struct_list.append(line.strip())
    f.close
    print (len(struct_list))
        
    source_paths = "..//..//data"
    result_paths = "..//..//data//selected" 
    PrepareData().select_FR3D_PDB_to_analyze(struct_list,source_paths,result_paths)  
    
    import src.data.simrna3
    path_PDB = "..//..//data//selected//PDB"
    path_res = "..//..//data//selected//PDB_SimRNA3"
    src.data.simrna3.PDB_to_simRNA3(path_PDB,path_res)
    """

    path_FR3D = "..//..//data//selected//fr3d"
    path_PDB = "..//..//data//selected//PDB_SimRNA3"
    mode = "pdb.sim3"
    path_res = "..//..//data//selected//distances_SimRNA3"
    
    from os import listdir
    struct2_list = [ f[:4] for f in listdir(path_FR3D) 
                         if os.path.isfile(os.path.join(path_FR3D,f)) 
                         and os.path.getsize(os.path.join(path_FR3D,f))>0 ]
    
    PrepareData().create_all_training_files(struct2_list,path_FR3D,path_PDB,path_res,mode)

        