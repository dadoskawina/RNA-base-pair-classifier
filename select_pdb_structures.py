# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Select needed structures from PDB of-line database """

import gzip
import os.path

def select_pdb_structures (filelist,pathPDB,path_res):
    """ Select needed structures from PDB of-line database
        
        :param list: path to file with ids of needed structures
        :param pathPDB: path to stored (off-line) complete PDB database (with compressed files)
        :param path_res: destination path
    """

    id_list = os.listdir(pathPDB)
      
    f = open(filelist,"r")
    for line in f:
        Id = line[:4]
        print(Id)
        for structure in id_list:
            if Id.lower() in structure:
                print("____struktura")
                filename_source = os.path.join(pathPDB, structure)
                filename_dest = os.path.join(path_res, Id+".pdb")
                #copyfile(filename_source, filename_dest)
                gz = gzip.open(filename_source, 'rb')
                with open(filename_dest, 'wb') as out:
                    out.writelines(gz)
                gz.close()
                


if __name__ == '__main__':
    filelist = "d:\\Aptana-python\\ImpeRNA\\data\\RNA_3188_IDs.txt"
    pathPDB = "e:\\PDB\\PDB_2014-10-26_(104371)-zipped-21GB"
    path_res = "d:\\Aptana-python\\ImpeRNA\\data\\pdb"
    select_pdb_structures(filelist,pathPDB,path_res)
    print('done')