# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Download PDB structures via ftp """

import gzip
import os

from Bio.PDB.PDBList import PDBList

#Importing this function with leading underscore as not intended for reuse
from Bio._py3k import urlretrieve as _urlretrieve

class PDBListExtension(PDBList):
    """ Functionality extension of PDBList class from biopython package """ 
    
    def retrieve_pdb_biounit_file(self, pdb_code, pdir=None):
        """ Retrieves a PDB biounit structure file from the PDB server. 
        
        If more than one biounit file exists, all are downloaded.
        
        This code is based on retrieve_pdb_file method, by Kristian Rother,
        from biopython package (http://biopython.org).
        
        The PDB structure's file name(s) is returned as a single string.

        :param pdir: put the file in this directory (default: create a PDB-style
                     directory tree)
        :type pdir: string

        :return: filename(s)
        :rtype: string
        """

        code = pdb_code.lower()

        # Where does the final PDB file get saved?
        if pdir is None:
            path = self.local_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)

        # If several files for one structure exists, alle of them are downloaded
        loop = True
        count = 1
        final_files = ""
        while loop: 
            archive_fn = "%s.pdb%d.gz" % (code, count) 
            url = (self.pdb_server +
                   '/pub/pdb/data/biounit/coordinates/divided/%s/%s' %
                   (code[1:3], archive_fn))

            filename = os.path.join(path, archive_fn)
            final_file = os.path.join(path, "%s.pdb%d.unit.pdb" % (code,count))  # (decompressed)

            # Skip download if the file already exists
            if not self.overwrite:
                if os.path.exists(final_file):
                    print("Structure exists: '%s' " % final_file)
                    break

            # Get the compressed PDB structure
            # Retrieve the file
            print("Downloading PDB structure '%s'..." % pdb_code)
            try:
                _urlretrieve(url, filename)
            except IOError:
                print("No such structure")
                loop = False
            else:
                print("OK")
            # Uncompress the archive, delete when done
            #Can't use context manager with gzip.open until Python 2.7
                gz = gzip.open(filename, 'rb')
                with open(final_file, 'wb') as out:
                    out.writelines(gz)
                gz.close()
                os.remove(filename)
                final_files += final_file
                count += 1

        return final_files
    
    
# Execution for all RNA structures
if __name__ == '__main__':
    retriveRNAstructures = PDBListExtension()

   
    RNAlist = open("D:\Aptana-Python\ImpeRNA\data\RNA_3188_IDs.txt","r")
    destination_dir = "D:\\Aptana-Python\\ImpeRNA\\data\\biounits"

    #retriveRNAstructures.retrieve_pdb_biounit_file("7msf",destination_dir)
    
    for line in RNAlist:
        pdb_code = line.strip()
        retriveRNAstructures.retrieve_pdb_biounit_file(pdb_code,destination_dir)