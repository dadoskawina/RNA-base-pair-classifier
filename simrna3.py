# This is a part of ImpeRNA package.
# @author: Jacek Smietanski
# http://jaceksmietanski.net

""" Compact PDB files to SimRNA3 model """

import os.path
from os import listdir

def PDB_to_simRNA3(path_PDB,path_res):
    """ Compact PDB files to SimRNA3 model

    Skrypt usuwa z plikow PDB część nagłówkową oraz zbędne łańcuchy, nie
    reprezentujące RNA. W przypadku struktur wielomodelowych, pozostawia tylko
    pierwszy model.

    Zostawiamy tylko atomy z modelu SimRNA3 (bez heteroatomów)

    :param path_PDB = folder ze źródłową bazą PDB
    :param path_res = folder, do którego będą zapisywane zmodyfikowane dane.
         Do oryginalnych nazw plików dodawane będzie rozszerzenie '.sim3' """

    for fname in listdir(path_PDB):
        PDBid = fname[:4].upper()
        print(fname)
        filename_PDB = os.path.join(path_PDB, fname)
        filename_res = os.path.join(path_res, fname+'.sim3')

        fr = open(filename_res,'w')
        f = open(filename_PDB,"r")

        model = False  #mówi czy jestesmy wewnatrz pierwszego modelu
        skip = False   #mówi czy mamy pominąć analizę dalszej części pliku

        for line in f:
            recordName = line[:6]
         
            if recordName=="MODEL ":
                if model==True:
                    skip = True #to drugie wystąpienie słowa model; koniec analizowania struktury
                else:
                    model = True
            elif recordName=="ENDMDL":
                skip = True

            if not skip:
                if recordName in ["ATOM  ", "HETATM", "TER   "]:
                    chainId = line[21]
                    if recordName=="ATOM  ":
                        resName = line[17:20].strip()
                        atomName = line[12:16].strip()
                        if resName in ['U','C']:
                            if atomName in ["P","C4'","N1","C2","C4"]:
                                fr.write(line)
                        elif resName in ['A','G']:
                            if atomName in ["P","C4'","N9","C2","C6"]:
                                fr.write(line)
        f.close()
        fr.close()
 

if __name__ == '__main__':
    path_PDB = "..//..//data//PDB"
    path_res = "..//..//data//PDB_SimRNA3"
    PDB_to_simRNA3(path_PDB,path_res)


