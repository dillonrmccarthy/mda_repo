#!/usr/bin/microconda/envs/mda/python
#calculate dihedral angles given four atoms, or given 4 atom selections (Not Implemented) of a given nucleic acid chain. 
#==========================================
#imports

import MDAnalysis as mda
import warnings
warnings.filterwarnings('ignore',category=DeprecationWarning) #need to do AFTER mda import because mda will reset warnings
from MDAnalysis.lib.distances import calc_dihedrals #function to calculate dihedral 
from MDAnalysis.analysis.dihedrals import Dihedral #general dihedrals
import numpy as np

#==========================================
#Classes 

class NAChainDihedral:

    def __init__(self,chain,universe): #ch = chain name you want to process... 
        self._resnums = list(set(universe.select_atoms(f"chainID {chain}").resnums))[1:-1] #strip the first and last residue!
        self._resnames = [list(set(universe.select_atoms(f"chainID {chain} and resnum {_}").resnames))[0] for _ in self._resnums] #need the residue names (AUCG) so I can select the right atoms of the base
        self.chain = chain
        self._universe = universe

    def makeAG(self,**kwargs): #should break this into two different functions. one that sets the kwargs and types, one that does the actual creation 
        """
        given the chain, creates a list of atom groups, where each atom group has 4 unique atoms in it
        """
        #SETTING DEFAULT ASL IF NONE IS GIVEN
        if 'mod_resnum' in kwargs.keys():
            mod_exists = True
            mod_resnum = kwargs['mod_resnum']
            if 'mod_type' in kwargs.keys(): #only check if mod_resnum is also given
                mod_type = kwargs['mod_type']
        else: mod_exists = False

        #create standard asl 
        if 'standard_asl' not in kwargs.keys():
            kwargs['standard_asl'] = ["P","C4'","N#","O3'"]
        #create the thermally destabalizing asl if applicable
        if 'mod_asl' not in kwargs.keys() and mod_exists: #if mod is given create asl if mod_asl was not defined
            #if the type is bbb its special, otherwise it should match the standard selection
            kwargs['mod_asl'] = ["P","C2'","N#","O2'"] if mod_type == 'bbb' else kwargs['standard_asl']

        #------------------------------------------------------------
        self.chain_atomgroup = [] #master list
        for resnum, resname in zip(self._resnums,self._resnames): 
            #begin looping through residues
            _list_of_atom_asl = kwargs['standard_asl']
            if mod_exists and resnum == mod_resnum:
                _list_of_atom_asl = kwargs['mod_asl']
            
            #now loop through each atom of the list 
            _res_asl = [] #will have list of asl's
            for _single_atom in _list_of_atom_asl:
                if _single_atom == "N#":
                    _single_atom = "N9" if resname in ["A","G"] else "N1" #n1 for C and U and everything else :-)
                _res_asl.append(str(f"chainID {self.chain} and resnum {resnum} and name {_single_atom}"))
            
            #-------------------------------------------------------------
            #generate the atomgroup containing 4 atoms
            residue_atomgroup = self._universe.select_atoms(_res_asl[0], _res_asl[1], _res_asl[2], _res_asl[3])
            #append to the master group
            self.chain_atomgroup.append(residue_atomgroup)
            #check to make sure you have the correct number of atoms :) 
            assert len(residue_atomgroup) == 4

        return self

    def compute_dihedrals(self,**kwargs):
        self.makeAG(**kwargs) #make the atom groups
        self.mda_dihedral = Dihedral(self.chain_atomgroup).run(start=0,stop=1) #compute the dihedrals. result will be array of shape (frames,dihedrals)
        self.angles = self.mda_dihedral.results.angles
        return self

#==========================================
if __name__ == '__main__':
    dms = sys.argv([1])
    dcd = sys.argv([2])
    universe = mda.Universe(dms,dcd)

    # NOTE: the use of 'N#' means pick the connective nitrogen (sugar--base). This value changes for A,U,C,G so this is easier 
    # Chain A
    # no atom selection given, use default values
    x = NAChainDihedral("A", universe).compute_dihedrals()
    
    #give custom ASL
    my_asl = ["OP2","C4'","N#","O3'"]
    x = NAChainDihedral("A", universe).compute_dihedrals(standard_asl=my_asl)
    print(x.angles)

    #give custom ASL, but its actually the default so it should match the first one
    standard_asl = ["P","C4'","N#","O3'"]
    x = NAChainDihedral("A", universe).compute_dihedrals(standard_asl=standard_asl)
    print(x.angles)


    #==========================================================
    #Chain B
    
    standard_asl = ["P","C4'","N#","O3'"]
    mod_asl = ["P","C2'","N#","O2'"] 
    #y = NAChainDihedral("B", universe).compute_dihedrals(mod_resnum=7,mod_type='bbb') #not as good
    y = NAChainDihedral("B", universe).compute_dihedrals(mod_resnum=7,standard_asl=standard_asl,mod_asl=mod_asl)
    print(y.angles)
    
    exit()
