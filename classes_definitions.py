from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
import atom_dict
import math

class StructureAnalysis:
    column_names = ['b_factor', 'sasa', 'aliphatic', 'aromatic', 'donor', 'acceptor',
       'don_acc', 'charge', 'CG', 'totalC', 'ND2', 'totalN', 'OE2', 'totalO',
       'CD', 'OG1', 'CG2', 'CB', 'CA', 'N', 'C', 'CD2', 'CE2', 'CD1', 'OE1',
       'O', 'CE', 'SD', 'totalS', 'environment_density', 'CZ', 'NZ', 'OH',
       'CE1', 'CG1', 'OD1', 'NE2', 'NH1', 'OD2', 'OG', 'NH2', 'NE', 'ND1',
       'SG', 'NE1', 'CH2', 'CZ2', 'CE3', 'CZ3']
    
        
    def __init__(self, pdb_file):
        
        parser = PDBParser()
        self.structure = parser.get_structure('structure', pdb_file)
        self.standard_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        self.ns = NeighborSearch(list(self.structure.get_atoms()))
        self.sr = ShrakeRupley()
        self.sasa = self.sr.compute(self.structure, level = 'A')



    def get_ligands_from_structure(self):
        '''Returns list of ligands in the structure'''
        ligand_residues = []

        for residue in self.structure.get_residues():
            resname = residue.get_resname()
            if resname not in self.standard_residues and resname != 'HOH':
                ligand_residues.append(resname)

        return ligand_residues
    


    def get_ligand_binding_site_atoms(self):
        '''Returns list of atoms in the ligand binding site'''
    
        list_of_ligands = self.get_ligands_from_structure() # get ligands from the structure
        ligand_binding_site_atoms = []
    

        # Iterate over the residues in the structure
        for residue in self.structure.get_residues():
            residue_name = residue.get_resname()
            # Check if the residue is a ligand
            if len(list_of_ligands) > 0:    
                if residue_name in list_of_ligands:
                    
                    # Iterate over the atoms in the ligand
                    for atom in residue.get_atoms():
                        # Find all atoms within 8A of the ligand atom
                        nearby = self.ns.search(atom.coord, 8)
                        # Filter atoms that belong to the protein
                        nearby_protein_atoms = [a for a in nearby if a.get_parent().get_resname() in self.standard_residues]
                        # Add nearby protein atoms to the ligand binding site atoms list
                        ligand_binding_site_atoms.extend([(nearby_atom.get_serial_number() , nearby_atom.get_name(), nearby_atom.coord, nearby_atom.get_parent().get_resname()) for nearby_atom in nearby_protein_atoms])

        return ligand_binding_site_atoms


    def get_structure_environments(self):
        '''Samples the atoms of the structure, computing the features of the environment of each atom. 
        Designed to create the training set, sampling as many atoms from non-ligand-binding residues as from ligand-binding residues.'''

        # Get the ligands from the structure
        ligand_names = self.get_ligands_from_structure()
        
        # Check if there are ligands in the structure
        if len(ligand_names) > 0:
            
            # Get the atoms in the ligand binding site
            ligand_binding_atoms = self.get_ligand_binding_site_atoms()
            # Get the serial numbers of the atoms in the ligand binding site
            serial_numbers_lbs = set([t[0] for t in ligand_binding_atoms])

            # create a list of non lbs atoms
            non_lbs_atoms = []

            for residue in self.structure.get_residues():
                for atom in residue.get_atoms():
                    if atom.get_serial_number() not in serial_numbers_lbs:
                        non_lbs_atoms.append(atom)


            # Initialize the dictionary of environments
            environments = {}

            # Define the radius of the sphere
            radio = 6.2
            # Calculate the volume of the sphere
            sphere_volume = 4/3 * math.pi * radio**3




            for residue in self.structure.get_residues():
            

                # iterate over the atoms of the protein 

                if residue.get_resname() not in ligand_names and residue.get_resname() != 'HOH':  # Exclude water residues and ligand residues

                    for selected_atom in residue.get_atoms():


                        nearby_atoms = self.ns.search(selected_atom.coord, radio)  # Change the radius if needed

                        # Filter nearby atoms that are not water and are not in residues of ligand_names
                        nearby_atoms = [atom for atom in nearby_atoms if atom.get_parent().get_resname() != 'HOH' and atom.get_parent().get_resname() not in ligand_names]
                        
                        # Initialize the atom counts and counter outside of the loop
                        atom_counts = {}
                        atom_counter = 0
                        
                        # set column names
                        
                        for i in self.column_names:
                            atom_counts[i] = 0

                        for nearby_atom in nearby_atoms:
                            atom_counter += 1
                            
                            # Count the number of atoms of each type
                            if nearby_atom.get_id() in atom_counts:
                                atom_counts[nearby_atom.get_id()] += 1
                            else:
                                atom_counts[nearby_atom.get_id()] = 1
                            

                            # Count the total number of atoms of each element

                            if str("total"+nearby_atom.element) in atom_counts:
                                atom_counts["total"+nearby_atom.element] += 1
                            else:
                                atom_counts["total"+nearby_atom.element] = 1
                        
                            # add b_factor and sasa
                        
                            atom_counts["b_factor"] += nearby_atom.get_bfactor()

                            atom_counts["sasa"] += nearby_atom.sasa

                            # Calculate the features of the atom
                            if nearby_atom.get_id() in list(atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()].keys()):
                                # This is a counter for each atom is it is aliphatic, aromatic, donor, acceptor or donor_acceptor.
                                atom_counts[atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()][nearby_atom.get_id()]] += 1

                            if nearby_atom.get_id() == "N":
                                atom_counts["donor"] +=1

                            if nearby_atom.get_id() == "O":
                                atom_counts["acceptor"] +=1

                            if nearby_atom.get_id() == "C":

                                atom_counts["aromatic"] +=1
                                
                            resname = nearby_atom.get_parent().get_resname().capitalize()
                            if resname in atom_dict.charges and nearby_atom.get_id() in atom_dict.charges[resname]:
                                atom_counts["charge"] += atom_dict.charges[resname][nearby_atom.get_id()]
                    
                    
                        
                        
                        atom_proportions = {atom_id: count / atom_counter for atom_id, count in atom_counts.items()}
                        
                        # create label of atom in lbs

                        atom_proportions["environment_density"] = atom_counter/sphere_volume

                        if selected_atom.get_serial_number() in serial_numbers_lbs:
                            atom_proportions["is_lbs"] = 1

                        else:
                            atom_proportions["is_lbs"] = 0


                        # use a chimera recognisable code as key 
                            
                        chimera_code = ":" + str(residue.get_id()[1]) + "@" + selected_atom.get_name()
                        environments[chimera_code] = atom_proportions

            return environments
            
        else:
            print("No ligand names provided")
            return None
        








    def get_input_environments(self):
        '''Returns dictionary of environments of each residue in the structure of the input PDB file'''
        environments = {}
        radio = 6.2
        sphere_volume = 4/3 * math.pi * radio**3

        for residue in self.structure.get_residues():
            
            if residue.get_resname() in self.standard_residues:  # Just include standard residues. 
                for selected_atom in residue.get_atoms():
                    nearby_atoms = self.ns.search(selected_atom.coord, radio)  # Change the radius if needed

                    # Filter nearby atoms that are standard residues.
                    nearby_atoms = [atom for atom in nearby_atoms if atom.get_parent().get_resname() in self.standard_residues]
                    
                    # Initialize the atom counts and counter outside of the loop
                    atom_counts = {}
                    atom_counter = 0
                    
                    for i in self.column_names:
                        atom_counts[i] = 0

                    for nearby_atom in nearby_atoms:
                        atom_counter += 1
                        
                        # Count the number of atoms of each type
                        if nearby_atom.get_id() in atom_counts:
                            atom_counts[nearby_atom.get_id()] += 1
                        else:
                            atom_counts[nearby_atom.get_id()] = 1
                        

                        # Count the total number of atoms of each element
                

                        if str("total"+nearby_atom.element) in atom_counts:
                            atom_counts["total"+nearby_atom.element] += 1
                        else:
                            atom_counts["total"+nearby_atom.element] = 1
                    
                        # Count the total number of atoms
                    
                        atom_counts["b_factor"] += nearby_atom.get_bfactor()

                        atom_counts["sasa"] += nearby_atom.sasa

                        # Calculate the features of the atom
                        if nearby_atom.get_id() in list(atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()].keys()):
                                # This is a counter for each atom is it is aliphatic, aromatic, donor, acceptor or donor_acceptor.
                                atom_counts[atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()][nearby_atom.get_id()]] += 1


                        if nearby_atom.get_id() == "N":
                            atom_counts["donor"] +=1

                        if nearby_atom.get_id() == "O":
                            atom_counts["acceptor"] +=1

                        if nearby_atom.get_id() == "C":
                            atom_counts["aromatic"] +=1

                    # Calculate charges:
                        
                        resname = nearby_atom.get_parent().get_resname().capitalize()
                        if resname in atom_dict.charges and nearby_atom.get_id() in atom_dict.charges[resname]:
                            atom_counts["charge"] += atom_dict.charges[resname][nearby_atom.get_id()]
                       
                       
                        atom_proportions = {atom_id: count / atom_counter for atom_id, count in atom_counts.items()}
                        
                        # create label of atom in lbs

                        atom_proportions["environment_density"] = atom_counter/sphere_volume

                        # create a label for the residues
                        
                        atom_proportions["residue"] = residue.get_resname() + str(residue.get_id()[1])

                        # use a chimera recognisable code as key 
                            
                        chimera_code = ":" + str(residue.get_id()[1]) + "@" + selected_atom.get_name()
                        environments[chimera_code] = atom_proportions

        return environments
            




