from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
import atom_dict
import math

class StructureAnalysis:
    column_names = ['b_factor', 'sasa', 'aliphatic', 'aromatic', 'donor', 'acceptor',
       'don_acc', 'CD', 'totalC', 'NZ', 'totalN', 'CD1', 'CG', 'OE2', 'totalO',
       'CD2', 'OE1', 'CB', 'N', 'C', 'O', 'CA', 'OG1', 'CG2',
       'environment_density', 'CZ', 'NE', 'NH1', 'CE1', 'CE', 'NH2', 'CE2',
       'OD1', 'ND2', 'CG1', 'SG', 'totalS', 'SD', 'CZ3', 'CH2', 'CZ2', 'CE3',
       'NE1', 'OD2', 'OH', 'OG', 'NE2', 'ND1']
        
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
        '''Returns dictionary of environments of each residue in the structure'''
    
        ligand_names = self.get_ligands_from_structure()
        ligand_binding_atoms = self.get_ligand_binding_site_atoms()
        serial_numbers_lbs = set([t[0] for t in ligand_binding_atoms])
        environments = {}
        radio = 8
        sphere_volume = 4/3 * math.pi * radio**3
        if len(ligand_names) > 0:
            for residue in self.structure.get_residues():
            
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
                        
                            # Count the total number of atoms
                        
                            atom_counts["b_factor"] += nearby_atom.get_bfactor()

                            atom_counts["sasa"] += nearby_atom.sasa

                            # Calculate the features of the atom
                            if nearby_atom.get_id() in list(atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()].keys()):
                                atom_counts[atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()][nearby_atom.get_id()]] += 1


                            if nearby_atom.get_id() == "N":
                                atom_counts["donor"] +=1

                            if nearby_atom.get_id() == "O":
                                atom_counts["acceptor"] +=1

                            if nearby_atom.get_id() == "C":
                                atom_counts["aromatic"] +=1
                            

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
        radio = 8
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
                            atom_counts[atom_dict.characteristics[nearby_atom.get_parent().get_resname().capitalize()][nearby_atom.get_id()]] += 1


                        if nearby_atom.get_id() == "N":
                            atom_counts["donor"] +=1

                        if nearby_atom.get_id() == "O":
                            atom_counts["acceptor"] +=1

                        if nearby_atom.get_id() == "C":
                            atom_counts["aromatic"] +=1
                        
                        atom_proportions = {atom_id: count / atom_counter for atom_id, count in atom_counts.items()}
                        
                        # create label of atom in lbs

                        atom_proportions["environment_density"] = atom_counter/sphere_volume

                        # use a chimera recognisable code as key 
                            
                        chimera_code = ":" + str(residue.get_id()[1]) + "@" + selected_atom.get_name()
                        environments[chimera_code] = atom_proportions

        return environments
            




