from Bio.PDB import PDBParser, NeighborSearch

class StructureAnalysis:
    def __init__(self, pdb_file):
        
        parser = PDBParser()
        self.structure = parser.get_structure('structure', pdb_file)
        self.standard_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        self.ns = NeighborSearch(list(self.structure.get_atoms()))

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

        if len(ligand_names) > 0:
            for residue in self.structure.get_residues():
            
                if residue.get_resname() not in ligand_names and residue.get_resname() != 'HOH':  # Exclude water residues and ligand residues
                    for selected_atom in residue.get_atoms():
                        nearby_atoms = self.ns.search(selected_atom.coord, 8)  # Change the radius if needed

                        # Filter nearby atoms that are not water and are not in residues of ligand_names
                        nearby_atoms = [atom for atom in nearby_atoms if atom.get_parent().get_resname() != 'HOH' and atom.get_parent().get_resname() not in ligand_names]
                        
                        # Initialize the atom counts and counter outside of the loop
                        atom_counts = {}
                        atom_counter = 0
                        
                        
                        atom_counts["b_factor"] = 0 # initialize b_factor

                        # iterate over the nearby atoms and count the atoms
                        for nearby_atom in nearby_atoms:
                            atom_counter += 1
                            
                            # Count the number of atoms of each type
                            if nearby_atom.get_id() in atom_counts:
                                atom_counts[nearby_atom.get_id()] += 1
                            else:
                                atom_counts[nearby_atom.get_id()] = 1

                        # Count the total number of atoms of each element
                        for nearby_atom in nearby_atoms:

                            if str("total"+nearby_atom.element) in atom_counts:
                                atom_counts["total"+nearby_atom.element] += 1
                            else:
                                atom_counts["total"+nearby_atom.element] = 1
                        
                        # Count the total number of atoms
                        for nearby_atom in nearby_atoms:
                            atom_counts["b_factor"] += nearby_atom.get_bfactor()

                        # Count the number of atoms in the environment that are considered lbs
                        atom_counts["ligand_binding_site"] = 0
                        for nearby_atom in nearby_atoms:
                            if nearby_atom.get_serial_number() in serial_numbers_lbs:
                                atom_counts["ligand_binding_site"] += 1
                            

                        atom_proportions = {atom_id: count / atom_counter for atom_id, count in atom_counts.items()}
                        
                        # create label of atom in lbs
                        if selected_atom.get_serial_number() in serial_numbers_lbs:
                            atom_proportions["is_lbs"] = "Yes"

                        else:
                            atom_proportions["is_lbs"] = "No"


                        # use a chimera recognisable code as key 
                            
                        chimera_code = ":" + str(residue.get_id()[1]) + "@" + selected_atom.get_name()
                        environments[chimera_code] = atom_proportions

            return environments
            
        else:
            print("No ligand names provided")
            return None
