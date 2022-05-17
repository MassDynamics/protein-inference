class ProcessProteinInferenceOutput():
    
    def get_major_proteins(self, protein_table):
        protein_table = self.convert_protein_table_to_uniprot(protein_table)
        protein_ids_set = set(protein_table.Group.to_list())
        protein_ids_set = {i for i in protein_ids_set if type(i) == str}
        return protein_ids_set

    def get_protein_group(self, protein_table, group):
        group_members = protein_table.protein_id[protein_table.Group == group]
        return set(list(group_members))

    def convert_protein_table_to_uniprot(self, protein_table):
        '''Uses regex to retrieve uniprot id's and
        return updated protein table'''

        protein_table.protein_id = protein_table.protein_id.str.extract('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', 
                                    expand=False)[0]

        protein_table.Group = protein_table.Group.str.extract('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', 
                                    expand=False)[0]
        return protein_table

    def get_protein_groups_dictionary(self, protein_table):
        '''for a protein table, creates a dictionary of the set
        of proteins contained in a protein group, labeled by
        major protein'''

        protein_groups = dict()
        for protein in self.get_major_proteins(protein_table):
            protein_groups[protein] = self.get_protein_group(protein_table,protein)

        return protein_groups

    def convert_peptide_table_to_uniprot(self, peptide_table):
        '''Uses regex to retrieve uniprot id's and
        return updated protein table'''

        peptide_table.major = peptide_table.major.str.extract('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', 
                                    expand=False)[0]


        return peptide_table
