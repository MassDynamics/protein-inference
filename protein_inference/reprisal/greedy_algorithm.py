from operator import itemgetter
from protein_inference.table_maker import TableMaker
from multiprocessing import cpu_count, Pool


class GreedyAlgorithm:

    def run(self, pn):  # should this be >1 function?

        # initialize allocation property
        for peptide in pn.get_peptides():
            pn.network.nodes[peptide]["allocated"] = 0

        # initialize major property
        for protein in pn.get_proteins():
            pn.network.nodes[protein]["major"] = 0

        n_allocated = 0
        circuit_breaker = 0
        max_iterations = 100

        while n_allocated < len(pn.get_peptides()):

            # Step 1. identify allocatee
            if len(self.get_unique_non_major_proteins(pn)) > 0:
                unique_only = 1
            else:
                unique_only = 0

            allocatee = self.get_highest_scoring_protein(pn, unique_only)

            # Step 2: Update Network

            # 2.0 Call Protein Major
            pn.network.nodes[allocatee]["major"] = allocatee

            # 2.1 Score
            pn.network.nodes[allocatee]["score"] = self.score_protein(
                pn, allocatee)

            # 2.2 Allocate Peptides:
            pn = self.peptide_allocator(pn, allocatee)  # dangerous line

            # 2.3 Alocate Proteins:
            pn = self.protein_allocator(pn, allocatee)

            # Step 3. Check if you can finish
            n_allocated = sum(
                [i != 0 for i in pn.get_node_attribute_dict("allocated").values()])

            # Circuit Breaker.
            circuit_breaker = circuit_breaker + 1
            if circuit_breaker > max_iterations:
                print("Warning: Problem Network including protein: ",
                      pn.get_proteins()[0], "couldn't be solved")
                break

        # set all remaining protein scores to 0
        for protein in pn.get_proteins():
            if "score" not in pn.network.nodes[protein].keys():
                pn.network.nodes[protein]["score"] = 0

        pn = self._tag_razor(pn)

        return pn

    def run_system_wide(self, pns):

        p = Pool(cpu_count())
        solved_pns = p.map(self.run, pns)

        return solved_pns

    def get_highest_scoring_protein(self, pn, unique_only=False):

        score_dict = self.score_all_proteins(pn, unique_only=unique_only)
        best_scoring_protein = max(score_dict.items(), key=itemgetter(1))[0]
        # deal with same score proteins! (id's should be unique even if scores aren't)
        best_scoring_proteins = [
            k for k, v in score_dict.items() if v == score_dict[best_scoring_protein]]

        return sorted(best_scoring_proteins)[0]

    def score_all_proteins(self, pn, unique_only=False):
        score_dict = {}
        if unique_only:
            proteins = pn.pick_nodes("unique_evidence", True)
        else:
            proteins = pn.get_proteins()

        for protein in proteins:
            score_dict[protein] = self.score_protein(pn, protein)

        return score_dict

    def score_protein(self, pn, protein):
        '''uses large arbitrary number to 
        ensure uniqueness dominates scores on non-unique proteins '''

        score = 0

        for peptide in pn.network.neighbors(protein):
            # otherwise only add score for non-allocated peptides
            if not pn.network.nodes[peptide]["allocated"]:
                # I've left architecture here for treating unique and non-unique
                # peptides differently.
                if pn.network.nodes[peptide]["unique"]:
                    score = score + pn.network.edges[peptide, protein]["score"]
                else:
                    score = score + pn.network.edges[peptide, protein]["score"]

        # this is just to ensure proteins that are allocated don't get picked.
        if protein in pn.get_node_attribute_dict("allocated").values():
            return -10

        return score

    def peptide_allocator(self, pn, protein):
        '''Allocates all peptide neighbours, not already
        claimed to this protein'''

        # peptides need to be allocated
        for peptide in pn.network.neighbors(protein):
            if not pn.network.nodes[peptide]["allocated"]:
                pn.network.nodes[peptide]["allocated"] = protein

        return pn

    def protein_allocator(self, pn, allocatee):

        # get list of proteins that are neighbours
        # of your proteins peptides and not allocated
        peptide_neighbours = set(pn.network.neighbors(allocatee))
        protein_over_neighbours = []
        for peptide in peptide_neighbours:
            for protein_neighbour in pn.network.neighbors(peptide):
                protein_over_neighbours.append(protein_neighbour)

        protein_over_neighbours = set(protein_over_neighbours)
        unallocated_proteins = set(pn.pick_nodes("major", 0))

        # for each protein in set of overneighbours
        # if their unallocated peptides are a subset of
        # the peptide list, set their group:

        unnallocated_peptides = pn.pick_nodes("allocated", 0)

        for protein_neighbour in protein_over_neighbours:
            if protein_neighbour in unallocated_proteins:
                peptide_set = set(pn.network.neighbors(protein_neighbour))
                peptide_set = peptide_set.intersection(unnallocated_peptides)
                if peptide_set.issubset(peptide_neighbours):
                    pn.network.nodes[protein_neighbour]["major"] = allocatee

        return pn

    def get_unique_non_major_proteins(self, pn):
        set_of_unique_evidenced = set(pn.pick_nodes("unique_evidence", True))
        set_of_allocated = set(pn.get_node_attribute_dict("major").values())
        set_unique_non_major = [
            p for p in set_of_unique_evidenced if p not in set_of_allocated]
        return set_unique_non_major

    def _tag_razor(self, pn):

        for peptide in pn.get_peptides():
            neighbour_groups = []
            for protein in pn.network.neighbors(peptide):
                neighbour_groups.append(pn.network.nodes[protein]["major"])
            if len(set(neighbour_groups)) > 1:  # could be dense but clearer this way
                pn.network.nodes[peptide]["razor"] = False
            else:
                pn.network.nodes[peptide]["razor"] = True

        return pn
