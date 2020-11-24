from pyvis.network import Network
from copy import deepcopy
import matplotlib.pyplot as plt 
from matplotlib.colors import to_hex

class NetworkGrapher():
    '''
    This class generates a pyvis graph representing a PSM network. 

    This class uses information encoded in the network component of the
    ProblemNetwork class. It creates an interactive network representation
    using the pyvis package. Several colouring options are available such
    as by status, groups or score. 

    '''

    def draw(self, problem_network, by="status", name="nx", size = [800,800]):
        '''
        Draw creates the interactive network visualization from the 
        ProblemNetwork object using pyvis and the colouring methods in 
        this classs. 

        Parameters
        ----------
        problem_network : ProblemNetwork
            An instance of the problem network object
        by : string
            A string indicating the colouring method ("status", "colour" or 
            "group")
        name : string
            A string that will be used to generate the file name when the 
            visualization is saved as an HTML file. 
        size: list
            A list of length two indicating the size in pixels of the 
            resulting visualization

        Returns
        -------
        A static HTML object that is also written to the current working directory.
        
        '''

        pn = deepcopy(problem_network)

        if by.lower() == "group":
            pn = self.colour_by_group(pn)
        elif by.lower() == "score":
            pn = self.colour_by_score(pn)
        else:
            pn = self.colour_by_status(pn)

        # make proteins bigger:
        pn.update_nodes(pn.pick_nodes("protein", 1), "size", "10")
        pn.update_nodes(pn.pick_nodes("protein", 0), "size", "5")

        # if new scores are available, highlight with scores
        score_dict = pn.get_node_attribute_dict("score")
        if score_dict:
            for protein in score_dict.keys():
                title = "Protein Score: {}".format(round(score_dict[protein], 2))
                pn.network.nodes[protein]["title"] = title

        # add peptide scores
        score_dict = pn.get_edge_attribute_dict("score")
        for psm in score_dict.keys():
            title = "PSM Score: {}".format(round(score_dict[psm], 2))
            pn.network.edges[psm]["title"] = title

        for edge in pn.network.edges():
            if edge[0] in pn.get_proteins():
                pn.network.edges[edge]["color"] = pn.network.nodes[edge[0]]["color"]
            else: 
                pn.network.edges[edge]["color"] = pn.network.nodes[edge[1]]["color"]


        nt = Network(str(size[0])+"px", str(size[1])+"px", notebook=True, heading = "")
        # populates the nodes and edges data structures
        nt.from_nx(pn.network)
        nt.save_graph(name + ".html")
        return nt.show(name + ".html")

    def colour_by_group(self, pn):
        """
        Updates the node colour attributes in a ProblemNetworks network attribute
        according to the "allocated" and "major" attributes of those nodes"

        :param: pn: the ProblemNetwork to update.
        """
        groups = list(set(pn.get_node_attribute_dict("allocated").values()))
        #colour_dict = get_colour_dict(groups)
        colormap = [
            "#a6cee3",
            "#1f78b4",
            "#b2df8a",
            "#33a02c",
            "#fb9a99",
            "#e31a1c",
            "#fdbf6f",
            "#ff7f00",
            "#cab2d6",
            "#6a3d9a",
            "#ffff99",
            "#b15928"]
        colormap = colormap[:len(groups)]
        colour_dict = dict(zip(groups, colormap))
        for group in groups:
            pn.update_nodes(pn.pick_nodes("allocated", group),
                            "color", (colour_dict[group]))
            pn.update_nodes(pn.pick_nodes("major", group),
                            "color", colour_dict[group])

        return pn

    def colour_by_status(self, pn):
        """"
        Updates the node colour attributes in a ProblemNetworks network attribute
        according to the status attributes of those nodes"

        :param: pn: the ProblemNetwork to update.
        """

        pn.update_nodes(pn.get_proteins(), "color", "red")
        pn.update_nodes(pn.get_peptides(), "color", "blue")
        pn.update_nodes(pn.pick_nodes("razor", 1), "color", "green")
        pn.update_nodes(pn.pick_nodes("unique", 1), "color", "purple")
        pn.update_nodes(pn.pick_nodes("unique_evidence", 1), "color", "orange")

        return pn

    def colour_by_score(self,pn):
        """
        Updates the node colour attributes in a ProblemNetworks network attribute
        according to the score attributes of those nodes"

        :param: pn: the ProblemNetwork to update.
        """

        #start with default status colouring
        pn = self.colour_by_status(pn)
        
        # colour proteins by score
        for protein in pn.get_proteins():
            score = pn.network.nodes[protein]["score"] 
            pn.network.nodes[protein]["color"] = to_hex(plt.get_cmap("Reds")(score))

        # colour peptides by score
        for edge in pn.network.edges():
            score = pn.network.edges[edge]["score"]
            if edge[0] in pn.get_proteins():
                pn.network.nodes[edge[1]]["color"] = to_hex(plt.get_cmap("Blues")(score))
            else: 
                pn.network.nodes[edge[0]]["color"] = to_hex(plt.get_cmap("Blues")(score))

        return pn


