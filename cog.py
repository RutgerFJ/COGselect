# Standard library imports
from statistics import mean, stdev

# Local imports
from protein import Protein


class Cog:
    def __init__(self, cog_id: str, go_annotation: str, function: str, proteins: list[Protein]):
        """
        A class representing a Cluster of Orthologous Groups (COG).

        :param cog_id: string identifier for the COG.
        :param function: string of the function of the proteins in the COG.
        :param proteins: a list of Protein objects in the COG.
        """
        self.cog_id = cog_id
        self.go_annotation = go_annotation
        self.function = function
        self.proteins = proteins

    def get_cog_id(self) -> str:
        """
        Get the string identifier for the COG.

        :return: string identifier for the COG.
        """
        return self.cog_id

    def set_cog_id(self, cog_id: str) -> None:
        """
        Set the string identifier of the COG.

        :param cog_id: new string identifier for the COG.
        """
        self.cog_id = cog_id

    def get_go_annotation(self) -> str:
        """
        Get the GO annotation of the COG.

        :return: string with the GO annotation of the COG
        """
        return self.go_annotation

    def set_go_annotation(self, go_annotation: str) -> None:
        """
        Set the GO annotation of the COG.

        :param go_annotation: new string with the GO annotation of the COG.
        """
        self.go_annotation = go_annotation

    def get_function(self) -> str:
        """
        Get the function of the COG.

        :return: string of the function of the COG
        """
        return self.function

    def set_function(self, function: str) -> None:
        """
        Set the function of the COG.

        :param function: new string of the function of the COG.
        """
        self.function = function

    def get_proteins(self) -> list[Protein]:
        """
        Get the list of proteins in the COG.

        :return: list of Protein objects.
        """
        return self.proteins

    def set_proteins(self, proteins: list[Protein]) -> None:
        """
        Set the list of proteins in the COG.

        :param proteins: new list of Protein objects.
        """
        self.proteins = proteins

    def get_avg_seq_len(self) -> float:
        """
        Get the average length of the amino acid sequences for the proteins in
        the COG.

        :return: float of the average length of the amino acid sequences.
        """
        sequences = [p.get_sequence() for p in self.proteins]
        return round(mean(map(len, sequences)), 2)

    def get_rsd(self) -> float:
        """
        Get the relative standard deviation (RSD) of the amino acid sequence
        lengths for the proteins in the COG.

        :return: float of the RSD of the amino acid sequence lengths.
        """
        sequences = [p.get_sequence() for p in self.proteins]
        return round(stdev(map(len, sequences)) / self.get_avg_seq_len(), 2)

    def get_methionines(self) -> str:
        """
        Get the amount of amino acid sequences that start with methionine in
        the COG. The string shows the count of sequences that start with 'M' in
        the COG, divided by the total number of sequences in the COG.

        :return: string showing methionine starts.
        """
        methionines = sum(1 for protein in self.proteins
                          if protein.get_sequence().startswith('M'))
        return f"[{methionines} / {len(self.proteins)}]"
