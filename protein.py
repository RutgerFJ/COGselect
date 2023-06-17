# Standard library imports
import textwrap


class Protein:
    def __init__(self, protein_id, organism, cog_id, sequence):
        """
        A class representing a Protein.

        :param protein_id: the protein identifier.
        :param organism: the Latin name for the organism the protein came from.
        :param cog_id: the identifier for the COG that contains the protein.
        :param sequence: the amino acid sequence for the protein.
        """
        self.protein_id = protein_id
        self.organism = organism
        self.cog_id = cog_id
        self.sequence = sequence

    def __repr__(self, line_length=70) -> str:
        """
        Get a FASTA formatted string displaying the protein_id, organism and
        amino acid sequence.

        :param line_length: integer for the line length to use in FASTA string.
        :return: string containing a FASTA formatted amino acid sequence.
        """
        fasta = f'>{self.protein_id}|{self.organism}\n'
        seq_lines = textwrap.wrap(self.sequence, line_length)
        fasta += '\n'.join(seq_lines) + '\n'
        return fasta

    def get_protein_id(self) -> str:
        """
        Get the identifier for the protein.

        :return: string identifier for the protein.
        """
        return self.protein_id

    def set_protein_id(self, protein_id) -> None:
        """
        Set the identifier for the protein.

        :param protein_id: string that contains the new identifier
        """
        self.protein_id = protein_id

    def get_sequence(self) -> str:
        """
        Get the amino acid sequence string that belongs to the protein.

        :return: string containing amino acid sequence of the protein.
        """
        return self.sequence

    def set_sequence(self, sequence) -> None:
        """
        Set the amino acid sequence for the protein.

        :param sequence: string that contains the new amino acid sequence.
        """
        self.sequence = sequence

    def get_length(self) -> int:
        """
        Get the integer length of the amino acid sequence.

        :return: integer length of the amino acid sequence.
        """
        return len(self.sequence)

    def get_organism(self) -> str:
        """
        Get the organism the protein came from.

        :return: string that contains the Latin name of the organism.
        """
        return self.organism

    def set_organism(self, organism) -> None:
        """
        Set the organism the protein came from.

        :param organism: string that contains the new name for the organism.
        """
        self.organism = organism

    def get_cog_id(self) -> str:
        """
        Get the string identifier for the COG.

        :return: string identifier for the COG.
        """
        return self.cog_id

    def set_cog_id(self, cog_id) -> None:
        """
        Set the string identifier for the COG the protein belongs to.

        :param cog_id: new string identifier for the COG.
        """
        self.cog_id = cog_id
