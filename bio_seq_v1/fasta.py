import warnings

class sequence():
    """
    Represents a biological sequence (DNA/RNA) with methods to calculate 
    length, base counts, GC content, and reverse complement.

    Attributes:
        id (str): Identifier for the sequence.
        sequence (str): Uppercase string of the sequence.
    """

    revcomp_dict = {
        "A":"T", "T":"A", "G":"C", "C":"G", "U":"A",
        "R":"Y", "Y":"R", "S":"S", "W":"W",
        "K":"M", "M":"K", "B":"V", "D":"H",
        "H":"D", "V":"B", "N":"N"
    }
    
    valid = "ACGTUNRYSWKMBDHV-."

    def __init__(self, id, sequence):
        """
        Initialize a Sequence object with an ID and sequence string.

        Args:
            id (str): Identifier for the sequence.
            sequence (str): Sequence string containing valid bases.

        Raises:
            ValueError: If the sequence is empty or contains invalid characters.
        """
        if not sequence:
             raise ValueError(f"Sequence for ID '{id}' is empty")
        invalid_chars = set(sequence.upper())-set(self.valid)
        if invalid_chars:
             raise ValueError(f"Sequence '{id}' contains invalid characters: {invalid_chars}")
        self.id = id 
        self.sequence = sequence.upper()
             
    def sequence_lengths(self):
        """
        Return the length of the sequence.

        Returns:
            int: Number of bases in the sequence.
        """
        return len(self.sequence)
    
    def base_count(self):
        """
        Count the occurrences of each valid base in the sequence.

        Returns:
            dict: Dictionary with bases as keys and counts as values.

        Notes:
            Issues a warning if invalid characters are found.
        """
        counts = {b:0 for b in self.valid}
        for b in self.sequence:
            if b in counts:
                counts[b] += 1
            else:
                 warnings.warn(f"Invalid character in id: '{self.id}' and sequence: '{self.sequence}")
        return counts
    
    def gc_content(self):
        """
        Calculate the GC content percentage of the sequence.

        Returns:
            float: GC content as a percentage.

        Raises:
            ValueError: If the sequence contains no valid bases.
        """
        if not self.sequence:
             return 0.0
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        total = sum(self.sequence.count(b) for b in self.valid)
        if total == 0:
             raise ValueError("Sum of total bases is 0 due to invalid bases as input")
        return ((g+c)/total)*100

    def rev_complement(self):
        """
        Compute the reverse complement of the sequence.

        Returns:
            str: Reverse complement string.

        Raises:
            ValueError: If the sequence contains invalid bases.
        """
        reverse = self.sequence[::-1]
        complement = ""
        for b in reverse:
            if b in self.revcomp_dict:
                  complement += self.revcomp_dict[b]
            else:
                 raise ValueError(f"Sequence has invalid characters, id: '{self.id}', sequence: '{self.sequence}'")
        return complement