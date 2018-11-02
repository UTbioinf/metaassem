# TODO: rename this scaffold?

import util

class Assembly:
    def __init__(self, assembly=[], left=[None], right=[None]):
        self.assembly = assembly
        self.sequence = []
        self.length = None
        self.left = left
        self.right = right
    def __repr__(self):
        return str(self.assembly)
    def __str__(self):
        return self.__repr__()
    def anchor(self, contig):
        if contig not in self.assembly:
            if contig.left == None: self.assembly.insert(0, contig)
            elif contig.right == None: self.assembly.append(contig)
    def insert(self, contig):
        if contig not in self.assembly:
            if contig.left in self.assembly:
                index = self.assembly.index(contig.left)
                self.assembly.insert(index + 1, contig)
            elif contig.right in self.assembly:
                index = self.assembly.index(contig.right)
                self.assembly.insert(index, contig)
            else:
                self.assembly.insert(1, contig) # just stick it in somewhere
                # assumes contig is definitely part of this assembly
    def get_sequence(self):
        # NEW APPROACH
        if len(self.assembly) >= 2:
            for contig1,contig2 in util.pairs(self.assembly):
                uniq1 = contig1.get_unique_sequence()
                uniq2 = contig2.get_unique_sequence()
                overlap = contig1.get_overlap_sequence(contig2)
                if uniq1 not in self.sequence:
                    self.sequence.append(uniq1)
                self.sequence.append(overlap)
                self.sequence.append(uniq2)
            # cleanup steps
            self.sequence = "".join(self.sequence)
            self.length = len(self.sequence)
        else:
            self.sequence = self.assembly[0].sequence