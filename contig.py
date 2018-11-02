class Contig:
    def __init__(self, name, origin, length=0, sequence=""):
        self.name = name
        self.length = int(length) # coerces into being an int
        self.overlaps = {}
        self.sequence = sequence
        self.origin = origin        # expects reference, query, or indel
        # TODO: add check such that only these values are allowed?
        self.unique_sequence = None
        self.left = [None]
        self.right = [None]
        self.reversed = False
    def __repr__(self):
        return "Contig %s" % self.name
#        return "<Contig Name:%s Length: %s>" % (self.name, self.length)
    def __str__(self):
        return Contig.__repr__(self)
#        return "From str method of Contig: name is %s, length is %s" % (self.name, self.length)
#    def __add__(self, other):
#        pass
#        TODO: implement this. two Contigs add to make a new Assembly; Contig + Assembly => appends Contig to Assembly
#    def __sub__(self, other):
#        pass
#        TODO: implement this; takes Assembly as "other" and removes self from Assembly if present, otherwise fails
    def overlap(self, other, self_span, other_span):
        try:
            self.overlaps[other] += [(self_span, other_span)]
        except:
            self.overlaps[other] = [(self_span, other_span)]
            
    def update_overlaps(self, original_name, new_name=None, self_span=None, other_span=None):
        other = original_name
        if new_name != None:
            print self.overlaps
            print self.overlaps[other]
            self.overlaps[new_name] = self.overlaps[other]
            del(self.overlaps[other])
        if self_span != None:
            pass
        if other_span != None:
            pass
        
    def get_unique_sequence(self):
        '''Returns unique sequence. If not already determined, computes and assigns it.'''
        if self.unique_sequence == None:
            if self.overlaps == {}:
                self.unique_sequence = self.sequence
            else:
                if self.left == [None]:
                    begin, end = self.overlaps.values()[0][0][0]
                    self.unique_sequence = self.sequence[:begin]
                elif self.right == [None]:
                    begin, end = self.overlaps.values()[0][0][0]
                    self.unique_sequence = self.sequence[end:]
                else:
                    # in the world of many neighbors, this only works for
                    # things with exactly one neighbor...should this change?
                    
                    # quick fix:
                    begin_r, end_r = self.overlaps[self.right[0]][0][0]
                    begin_l, end_l = self.overlaps[self.left[0]][0][0]
                    
                    self.unique_sequence = self.sequence[end_l:begin_r]
        return self.unique_sequence
    def get_overlap_sequence(self, other):
        '''Returns overlap sequence. Does not assign.'''
        try:
            begin, end = self.overlaps[other][0][0]
        except:
            # print "KeyError"
            pass
            # TODO: should this throw an error?
        else:
            return self.sequence[begin:end]
    
    def add_neighbor(self, other, side): # add , direction?
        if side == "right":
            self.add_right_neighbor(other)
        elif side == "left":
            self.add_left_neighbor(other)
        else:
            raise ContigError("specified side does not exist")

    def add_right_neighbor(self, other):
        if self.right[0] == None:
            self.right[0] = other
        else:
            if other not in self.right:
                self.right.append(other)
                
    def add_left_neighbor(self, other):
        if self.left[0] == None:
            self.left[0] = other
        else:
            if other not in self.left:
                self.left.append(other)
    # these methods should ensure that something's left or right is only [None]
    # if it truly has no neighbors on that side
    # also need to double-check they aren't already neighbors?
        
    def get_neighbors(self):
        for name, alignment in self.overlaps.iteritems():
            x, y = alignment[0][0]
            a, b = alignment[0][1]
            self_length     = self.length
            print name
#            other_length    = 
            if x == 1:
                self.left = name
            elif a == 1:
                self.right = name
# this method has been obsoleted
# will eventually need to adjust to include the 0 0 0
