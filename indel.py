import util

class Indel:
    def __init__(self, contig1, contig2, name=None):
        self.contig1 = contig1
        self.contig2 = contig2
        self.a, self.b = contig1.overlaps[contig2][0][0]     
        self.c, self.d = contig1.overlaps[contig2][1][0]
        self.x, self.y = contig1.overlaps[contig2][0][1]
        self.z, self.w = contig1.overlaps[contig2][1][1]
        self.name = str(self.contig1) + " indel " + str(self.contig2)        
        if self.x < self.w:
            self.direction = "forward"
        elif self.x > self.w:
            self.direction = "reverse"
        else:
            self.direction = "confused"
            
        if self.direction == "forward":
            self.R = contig1.sequence[self.c-1:self.b-1]
            self.I = contig2.sequence[self.y-1:self.z-1]
        elif self.direction == "reverse":
            self.R = util.reverse_complement(contig1.sequence[self.b-1:self.c-1])
            self.I = util.reverse_complement(contig2.sequence[self.z-1:self.y-1])
            
    def __repr__(self):
        return str(self.contig1) + " indel " + str(self.contig2)
