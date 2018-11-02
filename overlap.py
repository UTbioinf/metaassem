# Currently not being used

class Overlap:
    def __init__(self, contig1, contig2, contig1_span, contig2_span):
        self.contig1 = contig1
        self.contig2 = contig2
        self.contig1_span = contig1_span # 2-tuple
        self.contig2_span = contig2_span # 2-tuple
    def __repr__(self):
        return "%s and %s" % (self.contig1, self.contig2)
