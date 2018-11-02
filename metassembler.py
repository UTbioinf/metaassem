#!/usr/bin/env python26
"""Metassembler
put more usage info here>
"""

__author__ = "Paul Baranay (pbaranay@cshl.edu)"

import sys
import getopt
import subprocess
import util
# import argparse    # TODO: redo argparser
import math
from contig import Contig
from assembly import Assembly
from indel import Indel
from util import printwithtime

# TODO: if delta file/alignments is empty, abort

class Run:
    def __init__(self, argv):
        self.argv = argv
        self.input_delta = None
        self.delta_file = None
        self.delta_filter_file = None
        self.keep_delta = True
        self.give_output = True
        self.output_prefix = "out"
        # self.output_delta = None # TODO: ?
        self.reference = None
        self.query = None
        self.silent = True # for nucmer	# TODO: currently no way for user to turn this off outside of code-diving
        self.alignments = []
        self.contigs = {}
        self.assemblies = {} # of form {0: Assembly()}
        self.process_indels_bool = True
        self.perform_assembly_bool = True
    def main(self):
        printwithtime("Welcome to the metassembler.")
        self.process_arguments()
        self.check_input()
        self.get_delta()
        self.process_delta()
        self.contigs_from_fasta()
        self.contigs_from_alignments() # TODO: these two methods should be merged?

        if self.process_indels_bool: self.process_indels()

        self.get_all_neighbors()
        if self.perform_assembly_bool:
            self.assemble()
        else:
            for contig in self.contigs.values(): self.assemblies[contig.name] = Assembly([contig])
        
        for assembly in self.assemblies.values():
            assembly.get_sequence()
        self.cleanup()
        self.summary()
        self.write_output()
        printwithtime("Done.")
        
    def usage(self):
        # TODO
        print "Usage: metassembler [OPTIONS] <reference> <query>"
        print "Options: [-o output_prefix] [-d delta_file] [-n] [-i] [-a]"
        print "Providing -n causes delta files to not be saved."
        print "Providing -i causes indel processing to be skipped."
        print "Providing -a causes assembly to be skipped"
        
    def process_arguments(self):
        printwithtime("Processing arguments...")
        try:
            opts, args = getopt.getopt(self.argv, "hknaid:o:", ["help", "keep-delta", "no-keep-delta", "no-assembly", "no-indel", "delta-input=" "output="])
        except getopt.GetoptError:
            print "Unreocgnized option"
            self.usage()
            sys.exit(2)
        input_provided = False
    	try:
    	    self.reference, self.query = args
    	    input_provided = True
    	except:
    	    print "Must specify reference and query" # TODO: more informative
    	    self.usage()
    	    sys.exit(2)
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.usage()
                sys.exit()
            if opt in ("-k", "--keep-delta"): # TODO: default to keeping delta
                self.keep_delta = True
            if opt in ("-n", "--no-keep-delta"): # TODO: default to keeping delta
                self.keep_delta = False
            if opt in ("-o", "--output"):
                self.give_output = True
                self.output_prefix = arg
            if opt in ("-d", "--delta-input"):
                self.input_delta = arg
                self.delta_file = arg
            if opt in ("-a", "--no-assembly"):
                print "Will NOT perform assembly" # TODO: more descriptive/put this elsewhere
                self.perform_assembly_bool = False
            if opt in ("-i", "--no-indel"):
                print "Will NOT process indels" # TODO: more descriptive/put this elsewhere
                self.process_indels_bool = False
        if not input_provided:
            print "Must specify input file"
            self.usage()
            sys.exit(2)
        print("Reference input file: \t %s") % self.reference
        print("Query input file: \t %s") % self.query
        print("Input delta file: \t %s") % self.input_delta
        
    def check_input(self):
        # TODO: check that delta-file appears to be a delta file
        # exit gracefully if files don't exist
        try:
            open(self.reference)
        except IOError:
            self.error("Reference file '%s' does not exist." % self.reference)
        try:
            open(self.query)
        except IOError:
            self.error("Query file '%s' does not exist." % self.query)
        if self.input_delta:
            try:
                open(self.input_delta)        
            except IOError:
                self.error("Delta file '%s' does not exist." % self.input_delta)
                print "Error with input file(s)."
                self.usage()
                sys.exit(2)

    def error(self, message):
        print message
        sys.exit(2)

    def get_delta(self):
        if self.delta_file == None:
            printwithtime("Running nucmer...")
            self.nucmer()
        else:
            self.delta_file = self.input_delta
            self.keep_delta = True # prevent you from deleting your own delta file
        
        print("Delta file: \t\t %s") % self.delta_file
        self.filter()
    
    def nucmer(self):
        # TODO: deal gracefully with nucmer errors
        prefix = "metassembler"
        if self.silent:
            f = open("/dev/null")
            subprocess.call(["nucmer", "-prefix", prefix, self.reference, self.query], stderr=f)
            f.close()
        else:
            subprocess.call(["nucmer", "-prefix", prefix, self.reference, self.query])
            # TODO: add support for passing options to nucmer
        self.delta_file = prefix + ".delta"
        printwithtime("Nucmer finished.")
    
    def filter(self):
        printwithtime("Filtering reads...") # TODO: output specifications of filtering here
        self.delta_filter(self.delta_file, self.delta_file + ".filtered")
        if self.keep_delta:
            printwithtime("Delta file will be kept at end of run.")
        else:
            printwithtime("Delta file will be deleted at end of run.")
        
    def delta_filter(self, input, output):
        printwithtime("Filtering delta file...")
        self.delta_filter_file = output
        f = open(output, "w")        
        subprocess.call(["delta-filter", "-1", "-i", "95", input], stdout=f)
        f.close()
        print("Filtered delta file: \t %s" % self.delta_filter_file)
         
    def process_delta(self):
        printwithtime("Processing delta file...")
        printwithtime("Extracting alignments...")
        f = open(self.delta_filter_file, "r")
        lines = f.readlines()
        f.close()
        lines = lines[2:] # discard first two (useless) lines # TODO: add check that they are actually useless
        lines = [x.strip('\n') for x in lines] # strip \n off each line
        first_chars = [line[0] for line in lines]
        entries = first_chars.count(">")
        alignments = []
        for i in xrange(entries):
            try:
                breakpoint = first_chars.index(">", 1)
            except:
                breakpoint = len(lines)
            alignments.append(lines[0:breakpoint])
            lines = lines[breakpoint:]
            first_chars = first_chars[breakpoint:]
        new_alignments = []
        for alignment in alignments:
            try:
                zero_index = alignment.index("0")
                if zero_index + 1 != len(alignment):
                    breakpoint = zero_index
                    new_align1 = alignment[0:breakpoint+1]
                    new_align2 = [alignment[0]] + alignment[breakpoint+1:]
                    new_alignments.append(new_align1)
                    new_alignments.append(new_align2)
                    # TODO: this may fail horribly if needs to split twice?
                else:
                    new_alignments.append(alignment)
            except:
                printwithtime("fail") # TODO: better error message
                pass
        self.alignments = new_alignments
        # parse through alignments, find any with multiple zeros
        # do operations on those (split into two)
        # then chop out all zeroes

    def contigs_from_fasta(self):
        print "Creating contigs..."
        for file in [self.reference, self.query]:
            if file == self.reference:
                type = "reference"
            elif file == self.query:
                type = "query"            
            f = open(file)
            lines = f.readlines()
            for line in lines:
                if line[0] == ">":
                    name = line[1:-1]
                    name = name.split(" ")[0] # sanitize name so it matches with results from contigs_from_alignments
                    contig = Contig(name, type)
                    self.contigs[name] = contig
                else:
                    sequence = line.split("\n")[0]
                    self.contigs[name].sequence += sequence
            for name, contig in self.contigs.items():
                self.contigs[name].length = len(self.contigs[name].sequence)
            f.close()
        print("%i contigs created") % len(self.contigs)

    def contigs_from_alignments(self):
        print "Adding alignments to contigs..."
        for i in xrange(len(self.alignments)):
            parts = self.alignments[i][0].split(" ")
            c1 = Contig(parts[0][1:], "reference", int(parts[2]))   # ref contig   # [1:] cleans off ">" character
            c2 = Contig(parts[1], "query", int(parts[3]))           # query contig
            # sanitize name so it matches with results from contigs_from_fasta
            c1.name = c1.name.split(" ")[0]
            c2.name = c2.name.split(" ")[0]
            parts2 = self.alignments[i][1].split(" ")
            # checks if these contigs already exist
            if c1.name in self.contigs:
                c1 = self.contigs[c1.name]
            if c2.name in self.contigs:
                c2 = self.contigs[c2.name]
            # adds overlap information to contig
            c1.overlap(c2, (int(parts2[0]), int(parts2[1])), (int(parts2[2]), int(parts2[3])))
            c2.overlap(c1, (int(parts2[2]), int(parts2[3])), (int(parts2[0]), int(parts2[1])))
            # adds contigs to dictionary for storage
            self.contigs[c1.name] = c1
            self.contigs[c2.name] = c2

    def get_all_neighbors(self):
        for self_contig in self.contigs.values():
            for other_contig, alignment in self_contig.overlaps.iteritems():
                x, y = alignment[0][0]
                a, b = alignment[0][1]
                self_length     = self_contig.length
                other_length    = self.contigs[other_contig.name].length
                if a < b:       # same direction
                    if x == 1:
                        make_neighbors(other_contig, self_contig)
                        # case 2
                    elif a == 1:
                        # case 1
                        make_neighbors(self_contig, other_contig)
                elif a > b:    # 'other' is in reversed direction
                    if x == 1:
                        make_neighbors(self_contig, other_contig)
                        if not other_contig.reversed:
                            # TODO: this...probably needs to be more sensitive 
                            other_contig.sequence = util.reverse_complement(other_contig.sequence)
                            other_contig.reversed = True                      
                        # case 4
                    elif a == other_length:
                        make_neighbors(other_contig, self_contig)
                        if not other_contig.reversed: 
                            other_contig.sequence = util.reverse_complement(other_contig.sequence)
                            other_contig.reversed = True                      
                        # case 3                        

            
#        for contig in self.contigs.values():
#            if contig.left == None and contig.right == None:
#                contig.get_neighbors()

    def process_indels(self):
        # TODO: add "number of indels collapsed", "number of contigs deleted"
        printwithtime(len(self.contigs))
        printwithtime("Searching for indels...")
        self.indels = []
        for key, contig in self.contigs.items():
            for overlap in contig.overlaps.items():
                if len(overlap[1]) == 2:
                    other_name = overlap[0].name
                    print overlap[1], contig.name, other_name
                    other_contig = self.contigs[other_name]
                    new_indel = Indel(contig, other_contig)
                    if contig.length > other_contig.length:
                        self.indels.append(new_indel)
                    elif contig.length == other_contig.length:
                        printwithtime("congrats, you've found an edge case; please report this")
                        # TODO: improve this                
                elif len(overlap[1]) >= 3:
                    # TODO: fix reason why three overlaps aren't being stored
                    # e.g. the A9-B8 issue for test case 7
                    printwithtime("something weird is going on! please report this")
                    # TODO: improve this
        
        # this processes indels from right to left
        indel_dict = {}
        for item in set([x.contig1.name for x in self.indels]):
            indel_dict[item] = {}
        for indel in self.indels:
            indel_dict[indel.contig1.name][indel.d] = indel
        for contig_set in indel_dict.values():            
            indel_keys = contig_set.keys()
            indel_keys.sort(reverse=True)

            for key in indel_keys:
                # TODO:if contigs from ref & query have same names, badness will ensue
                indel = contig_set[key]
                print indel.name, indel.direction
                print indel.a, "\t", indel.b, "\t", indel.x, "\t", indel.y
                print indel.c, "\t", indel.d, "\t", indel.z, "\t", indel.w
                print "R: ", len(indel.R)
                print "I: ", len(indel.I)
                
                indel_sequence = indel.R + indel.I
                original = self.contigs[indel.contig1.name].sequence
                self.contigs[indel.contig1.name].sequence = original[:indel.c] + indel_sequence + original[indel.c:]
                self.contigs[indel.contig1.name].length   = len(self.contigs[indel.contig1.name].sequence)
                self.contigs[indel.contig1.name].origin   = "indel" # TODO: adjust this to keep a running trail of added contigs
                try:
                    del(self.contigs[indel.contig2.name].overlaps[indel.contig1])
                    if len(self.contigs[indel.contig2.name].overlaps) == 0:
                        print("Deleted contig %s\n") % indel.contig2.name
                        del(self.contigs[indel.contig2.name])
                    else:
                        pass
                        # update overlaps
                        # MAJOR TODO: also, need to handle updating overlap info for contig1
                    del(self.contigs[indel.contig1.name].overlaps[indel.contig2])
                except:
                    printwithtime("Hey, you've probably found an error, please report this")
                    # TODO: make this better              
                
#                if indel.signature == "basic":
#                    self.basic_indel(indel)
#                elif indel.signature == "collapse":
#                    self.collapse_indel(indel)
#                elif indel.signature == "skip":
#                    print("skipping indel %s") % indel.name

#####    OLD IMPLEMENTATION
#                my_sequence     = indel.contig1.sequence
#                other_sequence  = indel.contig2.sequence
#                indel_start     = my_sequence[:indel.b]
#                if indel.direction == "forward":
#                    indel_middle = other_sequence[indel.y:indel.z-1]
#                elif indel.direction == "reverse":
#                    indel_middle  = util.reverse_complement(other_sequence[indel.z:indel.y-1])
#                elif indel.direction == "skip":
#                    indel_middle = ""
#                    printwithtime("This contig is odd")
#                indel_end       = my_sequence[indel.c-1:]
#                indel_sequence  = indel_start + indel_middle + indel_end
#                indel_len       = indel.z - indel.y - 1
#                L = indel_len
#                self.contigs[indel.contig1.name].sequence = indel_sequence
#                self.contigs[indel.contig1.name].length   = len(indel_sequence)
#                self.contigs[indel.contig1.name].origin   = "indel" # TODO: adjust this to keep a running trail of added contigs
#                try:
#                    del(self.contigs[indel.contig2.name].overlaps[indel.contig1])
#                    if len(self.contigs[indel.contig2.name].overlaps) == 0:
#                        print("Deleted contig %s") % indel.contig2.name
#                        del(self.contigs[indel.contig2.name])
#                    else:
#                        pass
#                        # update overlaps
#                        # MAJOR TODO: also, need to handle updating overlap info for contig1
#                    del(self.contigs[indel.contig1.name].overlaps[indel.contig2])
#                except:
#                    printwithtime("Hey, you've probably found an error, please report this")
#                    # TODO: make this better
#####

#                print("Updated %s with indel from %s\n") % (indel.contig1.name, indel.contig2.name)
        printwithtime(len(self.contigs))

    def basic_indel(self, indel):
        my_sequence     = indel.contig1.sequence
        other_sequence  = indel.contig2.sequence
        indel_start     = my_sequence[:indel.b]
        if indel.direction == "forward":
            indel_middle = other_sequence[indel.y:indel.z-1]
        elif indel.direction == "reverse":
            indel_middle  = util.reverse_complement(other_sequence[indel.z:indel.y-1])
        indel_end       = my_sequence[indel.c-1:]
        indel_sequence  = indel_start + indel_middle + indel_end
        indel_len       = indel.z - indel.y - 1
        L = indel_len
        self.contigs[indel.contig1.name].sequence = indel_sequence
        self.contigs[indel.contig1.name].length   = len(indel_sequence)
        self.contigs[indel.contig1.name].origin   = "indel" # TODO: adjust this to keep a running trail of added contigs
        try:
            del(self.contigs[indel.contig2.name].overlaps[indel.contig1])
            if len(self.contigs[indel.contig2.name].overlaps) == 0:
                print("Deleted contig %s") % indel.contig2.name
                del(self.contigs[indel.contig2.name])
            else:
                pass
                # update overlaps
                # MAJOR TODO: also, need to handle updating overlap info for contig1
            del(self.contigs[indel.contig1.name].overlaps[indel.contig2])
        except:
            printwithtime("Hey, you've probably found an error, please report this")
            # TODO: make this better

    def collapse_with_insertion_indel(self, indel):
        try:
            del(self.contigs[indel.contig2.name].overlaps[indel.contig1])
            if len(self.contigs[indel.contig2.name].overlaps) == 0:
                print("Deleted contig %s") % indel.contig2.name
                del(self.contigs[indel.contig2.name])
            else:
                pass
                # update overlaps
                # MAJOR TODO: also, need to handle updating overlap info for contig1
            del(self.contigs[indel.contig1.name].overlaps[indel.contig2])
        except:
            printwithtime("Hey, you've probably found an error, please report this")
            # TODO: make this better        

    def assemble(self):
        printwithtime("Assembling contigs...")
        
        contigs_copy = self.contigs.copy()
        
        # greedily finds "anchor" contigs (with no neighbors on one side)
        for key, contig in self.contigs.items():
            if contig.left == [None] or contig.right == [None]:
                index = len(self.assemblies)
                assembly = Assembly([contig])
                assembly.left  = contig.left
                assembly.right = contig.right
                self.assemblies[index] = assembly
                del(contigs_copy[key])
                print("Assigning %s to assembly %i // %i contigs remaining to assign") % (key, index, len(contigs_copy))
        
        self.contigs = contigs_copy

        num_contigs_list = [None, None]
        while contigs_copy != {}:
            if num_contigs_list[0] == num_contigs_list[1] != None:
                break
            else:
                for index, assembly in self.assemblies.items():
                    for key, contig in self.contigs.items():
                        if assembly.assembly[-1] in contig.left:
                            if len(assembly.right) == 1 and len(contig.left) == 1:                        
                                assembly.assembly.append(contig)
                                assembly.right = contig.right 
                                del(contigs_copy[key])
                                print("Assigning %s to assembly %i // %i contigs remaining to assign") % (key, index, len(contigs_copy))
                        elif assembly.assembly[0] in contig.right:
                            if len(assembly.left) == 1 and len(contig.right) == 1:                            
                                assembly.assembly.insert(0, contig)
                                assembly.left = contig.left
                                del(contigs_copy[key])
                                print("Assigning %s to assembly %i // %i contigs remaining to assign") % (key, index, len(contigs_copy))
                num_contigs_list[0] = num_contigs_list[1]
                num_contigs_list[1] = len(contigs_copy)

# TODO: list value of xrange(100) [1st pass, 2nd pass...]
# TODO: make it more clear when it enters this phase

####### handle any unassigned contigs
        if contigs_copy != {}:
            printwithtime("Assigning unassignable contigs")
            index = int(max(self.assemblies.keys())) # ensures that nothing is being overwritten
            # TODO: add a try/except loop that double-checks nothing is being overwritten
            for key, contig in self.contigs.items():
                index += 1
                self.assemblies[index] = Assembly([contig], contig.left, contig.right)
                del(contigs_copy[key])
                print("Assigning %s to assembly %i // %i contigs remaining to assign") % (key, index, len(contigs_copy))

####### merging assemblies
# TODO: current issue, sometimes combines one assembly into multiple places...
        printwithtime("Merging assemblies...")
        num_merges = 0
        assemblies_copy = self.assemblies.copy()
        for asmb1 in assemblies_copy.items():
            for asmb2 in assemblies_copy.items():
                if asmb1 != asmb2 and asmb1 in self.assemblies.items():
                    # the and clause helps ensure that assemblies are not merged into
                    # assemblies that no longer exist
                    asmb1_index = asmb1[0]
                    asmb1_asmb  = asmb1[1]
                    asmb2_index = asmb2[0]
                    asmb2_asmb  = asmb2[1]

                    if (asmb2_asmb.assembly[0] in asmb1_asmb.right and asmb1_asmb.right != [None] and asmb2_asmb.left != [None]) or \
                        (asmb1_asmb.assembly[-1] in asmb2_asmb.left and asmb2_asmb.left != [None] and asmb1_asmb.right != [None]):
                        if len(asmb1_asmb.right) == 1 and len(asmb2_asmb.left) == 1:
                            print("Merging assembly %s into assembly %s") % (asmb2_index, asmb1_index)
  #                          print("Old assemblies: %s and %s") % (asmb1_asmb, asmb2_asmb)
                            asmb1_asmb.assembly.extend(asmb2_asmb.assembly)
  #                          print("New assembly: %i  %s") % (asmb1_index, asmb1_asmb)
                            asmb1_asmb.right = asmb2_asmb.right
                            num_merges += 1
                            try:
                                del(assemblies_copy[asmb2_index])
                                del(self.assemblies[asmb2_index])
                            except:
                                printwithtime("exception!")
   #                         printwithtime("---")                            


                    elif (asmb1_asmb.assembly[0] in asmb2_asmb.right and asmb2_asmb.right != [None] and asmb1_asmb.left != [None]) or \
                        (asmb2_asmb.assembly[-1] in asmb1_asmb.left and asmb1_asmb.left != [None] and asmb2_asmb.right != [None]):
                        if len(asmb2_asmb.right) == 1 and len(asmb1_asmb.left) == 1:
                            print("Merging assembly %s into assembly %s") % (asmb1_index, asmb2_index)
 #                           print("Old assemblies: %s and %s") % (asmb1_asmb, asmb2_asmb)
                            asmb2_asmb.assembly.extend(asmb1_asmb.assembly)
 #                           print("New assembly: %i  %s") % (asmb1_index, asmb2_asmb)
                            asmb2_asmb.right = asmb1_asmb.right
                            num_merges += 1
                            try:
                                del(assemblies_copy[asmb1_index])
                                del(self.assemblies[asmb1_index])
                            except:
                                printwithtime("exception!")
    #                            printwithtime(asmb1_asmb, asmb1_index, asmb1_asmb.left, asmb1_asmb.right
    #                            printwithtime(asmb2, asmb2_index, asmb2_asmb.left, asmb2_asmb.right
    #                        printwithtime("---"
        print("%i merges made.") % (num_merges)


        
        # TODO:
        # clean up irrelevant code
        # redirect printwithtime(statements somewhere else...log file?
            
    def cleanup(self):
        printwithtime("Cleaning up...")
        if not self.keep_delta:
            subprocess.call(["rm", self.delta_file])
            subprocess.call(["rm", self.delta_filter_file])

    def summary(self):
        printwithtime("Writing summary...")
        f = open(self.output_prefix + ".summary.txt", "w")
        num_assemblies = len(self.assemblies)
        asmb_lengths = [len(x.assembly) for x in self.assemblies.values()]
        num_contigs = sum(asmb_lengths)
        max_length  = max(asmb_lengths)
        min_length  = min(asmb_lengths)
        avg_length  = num_contigs/num_assemblies
        # add mean?
        n50         = N50(asmb_lengths)
        f.write("%i assemblies incorporating %i contigs" % (num_assemblies, num_contigs) + "\n")
        f.write("MIN: %i \t MAX: %i \t AVG: %.2f \t N50: %.2f" % (min_length, max_length, avg_length, n50) + "\n" * 2)
        for index, assembly in self.assemblies.items():
            f.write(str(index) + " " + str(assembly) + "\n")
        f.close()
            
    def write_output(self):
        if self.give_output:
            printwithtime("Writing output...")
            f = open(self.output_prefix + ".fa", "w")
            for index, assembly in self.assemblies.items():    
                f.write(">ma_" + str(index) + " " + str(assembly) + "\n")
                f.write(assembly.sequence + "\n")
            f.close()

####################################################

def N50(values):
    '''Computes N50 using the Broad's definition.'''
    values.sort()
    prime = []
    for value in values:
        prime += [value]*value
    return median(prime)

def median(values):
    '''Calculates median without using numpy.'''
    if len(values) == 1:
        return values[0]
    middle = (1+len(values))/2.
    if math.floor(middle) == math.ceil(middle):
        return values[int(middle)]
    else:
        return (values[int(math.floor(middle))-1] + values[int(math.ceil(middle))-1]) / 2.0

####################################################            

def make_neighbors(left_contig, right_contig):
    left_contig.add_neighbor(right_contig, "right")
    right_contig.add_neighbor(left_contig, "left")

####################################################
        
if __name__ == "__main__":
    r = Run(sys.argv[1:])
    r.main()
