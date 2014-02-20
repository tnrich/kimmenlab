#!/usr/local/bin/python2.4
#should run with python2.4



#To put this file on Ohana:
# scp /Users/Tomaso/Documents/kimmenLab/Annotate_PDB_structure_with_CSA.py tnrich@ohana.berkeley.edu:/home/tnrich/blastp_test_dir

#To run on ohana terminal:
# python2.4 Annotate_PDB_structure_with_CSA.py -i '1a0i|A'

#Overall idea is to build an annotated dataset composed of proteins from PDB

#To allow for command line access:
import sys
import os
from pfacts003.intrepid.models import CSA
from pfacts003.intrepid.utils import read_sequence_from_pdb_file
from subprocess import call
from Bio.Blast import NCBIXML
import optparse
import unicodedata


#parse args with optparse
parser = optparse.OptionParser()

parser.add_option("-i", "--CSAInput", dest="CSA",
    help="Should be a PDB ID that is part of CSA")

opts, args = parser.parse_args()

#Input: CSA PDB ID. (We'd run this for every CSA entry, one after the next)
#Output: all PDB structures with detectable homology (based on BLAST or an HMM approach) are identified and residues in the PDB structure that agree with the CSA sequence are labelled as ?presumed_active? with the corresponding supporting evidence (BLAST E-value, %ID over the region aligned).

#STEPS:

#print 'foolll'
#print opts.CSA
csaInputRaw = opts.CSA
csaInputRHS = None;
csaInputLHS = None;

if len(csaInputRaw) == 4:
    csaInput = csaInputRaw
else:
    csaInputLHS, csaInputRHS = csaInputRaw.split("|")
    csaInput = csaInputLHS +"."+csaInputRHS
    #print csaInputLHS +"."+csaInputRHS + '   csaInputLHS . csaInputRHS'

# checks PDB structure files located at:
pdbLocation = "/clusterfs/ohana/external/pdb/"
pdbId = csaInputLHS
chainId = csaInputRHS

#get PDB sequence
sequence, numberedResidues = read_sequence_from_pdb_file(pdbId, chainId, pdbLocation)
numberedResiduesDict = dict(numberedResidues)
print sequence
print numberedResidues
print numberedResiduesDict


#parse sequence to get it into the fasta standard
n = 80
sequenceFasta=[sequence[i:i+n] for i in range(0, len(sequence), n)]
#write to a pdb.fasta file
file = open("%s.fasta" % (csaInput),'w') 
file.write('>pdb|%s mol:protein length:%s\n' % (csaInputRaw,str(len(sequence))))
file.writelines(sequenceFasta)
    


#log = open('debug.txt','w')
#log.write(csaInput + '\n')
#command = "blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry ' " + str(csaInputRaw) + "' > "+ str(csaInput) + '.fasta'
        #eg: blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry 12as > 12as.fasta
        #eg: blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry "12as|A" > 12as.A.fasta
#log.write(command + '\n')
#os.system(command)


#for each protein in csa:
#run blastp
command = "blastp -db /clusterfs/ohana/external/pdb/blastdbs/pdb -query "+ str(csaInput) + ".fasta -out "+ str(csaInput) + ".xml -outfmt 5"
log.write(command + '\n')
os.system(command)


#translation stuff:
trans = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



#grab the CSA data to have it on hand for the query sequence
#print csaInputLHS +"."+csaInputRHS + '   csaInputLHS . csaInputRHS'
csa = CSA.objects.filter(pdb_id = csaInputLHS, chain_id = csaInputRHS)
site_list=[]
csaResDict = {}
for site in csa:
    site_list.extend([site.residue_number])
    
    #print 'site.residue_number'
    #print site.residue_number
    
    resString = unicodedata.normalize('NFKD', site.residue).encode('ascii','ignore')
    resOneLetter = trans[resString.upper()]
    #print resOneLetter
    #print site.residue
    csaResDict[site.residue_number] = resOneLetter


#should check for hits 
def hitChecker(csaRes_dict,csaSeq,csaStart,pdbSeq,pdbStart):
#make sure that the csaRes is in the hsp
    csaSeq = unicodedata.normalize('NFKD', csaSeq).encode('ascii','ignore')
    pdbSeq = unicodedata.normalize('NFKD', pdbSeq).encode('ascii','ignore')
    posCounter = 0;
    csaResCounter = int(csaStart) - 1 
    pdbResCounter = int(pdbStart) - 1
    #print '  csaRes_list   '
    #print csaRes_list 
    #print '  csaSeq'
    #print csaSeq 
    #print '  csaStart'
    #print csaStart 
    #print  '  pdbSeq'
    #print pdbSeq 
    #print  '  pdbStart'
    #print pdbStart 
    results = {}
    for char in csaSeq:
        #if char == '-':  #was used to make sure the unicode dashes were being converted correctly
            #print "watchout"
            #pass
        if char != '-':
            for resNum in csaRes_dict:
                if int(csaResCounter) == (int(resNum) -1):
                    #make sure that the CSA res (as applied to the PDB by the resnum)
                    #actually matches with CSA res from CSA
                    if str(char) == csaResDict[resNum]:
                        print 'true'
                        print 'resNum'
                        print resNum
                        print 'char number'
                        print csaResCounter
                        print '  csaStart'
                        print csaStart
                        print 'csa seq'
                        print csaSeq
                        print 'char'
                        print char
                        print 'csaResDict'
                        print csaResDict[resNum]
                    else:
                        print 'false'
                        print 'resNum'
                        print resNum
                        print 'char number'
                        print csaResCounter
                        print '  csaStart'
                        print csaStart
                        print 'csa seq'
                        print csaSeq
                        print 'char'
                        print char
                        print 'csaResDict'
                        print csaResDict[resNum]
                        
                    #print '  resNum'
                    ##print resNum
                    ##print 'csaResCounter'
                    ##print csaResCounter
                    ##print "little success"
                    #print char + '   char'
                    #print str(pdbSeq[posCounter]) + '   str(pdbSeq[posCounter]) '
                    if str(char) == str(pdbSeq[posCounter]):
                        #print str(pdbResCounter) + 'pdbResCounter'
                        #print str(posCounter) + '  posCounter'
                        #print 'great success!'
                        results[pdbResCounter] = char
            csaResCounter +=1
        if pdbSeq[posCounter] != '-':
            pdbResCounter += 1
        posCounter += 1
    #print 'new'
    for i in results:
        #print i
        pass
    return results



#open the .xml file made by blastp using biopython blast_record module
result_handle = open(str(csaInput) + ".xml")
blast_record = NCBIXML.read(result_handle)
#if E-value < .001
E_VALUE_THRESH = 100 #change this when actually running it
counter = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            counter += 1
#call hitChecker(csaRes_list,csaSeq,csaStart,pdbSeq,pdbStart)     
            pdbHits = hitChecker(csaResDict,hsp.query,hsp.query_start,hsp.sbjct,hsp.sbjct_start)
            for i in pdbHits:
                #print i
                pass
            
            #print 'test'    
            #print pdbHits
#            for resnum in site_list:
#                queryRes = hsp.query
            
            ##print('****Alignment****')
            ##print('sequence:', alignment.title)
            ##print('length:', alignment.length)
            ##print('e value:', hsp.expect)
            ##print(hsp.query[0:] + '...')
            ##print(hsp.match[0:75] + '...')
            ##print(hsp.sbjct[0:75] + '...')
            ##print(dir(hsp))
if counter == 0:
    sys.exit()


   


#identify hits with significant E-values (.001)

#retrieve these sequences using fastacmd (blast), add the query to create ?hits.query.fasta?

#build HMM using SAM for the query (modelfromalign) -> query.mod

#align the PDB hits to the HMM using align2model (align2model query.hits -i query.mod -db hits.query.fasta); this will produce query.hits.a2m