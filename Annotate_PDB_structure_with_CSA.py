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
from time import sleep


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
csaInput = None;
#parse input (example input should look like 12as|A)

if len(csaInputRaw) == 4: #if no chain ID is given
    csaInput = csaInputRaw
else: #a chain ID is given
    csaInputLHS, csaInputRHS = csaInputRaw.split("|")
    csaInput = csaInputLHS +"."+csaInputRHS #csaInput looks like 12as.A
    #print csaInputLHS +"."+csaInputRHS + '   csaInputLHS . csaInputRHS'

#assign variables for read_sequence_from_pdb_file:
#checks PDB structure files located at:
pdbLocation = "/clusterfs/ohana/external/pdb/"
pdbId = csaInputLHS
chainId = csaInputRHS

#get PDB sequence using read_sequence_from_pdb_file; this script was written by Cyrus Afrasiabi
#and extracts the sequence and residue numbers from a PDB structure file, thus giving the sequence from structure
sequence, numberedResidues = read_sequence_from_pdb_file(pdbId, chainId, pdbLocation)
numberedResiduesDict = dict(numberedResidues)
#print sequence
#print numberedResidues
#print numberedResiduesDict


#parse pdb sequence from structure to a .fasta file
n = 80 #how long of lines to make
sequenceFasta = [sequence[i:i+n] for i in range(0, len(sequence), n)]
#write to a pdb.fasta file
file = open("/home/tnrich/blastp_test_dir/%s.fasta" % (csaInput),'w') 
try:
    file.write('>pdb|%s mol:protein length:%s\n' % (csaInputRaw,str(len(sequence))))
    #for each n length segment of sequenceFasta, append it to the pdb.chain.fasta file with a new line
    for sequencePart in sequenceFasta:
        file.writelines(sequencePart + '\n')
finally:
    file.close()
    

#This is the outdated and incorrect way of obtaining the PDB sequence:
#print csaInput
#command = "blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry ' " + str(csaInputRaw) + "' > "+ str(csaInput) + '.fasta'
        #eg: blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry 12as > 12as.fasta
        #eg: blastdbcmd -db /clusterfs/ohana/external/pdb/blastdbs/pdb -entry "12as|A" > 12as.A.fasta
#print command
#os.system(command)


#for each protein in csa:
#Run blastp on the pdb.chain.fasta file we just created
#command = "cat /home/tnrich/blastp_test_dir/%s.fasta" % (csaInput)
#print command
#print 'here it comes'
#os.system(command)
command = "blastp -db /clusterfs/ohana/external/pdb/blastdbs/pdb -query /home/tnrich/blastp_test_dir/%s.fasta -out /home/tnrich/blastp_test_dir/%s.xml -outfmt 5" % (str(csaInput), str(csaInput))
#print command
os.system(command)
#os.system(command)
#command = 'cat /home/tnrich/blastp_test_dir/%s.xml' %(csaInput)
#print command
#os.system(command)


#translation table:
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
def hit_checker(csaRes_dict,csaSeq,csaStart,csaNumberedResidues,csaNumberedResiduesDict,pdbSeq,pdbStart):
    
    #make sure that the csaRes is in the hsp
    
    csaSeq = unicodedata.normalize('NFKD', csaSeq).encode('ascii','ignore') #convert the unicode csaSeq from blastp to ascii string
    pdbSeq = unicodedata.normalize('NFKD', pdbSeq).encode('ascii','ignore')
    posCounter = 0;
    csaResCounter = int(csaStart) - 1 
    pdbResCounter = int(pdbStart) 
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
    #Step 1:
    #we want to parse through every char in the CSA seq and first check if it is one of the catalytic residues in
    #the CSA list csaRes_Dict. Every time a characer is a -, we discard it. If it is a character, we increase the
    #csaResCounter by 1 (which starts at csaStart-1 (for example if the csaStart=1, our first time through, our csaResCounter would =0 ))
    #This allows us to know the number of the character we're on in the csaSeq, which we want to compare to the CSA site list
    #(do this by comparing csaNumberedResidues[csaResCounter][0] == resNum in csaRes_dict )
    #for example, if resnum = 1 (aka the first residue in the sequence is catalytic according to CSA),
    #then csaNumberedResidues[csaResCounter][0] should be 1 also.
    #it is important to pull this number from the csaNumberedResidues because there can be breaks in the sequence that the csaResCounter alone will not account for
    
    #Step 2:
    #then check if the residue is the same as the one in the pdb sequence
    #eg if str(char) == str(pdbSeq[posCounter]):
    #posCounter starts at 0 and is updated every char, thus by the 30th time through the char loop, poscounter = 29
    #so pdbSeq[29] is asking about the 30th character in the pdb seq which, if it is aligned to the catalytic char in csaSeq,
    #should make the above relationship hold
    
    #Step 3
    #if the pdb residue and the CSA residues are aligned, we then want to transfer that annotation to the correct residue number of the pdbseq
    #we can extract that number using pdbResCounter which should be set to the residue number (it starts at pdbStart and increments each time its char!='-')
    #eg pdbseq residue 1 is a match, pdbResCounter will also equal 1
    
    for char in csaSeq:     
        #if char == '-':  #was used to make sure the unicode dashes were being converted correctly
            #print "watchout"
            #pass
        if char != '-':
            #test that the numbering is correct:
            #CUR
            #if csaNumberedResidues[csaResCounter][1] == char: #This test was passed!!
            #    print "true"
            #else:
            #    print "falso!"
            
                
            for resNum in csaRes_dict:
                if csaNumberedResidues[csaResCounter][0] == resNum:
                        
                #if int(csaResCounter) == (int(resNum) -1):  
                    #make sure that the CSA res (as applied to the PDB by the resnum)
                    #actually matches with CSA res from CSA
                    #if char == csaNumberedResidues[csaResCounter][1]: #make sure char and csaNumberedResidues are matching up
                    #    print 'char == csaNumberedResidues[csaResCounter][1]'
                    #else:
                    #    print "no char =="
                    #if str(char) == csaResDict[resNum]: #test
                    #    print 'true'
                    #    print 'resNum'
                    #    print resNum
                    #    print 'char number'
                    #    print csaResCounter
                    #    print '  csaStart'
                    #    print csaStart
                    #    print 'csa seq'
                    #    print csaSeq
                    #    print 'char'
                    #    print char
                    #    print 'csaResDict'
                    #    print csaResDict[resNum]
                    #else:
                    #    print 'false'
                    #    print 'resNum'
                    #    print resNum
                    #    print 'char number'
                    #    print csaResCounter
                    #    print '  csaStart'
                    #    print csaStart
                    #    print 'csa seq'
                    #    print csaSeq
                    #    print 'char'
                    #    print char
                    #    print 'csaResDict'
                    #    print csaResDict[resNum]
                        
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

#TODO
def populate_database(): #
    pass



#open the .xml file made by blastp using biopython blast_record module
result_handle = open("/home/tnrich/blastp_test_dir/%s.xml" % (str(csaInput)))
blast_record = NCBIXML.read(result_handle)
#if E-value < .001
E_VALUE_THRESH = 100 #change this when actually running it
counter = 0

##figure out what is in blast_record
#for item in blast_record.__dict__.items():
        #print item

for alignment in blast_record.alignments:
    ##figure out what is in alignment
    #for item in alignment.__dict__.items():
    #    print item
    
    
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            counter += 1
#call hitChecker(csaRes_list,csaSeq,csaStart,pdbSeq,pdbStart)     
            pdbHits = hit_checker(csaResDict,hsp.query,hsp.query_start,numberedResidues,numberedResiduesDict,hsp.sbjct,hsp.sbjct_start)
            
            if pdbHits != {}:
                #for each non-empty hit:
                print '----------------------------------------------------------'
                #print query info
                print 'Query:'
                print vars(blast_record)['query']
                print 'Annotation Type:'
                print 'CsaToPdbBlastP'
                #print alignment info:
                print 'Hit:'
                print vars(alignment)['title']
                print 'Hit Length:'
                print vars(alignment)['length']
                #print hit info:
                print 'Residue Hits:'
                print pdbHits
                for item in hsp.__dict__.items():
                    print item
                
            #
            #for key, value in pdbHits:
            #    catalyticResidueIdentifier = {}
            #
            #
            #print 'pdbHits'
            #print pdbHits
            #print "end"
            #print(dir(hsp))
            #pdbHits.update(hsp.__dict__.items())
            #
            ##get e value:
            #print hsp.expect
            #
            ##get Bit score:
            #print hsp.score
            
            
                
                
            
            
            
        
            #Here's a list of the attributes we want to send to our database table:
            #identifier for the residue being annotated (e.g., D185)
            #EvidenceCode (alphanumeric): foreign key to another table where we explain what these evidence codes mean; BCSA (BLAST from CSA), (BSP: BLAST from SwissProt).
            #Source: UniProt sequence accession/identifier or PDB identifier; 
            #Bit score
            #E-value
            #Percent identity
            

            
            
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
            
if counter == 0:
    sys.exit()


   


#identify hits with significant E-values (.001)

#retrieve these sequences using fastacmd (blast), add the query to create ?hits.query.fasta?

#build HMM using SAM for the query (modelfromalign) -> query.mod

#align the PDB hits to the HMM using align2model (align2model query.hits -i query.mod -db hits.query.fasta); this will produce query.hits.a2m