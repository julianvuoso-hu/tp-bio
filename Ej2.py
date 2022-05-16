from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

fasta_string = open("out1/sequence-1.fasta").read()
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string)
blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.04

output = ""
for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
          if hsp.expect < E_VALUE_THRESH:
               output += "****Alignment****\n"
               output += "sequence: %s\n" % alignment.hit_def.split(' >')[0]
               output += "accession: %s\n" % alignment.hit_id.split('|')[1]
               output += "length: %d\n" % alignment.length
               output += "score: %s\n" % str(hsp.score)
               output += "identity: %d/%d(%.2f%%)\n" % (hsp.identities, hsp.align_length, (100 * hsp.identities / hsp.align_length))
               output += "E-value: %f\n" % hsp.expect
               output += "query: %s\n" % hsp.query
               output += "match: %s\n" % hsp.match
               output += "sbjct: %s\n\n" % hsp.sbjct

f = open("out2/blast.out", "w")
f.write(output)
f.close()