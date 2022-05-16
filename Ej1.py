from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


input_seq = "sequence.gb"
record = SeqIO.read(input_seq, "genbank")
table = 1
min_pro_len = 100

i = 1
for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
     for frame in range(3):
         length = 3 * ((len(record)-frame) // 3) #Multiple of three
         for pro in nuc[frame:frame+length].translate(table).split("*"):
             if len(pro) >= min_pro_len:
                 print("%s...%s - length %i, strand %i, frame %i" \
                       % (pro[:30], pro[-3:], len(pro), strand, frame))
                 orf_seq = SeqRecord(pro)
                 orf_seq.id = "lcl"
                 orf_seq.description = "ORF%d" % i
                 SeqIO.write(orf_seq, "out1/%s-%d.fasta" % (input_seq[:-3], i), 'fasta')
                 i += 1
