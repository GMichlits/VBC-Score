INSTRUCTIONS for use of generate_VBCscore scripts:

current version:
for python2: generate_VBCscore_v0.0_py2.py
for python3: generate_VBCscore_v0.0_py3.py

execution:
copy the files 
generate_VBCscore_v0.0_py3.py
generate_VBCscore_v0.0_py2.py
demo_input.tsv
to a folder on your computer
Navigate to the folder
execute the script (e.g. on mac in terminal by typing “python generate_VBCscore_v0.0_py2.py”)
an input prompt will open to give the path of the input tsv file (that file must be in the same folder). 
Just pressing enter will automatically load the demo_input.tsv file.

If you have any queries do not hestitate to contact me. georg.michlits@gmail.com

################################################################
Instructions for input file:
needs to be a tab separated value file with 18 columns the first row will be skipped (contains table titles).
See demo_input.tsv as reference.

column [0] - sgRNA_ID: optional for your own reference

column [1] - ntseq: nucleotide sequence, optional to be used for your own reference. (doench2016 and tracrv2 sgRNA activity scores require 30bp to calculate on-target activity scores.  4bp - 20bp + NGG + 3bp = 30bp
If you do not have tracrv2 scores available you can leave column 18 empty (it needs to be present but empty) and tracrv2 scores
will be calculated from the 30bp nucleotide sequence in column[1] not that for the script to work it must be 30bp and must contain
GG at positions 25 and 26.

column [2]  - phylo21: The average phylo score of a window of 21 bp around the CRISPR cut site. the CRISPR cut site is defined as 17th basepair of the sgRNA spaver sequence. 
e.g. for spacer: CTGGGCAGAGTTCTGGAAATT CAGTGGAAT the expected cut site is 
CTGGGCAGAGTTCTGGAAATT / CAGTGGAAT so the nucleotide C is considered the “cut site”, The 10bp upstream and downstream are added to give the 21bp range. Phylogenetic scores can be accessed via the UCSC genome browser.

column [3]  - phast21: same as for phylo21. 

column [4] - Pfam_domain: If the cut site is in an annotated Pfam domain the feature is 1 else 0.

column [5] - Pfam_family: If the cut site is in an annotated Pfam family the feature is 1 else 0.

column [6] - distance: is the relative position (ranging from 0 to 1) of the cute site along the mRNA sequence (5’ to 3’ oritenation with 0 being the start  codon and 1 the stop codon).

column [7] - mod3: modal 3. If the exon length is a multiple of 3 the value is 1 else the value is 0.

column [8] - exlen: the exon length in [bp]

column [9] - SD_distance: The distance to the Splice donor (exon boundary) in bp. (the script automatically calculates with a value of 21 if the splice donor is further away than 21 bp.)

column [10] - SA_distance: The distance to the Splice acceptor (exon boundary) in bp. (the script automatically calculates with a value of 21 if the splice donor is further away than 21 bp.)

column [11] - AAcons7: If a score for phylogenetic conservation that is based on amino acid alignments. The score is to a large extend redundant to phylo21 and phast21 but has shown higher predictive power in our hands. Since derivation of the score can be challenging (for an instruction see Online methods “Generation of aa conservation score” of associated publication) it is also possible to calculate VBC-Scores without AAcons7 input, in that case the median AAcons7 value of sgRNAs will be used for the purpose of calculating a VBC-Score.

column [12] - AA_before_cut. The 6 amino acids upstream of the CRISPR cut site. (or less amino acids if there are for example less aa upstream because the cut is close to the start codon). For the codon “CAT” in the sequence NNNCATNNN the cuts N / C, C / A and C / T belong to the amino acid histidin codon “CAT”. A CRISPR Cut between T / N would already refer to a cut “at” the next amino acid. Use the one letter amino acid code. (sequence N to C terminal - or 5’ to 3’)

column [13] - AA_at_cut. see definition before (always 1 amino acid) (sequence N to C terminal - or 5’ to 3’)

column [14] - AA_after_cut. see definition before give 6 amino acids at most (sequence N to C terminal - or 5’ to 3’)

column[15] - Framshift_ratio_inDelphi: The predicted ratio of frameshifts generated by a CRISPR cut using Cas9 at this site. Use prediction tool inDelphi. (rang 0 to 1)

column [16] - Doench 2016: The Doench2016 sgRNA score “Rule Set 2”.

column [17] - tracrv2: The tracrv2 sgRNA score - can be derived using the table presented in the associated publication. Alternatively this column can be left empty. In that case tracrv2 score will be calculated from ntseq (if it is available in column 2).

output:
The script will add 3 columns to the end of the manuscript containing
(1) combined sgRNA score (this is roughly a 50:50 merge of the Doench2016 score and the tracrv2 score.)
(2) Bioscore (this is a merge of conservation, Pfam, exon structure and amino acid identity)
(3) VBC (this merges the frameshift-prediction with combined sgRNA score and Bioscore)

If you have any queries do not hestitate to contact me. 
georg.michlits@gmail.com
