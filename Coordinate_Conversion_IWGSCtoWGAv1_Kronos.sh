#!/bin/sh
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=15000
#SBATCH -J BLAST_database_v1.0tetraploid
#SBATCH -o BLAST_database_v1.0tetraploid.%N.%j.out
#SBATCH -e BLAST_database_v1.0tetraploid.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.harrington@jic.ac.uk

#Conversion of coordinates for the original MAPS-derived SNPs (in this case from line K2282)
#the MAPS SNP file for each line can be downloaded from www.wheat-tilling.com

# To do this, we want to BLAST the SNP, including 150 bp of flanking sequence on either side, against the NRGene assembly. This contig size should give us enough confidence in the
# position of the alignment and yet be small enough to avoid InDels that might disrupt the alignment.

source blast+-2.2.30
source bedtools-2.24.0
source samtools-1.3

Input='K2282_maps_main_set.txt'
filename='K2282_maps_main_set_KronosRef'

#Make a directory to hold the coordinate conversion, and change into that directory
cd /nbi/group-data/ifs/NBI/Cristobal-Uauy/Sophie/ENQ2022/CoordinateConversion

#Construct a BED file with all SNP positions
#Extract first and second column, with Chrom/Scaffold and position; exclude first row for header
sed -n '2,$ p' $Input | cut -f 1 /dev/stdin > $filename-Chrom.tab
sed -n '2,$ p' $Input | cut -f 2 /dev/stdin > $filename-Pos.tab

#create a third fle with every entry changed to the '_' symbol, same number of lines
sed 's/[0-9]*/_/ g' $filename-Pos.tab > $filename-underscore.tab

#Combine the three files together to give each SNP a unique ID (-d '\0' means no delimiter included when pasting)
paste -d '\0' \
		$filename-Chrom.tab \
		$filename-underscore.tab \
		$filename-Pos.tab \
		> $filename-uniqueID.tab


#Now, to make coordinates base 0, need to subtract 1 from each entry, by creating one file where each line entry is '1', then combine with the position file
sed 's/[0-9]*/1/ g' $filename-Pos.tab > $filename-1.tab
paste \
		$filename-Pos.tab \
		$filename-1.tab \
		> $filename-Pos-1.tab

#This allows the two columns to be subtracted from each other to give the new, base 0 position
awk '{ $3 = $1 - $2; print $3 }' $filename-Pos-1.tab > $filename-substracted.tab

#Now combine all of these files to make the correct BED file
paste \
		$filename-Chrom.tab \
		$filename-substracted.tab \
		$filename-Pos.tab \
		$filename-uniqueID.tab \
		> $filename.bed


#Next step: Extend coordinates to 150bp either side of the SNP, but use fasta index to prevent extension past the start/end of scaffold.
#Use the command "slop" from BEDTools --> might be better just to have 300 bo

bedtools slop \
		-i $filename.bed \
		-g /nbi/group-data/ifs/NBI/Cristobal-Uauy/IWGSC_data/IWGSC_all/IWGSC-all_UCW-Kronos_U.fa.fai \
		-l 300 \
		-r 0 \
		> $filename-slop-left-300.bed

#Using the modified BED file, extract SNP and flanking sequence using "getfasta" in BEDTools; -name give the fasta file the unique ID from the BED file
bedtools getfasta \
		-fi /nbi/group-data/ifs/NBI/Cristobal-Uauy/IWGSC_data/IWGSC_all/IWGSC-all_UCW-Kronos_U.fa \
		-bed $filename-slop-left-300.bed \
		-fo $filename-slop-left-300.fa \
		-name

#Finally, index the fasta file to get the length and name of each entry
samtools faidx $filename-slop-left-300.fa

#Preparation is complete- now time to BLAST against the NRGene assembly

#To prevent odd mismatching that has been seen when BLASTing large numbers of fasta files at once, separate each FASTA entry into a separate file (split every other lines to get the two lines of the fasta entries out)
mkdir -p split-files-300
split -l 2 --suffix-length=4 --additional-suffix=.txt $filename-slop-left-300.fa ./split-files-300/$filename-slop-left-300_split-

#Now cycle through all of the fasta files and BLAST against the NRGene, creating a separate output file for each. Use stringent parameters (single hit and single HSP)

cd ./split-files-300/
mkdir ./split-BLAST-output-slop-left-300-Kronos

for file in $filename-slop-left-300_split-*; do
	blastn \
		-task megablast \
		-db /nbi/group-data/ifs/NBI/Cristobal-Uauy/WGAv1.0/161010_Chinese_Spring_v1.0_pseudomolecules_parts_tetraploid \
		-query $file \
		-max_target_seqs 1 \
		-max_hsps 1 \
		-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop" \
		-parse_deflines \
		-out ./split-BLAST-output-slop-left-300-Kronos/WGAv1.0_parts-megablast-$file
	done

#Append all of the output files together; can't use cat due to number of files opened, so use for-loop to append each file one at a time
cd ./split-BLAST-output-slop-left-300-Kronos
for file in WGAv1.0_parts-megablast-$filename-slop-left-300_split-*; do \
	cat $file >> /nbi/group-data/ifs/NBI/Cristobal-Uauy/Sophie/ENQ2022/CoordinateConversion/$filename-slop-left-300_vs_WGAv1.0_parts_tetraploid-megablast.txt \
	; done

#for neatness, remove empty spaces from Bitscore column
sed 's/  // g' \
		/nbi/group-data/ifs/NBI/Cristobal-Uauy/Sophie/ENQ2022/CoordinateConversion/$filename-slop-left-300_vs_WGAv1.0_parts_tetraploid-megablast.txt \
		> /nbi/group-data/ifs/NBI/Cristobal-Uauy/Sophie/ENQ2022/CoordinateConversion/$filename-slop-left-300_vs_WGAv1.0_parts_tetraploid-megablast-sed.txt

#Now convert coordinates -> go into Excel and use "VLOOKUP" options as detailed below to finalise the data:

## 1) Add columns for the positions flanking the SNP (i.e. the left and right positions flanking the BLAST result; should be 150bp or less for each, based on the BLAST requirements)
## 		To do this, add the left and right positions from the -150.bed file created earlier as the "left position" and "right position" columns. Then subtract the "left position" from the SNP position to get "bp left of SNP", and similarly to obtain the right side distance.

## 2) Subtract the send and sstart columns (from the BLAST output) to determine the orientation of the BLAST sequence (i.e. a negative number is a (-) orientation region)

## 3) Get the correct SNP position for the RefSeq v1.0 assembly 
##		3a) if the sequence is in the + orientation (See #2 above), use the following formula: "subject-start" + "bp left of SNP" - "query-start" + 1
##		3b) if the sequence is in the - orientation (see #2 above), use the following formula: "subject-start" - "bp left of SNP" + "query-start" - 1
#			3.1)

## other columns can also be added to the file as desired from i.e. the original SNP set (such as the WT and Mut allele etc.)

