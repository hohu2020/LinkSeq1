#!/usr/bin/bash
# 22/01, 2021
# hongxing
# contact: hhxing@gmail.com

# BCR fastq file can be downloaded from ENA with accession number: PRJEB57285
# fill the contents in ""
cd "fill_with_the_path_to_TRUST4_dir"

./run-trust4 \
-t 7 \
-f hg38_bcrtcr.fa \
--ref human_IMGT+C.fa \
--od "output_dir" \
-u "fill_with_the_BCR_read-1_file"  \
--barcode "fill_with_the_BCR_Index_file(with cell barcode and UMI sequences)" \
--barcodeRange 0 11 +  \
--UMI "fill_with_the_BCR_Index_file(with cell barcode and UMI sequences)"  \
--umiRange 12 19 + \
-o  "prefix_to_output_files"


#================================================#
# example in my laptop:
#================================================#

cd /home/hohu/TRUST4

./run-trust4 \
-t 7 \
-f hg38_bcrtcr.fa \
--ref human_IMGT+C.fa \
--od "/home/hohu/Desktop/LinSeq/result/summary/no1/" \
-u "/home/hohu/Desktop/LinSeq/data/000000000-J95WL_HuAbSeq_02_20s003655-1-1_Hu_lane1EM960HL_2_sequence.txt.gz"  \
--barcode "/home/hohu/Desktop/LinSeq/data/000000000-J95WL_HuAbSeq_02_20s003655-1-1_Hu_lane1EM960HL_1_sequence.txt.gz" \
--barcodeRange 0 11 +  \
--UMI "/home/hohu/Desktop/LinSeq/data/000000000-J95WL_HuAbSeq_02_20s003655-1-1_Hu_lane1EM960HL_1_sequence.txt.gz"  \
--umiRange 12 19 + \
-o  "HD1_stD2HL"
