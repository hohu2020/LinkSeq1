#!/usr/bin/bash
# 22/01, 2021
# hongxing


# This code is to:
#   - extract the consensus BCR sequences from TRUST4 output file: "*_barcode_airr.tsv"
#   -


###########################################################################
#  1. extract seqence_id and sequence from "*_barcode_airr.tsv" for ChangeO-Igblast
############################################################################

# "file1" "file2"  "file3" "fileN" are dummy file names which should be replaced with YOUR file names

for fileName in  "file1" "file2"  "file3" "fileN"
do
  # input
  main_inputPath="fill_with_TRUST4_output_path"

  inputFile="${main_inputPath}${fileName}HL_barcode_airr.tsv"
  # output
  outputPath="fill_with_your_output_path"
  mkdir -p $outputPath

  # extract seqence_id and sequence  for igblast
  awk '{ print ">"$1 "\n" $2}' $inputFile \
  > "${outputPath}${fileName}_noPear_barcode_airr.fa"


done


#================================================#
# example in my laptop:
#================================================#

for fileName in  "HD1_stD2"
do
  # input
  main_inputPath="/home/hohu/Desktop/LinSeq/result/summary/no1/"

  inputFile="${main_inputPath}${fileName}HL_barcode_airr.tsv"
  # output
  outputPath="/home/hohu/Desktop/LinSeq/result/summary/no2/"
  mkdir -p $outputPath

  # extract seqence_id and sequence  for igblast
  awk '{ print ">"$1 "\n" $2}' $inputFile \
  > "${outputPath}${fileName}_noPear_barcode_airr.fa"


done
