#!/usr/bin/bash
# 22/01, 2021
# hongxing



#############
#  it runs in docker
# This code is to:
#   - run IgBlast in ChangeO
#   - align bcr to Germline sequence
#############

# activate docker and mount "Linkseq" folder to the data dir in docker

 sudo docker run -it -v /home/hohu/Desktop:/data:z immcantation/suite:4.3.0 bash
# password

######################################################################################
# IgBlast in immcantation docker
# NOTE : mannually  copy paste the code bellow in docker bash after logging into Docker
#######################################################################################

# Directory containing IMGT-gapped reference germlines. Defaults to /usr/local/share/germlines/imgt/[species name]/vdj.
# no PEAR use "*_noPear_barcode_airr.fa"  file
# already collapsed at cell barcode level


# "file1" "file2"  "file3" "fileN" are dummy file names which should be replaced with YOUR file names

for fileName in "file1" "file2"  "file3" "fileN"
do
  # input of no3 is from the output of no2
  main_inputPath="fill_with_your_output_path_from_no2"

  # use the "_noPear_barcode_airr.fa"
   READS="${main_inputPath}${fileName}_noPear_barcode_airr.fa"


  # output
  OUT_DIR="fill_with_your_output_path_for_no3"
  mkdir -p $OUT_DIR

  SAMPLE_NAME="${fileName}_noPearIgblastBarcodeAirr"
  NPROC="8"
  #-------------------------------------------------------------
  cd / # / (back to the root level)

  changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

  # Define clone
  cd $OUT_DIR
  #-------------------------------------------------------------
  inputFileName1="${fileName}_noPearIgblastBarcodeAirr_db-pass.tsv"
  #-------------------------------------------------------------
  DefineClones.py -d $inputFileName1 --act set --model ham \
      --norm len --dist 0.16

  # repair clone
  #-------------------------------------------------------------
  inputFileName2="${fileName}_noPearIgblastBarcodeAirr_db-pass_clone-pass.tsv"
  #-------------------------------------------------------------
  CreateGermlines.py -d $inputFileName2 -g dmask --cloned \
          -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
          /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
          /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta

done






#================================================#
# example in my laptop:
#================================================#

# activate docker and mount "Linkseq" folder to the data dir in docker

 sudo docker run -it -v /home/hohu/Desktop:/data:z immcantation/suite:4.3.0 bash
# password

######################################################################################
# IgBlast in immcantation docker
# NOTE : mannually  copy paste the code bellow in docker bash after logging into Docker
#######################################################################################

for fileName in "HD1_stD2"
do
  # input
  main_inputPath="/data/LinSeq/result/summary/no2/"

  # use the "_noPear_barcode_airr.fa" instead of "annot.fa"
   READS="${main_inputPath}${fileName}_noPear_barcode_airr.fa"


  # output
  OUT_DIR="/data/LinSeq/result/summary/no3/"
  mkdir -p $OUT_DIR

  SAMPLE_NAME="${fileName}_noPearIgblastBarcodeAirr"
  NPROC="8"
  #-------------------------------------------------------------
  cd / # / (back to the root level)

  changeo-igblast -s $READS -n $SAMPLE_NAME -o $OUT_DIR -p $NPROC

  # Define clone
  cd $OUT_DIR
  #-------------------------------------------------------------
  inputFileName1="${fileName}_noPearIgblastBarcodeAirr_db-pass.tsv"
  #-------------------------------------------------------------
  DefineClones.py -d $inputFileName1 --act set --model ham \
      --norm len --dist 0.16

  # repair clone
  #-------------------------------------------------------------
  inputFileName2="${fileName}_noPearIgblastBarcodeAirr_db-pass_clone-pass.tsv"
  #-------------------------------------------------------------
  CreateGermlines.py -d $inputFileName2 -g dmask --cloned \
          -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
          /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
          /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta

done
