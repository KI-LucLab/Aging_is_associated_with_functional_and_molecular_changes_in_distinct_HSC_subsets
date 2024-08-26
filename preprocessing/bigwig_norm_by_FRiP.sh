#! /bin/bash -l

#SBATCH -A naiss2023-22-550
#SBATCH -o bigwig_norm_FRiP.out
#SBATCH -e bigwig_norm_FRiP.err
#SBATCH -J bigwig_norm_FRiP.job
#SBATCH -n 1
#SBATCH -t 24:00:00

# LOAD REQUIRED MODULES
module load bioinfo-tools
module load ucsc-utilities

# GET USER ARGUMENTS
while getopts 'f:' flag; do
  case "${flag}" in
    f) fileSamples="${OPTARG}" ;;
    *) echo "Unexpected option ${flag}" ;;
  esac
done

########## CHANGE THE FOLLOWING IF NEEDED #######
scriptPath="/proj/snic2020-16-175/scripts"	# location of this script, should be existing
fileSamplesPath="/proj/snic2020-16-175/analysisFileLists"	# location of fileSamples, should be existing
outputPath="/proj/snic2020-16-175/webexport/mm10_atacseq_nfcore_pipeline_bigwigs_FRiP_normalised/" # location of result files
chromsizes="/proj/snic2020-16-175//genomes/mouse.mm10.genome"
#################################################

########## MAKE output dir if needed ##########
mkdir -p $outputPath

########## READ INPUT FILES FROM $fileSamples ##########
readarray -t filesToAnalyse < $fileSamplesPath/$fileSamples # read $fileSamples file line by line, remove newlines (-t) and read into array (filesToAnalyse)
if [ ! -f $fileSamplesPath/$fileSamples ]; then
  echo "ERROR: Can't find the input file list, run terminated!"
  exit
fi

########## PROCESS INPUT FILE LIST ##########
bigwigList=("none")
FRiPList=("none")
i=0 # line index
for line in "${filesToAnalyse[@]}"; do
  IFS=' ' read -r -a splitLine <<< "${line}" # split line on space
  bigwigList[i]=${splitLine[0]}
  FRiPList[i]=${splitLine[1]}
  i=$((i+1))
done

# Normalise by FRiP
for (( i=0; i<${#bigwigList[@]}; i++ )); do # for each tagdir in list, find peaks
  bigwig_basename="$(basename "${bigwigList[i]}")"
  bigWigToBedGraph "${bigwigList[i]}" ${outputPath}/${bigwig_basename%.bigWig}.bedGraph
  awk -v var=${FRiPList[i]} '{print $1"\t"$2"\t"$3"\t"($4/var)}' ${outputPath}/${bigwig_basename%.bigWig}.bedGraph > ${outputPath}/${bigwig_basename%.bigWig}.norm.bedGraph
  bedGraphToBigWig ${outputPath}/${bigwig_basename%.bigWig}.norm.bedGraph $chromsizes ${outputPath}/${bigwig_basename%.bigWig}.norm.bigWig
done

wait

# Remove .bedGraph files
for (( i=0; i<${#bigwigList[@]}; i++ )); do # for each tagdir in list, find peaks
  bigwig_basename="$(basename "${bigwigList[i]}")"
  rm ${outputPath}/${bigwig_basename%.bigWig}.bedGraph
  rm ${outputPath}/${bigwig_basename%.bigWig}.norm.bedGraph
done
