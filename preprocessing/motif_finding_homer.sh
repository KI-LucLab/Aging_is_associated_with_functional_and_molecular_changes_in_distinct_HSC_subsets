#! /bin/bash -l

#SBATCH -A naiss2023-22-550
#SBATCH -o Motif_finding.out
#SBATCH -e Motif_finding.err
#SBATCH -J Motif_finding.job
#SBATCH -n 1
#SBATCH -t 24:00:00

# LOAD REQUIRED MODULES
module load bioinfo-tools
module load HOMER/4.11

# GET USER ARGUMENTS
while getopts 'f:' flag; do
  case "${flag}" in
    f) inputfile="${OPTARG}" ;;
    *) echo "Unexpected option ${flag}" ;;
  esac
done

inputpath="/proj/snic2020-16-175/nobackup/processedData/atacseq/motif_enrichment/input/" # loaction of input files
outputpath="/proj/snic2020-16-175/nobackup/processedData/atacseq/motif_enrichment/output/" # location of output files

mkdir $outputpath/${inputfile%.bed}

findMotifsGenome.pl $inputpath/$inputfile mm10 $outputpath/${inputfile%.bed} -size given
