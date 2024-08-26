#! /bin/bash -l

#SBATCH -A naiss2023-22-550
#SBATCH -o bamToTagdir.out
#SBATCH -e bamToTagdir.err
#SBATCH -J bamToTagdir.job
#SBATCH -n 1
#SBATCH -t 24:00:00

# GET USER ARGUMENTS
while getopts 'f:' flag; do
  case "${flag}" in
    f) fileSamples="${OPTARG}" ;;
    *) echo "Unexpected option ${flag}" ;;
  esac
done

# LOAD REQUIRED MODULES
module load bioinfo-tools
module load BEDTools
module load bamtools
module load samtools
module load HOMER

# Location of input files
REFERENCE="/proj/snic2020-16-175/genomes/mouse.mm10.genome"
run_path="/proj/snic2020-16-175/nobackup/processedData/atacseq/TagDirs_shifted_from_bam/"
fileSamplesPath="/proj/snic2020-16-175/analysisFileLists"

########## READ INPUT FILES FROM $fileSamples ##########
readarray -t filesToAnalyse < $fileSamplesPath/$fileSamples # read $fileSamples file line by line, remove newlines (-t) and read into array (filesToAnalyse)
if [ ! -f $fileSamplesPath/$fileSamples ]; then
  echo "ERROR: Can't find the input file list, run terminated\!"
  exit
fi

# Go to the run path
cd ${run_path}

for input_file_whole in ${filesToAnalyse[@]}; do

	input_file=$(basename "${input_file_whole}")
	cp -r ${input_file_whole} ${run_path}/${input_file}

# Shift reads
	# 1. convert bam to bed
	# 2. shift by +4bp on the "+" strand and by -5bp on the "-" strand. Total shift has to be 9bp.
	# 3. convert bed to bam
	# 4. sort bam

	bedtools bamtobed -i $input_file | \
	awk -v OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4,$4,$5,$6}else if($6=="-" && $2 > 4){print $1,$2-5,$3-5,$4,$5,$6}}' > ${input_file}.shifted.bed
	bedtools bedtobam -i ${input_file}.shifted.bed -g ${REFERENCE} | \
	bamtools sort -out ${input_file}.shifted.sorted.bam

# Convert to tagdir
	tagDirName=${input_file}.shifted.sorted.bam_TagDir
	samtools view -h ${input_file}.shifted.sorted.bam > ${input_file}.shifted.sorted.sam
	file_sam=${input_file}.shifted.sorted.sam
	makeTagDirectory ${tagDirName} $file_sam -format sam
	rm $input_file

done
