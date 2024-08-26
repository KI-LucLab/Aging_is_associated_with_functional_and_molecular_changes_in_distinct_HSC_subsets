#!/bin/bash -l

#######  
## There is no need to submit a job for this: do not run "sbatch"
## To run: 
## bash nfcore.sh -f [file_list.csv] -r [raw_data_folder]

# SET DEFAULT ARGUMENTS
resume=false
# GET USER ARGUMENTS
while getopts 'f:r:x:' flag; do
  case "${flag}" in
    f) fileSamples="${OPTARG}" ;;
    r) rawDataFolder="${OPTARG}" ;;
    x) resume="${OPTARG}" ;;
    *) echo "Unexpected option ${flag}" ;;
  esac
done

########## CHANGE THE FOLLOWING IF NEEDED #######
scriptPath="/proj/snic2020-16-175/scripts"    # location of this script
fileSamplesPath="/proj/snic2020-16-175/analysisFileLists"    # location of file_list.csv
rawDataPath="/proj/snic2020-16-175/rawData/atacseq/$rawDataFolder" # The whole path needs to be specified because several raw data files have the same name
runPath="/proj/snic2020-16-175/nobackup/processedData/atacseq" # location of output files for this run, will be created if needed
#################################################

########## SET NAME FOR THIS RUN ##########
runName=${fileSamples%.txt} # the name of the run, will be based on the name of the input file with samples

########## READ INPUT FILES FROM $fileSamples ##########
if [ ! -f $fileSamplesPath/$fileSamples ]; then
  echo "ERROR: Can't find the input file list, run terminated!"
  exit
fi

########## SET PROCESSED AND EXPORT DIRS ##########
filesToAnalyse=($(awk -F"," -v OFS=" " '{print $3; print $4}' $fileSamplesPath/$fileSamples)) #reads filenames into array

rawDataDir=`find $rawDataPath -name ${filesToAnalyse[0]}` # full path to the input files
rawDataDir=$(dirname "${rawDataDir}")

########## SET RUN DIRS AND LOG FILE ##########
runDir=$runPath/$runName # run specific output data will be in a folder with same name as the $fileSamples input file, but in runPath location
mkdir -p $runDir
log_file_name=$runDir/${runName}"_nextflow_atacseq_log.txt"
echo "Log of this run: ${log_file_name}"
exec &> $log_file_name # write all output and error to this file from here on

# Format paths of files to analyze in Design
formatted_design=$runDir/$runName"_design.csv"
if [ -f $formatted_design ]; then
	rm $formatted_design
fi

# Adds header to design file
awk -v rawdir=$rawDataDir"/" 'BEGIN{FS=OFS=","}{$3 = rawdir $3; $4 = rawdir $4; print}' $fileSamplesPath/$fileSamples > $formatted_design  # adds path to each of the files
printf '%s\n%s\n' "group,replicate,fastq_1,fastq_2" "$(cat $formatted_design)" > $formatted_design # adds header to design file

# if raw files not found, exit
for file in ${filesToAnalyse[@]}; do 
  if [ ! -f $rawDataDir/${file} ]; then
    echo "ERROR: Can't find the input file ${file} (and possibly more after that). Run terminated!"
    exit
  fi
done

########### GENERAL FEEDBACK BEFORE START ##########
echo "########## GENERAL FEEDBACK BEFORE START ##########"
printf "\nDefaults:\n"
echo "Location of the script that does this analysis = $scriptPath"
echo "Location of the input file list with samples to analyse = $fileSamplesPath"
echo "Location of the formatted design with input list to analyse = $formatted_design"
echo "Location of the raw data files = $rawDataDir"
echo "Location of the output files to this run = $runDir (will be created if not existing before run start)"
echo "Location of the log/error file for this run = $log_file_name"
printf "\nInputs:\n"
echo "Text file with samples to analyse = $fileSamples (if not located in same folder as script, the full path should be included)"
echo "Samples were fould in this location: $rawDataDir/"
echo "Samples found in this file:"
printf '\t%s\n' "${filesToAnalyse[@]}"
echo ""

# SETTING CACHE DIR FOR SINGULARITY
export NXF_SINGULARITY_CACHEDIR=/proj/snic2020-16-175/tools/nextflow/work/singularity/

# RUNNING THE PIPELINES
cd $runDir

if [ $resume == true ]; then
	echo "RESUMING......"
	nohup nextflow run nf-core/atacseq -r 1.0.0 -resume $runName \
	--project naiss2023-22-550 --email $email \
	--genome mm10 -profile uppmax \
	--design $formatted_design \
	--outdir $runDir \
	-c /proj/snic2020-16-175/tools/nextflow/atacseq_config \
	-with-report -with-trace -with-timeline \
	> $log_file_name 2>&1 &
fi
if [ $resume == false ]; then	
	echo "Submitting job"
	nohup nextflow run nf-core/atacseq -r 1.0.0 \
	--project naiss2023-22-550 --email $email -name $runName \
	--genome mm10 -profile uppmax \
	--design $formatted_design \
	--outdir $runDir \
	-c /proj/snic2020-16-175/tools/nextflow/atacseq_config \
	-with-report -with-trace -with-timeline \
	> $log_file_name 2>&1 &
fi
