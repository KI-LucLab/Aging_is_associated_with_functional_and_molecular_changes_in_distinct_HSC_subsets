#! /bin/bash -l

#SBATCH -A naiss2023-22-550
#SBATCH -o tagdirToCounts.out
#SBATCH -e tagdirToCounts.err
#SBATCH -J tagdirToCounts.job
#SBATCH -n 2
#SBATCH -t 24:00:00

module load bioinfo-tools
module load samtools
module load R/3.4.3
module load R_packages/3.4.3
module load HOMER/4.11

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
processedPath="/proj/snic2020-16-175/nobackup/processedData/" # location of tagDirs, should be existing
runPath="/proj/snic2020-16-175/nobackup/processedData/" # location of output files for this run, should be existsing from basicAnalysis run
#################################################

########## SET NAME FOR THIS RUN ##########
runName=${fileSamples%.txt} # the name of the run, will be based on the name of the input file with samples

########## SET RUN DIRS AND LOG FILE ##########
runDir=$runPath/atacseq/$runName # run specific output data will be in a folder with same name as the $fileSamples input file, but in runPath location
mkdir -p $runDir
cat ${scriptPath}/tagdirToCounts.out ${scriptPath}/tagdirToCounts.err > $runDir/${runName}_analysisLog.txt # paste everything that is in output or error to one log file
exec &> $runDir/${runName}_analysisLog.txt # write all output and error to this file from here on

########## READ INPUT FILES FROM $fileSamples ##########
readarray -t filesToAnalyse < $fileSamplesPath/$fileSamples # read $fileSamples file line by line, remove newlines (-t) and read into array (filesToAnalyse)
if [ ! -f $fileSamplesPath/$fileSamples ]; then
  echo "ERROR: Can't find the input file list, run terminated!"
  exit
fi

########## PROCESS INPUT FILE LIST ##########
tagdirList=("none")
processedDir=""
i=0 # line index
for line in "${filesToAnalyse[@]}"; do
  IFS=' ' read -r -a splitLine <<< "${line}" # split line on space
  tagdirList[i]=${splitLine[0]}

  if [ $(basename "${tagdirList[i]}") != ${tagdirList[i]} ]; then # if a path is given, use it
    directory=${tagdirList[i]##$processedPath/} # directory will be the path between the set $processedPath and the actual location of the files
    processedDir=${runDir}/${runName}_homerPeakFiles # peak files will be in run output folder
    mkdir $processedDir
    exportDir=${runDir}/${runName}_peakBeds # peak bed files will be in run output folder
    mkdir $exportDir
  fi
  if [ $(basename "${tagdirList[i]}") == ${tagdirList[i]} ]; then # if no path is given, find it in processedPath
    tagdirList[i]=`find $processedPath -name ${tagdirList[i]}` # use full path to the file from now on
    directory=${tagdirList[i]##$processedPath/} # directory will be the path between the set $processedPath and the actual location of the files
    processedDir=${runDir}/${runName}_homerPeakFiles # peak files will be in run output folder
    mkdir $processedDir
    exportDir=${runDir}/${runName}_peakBeds # peak bed files will be in run output folder
    mkdir $exportDir
  fi
  i=$((i+1))
done

########## Set fragmentLength and peak type ###########
fragmentLength=50
type="-region -size 100 -minDist 200 -fragLength $fragmentLength"
typeName="region-size100-minDist200-fragLength$fragmentLength"

########### GENERAL FEEDBACK BEFORE START ##########
echo "########## GENERAL FEEDBACK BEFORE START ##########"
printf "\nDefaults:\n"
echo "Location of the script that does this analysis = $scriptPath"
echo "Location of the input file list with samples to analyse = $fileSamplesPath"
echo "Location of the tagDirs = $processedDir "
echo "Location of the output files to this run = $runDir (will be created if not existing before run start)"
printf "\nInputs:\n"
echo "Text file with samples to analyse = $fileSamples"
echo "Lines found in this file:"
printf '\t%s\n' "${filesToAnalyse[@]}"
echo ""

peakNames=""
for (( i=0; i<${#tagdirList[@]}; i++ )); do # for each tagdir in list, find peaks
  peakFileName="$(basename "${tagdirList[i]}")_peaks_${typeName}.txt"
  peakNames="$peakNames $peakFileName"

  findPeaks ${tagdirList[i]} -o $processedDir/$peakFileName $type}
  pos2bed.pl $processedDir/$peakFileName > $exportDir/${peakFileName%.txt}.bed
done

cd $processedDir
mergePeaks -d given $peakNames -code > $runDir/${runName}_${typeName}_mergedPeaks.txt
pos2bed.pl $runDir/${runName}_${typeName}_mergedPeaks.txt > $runDir/${runName}_${typeName}_mergedPeaks.bed
annotatePeaks.pl $runDir/${runName}_${typeName}_mergedPeaks.txt mm10 -d ${tagdirList[@]} -noadj -fragLength $fragmentLength > $runDir/${runName}_${typeName}_mergedPeaks_noadj.txt

echo "DONE"
