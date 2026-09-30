WorkDir=$1
File=$2
Output=$3
TrigMode=$4


cd $WorkDir
cd ../..
eval `scramv1 runtime -sh`
cd -

echo Input files are: $File
echo output dir is: $Output
echo Trigger mode: $TrigMode

runAsymmetry $File $Output $TrigMode 
