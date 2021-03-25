#$ -N ping_unzip_anon
#$ -o $HOME/anonymize/pbsout
#$ -j y
#$ -V
echo
echo "Starting ping_unzip_anon.csh at"
date

set subjid=\'$1\'
set forceflag=$2

set MAT_PING_PATH=\'/home/cooper/MMPS_mat/jcrpub/PING\'
set PING_SETUP_MFILE=\'$HOME/anonymize/PING_Setup.m\'
#set PING_SETUP_MFILE=\'$HOME/space_md8_5/PING/PING_Setup.m\'

echo "matlab -r PING_Import_and_Anon_Exams($subjid,$PING_SETUP_MFILE,$forceflag)"

/usr/pubsw/packages/matlab/R2007a/bin/matlab -nosplash -nojvm -nodisplay -r  "addpath($MAT_PING_PATH); PING_Import_and_Anon_Exams($subjid,$PING_SETUP_MFILE,$forceflag); exit;"

echo
echo "Script stopped at"
date

