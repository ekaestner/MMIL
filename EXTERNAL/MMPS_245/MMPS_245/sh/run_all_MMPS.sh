#!/bin/bash
#
# Run multiple MMPS processing and analysis steps
# The script waits for jobs to be finished.
# Usage:
#    /usr/bin/nohup run_all_MMPS.sh -p PROJID > log.txt &
#    run_all_MMPS.sh
#
# Created:  07/22/13 by Hauke Bartsch
# Prev Mod: 08/06/14 by Don Hagler
# Last Mod: 01/26/17 by Don Hagler
#

usage()
{
cat <<EOF
usage: $0 options

This script will run a series of processing steps as defined in a configuration script.
In order to work this script needs to run until the end of all processing steps. Use
nohup (or screen) to be able to disconnect from the shell (see example).

OPTIONS:
   -p      project name (e.g. PING)
   -s      suffix (e.g. proc)
   -m      cluster to use (optional, default is mmilcluster4)
   -f      force restart (removes the lock file)
   -c      configuration file to use (optional)
   -l      do not run anything, just print what would be done

Example:
  /usr/bin/nohup `basename $0` -p PING -m mmilcluster4 > log.txt &

EOF
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
getValue=$DIR/getValue.sh

proj=
suffix=
cluster=mmilcluster4
force=0
testrun=0
configuration=
defaultconfiguration=$DIR/../parms/default_ProcSteps.csv

while getopts "lfp:s:c:m:" OPTION
do
     case $OPTION in
         p)
             proj=$OPTARG
             ;;
         s)
             suffix=$OPTARG
             ;;
         m)
             cluster=$OPTARG
             ;;
         c)
             configuration=$OPTARG
             ;;
         f)
             force=1
             ;;
         l)
             testrun=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $proj ]]; then
     usage
     exit 1
fi

if [ ! -f "$configuration" ]; then
     conf=~/ProjInfo/${proj}/${proj}_ProcSteps.csv
     if [ -e "$conf" ]; then
        configuration=$conf
     else
        configuration=$defaultconfiguration
     fi
     if [ -r "$configuration" ]; then
       echo "Using: $configuration"
     else
       echo "Error: could not find configuration file in $conf or $defaultconfiguration"
       exit 1
     fi
fi

if [ ! -r "$configuration" ]; then
   echo "Error: configuration file not readable ($configuration)"
   exit 1
fi

lockfile=~/.lock-${proj}
if [ ! -z $suffix ]; then
  lockfile=${lockfile}-${suffix}
fi

if [ "$force" == "1" ]; then
  if [ -e "${lockfile}" ]; then
     echo "Warning: a lock file has been removed..."
     rm -f "${lockfile}"
  fi
fi

if [ -e "${lockfile}" ]
then
   echo "Error: log file (${lockfile}) exists, last job might not be finished yet."
   echo "  If you are sure this is an error delete the lock file and try again."
   exit 0
else
   touch "${lockfile}"
fi

date

#
# define proj, batchname
# returns jobids string
#
function getJobIDs()
{
  d=`pwd`
  cd ~/batchdirs/${proj}_${batchname}/
  joblist=`cat scriptlist.txt`
  cd $d
  jobids=''
  for job in $joblist
  do
    jobids="$jobids ${job}"
  done
}

#
# define proj, exam, jobids
# returns once jobs are done
function waitForJobs()
{
  while [ 1 = 1 ]
  do
    sleep 10
    stillworking=0
    valSum=`ssh ${cluster} qstat -r`
    if [ $? -ne 0 ]
    then
       sleep 5
       /bin/echo -ne "[try again]";
       valSum=`ssh ${cluster} qstat -r`
    fi
    for job in $jobids
    do
      #val=`ssh ${cluster} qstat -j ${proj}_${exam}_$job`
      val=`echo $valSum | grep ${proj}_${batchname}_$job`
      if [ "$val" != "" ]
      then
        stillworking=1
      fi
    done
    if [ "$stillworking" == 0 ]
    then
       printf "\n---------step is done---------\n"
       break
    else
       /bin/echo -ne ".";
    fi
  done
}

# now loop through all the jobs in [project name].csv
step=0
while [ 1 ]; do
  T="$(date +%s)"
  parms=`$getValue ${configuration} parms $step`
  if [ $? -ne 0 ]
  then
     break;
  fi

  command=`$getValue ${configuration} command $step`
  parms=$(echo $parms | sed -e "s/\${proj}/$stepname/g")
  #parms=$(echo $parms | sed -e "s/{/\\{/g")
  #parms=$(echo $parms | sed -e "s/}/\\}/g")
  parms=$(echo $parms | sed -e "s/\.\.\./\ /g")
  proc=`$getValue ${configuration} type $step`
  batchname=`$getValue ${configuration} batchname $step`
  clustercmd=`$getValue ${configuration} cluster $step`

  echo "#"
  echo "# STEP $step ($command)"
  echo "#"

  if [ "$proc" == "matlab" ]; then
    if [ "${batchname}" != "" ]; then
       parms="${parms}; parms.batchname = [ '$batchname' ]"
    fi
    if [ "${parms}" == "" ]; then
       cmd="${command}('$proj');"
    else
       cmd="parms = []; ${parms}; args=mmil_parms2args(parms); ${command}('$proj', args{:});"
    fi

    if [ $testrun == "1" ]; then
       echo "matlab -nosplash -nojvm -r \"try $cmd exit; catch e; fprintf('ERROR in matlab: %s\n',e.message); exit; end;\""
    else
       matlab -nosplash -nojvm -r "try $cmd exit; catch e; fprintf('ERROR in matlab: %s\n',e.message); exit; end;"
    fi
    if [ "${clustercmd}" == "" ]; then
       echo "no cluster run required..."
    else
       if [ "${batchname}" == "" ]; then
          cmd="ssh ${cluster} ${clustercmd} $command"
       else
          cmd="ssh ${cluster} ${clustercmd} ${proj}_${batchname}"
       fi
       if [ "$testrun" == "1" ]; then
          echo $cmd
       else
          echo "now run $cmd"
          eval $cmd
          echo "did run $cmd"
       fi
    fi
    if [ "$testrun" == "0" -a "${clustercmd}" != "" ]; then
       if [ "$batchname" == "" ]; then
          batchname=$command
       fi
       getJobIDs
       waitForJobs

       for u in `find ~/batchdirs/${proj}_${batchname}/pbsout/ -type f -name "*.err" -not -size 0`
       do
          echo "Check error log: $u"
       done
    fi
  else  # in case we do not have matlab run this using the 
    if [ "$testrun" == "1" ]; then
      echo "/usr/bin/env ${proc} ${command}"
    else
      /usr/bin/env ${proc} ${command}
    fi
  fi

  T="$(($(date +%s)-T))"
  echo "Timing (step $step): $T seconds"
  let step=step+1
done

if [ -e "${lockfile}" ]
then
   rm -f "${lockfile}"
else
   echo "Error: lock file (${lockfile}) could not be found at end of processing"
   echo " A lock file is created at the beginning of the process and needs to be"
   echo " removed at the end."
fi

date
