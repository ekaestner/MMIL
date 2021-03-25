#!/bin/csh
# ping_unzip_anon_all.csh

set PING_PATH_IN='/space/md10/10/data/MMILDB/PING/rsync'
cd $PING_PATH_IN

#foreach d (`ls -d P041*`)
foreach d (P0319)
   set subjid = $d
   if (-d $subjid) then
      qsub $HOME/bin/ping_unzip_anon.csh $subjid
   endif
end
