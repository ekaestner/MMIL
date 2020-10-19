#!/bin/tcsh -f

set addpathlist = ( \
/opt/neuromag/bin/vue \
/opt/neuromag/bin/util \
/opt/neuromag/bin/X11 \
)
#$PUBSW/packages/mne/bin/mne \

foreach dir ( $addpathlist )
  if ("$path" !~ *$dir*) set path = ($dir $path)
end
unset addpathlist dir

### Remove duplicates from path WITHOUT reording
set path=`echo $path | awk '{for(i=1;i<=NF;i++){if(!($i in a)){a[$i];printf s$i;s=" "}}}'`

