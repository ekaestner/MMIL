MMPS_extmat = getenv('MMPS_EXTMAT');
fieldtrip_dir = [MMPS_extmat '/fieldtrip-20110912'];
if exist(fieldtrip_dir,'dir')
  % remove redundant spm path
  spmdir = '/usr/pubsw/packages/spm/spm5b';
  if exist(spmdir,'dir'), rmpath(genpath(spmdir)); end;
  % remove old fieldtrip versions
  rmpath(genpath([MMPS_extmat '/fieldtrip-0.9.7']));
  rmpath(genpath([MMPS_extmat '/fieldtrip-0.9.7_private']));
  rmpath(genpath([MMPS_extmat '/fieldtrip-20080624']));
%  rmpath(genpath([MMPS_extmat '/fieldtrip-20080624_private']));
  % add new fieldtrip
  addpath(genpath(fieldtrip_dir));
end;

