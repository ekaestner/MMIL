% startup_MMPS
%
% early mod: 06/08/11 by Don Hagler
% last mod:  10/28/14 by Don Hagler
%

fprintf('%s: starting up MMPS...',mfilename);

% for processing fif format MEG data
neuromag_tb = '/opt/neuromag/meg_pd_1.2';
if exist(neuromag_tb,'dir'), path(path,neuromag_tb); end;

MMPS_extmat = getenv('MMPS_EXTMAT');
if ~isempty(MMPS_extmat), addpath(genpath(MMPS_extmat)); end;

MMPS_matlab = getenv('MMPS_MATLAB');
if ~isempty(MMPS_matlab), addpath(genpath(MMPS_matlab)); end;

% remove new fieldtrip (only used for OpenMEEG and it causes problems otherwise)
fieldtrip_dir = [MMPS_extmat '/fieldtrip-20110912'];
if exist(fieldtrip_dir,'dir')
  rmpath(genpath(fieldtrip_dir));
end;

% use some functions for reading dicoms and PET
spmdir = '/usr/pubsw/packages/spm/spm5b';
if exist(spmdir,'dir'), addpath(genpath(spmdir)); end;

user = getenv('USER');
nonsvn = ['/home/' user '/matlab/MMPS_nonsvn'];
if exist(nonsvn,'dir')
	dlist = dir([nonsvn '/*.m']);
	if ~isempty(dlist)
		addpath(nonsvn);
	end;
end;

clear MMPS_extmat MMPS_matlab neuromag_tb fieldtrip_dir spmdir user nonsvn dlist

fprintf(' finished\n');

