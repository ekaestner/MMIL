function MEG_MMIL_Browse_Raw(varargin)
%function MEG_MMIL_Browse_Raw([options])
%
% Purpose: view raw MEG data using ts_browseraw
%
% Optional Input:
%   'prefix': processing prefix
%     {default = 'proc'}
%   'ContainerPath': full path of processed MEG dir
%     {default = pwd}
%   'filenum': file number of MEG data
%     {default = 1}
%
% Created:  02/22/11 by Don Hagler
% Last Mod: 10/25/11 by Don Hagler
%

parms = mmil_args2parms(varargin, { ...
  'prefix','proc',[],...
  'ContainerPath',pwd,[],...
  'filenum',1,[],...
},false); % not strict -- allows anything to pass

fname = sprintf('%s/matfiles/%s_parms.mat',parms.ContainerPath,parms.prefix);
if ~exist(fname)
  if isempty(varargin), help(mfilename); end;
  error('file %s does not exist',fname);
end;
proc = load(fname);
proc.parms.browseraw = parms.filenum;
input_data_files = proc.parms.datafile;
args = MMIL_Args(proc.parms,'ts_process_fif_data');
ts_process_fif_data(input_data_files,args{:});

