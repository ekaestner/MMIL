function ts_avg2fif(data,varargin)
% Purpose: saves each condition in avg_data structure as fif file
%
% Usage: ts_avg2fif(data,'template_fif',val1,'key2', val2, ... );
%
% Example: ts_avg2fif(avg_data,...
%                       'template_fif',    '/rawdata/template.fif',...
%                       'outstem,          'home/user/fifs/template_',...
%                       'evcode_flag'      ,1,...
%                       'forceflag'        ,1)
%
% Inputs:
%      data: timesurfer average data structure
%
% Parameters: default : options :  description
%     template_fif: [] :: full pathname of fif file to be used as template
%                         (commonly the raw data .fif used to generate the
%                         average data) Note: This parameter is REQUIRED.
%     outstem: [] :: path to output files up to and including prefix
%     evcode_flag: 1 : 0,1 : whether to append filenames with event codes (1)
%                            or condition numbers (0)
%     forceflag: 1 : 0,1 : whether to overwrite existing fif files
%
% created:  11/24/05 by Don Hagler
% last mod: 03/26/09 by Jason Sherfey
% last mod: 08/05/09 by Don Hagler
% last mod: 02/09/11 by Burke Rosen
data  = ts_checkdata_header(data);
% backcompatability
if ~sum(ismember(varargin(cellfun(@ischar,varargin)),'template_fif'))
    parms.template_fif = varargin{1};
    if length(varargin) == 4
        parms.outstem = varargin{2};
        parms.evcode_flag = varargin{3};
        parms.forceflag = varargin{4};
    elseif length(varargin) == 3
        parms.outstem = varargin{2};
        parms.evcode_flag = varargin{3};
        parms.forceflag = 1;
    elseif length(varargin) == 2
        parms.outstem = varargin{2};
        parms.evcode_flag = 0;
        parms.forceflag = 1;
    elseif length(varargin) == 1
        parms.outstem = [pwd '/'];
        parms.evcode_flag = 0;
        parms.forceflag = 1;
    end
else
    parms = mmil_args2parms(varargin,...
                            {'outstem',[pwd '/'],[],...
                             'evcode_flag',[0],{0,1},...
                             'forceflag',[1],{0,1},... 
                             'template_fif', [],[],... 
                            },false);
end


% if (~mmil_check_nargs(nargin,3)) return; end;
% if ~exist('evcode_flag','var') || isempty(evcode_flag), evcode_flag=0; end;
% if ~exist('forceflag','var') || isempty(forceflag), forceflag=1; end;

[fpath,fstem]=fileparts(parms.outstem);
if ~isempty(fpath) && ~exist(fpath,'dir')
  [s,m] = mkdir(fpath);
  if ~s
    error('failed to create output directory %s:\n%s',fpath,m);
  end;
end;

% get info from online average
if ~exist(parms.template_fif,'file')
  error('template fif file %s not found',parms.template_fif);
end;

megmodel('head',[0,0,0],parms.template_fif);
[na,ki,nu,ct,t]=channames(parms.template_fif);

% save each condition as a fif file
nconds = length(data.averages);
for j=1:nconds
  if parms.evcode_flag
    evcode = data.averages(j).event_code;
    outfile = sprintf('%s_event%d.fif',parms.outstem,evcode);
  else
    outfile = sprintf('%s_cond%d.fif',parms.outstem,j);
  end;
  if ~exist(outfile,'file') || parms.forceflag
    B = double(data.averages(j).data);
    sfreq = data.sfreq;
    t0 = data.averages(j).time(1);
    % added by JSS on 26-Mar-2009
    if size(B,1) < length(na)
      [sel1 sel2] = match_str({data.sensor_info.label},na);
      na = na(sel2);
    end
    savefif(outfile,B,na,sfreq,t0);
  end;
end
