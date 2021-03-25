function data = ts_load_data(datafile,varargin)
% Purpose: input a single data file and return continuous data
%
% Usage: data = ts_load_data(datafile)
%
% Example: cont_data = ts_load_data('/home/user/datafile.fif')
%
% Inputs: 
%     datafile: the full path to a datafile in one of the supported formats:
%     .fif, .edf, .set, .eeg, .cnt, .avg, .dat, or .mat
%
% Outputs:
%     data: timesurfer data struct of continuous data
%   
% Parameters:
%     If loading .fif data, see ts_MNE_loadfif
%     If loading .edf data, see ts_loadedf
%     If loading .set data, see ts_iEEG_eeg2epoch
%     If loading .vhdr data,see ts_iEEG_eeg2epoch
%     If loading .eeg data, see ts_iEEG_eeg2epoch
%     If loading .cnt data, see ts_loadcnt
%     If loading .avg data, see ts_iEEG_neuroscan2avg
%     If loading .nsX data, see ts_load_nsX
%     If loading .mat data, there are no additional parameters
%     If loading .dat data
%         channels :: indices of channels to return
%         is_continuous : 1 : 0,1 : whether to return continuous data or epochs
%         dio_file :: file containing trigger info from NSpike system
%         exp_file :: file containing experiment structure associated with
%                     NSpike recording
%         Note: dio_file and exp_file are necessary only when returning
%           epoched data (i.e., is_continuous = 0).
%         Undocumented options: max_epochs, coords.
%         Also see ntools_gen_epoch_data
%
% Created:  05-Jan-2011, Jason Sherfey
% Last mod: 17-May-2011, Jason Sherfey

data = [];
if ~ischar(datafile)
  error('Input must be a single data file.\nTip: use ts_process_data to process multiple files at once.');
  return;
elseif ~exist(datafile,'file')
  error('Data file does not exist: %s',datafile);
  return;
end

[a,b,ext] = fileparts(datafile); 
switch lower(ext)
  case '.fif'                                           % Neuromag
    data = ts_MNE_loadfif(datafile,varargin{:});        
  case '.edf'                                           % European data format
    data = ts_loadedf(datafile,varargin{:});            
  case '.vhdr'                                          % NKI BrainVision
    data = ts_iEEG_eeg2epoch(datafile,varargin{:});     
  case '.eeg'                                           % Neuroscan (epochs)
    data = ts_iEEG_eeg2epoch(datafile,varargin{:});     
  case '.cnt'                                           % Neuroscan (continuous)
    data = ts_loadcnt(datafile,varargin{:});            
  case '.avg'                                           % Neuroscan (average)
    data = ts_iEEG_neuroscan2avg(datafile,varargin{:}); 
  case '.set'                                           % EEGLAB
    data = ts_iEEG_eeg2epoch(datafile,varargin{:});    
  case {'.ns1','.ns2','.ns3','.ns4','.ns5'}             % Blackrock
    data = ts_load_nsX(datafile,varargin{:});
  case '.mat'                                           % Matlab
    s = load(datafile);
    f = fieldnames(s);
    x = s.(f{1});
    clear s f
    if ts_object_info(x)                                    % this is a TimeSurfer structure
      data = x;
    elseif isnumeric(x)                                     % this is a data matrix
      data = ts_matrix2data(x,varargin{:});
    else
      fprintf('Error: MAT file data not recognized.\n');
      return;
    end
  case '.dat'                                           % NSpike
    % use ntools to convert *.dat into TimeSurfer structure
    parms = mmil_args2parms( varargin, ...
                             {  'dio_file',[],[],...
                                'exp_file',[],[],...
                                'channels',[],[],...
                                'max_epochs',[],[],...
                                'coords',[],[],...
                                'is_continuous',1,[],...
                             }, ...
                             false );
    if parms.is_continuous == 0
      if isempty(parms.dio_file) 
          error('you must specify a dio file to epoch NSpike data.');
      elseif ~exist(parms.dio_file,'file')
        error('dio file does not exist: %s',parms.dio_file);
      elseif isempty(parms.exp_file)
          error('you must specify an experiment file to epoch NSpike data.');
      elseif ~exist(parms.exp_file,'file')
          error('experiment file does not exist: %s',parms.exp_file);
      end
      load(parms.exp_file,'experiment');
    else
      experiment = [];
    end
    tags = fieldnames(rmfield(parms,{'dio_file','exp_file'}));
    tags = tags(~cellfun(@(x)isempty(parms.(x)),tags));
    args = mmil_parms2args(parms,tags);
    data = ntools_gen_epoch_data(datafile,parms.dio_file,experiment,args{:});
  otherwise
    fprintf('Error: data type not recognized.\n');
    return;
end