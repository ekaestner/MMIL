function varargout = fs_resample_frames(mghfile, varargin)
%function varargout fs_resample_frames(mghfile, varargin)
%
% Purpose:
%   This function does temporal resampling on an mgh timeseries,
%   allowing users to add, remove, or simply change
%   frames with a variety of resampling algorithms
%
% Input Parameters:
%   mghfile  : source mghfile to resample from
%   varargin : optional arguments as 'name' value pairs
%
% Named Parameters:
%   window_size : # of input frames to use to create a single output frame.
%                    window_size = 1 means output has same # of frames as input
%                    window_size = 5 means output has 1/5 of frames as input
%                    NOTE: If window_size is not an integer, then we need to 
%                      interpolate frames, which limits the algorthims that can
%                      be used.
%      OR
%   frames_out : number of frames you want to have in your output mgh file
%   
% Optional Parameters:
%   tfirst       : [A] first timepoint to render
%                  {default: 1}
%   tlast        : [A] last timepoint to render.
%                  {default: Inf (use all)}
%   algorithm  : resampling algorithm
%                   FOR RESAMPLING INVOLVING INTERPOLATION:
%                       linear - use linear algorithm to either interpolate a frame, or
%                                to derive a new frame from multiple values
%
%                   FOR RESAMPLING VIA SELECTING FRAMES:
%                       average - average vertex f-values over frames within each resampling window
%                       max     - find the max vertex f-value over frames within each resampling window 
%                       first   - use first frame from each resampling window
%                       mid     - use midpoint frame from each resampling window
%                       last    - use last frame from each resampling window
%
%                {default: linear}
%   outfile    : output filename (if not a full path, then file will be written
%                to the current directory.
%                {default: {mghfile}-resampled-{algorithm}-fr{framesOut}}
%
%
%
%  Created:   07/10/07 by Ben Cipollini
%  Last Mod:  05/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, {...
  'window_size',[],[0,100],...
  'frames_out',[],[1:10000],...
  'algorithm','linear',{'linear','average', 'max', 'first','mid','last'},...
  'tfirst',1,[0 Inf],...
  'tlast',Inf,[0 Inf],...
... %I/O args
  'outfile',[],[],...
... %user prefs
  'force',false,[false true],... 
  'verbose',true,[false true], ...
  'logfile',[],[],...
  'logfid',1,[],...
},false); %% todo: why false?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist(mghfile,'file'))
  mmil_error(parms, 'file %s not found', mghfile);
end;

% get the mgh
mmil_logstr(parms, 'loading mgh file %s...',mghfile);
[vol, M, mr_parms, volsz] = fs_load_mgh(mghfile);

if (isinf(parms.tlast))
  parms.tlast = volsz(4)-1;
end;

% check tfirst and tlast
if (parms.tfirst < 1) 
  mmil_error(parms, 'tfirst must be >=1');

elseif (parms.tlast > volsz(4))
  mmil_error(parms, 'tlast must be <= the last frame index (%d)', volsz(4));

elseif (parms.tfirst > parms.tlast)
  mmil_error(parms, 'first frame (%d) > last frame (%d)',...
                     parms.tfirst, parms.tlast);

% to make calculations easy, let's just choose 
% the first & last frames right here & now
else
  vol   = vol(:,:,:,parms.tfirst:parms.tlast);
  volsz = size(vol);
end;

% calculate "other" param: window_size vs. frames_out  
if (~isempty(parms.window_size))
  parms.frames_out = ceil((volsz(4))/parms.window_size);

elseif (~isempty(parms.frames_out))
  parms.window_size = (volsz(4)) / parms.frames_out;

else
  mmil_error(parms, 'must specify either window_size or frames_out parameter.');
end;

% construct output filename
if (isempty(parms.outfile))
  parms.outfile = sprintf('%s-resampled-%s-t%d-t%d-win%.0f.mgh',...
    regexprep(mghfile, '\.mgh$', ''),parms.algorithm, ...
    parms.tfirst,parms.tlast,parms.window_size);
end;

% No-op if output file exists.
if (~parms.force && exist(parms.outfile,'file'))
  mmil_logstr(parms, 'NOTICE: not recomputing resampling: output file %s exists', parms.outfile);
  return;
% Just copy if there's no calculation necessary.
elseif (parms.window_size==1 || parms.frames_out==volsz(4))
  mmil_logstr(parms, 'NOTICE: no resampling needs to be done; input mgh copied to output destination.')
  copyfile(mghfile, parms.outfile);
end;

% No interpolation needed; we can use our
% custom 'picking & choosing' algorithm.
if (mmil_isint(parms.window_size))

  % We have a window from which we will choose frames.
  % We'll loop over frames, grabbing a window-full of
  % frames to operate on.  After processing that window-full
  % of frames into a single frame via our "algorithm", we'll
  % store the result as a frame in a new mgh matrix, and
  % advancing the window to the next set of frames.

  window            = 1:parms.window_size;
  num_output_frames = ceil(volsz(4)/parms.window_size);
  new_mgh           = zeros( [volsz(1:3) num_output_frames] );
  new_frame_idx     = 1;

  while (~isempty(window))
    mmil_logstr(parms, 'Processing window [%d] of [%d]', new_frame_idx, num_output_frames);

    % Check to see if final window has fewer frames than other windows;
    % if so, it's something the user should know.
    if (length(window) >0 && length(window)~=parms.window_size)
      mmil_logstr(parms, 'WARNING: last window has only %d frame(s) in it.', length(window));
      end;


    curSamples = vol(:,:,:,window);

    switch (parms.algorithm)
      case {'linear','average'}
        new_frame = mean(curSamples,4);

      case 'max'
        new_frame = max(curSamples,[],4);

      case 'first'
          new_frame = curSamples(:,:,:,1);

      case 'mid'
          new_frame = curSamples(:,:,:,ceil(length(windows)/2));

      case 'last'
          new_frame = curSamples(:,:,:,end);

      otherwise
        mmil_error(parms, 'algorithm %s does not apply to resampling via selection',...
                   parms.algorithm);
    end;


    new_vol(:,:,:,new_frame_idx) = reshape(new_frame, [volsz(1:3) 1]);
    new_frame_idx = new_frame_idx + 1;

    % By choosing a window in this way, the last frame
    % may not have the same of frames in it
    % as the others.  Depending on who you ask, this
    % may be good, this may be bad.
    window = window(end)+1:window(end)+length(window);
    window = window(find(window<=volsz(4)));
  end; % while

% Need interpolation
else
  %default value specific to interpolating
  if (isempty(parms.algorithm))
  parms.algorithm = 'linear';

elseif (~ismember(parms.algorithm, {'linear', 'zod'}))
    mmil_error(parms, 'algorithm %s does not apply to resampling that uses interpolation.',...
               parms.algorithm);
end;
  % Create a "timeseries" object (first dimension is time)
  ts_dim = [ length(volsz) 2:length(volsz)-1 1 ];
  ts = timeseries(permute(vol,ts_dim));

  % resample ("time" is 0-based, so we have to subtract 1 from MATLAB indices)
  ts_new = resample(ts, 0:parms.window_size:(volsz(4)-1), parms.algorithm);

  % switch our time back to the last dimension to get mgh format.
  new_vol = permute(ts_new.data, ts_dim);
end;

% save out the mgh
mmil_logstr(parms, 'saving mgh file %s', parms.outfile);
fs_save_mgh(new_vol, parms.outfile, M, mr_parms);

mmil_logstr(parms, 'Done.', parms.outfile);
