function ts_surf2movie(subject, srcfile, varargin)
%
%  ts_mgh2movie(subject, srcfile, varargin)
%
% Purpose:
%   make movie from freesurfer time-course files.  Cache rendered images to
%   disk for easy reference or re-use.
%
%   NOTE: since a movie can be created from a set of images (1 or more),
%     many arguments can take either a single value (to apply to all
%     images), or multiple values (one value per image in the movie).
%     SUCH PARAMETERS ARE DENOTED WITH AN [A] IN THEIR DESCRIPTION.
%
%
% Input parameters:
%   subject  : [A] string, or cell array of strings, specifying subjects on
%              whose structural data time-course overlays will be rendered
%   srcfile  : [A] string, or cell array of strings, specifying time-course
%              overlay files to use in rendering (STC and MGH formats OK)
%   varargin : optional arguments, in standard form ('name', value pairs)
%
% Optional Arguments: 
%
%   [Program Flow]
%
%   force      : regenerate all images regardless if they exist in the image
%                cache
%                {default: false}
%   nocache    : remove all image files after movie has been created
%                {default: false}
%   hashnames  : insert unique identifier string into cache image file
%                names (necessary if some input srcfile file names are
%                identical except for different paths)
%                {default: true}
%
%   [I/O args:]
%   indir        : base directory for reading src files
%                  {default: current directory}
%   outdir       : base directory for creating output files & image cache
%                  {default: current directory}
%   outstem      : prefix to prepend to output movie
%                  {default: surf2movie}
%   SUBJECTS_DIR : see freesurfer documentation 
%                  {default: environment variable SUBJECTS_DIR}
%
%   [Data args:]
%
%   hemi         : [A] hemisphere to render (lh or rh)
%                  {default: lh}
%   polar        : polar angle retinotopy data (requires _r and _i complex data)
%                  {default: false}
%   tfirst       : [A] 0-indexed value for first timepoint to render
%                  {default: 0}
%   tlast        : [A] 0-indexed value for last timepoint to render.
%                  {default: Inf (use all)}
%   [Image args:]
%
%   surface      : [A] surface to render (inflated, sphere, etc; see
%                  freesurfer docs for more info)
%                  {default: inflated}
%   fthresh      : [A] F-stat threshhold value (for colorization)
%                  {default: 0}
%   fmid         : [A] F-stat midpoint value of sigmoid (for colorization)
%                  {default: 1}
%   fslope       : [A] F-stat slope of sigmoid at fmid (for colorization)
%                  {default: 1}
%   sparsesmooth : [A] # of sparse smoothing steps to apply
%                  {default: 0}
%   postsmooth   : [A] # of post smoothing steps to apply
%                  {default: 0}
%   offset       : Intensity offset for surface
%                  {default: 0.2}
%   cvfact       : Curvature contrast factor
%                  {default: 1.5}
%   scale        : [A] 'zoom' factor
%                  {default: 1.0}
%   view         : [A] perspective of brain (lat, med, dor, ven)
%                  {default: lat}
%   rotx         : [A] custom rotation about x axis (WRT current view)
%                  {default: 0}
%   roty         : [A] custom rotation about y axis (WRT current view)
%                  {default: 0}
%   rotz         : [A] custom rotation about z axis (WRT current view)
%                  {default: 0}
%   label        : [A] label to draw under image
%                  {default: none}
%   image_size   : [A] width x height in pixels; height of each image
%                  NOTE: if image size is too large (e.g. 400x400),
%                    resulting in very large montages,
%                    movies may get bad sequence header and be unplayable
%                  {default: 200x200}
%   scalebar     : whether to show the scale or not
%                  {default: false}
%   colorbar     : whether to show the mapping of F-values to color or not
%                  {default: false}
%
%   [Montage args:]
%
%   nrows        : # of rows in your montage
%                  {default: minimum to fit all images in a 4:3 aspect grid}
%   ncols        : # of columns in your montage
%                  {default: minimum to fit all images in a 4:3 aspect grid}
%   position     : [A] index of position to place your image
%                  {default {1 ... nImages} }
%   bgcolor      : background color of movie, in HTML format
%                  {default: #000000}
%   trim         : trim edges of each image
%                  {default: true}
%   timestamp    : whether to display the current timestamp or not and, if
%                  so, where to display it. Locations are with reference to
%                  North=top-center.
%                  Possible values: {'none', 'NorthWest', 'North', 'NorthEast', 
%                                    'East','SouthEast', 'South',
%                                    'SouthWest', 'West','NorthWest'}
%                  {default: South}
%   timestamp_format : string for formatting the timestamp; must include
%                  formatting specifier for current time and final time.
%                  {default: Frame [%4d] of [%4d]}
%   timestamp_fontsize : font size in points
%                  {default: 16}
%   timestamp_fontname : font name
%                  {default: 'helvetica'}
%   timestamp_framedur : duration of one frame in msec
%                  {default: 1}
%   timestamp_prestim : duration of prestimulus period (msec)
%                  {default: 0}
%   timestamp_textcolor : color of time stamp text, in HTML format
%                  {default: #f00} (white)
%   timestamp_undercolor : color of filled area around time stamp, in HTML format
%                  {default: #808080} (gray)
%   image_type   : montaged image output format (gif, tiff, png, jpg, bmp)
%                  {default: gif}
%   delay        : # milliseconds to delay between frames
%                  {default: 5}
%   montage      : name of a pre-defined montage; will set all other optional
%                  params as necessary. 
%                  Current pre-defined montages: none.
%                  {default: 'custom'}
%  
% More documentation and examples at:
%  neuroimaging.pbwiki.com/Ts%20Movie-Maker
%
% Created:  06/01/07 by Ben Cipollini
% Rcnt Mod: 05/21/12 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

%% todo: change tfirst and tlast to be 1-based
%% todo: add documentation to MMPS documentation
%% todo: test for polar angle retinotopy data (polar flag)
%% todo:
%   frames       : [A] vector of 0-indexed values for timepoints to render
%                  supercedes tfirst and tlast if not []
%                  {default: []} (all)

if ~mmil_check_nargs(nargin,2), return; end;

%%%%%%%%%%%%%%%%
parms = mmil_args2parms(varargin,{...
  'nocache',false,[false,true],...
  'force',false,[false,true],...
  'hashnames',true,[false,true],...
  'IMAGEMAGICK_DIR',getImageMagickDir(),[],...
...
  'outdir',pwd,[],...
  'outstem','surf2movie',[],...
  'image_type','gif',{'gif','tif','tiff','png','jpg','bmp'}...
...
  'montage','custom',{'custom','onesub','comparesubs','groupsummary'},...
  'nrows',[],[1,50],...
  'ncols',[],[1,50],...
  'bgcolor','#000000',[],...
  'trim',true,[false,true],...
  'timestamp','South',{'none','NorthWest','North','NorthEast',...
                       'East','SouthEast','South','SouthWest','West','NorthWest'},...
  'timestamp_fontsize',16,[1,Inf],...
  'timestamp_fontname','helvetica',[],...
  'timestamp_format','Frame [%4d] of [%4d]',[],...
  'timestamp_framedur',1,[0.1,Inf],...
  'timestamp_prestim',0,[0,Inf],...
  'timestamp_textcolor','#f00',[],...
  'timestamp_undercolor','#808080',[],...
  'gravity','NorthEast',[],...
  ...
  'fgcolor','#ffffff',[],...%unpublishedfeature
  'transparent','#000000',[],...%unpublishedfeature
  'delay',2,[1,Inf],...
  ...
  'SUBJECTS_DIR',getenv('SUBJECTS_DIR'),[],...
  'srcfile',[],[],...
...
  'image_size','200x200',[],...
  'hemi','lh',[],...
  'surface','inflated',[],...
  'patch',[],[],...
  'fthresh',0,[0,Inf],...
  'fslope',1,[0,Inf],....
  'fmid',1,[0,Inf],...
  'sparsesmooth',0,[0,100],...
  'postsmooth',0,[0,100],...
  'scale',1.0,[0,100],...
  'view','lat',{'lat','ven','med','pos','dor','cus'},...
  'rotx',0,[-Inf,Inf],...
  'roty',0,[-Inf,Inf],...
  'rotz',0,[-Inf,Inf],...
  'label','',[],...
  'offset',0.2,[],...
  'cvfact',1.5,[],...
  'indir',pwd,[],...
  'frames',[],[],...
  'tfirst',0,[0,Inf],...
  'tlast',Inf,[0,Inf],...
  'tstep',1,[1,Inf],...%
  'step_type','average',[],...
  'position',[],[]...
  'scalebar',false,[false,true],...
  'colscalebar',false,[false,true],...
  'polar',false,[false,true],...
  'checksum','',[],...
...
  'verbose',true,[false,true],...
  'logfile',[],[],...
  'logfid',1,[]...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  check basic variables 
%    and set computed defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     (parms.bgcolor(1)     ~= '#' || length(parms.bgcolor)~=7)
  error('%s: bgcolor must have format #HHHHHH; passed value=%s', mfilename, parms.bgcolor);
elseif (parms.bgcolor(1)     ~= '#' || length(parms.bgcolor)~=7)
  error('%s: fgcolor must have format #HHHHHH; passed value=%s', mfilename, parms.fgcolor);
elseif (parms.bgcolor(1)     ~= '#' || length(parms.bgcolor)~=7)
  error('%s: transparent must have format #HHHHHH; passed value=%s', mfilename, parms.transparent);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  check dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_logstr(parms,'Checking dependencies...');

%	Just because things are set up doesn't mean
%	the binaries we want exist.  Do a cursory scan
%	to make sure all our dependencies are available.

%montage
if (~exist(fullfile(parms.IMAGEMAGICK_DIR, 'montage'), 'file'))
	error('%s: IMAGEMAGICK "montage" program not found at %s', mfilename, ...
	      fullfile(parms.IMAGEMAGICK_DIR, 'montage'));
%convert
elseif (~exist(fullfile(parms.IMAGEMAGICK_DIR, 'convert'), 'file'))
  error('%s: IMAGEMAGICK "convert" program not found at %s', mfilename, ...
	      fullfile(parms.IMAGEMAGICK_DIR, 'convert'));
%surmfgh      
elseif (mmil_unix('which fs_surfmgh.csh'))
  error('%s: fs_surfmgh.csh not found in path', mfilename);

%freesurfer
elseif (isempty(getenv('FREESURFER_HOME')))
  error('%s: FREESURFER_HOME is not set; you need to setup your freesurfer installation.', mfilename);
end;

%  if ( "`which $TKSURFER_PROG | grep /`" == "" ) goto no_tksurfer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  distribute variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_logstr(parms,'Setting up image properties...');

% generic vars first
parms.subject = subject;
parms.srcfile = srcfile;

allarrs = { 'SUBJECTS_DIR' 'srcfile' 'subject' 'checksum' 'hemi' 'surface' 'patch' 'polar' 'scalebar' 'colscalebar' ...
            'fthresh' 'fslope' 'fmid' ...
            'sparsesmooth' 'postsmooth' 'scale' 'view' 'rotx' 'roty' 'rotz' 'label' 'bgcolor' 'fgcolor' 'transparent' ...
            'offset' 'cvfact' 'indir' 'frames' 'tfirst' 'tlast' };

% get the length of each array			  
arrlens = [];
for i=1:length(allarrs)
  curarr = getfield(parms, allarrs{i});
  if iscell(curarr)
  	arrlens(i) = length(curarr);
  else
	  arrlens(i) = 1;
  end;
end;
num_images = max(arrlens);

% special vars second
if isempty(parms.position)
  parms.position=num2cell([1:num_images]);
end;

% check for inconsistencies in lengths
%   - check for all arrays with length ~= 1
%   - if there are multiple lengths, then we have a problem.
if (length(unique(arrlens(find(arrlens~=1)))) > 1)
  parms_length_below_max = allarrs(find(arrlens>1 & arrlens<num_images));
	parms_length_at_max    = allarrs(find(arrlens==num_images));
  error('Parameters had different lengths: %s', ...
	      sprintf('''%s'' ', parms_length_at_max{:}, parms_length_below_max{:}));
end;

% distribute values of non-cell arguments.
for i=1:length(allarrs)
  curarr = getfield(parms, allarrs{i});
  if ~iscell(curarr)
  	curarr={curarr};
	end;

  % length MUST be 1 if we're here
	if (length(curarr) < num_images)
	  parms = setfield(parms, allarrs{i}, repmat(curarr, [1 num_images]));
  else
    parms = setfield(parms, allarrs{i}, curarr);
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  check variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check output file
%if (exist(fullfile(parms.outdir, [parms.outstem '.mpg']), 'file'))
%  fprintf('%s: WARNING: Output file exists; exiting.\n', mfilename);
%end;

% check subjects dir & subject
for i=1:length(parms.SUBJECTS_DIR)
  if ~exist(parms.SUBJECTS_DIR{i},'file')
    error('%s: SUBJECTS_DIR ''%s'' not found.', mfilename, parms.SUBJECTS_DIR{i});
  elseif ~exist(fullfile(parms.SUBJECTS_DIR{i}, parms.subject{i}), 'file')
    error('%s: Subject''s structurals directory not found: %s', mfilename, ...
      fullfile(parms.SUBJECTS_DIR{i}, parms.subject{i}) );
  end;
end;

% Make sure the # of rows & columns can 
% support the total # of images displayed.
%
%  Calculate a grid if no grid was specified.
if ( isempty(parms.nrows) && isempty(parms.ncols) )
  parms.nrows = floor(sqrt(num_images));
end;
if (isempty(parms.nrows))
  parms.nrows = ceil(num_images/parms.ncols);
elseif (isempty(parms.ncols))
  parms.ncols = ceil(num_images/parms.nrows);
end

if ((parms.nrows*parms.ncols) < num_images)
  error(['grid of size %d x %d cannot support %d images; \n '...
	    'please increase your nrows and/or ncols, ' ...
		'or remove them as parameters.'], ...
		parms.nrows, parms.ncols, num_images );
end;

% Make sure all the grid positions 
% will fit within the current grid
if (~isempty(find([parms.position{:}]>(parms.nrows*parms.ncols))))
	error('A %d x %d grid cannot support the following grid positions:\n',...
	      sprintf('\t%d\n', parms.position{find([parms.position{:}]>(parms.nrows*parms.ncols))}));
end;

% Check datafile consistency with params and,
% if necessary, read the # of frames to use
%
% Grab all the necessary params
for i=1:length(parms.srcfile)
  if (isempty(parms.srcfile{i}))
    error('%s: empty source file in input.', mfilename);
  elseif (parms.srcfile{i}(1) ~= '/')
    parms.srcfile{i} = fullfile(parms.indir{i}, parms.srcfile{i});
  end;

  if (~exist(parms.srcfile{i}, 'file'))
    computed_file = [parms.srcfile{i} '-' parms.hemi{i} '.mgh'];
    if (exist(computed_file, 'file'))
      parms.srcfile{i} = computed_file;
    else
      error('%s: source file not found: %s', mfilename, parms.srcfile{i});
    end;
  end;

  setenv('SUBJECTS_DIR', parms.SUBJECTS_DIR{i});

	% Convert all input files to mgh
	if (strcmp(mmil_fileparts(parms.srcfile{i},'ext'),'.stc'))
	  mghfile	= regexrep(parms.srcfile{i}, '\.stc$', '.mgh');
    ts_stc2mgh(parms.srcfile{i}, mghfile);
	  parms.srcfile{i} = mghfile;
	end;

	% Resample files
	if (parms.tstep ~= 1)
	  resampled_mghfile = regexrep(parms.srcfile{i}, '\.mgh$', sprintf('-tstep%d.mgh', parms.tstep));
	  fs_resample_frames(parms.srcfile{i},...
      'window_size',parms.tstep,...
      'outfile',resampled_mghfile);
	  parms.srcfile{i} = resampled_mghfile;
	end;

  % Pull metadata
  [v,M,mr,vs] = fs_load_mgh(parms.srcfile{i}, [], [], true);
  nFrames = vs(4);

  % save out values
	if (isinf(parms.tlast{i}))
    parms.tlast{i} = nFrames-1;
  end;

	%check params
  if (parms.tlast{i} > nFrames)
	  error('%s: last frame (%d) > total # frames (%d): %s', ...
	        mfilename, parms.tlast{i}, nFrames, parms.srcfile{i});
  end; 
end;

% check tfirst and tlast
%
% Note: tlast checked within bounds previously
if (~isempty(find([parms.tfirst{:}]<0)))
  error('%s: all tfirst must be >=0', mfilename);

elseif (~isempty(find([parms.tfirst{:}]>[parms.tlast{:}])))
  error('%s: first frame (%d) > last frame (%d)', mfilename,...
         parms.tfirst{find([parms.tfirst{:}]>[parms.tlast])},...
	  	   parms.tlast{find([parms.tfirst{:}]>[parms.tlast])});
end;

% check total # of frames
if (length(unique([parms.tlast{:}] - [parms.tfirst{:}])) > 1)
  error('%s: Frames per subject don''t match:\n', mfilename,...
	    sprintf('\t%d\n', [parms.tlast{:}] - [parms.tfirst{:}]));
end;

if isempty(parms.frames)
  parms.frames = [parms.tfirst:parms.tlast];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERATE THE IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we're going to create all necessary
% images, calling fs_surfmgh.csh only if/when necessary.

%	Loop over each image set
for i=1:num_images
  mmil_logstr(parms, 'Rendering images for surface [%d] of [%d], frames [%d] - [%d]',...
                 i, num_images, parms.tfirst{i}, parms.tlast{i});

  if (isempty(parms.checksum{i}))
    parms.checksum{i} = mmil_fileparts(parms.srcfile{i},'name');
    if parms.hashnames
      parms.checksum{i} = sprintf('%s-hash%d',...
                                  parms.checksum{i},...
                                  bogushash(mmil_fileparts(parms.srcfile{i},'path')));
    end;      
    parms.checksum{i} = sprintf('%s-ft%0.3f-fm%0.3f-fs%0.3f',...
                                parms.checksum{i},...
                                parms.fthresh{i},parms.fmid{i},parms.fslope{i});
    if parms.rotx{i} | parms.roty{i} | parms.rotz{i}
      parms.checksum{i} = sprintf('%s-rx%d-ry%d-rz%d',...
        parms.checksum{i},parms.rotx{i},parms.roty{i},parms.rotz{i});
    end;
  end;

  parms.image_stem{i} = sprintf('%s-%s', parms.subject{i}, parms.checksum{i});
%    rgbdir              = fullfile(parms.outdir,'rgb');
  rgbdir              = sprintf('%s/rgb/%s',parms.outdir,parms.image_stem{i});

  if (~exist(rgbdir,'dir'))
    mkdir(rgbdir);
  end;

  rargs = {};
  if (parms.force)          rargs = {rargs{:} '-force'}; end;
  if (parms.patch{i})       rargs = {rargs{:} '-flat' '-patch' parms.patch{i}}; end;
  if (parms.scalebar{i})    rargs = {rargs{:} '-scalebar'}; end;
  if (parms.colscalebar{i}) rargs = {rargs{:} '-colscalebar'}; end;
  if (parms.polar{i})       rargs = {rargs{:} '-polar'}; end;

%    setenv('SUBJECTS_DIR', parms.SUBJECTS_DIR{i});  %% this doesn't necessarily work for some reason

  % cus view is lateral view plus rotations ('cus' is just for the file name)
  tmp_view = parms.view{i};
  if strcmp(tmp_view,'cus')
    tmp_view = 'lat';
  end;

  cmd = ['fs_surfmgh.csh ' parms.subject{i} ...
                       ' ' parms.srcfile{i} ...
                       ' ' parms.hemi{i} ...
                 ' -sdir ' parms.SUBJECTS_DIR{i} ...
                 ' -surf ' parms.surface{i} ...
               ' -tfirst ' sprintf('%d',parms.tfirst{i}) ...
                ' -tlast ' sprintf('%d',parms.tlast{i}) ...
              ' -fthresh ' sprintf('%0.3f',parms.fthresh{i}) ...
                 ' -fmid ' sprintf('%0.3f',parms.fmid{i}) ...
               ' -fslope ' sprintf('%0.3f',parms.fslope{i}) ...
               ' -offset ' sprintf('%3.1f',parms.offset{i}) ...
               ' -cvfact ' sprintf('%3.1f',parms.cvfact{i}) ...
                ' -scale ' sprintf('%3.1f',parms.scale{i}) ...
                 ' -view ' tmp_view ...
                 ' -rotx ' sprintf('%d',parms.rotx{i}) ...
                 ' -roty ' sprintf('%d',parms.roty{i}) ...
                 ' -rotz ' sprintf('%d',parms.rotz{i}) ...
               ' -smooth ' sprintf('%d',parms.postsmooth{i}) ...
              ' -outstem ' parms.image_stem{i} ...
               ' -outdir ' rgbdir ...
               sprintf(' %s', rargs{:}), ...
           ' -rmintermed'...
              ' -savergb'...
            ' -offscreen'];%...
  if parms.force
    cmd = [cmd ' -force'];
  end;

%                 ' -silent'];
  [status,result]=mmil_unix(cmd);
  if status
    error('%s: Failed to run fs_surfmgh.csh: %s\n\n%s',mfilename,result,cmd);
  end;
end; %for each image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERATE THE MONTAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output montage dir
montage_dir = fullfile(parms.outdir, 'montages');
if (~exist(montage_dir, 'dir'))
  mkdir(montage_dir);
end;

montaged_files = {};
nslices = 1 + parms.tlast{1} - parms.tfirst{1};

mmil_logstr(parms, 'Creating %d montages ...', nslices);

% montage each frame
for i=1:nslices
  montage_file          = sprintf('%s-%s-%s.%s', parms.outstem, ...
                                  mmil_padstr(i,1+floor(log10(nFrames))),...
                                  parms.image_size, parms.image_type);
  montaged_files{end+1} = fullfile(montage_dir, montage_file);
  mmil_logstr(parms, 'Montage %d of %d...', i, nslices );
  if ~exist(montaged_files{end},'file') || parms.force
    %	Loop over each image set to add the image to the list that
    % needs to be montaged
    image_files = {};
    for j=1:num_images
      t = i + parms.tfirst{j} - 1;
      rgbdir = sprintf('%s/rgb/%s',parms.outdir,parms.image_stem{j});
      image_file = [rgbdir '/' parms.image_stem{j}];
      if nFrames>1
        image_file = [image_file '-' mmil_padstr(t,1+floor(log10(nFrames)))];
      end;
      image_file = [image_file '-' parms.hemi{j} '-' parms.surface{j} '-' parms.view{j} '.rgb'];
      if parms.trim
        fname = sprintf('%s/%s-trim.rgb',rgbdir,mmil_fileparts(image_file,'name'));
        if ~exist(fname,'file') || parms.force
          cmd = sprintf('%s/convert -trim %s %s', ...
                        parms.IMAGEMAGICK_DIR, image_file,fname);
          [status,result] = mmil_unix(cmd);
          if status
            disp(cmd);
            error('%s: failed to trim image:\n%s', mfilename,result);
          end;
        end;
        image_file = fname;
      end;
      image_files{end+1} = image_file;
    end

    %	Now we have a set of $nimages images and a grid of $positions.
    %   We need to reorder images to follow the $positions grid, and to
    %   insert blank images for 'skipped' grid positions.
    gridded_images = {};

    % loop through every grid position
    for j=1:(parms.nrows*parms.ncols)

     % search for an image specified for this grid location
     for im=1:length(parms.position)
       if (parms.position{im} == j)
         break;
       end;
     end;

      % no image specified for this grid location, so use "blank"
      if (parms.position{im} == j)
        if (~exist(image_files{im},'file'))
          error('%s: can''t find surface image (error running fs_surfmgh.csh?): %s\n', ...
                mfilename, image_files{im});
        end;
        gridded_images{end+1} = [ ' -label "'       parms.label{im} '"' ...
                                  ' -background \'  parms.bgcolor{im} ...
                                  ' -geometry '    parms.image_size ...
                                  ' -stroke \'      parms.fgcolor{im} ...
                                  ' -transparent \' parms.transparent{im} ...
                                  ' ' image_files{im} ];
      else
        if (~exist(fullfile(parms.outdir, 'rgb', '_blank.gif')))
          cmd = sprintf('echo "" | "%s/convert" -background \\%s -page 1x1 text:- %s', ...
                        parms.IMAGEMAGICK_DIR, parms.bgcolor{im}, fullfile(parms.outdir, 'rgb', '_blank.gif'));
          [status,result] = mmil_unix(cmd);
          if status
            disp(cmd);
            error('%s: failed to create dummy spacer gif:\n%s', mfilename,result);
          end;
        end;
        gridded_images{end+1} = [ '+label ' fullfile(parms.outdir, 'rgb', '_blank.gif')];
      end;
    end;

    % Make the montage
    cmd = sprintf('%s -tile %dx%d %s %s', ...
                   fullfile(parms.IMAGEMAGICK_DIR,'montage'), ...
                   parms.ncols, parms.nrows,... % NOTE that imagemagick is opposite logic
                   sprintf('%s ',gridded_images{:}), ...
                   montaged_files{end});
    %% so why are we using logic opposite to imagemagick? DH

    [status,result] = mmil_unix(cmd);
    if status
      disp(cmd);
      error('%s: failed to montage rgb files:\n%s', mfilename,result);
    end;

    % Annotate the montage
	  if (~strcmp('none',parms.timestamp))
	    switch (parms.timestamp)
        case {'NorthWest','SouthWest', 'West'}
  		    timestamp_relative_left = 5;
        case {'NorthEast', 'SouthEast', 'East'}
  		    timestamp_relative_left = -50; % approx
        otherwise
  		    timestamp_relative_left = 0; % approx
	    end;

	    switch (parms.timestamp)
        case {'NorthWest','North','NorthEast'}
          timestamp_relative_top = 0;
        case {'SouthWest', 'South', 'SouthEast'}
          timestamp_relative_top = 25;
        otherwise
          timestamp_relative_top = 0; % approx
      end;

      fcur = i-1+parms.tfirst{1};
      flast = parms.tlast{1};

      tcur = round(fcur*parms.timestamp_framedur - parms.timestamp_prestim);
      tlast = round(flast*parms.timestamp_framedur - parms.timestamp_prestim);

	    parms.timestamp_format(findstr('''',parms.timestamp_format)) = '`';
      cmd = sprintf(['%s %s -undercolor \\%s -fill \\%s -pointsize %d -font %s -gravity %s -draw "text %d,%d ''' parms.timestamp_format '''" %s'], ...
                     fullfile(parms.IMAGEMAGICK_DIR,'convert'), ...
                     montaged_files{end},...
                     parms.timestamp_undercolor,...
                     parms.timestamp_textcolor,...
                     parms.timestamp_fontsize,...
                     parms.timestamp_fontname,...
                     parms.gravity,...
                     timestamp_relative_left, ...
                     timestamp_relative_top, ...
                     tcur, tlast,...
                     montaged_files{end});

      [status,result] = mmil_unix(cmd);
      if status
        disp(cmd);
        error('%s: failed to montage rgb files:\n%s', mfilename,result);
      end;
    end;
  end;
end; %for each slice


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERATE THE MOVIE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_logstr(parms, 'Concatenating images to create mpg...'); 

mpg_file  = fullfile(parms.outdir, [parms.outstem '.mpg']);
if (exist(mpg_file, 'file'))
  mmil_logstr(parms, 'NOTICE: overwriting existing movie: %s', mpg_file);
  delete(mpg_file);
end;

cmd = sprintf('%s -delay %d -adjoin %s mpg:%s', fullfile(parms.IMAGEMAGICK_DIR, 'convert'),...
              parms.delay, sprintf('%s ', montaged_files{:}), mpg_file);
[status,result] = mmil_unix(cmd);
disp(cmd);
disp(result);
if status
  if (exist(mpg_file,'file'))
    mmil_logstr(parms, 'WARNING: movie was created, but with warnings:\n\n%s',result);
	else
    disp(cmd);
    error('failed to create mpg from montages:\n%s',result,cmd);
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove intermediates
if ( parms.nocache )
  mmil_logstr(parms, 'Removing intermediate cache files');
  rmdir(fullfile(parms.outdir,'rgb'), 's');
  rmdir(fullfile(parms.outdir,'montages'), 's');
end;

% ALWAYS clean up after fs_surfmgh.csh
if (length(dir('tempsurf.tcl.*')) > 0)
  delete('tempsurf.tcl.*');
end;
if (exist('error.log','file'))
  delete('error.log');
end;

mmil_logstr(parms, 'DONE.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getImageMagickDir() %getMontageParms(montageName)
%
%
  [status,result] = unix('which montage');
  d = mmil_fileparts(result, 'path');
  if (status || ~exist(fullfile(d, 'convert'), 'file'))
    error('%s: Couldn''t determine IMAGEMAGICK_DIR.', mfilename);
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = getMontageParms(montageName)
%
% Define default array values for a named montage
%
  parms.montage=montageName;

  switch montageName
    case 'custom'
  	otherwise, error('Montage type NYI: %s', montageName);
  end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hash = bogushash(str)
  str = uint32(str(:))';
  hash = uint32(5381);
  for i=length(str):-1:2
      hash = hash + i*mod(str(i),bitshift(str(i-1),3));
  end
  hash = hash + str(1);

