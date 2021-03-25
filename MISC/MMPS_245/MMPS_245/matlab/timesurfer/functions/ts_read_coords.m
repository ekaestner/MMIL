function varargout = ts_read_coords(varargin)
% Use to read a 4-column ascii file containing label & xyz coordinates
%
% useage:
%   ts_read_coords
%   ts_read_coords('coordfile')
%   ts_read_coords('coordfile','matfile')
%   ts_read_coords('coordfile','matfile','fieldname')
%   ts_read_coords('coordfile',STRUCT)
%   ts_read_coords('coordfile',STRUCT,'fieldname')
%   elec = ts_read_coords;
%   elec = ts_read_coords('coordfile')
%   data = ts_read_coords('coordfile','matfile')
%   data = ts_read_coords('coordfile','matfile','fieldname') 
%   data = ts_read_coords('coordfile',STRUCT)
%   data = ts_read_coords('coordfile',STRUCT,'fieldname')
%
% inputs:
%   coordfile - text file w/ coordinates (will be prompted if not supplied)
%   matfile   - mat file w/ TS structure to which coords should be added
%   STRUCT    - TS structure to which coords should be added
%   fieldname - the field name to use for coordinates (default: 'coords')
%     note: this is useful when adding multiple sets of coordinates to the
%     same TS structure (ex. 'mni_coords','talairach_coords','mri_coords')
%
% outputs: 
%   - results are saved if no output is given.
%   - returns/saves 'elec' structure if no matfile or struct is given.
%   - returns/saves STRUCT structure if matfile or struct is given.
%
% Created by JSS on 13-May-2009

tsflag = 0;

% get the ascii coordinate file (required)
if nargin == 0
  [filename,pathname] = uigetfile('*.txt','Pick a text file with labels and xyz coordinates');
  coordfile = fullfile(pathname,filename);
else
  if ischar(varargin{1})
    coordfile = varargin{1};
  else
    error('The first parameter must be a text file with xyz coordinates.');
  end
end
elec = readcoords(coordfile);

% get timesurfer structure (if appropriate)
if nargin > 1
  if ischar(varargin{2})
    % load data from user-supplied matfile
    matfile = varargin{2};
    if exist(matfile,'file')
      S = load(matfile);
      f = fieldnames(S);
      if ts_object_info(S.(f{1}))
        data = S.(f{1});
      else
        error('%s does not contain timesurfer data.',matfile);
      end
      clear S f
    else
      error('%s does not exist.',matfile);
    end
  elseif isstruct(varargin{2}) && any(ts_object_info(varargin{2}))
    % use the user-supplied timesurfer structure
    data   = varargin{2};
  else
    error('The second parameter must be a mat file containing timesurfer data.');
  end
  tsflag   = 1;
  datatype = ts_object_info(data);
  if nargin > 2 && ischar(varargin{3})
    coorfield = varargin{3};
  else
    coorfield = 'coords';
  end
end

if tsflag
  % add coords to the sensor_info field
  [sel1 sel2] = match_str({data.sensor_info.label},{elec.label});
  for k = 1:length(sel1)
    data.sensor_info(sel1(k)).(coorfield) = elec(sel2(k)).coords;
  end
  % report differences between sensors in the structure and file
  dextra = {data.sensor_info(setdiff(1:length(data.sensor_info),sel1)).label};
  fextra = {elec(setdiff(1:length(elec),sel2)).label};
  if ~isempty(fextra)
    fprintf('Labels in the coordinate file not found in the timesurfer data:\n');
    for k = 1:length(fextra), fprintf('%s ',fextra{k}); end
    fprintf('\n');
  end
  if ~isempty(dextra)
    fprintf('Labels in the timesurfer data not found in the coordinate file:\n');
    for k = 1:length(dextra), fprintf('%s ',dextra{k}); end
    fprintf('\n');    
  end  
end

if nargout == 0
  % save data
  if tsflag
    eval(sprintf('%s = data;',datatype));
    save(matfile,datatype);
  else
    save('elec.mat','elec');
  end
elseif nargout == 1
  % return data
  if tsflag
    varargout{1} = data;
  else
    varargout{1} = elec;
  end
else
  error('type ''help ts_read_coords'' for info.');
end


% SUBFUNCTION
% read the coordinate information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elec = readcoords(filename)
if ~exist(filename, 'file')
  error('could not open coordinate file: %s', filename);
end
[Lbl,X,Y,Z] = textread(filename,'%q %f %f %f');
for k = 1:length(X)
  elec(k).label = Lbl{k};
  elec(k).coords = [X(k) Y(k) Z(k)];
end