function [X,M,info] = interfileread_ADNI(varargin)
%INTERFILEREAD Read images in Interfile 3.3 format.
%   A = INTERFILEREAD(FILENAME) reads the images in the first energy window
%   of FILENAME into A, where A is an M-by-N array for a single image and
%   an M-by-N-by-P array for multiple images.  The file must be in the
%   current directory or in a directory on the MATLAB path.
%
%   A = INTERFILEREAD(FILENAME,WINDOW) reads the images in energy window
%   WINDOW of FILENAME into A.
%
%   The images in the energy window must be of the same size.
%
%   Examples
%   --------
%
%       img  = interfileread('dyna3.hdr');
%
%   The sample header and image files in the example above can be found at:
%
%   http://www.keston.com/Phantoms/
%
%   See also INTERFILEINFO.

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $   $Date: 2005/06/20 03:10:32 $

% check number of inputs
error(nargchk(1, 2, nargin));
% get file parameters
[filename, window] = parse_inputs(varargin{:});
info = interfileinfo_ADNI(filename);
params = get_params(info);
if (window > params.num_windows)
    err_id = 'Images:interfileread:invalidInput';
    err_msg = 'Invalid energy window.';
    error(err_id, err_msg);
end

% Handle LONIfied files (AMD added 8/20/06)
%if strcmp(params.data_file,'file name removed during de-identification')
if 1 % Always ignore filename in file
  if strcmp(filename((length(filename)-1):end),'.i')
    params.data_file = filename;
  else
    params.data_file = filename(1:(end-4));
  end
end

% open the image file for reading
fid = fopen(params.data_file, 'r', params.byte_order);
if (fid == -1)
    err_id = 'Images:interfileread:invalidValue';
    err_msg = 'Invalid value found at key ''name of data file''.';
    error(err_id, err_msg);
end

fseek(fid, params.offset, 'cof');

% read images according to data type
switch (lower(params.data_type))
    case {'static', 'roi'}
        X = get_stat(info, window, fid);
        
    case ('dynamic')
        X = get_dyna(info, window, fid);
        
    case ('gated')
        X = get_gated(info, window, fid);
        
    case ('tomographic')
        X = get_tomo(info, window, fid);
        
    case ('curve')
        X = get_curve(info, fid);
        
    case ('other')
        fclose(fid);
        err_id = 'Images:interfileread:unsupportedFeature';
        err_msg = 'Data of type ''other'' is not supported.';
        error(err_id, err_msg);
        
    otherwise
        fclose(fid);
        err_id = 'Images:interfileread:invalidValue';
        err_msg = 'Invalid value found at key ''type of data''.';
        error(err_id, err_msg);
end

fclose(fid);

% Specify geometry (NMvox2ras) (AMD added 8/20/06)
M = [ 0 -info.ScalingFactorMmPixel1 0 0; -info.ScalingFactorMmPixel2 0 0 0; 0 0 -info.ScalingFactorMmPixel3 0; 0 0 0 1];
M(1:3,4) = M(1:3,4)-M(1:3,:)*[(1+size(X))/2 1]'; % Assume axial acq, S-I slice ordering, radiological conventions 

return


%%%
%%% Function parse_inputs
%%%
function [filename, window] = parse_inputs(varargin)

filename = varargin{1};

if (nargin == 1)
    window = 1;
    
else
    window = varargin{2};
end

if window < 1
    err_id = 'Images:interfileread:invalidInput';
    err_msg = 'Invalid energy window.';
    error(err_id, err_msg);
end



%%%
%%% Function get_params
%%%
function params = get_params(info)

% get version of keys
%params.version = get_val(info, 'VersionOfKeys', 'version of keys', 1);
params.version = get_val(info, 'VersionOfKeys', 'version of keys', 0); % Ignore missing (AMD added 8/20/06)
if (isempty(params.version))
    params.version = 3.3;
end

% get data offset
try
    offset = 2048*get_val(info, 'DataStartingBlock', 'data starting block', 1);
    
catch
    try
        offset = get_val(info, 'DataOffsetInBytes', 'data offset in bytes', 1);
        
    catch
        offset = 0; % Assume zeros offset, if fields missing (to work with ADNI data) (added AMD 8/20/06)
%         err_id = 'Images:interfileread:missingKey';
%         err_msg = 'Missing required key ''data starting block'' or ''data offset in bytes''.';
%         error(err_id, err_msg);
    end
end

if (isempty(offset))
    params.offset = 0;
else
    params.offset = offset;
end

% get image file name
params.data_file = get_val(info, 'NameOfDataFile', 'name of data file', 1);
if (isempty(params.data_file))
    err_id = 'Images:interfileread:noDataFile';
    err_msg = 'No data file.';
    error(err_id, err_msg);
end

% get image type
%params.data_type = get_val(info, 'TypeOfData', 'type of data', 1);
params.data_type = get_val(info, 'TypeOfData', 'type of data', 0); % (AMD added)
if (isempty(params.data_type))
%    params.data_type = 'Other';
    params.data_type = 'Static'; % (AMD added)
end

% get endianness
byte_order = get_val(info, 'ImagedataByteOrder', 'imagedata byte order', 0);
%if strcmpi(byte_order, 'bigendian') || isempty(byte_order)
if strcmpi(byte_order, 'bigendian') % (AMD make default littleendian)
    params.byte_order = 'ieee-be';
    
elseif strcmpi(byte_order, 'littleendian') || isempty(byte_order) % (AMD make default littleendian)
    params.byte_order = 'ieee-le';

else
    err_id = 'Images:interfileread:invalidValue';
    err_msg = 'Invalid value found at key ''imagedata byte order''.';
    error(err_id, err_msg);
end

% get number of energy windows
params.num_windows = get_val(info, 'NumberOfEnergyWindows', 'number of energy windows', 0);
if (isempty(params.num_windows))
    params.num_windows = 1;
end



%%%
%%% Function get_stat
%%%
function X = get_stat(info, window, fid)

X = [];
i = 1;

while (i <= window)
    num_img = get_val(info, 'NumberOfImagesEnergyWindow', 'number of images/energy window', 0, fid, i);
    if (isempty(num_img))
%        num_img = 1;
        num_img = info.MatrixSize3; % (AMD added to deal with ADNI data)
    end
    
    % go through images
    for j = 1:num_img
%        num_cols = get_val(info, 'MatrixSize1', 'matrix size [1]', 2, fid, (i-1)*num_img+j);
%        num_rows = get_val(info, 'MatrixSize2', 'matrix size [2]', 2, fid, (i-1)*num_img+j);
%        num_bytes = get_val(info, 'NumberOfBytesPerPixel', 'number of bytes per pixel', 2, fid, (i-1)*num_img+j);
        num_cols = info.MatrixSize1;
        num_rows = info.MatrixSize2;
        num_bytes = info.NumberOfBytesPerPixel;
        if (i == window)
%            num_format = get_val(info, 'NumberFormat', 'number format', 1, fid, (i-1)*num_img+j);
            num_format = info.NumberFormat;
            if (isempty(num_format))
                num_format = 'unsigned integer';
            end
           
            % read the image
            img = read_img(fid, num_rows, num_cols, num_bytes, num_format);
            if isempty(X)
                X = img;
                
            % add image to image array
            else
                check_image_size(size(X), size(img));
                X(:,:,end+1) = img;
            end
            
        else
            fseek(fid, num_cols*num_rows*num_bytes, 'cof');
        end
    end
    
    i = i+1;
end



%%%
%%% Function get_dyna
%%%
function X = get_dyna(info, window, fid)

X = {};
i = 1;

while (i <= window)
    % get number of frame groups
    num_frames = get_val(info, 'NumberOfFrameGroups', 'number of frame groups', 1, fid, i);
    if (isempty(num_frames))
        num_frames = 1;
    end
    
    % read images in each frame group
    for j = 1:num_frames
        num_cols = get_val(info, 'MatrixSize1', 'matrix size [1]', 2, fid, (i-1)*num_frames+j);
        num_rows = get_val(info, 'MatrixSize2', 'matrix size [2]', 2, fid, (i-1)*num_frames+j);
        num_bytes = get_val(info, 'NumberOfBytesPerPixel', 'number of bytes per pixel', 2, fid, (i-1)*num_frames+j);
        num_img = get_val(info, 'NumberOfImagesThisFrameGroup', 'number of images this frame group', 2, fid, (i-1)*num_frames+j);
        % read images if input energy window
        if (i == window)
            num_format = get_val(info, 'NumberFormat', 'number format', 1, fid, (i-1)*num_frames+j);
            if (isempty(num_format))
                num_format = 'unsigned integer';
            end
            
            for k = 1:num_img
                % read the image
                img = read_img(fid, num_rows, num_cols, num_bytes, num_format);
                if isempty(X)
                    X = img;

                % add image to image array
                else
                    check_image_size(size(X), size(img));
                    X(:,:,end+1) = img;
                end
            end
            
        % otherwise increase offset
        else
            fseek(fid, num_img*num_cols*num_rows*num_bytes, 'cof');
        end
    end
    
    i = i+1;
end



%%%
%%% Function get_gated
%%%
function X = get_gated(info, window, fid)

X = {};
i = 1;

while (i<=window)
    num_cols = get_val(info, 'MatrixSize1', 'matrix size [1]', 2, fid, i);
    num_rows = get_val(info, 'MatrixSize2', 'matrix size [2]', 2, fid, i);
    num_bytes = get_val(info, 'NumberOfBytesPerPixel', 'number of bytes per pixel', 2, fid, i);
    num_format = get_val(info, 'NumberFormat', 'number format', 1, fid, i);
    if (isempty(num_format))
        num_format = 'unsigned integer';
    end
    
    num_time = get_val(info, 'NumberOfTimeWindows', 'number of time windows', 0, fid, i);
    if (isempty(num_time))
        num_time = 1;
    end

    % read images in time window or increase offset
    for j = 1:num_time
        num_img = get_val(info, 'NumberOfImagesInWindow', 'number of images in window', 2, fid, (i-1)*num_time+j);
        if (i == window)
            for k = 1:num_img
                % read the image
                img = read_img(fid, num_rows, num_cols, num_bytes, num_format);
                if isempty(X)
                    X = img;

                % add image to image array
                else
                    check_image_size(size(X), size(img));
                    X(:,:,end+1) = img;
                end
            end
            
        else
            fseek(fid, num_img*num_cols*num_rows*num_bytes, 'cof');
        end
    end
    
    i = i+1;
end



%%%
%%% Function get_tomo
%%%
function X = get_tomo(info, window, fid)

X = {};
i = 1;

while (i <= window)
    num_heads = get_val(info, 'NumberOfDetectorHeads', 'number of detector heads', 0, fid, i);
    if (isempty(num_heads))
        num_heads = 1;
    end

    % go through each detector head
    for j = 1:num_heads
        num_cols = get_val(info, 'MatrixSize1', 'matrix size [1]', 2, fid, (i-1)*num_heads+j);
        num_rows = get_val(info, 'MatrixSize2', 'matrix size [2]', 2, fid, (i-1)*num_heads+j);
        num_bytes = get_val(info, 'NumberOfBytesPerPixel', 'number of bytes per pixel', 2, fid, (i-1)*num_heads+j);
        % get process status
        status = get_val(info, 'ProcessStatus', 'process status', 1, fid, (i-1)*num_heads+j);
        if ((isempty(status)) || (strcmpi(status, 'reconstructed')))
            num_img = get_val(info, 'NumberOfSlices', 'number of slices', 2, fid, (i-1)*num_heads+j);
            
        elseif (strcmp(status, 'acquired'))
            num_img = get_val(info, 'NumberOfProjections', 'number of projections', 2, fid, (i-1)*num_heads+j);
            
        else
            fclose(fid);
            err_id = 'Images:interfileread:invalidValue';
            err_msg = 'Invalid value found at key ''process status'''.';
            error(err_id, err_msg);
        end
        
        % read the images based on process status or increase offset
        if (i == window)
            num_format = get_val(info, 'NumberFormat', 'number format', 1, fid, (i-1)*num_heads+j);
            if (isempty(num_format))
                num_format = 'unsigned integer';
            end
        
            for k = 1:num_img
                % read the image
                img = read_img(fid, num_rows, num_cols, num_bytes, num_format);
                if isempty(X)
                    X = img;

                % add image to image array
                else
                    check_image_size(size(X), size(img));
                    X(:,:,end+1) = img;
                end
            end
            
        else
            fseek(fid, num_img*num_cols*num_rows*num_bytes, 'cof');
        end
    end
    
    i = i+1;
end



%%%
%%% Function get_curve
%%%
function X = get_curve(info, fid)

num_cols = get_val(info, 'MatrixSize1', 'matrix size [1]', 2, fid);
num_rows = get_val(info, 'MatrixSize2', 'matrix size [2]', 2, fid);
num_bytes = get_val(info, 'NumberOfBytesPerPixel', 'number of bytes per pixel', 2, fid);
num_format = get_val(info, 'NumberFormat', 'number format', 1, fid);
if (isempty(num_format))
    num_format = 'unsigned integer';
end

X = read_img(fid, num_rows, num_cols, num_bytes, num_format);



%%%
%%% Function read_img
%%%
function img = read_img(fid, num_rows, num_cols, num_bytes, num_format)

% convert from Interfile to MATLAB data type
switch (lower(num_format))
    case ('unsigned integer')
        switch (num_bytes)
            case (1)
                num_type = 'uint8';
                
            case (2)
                num_type = 'uint16';
                
            case (4)
                num_type = 'uint32';
                
            case (8)
                num_type = 'uint64';
                
            otherwise
                fclose(fid);
                err_id = 'Images:interfileread:unsupportedFeature';
                err_msg = 'Number format ''ubit'' is not supported.';
                error(err_id, err_msg);
        end
        
    case ('signed integer')
        switch (num_bytes)
            case (1)
                num_type = 'int8';
                
            case (2)
                num_type = 'int16';
                
            case (3)
                num_type = 'int32';
                
            case (4)
                num_type = 'int64';
                
            otherwise
                fclose(fid);
                err_id = 'Images:interfileread:unsupportedFeature';
                err_msg = 'Number format ''bit'' is not supported.';
                error(err_id, err_msg);
        end
                
    case ('long float')
        num_type = 'float64';
        
%    case ('short float')
    case {'float', 'short float'} % (AMD added to deal with ADNI data)
        num_type = 'float32';
        
    case ('bit')
        num_type = 'ubit1';
        
    case ('ascii')
        num_type = 'char';
        
    otherwise
        fclose(fid);
        err_id = 'Images:interfileread:invalidValue';
        err_msg = 'Invalid value found at key ''number format''.';
        error(err_id, err_msg);
end

% read the image from file
img = reshape(fread(fid, num_rows*num_cols, num_type), num_cols, num_rows)';



%%%
%%% Function get_val
%%%
function val = get_val(varargin)

info = varargin{1};
key = varargin{2};
err_key = varargin{3};
required = varargin{4};
if (nargin >= 5)
    fid = varargin{5};
    if (nargin == 6)
        num = varargin{6};
    end

else
    num = [];
end

val = [];
err = 0;

% check to see if is a field
if (isfield(info,key))
    if (isempty(num))
        val = info.(key);
        
    else
        % check for over-indexing
        if (num > length(info.(key)))
            err = 2;
            
        else
            % read from cell
            if (iscell(info.(key)))
                val = info.(key){num};

            else
                % if string
                if (ischar(info.(key)))
                    if (num > 1)
                        err = 2;
                        
                    else
                        val = info.(key);
                    end
                    
                % if numeric
                else
                    val = info.(key)(num);
                end
            end
        end
    end
    
elseif (required)
    err = 1;
end

% set error if required but not specified
if ((required == 2) && (isempty(val)))
    err = 2;
end

if (err == 1)
    if (nargin >= 5)
        fclose(fid);
    end
    
    err_id = 'Images:interfileread:missingKey';
    err_msg = sprintf('Missing required key ''%s''.', err_key);
    error(err_id, err_msg);
    
elseif (err == 2)
    if (nargin >= 5)
        fclose(fid);
    end

    err_id = 'Images:interfileread:invalidValue';
    err_msg = sprintf('Invalid value found at key ''%s''.', err_key);
    error(err_id, err_msg);
end



%%%
%%% Function check_image_sizes
%%%
function check_image_size (orig_size, new_size)

% error if images are different sizes
if ~isequal(new_size, orig_size(1:2))
    err_id = 'Images:interfileread:diffImageSizes';
    err_msg = 'Images of different sizes in a window are not supported.';
    error(err_id, err_msg);
end
