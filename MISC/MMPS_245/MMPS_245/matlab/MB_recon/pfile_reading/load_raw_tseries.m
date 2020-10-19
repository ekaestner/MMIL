function [data, params] = load_raw_tseries(fname,params,slices,passes)
%
% Modified from rawloadX.m, to load k-space data of a time series conntained
% in a p-file.
%
% function [data, params] = load_raw_tseries(fname,params,slices,passes);
%
%	Function reads selected (or all) data from a P-file,
%	using given lists of slices and/or passes (i.e., timepoints).
%	Reads 12.x, 14.x and 20.x EPIC files.
%
%	INPUT:
%		fname   = path and file name.
%       params  = The structure with all parameters of the scan and
%                 for the reconstruction. Please refer to mux_epi_params.m
%                 for the definition of each field in this structure.
%                 The 'pfile_header' field in 'params' is used in
%                 this function. If 'params.pfile_header' is empty, the
%                 header information will be read using function 'rawheadX.m'.
%		slices  = slice numbers to read (1...).
%       passes  = pass numbers to read (1...). For the fMRI scan,
%                 this must be continuous from 1 to the last pass
%                 number to read. For an external calibration scan,
%                 this should correspond to the last group of slice
%                 phase cycling.

%
%	OUTPUT:
%		data = data array (up to 6 dimensions, including all
%			specified frames, echoes, slices, coils, passes.
%               params = The output parameter structure, having the same
%                       fields as the input 'params'. The following fields might
%                       have been changed inside this function if the scan was aborted
%                       midway and the p-file was only partially collected: num_mux_cycle,
%                       num_passes, nt_to_recon, tpoints_to_load.
%
%
%	NOTES:
%	1) If ranges are beyond those of the header, then they are
%		limited by ranges in the header.
%	2) If not all expected frames are read, then the data are returned
%		as a 2D array, with all frames read.
%
%
%	B.Hargreaves -- April 2008.
%       Modified by Kangrong Zhu,   Stanford University       July 2012

% -- Header information
if isempty(params.pfile_header)
    header = rawheadX(fname);
else
    header = params.pfile_header;
end

% -- Set defaults, if not provided.
if (~exist('slices','var') || isempty(slices))
    slices = [1:header.nslices/header.npasses];
end;
if (~exist('passes','var') || isempty(passes))
    passes = [1:header.npasses];
end;

rds_data_order = params.rds_data_order;
inplane_R = params.inplane_R;

% -- Open file.
d = dir(fname);
file_bytes = d.bytes;
nframes = header.nframes+header.hnover;
nechos = header.nechoes;
ncoils  = header.ncoils;
num_acquired_slices = params.num_slices;

slices = checkrange(slices,0,num_acquired_slices,'Slices');
passes = checkrange(passes,0,header.npasses,'Passes');

if (passes(1)==0) && (passes(end) ~= length(passes)-1) % Assuming num_mux_cycle is always larger than 1, so the p-file must be from an fMRI scan, not from an external calibration scan, when passes(1)==1
    error('The pass numbers to read must be continuous from 1 to the last pass number to read.');
end

% --- Allocate array for data, based on passed arguments or defaults.
data = zeros(header.frsize, nframes, nechos, length(slices), ncoils, length(passes));

ptsize = header.ptsize;			% Sample size (bytes)
rawhdrsize = header.rawhdrsize;

aborted_scan = false;                   % true: the scan was aborted midway and the k-space data was collected only partially.

dattypes = {'int16','int32'};

if(rds_data_order)
    % It took a year, but I think I finally cracked the flip even/odd thing! Time will tell as the scan onslaught rages on...
    % The important rule seems to be: if the number of acquired echos is even, flip even frames. If it's odd, flip odd frames.
    % And note that the number of acquired echos is nframes/inplane_R.
    if(mod(nframes/inplane_R, 2))
        flip_even_frames = false;
    else
        flip_even_frames = true;
    end
    fip = fopen(fname,'r','l');
    if fip == -1
      error('File %s not found\n',fn);
    end
    % The following hack is no longer needed-- the rdsclient ensures that
    % the header size is set correctly.
    % if params.pfile_header.version<20.0065
    %     hdr_size = 149800;
    % else
    %     hdr_size = 149808;
    % end
    % if(rawhdrsize ~= hdr_size) %file reports 149788, ese22 was 149800, file size suggests 149808 for ese23
    %     fprintf('  Fixing raw header size.\n');
    %     rawhdrsize = hdr_size;
    % end
    
    % Data are arranged (from fastest to slowest moving):
    % coil - echo - slice - pass (timeframe)
    % note that in here a "frame" refers to a line of kspace (freq-encode line).
    framesize = header.frsize;
    viewsize = ncoils * nechos * framesize;
    viewbytes = viewsize * 2 * ptsize;
    % Size of one slice (bytes)
    slicebytes = viewbytes * (header.nframes+header.hnover)/inplane_R;
    passbytes = slicebytes * num_acquired_slices;      % Size of one pass (bytes)
    expected_num_bytes = params.pfile_header.npasses * passbytes + params.mux*(inplane_R-1)*params.num_mux_cycle * passbytes + rawhdrsize;
    missing_bytes = expected_num_bytes - file_bytes;
    start_pass = 1;
    if(missing_bytes > 0)
        fprintf('  File is missing %d bytes! Will assume they are missing from the end\n', missing_bytes);
    elseif(missing_bytes<0)
        % This really shouldn't happen. Maybe we should error out?
        fprintf('  File has %d extra bytes! Not sure what to do. Maybe this is a regular p-file?\n', abs(missing_bytes));
        fclose(fip);
        rds_data_order = false;
    end
end
% we might have decided that this isn't an rds file after counting bytes...
if(rds_data_order)
    passes = passes - 1;
    slices = slices - 1;
    for(passindex = start_pass:length(passes))
        if(passes(passindex) < params.mux*params.num_mux_cycle) && (inplane_R > 1)
            % These are calibration scans.
            acquired_pass = passes(passindex) * inplane_R;
            for(shot = 0:inplane_R-1)
                for(sliceindex = 1:length(slices))
                    slices_act = slices;
                    if (params.slice_shuffling == 1)
                        slices_act = mod(slices - passindex + 1, params.num_slices);
                    end
                    fseek(fip, rawhdrsize + (acquired_pass + shot)*passbytes + slices_act(sliceindex)*slicebytes, -1);
                    [d,nread] = fread(fip, slicebytes/ptsize, dattypes{ptsize/2});
                    if (nread == slicebytes/ptsize)
                        d = d(1:2:end) + 1i*d(2:2:end);
                        for frame = 0:inplane_R:nframes-1
                            for echo = 0:nechos-1
                                for coil = 0:ncoils-1
                                    startind = coil*framesize + echo*ncoils*framesize + frame/inplane_R*viewsize;
                                    endind = startind + framesize;
                                    even = mod(floor(frame/inplane_R),2)==0;
                                    if((even && flip_even_frames) || (~even && ~flip_even_frames))
                                        data(1:end, nframes-frame-shot, echo+1, sliceindex, coil+1, passindex) = d(endind:-1:startind+1);
                                    else
                                        data(1:end, nframes-frame-shot, echo+1, sliceindex, coil+1, passindex) = d(startind+1:endind);
                                    end
                                end
                            end
                        end
                    else
                        aborted_scan = true;
                        acquired_num_passes = passindex - 1;    % Assuming the data until the previous pass has been collected correctly.
                        break;                                  % Once we know the p-file was only collected partially, we don't need to read in any more data.
                    end;
                end;  % --Slice loop.
                if aborted_scan
                    break;
                end
            end % -- Shot loop
            if aborted_scan
                break;
            end
        else
            % Done with calibration scans-- load the normal scans here.
            acquired_pass = passes(passindex) + params.mux*(inplane_R-1)*params.num_mux_cycle;
            for(sliceindex = 1:length(slices))
                slices_act = slices;
                if (params.slice_shuffling == 1)
                    slices_act = mod(slices - passindex + 1, params.num_slices);
                end
                % Seek to the correct slice.
                fseek(fip, rawhdrsize + acquired_pass*passbytes + slices_act(sliceindex)*slicebytes, -1);
                % Read the entire slice (all coils, all echos, all frames)
                [d,nread] = fread(fip, slicebytes/ptsize, dattypes{ptsize/2});
                if (nread == slicebytes/ptsize)
                    % convert to complex
                    d = d(1:2:end) + 1i*d(2:2:end);
                    % We only get the actual acquired data. We'll fill the acquired lines and leave the rest zero-padded.
                    for frame = 0:inplane_R:nframes-1
                        for echo = 0:nechos-1
                            for coil = 0:ncoils-1
                                % EPI acquires data in a raster pattern, so we need to flip every other FE line.
                                startind = coil*framesize + echo*ncoils*framesize + frame/inplane_R*viewsize;
                                endind = startind + framesize;
                                even = mod(floor(frame/inplane_R),2)==0;
                                if((even && flip_even_frames) || (~even && ~flip_even_frames))
                                    data(:, nframes-frame-(inplane_R-1), echo+1, sliceindex, coil+1, passindex) = d(endind:-1:startind+1);
                                else
                                    data(:, nframes-frame-(inplane_R-1), echo+1, sliceindex, coil+1, passindex) = d(startind+1:endind);
                                end
                            end
                        end
                    end
                else
                    aborted_scan = true;
                    acquired_num_passes = passindex - 1;    % Assuming the data untill the previous pass has been collected correctly.
                    break;                                  % Once we know the p-file was only collected partially, we don't need to read in any more data.
                end;
            end;  % --Slice loop.
            if aborted_scan
                break;
            end
        end % --if(passes...
    end % --Pass loop
    % -- No need to skip remaining passes, as no more data to read!
    fclose(fip);	% -- Close file.
    if params.debug
        fprintf('  RDS-format data loaded.\n');
    end
else
    % TODO: rewrite this code here. It can be simplified. Also, we should load the data into a more compute-friendly
    % array order so that when we start the calculations we aren't thrashing all over RAM.
    [data, params] = rawloadX_tseries(fname,[],[],slices,[],passes,params);
end


% -- Deal with partially acquired p-files.
if aborted_scan
    if params.debug
        fprintf('  The scan was aborted midway and the p-file was only collected partially.\n');
    end

    if ~exist('params', 'var')
        error('The parameter structure doesn''t exist.');
    end

    if acquired_num_passes < 1          % p-file contains no data.
        error('The p-file doesn''t contain any data.');
    end

    if acquired_num_passes < params.mux % Data in p-file is less than one mux phase cycling.
        error(['The number of temporal frames contained in the p-file is %d', ...
            'mux is %d. Can''t conduct any reconstruction in this case.'], acquired_num_passes, params.mux);
    end

    if acquired_num_passes <= params.mux*params.num_mux_cycle; % Data in p-file are all mux phase cycling frames.
        params.num_mux_cycle = floor(acquired_num_passes/params.mux);
        params.num_passes = params.mux*params.num_mux_cycle;
        params.nt_to_recon = 0;
        params.tpoints_to_load = 1 : params.num_passes;
        data = data(:, :, :, :, :, params.tpoints_to_load);
        if params.debug
            fprintf(['  Number of frames contained in the p-file is %d, mux is %d, num_mux_cycle is %d.', ...
                'Only %d calibration frames will be reconstructed.\n'], ...
                acquired_num_passes, params.mux, params.num_mux_cycle, params.num_passes);
        end
    end

    if acquired_num_passes > params.mux*params.num_mux_cycle; % Some accelerated frames are contained in the p-file.
        params.num_passes = acquired_num_passes;
        params.nt_to_recon = acquired_num_passes - params.mux*params.num_mux_cycle;
        params.tpoints_to_load = 1 : min(params.num_passes, size(data,6));
        data = data(:, :, :, :, :, params.tpoints_to_load);
        if params.debug
            fprintf('  Only the acquired %d frames will be loaded and reconstructed.\n', params.num_passes);
        end
    end
end

params_out = params;


return;		% -- End of main function.



function arrout = checkrange(arr,amin,amax,atype)
%%%% -- Check values in array are within range, and remove ones that
	%%%are not!
f = find(arr >= amin);
arrout = arr(f);
f = find(arrout <= amax);
arrout = arrout(f);

if (length(arrout) < length(arr))
  arr = arr(:);
  f1 = find(arr < amin);
  f2 = find(arr > amax);
  tt = sprintf('%d %s out of range: ',length(f1)+length(f2),atype); disp(tt);
  disp([arr(f1); arr(f2)].');
end;
return;
%%%% --------


function [fileptr,counter] = skipsection(fileptr,counter,target,secsize)
% Skip ahead so that counter matches target, in chunks of secsize bytes.

secskip = target-counter;
fseek(fileptr,secskip*secsize,0);
counter = target;
% ---- End;


