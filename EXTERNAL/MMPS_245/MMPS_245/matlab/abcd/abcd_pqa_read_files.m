function [vol4D, meta, fwhm] = abcd_pqa_read_files(input, output)
%PHANTOM_FMRI Summary of this function goes here
%
%   Reads dicoms from fMRI phantom QC
%
%   Detailed explanation goes here

filelist = dir(input);
filelist = filelist(~[filelist.isdir] & ~strncmpi('.', {filelist.name}, 1));
fullfilename = fullfile(input,filelist(1).name);
if file_is_dicom(fullfilename)
    info=dicominfo(fullfilename);
else
    fprintf('Only dicoms should be included in this folder. Aborting...\n');
    return;
end

if (strfind(info.Manufacturer, 'GE'))
    [meta] = read_ge_phantom(input,filelist);
elseif  (strfind(info.Manufacturer, 'SIEMENS'))
    [meta] = read_siemens_phantom(input,filelist);
elseif (strfind(info.Manufacturer, 'Philips'))
    [meta] = read_philips_phantom(input,filelist);
else
    return;
end

[vol4D, fwhm] = getFWHM(input,output, meta);

end    

function [meta]=read_siemens_phantom(input,filelist)
osLevel = '';
json_fname = strcat(input,'.json');
if exist(json_fname,'file') 
    json_data=loadjson(json_fname);
    if isfield(json_data,'OSLevel'), osLevel = json_data.OSLevel; end   
end

fullfilename = fullfile(input,filelist(1).name);
if file_is_dicom(fullfilename)
    info=dicominfo(fullfilename);
    imatrix = nonzeros(info.AcquisitionMatrix);   
    nVx = double(imatrix(1));
    nVy = double(imatrix(2));
    nSlices = double(info.Private_0019_100a);
    nFrames = length(filelist);
    
    imageFreq = info.ImagingFrequency;
    transmitGain = 0;
    aRecGain = 0;
    
    sx = info.PixelSpacing(1);
    sy = info.PixelSpacing(2);
    sz = info.SliceThickness;
    
    TR = info.RepetitionTime;
    FA = info.FlipAngle;
    TE = info.EchoTime;
    
    s_date = info.StudyDate;
    s_time = info.StudyTime;
    si_UID = info.StudyInstanceUID;
    se_time = info.SeriesTime;
    se_number = info.SeriesNumber; 
    manufact = info.Manufacturer;
    model = info.ManufacturerModelName;
    sDes = info.SeriesDescription;
    
    if (isfield(info,'DeviceSerialNumber')) serialNumber = info.DeviceSerialNumber; else serialNumber=''; end;   
    if (isfield(info,'SoftwareVersion')) softVersion = info.SoftwareVersion; else softVersion=''; end;
    if (isfield(info,'ImageComments')) imComments = info.ImageComments; else imComments=''; end;
    
    if (isfield(info, 'Private_0051_100f'))
        coilTypes = info.Private_0051_100f;
        if any(strfind(coilTypes, 'HEA;HEP'))
            coil = '32Ch';
        elseif any(strfind(coilTypes, 'HC1-7'))
            coil = '64Ch';
        else
            coil = 'Coil Not Recognized';
        end
    else
        coil = 'Error reading coil dicom tag';
    end
else
    fprintf('Only dicoms should be included in this folder. Aborting...');
    return;
end

meta = struct('TR',TR, 'FA',FA, 'TE', TE, 'imageFreq', imageFreq, 'transmitGain', transmitGain, 'aRecGain', aRecGain, 'sx', sx, 'sy', sy, 'sz', sz,...,
    's_date', s_date, 's_time', s_time, 'si_UID', si_UID, 'se_time', se_time, 'se_number', se_number, 'manufact', manufact, 'model', model, 'sDes', sDes, 'coil', coil,...,
    'serialNumber', serialNumber, 'softVersion', softVersion, 'osLevel', osLevel, 'imComments', imComments);
end


function [meta]=read_ge_phantom(input,filelist)
osLevel = '';
json_fname = strcat(input,'.json');
if exist(json_fname,'file') 
    json_data=loadjson(json_fname);
    if isfield(json_data,'OSLevel'), osLevel = json_data.OSLevel; end   
end
fullfilename = fullfile(input,filelist(1).name);
if(file_is_dicom(fullfilename))
    info=dicominfo(fullfilename);
    nVy = info.Rows;
    nVx = info.Columns;
    nImages = info.ImagesInAcquisition;
    nFrames = info.NumberOfTemporalPositions;
    
    sx = info.PixelSpacing(1);
    sy = info.PixelSpacing(2);
    sz = info.SliceThickness;
    
    imageFreq = info.Private_0019_1093;
    transmitGain = info.Private_0019_1094;
    aRecGain = info.Private_0019_1095;
    if (nImages > nFrames)
        while(mod(nImages,nFrames)~=0)
            nFrames = nFrames-1; 
        end
        nSlices = nImages/nFrames;
    else
        nSlices = nImages;
    end
    TR = info.RepetitionTime;
    FA = info.FlipAngle;
    TE = info.EchoTime;
    
    s_date = info.StudyDate;
    s_time = info.StudyTime;
    si_UID = info.StudyInstanceUID;
    se_time = info.SeriesTime;
    se_number = info.SeriesNumber; 
    manufact = info.Manufacturer;
    model = info.ManufacturerModelName;
    sDes = info.SeriesDescription;
    
    if (isfield(info,'DeviceSerialNumber')) serialNumber = info.DeviceSerialNumber; else serialNumber=''; end;   
    if (isfield(info,'SoftwareVersion')) softVersion = info.SoftwareVersion; else softVersion=''; end;
    if (isfield(info,'ImageComments')) imComments = info.ImageComments; else imComments=''; end;
       
    if (isfield(info, 'ReceiveCoilName'))
        coil = info.ReceiveCoilName;
    else
        coil = '';
    end
    
else
    fprintf('Only dicoms should be included in this folder. Aborting...');
    return;
end

meta = struct('TR',TR, 'FA',FA, 'TE', TE, 'imageFreq', imageFreq, 'transmitGain', transmitGain, 'aRecGain', aRecGain, 'sx', sx, 'sy', sy, 'sz', sz,...,
    's_date', s_date, 's_time', s_time, 'si_UID', si_UID, 'se_time', se_time, 'se_number', se_number, 'manufact', manufact, 'model', model, 'sDes', sDes, 'coil', coil,...,
    'serialNumber', serialNumber, 'softVersion', softVersion, 'osLevel', osLevel, 'imComments', imComments);

end


function [meta]=read_philips_phantom(input,filelist)
osLevel = '';
json_fname = strcat(input,'.json');
if exist(json_fname,'file') 
    json_data=loadjson(json_fname);
    if isfield(json_data,'OSLevel'), osLevel = json_data.OSLevel; end   
end
fullfilename = fullfile(input,filelist(1).name);
if(file_is_dicom(fullfilename))
    info=dicominfo(fullfilename);
    nVy = info.Rows;
    nVx = info.Columns;
    nImages = length(filelist);
    nFrames = double(info.NumberOfTemporalPositions);
    
    sx = info.PixelSpacing(1);
    sy = info.PixelSpacing(2);
    sz = info.SliceThickness;
    
    imageFreq = info.ImagingFrequency;
    transmitGain = 0;
    aRecGain = 0;
    
    nSlices = nImages/nFrames;
    TR = info.RepetitionTime;
    FA = info.FlipAngle;
    TE = info.EchoTime;
    
    
    s_date = info.StudyDate;
    s_time = info.StudyTime;
    si_UID = info.StudyInstanceUID;
    se_time = info.SeriesTime;
    se_number = info.SeriesNumber; 
    manufact = info.Manufacturer;
    model = info.ManufacturerModelName;
    sDes = info.SeriesDescription;

    if (isfield(info,'DeviceSerialNumber')) serialNumber = info.DeviceSerialNumber; else serialNumber=''; end;   
    if (isfield(info,'SoftwareVersion')) softVersion = info.SoftwareVersion; else softVersion=''; end;
    if (isfield(info,'ImageComments')) imComments = info.ImageComments; else imComments=''; end;
    
    if (isfield(info, 'ReceiveCoilName'))
        coil = info.ReceiveCoilName;
    else
        coil = '';
    end
    
else
    fprintf('Only dicoms should be included in this folder. Aborting...');
    return;
end

meta = struct('TR',TR, 'FA',FA, 'TE', TE, 'imageFreq', imageFreq, 'transmitGain', transmitGain, 'aRecGain', aRecGain, 'sx', sx, 'sy', sy, 'sz', sz,...,
    's_date', s_date, 's_time', s_time, 'si_UID', si_UID, 'se_time', se_time, 'se_number', se_number, 'manufact', manufact, 'model', model, 'sDes', sDes, 'coil', coil,...,
    'serialNumber', serialNumber, 'softVersion', softVersion, 'osLevel', osLevel, 'imComments', imComments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol4D, fwhm] = getFWHM(input, output, meta)

filelist = dir(input);
filelist = filelist(~[filelist.isdir] & ~strncmpi('.', {filelist.name}, 1));
i=1;
fullfilename = fullfile(input,filelist(i).name);
while (~file_is_dicom(fullfilename)) && (i<length(filelist))
    fullfilename = fullfile(input,filelist(i).name);
end

fname=fullfile(output,'original_vol.nii');

if ~exist(fname, 'file')
    if any(strfind(meta.manufact, 'SIEMENS'))
        cmd = sprintf('mri_convert -it siemens_dicom -ot nii %s %s/original_vol.nii', fullfilename, output);
    else
        cmd = sprintf('mri_convert -it dicom -ot nii %s %s/original_vol.nii', fullfilename, output);
    end
    unix(cmd);
end

nifti_image = load_nii(fname);
vol4D = rot90(flip(double(nifti_image.img),1),3);

cmd = sprintf('mkdir %s/AFNI', output);
unix(cmd);
cmd = sprintf('3dcopy %s %s/AFNI/dset+orig', fname, output);
unix(cmd);

afni_output = [output, '/AFNI'];
cmd = sprintf('3dvolreg -prefix %s/volreg %s/dset+orig',afni_output, afni_output);
unix(cmd);
cmd = sprintf('3dDetrend -polort 2 -prefix %s/voldetrend %s/volreg+orig', afni_output, afni_output);
unix(cmd);
cmd = sprintf('3dTstat -mean -prefix %s/volmean %s/volreg+orig', afni_output, afni_output);
unix(cmd);
cmd = sprintf('3dAutomask -q -prefix %s/volmask %s/volmean+orig', afni_output, afni_output);
unix(cmd);
cmd = sprintf('3dFWHMx -dset %s/voldetrend+orig -mask %s/volmask+orig -out %s/FWHMVALS', afni_output, afni_output, afni_output);
unix(cmd);

fname = fullfile(afni_output,'FWHMVALS');
fileID = fopen(fname,'r');
formatSpec = '%f';
sizeA = [3 size(vol4D,4)];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);

fwhm_x = A(1,:);
fwhm_x(fwhm_x==-1)=0;
fwhm(1)=mean(nonzeros(fwhm_x));

fwhm_y = A(2,:);
fwhm_y(fwhm_y==-1)=0;
fwhm(2)=mean(nonzeros(fwhm_y));

fwhm_z = A(3,:);
fwhm_z(fwhm_z==-1)=0;
fwhm(3)=mean(nonzeros(fwhm_z));

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isdicom=file_is_dicom(filename)
isdicom=false;
try
    fid = fopen(filename, 'r');
    status=fseek(fid,128,-1);  
    if(status==0)
        tag = fread(fid, 4, 'uint8=>char')';
        isdicom=strcmpi(tag,'DICM');
    end
    fclose(fid);
catch me
end

end