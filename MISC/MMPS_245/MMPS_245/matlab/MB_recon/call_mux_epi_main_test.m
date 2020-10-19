function [  ] = call_mux_epi_main(pfile, dcmfile, totalN)
%UNTITLED Summary of this function goes here

% convert to number
totalN = str2num(totalN);

dataFileHeader = read_MR_headers(pfile, 'all', 'raw');
nslices = dataFileHeader.rdb_hdr.nslices/dataFileHeader.rdb_hdr.reps;
ncoil = (dataFileHeader.rdb_hdr.dab(2)-dataFileHeader.rdb_hdr.dab(1))+1;
n_vcoils = ncoil/2; % Coil compression - number of virtual coils

% Call once for first timepoint to set scaling factor for dicomwrite
d = mux_epi_main(pfile, [], [], [], 1, n_vcoils, 1, '1Dgrappa', 1, 1, 1, dcmfile, dataFileHeader, 0, 0);
dataFileName = 'dataFile.mat';
save(sprintf('%s', dataFileName), 'dataFileHeader')
T2 = d(:,:,:,2);
maxim = max(T2(:));
clear T2 d
% Get new series UID for dicom write
new_se_uid = get_new_se_uid(dcmfile);

% % Set up parallel processing
% if matlabpool('size') == 0
%     matlabpool('open',4);  %    matlabpool
% end
% parfor sl = 1:nslices
tic
%nvols = dataFileHeader.rdb_hdr.pasframe;
% outside parallelism using totalN instances, this one is instance thisN
chunk = floor(nslices/totalN);
for z = 1:totalN
startSlice(z) = chunk * (z-1) + 1;
endSlice(z)   = chunk * ((z-1) + 1);
  if z == totalN,
      endSlice(z) = nslices;
  end;
  %cmd = sprintf('matlab -nodesktop  -nosplash -r "sliceProcess(%d, %d, %d, %s, [dataFileHeader], %d, %s)" &', startSlice(z), endSlice(z), n_vcoils, dcmfile, maxim, new_se_uid);
  cmd = sprintf('matlab -nodesktop -nosplash -r "sliceProcess(%d, %d, %d, ''%s'', ''%s'', ''%s'', %d, ''%s'')" &', startSlice(z), endSlice(z), n_vcoils, pfile, dcmfile,dataFileName, maxim, new_se_uid);
  %cmd
  system([cmd]);
  % matlabpool close

end


toc
