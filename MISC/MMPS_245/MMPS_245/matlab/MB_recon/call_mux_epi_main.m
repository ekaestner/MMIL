function [  ] = call_mux_epi_main(pfile, dcmfile, outdir)
%UNTITLED Summary of this function goes here

dataFileHeader = read_MR_headers(pfile, 'all', 'raw');
nslices = dataFileHeader.rdb_hdr.nslices/dataFileHeader.rdb_hdr.reps;
ncoil = (dataFileHeader.rdb_hdr.dab(2)-dataFileHeader.rdb_hdr.dab(1))+1;
n_vcoils = ncoil/2; % Coil compression - number of virtual coils

p = mux_epi_params(pfile, [], [], n_vcoils, 1, '1Dgrappa', 1, 1, 1, []);

% -- Fix phase-encode direction for pepolar scans
if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
    pepolar = 1;
else
    pepolar = 0;
end

if pepolar == 1,
    d = mux_epi_main(pfile, [], [], [], [], n_vcoils, 1, '1Dgrappa', 1, 1, 1, dcmfile, dataFileHeader, 0, 0);
else
    d = mux_epi_main(pfile, [], [], [], 1, n_vcoils, 1, '1Dgrappa', 1, 1, 1, dcmfile, dataFileHeader, 0, 0);    
end

T2 = d(:,:,:,2);
maxim = real(max(T2(:)));
clear T2 d

% Get new series UID for dicom write
new_se_uid = get_new_se_uid(dcmfile);

% % Set up parallel processing
% if matlabpool('size') == 0
%     matlabpool('open',4);  %    matlabpool
% end
% parfor sl = 1:nslices
warning off;
for sl = 1:nslices % Non-parallel processing - use this for debugging in GUI
    fprintf('Processing slice #%d/%d...\n',sl,nslices);
    mux_epi_main(pfile, [], [], sl, [], n_vcoils, 1, '1Dgrappa', 1, 1, 1, dcmfile, dataFileHeader, maxim, new_se_uid, outdir);
end

% matlabpool close

end


