clear; clc;

prj_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/';
    red_cap_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'covariates' '/'];
    t1w_fsr_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'freesurfer' '/' 't1'];
    t1w_mca_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'micapipe' '/' 't1' '/'];
    t1w_c12_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'cat12' '/' 't1' '/'];

eng_dir = '/space/mcdonald-syn01/1/BIDS/enigma_conglom/';
    
mca_drv_dir = [ eng_dir '/' 'derivatives' '/' 'micapipe_v0.2.0' '/' ];
    t1w_mca_drv_txt = [  'ses-01' '/' 'anat' '/'];
    t1w_mca_drv_txt_no0 = [  'ses-1' '/' 'anat' '/'];
    t1w_mca_drv_txt_nse = [  'anat' '/'];
    t1w_mca_drv_fle     = '_ses-01_space-fsnative_T1w.nii.gz';
    t1w_mca_drv_fle_no0 = '_ses-1_space-fsnative_T1w.nii.gz';
    t1w_mca_drv_fle_nse = '_space-fsnative_T1w.nii.gz';

c12_drv_ver_one_dir = [ eng_dir '/' 'derivatives' '/' 'cat12' '/' ];
c12_drv_ver_two_dir = [ eng_dir '/' 'derivatives' '/' 'cat12_2' '/' ];

ses_hld = mmil_readtext([ prj_dir '/' 'enigma_session.csv']);

cnt_loc = '/space/mcdonald-syn01/1/software/containers/';
mni_loc = '/opt/micapipe/MNI152Volumes/MNI152_T1_0.8mm.nii.gz';

red_cap_fle = 'ENIGMAEpilepsy-CNN_DATA_2024-04-10_1209_comma.csv';

%% Redcap Codes
cde_bok.dx = { 0 'HC' ; ...
               1 'EPD' ; ...
               9 'other' };

cde_bok.sex = { 1 'Male' ; ...
                2 'Female' ; ...
                3 'other' ; ...
                4 'unknown'};

cde_bok.handedness = { 1 'Right' ; ...
                       2 'Left' ; ...
                       3 'Ambi' ; ...
                       4 'unknown'};

cde_bok.localization = { 1 'temporal' ; ...
                         2 'temporal_plus' ; ...
                         3 'ex-tle-frontal' ; ...
                         4 'unknown'};

cde_bok.lateralization = { 1 'left' ; ...
                           2 'right' ; ...
                           3 'bilateral' ; ...
                           4 'multifocal' ; ...
                           9 'unknown' };

cde_bok.mts = { 0 'noHS' ; ...
                1 'HS_mri' ; ...
                2 'HS_pathology' ; ...
                9 'unknown' };

%% Micapipe (fspro)
% eng_sbj_lst = dir(mca_drv_dir); eng_sbj_lst = { eng_sbj_lst(:).name }; eng_sbj_lst = eng_sbj_lst(string_find(eng_sbj_lst,'sub-'));
% err_lst = {};
% for iS = 1:numel(eng_sbj_lst)
%     inp_fle = [ mca_drv_dir '/' eng_sbj_lst{iS} '/' t1w_mca_drv_txt '/' eng_sbj_lst{iS} t1w_mca_drv_fle];
%     out_fle = [ t1w_mca_dir '/' eng_sbj_lst{iS} t1w_mca_drv_fle ];
%     if ~exist(out_fle,'file')
%         try
%             copyfile(inp_fle,out_fle);
%         catch
%             try
%                 inp_fle = [ mca_drv_dir '/' eng_sbj_lst{iS} '/' t1w_mca_drv_txt_nse '/' eng_sbj_lst{iS} t1w_mca_drv_fle_nse];
%                 copyfile(inp_fle,out_fle);
%             catch
%                 try
%                     inp_fle = [ mca_drv_dir '/' eng_sbj_lst{iS} '/' t1w_mca_drv_txt_no0 '/' eng_sbj_lst{iS} t1w_mca_drv_fle_no0];
%                     copyfile(inp_fle,out_fle);
%                 catch
%                     try
%                         ses_fld = dir([ mca_drv_dir '/' eng_sbj_lst{iS} '/']); ses_fld = ses_fld(3).name;
%                         inp_fle = [ mca_drv_dir '/' eng_sbj_lst{iS} '/' ses_fld '/' t1w_mca_drv_txt_nse '/' eng_sbj_lst{iS} '_' ses_fld t1w_mca_drv_fle_nse];
%                         copyfile(inp_fle,out_fle);
%                     catch
%                         err_lst{end+1} = sprintf('Missing: %s',inp_fle);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% tt1 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{900} t1w_mca_drv_fle ] );
% tt2 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{1700} t1w_mca_drv_fle ] );
% tt3 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{800} t1w_mca_drv_fle ] );
% tt4 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{2000} t1w_mca_drv_fle ] );
% subplot(2,2,1)
% imagesc( rot90(squeeze(tt1(:,:,70)),3) );
% subplot(2,2,2)
% imagesc( rot90(squeeze(tt2(:,:,70)),3) );
% subplot(2,2,3)
% imagesc( rot90(squeeze(tt3(:,:,70)),3) );
% subplot(2,2,4)
% imagesc( rot90(squeeze(tt4(:,:,70)),3) );

%% Micapipe (nativepro)
% Computers w/ singularity installed bigmac1, ip15
eng_sbj_lst = dir(mca_drv_dir); eng_sbj_lst = { eng_sbj_lst(:).name }; eng_sbj_lst = eng_sbj_lst(string_find(eng_sbj_lst,'sub-'));
err_lst = {};

cmd = sprintf(['singularity' ' ' 'exec' ' ' '-B /space' ' ' '%s' '/' '' 'micapipe_v0.2.2.sif' ' ' 'antsApplyTransforms'],cnt_loc);

for iS = 1:numel(eng_sbj_lst)

    sbj_cmd = sprintf([ '%s' ' ' '-d 3' ' ' '-i %s/%s/%s/anat/%s_%s_space-nativepro_T1w_brain.nii.gz' ],cmd,mca_drv_dir,eng_sbj_lst{iS},ses_hld{iS,2},eng_sbj_lst{iS},ses_hld{iS,2});
    sbj_cmd = sprintf([ '%s' ' ' '-r %s' ],sbj_cmd,mni_loc);
    sbj_cmd = sprintf([ '%s' ' ' '-t %s/%s/%s/xfm/%s_%s_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat'],sbj_cmd,mca_drv_dir,eng_sbj_lst{iS},ses_hld{iS,2},eng_sbj_lst{iS},ses_hld{iS,2});
    sbj_cmd = sprintf([ '%s' ' ' '-t %s/%s/%s/xfm/%s_%s_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz'],sbj_cmd,mca_drv_dir,eng_sbj_lst{iS},ses_hld{iS,2},eng_sbj_lst{iS},ses_hld{iS,2});
    sbj_cmd = sprintf([ '%s' ' ' '-o %s/%s_micapipe_mni152space.nii.gz'],sbj_cmd,t1w_mca_dir,eng_sbj_lst{iS});
    sbj_cmd = sprintf([ '%s' ' ' '-v'],sbj_cmd);

    unix(sbj_cmd)

end

tt1 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{1} '_micapipe_mni152space.nii.gz' ] );
tt2 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{2} '_micapipe_mni152space.nii.gz' ] );
tt3 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{3} '_micapipe_mni152space.nii.gz' ] );
tt4 = niftiread( [ t1w_mca_dir '/' eng_sbj_lst{4} '_micapipe_mni152space.nii.gz' ] );
subplot(2,2,1)
imagesc( rot90(squeeze(tt1(:,150,:))),  [0 40] ); colorbar;
subplot(2,2,2)
imagesc( rot90(squeeze(tt2(:,150,:))),  [0 40] ); colorbar;
subplot(2,2,3)
imagesc( rot90(squeeze(tt3(:,150,:))),  [0 40] ); colorbar;
subplot(2,2,4)
imagesc( rot90(squeeze(tt4(:,150,:))),  [0 40] ); colorbar;

%% CAT12 (Donatello)
c12_ver_one_sbj_lst = dir(c12_drv_ver_one_dir); c12_ver_one_sbj_lst = { c12_ver_one_sbj_lst(:).name }; c12_ver_one_sbj_lst = c12_ver_one_sbj_lst(string_find(c12_ver_one_sbj_lst,'sub-'));
c12_ver_two_sbj_lst = dir(c12_drv_ver_two_dir); c12_ver_two_sbj_lst = { c12_ver_two_sbj_lst(:).name }; c12_ver_two_sbj_lst = c12_ver_two_sbj_lst(string_find(c12_ver_two_sbj_lst,'sub-'));
err_lst = {};

tot_sbj_lst = unique( [ c12_ver_one_sbj_lst' ; c12_ver_two_sbj_lst' ]);

c12_dta_out = [ tot_sbj_lst cell(numel(tot_sbj_lst),3) ];

for iS = 1:numel(tot_sbj_lst)

    % c12_ver_one_sbj_lst %%%
    ses_dir = dir([ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} ]); ses_dir = {ses_dir(:).name};
    if any(string_find(ses_dir,'ses'))
        ses_dir = ses_dir{string_find(ses_dir,'ses')};
        mwp_pth = [ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} '/' ses_dir '/' 'anat' '/'];
        mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
        try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        
    else
        mwp_pth = [ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} '/' 'anat' '/'];
        mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
        try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
    end
    if ~isempty(mwp_fle)
        c12_dta_out{iS,2} = mwp_fle;
        c12_dta_out{iS,4} = [ mwp_pth '/' mwp_fle];
    else
        c12_dta_out{iS,2} = 'MISSING_mwp_fle';
    end

    % c12_ver_two_sbj_lst %%%
    if any(strcmpi(c12_ver_two_sbj_lst,tot_sbj_lst{iS}))
        ses_dir = dir([ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} ]); ses_dir = {ses_dir(:).name};
        if any(string_find(ses_dir,'ses'))
            ses_dir = ses_dir{string_find(ses_dir,'ses')};
            mwp_pth = [ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} '/' ses_dir '/' 'anat' '/'];
            mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
            try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        else
            mwp_pth = [ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} '/' 'anat' '/'];
            mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
            try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        end
        if ~isempty(mwp_fle)
            c12_dta_out{iS,3} = mwp_fle;
            if isempty(c12_dta_out{iS,4})
                c12_dta_out{iS,4} = [ mwp_pth '/' mwp_fle];
            else
                c12_dta_out{iS,5} = 'Double Files';
            end
        else
            c12_dta_out{iS,3} = 'MISSING_mwp_fle';
        end
    end
end

% Check the outputs %%%
clear sum_tbl
sum_tbl{1,2} = numel(string_find(c12_dta_out(:,2),'.nii')); sum_tbl{1,1} = 'derivatives/cat12 has mwp file';
sum_tbl{2,2} = numel(string_find(c12_dta_out(:,3),'.nii')); sum_tbl{2,1} = 'derivatives/cat12_2 has mwp file';

sum_tbl{3,2} = sum(strcmpi(c12_dta_out(:,2),'MISSING_mwp_fle') & cellfun(@isempty,c12_dta_out(:,3))); sum_tbl{3,1} = 'derivatives/cat12 missing mwp file';
sum_tbl{4,2} = sum(strcmpi(c12_dta_out(:,3),'MISSING_mwp_fle')); sum_tbl{4,1} = 'derivatives/cat12_2 missing mwp file';

% Save Out %%%
cell2csv([ prj_dir '/' 'CAT12_wmp1_file_locations_2024_05_22.csv'],[{'subject'} {'derivatives/cat12'} {'derivatives/cat12_2'} {'file_location'} {'notes'}; c12_dta_out])
cell2csv([ prj_dir '/' 'CAT12_wmp1_file_totals_2024_05_22.csv'],sum_tbl);

% Copy Files over %%%
for iS = 1:size(c12_dta_out,1)
    if ~isempty(c12_dta_out{iS,4}) && ~exist([ t1w_c12_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii'],'file') 
        copyfile(c12_dta_out{iS,4}, [ t1w_c12_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii']);
    end
end

%% Load ENIGMA Redcap
%
cov_var = { 'site' 'dx' 'sex' 'age' 'handedness' 'localization' 'lateralization' 'mts'};
sbj_nme_var = 'subjid';
sbj_nme_org = 'old_subjid';

%
fcfg = [];
fcfg.dta_loc = [red_cap_dir '/' red_cap_fle];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

reg_dta_row = [ 3:9 14:19];
img_dta_row = [ 1:2 10:13 ];

cmb_sbj = unique(red_cap_sbj);
cmb_col = red_cap_col;
cmb_dta = cell( numel(cmb_sbj), numel(red_cap_col) );
for iS = 1:numel(cmb_sbj)

    sbj_row = find(strcmpi(red_cap_sbj,cmb_sbj{iS}));
    if numel(sbj_row)>2; error('rows greater than 2'); end
    reg_row = sbj_row( ~strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );
    img_row = sbj_row( strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );

    cmb_dta(iS,reg_dta_row) = red_cap_dta(reg_row,reg_dta_row);
    if ~isempty(img_row); cmb_dta(iS,img_dta_row) = red_cap_dta(img_row,img_dta_row); end

    % fix ses BS
    ses_pos = string_find( cmb_dta(iS,strcmpi(red_cap_col,'subjid')), 'ses' );
    if ~isempty(ses_pos)
        ses_pos = strfind( cmb_dta{iS,strcmpi(red_cap_col,'subjid')}, 'ses' );
        cmb_dta{iS,strcmpi(red_cap_col,'subjid')}(ses_pos-1:end) = [];
    end

end

% Put together filenames & characteristics
cov_img_dta = cell(size(cmb_dta,1),4+numel(cov_var));
cov_img_col = [ {'sbj_nme'} {'img_nme'} {'t1w_c12_fle'} {'t1w_mca_fle'} cov_var ];

t1w_c12_fle = dir(t1w_c12_dir); t1w_c12_fle = {t1w_c12_fle(:).name}; t1w_c12_fle = t1w_c12_fle(string_find(t1w_c12_fle,".nii"));
t1w_mca_fle = dir(t1w_mca_dir); t1w_mca_fle = {t1w_mca_fle(:).name}; t1w_mca_fle = t1w_mca_fle(string_find(t1w_mca_fle,".nii"));

[~, cov_col ] = ismember(cov_var,cmb_col);

for iS = 1:size(cmb_dta,1)

    cov_img_dta{iS,1} = cmb_dta{iS,strcmpi(cmb_col,sbj_nme_org)};
    cov_img_dta{iS,2} = cmb_dta{iS,strcmpi(cmb_col,sbj_nme_var)};
    try cov_img_dta{iS,3} = t1w_c12_fle{string_find( t1w_c12_fle, cov_img_dta{iS,2} )}; catch cov_img_dta{iS,3} = ''; end
    try cov_img_dta{iS,4} = t1w_mca_fle{string_find( t1w_mca_fle, cov_img_dta{iS,2} )}; catch cov_img_dta{iS,4} = ''; end
    cov_img_dta(iS,5:end) = cmb_dta(iS,cov_col);

end

% Find Missing
mss_fle = sort(setxor(t1w_c12_fle, unique(cov_img_dta_rcd(:,3)) ));
cell2csv( [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'missing_enigma_imaging_and_covariates_2024_05_22.csv'],mss_fle);

% Fix covariates
cov_nme = fieldnames(cde_bok);
cov_val = cell(1,numel(cov_nme));
for iCV = 1:numel(cov_nme)
    cov_val{iCV} = cde_bok.(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
end

fcfg = [];

fcfg.dta     = cov_img_dta;
fcfg.dta_col = cov_img_col;

fcfg.rcd_nme = cov_nme;
fcfg.rcd_val = cov_val;

fcfg.swt_val = 1;

cov_img_dta_rcd = ejk_recode(fcfg);

%
[ sum(~strcmpi(cov_img_dta_rcd(:,3),'')) numel(t1w_c12_fle); 
  sum(~strcmpi(cov_img_dta_rcd(:,4),'')) numel(t1w_mca_fle) ]

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'dx'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'localization'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'lateralization'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'mts'))))

cell2csv( [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'enigma_imaging_and_covariates_2024_05_22.csv'],[ {'sbj_nme'} cov_img_col ; cmb_sbj cov_img_dta ])

%% Preliminary data
typ_dta = cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'dx')));
    typ_dta(cellfun(@isempty,typ_dta)) = {''};
loc_dta = cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'localization')));
    loc_dta(cellfun(@isempty,loc_dta)) = {''};
lat_dta = cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'lateralization')));
    lat_dta(cellfun(@isempty,lat_dta)) = {''};
c12_dta = ~cellfun(@isempty,cov_img_dta_rcd(:,3));
mca_dta = ~cellfun(@isempty,cov_img_dta_rcd(:,4));

% type
grp_nme = 'typ';
grp.(grp_nme).hc      = find( strcmpi(typ_dta,'HC')  & c12_dta==1 );
grp.(grp_nme).tle     = find( strcmpi(typ_dta,'EPD') & strcmpi(loc_dta,'temporal') & c12_dta==1);
grp.(grp_nme).non_tle = find( strcmpi(typ_dta,'EPD') & ~strcmpi(loc_dta,'temporal') & c12_dta==1);

% lateralization
grp_nme = 'lat';
grp.(grp_nme).lft_tle = find( strcmpi(typ_dta,'EPD') & strcmpi(loc_dta,'temporal') & strcmpi(lat_dta,'left') & c12_dta==1);
grp.(grp_nme).rgh_tle = find( strcmpi(typ_dta,'EPD') & strcmpi(loc_dta,'temporal') & strcmpi(lat_dta,'right') & c12_dta==1);
grp.(grp_nme).bil_tle = find( strcmpi(typ_dta,'EPD') & strcmpi(loc_dta,'temporal') & strcmpi(lat_dta,'bilateral') & c12_dta==1);

% Load Data


% Run VBM

%% Freesurfer
% # antsApplyTransforms \
% #     -d dimension \
% #     -i input_image \
% #     -r reference_image \
% #     -t freesurfer_to_nativepro_matrix \
% #     -t nativepro_to_MNI_warpfield.nii.gz \
% #     -t nativepro_to_MNI_affine_matrix \
% #     -o output \
% #     -v # verbose flag
% 
% singularity \
%     exec \
%     -B /space
%     ~/DockerImages/micapipe_v0.2.3.sif \
%     antsApplyTransforms \
%         -d 3 \
%         -i /space/mcdonald_syn01/1/BIDS/enigma_conglom/derivatives/fastsurfer/sub-123_ses-01/mri/T1.mgz \
%         -r /opt/micapipe/MNI152Volumes/MNI152_T1_0.8mm.nii.gz \
%         -t /space/mcdonald_syn01/1/BIDS/enigma_conglom/derivatives/micapipe_v0.2.0/sub-TESTUNAM37348/ses-01/xfm/sub-TESTUNAM37348_ses-01_from-fsnative_to_nativepro_T1w_0GenericAffine.mat \
%         -t /space/mcdonald_syn01/1/BIDS/enigma_conglom/derivatives/micapipe_v0.2.0/sub-TESTUNAM37348/ses-01/xfm/sub-TESTUNAM37348_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz \
%         -t /space/mcdonald_syn01/1/BIDS/enigma_conglom/derivatives/micapipe_v0.2.0/sub-TESTUNAM37348/ses-01/xfm/sub-TESTUNAM37348_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat \
%         -o /path/to/save/dir/T1_warp.nii.gz
%         -v



%% T1.mgz apply transforms





% 1st pass %%%
% tt1 = niftiread( [ t1w_c12_dir '/' c12_dta_out{100,1} '_' 'cat12' '_' 'mwp1.nii'] );
%   
% subplot(2,2,1)
% imagesc( rot90(squeeze(tt1(:,79,:)))); colorbar;
% subplot(2,2,2)
% imagesc( rot90(squeeze(tt2(:,79,:)))); colorbar;
% subplot(2,2,3)
% imagesc( rot90(squeeze(tt3(:,79,:)))); colorbar;
% subplot(2,2,4)
% imagesc( rot90(squeeze(tt4(:,79,:)))); colorbar;
