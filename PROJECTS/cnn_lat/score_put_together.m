clear; clc;

ovr_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/';

dem_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/';

dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/performance_final';

lft_dir = [ ovr_dir '/' 'JohnnyData' '/' 'hippocampal_volume' '/' 'ltle hippocampal vol'];
rgh_dir = [ ovr_dir '/' 'JohnnyData' '/' 'hippocampal_volume' '/' 'rtle hippocampal vol'];

seg_dir = [ ovr_dir '/' 'JohnnyData' '/' 'automatedQC' '/'];

%% Load Data
lft_sbj_nme = mmil_readtext([ ovr_dir '/' 'JohnnyData' '/' 'performance_final' '/' 'fn_left.csv' ]);
rgh_sbj_nme = mmil_readtext([ ovr_dir '/' 'JohnnyData' '/' 'performance_final' '/' 'fn_right.csv' ]);

fcfg = [];
fcfg.dta_loc = [ dem_loc '/' 'demographic_table.csv'];
[ dem_dta, dem_dta_sbj, dem_dta_col ] = ejk_dta_frm(fcfg);
dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')) = strcat('Rush',cellfun(@num2str,dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')),'uni',0));

seg_dta = mmil_readtext([ seg_dir '/' 'segmentation_scores_combined.csv' ]);
seg_dta(cellfun(@isnumeric,seg_dta(:,1)),1) = strcat('Rush',cellfun(@num2str,seg_dta(cellfun(@isnumeric,seg_dta(:,1)),1),'uni',0));

%% hippocampal volume
lft_fle = dir(lft_dir);
    lft_fle = {lft_fle(3:end).name};
rgh_fle = dir(rgh_dir);
    rgh_fle = {rgh_fle(3:end).name};
    rgh_fle(string_find(rgh_fle,'icv')) = [];

ord_hip = [ 1 2 3 4 6 5 ];
    
lft_hip_vol = cell(0); 
rgh_hip_vol = cell(0); 
for iF = 1:numel(ord_hip)
    lft_dta_hld = mmil_readtext( [ lft_dir '/' lft_fle{ord_hip(iF)}] );
    lft_hip_vol = [ lft_hip_vol ; repmat({lft_fle{ord_hip(iF)}},size(lft_dta_hld,1),1) lft_dta_hld ];
    
    rgh_dta_hld = mmil_readtext( [ rgh_dir '/' rgh_fle{ord_hip(iF)}] );
    rgh_hip_vol = [ rgh_hip_vol ; repmat({rgh_fle{ord_hip(iF)}},size(rgh_dta_hld,1),1) rgh_dta_hld ];
end

hip_vol = [ lft_hip_vol ; rgh_hip_vol ];
cell2csv( [ dta_loc '/' 'hippocamapal_volume.csv'], [ {'sbj_nme' 'file' 'LeftHippocampalVolume' 'RighHippocampaltVolume'}  ; [ lft_sbj_nme ; rgh_sbj_nme] hip_vol ] )

% Re-order to match demographics\
hip_sbj_nme = [ lft_sbj_nme ; rgh_sbj_nme];
hip_vol_ord = cell(size(dem_dta_sbj,1),3);
for iS = 1:size(hip_vol_ord,1)
    hip_ind = strcmpi(hip_sbj_nme,dem_dta_sbj{iS});
    if sum(hip_ind)>0
        hip_vol_ord(iS,:) = hip_vol(hip_ind,:);
    else
        hip_vol_ord(iS,:) = {'' NaN NaN};
    end
end

cell2csv( [ dta_loc '/' 'hippocamapal_volume_ordered.csv'], [ {'sbj_nme' 'file' 'LeftHippocampalVolume' 'RighHippocampaltVolume'}  ; dem_dta_sbj hip_vol_ord ] )

%% AutomatedQC
seg_ord = cell(size(dem_dta_sbj,1),2);
for iS = 1:size(seg_ord,1)
    seg_ind = strcmpi(seg_dta(:,1),dem_dta_sbj{iS});
    if sum(seg_ind)>0
        seg_ord(iS,:) = seg_dta(seg_ind,:);
    else
        seg_ord(iS,:) = {'' NaN};
    end
end

cell2csv( [ dta_loc '/' 'automatedQC_ordered.csv'], [ {'sbj_nme' 'QC_Score'}  ; dem_dta_sbj seg_ord(:,2) ] );





