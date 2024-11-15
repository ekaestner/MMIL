clear; clc;

prj_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/';

c12_vbm_out_dir = [ prj_dir '/' 'vbm' '/' 'cat12' '/'];
mca_vbm_out_dir = [ prj_dir '/' 'vbm' '/' 'micapipe' '/'];

t1w_c12_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'cat12' '/' 't1' '/'];
t1w_mca_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'micapipe' '/' 't1' '/'];

%% Load Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'enigma_imaging_and_covariates_colnames.csv'];
[ cov_dta, cov_sbj, cov_col ] = ejk_dta_frm( fcfg );

%% Create Groups
typ_dta = cov_dta(:,find(strcmpi(cov_col,'dx')));
    typ_dta(cellfun(@isempty,typ_dta)) = {NaN};
    typ_dta = cell2mat(typ_dta);
loc_dta = cov_dta(:,find(strcmpi(cov_col,'localization')));
    loc_dta(cellfun(@isempty,loc_dta)) = {NaN};
    loc_dta = cell2mat(loc_dta);
lat_dta = cov_dta(:,find(strcmpi(cov_col,'lateralization')));
    lat_dta(cellfun(@isempty,lat_dta)) = {NaN};
    lat_dta = cell2mat(lat_dta);
c12_dta = ~cellfun(@isempty,cov_dta(:,find(strcmpi(cov_col,'c12_fle'))));
mca_dta = ~cellfun(@isempty,cov_dta(:,find(strcmpi(cov_col,'mca_fle'))));

% CAT12
grp_nme = 'c12_epd_con';
grp.(grp_nme).hc = find( typ_dta==0  & c12_dta==1 );
grp.(grp_nme).tle = find( typ_dta==1 & loc_dta==1 & c12_dta==1);
grp.(grp_nme).lft_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==2 & c12_dta==1);
grp.(grp_nme).rgh_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==1 & c12_dta==1);
grp.(grp_nme).bil_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==3 & c12_dta==1);

% micapipe
grp_nme = 'mca_epd_con';
grp.(grp_nme).hc = find( typ_dta==0  & mca_dta==1 );
grp.(grp_nme).tle = find( typ_dta==1 & loc_dta==1 & mca_dta==1);
grp.(grp_nme).lft_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==2 & mca_dta==1);
grp.(grp_nme).rgh_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==1 & mca_dta==1);
grp.(grp_nme).bil_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==3 & mca_dta==1);

%% Load CAT12 Data
grp_nme = 'c12_epd_con';
c12_dta_ind = unique( [ grp.(grp_nme).hc ; grp.(grp_nme).tle] );

c12_col = find(strcmpi(cov_col,'c12_fle'));

tic
c12_dta_mtx = nan(113,137,113,numel(c12_dta_ind));
for iS = 1:numel(c12_dta_ind)
    c12_dta_mtx(:,:,:,iS) = niftiread( [ t1w_c12_dir '/' cov_dta{ c12_dta_ind(iS), c12_col } ] );
end
toc

%% Run CAT12 VBM
% Update covariates & groups
c12_cov_dta = cov_dta(c12_dta_ind,:);
c12_cov_sbj = cov_sbj(c12_dta_ind,:);

% Update groups
typ_dta = c12_cov_dta(:,find(strcmpi(cov_col,'dx')));
    typ_dta(cellfun(@isempty,typ_dta)) = {NaN};
    typ_dta = cell2mat(typ_dta);
loc_dta = c12_cov_dta(:,find(strcmpi(cov_col,'localization')));
    loc_dta(cellfun(@isempty,loc_dta)) = {NaN};
    loc_dta = cell2mat(loc_dta);
lat_dta = c12_cov_dta(:,find(strcmpi(cov_col,'lateralization')));
    lat_dta(cellfun(@isempty,lat_dta)) = {NaN};
    lat_dta = cell2mat(lat_dta);

grp_nme = 'c12_epd_con';
grp_c12.(grp_nme).hc      = find( typ_dta==0  );
grp_c12.(grp_nme).tle     = find( typ_dta==1 & loc_dta==1 );
grp_c12.(grp_nme).lft_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==2 );
grp_c12.(grp_nme).rgh_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==1 );
grp_c12.(grp_nme).bil_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==3 );    

% Setup analyses
grp_nme = { 'c12_epd_con' };

grp_vbm = { { 'hc' 'lft_tle'} { 'hc' 'rgh_tle'} { 'hc' 'bil_tle'} };

% Run VBM
iG = 1;
iA = 3;

fcfg = [];
fcfg.dta_one = c12_dta_mtx(:,:,:,grp_c12.(grp_nme{iG}).(grp_vbm{iA}{1}));
fcfg.dta_two = c12_dta_mtx(:,:,:,grp_c12.(grp_nme{iG}).(grp_vbm{iA}{2}));
vbm_out = ejk_vbm(fcfg);

save([ c12_vbm_out_dir '/' 'vbm_' grp_nme{iG} '_' grp_vbm{iA}{1} '_VS_' grp_vbm{iA}{2} '.mat'],'vbm_out')

% Montage
avg_brn = nanmean(c12_dta_mtx(:,:,:,grp_c12.(grp_nme{iG}).hc),4);
avg_brn_cut = 0.15;

slc = 1+32:4:size(vbm_out.tvl,2)-44;

ylm_low = -8;
ylm_hgh = 8;

pvl_msk = .005;

fcfg = [];

fcfg.plt_dta = vbm_out.tvl;

fcfg.brn_dta = avg_brn;
fcfg.brn_cut_off = avg_brn_cut;

fcfg.pvl_dta = vbm_out.pvl;
fcfg.pvl_msk = pvl_msk;

fcfg.slc      = slc;

fcfg.ylm_low = ylm_low;
fcfg.ylm_hgh = ylm_hgh;

fcfg.plt_dir = c12_vbm_out_dir;
fcfg.plt_nme = [ 'montage' '_' grp_nme{iG} '_' grp_vbm{iA}{1} '_VS_' grp_vbm{iA}{2} ];

ejk_montage(fcfg)


%% Load Micapipe Data
grp_nme = 'mca_epd_con';
mca_dta_ind = unique( [ grp.(grp_nme).hc ; grp.(grp_nme).tle] );

mca_col = find(strcmpi(cov_col,'mca_fle'));

tic
mca_dta_mtx = nan(227,272,227,numel(mca_dta_ind));
for iS = 1714:numel(mca_dta_ind)
    mca_dta_mtx(:,:,:,iS) = niftiread( [ t1w_mca_dir '/' cov_dta{ mca_dta_ind(iS), mca_col } ] );
end
toc

%%
% Update covariates & groups
mca_cov_dta = cov_dta(mca_dta_ind,:);
mca_cov_sbj = cov_sbj(mca_dta_ind,:);

% Update groups
typ_dta = mca_cov_dta(:,find(strcmpi(cov_col,'dx')));
    typ_dta(cellfun(@isempty,typ_dta)) = {NaN};
    typ_dta = cell2mat(typ_dta);
loc_dta = mca_cov_dta(:,find(strcmpi(cov_col,'localization')));
    loc_dta(cellfun(@isempty,loc_dta)) = {NaN};
    loc_dta = cell2mat(loc_dta);
lat_dta = mca_cov_dta(:,find(strcmpi(cov_col,'lateralization')));
    lat_dta(cellfun(@isempty,lat_dta)) = {NaN};
    lat_dta = cell2mat(lat_dta);

grp_nme = 'mca_epd_con';
grp_mca.(grp_nme).hc      = find( typ_dta==0  );
grp_mca.(grp_nme).tle     = find( typ_dta==1 & loc_dta==1 );
grp_mca.(grp_nme).lft_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==2 );
grp_mca.(grp_nme).rgh_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==1 );
grp_mca.(grp_nme).bil_tle = find( typ_dta==1 & loc_dta==1 & lat_dta==3 );    

% Setup analyses
grp_nme = { 'mca_epd_con' };

grp_vbm = { { 'hc' 'lft_tle'} { 'hc' 'rgh_tle'} { 'hc' 'bil_tle'} };

% Run VBM
iG = 1;
iA = 2;

fcfg = [];
fcfg.dta_one = mca_dta_mtx(:,:,:,grp_mca.(grp_nme{iG}).(grp_vbm{iA}{1}));
fcfg.dta_two = mca_dta_mtx(:,:,:,grp_mca.(grp_nme{iG}).(grp_vbm{iA}{2}));
vbm_out = ejk_vbm(fcfg);

save([ mca_vbm_out_dir '/' 'vbm_' grp_nme{iG} '_' grp_vbm{iA}{1} '_VS_' grp_vbm{iA}{2} '.mat'],'vbm_out')

% Montage
avg_brn = nanmean(mca_dta_mtx(:,:,:,grp_mca.(grp_nme{iG}).hc),4);
avg_brn_cut = 0.15;

slc = 1+32:16:size(vbm_out.tvl,2)-44;

ylm_low = -8;
ylm_hgh = 8;

pvl_msk = .01;

fcfg = [];

fcfg.plt_dta = vbm_out.tvl;

fcfg.brn_dta = avg_brn;
fcfg.brn_cut_off = avg_brn_cut;

fcfg.pvl_dta = vbm_out.pvl;
fcfg.pvl_msk = pvl_msk;

fcfg.slc      = slc;

fcfg.ylm_low = ylm_low;
fcfg.ylm_hgh = ylm_hgh;

fcfg.plt_dir = mca_vbm_out_dir;
fcfg.plt_nme = [ 'montage' '_' grp_nme{iG} '_' grp_vbm{iA}{1} '_VS_' grp_vbm{iA}{2} ];

ejk_montage(fcfg)




