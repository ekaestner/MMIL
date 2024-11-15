clear; clc;

dta_dir = '/space/mcdonald-syn01/1/projects/ekaestner/Diffusion_T1_forErik_12.15.23/';

plt_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/';

grp_dir = { 'Processed-Controls' 'Processed-Patients' };
grp_nme = { 'HC'                 'TLE' };

%% Load Data
tic
for iG = 1:numel(grp_nme)
    
    sbj_dir = dir([ dta_dir '/' grp_dir{iG}]); sbj_dir = {sbj_dir(3:end).name};
    
    dum_dta = load([ dta_dir '/' grp_dir{iG} '/' sbj_dir{1}]);
    ses_nme = fieldnames(dum_dta);
    dum_dta_sze = size(dum_dta.(ses_nme{1}).vbm_gm.dat);
    raw_dta.(grp_nme{iG}) = nan(dum_dta_sze(1),dum_dta_sze(2),dum_dta_sze(3),numel(sbj_dir));
    raw_dta.([grp_nme{iG} '_' 'nme']) = cell(numel(sbj_dir),1);
    
    for iS = 1:numel(sbj_dir)
        dta_hld = load([ dta_dir '/' grp_dir{iG} '/' sbj_dir{iS}]);
        ses_nme = fieldnames(dta_hld);
        raw_dta.(grp_nme{iG})(:,:,:,iS) = dta_hld.(ses_nme{1}).vbm_gm.dat;
        raw_dta.([grp_nme{iG} '_' 'nme']){iS,1} = sbj_dir{iS}(1:end-4);
        if rem(iS,10)==0 
            sprintf('Group %i: %i out of %i',iG,iS,numel(sbj_dir))
        end
    end   
end

save([ dta_dir '/' 'data.mat'],'raw_dta','-v7.3');

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'covariates_fixed.csv'];
[ cov_out, cov_sbj, cov_col ] = ejk_dta_frm(fcfg);

cov_sbj_mtc = strrep(cov_sbj,'_','');

img_nme = [ raw_dta.HC_nme(:,1) ; raw_dta.TLE_nme(:,1) ];
img_dta = cat(4,raw_dta.HC,raw_dta.TLE);
[ raw_dta_mbm, cov_dta_ind ] = ismember(img_nme,cov_sbj_mtc);
[ img_nme(logical(raw_dta_mbm)) cov_sbj_mtc(cov_dta_ind(cov_dta_ind~=0))]

img_nme = img_nme(logical(raw_dta_mbm));
img_dta = img_dta(:,:,:,logical(raw_dta_mbm));

cov_sbj_mtc = cov_sbj_mtc(cov_dta_ind(cov_dta_ind~=0));
cov_out     = cov_out(cov_dta_ind(cov_dta_ind~=0),:);

tme_out = toc;
sprintf('Data Loaded: %f seconds',tme_out)

%% Find Reihaneh HC
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'predictions.csv'];
[ rei_dta, rei_sbj, rei_col ] = ejk_dta_frm(fcfg);

rei_sbj_mtc = strrep(rei_sbj,'_','');

rei_inc = ismember(cov_sbj_mtc,rei_sbj_mtc);
rei_inc = num2cell(rei_inc);
yes_ind = cellfun(@(x) x==1,rei_inc);
no_ind = cellfun(@(x) x==0,rei_inc);
rei_inc(yes_ind) = {'yes'};
rei_inc(no_ind) = {'no'};

cov_col = [cov_col {'Reihaneh_data'}];
cov_out = [cov_out rei_inc];

%% Investigate Groups
cls_col = find(strcmpi(cov_col,'class'));
loc_col = find(strcmpi(cov_col,'Location of Epilepsy'));
sde_col = find(strcmpi(cov_col,'Side of Epilepsy'));
rei_col = find(strcmpi(cov_col,'Reihaneh_data'));

% unique(cov_out(:,sde_col))

% EPD vs HC ALL
grp_nme = 'epd_con_all';
grp.(grp_nme).epd = find(  strcmpi(cov_out(:,cls_col),'EPD')  );
grp.(grp_nme).con = find(  strcmpi(cov_out(:,cls_col),'HC')  );

% TLE vs HC ALL
grp_nme = 'tle_con_all';
grp.(grp_nme).tle = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') );
grp.(grp_nme).con = find(  strcmpi(cov_out(:,cls_col),'HC')  );

% L-TLE vs HC 
grp_nme = 'lft_tle_con_all';
grp.(grp_nme).lft_tle = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'L') );
grp.(grp_nme).con = find(  strcmpi(cov_out(:,cls_col),'HC')  );

% R-TLE vs HC 
grp_nme = 'rgh_tle_con_all';
grp.(grp_nme).rgh_tle = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'R') );
grp.(grp_nme).con = find(  strcmpi(cov_out(:,cls_col),'HC')  );

% R-TLE vs L-TLE 
grp_nme = 'rgh_tle_lft_tle_all';
grp.(grp_nme).rgh_tle = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'R') );
grp.(grp_nme).lft_tle = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'L') );

% L-TLE vs HC 
grp_nme = 'lft_tle_con_rei';
grp.(grp_nme).lft_tle_rei = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'L')  & strcmpi(cov_out(:,rei_col),'yes') );
grp.(grp_nme).con_rei     = find(  strcmpi(cov_out(:,cls_col),'HC')  & strcmpi(cov_out(:,rei_col),'yes') );

% L-TLE vs nonHCP-HC 
grp_nme = 'lft_tle_non_humc';
grp.(grp_nme).lft_tle     = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'L') );
grp.(grp_nme).con_non_hmc = setxor(find(  strcmpi(cov_out(:,cls_col),'HC') ) , string_find(cov_sbj_mtc,'HUM') );

% L-TLE vs nonHCP-HC 
grp_nme = 'lft_tle_yes_humc';
grp.(grp_nme).lft_tle     = find(  strcmpi(cov_out(:,cls_col),'EPD') & strcmpi(cov_out(:,loc_col),'TL') & strcmpi(cov_out(:,sde_col),'L') );
grp.(grp_nme).con_yes_hmc = intersect(find(  strcmpi(cov_out(:,cls_col),'HC') ) , string_find(cov_sbj_mtc,'HUM') );

%% Run VBM %%%%%%%%%%%%%%%

grp_nme = fieldnames(grp);

for iG = 1:numel(grp_nme)

    grp_vbm = fieldnames(grp.(grp_nme{iG}));

    %% VBM
    fcfg = [];
    fcfg.dta_one = img_dta(:,:,:,grp.(grp_nme{iG}).(grp_vbm{1}));
    fcfg.dta_two = img_dta(:,:,:,grp.(grp_nme{iG}).(grp_vbm{2}));
    vbm_out = ejk_vbm(fcfg);

    save([ plt_dir '/' 'vbm' '/' 'vbm_' grp_nme{iG} '_' grp_vbm{1} '_VS_' grp_vbm{2} '.mat'],'vbm_out')

    %% Montage
    avg_brn = nanmean(raw_dta.HC,4);
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

    fcfg.plt_dir = [ plt_dir '/' 'vbm' '/'];
    fcfg.plt_nme = [ 'montage' '_' grp_nme{iG} '_' grp_vbm{1} '_VS_' grp_vbm{2} ];

    ejk_montage(fcfg)

end




