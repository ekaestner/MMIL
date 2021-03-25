clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
red_fle = 'mmil_redcap_2020_05_25.csv';

%%
% Redcap Data %%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[ sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo] = mmil_load_redcap(fcfg);

sbj_ind_bdi = find(~isnan(sbj_emo.bdi));
sbj_ind_bai = find(~isnan(sbj_emo.bai));
sbj_ind_bth = intersect( sbj_ind_bdi, sbj_ind_bai);
sbj_ind_all = unique( [ sbj_ind_bdi; sbj_ind_bai ]);

gry_thk = mmil_readtext([ prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'MRI_thickness_aparc_xsplit.csv']);
    [ ~, use_ind ] = intersect( gry_thk(:,1), sbj_emo.sbj_nme(sbj_ind_all, 1));
    gry_thk = gry_thk( [ 1 ; use_ind ], :);
    
sub_vol = mmil_readtext([ prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Volumes_aparc.csv']);    
    [ ~, use_ind ] = intersect( sub_vol(:,1), sbj_emo.sbj_nme(sbj_ind_all, 1));
    sub_vol = sub_vol( [ 1 ; use_ind ], :);

fib_tfa = mmil_readtext([ prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Fibers_aparc_FA.csv']);
    [ ~, use_ind ] = intersect( fib_tfa(:,1), sbj_emo.sbj_nme(sbj_ind_all, 1));
    fib_tfa = fib_tfa( [ 1 ; use_ind ], :);

fib_tmd = mmil_readtext([ prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Fibers_aparc_MD.csv']);
    [ ~, use_ind ] = intersect( fib_tmd(:,1), sbj_emo.sbj_nme(sbj_ind_all, 1));
    fib_tmd = fib_tmd( [ 1 ; use_ind ], :);

fib_wmd = mmil_readtext([ prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'WMParc_aparc__MD.csv']);
    [ ~, use_ind ] = intersect( fib_wmd(:,1), sbj_emo.sbj_nme(sbj_ind_all, 1));
    fib_wmd = fib_wmd( [ 1 ; use_ind ], :);
    
%% Characterize Spread
% Number
num_con_bai = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bai,1),{'fc'}) );
num_epd_bai = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bai,1),{'epd'}) );

num_con_bdi = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bdi,1),{'fc'}) );
num_epd_bdi = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bdi,1),{'epd'}) );

num_con_bth = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bth,1),{'fc'}) );
num_epd_bth = numel( string_find(sbj_emo.sbj_nme(sbj_ind_bth,1),{'epd'}) );

sbj_cat = sbj_emo.sbj_nme(:,1);
for iS = 1:size(sbj_cat, 1)

    % Categorize subjects
    if ~isempty( string_find( sbj_cat(iS, 1), {'epd'} ) )
        sbj_cat{iS, 2} = 'EPD';
    elseif ~isempty( string_find( sbj_cat(iS, 1), {'fc'} ) )
        sbj_cat{iS, 2} = 'HC';
    end
    
    % Categorize BDI
    if ~isnan(sbj_emo.bdi(iS)) && sbj_emo.bdi(iS)<=13
        sbj_cat{iS, 3} = 'Minimal';
    elseif ~isnan(sbj_emo.bdi(iS)) && sbj_emo.bdi(iS)>13 && sbj_emo.bdi(iS)<=19
        sbj_cat{iS, 3} = 'Mild';
    elseif ~isnan(sbj_emo.bdi(iS)) && sbj_emo.bdi(iS)>19 && sbj_emo.bdi(iS)<=28
        sbj_cat{iS, 3} = 'Moderate';
    elseif ~isnan(sbj_emo.bdi(iS)) && sbj_emo.bdi(iS)>28
        sbj_cat{iS, 3} = 'Severe';
    elseif isnan(sbj_emo.bdi(iS))
        sbj_cat{iS, 3} = 'N/A';
    end
    sbj_cat{iS, 4} = sbj_emo.bdi(iS);
    
    % Categorize BAI
    if ~isnan(sbj_emo.bai(iS)) && sbj_emo.bai(iS)<=9
        sbj_cat{iS, 5} = 'Minimal';
    elseif ~isnan(sbj_emo.bai(iS)) && sbj_emo.bai(iS)>9 && sbj_emo.bai(iS)<=16
        sbj_cat{iS, 5} = 'Mild';
    elseif ~isnan(sbj_emo.bai(iS)) && sbj_emo.bai(iS)>16 && sbj_emo.bai(iS)<=29
        sbj_cat{iS, 5} = 'Moderate';
    elseif ~isnan(sbj_emo.bai(iS)) && sbj_emo.bai(iS)>29
        sbj_cat{iS, 5} = 'Severe';
    elseif isnan(sbj_emo.bai(iS))
        sbj_cat{iS, 5} = 'N/A';
    end
    sbj_cat{iS, 6} = sbj_emo.bai(iS);
    
end

%% Correlation of BDI & BAI
ind_con = intersect( sbj_ind_bth, find(strcmpi(sbj_cat(:,2), 'HC')));
ind_epd = intersect( sbj_ind_bth, find(strcmpi(sbj_cat(:,2), 'EPD')));

subplot(1,2,1)
scatter( sbj_emo.bdi(ind_con) , sbj_emo.bai(ind_con) );
    ylabel('Beck''s Anxiety');xlabel('Beck''s Depression'); xlim([0 63]); ylim([0 63]);
    line([ 0 63], [0 63])
subplot(1,2,2)
scatter( sbj_emo.bdi(ind_epd) , sbj_emo.bai(ind_epd) );
    ylabel('Beck''s Anxiety');xlabel('Beck''s Depression'); xlim([0 63]); ylim([0 63]);
    line([ 0 63], [0 63])
   
%% T-TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Create IV
sbj_nme_bdi = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 1 );

sbj_grp_bdi(:,1) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 2 );
sbj_grp_bdi(:,2) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 3 );
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,1),'HC'), 2)       = {'N/A'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Mild'), 2)     = {'N/A'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Severe'), 2)   = {'Depressed'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Moderate'), 2) = {'Depressed'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Minimal'), 2)  = {'Not Depressed'};

%% Cognitive Categories   
% Create DV
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;
fcfg.dta_fld = { 'log_mem_nor_scr_one' 'log_mem_nor_scr_two' ...
                 'cvl_lfr_nor_scr' ...
                 'vp1_nor_scr' 'vp2_nor_scr' ...
                 'bnt_nor_scr' 'cat_flu_nor_scr' ...
                 'ltr_tot_raw_scr' 'swt_cor_nor_scr' 'swt_acc_nor_scr' };

cog_dta_frm = ejk_struct_2_dataframe( fcfg, sbj_cog );

 
% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = cog_dta_frm;
fcfg.dta_nme = { 'LM1' 'LM2' 'CVLT' 'VP1' 'VP2' 'BNT' 'CF' 'LTR' 'SWITCH_COR' 'SWITCH_ACC' };

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/cognitive';

ejk_ttest2_independent(fcfg)
 
%% Thickness Categories
% Create DV
gry_thk(:,find(isnan(cell2mat(gry_thk(2,2:end))))+1) = [];
gry_thk_use = nan( size(sbj_nme_bdi,1), size(gry_thk,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( gry_thk(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        gry_thk_use(iS,:) = cell2mat(gry_thk( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = gry_thk_use;
fcfg.dta_nme = gry_thk(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/corticalthickness';

ejk_ttest2_independent(fcfg) 
 
%% tract FA Categories - DEPRESSION
% Create DV
fib_tfa(:,find(isnan(cell2mat(fib_tfa(2,2:end))))+1) = [];
fib_tfa_use = nan( size(sbj_nme_bdi,1), size(fib_tfa,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( fib_tfa(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        fib_tfa_use(iS,:) = cell2mat(fib_tfa( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = fib_tfa_use;
fcfg.dta_nme = fib_tfa(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA';

ejk_ttest2_independent(fcfg) 

%% tract FA Categories - ANXIETY
% Create IV
sbj_nme_bai = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 6))), 1 );

sbj_grp_bai(:,1) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 6))), 2 );
sbj_grp_bai(:,2) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 6))), 5 );
    sbj_grp_bai( strcmpi(sbj_grp_bai(:,1),'HC'), 2)       = {'N/A'};
    sbj_grp_bai( strcmpi(sbj_grp_bai(:,2),'Mild'), 2)     = {'N/A'};
    sbj_grp_bai( strcmpi(sbj_grp_bai(:,2),'Severe'), 2)   = {'Depressed'};
    sbj_grp_bai( strcmpi(sbj_grp_bai(:,2),'Moderate'), 2) = {'Depressed'};
    sbj_grp_bai( strcmpi(sbj_grp_bai(:,2),'Minimal'), 2)  = {'Not Depressed'};

% Create DV
fib_tfa(:,find(isnan(cell2mat(fib_tfa(2,2:end))))+1) = [];
fib_tfa_use = nan( size(sbj_nme_bai,1), size(fib_tfa,2)-1 );
for iS = 1:size(sbj_nme_bai,1)
    use_ind = find(strcmpi( fib_tfa(:,1), sbj_nme_bai{iS,1}));
    if ~isempty(use_ind)
        fib_tfa_use(iS,:) = cell2mat(fib_tfa( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bai;

fcfg.dta     = fib_tfa_use;
fcfg.dta_nme = fib_tfa(1, 2:end);

fcfg.grp     = sbj_grp_bai;
fcfg.grp_nme = { 'Patient' 'PatientAnxiety' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA_anxiety';

ejk_ttest2_independent(fcfg) 

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Create IV
sbj_nme_bdi = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 1 );

sbj_grp_bdi(:,1) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 2 );
sbj_grp_bdi(:,2) = sbj_cat( ~isnan(cell2mat(sbj_cat(:, 4))), 3 );
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,1),'HC'), 2)       = {'HC'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Mild'), 2)     = {'N/A'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Severe'), 2)   = {'Depressed'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Moderate'), 2) = {'Depressed'};
    sbj_grp_bdi( strcmpi(sbj_grp_bdi(:,2),'Minimal'), 2)  = {'Not Depressed'};

%% Cognitive - Depression
% Create DV
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;
fcfg.dta_fld = { 'log_mem_nor_scr_one' 'log_mem_nor_scr_two' ...
                 'cvl_lfr_nor_scr' ...
                 'vp1_nor_scr' 'vp2_nor_scr' ...
                 'bnt_nor_scr' 'cat_flu_nor_scr' ...
                 'ltr_tot_raw_scr' 'swt_cor_nor_scr' 'swt_acc_nor_scr' };

cog_dta_frm = ejk_struct_2_dataframe( fcfg, sbj_cog );

 
% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = cog_dta_frm;
fcfg.dta_nme = { 'LM1' 'LM2' 'CVLT' 'VP1' 'VP2' 'BNT' 'CF' 'LTR' 'SWITCH_COR' 'SWITCH_ACC' };

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/cognitive/ANOVA';

ejk_1way_anova(fcfg)
    
%% Thickness - Depression
% Create DV
gry_thk(:,find(isnan(cell2mat(gry_thk(2,2:end))))+1) = [];
gry_thk_use = nan( size(sbj_nme_bdi,1), size(gry_thk,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( gry_thk(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        gry_thk_use(iS,:) = cell2mat(gry_thk( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = gry_thk_use;
fcfg.dta_nme = gry_thk(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/corticalthickness/ANOVA/';

ejk_1way_anova(fcfg) 

%% Subcortical Volume - Depression
% Create DV
sub_vol(:,find(isnan(cell2mat(sub_vol(2,2:end))))+1) = [];
sub_vol_use = nan( size(sbj_nme_bdi,1), size(sub_vol,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( sub_vol(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        sub_vol_use(iS,:) = cell2mat(sub_vol( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = sub_vol_use(:,1:32);
fcfg.dta_nme = strcat('x', cellfun( @(x) mmil_spec_char(x,{'-'}),sub_vol(1, 2:33),'UniformOutput',false));

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/Volume/ANOVA/';

ejk_1way_anova(fcfg) 

%% tFA - Depression
% Create DV
fib_tfa(:,find(isnan(cell2mat(fib_tfa(2,2:end))))+1) = [];
fib_tfa_use = nan( size(sbj_nme_bdi,1), size(fib_tfa,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( fib_tfa(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        fib_tfa_use(iS,:) = cell2mat(fib_tfa( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = fib_tfa_use;
fcfg.dta_nme = fib_tfa(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA/ANOVA/';

ejk_1way_anova(fcfg) 

%% tMD - Depression
% Create DV
fib_tmd(:,find(isnan(cell2mat(fib_tmd(2,2:end))))+1) = [];
fib_tmd_use = nan( size(sbj_nme_bdi,1), size(fib_tmd,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( fib_tmd(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        fib_tmd_use(iS,:) = cell2mat(fib_tmd( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = fib_tmd_use;
fcfg.dta_nme = fib_tmd(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractMD/ANOVA/';

ejk_1way_anova(fcfg) 

%% wMD - Depression
% Create DV
fib_wmd(:,find(isnan(cell2mat(fib_wmd(5,2:end))))+1) = [];
fib_wmd_use = nan( size(sbj_nme_bdi,1), size(fib_wmd,2)-1 );
for iS = 1:size(sbj_nme_bdi,1)
    use_ind = find(strcmpi( fib_wmd(:,1), sbj_nme_bdi{iS,1}));
    if ~isempty(use_ind)
        fib_wmd_use(iS,:) = cell2mat(fib_wmd( use_ind, 2:end));
    end
end

% Run ttest
fcfg = [];

fcfg.sbj_nme = sbj_nme_bdi;

fcfg.dta     = fib_wmd_use;
fcfg.dta_nme = fib_wmd(1, 2:end);

fcfg.grp     = sbj_grp_bdi;
fcfg.grp_nme = { 'Patient' 'PatientDepression' };

fcfg.out_dir  = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/wmparcMD/ANOVA/';

ejk_1way_anova(fcfg) 

%% DFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
out_dta.sbj_nme = sbj_nme_bdi;

out_dta.sbj_grp = sbj_grp_bdi(:,2);

out_dta.lhs_ist_cng = gry_thk_use(:,7); % L - Isthmus Cingulate
out_dta.rhs_ist_cng = gry_thk_use(:,56); % R - Isthmus Cingulate

out_dta.lhs_cng_pra = fib_tmd_use(:,5); % fib_tfa_use(:,5); % L - Cingulum (ParaHippocampal)
out_dta.rhs_cng_pra = fib_tmd_use(:,6); % fib_tfa_use(:,6); % R - Cingulum (ParaHippocampal)

%% epd_dep vs HC
chs_ind = logical(sum([ strcmpi(out_dta.sbj_grp,'HC') strcmpi(out_dta.sbj_grp,'Depressed') ], 2));

out_dta_frs.sbj_nme = out_dta.sbj_nme(chs_ind,1);

out_dta_frs.sbj_grp = out_dta.sbj_grp(chs_ind,1);

out_dta_frs.lhs_ist_cng = out_dta.lhs_ist_cng(chs_ind,1); % L - Isthmus Cingulate
out_dta_frs.rhs_ist_cng = out_dta.rhs_ist_cng(chs_ind,1); % R - Isthmus Cingulate
out_dta_frs.lhs_cng_pra = out_dta.lhs_cng_pra(chs_ind,1); % L - Cingulum (ParaHippocampal)
out_dta_frs.rhs_cng_pra = out_dta.rhs_cng_pra(chs_ind,1); % R - Cingulum (ParaHippocampal)

save( [ '/home/ekaestner/Downloads' '/' 'out_dta_frs.mat' ], 'out_dta_frs' )

%% epd_dep vs non
chs_ind = logical(sum([ strcmpi(out_dta.sbj_grp,'Not Depressed') strcmpi(out_dta.sbj_grp,'Depressed') ], 2));

out_dta_frs.sbj_nme = out_dta.sbj_nme(chs_ind,1);

out_dta_frs.sbj_grp = out_dta.sbj_grp(chs_ind,1);

out_dta_frs.lhs_ist_cng = out_dta.lhs_ist_cng(chs_ind,1); % L - Isthmus Cingulate
out_dta_frs.rhs_ist_cng = out_dta.rhs_ist_cng(chs_ind,1); % R - Isthmus Cingulate
out_dta_frs.lhs_cng_pra = out_dta.lhs_cng_pra(chs_ind,1); % L - Cingulum (ParaHippocampal)
out_dta_frs.rhs_cng_pra = out_dta.rhs_cng_pra(chs_ind,1); % R - Cingulum (ParaHippocampal)

save( [ '/home/ekaestner/Downloads' '/' 'out_dta_frs.mat' ], 'out_dta_frs' )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % Correlation of BAI with Memory
% subplot(1,2,1)
% scatter( sbj_cog.log_mem_nor_scr_one(ind_con) , sbj_emo.bai(ind_con) );
%     ylabel('Beck''s Anxiety');xlabel('Log Mem 1'); xlim([0 27]); ylim([0 63]);
% subplot(1,2,2)
% scatter( sbj_cog.log_mem_nor_scr_one(ind_epd) , sbj_emo.bai(ind_epd) );
%     ylabel('Beck''s Anxiety');xlabel('Log Mem 1'); xlim([0 27]); ylim([0 63]);
% 
% %
% subplot(1,2,1)
% scatter( sbj_cog.log_mem_nor_scr_two(ind_con) , sbj_emo.bai(ind_con) );
%     ylabel('Beck''s Anxiety');xlabel('Log Mem 2'); xlim([0 27]); ylim([0 63]);
% subplot(1,2,2)
% scatter( sbj_cog.log_mem_nor_scr_two(ind_epd) , sbj_emo.bai(ind_epd) );
%     ylabel('Beck''s Anxiety');xlabel('Log Mem 2'); xlim([0 27]); ylim([0 63]);
% 
% %
% subplot(1,2,1)
% scatter( sbj_cog.cvl_lfr_nor_scr(ind_con) , sbj_emo.bai(ind_con) );
%     ylabel('Beck''s Anxiety');xlabel('CVLT LFR'); xlim([0 8]); ylim([0 63]);
% subplot(1,2,2)
% scatter( sbj_cog.cvl_lfr_nor_scr(ind_epd) , sbj_emo.bai(ind_epd) );
%     ylabel('Beck''s Anxiety');xlabel('CVLT LFR'); xlim([0 8]); ylim([0 63]);
% 
% %
% subplot(1,2,1)
% scatter( sbj_cog.vp1_nor_scr_pst(ind_con) , sbj_emo.bai(ind_con) );
%     ylabel('Beck''s Anxiety');xlabel('VP1'); xlim([0 15]); ylim([0 63]);
% subplot(1,2,2)
% scatter( sbj_cog.vp1_nor_scr_pst(ind_epd) , sbj_emo.bai(ind_epd) );
%     ylabel('Beck''s Anxiety');xlabel('VP1'); xlim([0 15]); ylim([0 63]);
% 
% %
% subplot(1,2,1)
% scatter( sbj_cog.vp2_nor_scr_pst(ind_con) , sbj_emo.bai(ind_con) );
%     ylabel('Beck''s Anxiety');xlabel('VP2'); xlim([0 15]); ylim([0 63]);
% subplot(1,2,2)
% scatter( sbj_cog.vp2_nor_scr_pst(ind_epd) , sbj_emo.bai(ind_epd) );
%     ylabel('Beck''s Anxiety');xlabel('VP2'); xlim([0 15]); ylim([0 63]);








