clear; clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/04_Developement/ict_cog/data/psych_info';
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/04_Developement/ict_cog/data/ejk';
dem_fle = 'Ictal2_Demographics.csv';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_ltr = 'S';
% par_ltr = 'E';
    shm_enc = 'sham_encoding';
    shm_imd = 'sham_immediate_recall';
    shm_dly = 'sham_delayed_recall';
    sze_enc = 'seizure_encoding';
    sze_imd = 'seizure_immediate_recall';
    sze_dly = 'seizure_delayed_recall';
    sod_ata = 'sodata';
cmb_fld = 'behav_seizure_delayed_recall';
    cmb_fle = 'seizure_delayed_recall';
    cmb_trl = '';
    
%% Define extraction
enc_fld_nme = { 'trial_number'          'correct_response' 'reaction_time'                                          'seizure_str' 'word' 'word_type' 'word_length' 'frequency' 'concrete' };
rcl_fld_nme = { 'trial_number' 'target' 'correct_response' 'reaction_time' 'confidence' 'confidence_reaction_time'  'seizure_str' 'word' 'word_type' 'word_length' 'frequency' 'concrete' };
    
%% Load Data
sbj_dir = dir(dta_dir); sbj_dir = {sbj_dir([sbj_dir(:).isdir]).name}; sbj_dir = sbj_dir(cellfun(@(x) strcmpi(x(1),sbj_ltr),sbj_dir));

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dem_fle];
[dem_dta, dem_sbj, dem_col] = ejk_dta_frm(fcfg);

%% File Exploration
% Subject 1
iS = 1;
tot_cat = load([ dta_dir '/' sbj_dir{iS} '/' cmb_fld '/' 'S21_145_seizure_delayed_recall.mat']); % Concatenated total events
tot_srt = load([ dta_dir '/' sbj_dir{iS} '/' cmb_fld '/' 'trialinfo_behav_seizure_delayed_recall.mat']); % Sorted total events

beh_fld = dir([dta_dir '/'  sbj_dir{iS} '/']); beh_fld = {beh_fld([beh_fld(:).isdir]).name}; beh_fld = beh_fld(cellfun(@(x) strcmpi(x(1),par_ltr),beh_fld));
    
mat_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{2} '/' 'S21_145_sham_immediate_recall.mat']);
trl_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{2} '/' 'trialinfo_E21-273_0023.mat']);

% Subject 22
iS = 22;
beh_fld = dir([dta_dir '/'  sbj_dir{iS} '/']); beh_fld = {beh_fld([beh_fld(:).isdir]).name}; beh_fld = beh_fld(cellfun(@(x) strcmpi(x(1),par_ltr),beh_fld));
    
mat_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{2} '/' 'sodata.S23_199.01.04.2023.10.38.mat']);
trl_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{2} '/' 'trialinfo_E23-457_0006.mat']);

% Subject 24
iS = 24;
beh_fld = dir([dta_dir '/'  sbj_dir{iS} '/']); beh_fld = {beh_fld([beh_fld(:).isdir]).name}; beh_fld = beh_fld(cellfun(@(x) strcmpi(x(1),par_ltr),beh_fld));
    
mat_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{1} '/' 'sodata.S23_202.02.05.2023.16.10.mat']);
trl_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{1} '/' 'trialinfo_E23-610_0003.mat']);

% Seizure trialinfo
ttt = load('/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/04_Developement/ict_cog/data/psych_info/ict_Seizure_trialinfos_complete/S21_145_PCb_E21-273_0056.mat');

ttt = load('/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/04_Developement/ict_cog/data/psych_info/ict_Seizure_trialinfos_complete/S22_178_AF_E22-177_0139.mat');

ttt = load('/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/04_Developement/ict_cog/data/psych_info/ict_Seizure_trialinfos_complete/S21_165_WN_E21-465_0093.mat');


%% File Presence
out_fle = cell(numel(sbj_dir),8);

for iS = 1:numel(sbj_dir)
    
    shm_enc_num = 0;
    shm_imd_num = 0;
    shm_dly_num = 0;
    sze_enc_num = 0;
    sze_imd_num = 0;
    sze_dly_num = 0;
    sod_ata_num = 0;
    
    beh_fld = dir([dta_dir '/'  sbj_dir{iS} '/']); beh_fld = {beh_fld([beh_fld(:).isdir]).name};
    for iB = 1:numel(beh_fld)
        
        prs_fle = dir([dta_dir '/'  sbj_dir{iS} '/' beh_fld{iB} '/']); prs_fle = {prs_fle(~[prs_fle(:).isdir]).name};
        
        if     any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_enc),prs_fle,'uni',0))); shm_enc_num = shm_enc_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_imd),prs_fle,'uni',0))); shm_imd_num = shm_imd_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_dly),prs_fle,'uni',0))); shm_dly_num = shm_dly_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_enc),prs_fle,'uni',0))); sze_enc_num = sze_enc_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_imd),prs_fle,'uni',0))); sze_imd_num = sze_imd_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_dly),prs_fle,'uni',0))); sze_dly_num = sze_dly_num + 1;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sod_ata),prs_fle,'uni',0))); sod_ata_num = sod_ata_num + 1;
        end        
    end    
    
    out_fle(iS,:) = [ sbj_dir{iS} {shm_enc_num} {shm_imd_num} {shm_dly_num} {sze_enc_num} {sze_imd_num} {sze_dly_num} {sod_ata_num} ];
    
end

cell2csv([ out_dir '/' 'file_check.csv'],out_fle)

%% Load Subject Data
clear out_dta
for iS = 1:numel(sbj_dir)
    
    beh_fld = dir([dta_dir '/'  sbj_dir{iS} '/']); beh_fld = {beh_fld([beh_fld(:).isdir]).name};
    
    for iB = 1:numel(beh_fld)
        
        prs_fle = dir([dta_dir '/'  sbj_dir{iS} '/' beh_fld{iB} '/']); prs_fle = {prs_fle(~[prs_fle(:).isdir]).name};
        
        % Find folder type
        if     any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_enc),prs_fle,'uni',0)))
            typ_nme = shm_enc; dta_pul = enc_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_imd),prs_fle,'uni',0)))
            typ_nme = shm_imd; dta_pul = rcl_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,shm_dly),prs_fle,'uni',0)))
            typ_nme = shm_dly; dta_pul = rcl_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_enc),prs_fle,'uni',0)))
            typ_nme = sze_enc; dta_pul = enc_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_imd),prs_fle,'uni',0)))
            typ_nme = sze_imd; dta_pul = rcl_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sze_dly),prs_fle,'uni',0)))
            typ_nme = sze_dly; dta_pul = rcl_fld_nme;
        elseif any(~cellfun(@isempty,cellfun(@(x) strfind(x,sod_ata),prs_fle,'uni',0)))
            typ_nme = sod_ata; dta_pul = [];
        else
            typ_nme = []; dta_pul = [];
        end
        
        % Extract info
        if ~isempty(dta_pul)
        
            % Load file
            trl_ind = string_find(prs_fle,'trialinfo');
            trl_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{iB} '/' prs_fle{trl_ind} ]);
            mat_ind_one = string_find(prs_fle,'S2'); mat_ind_two = string_find(prs_fle,typ_nme); mat_ind = intersect(mat_ind_one,mat_ind_two);
            mat_fle = load([ dta_dir '/' sbj_dir{iS} '/' beh_fld{iB} '/' prs_fle{mat_ind} ]);
           
            trl_fld = fieldnames(trl_fle.trialinfo);
            [ ~, trl_fld_lap, dta_pul_lap ] = intersect(trl_fld,dta_pul);
            
            out_dta.(sbj_dir{iS}).(typ_nme).fieldnames = setxor(trl_fld,dta_pul);
            for iOL = 1:numel(trl_fld_lap)
                out_dta.(sbj_dir{iS}).(typ_nme).(trl_fld{trl_fld_lap(iOL)}) = table2array(trl_fle.trialinfo(:,trl_fld_lap(iOL)));
                if isnumeric(out_dta.(sbj_dir{iS}).(typ_nme).(trl_fld{trl_fld_lap(iOL)}))
                out_dta.(sbj_dir{iS}).(typ_nme).(trl_fld{trl_fld_lap(iOL)})(isnan(out_dta.(sbj_dir{iS}).(typ_nme).(trl_fld{trl_fld_lap(iOL)}))) = 0;
                end
            end            
        end
    end
end

%% Code subject performance
clear clc_dta clc_dta

typ_cdn = { 'sham_immediate_recall' 'sham_delayed_recall' };
clc_dta = NaN(numel(sbj_dir),0); clc_col = cell(0);

%
trg_col = 'target';
cor_col = 'correct_response';
rea_col = 'reaction_time';

% Hit Rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col_nme = 'hit_rate';
use_ind = size(clc_dta,2);
clc_dta = [ clc_dta NaN(numel(sbj_dir),2)];
for iC = 1:numel(typ_cdn)
    clc_col = [ clc_col [col_nme '_' typ_cdn{iC}] ];
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            dta_hld = out_dta.(sbj_dir{iS}).(typ_cdn{iC});
            clc_dta(iS,use_ind+iC) = (sum(dta_hld.(trg_col) & dta_hld.(cor_col)) / sum(dta_hld.(trg_col))) * 100;
        else
            clc_dta(iS,use_ind+iC) = NaN;
        end
    end
end

% False-Alarm Rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col_nme = 'false_alarm';
use_ind = size(clc_dta,2);
clc_dta = [ clc_dta NaN(numel(sbj_dir),2)];
for iC = 1:numel(typ_cdn)
    clc_col = [ clc_col [col_nme '_' typ_cdn{iC}] ];
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            dta_hld = out_dta.(sbj_dir{iS}).(typ_cdn{iC});
            clc_dta(iS,use_ind+iC) = (sum(~dta_hld.(trg_col) & ~dta_hld.(cor_col)) / sum(~dta_hld.(trg_col))) * 100;
        else
            clc_dta(iS,use_ind+iC) = NaN;
        end
    end
end

% d'prime rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col_nme = 'd_prime';
use_ind = size(clc_dta,2);
clc_dta = [ clc_dta NaN(numel(sbj_dir),2)];
for iC = 1:numel(typ_cdn)
    clc_col = [ clc_col [col_nme '_' typ_cdn{iC}] ];
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            hit_col = strcmpi(clc_col,['hit_rate' '_' typ_cdn{iC}]);
            fls_col = strcmpi(clc_col,['false_alarm' '_' typ_cdn{iC}]);
            clc_dta(iS,use_ind+iC) = norminv(clc_dta(iS,hit_col)/100) - norminv(clc_dta(iS,fls_col)/100);
        else
            clc_dta(iS,use_ind+iC) = NaN;
        end
    end
end

% correct response RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col_nme = 'hit_rate_RT';
use_ind = size(clc_dta,2);
clc_dta = [ clc_dta NaN(numel(sbj_dir),2)];
for iC = 1:numel(typ_cdn)
    clc_col = [ clc_col [col_nme '_' typ_cdn{iC}] ];
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            dta_hld = out_dta.(sbj_dir{iS}).(typ_cdn{iC});
            clc_dta(iS,use_ind+iC) = (sum(dta_hld.(rea_col)(dta_hld.(cor_col) & dta_hld.(trg_col)))*1000) / sum(dta_hld.(cor_col) & dta_hld.(trg_col));
        else
            clc_dta(iS,use_ind+iC) = NaN;
        end
    end
end

% correct response RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col_nme = 'hit_rate_RT_std';
use_ind = size(clc_dta,2);
clc_dta = [ clc_dta NaN(numel(sbj_dir),2)];
for iC = 1:numel(typ_cdn)
    clc_col = [ clc_col [col_nme '_' typ_cdn{iC}] ];
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            dta_hld = out_dta.(sbj_dir{iS}).(typ_cdn{iC});
            clc_dta(iS,use_ind+iC) = std((dta_hld.(rea_col)(dta_hld.(cor_col) & dta_hld.(trg_col)))*1000);
        else
            clc_dta(iS,use_ind+iC) = NaN;
        end
    end
end

% correct response RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell2csv([ out_dir '/' 'performance.csv' ], [ {'sbj_nme'}  clc_col;  sbj_dir' num2cell(clc_dta)] )

%% Hack together ROCs
typ_cdn = { 'sham_immediate_recall' 'sham_delayed_recall' };

trg_col = 'target';
cor_col = 'correct_response';
cnf_col = 'confidence';

% Recode into bins
for iC = 1:numel(typ_cdn)
    bin_cnt.(typ_cdn{iC}) = zeros(numel(sbj_dir),8);
    for iS = 1:numel(sbj_dir)
        if any(ismember(fieldnames(out_dta),sbj_dir{iS})) && any(ismember(fieldnames(out_dta.(sbj_dir{iS})),typ_cdn{iC}))
            dta_hld = out_dta.(sbj_dir{iS}).(typ_cdn{iC});
            for iT = 1:numel(dta_hld.(trg_col))
                if     dta_hld.(trg_col)(iT)==1 && dta_hld.(cor_col)(iT)==0 && dta_hld.(cnf_col)(iT)==2 % Target, New, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,1) = bin_cnt.(typ_cdn{iC})(iS,1) + 1;
                elseif dta_hld.(trg_col)(iT)==1 && dta_hld.(cor_col)(iT)==0 && dta_hld.(cnf_col)(iT)==1 % Target, New, Low Confidence
                    bin_cnt.(typ_cdn{iC})(iS,2) = bin_cnt.(typ_cdn{iC})(iS,2) + 1;
                elseif dta_hld.(trg_col)(iT)==1 && dta_hld.(cor_col)(iT)==1 && dta_hld.(cnf_col)(iT)==1 % Target, Old, Low Confidence
                    bin_cnt.(typ_cdn{iC})(iS,3) = bin_cnt.(typ_cdn{iC})(iS,3) + 1;
                elseif dta_hld.(trg_col)(iT)==1 && dta_hld.(cor_col)(iT)==1 && dta_hld.(cnf_col)(iT)==2 % Target, Old, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,4) = bin_cnt.(typ_cdn{iC})(iS,4) + 1;
                elseif dta_hld.(trg_col)(iT)==0 && dta_hld.(cor_col)(iT)==1 && dta_hld.(cnf_col)(iT)==2 % Foil, New, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,5) = bin_cnt.(typ_cdn{iC})(iS,5) + 1;
                elseif dta_hld.(trg_col)(iT)==0 && dta_hld.(cor_col)(iT)==1 && dta_hld.(cnf_col)(iT)==1 % Foil, New, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,6) = bin_cnt.(typ_cdn{iC})(iS,6) + 1;
                elseif dta_hld.(trg_col)(iT)==0 && dta_hld.(cor_col)(iT)==0 && dta_hld.(cnf_col)(iT)==1 % Foil, Old, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,7) = bin_cnt.(typ_cdn{iC})(iS,7) + 1;
                elseif dta_hld.(trg_col)(iT)==0 && dta_hld.(cor_col)(iT)==0 && dta_hld.(cnf_col)(iT)==2 % Foil, Old, High Confidence
                    bin_cnt.(typ_cdn{iC})(iS,8) = bin_cnt.(typ_cdn{iC})(iS,8) + 1;
                end
            end
        else
            bin_cnt.(typ_cdn{iC})(iS,:) = NaN;
        end
    end
end

% Calculate 
bin_num = size(bin_cnt.(typ_cdn{iC}),2)/2;
bin_seq = bin_num:-1:2;
for iC = 1:numel(typ_cdn)
    bin_pct.(typ_cdn{iC}) = zeros(numel(sbj_dir),bin_num-1);
    for iS = 1:numel(sbj_dir)
        dta_hld = bin_cnt.(typ_cdn{iC})(iS,:);
        tot_trg = sum(dta_hld(1:bin_num));
        tot_foi = sum(dta_hld(1+bin_num:bin_num*2));
        for iCO = 1:numel(bin_seq)
            bin_pct.(typ_cdn{iC})(iS,iCO)           = (sum(dta_hld(bin_seq(iCO):bin_num))/tot_trg)*100;
            bin_pct.(typ_cdn{iC})(iS,iCO+bin_num-1) = (sum(dta_hld(bin_seq(iCO)+bin_num:bin_num*2))/tot_foi)*100;
        end
    end
end

%% Sham analyses (individual)
sbj_nme = fieldnames(out_dta);
sbj_col = distinguishable_colors(numel(sbj_nme));

% Overall Performance %%%%%%%%%%%%%%%%%%%
plt_fld_nme = {'sham_encoding' ...
               'sham_immediate_recall' ...
               'sham_delayed_recall' };
prf_fld_nme = { {'correct_response' 'reaction_time'} ...
                {'correct_response' 'reaction_time'} ...
                {'correct_response' 'reaction_time'} };
prf_fld_mod = { [100                1000] ...
                [100                1000] ...
                [100                1000] };

for iP = 1:numel(plt_fld_nme)
    
    figure('Visible','off');
    
    for iSB = 1:numel(prf_fld_nme{iP})
        
        sbp_hld = subplot(1,numel(prf_fld_nme{iP}),iSB);
        
        fcfg = [];
        
        sbj_ind = 1;
        for iS = 1:numel(sbj_nme)            
            try
                fcfg.ydt{sbj_ind}     = nansum(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}).(prf_fld_nme{iP}{iSB}))) / numel(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}).(prf_fld_nme{iP}{iSB}))) * prf_fld_mod{iP}(iSB);
                fcfg.xdt{sbj_ind}     = 1;
                fcfg.fce_col{sbj_ind} = sbj_col(iS,:);
                sbj_ind = sbj_ind + 1;
            catch; end
        end
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
        
        fcfg.xlb = {'Patients'};
        fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
        fcfg.ylb = prf_fld_nme{iP}(iSB);
        
        fcfg.jtr = 1;
        fcfg.jtr_wdt = 0.075;
            
        fcfg.sbp = sbp_hld;
        
        ejk_scatter(fcfg)
        
    end
    
    print( [ out_dir '/' 'explore' '/' 'raw_performance' '_' plt_fld_nme{iP} '.png'], '-dpng')
    close all;
end

%% Sham analyses (across)
sbj_nme = fieldnames(out_dta);
sbj_col = distinguishable_colors(numel(sbj_nme));

% Overall Performance %%%%%%%%%%%%%%%%%%%
plt_fld_nme = {'sham_encoding' ...
               'sham_immediate_recall' ...
               'sham_delayed_recall' };
prf_fld_nme = { {'correct_response' 'reaction_time'} ...
                {'correct_response' 'reaction_time'} ...
                {'correct_response' 'reaction_time'} };
prf_fld_mod = { [100                1000] ...
                [100                1000] ...
                [100                1000] };

figure('Visible','off');

for iSB = 1:numel(prf_fld_nme{iP})
    
    sbp_hld = subplot(1,numel(prf_fld_nme{iP}),iSB);
    
    fcfg = [];
    
    for iS = 1:numel(sbj_nme)
        tot_fld_nme = fieldnames(out_dta.(sbj_nme{iS}));
        for iP = 1:numel(plt_fld_nme)
            if any(ismember(tot_fld_nme,plt_fld_nme{iP}))
                fcfg.ydt{iS}(iP) = nansum(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}).(prf_fld_nme{iP}{iSB}))) / numel(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}).(prf_fld_nme{iP}{iSB}))) * prf_fld_mod{iP}(iSB);
            else
                fcfg.ydt{iS}(iP) = NaN;    
            end
        end
        fcfg.xdt{iS} = 1:numel(plt_fld_nme);
        fcfg.fce_col{iS} = sbj_col(iS,:);
    end
    
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    
    fcfg.xlb = mmil_spec_char(plt_fld_nme,{'_'},{' '});
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = prf_fld_nme{iP}(iSB);
    
    fcfg.jtr = 1;
    fcfg.jtr_wdt = 0.125;
    
    fcfg.sbp = sbp_hld;
    
    ejk_scatter(fcfg)
    
end

print( [ out_dir '/' 'explore' '/' 'across_sham_performance.png'], '-dpng')
close all;

%% Subtraction Analyses (individual)
sbj_nme = fieldnames(out_dta);
sbj_col = distinguishable_colors(numel(sbj_nme));

% Overall Performance %%%%%%%%%%%%%%%%%%%
plt_out_nme = { 'sham: enc - imd' ...
                'sham: imd - del' };
plt_fld_nme = { { 'sham_encoding' 'sham_immediate_recall'} ...
                { 'sham_immediate_recall' 'sham_delayed_recall'} };
prf_fld_nme = { {'correct_response' 'reaction_time'} ...
                {'correct_response' 'reaction_time'} };
prf_fld_mod = { [100                1000] ...
                [100                1000] };

figure('Visible','off');

for iSB = 1:numel(prf_fld_nme{iP})
    
    sbp_hld = subplot(1,numel(prf_fld_nme{iP}),iSB);
    
    fcfg = [];
    
    for iS = 1:numel(sbj_nme)
        tot_fld_nme = fieldnames(out_dta.(sbj_nme{iS}));
        for iP = 1:numel(plt_fld_nme)
            if any(ismember(tot_fld_nme,plt_fld_nme{iP}))
                fcfg.ydt{iS}(iP) = ((nansum(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}{1}).(prf_fld_nme{iP}{iSB}))) / numel(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}{1}).(prf_fld_nme{iP}{iSB})))) - ...
                                   (nansum(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}{2}).(prf_fld_nme{iP}{iSB}))) / numel(table2array(out_dta.(sbj_nme{iS}).(plt_fld_nme{iP}{2}).(prf_fld_nme{iP}{iSB}))))) ...
                                    * prf_fld_mod{iP}(iSB);
            else
                fcfg.ydt{iS}(iP) = NaN;    
            end
        end
        fcfg.xdt{iS} = 1:numel(plt_fld_nme);
        fcfg.fce_col{iS} = sbj_col(iS,:);
    end
    
    fcfg.hln = 0;
    
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    
    fcfg.xlb = mmil_spec_char(plt_out_nme,{'_'},{' '});
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = prf_fld_nme{iP}(iSB);
    
    fcfg.jtr = 1;
    fcfg.jtr_wdt = 0.125;
    
    fcfg.sbp = sbp_hld;
    
    ejk_scatter(fcfg)
    
end

print( [ out_dir '/' 'explore' '/' 'across_sham_subtraction.png'], '-dpng')
close all;

%% Accelerated forgetting plots
sbj_col = distinguishable_colors(numel(sbj_dir));
plt_out_dir = [ out_dir '/' 'measures_plot' ];
typ_cdn = { 'sham_immediate_recall' 'sham_delayed_recall' };
mes_plt = { 'hit_rate' 'false_alarm' 'd_prime' 'hit_rate_RT' 'hit_rate_RT_std' };

for iM = 1:numel(mes_plt)
    
    figure('Visible','off');
    
    num_str = num2str(iM);
    if length(num_str)==1; num_str = ['0' num_str]; end
    
    % Measures plot
    sbp_hld = subplot(1,2,1);
    
    fcfg = [];
    
    for iS = 1:numel(sbj_dir)
        for iC = 1:numel(typ_cdn)
            fcfg.ydt{iS}(iC) = clc_dta(iS,strcmpi(clc_col,[mes_plt{iM} '_' typ_cdn{iC}]));
        end
        fcfg.xdt{iS} = 1:numel(typ_cdn);
        fcfg.fce_col{iS} = sbj_col(iS,:);
    end
    
    fcfg.hln = 0;
    
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    
    fcfg.xlb = mmil_spec_char(typ_cdn,{'_'},{' '});
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = mes_plt(iM);
    
    fcfg.jtr = 1;
    fcfg.jtr_wdt = 0.125;
    
    fcfg.sbp = sbp_hld;
    
    ejk_scatter(fcfg)
    
    % Subtraction plot
    sbp_hld = subplot(1,2,2);
    
    [~,pvl] = ttest(clc_dta(:,strcmpi(clc_col,[mes_plt{iM} '_' typ_cdn{1}])), clc_dta(:,strcmpi(clc_col,[mes_plt{iM} '_' typ_cdn{2}])));
    
    fcfg = [];
    
    for iS = 1:numel(sbj_dir)
        fcfg.ydt{iS} = clc_dta(iS,strcmpi(clc_col,[mes_plt{iM} '_' typ_cdn{1}])) - clc_dta(iS,strcmpi(clc_col,[mes_plt{iM} '_' typ_cdn{2}]));
        fcfg.xdt{iS} = 1;
        fcfg.fce_col{iS} = sbj_col(iS,:);
    end
    
    fcfg.hln = 0;
    
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    
    fcfg.xlb = {['Subtraction' ', ' 'p=' num2str(roundsd(pvl,3))]};
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = mes_plt(iM);
    
    fcfg.jtr = 1;
    fcfg.jtr_wdt = 0.125;
    fcfg.hln = 0;
     
    fcfg.sbp = sbp_hld;
    
    ejk_scatter(fcfg)
    
    %
    print( [ plt_out_dir '/' num_str '_' mes_plt{iM} '.png'], '-dpng')
    close all;
    
end

%% Accelerated forgetting by combination
fcfg = [];
fcfg.dta_loc = [out_dir '/' 'performance_factors_and_subject_factors.csv'];
[ dta_plt, sbj_plt, col_plt] = ejk_dta_frm(fcfg);

plt_out_dir = [ out_dir '/' 'measures_plot_split' ];
typ_cdn = { 'sham_immediate_recall' 'sham_delayed_recall' };
typ_cdn_ind = {[1 1.5 2] [3 3.5 4]};
mes_plt = { 'hit_rate' 'false_alarm' 'd_prime' 'hit_rate_RT' 'hit_rate_RT_std' };

grp_plt     = { 'sde' 'sde_tmp' };
sub_grp_plt = { 'lft' 'bil' 'rgh' };
sub_grp_col = { rgb('blue') rgb('purple') rgb('red') };

% Groups
lft_ind = find(strcmpi(dta_plt(:,strcmpi(col_plt,'Seizure_Hemisphere')),'L'));
rgh_ind = find(strcmpi(dta_plt(:,strcmpi(col_plt,'Seizure_Hemisphere')),'R'));
bil_ind = find(strcmpi(dta_plt(:,strcmpi(col_plt,'Seizure_Hemisphere')),'B'));
tmp_ind = find(strcmpi(dta_plt(:,strcmpi(col_plt,'SOZ_simple')),'Temp'));

grp.sde.lft = lft_ind;
grp.sde.rgh = rgh_ind;
grp.sde.bil = bil_ind;

grp.sde_tmp.lft = intersect(lft_ind,tmp_ind);
grp.sde_tmp.rgh = intersect(rgh_ind,tmp_ind);
grp.sde_tmp.bil = intersect(bil_ind,tmp_ind);

% Plot
for iM = 1:numel(mes_plt)
    
    figure('Visible','off');
    
    num_str = num2str(iM);
    if length(num_str)==1; num_str = ['0' num_str]; end
        
    for iG = 1:numel(grp_plt)       
        
        % Measures plot
        sbp_hld = subplot(2,2,((iG-1)*2)+1);
        
        fcfg = [];
        
        ydt_ind = 1;
        for iC = 1:numel(typ_cdn)            
            dta_col = find(strcmpi(col_plt,[ mes_plt{iM} '_' typ_cdn{iC} ]));
            for iSG = 1:numel(sub_grp_plt)
                fcfg.ydt{ydt_ind} = [dta_plt{grp.(grp_plt{iG}).(sub_grp_plt{iSG}),dta_col}]';
                fcfg.xdt{ydt_ind} = typ_cdn_ind{iC}(iSG);
                fcfg.fce_col{ydt_ind} = sub_grp_col{iSG};
                fcfg.xlb{ydt_ind} = sub_grp_plt{iSG};
                fcfg.box_plt(ydt_ind) = 1;
                fcfg.box_plt_col{ydt_ind} = sub_grp_col{iSG};
                ydt_ind = ydt_ind + 1;
            end
        end
                        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
        
        fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
        fcfg.ylb = mes_plt(iM);
        
        fcfg.jtr = 1;
        fcfg.jtr_wdt = 0.125;
        
        fcfg.sbp = sbp_hld;
        
        ejk_scatter(fcfg)
        
        % Subtraction plot
        sbp_hld = subplot(2,2,((iG-1)*2)+2);
        
        fcfg = [];
                
        dta_col_one = find(strcmpi(col_plt,[ mes_plt{iM} '_' typ_cdn{1} ]));
        dta_col_two = find(strcmpi(col_plt,[ mes_plt{iM} '_' typ_cdn{2} ]));
        for iSG = 1:numel(sub_grp_plt)
            fcfg.ydt{iSG}     = [dta_plt{grp.(grp_plt{iG}).(sub_grp_plt{iSG}),dta_col_one}]' - [dta_plt{grp.(grp_plt{iG}).(sub_grp_plt{iSG}),dta_col_two}]';
            fcfg.xdt{iSG}     = typ_cdn_ind{1}(iSG);
            fcfg.fce_col{iSG} = sub_grp_col{iSG};
            fcfg.xlb{iSG}     = sub_grp_plt{iSG};
            fcfg.box_plt(iSG) = 1;
            fcfg.box_plt_col{iSG} = sub_grp_col{iSG};
        end
        
        fcfg.hln = 0;
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
        
        fcfg.xlb = grp_plt(iG);
        fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
        fcfg.ylb = {'Subtraction'};
        
        fcfg.jtr = 1;
        fcfg.jtr_wdt = 0.125;
        fcfg.hln = 0;
        
        fcfg.sbp = sbp_hld;
        
        ejk_scatter(fcfg)
        
    end
    
    %
    print( [ plt_out_dir '/' num_str '_' mes_plt{iM} '.png'], '-dpng')
    close all;
end

[~,pvl] = ttest2(fcfg.ydt{1},[fcfg.ydt{2} ; fcfg.ydt{3}])
[~,pvl] = ttest2(fcfg.ydt{1},[fcfg.ydt{2} ; fcfg.ydt{3}],'Vartype','unequal')
[~,pvl] = ttest2(fcfg.ydt{1},fcfg.ydt{3})
[~,pvl] = ttest2(fcfg.ydt{1},fcfg.ydt{3},'Vartype','unequal')

%% Accelerated forgetting ROCs
sbj_col = distinguishable_colors(numel(sbj_dir));
typ_cdn = { 'sham_immediate_recall' 'sham_delayed_recall' };
plt_out_dir = [ out_dir '/' 'roc' ];

bin_num = size(bin_pct.(typ_cdn{iC}),2)/2;

% Total
fig_one = figure('Visible','off');
for iC = 1:numel(typ_cdn)
    sbp = subplot(1,2,iC);
    hold on;
    for iS = 1:numel(sbj_dir)
        plot(bin_pct.(typ_cdn{iC})(iS,1+bin_num:bin_num*2),bin_pct.(typ_cdn{iC})(iS,1:bin_num),'Color',sbj_col(iS,:), ...
            'Marker','diamond','MarkerFaceColor',sbj_col(iS,:),'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'LineWidth',3)
        
    end
    xlim([-5 105]);
    ylim([-5 105]);
    ylabel('Hit Rate');
    xlabel('False-Alarm Rate');
    title(mmil_spec_char(typ_cdn{iC},{'_'},{' '}))
    
    cell2csv([plt_out_dir '/' typ_cdn{iC} '_pct.csv'],  [ sbj_dir' num2cell(bin_pct.(typ_cdn{iC}))])
    cell2csv([plt_out_dir '/' typ_cdn{iC} '_count.csv'],[ sbj_dir' num2cell(bin_cnt.(typ_cdn{iC}))])
end
tightfig();
print(fig_one,[ plt_out_dir '/' 'total.png'],'-dpng')
close(fig_one)

% Individual Subject
for iS = 1:numel(sbj_dir)
    fig_one = figure('Visible','off');
    for iC = 1:numel(typ_cdn)
        sbp = subplot(1,2,iC);
        hold on;
        
        plot(bin_pct.(typ_cdn{iC})(iS,1+bin_num:bin_num*2),bin_pct.(typ_cdn{iC})(iS,1:bin_num),'Color',sbj_col(iS,:), ...
            'Marker','diamond','MarkerFaceColor',sbj_col(iS,:),'MarkerEdgeColor',[0 0 0],'MarkerSize',10,'LineWidth',3)
        
        xlim([-5 105]);
        ylim([-5 105]);
        ylabel('Hit Rate');
        xlabel('False-Alarm Rate');
        title(mmil_spec_char(typ_cdn{iC},{'_'},{' '}))
        
    end
    
    tightfig();
    print(fig_one,[ plt_out_dir '/' sbj_dir{iS} '.png'],'-dpng')
    close(fig_one)
    
end





