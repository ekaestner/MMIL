%% Grab Data
load(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/manuscript/figure7/stimulation' '/' 'stm_rsl.mat'])

sel_ele = stm_rsl.sel_ele;
sel_txt = stm_rsl.ele_txt;

ecg_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical' '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total' '_' 'split' '_' 'lhs' '_' 'ecog'],'[, \t]');
ecg_loc{1} = ecg_hld(:,[1 5:end]);
for iR = 1:size(ecg_loc{1},1)
    if isempty(ecg_loc{1}{iR,2})
        ecg_loc{1}{iR,2} = 'unknown';
    else
        num_hld = cellfun(@(x) str2num(x(1:end-1)),ecg_loc{1}(iR,string_find(ecg_loc{1}(iR,~cellfun(@isempty,ecg_loc{1}(iR,:))),'%')),'uni',0);
        [~,max_hld] = max([num_hld{:}]);
        ecg_loc{1}{iR,2} = ecg_loc{1}{iR,max_hld+1+(1*max_hld)};
    end
end
ecg_loc{1} = ecg_loc{1}(:,1:2);

fcfg.cmb_reg = { { 'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    { 'lateraloccipital' } ...
    { 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
    { 'caudal-MTG' 'middle-MTG' 'rostral-MTG' } ...
    { 'caudal-STG' 'middle-STG' 'rostral-STG' } ...
    { 'inferior-precentral' 'middle-precentral' } ...
    { 'parstriangularis' } ...
    { 'parsopercularis' } };
fcfg.reg_nme = { 'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' };

%% Get Location Data
for iS = 1:numel(sel_ele)
    sel_loc{iS} = cell(numel(sel_ele{iS}),1);
    rmv_ind = [];
    for iE = 1:numel(sel_ele{iS})
        if ~isempty(find(strcmpi(ecg_loc{1}(:,1),sel_ele{iS}{iE})))
            sel_loc{iS}(iE,1) = ecg_loc{1}(strcmpi(ecg_loc{1}(:,1),sel_ele{iS}{iE}),2);
        else
            rmv_ind = [rmv_ind iE];
        end
    end
    sel_ele{iS}(rmv_ind) = [];
    sel_txt{iS}(rmv_ind) = [];
    sel_loc{iS}(rmv_ind) = [];
end

for iS = 1:numel(sel_ele)
    for iC = 1:numel(fcfg.reg_nme)
        for iE = 1:numel(sel_loc{iS})
            if any(ismember(sel_loc{iS}{iE},fcfg.cmb_reg{iC})); sel_loc{iS}{iE} = fcfg.reg_nme{iC}; end
        end
    end
end

tabulate(sel_loc{1});

tabulate(sel_loc{2});

%% Effects in these patients
% Motor Movement
lng_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_lng_950/subjects/total/pap_lng_950_plt');
    [~,res_ind,~] = intersect(lng_ovr_lap(:,2),sel_ele{2});
    lng_ovr_lap = lng_ovr_lap(res_ind,1:4);
con_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_con_950/subjects/total/pap_con_950_plt');
    con_ovr_lap = con_ovr_lap(res_ind,3:4);
phn_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_phn_950/subjects/total/pap_phn_950_plt');
    phn_ovr_lap = phn_ovr_lap(res_ind,3:4);
mtc_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_mtc_1450/subjects/total/pap_mtc_1450_plt');
    mtc_ovr_lap = mtc_ovr_lap(res_ind,3);

% for iE = 1:numel(sel_ele{2})
%     hld{iE,1} = sel_ele{2}{iE};
%     hld{iE,2} = sel_loc{2}{iE}; 
%     hld(iE,3:4) = lng_ovr_lap(strcmpi(sel_ele{2}{iE},lng_ovr_lap(:,2)),3:4);
% end

[~,~,res_ind] = intersect(lng_ovr_lap(:,2),sel_ele{2});
[sel_ele{2}(res_ind) sel_loc{2}(res_ind) lng_ovr_lap(:,3:4) con_ovr_lap phn_ovr_lap mtc_ovr_lap]

% Language Sensation
lng_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_lng_950/subjects/total/pap_lng_950_plt');
    [~,res_ind,~] = intersect(lng_ovr_lap(:,2),sel_ele{1});
    lng_ovr_lap = lng_ovr_lap(res_ind,1:4);
con_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_con_950/subjects/total/pap_con_950_plt');
    con_ovr_lap = con_ovr_lap(res_ind,3:4);
phn_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_phn_950/subjects/total/pap_phn_950_plt');
    phn_ovr_lap = phn_ovr_lap(res_ind,3:4);
mtc_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_mtc_1450/subjects/total/pap_mtc_1450_plt');
    mtc_ovr_lap = mtc_ovr_lap(res_ind,3);

% for iE = 1:numel(sel_ele{1})
%     hld{iE,1} = sel_ele{1}{iE};
%     hld{iE,2} = sel_loc{1}{iE}; 
%     hld(iE,3:4) = lng_ovr_lap(strcmpi(sel_ele{1}{iE},lng_ovr_lap(:,2)),3:4);
% end

[~,~,res_ind] = intersect(lng_ovr_lap(:,2),sel_ele{1});
[sel_ele{1}(res_ind) sel_loc{1}(res_ind) lng_ovr_lap(:,3:4) con_ovr_lap phn_ovr_lap mtc_ovr_lap]

             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             