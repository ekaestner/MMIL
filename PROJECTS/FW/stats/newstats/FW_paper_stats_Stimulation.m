%% Grab Data
load(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure6/stimulation' '/' 'stm_rsl.mat'])

sel_ele = stm_rsl.sel_ele;
sel_txt = stm_rsl.ele_txt;

ecg_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total' '_' 'split' '_' 'lhs' '_' 'ecog'],'[, \t]');
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



%% Overlap











