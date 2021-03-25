%% Mapping
map = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/clinical/Mapping/mapping.csv');
ecg_loc{1} = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/electrode_location_files/total/output/total_split_lhs_ecog','[, \t]');
ecg_loc{1} = ecg_loc{1}(:,[1 5:end]);
ecg_loc{2} = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/electrode_location_files/total/output/total_split_rhs_ecog','[, \t]');
ecg_loc{2} = ecg_loc{2}(:,[1 5:end]);

for iH = 1:numel(ecg_hms)
    if ~isempty(ecg_loc{iH})
        for iR = 1:size(ecg_loc{iH},1)
            if isempty(ecg_loc{iH}{iR,2})
                ecg_loc{iH}{iR,2} = 'unknown';
            else
                num_hld = cellfun(@(x) str2num(x(1:end-1)),ecg_loc{iH}(iR,string_find(ecg_loc{iH}(iR,~cellfun(@isempty,ecg_loc{iH}(iR,:))),'%')),'uni',0);
                [~,max_hld] = max([num_hld{:}]);
                ecg_loc{iH}{iR,2} = ecg_loc{iH}{iR,max_hld+1+(1*max_hld)};
            end
        end
        ecg_loc{iH} = ecg_loc{iH}(:,1:2);
    end
    
end

ecg_loc{1} = [ecg_loc{1} ; ecg_loc{2}];

           %%
loc_mid_pre = find(strcmpi(ecg_loc{1}(:,2),'middle-precentral'));
map_ele_hld = strcat(map(:,1),'_XS_',map(:,2));
[~,map_ind,~] = intersect(map_ele_hld,ecg_loc{1}(loc_mid_pre,1));
[map(map_ind,1) map(map_ind,2) map(map_ind,5) map(map_ind,31)] 

loc_mid_pre = find(strcmpi(ecg_loc{1}(:,2),'superior-precentral'));
map_ele_hld = strcat(map(:,1),'_XS_',map(:,2));
[~,map_ind,~] = intersect(map_ele_hld,ecg_loc{1}(loc_mid_pre,1));
[map(map_ind,1) map(map_ind,2) map(map_ind,5) map(map_ind,31)] 

loc_mid_pre = find(strcmpi(ecg_loc{1}(:,2),'inferior-precentral'));
map_ele_hld = strcat(map(:,1),'_XS_',map(:,2));
[~,map_ind,~] = intersect(map_ele_hld,ecg_loc{1}(loc_mid_pre,1));
[map(map_ind,1) map(map_ind,2) map(map_ind,5) map(map_ind,31)] 


