stt = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_stimuli_stats.csv');

%% fw_1
lst = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_stimuli_1.csv');
lst(cell2mat(lst(:,3))==7,:) = [];
lst_out_one = cell(size(lst,1),5);

for iS = 1:size(lst,1)
    ind = find(strcmpi(lst{iS,2},stt(:,1)));
    ind(ind==1) = [];
    if ~isempty(ind)
        lst_out_one(iS,:) = stt(ind,[1 2 3 4 7]);
    elseif isempty(ind) && strcmpi(lst{iS,2},'Trigger=6')
        lst_out_one(iS,:) = {lst{iS,2} nan nan nan nan};
    else
        error('')
    end
end

lst_out_one(strcmpi(lst_out_one,'NA')) = cellfun(@(x) strrep(x,'NA','NaN'),lst_out_one(strcmpi(lst_out_one,'NA')),'uni',0);
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_1.csv',lst_out_one)

%% fw_2
lst = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_stimuli_2.csv');
lst(cell2mat(lst(:,3))==7,:) = [];
lst_out_two = cell(size(lst,1),5);

for iS = 1:size(lst,1)
    ind = find(strcmpi(lst{iS,4},stt(:,1)));
    ind(ind==1) = [];
    if ~isempty(ind)
        lst_out_two(iS,:) = stt(ind,[1 2 3 4 7]);
    elseif isempty(ind) && strcmpi(lst{iS,4},'Trigger=6')
        lst_out_two(iS,:) = {lst{iS,4} nan nan nan nan};
    else
        error('')
    end
end

lst_out_two(strcmpi(lst_out_two,'NA')) = cellfun(@(x) strrep(x,'NA','NaN'),lst_out_two(strcmpi(lst_out_two,'NA')),'uni',0);
cell2csv('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_2.csv',lst_out_two)






