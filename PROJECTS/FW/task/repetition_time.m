ttt = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/task/fw_stimuli_1.csv');

rep_mat = unique(ttt(cell2mat(ttt(:,3))==4,2));

for iW = 1:numel(rep_mat)
    rep_num{iW} = diff(find(strcmpi(rep_mat{iW},ttt(:,2))));
end

mean(cat(1,rep_num{:}))

min(cat(1,rep_num{:}))
max(cat(1,rep_num{:}))

