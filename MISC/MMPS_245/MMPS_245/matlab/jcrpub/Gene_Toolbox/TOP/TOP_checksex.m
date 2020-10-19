function [sex_prob_subjID] = TOP_checksex(sexchkfile)

if ~exist(sexchkfile,'file'),
   error(sprintf('%s does not exist!',sexchkfile));
%   cmd = sprintf('plink --bfile %s --check-sex --out %s', basefilename);
%   [cmd_stat,result] = system(cmd);
end

Sex_chk = importdata(sexchkfile,' ');

stat_col = find(find_cell(strfind(Sex_chk.textdata(1,:),'STATUS')));

prob_subj_ind = find(find_cell(strfind(Sex_chk.textdata(:,stat_col),'PROBLEM')));

subjID_col = find(find_cell(strfind(Sex_chk.textdata(1,:),'IID')));

sex_prob_subjID = str2double(Sex_chk.textdata(prob_subj_ind,subjID_col));
