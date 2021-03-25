clear; clc;

clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

lsa_lst = mmil_readtext([clr_fld '/' 'task' '/' 'SAi-actual.csv']);
lsz_lst = mmil_readtext([clr_fld '/' 'task' '/' 'SZi-actual.csv']);

lsa_lst_wrd = tabulate(lsa_lst(:,2)); lsa_lst_wrd = lsa_lst_wrd(cell2mat(lsa_lst_wrd(:,2))==40,1); 
for iW = 1:numel(lsa_lst_wrd); 
    num_rep.(lsa_lst_wrd{iW}) = diff(find(strcmpi(lsa_lst(:,2),lsa_lst_wrd{iW})));
    num_avg.(lsa_lst_wrd{iW}) = mean(num_rep.(lsa_lst_wrd{iW}));
    avg_hld(iW) = num_avg.(lsa_lst_wrd{iW});
end

lsz_lst_wrd = tabulate(lsz_lst(:,2)); lsz_lst_wrd = lsz_lst_wrd(cell2mat(lsz_lst_wrd(:,2))==20,1);
for iW = 1:numel(lsz_lst_wrd);
    num_rep.(lsz_lst_wrd{iW}) = diff(find(strcmpi(lsz_lst(:,2),lsz_lst_wrd{iW})));
    num_avg.(lsz_lst_wrd{iW}) = mean(num_rep.(lsz_lst_wrd{iW}));
    avg_hld(iW+10) = num_avg.(lsz_lst_wrd{iW});
end

