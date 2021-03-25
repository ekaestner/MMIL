clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

sbj_nme = mmil_readtext([clr_fld '/' 'subjects']);

%%
for iS = 1:numel(sbj_nme)
  
    ele_idn = mmil_readtext([ clr_fld '/' 'ele_idn' '/' sbj_nme{iS} '/' sbj_nme{iS} '.csv']);
    cell2csv([ clr_fld '/' 'ele_idn' '/' sbj_nme{iS} '/' sbj_nme{iS} '_include.csv'],[ele_idn num2cell(ones(size(ele_idn,1),1))]);    
    
end
