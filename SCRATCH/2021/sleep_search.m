clear; clc;

ttt = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Complete/slp_emo_mem/AffNap_Data_dprime.csv');

hld_dta = cell2mat(ttt(2:end,[ 6 7 8]));
hld_dta(:,4) = mean(hld_dta(:,1:2),2);

tst_num = 1;
hld_out = [];
for iI = 1:1000000
    
    smp_num = randsample(28,24);
       
    neu_men = mean(hld_dta(smp_num,3));
    neu_std = std(hld_dta(smp_num,3));
    
    pos_men = mean(hld_dta(smp_num,4));
    pos_std = std(hld_dta(smp_num,4));
    
    
    if (roundsd(neu_men,3)==1.84 && roundsd(neu_std,2)==0.68)
        hld_out = [ hld_out setxor(1:28,smp_num) ];
    end

end

tabulate(hld_out(:))

% 1, 5, 15, 20