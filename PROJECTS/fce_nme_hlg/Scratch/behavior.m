clear; clc

dta_loc = '/home/ekaestner/';

dta_loc_sub = 'Dropbox/McDonald Lab/Erik/Projects/Ideas/fce_nme_stm/PreliminaryData/HalgLab_FN';

out_plt_dir = [ dta_loc '/' dta_loc_sub '/' ];    

out_sht_col = { 'sbj_nme' 'Rp0_acc' 'Rp0_rsp_tme' 'Rp1_acc' 'Rp1_rsp_tme' 'Rp2_acc' 'Rp2_rsp_tme' 'Rp3_acc' 'Rp3_rsp_tme' ...
                'trl_num' 'trl_num_rsp_Rp0' 'trl_num_rsp_Rp1' 'trl_num_rsp_Rp2' 'trl_num_rsp_Rp3' };

sbj_nme = dir([ dta_loc '/' dta_loc_sub]);
    sbj_nme = {sbj_nme(:).name};
    sbj_nme = sbj_nme(string_find(sbj_nme,'DY1'));
    
%% Get data
sbj_nme_col = strcmpi(out_sht_col,'sbj_nme');
Rp0_acc_col = strcmpi(out_sht_col,'Rp0_acc');
Rp0_rsp_tme_col = strcmpi(out_sht_col,'Rp0_rsp_tme');
Rp1_acc_col = strcmpi(out_sht_col,'Rp1_acc');
Rp1_rsp_tme_col = strcmpi(out_sht_col,'Rp1_rsp_tme');
Rp2_acc_col = strcmpi(out_sht_col,'Rp2_acc');
Rp2_rsp_tme_col = strcmpi(out_sht_col,'Rp2_rsp_tme');
Rp3_acc_col = strcmpi(out_sht_col,'Rp3_acc');
Rp3_rsp_tme_col = strcmpi(out_sht_col,'Rp3_rsp_tme');

trl_num_col = strcmpi(out_sht_col,'trl_num');
trl_num_rsp_Rp0_col = strcmpi(out_sht_col,'trl_num_rsp_Rp0');
trl_num_rsp_Rp1_col = strcmpi(out_sht_col,'trl_num_rsp_Rp1');
trl_num_rsp_Rp2_col = strcmpi(out_sht_col,'trl_num_rsp_Rp2');
trl_num_rsp_Rp3_col = strcmpi(out_sht_col,'trl_num_rsp_Rp3');

out_hld = cell(numel(sbj_nme),numel(out_sht_col));

for iS = 1:numel(sbj_nme)
    
    dta_hld = mmil_readtext([ dta_loc '/' dta_loc_sub '/' sbj_nme{iS}]);
    
    out_hld{iS,sbj_nme_col} = sbj_nme{iS}(1:end-4);
    
    tot_ind = strcmpi(dta_hld(:,1),'Rep0') & ~(cell2mat(dta_hld(:,3))==2);
    cor_ind = strcmpi(dta_hld(:,1),'Rep0') & cell2mat(dta_hld(:,3))==1;
    out_hld{iS,Rp0_acc_col}     = round((sum(cor_ind) / sum(tot_ind))*100);
    out_hld{iS,Rp0_rsp_tme_col} = round((sum(cell2mat(dta_hld(cor_ind,4))) / sum(cor_ind))*1000);
    
    tot_ind = strcmpi(dta_hld(:,1),'Rep1') & ~(cell2mat(dta_hld(:,3))==2);
    cor_ind = strcmpi(dta_hld(:,1),'Rep1') & cell2mat(dta_hld(:,3))==1;
    out_hld{iS,Rp1_acc_col}     = round((sum(cor_ind) / sum(tot_ind))*100);
    out_hld{iS,Rp1_rsp_tme_col} = round((sum(cell2mat(dta_hld(cor_ind,4))) / sum(cor_ind))*1000);
    
    tot_ind = strcmpi(dta_hld(:,1),'Rep2') & ~(cell2mat(dta_hld(:,3))==2);
    cor_ind = strcmpi(dta_hld(:,1),'Rep2') & cell2mat(dta_hld(:,3))==1;
    out_hld{iS,Rp2_acc_col}     = round((sum(cor_ind) / sum(tot_ind))*100);
    out_hld{iS,Rp2_rsp_tme_col} = round((sum(cell2mat(dta_hld(cor_ind,4))) / sum(cor_ind))*1000);
    
    tot_ind = strcmpi(dta_hld(:,1),'Rep3') & ~(cell2mat(dta_hld(:,3))==2);
    cor_ind = strcmpi(dta_hld(:,1),'Rep3') & cell2mat(dta_hld(:,3))==1;
    out_hld{iS,Rp3_acc_col}     = round((sum(cor_ind) / sum(tot_ind))*100);
    out_hld{iS,Rp3_rsp_tme_col} = round((sum(cell2mat(dta_hld(cor_ind,4))) / sum(cor_ind))*1000);
    
    out_hld{iS,trl_num_col}     = sum(strcmpi(dta_hld(:,1),'Rep0'));
    out_hld{iS,trl_num_rsp_Rp0_col} = sum(strcmpi(dta_hld(:,1),'Rep0') & ~(cell2mat(dta_hld(:,3))==2));
    out_hld{iS,trl_num_rsp_Rp1_col} = sum(strcmpi(dta_hld(:,1),'Rep1') & ~(cell2mat(dta_hld(:,3))==2));
    out_hld{iS,trl_num_rsp_Rp2_col} = sum(strcmpi(dta_hld(:,1),'Rep2') & ~(cell2mat(dta_hld(:,3))==2));
    out_hld{iS,trl_num_rsp_Rp3_col} = sum(strcmpi(dta_hld(:,1),'Rep3') & ~(cell2mat(dta_hld(:,3))==2));
    
end

%% Make figures
dst_col = distinguishable_colors(numel(sbj_nme));

% Accuracy plot
fcfg = [];

for iS = 1:numel(sbj_nme)
    fcfg.xdt{iS}     = [1:4]+rand(1)/10;
    fcfg.ydt{iS}     = [out_hld{iS,[find(Rp0_acc_col) find(Rp1_acc_col) find(Rp2_acc_col) find(Rp3_acc_col)]}];
    fcfg.fce_col{iS} = dst_col(iS,:);
end

fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { '~5s ISI' '~45s ISI' '~45s ISI' '~45s ISI' };
fcfg.ylb = { '% Correct'  };

fcfg.hln     = 0;
fcfg.hln_col = rgb('black');

fcfg.xlm = [ 0 5 ];
fcfg.ylm = [ -10 110 ];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Accuracy';

ejk_scatter(fcfg)

% Response time plot
fcfg = [];

for iS = 1:numel(sbj_nme)
    fcfg.xdt{iS}     = [1:4]+rand(1)/10;
    fcfg.ydt{iS}     = [out_hld{iS,[find(Rp0_rsp_tme_col) find(Rp1_rsp_tme_col) find(Rp2_rsp_tme_col) find(Rp3_rsp_tme_col)]}];
    fcfg.fce_col{iS} = dst_col(iS,:);
end

fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'E0 (~5s ISI)' 'E1 (~45s ISI)' 'E2 (~45s ISI)' 'E3 (~45s ISI)' };
fcfg.ylb = { 'ResponseTime'  };

fcfg.hln     = 0;
fcfg.hln_col = rgb('black');

fcfg.xlm = [ 0 5 ];
fcfg.ylm = [ -10 4500 ];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'ResponseTime';

ejk_scatter(fcfg)

% Trial numbers
fcfg = [];

for iS = 1:numel(sbj_nme)
    fcfg.xdt{iS}     = [1:5]+rand(1)/10;
    fcfg.ydt{iS}     = [out_hld{iS,[find(trl_num_col) find(trl_num_rsp_Rp0_col) find(trl_num_rsp_Rp1_col) find(trl_num_rsp_Rp2_col) find(trl_num_rsp_Rp3_col)]}];
    fcfg.fce_col{iS} = dst_col(iS,:);
end

fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'Total Trials' 'E0 Responses' 'E1 Responses' 'E2 Responses' 'E3 Responses' };
fcfg.ylb = { 'ResponseTime'  };

fcfg.hln     = 0;
fcfg.hln_col = rgb('black');

fcfg.xlm = [ 0 6 ];
fcfg.ylm = [ -5 70 ];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'TrialNumbers';

ejk_scatter(fcfg)







