clear; clc;

dta_loc = '/home/ekaestne/PROJECTS/OUTPUT/SL/sig_chn/hgp/ecog/split';

sbj_nme = { 'NY439_SL' 'NY523_SL' 'NY590_SL' 'NY591_SL' 'NY598_SL' };

eff_nme = {   'pap_anv_1500'                              'pap_lng_950'                           'pap_phn_950'                        'pap_mtc_1450' };
eff_imp = { { 'VisualSelective' 'AuditorySelective' }   { 'VisualLanguage' 'AuditoryLanguage' } { 'VisualLetter' 'AuditoryPhoneme' } { 'Mismatch_sensory' } };


lbl_use = {'Fusiform' ...
           'STG' ...
           'Supramarginal' ...
           'Precentral' ...                    
           'Pars Opercularis' };
lbl_nme = {'Fusiform' ...
           'Lateral Occipital' ...
           'Precentral' ...
           'STG' ...
           'MTG' ...
           'Supramarginal' ...
           'Pars Opercularis' };
lbl_inc = { { 'lhs_caudal-fusiform' 'lhs_middle-fusiform' }...
            { 'lhs_lateraloccipital' } ...
            { 'lhs_inferior-precentral' 'lhs_middle-precentral' } ...
            { 'lhs_caudal-STG' 'lhs_middle-STG' } ...
            { 'lhs_caudal-MTG' 'lhs_middle-MTG' } ...
            { 'lhs_supramarginal' } ...
            { 'lhs_parsopercularis' } };

sbj_col = distinguishable_colors(9);
    sbj_col = flipud(sbj_col([1:7 9],:));
    sbj_col = sbj_col([1 2 6 7 8], :);
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iE = 1:numel(eff_nme)
    tbl_out.(eff_nme{iE})= [];
    tbl_tot.(eff_nme{iE})= [];
    
    for iE2 = 1:numel(eff_imp{iE})        
        tbl_out.(eff_nme{iE}).(eff_imp{iE}{iE2}) = [];
        tbl_tot.(eff_nme{iE}).(eff_imp{iE}{iE2}) = [];
                
        for  iL = 1:numel(lbl_nme)
            tbl_out.(eff_nme{iE}).(eff_imp{iE}{iE2}){1,iL+1} = lbl_nme{iL};
            tbl_tot.(eff_nme{iE}).(eff_imp{iE}{iE2}){1,iL+1} = lbl_nme{iL};
        
            for iS = 1:numel(sbj_nme)
                sbj_dta_hld = mmil_readtext( [ dta_loc '/' eff_nme{iE} '/' 'table' '/' 'lhs' '/' sbj_nme{iS} '_' eff_nme{iE} '_' 'lhs_table' ]);
                tbl_out.(eff_nme{iE}).(eff_imp{iE}{iE2}){iS+1,1} = sbj_nme{iS};
                tbl_tot.(eff_nme{iE}).(eff_imp{iE}{iE2}){iS+1,1} = sbj_nme{iS};
                out_eff_hld = 0;
                out_tot_hld = 0;
                
                for iL2 = 1:numel(lbl_inc{iL})                    
                    eff_col = find(strcmpi( sbj_dta_hld(1,:), eff_imp{iE}{iE2}));
                    row_ind = find(strcmpi( sbj_dta_hld(:,1), lbl_inc{iL}{iL2}));
                    
                    if ~isempty(row_ind)
                        out_eff_hld = out_eff_hld + sbj_dta_hld{ row_ind, eff_col};
                        out_tot_hld = out_tot_hld + sbj_dta_hld{ row_ind, eff_col+1};
                    else
                        out_eff_hld = out_eff_hld + 0;
                        out_tot_hld = out_tot_hld + 0;
                    end
                    
                end
                
                tbl_out.(eff_nme{iE}).(eff_imp{iE}{iE2}){iS+1,iL+1} = out_eff_hld;
                tbl_tot.(eff_nme{iE}).(eff_imp{iE}{iE2}){iS+1,iL+1} = out_tot_hld;
                
            end
            
        end
        
        
    end
end

%% Make Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eff_imp_hld = fieldnames(tbl_out);

for iEI = 1:numel(eff_imp_hld)
    eff_imp_use = eff_imp_hld{iEI};
    eff_nme_hld = fieldnames(tbl_out.(eff_imp_use));
    
    for iEN = 1:numel(eff_nme_hld)       
        eff_nme_use = eff_nme_hld{iEN};
        
        use_hld     = tbl_out.(eff_imp_use).(eff_nme_use);
        tot_hld     = tbl_tot.(eff_imp_use).(eff_nme_use);
        cnt = 1;
        for iL = 1:numel(lbl_use)
            for iS = 1:numel(sbj_nme)
                
                if ~(tot_hld{ strcmpi(tot_hld(:,1), sbj_nme{iS}), strcmpi(tot_hld(1,:), lbl_use{iL})}==0)
                    ydt{cnt} = use_hld{ strcmpi(use_hld(:,1), sbj_nme{iS}), strcmpi(use_hld(1,:), lbl_use{iL}) } / tot_hld{ strcmpi(tot_hld(:,1), sbj_nme{iS}), strcmpi(tot_hld(1,:), lbl_use{iL}) };
                    xdt{cnt} = iL;
                    
                    fce_col{cnt} = sbj_col(iS,:);
                    edg_col{cnt} = rgb('black');
                    
                    cnt = cnt+1;
                end
                
            end
        end
        
        %
        fcfg = [];
        
        fcfg.ydt     = ydt;
        fcfg.xdt     = xdt;
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = edg_col;
        
        fcfg.xlb = lbl_use;
        fcfg.ylb = {[ eff_nme_use ' (Proportion)']};
        
        fcfg.ylm = [0 1];
        
        fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/SL/revision/proportions';
        fcfg.out_nme = [ eff_imp_use '_' eff_nme_use];
        
        ejk_scatter(fcfg)
        
        clear ydt xdt fce_col edg_col
        
    end
end

%% Create Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt_col = 1;
for iEI = 1:numel(eff_imp_hld)
    eff_imp_use = eff_imp_hld{iEI};
    eff_nme_hld = fieldnames(tbl_out.(eff_imp_use));
    
    for iEN = 1:numel(eff_nme_hld)
        eff_nme_use = eff_nme_hld{iEN};
        tbl_prp{1, cnt_col+1} = eff_nme_use;
        tbl_sbj{1, cnt_col+1} = eff_nme_use;
        
        use_hld     = tbl_out.(eff_imp_use).(eff_nme_use);
        tot_hld     = tbl_tot.(eff_imp_use).(eff_nme_use);
        
        for iL = 1:numel(lbl_use)
            
            tbl_prp{iL+1, 1} = lbl_use{iL};
            tbl_sbj{iL+1, 1} = lbl_use{iL};
            
            fst = num2str(sum([use_hld{ 2:end, strcmpi(use_hld(1,:), lbl_use{iL}) }] > 0));
            scd = num2str(sum([tot_hld{ 2:end, strcmpi(tot_hld(1,:), lbl_use{iL}) }] > 0));
            
            tbl_prp{iL+1, cnt_col+1} = sum([use_hld{ 2:end, strcmpi(use_hld(1,:), lbl_use{iL}) }]) / sum([tot_hld{ 2:end, strcmpi(tot_hld(1,:), lbl_use{iL}) }]);
            tbl_sbj{iL+1, cnt_col+1} = [ fst ' of ' scd];
            
        end
        cnt_col = cnt_col + 1;
    end
end

cell2csv('/home/ekaestne/PROJECTS/OUTPUT/SL/revision/proportions/proportion_table.csv',tbl_prp);
cell2csv('/home/ekaestne/PROJECTS/OUTPUT/SL/revision/proportions/subject_table.csv',tbl_sbj);









