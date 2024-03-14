
%% Definitions
ini_sbj_fle = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data/subject_check/Total_Check.csv';

ejk_chk_dir(qal_dir);

qal_col_nme = { 'dti_pst_prc_qal' 'initial_dti_post_proc_quality' ; ...
                'low_pre_cog'     'initial_low_pre_cognitive' };

qal_fld_nme = { 'srg_sde_ana' 'initial' ; ...
                'srg_sde_ana' 'eng_one' ; ...
                'srg_sde_ana' 'eng_two' ; ...
                'srg_sde_ana' 'lng_lft' ; ...
                'sde_ana'     'initial' ; ...
                'sde_ana'     'eng_one' ; ...
                'sde_ana'     'eng_two' ; ...
                'sde_ana'     'lng_lft' };          
            
% %% Find initial subjects
% ini_sbj_hld = mmil_readtext(ini_sbj_fle);
% ini_sbj = ini_sbj_hld( strcmpi(ini_sbj_hld(:,2),'INCLUDED'), 1); ini_sbj(1,:) = [];
% 
% load([ dta_dir '/' 'group_quality_check.mat']);... 'group_initial.mat' ])
% 
% fcfg = [];
% fcfg.dta_loc = [ qal_dir '/' 'initial_QC_overall_ejk.csv'];
% fcfg.dta_col = 2;
% [ ini_qal_dta, ini_qal_sbj, ini_qal_col] = ejk_dta_frm( fcfg );
% 
% fcfg = [];
% fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
% fcfg.dta_col = 2;
% [ all_dta, all_sbj, all_col] = ejk_dta_frm( fcfg );
% 
% new_ind = [ grp.srg_sde.non_qc.lft_atl ; grp.srg_sde.non_qc.rgh_atl ; grp.srg_sde.non_qc.lft_slh ; grp.srg_sde.non_qc.rgh_slh ];
% new_dta = all_dta(new_ind,:);
% new_sbj = all_sbj(new_ind,:);
% new_col = all_col;
% 
% %% List of subjects
% tot_sbj = unique([ ini_sbj ; new_sbj ]);
% 
% out_csv = cell(numel(tot_sbj)+1,3+size(qal_col_nme,1));
% out_csv(:,4:3+size(qal_col_nme,1)) = {NaN};
% for iS = 1:numel(tot_sbj)
%     ini_exs = any(strcmpi(ini_sbj,tot_sbj{iS}));
%     new_exs = any(strcmpi(new_sbj,tot_sbj{iS}));
%     if ini_exs; out_csv{iS,1} = tot_sbj{iS}; end
%     if new_exs; out_csv{iS,2} = tot_sbj{iS}; end
%     if new_exs; out_csv{iS,3} = 1; else out_csv{iS,3} = 0; end
%     
%     ini_qal_ind = strcmpi(ini_qal_sbj,tot_sbj{iS});
%     if any(ini_qal_ind)
%         for iI = 1:size(qal_col_nme,1)
%             ini_col_ind = strcmpi(ini_qal_col,qal_col_nme{iI,1});
%             out_csv{iS,3+iI} = ini_qal_dta{ini_qal_ind,ini_col_ind};
%         end
%     end
%     
% end
% 
% tot_sbj_num = size(tot_sbj,1);
% new_sbj_num = sum(cellfun(@isempty,out_csv(:,1)));
% old_sbj_num = sum(cellfun(@isempty,out_csv(:,2)));
% inc_sbj_num = size(new_sbj,1);
% ini_sbj_num = size(ini_sbj,1);
% 
% out_csv{end,1} = sprintf('%i subjects; %i included; %i new; %i removed',tot_sbj_num, inc_sbj_num, new_sbj_num, old_sbj_num);
% cell2csv( [ qal_dir '/' 'subject_list.csv' ] , [ {'InitialSubjects'} {'CurrentSubjects'} {'Include'} qal_col_nme(:,2)' ; out_csv ])
% 
% %% Run DTI QC
% ejk_qc_roi

%% Apply QC
load([ dta_dir '/' 'group_initial.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ new_dta, new_sbj, new_col] = ejk_dta_frm( fcfg );

fcfg         = [];
fcfg.dta_loc = [ qal_dir '/' 'subject_list_ejk.csv' ];
fcfg.dta_col = 2;
[ app_qal_dta, app_qal_sbj, app_qal_col] = ejk_dta_frm( fcfg );

rmv_nme = app_qal_sbj(strcmpi(app_qal_dta( :, strcmpi(app_qal_col,'Plan') ),'Remove'));

for iQ = 1:size(qal_fld_nme,1)
    grp.(qal_fld_nme{iQ,1}).([qal_fld_nme{iQ,2} '_' 'QC']) = grp.(qal_fld_nme{iQ,1}).(qal_fld_nme{iQ,2});
    fld_nme = fieldnames(grp.(qal_fld_nme{iQ,1}).([qal_fld_nme{iQ,2} '_' 'QC']));
    
    for iF = 1:numel(fld_nme)
        cur_sbj = new_sbj(grp.(qal_fld_nme{iQ,1}).([qal_fld_nme{iQ,2} '_' 'QC']).(fld_nme{iF}));
        [~, rmv_sbj_ind ] = intersect(cur_sbj,rmv_nme);
        grp.(qal_fld_nme{iQ,1}).([qal_fld_nme{iQ,2} '_' 'QC']).(fld_nme{iF})(rmv_sbj_ind,:) = [];
    end
    
end

save([ dta_dir '/' 'group_quality_check.mat' ],'grp')

%%
ord_out = [ grp.srg_sde_ana.initial_QC.lft_atl ; grp.srg_sde_ana.initial_QC.rgh_atl ; grp.srg_sde_ana.initial_QC.lft_slh ; grp.srg_sde_ana.initial_QC.rgh_slh ];

cell2csv([ dta_dir '/' 'QC_Sample_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ] , [ 'sbj_nme' new_col ; new_sbj(ord_out,1) new_dta(ord_out,:)]);

