%% Check tasks for each participant
sbj_inf_nme = dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sbj_inf');
sbj_inf_nme = {sbj_inf_nme(3:end-1).name};
sbj_inf_nme(string_find(sbj_inf_nme,'~')) = [];

for iS = 1:numel(sbj_inf_nme)
    tsk         = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/' '/' 'sbj_inf' '/' sbj_inf_nme{iS}],'tsk');
    if ~isempty(tsk{1})
        sbj_nme{iS,1} = sbj_inf_nme{iS};
        tsk = unique(tsk);
        sbj_nme{iS,2} = [tsk{:}];
    end
end

%% Fix fsaverage locations
fix_ele_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/electrode_location_files/total/output_fixed/' '/' 'total_lhs_ecog']);
new_ele_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/electrode_location_files/total/output/' '/' 'total_lhs_ecog']);

for iE = 1:size(new_ele_loc,1)
    
    old_ele_loc = find(strcmpi(fix_ele_loc,new_ele_loc{iE,1}));
    if new_ele_loc{iE,2}~=fix_ele_loc{old_ele_loc,2} && ...
       new_ele_loc{iE,3}~=fix_ele_loc{old_ele_loc,3} && ...
       new_ele_loc{iE,4}~=fix_ele_loc{old_ele_loc,4}
   
        new_ele_loc{iE,2} = fix_ele_loc{old_ele_loc,2};
        new_ele_loc{iE,3} = fix_ele_loc{old_ele_loc,3};
        new_ele_loc{iE,4} = fix_ele_loc{old_ele_loc,4};
        fprintf('%s\n', new_ele_loc{iE,1})
       
    end
end

cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/electrode_location_files/total/output/' '/' 'total_lhs_ecog'],new_ele_loc)

%% Make Subtraction Times & Amplitude
% Times %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HGP
% Act
aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_act/subjects/total' '/' 'pap_aud_act' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_act/subjects/total' '/' 'pap_vis_act' '_' 'fst' ]);

bim_fst(:,1:2) = aud_fst(:,1:2);
bim_fst(:,3)   = {'bim'};
for iR = 2:size(bim_fst,1)
    if ~isnan(aud_fst{iR,3}) && ~isnan(vis_fst{iR,3})
        bim_fst{iR,3} = vis_fst{iR,3}-aud_fst{iR,3};
    else
        bim_fst{iR,3} = NaN;
    end
end

cell2csv([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_bim_act/subjects/total' '/' 'pap_bim_act' '_' 'fst' ],bim_fst)

% Rep
aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_rep/subjects/total' '/' 'pap_aud_rep' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_rep/subjects/total' '/' 'pap_vis_rep' '_' 'fst' ]);

bim_fst(:,1:2) = aud_fst(:,1:2);
bim_fst(:,3)   = {'bim'};
for iR = 2:size(bim_fst,1)
    if ~isnan(aud_fst{iR,3}) && ~isnan(vis_fst{iR,3})
        bim_fst{iR,3} = aud_fst{iR,3}-vis_fst{iR,3};
        fprintf(['%s : %3.3f\n'],bim_fst{iR,2},bim_fst{iR,3}*100)
    else
        bim_fst{iR,3} = NaN;
    end
end

cell2csv([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_bim_rep/subjects/total' '/' 'pap_bim_rep' '_' 'fst' ],bim_fst)

% LFP
% Act
aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_act/subjects/total' '/' 'pap_aud_act' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_act/subjects/total' '/' 'pap_vis_act' '_' 'fst' ]);

bim_fst(:,1:2) = aud_fst(:,1:2);
bim_fst(:,3)   = {'bim'};
for iR = 2:size(bim_fst,1)
    if ~isnan(aud_fst{iR,3}) && ~isnan(vis_fst{iR,3})
        bim_fst{iR,3} = vis_fst{iR,3}-aud_fst{iR,3};
    else
        bim_fst{iR,3} = NaN;
    end
end

cell2csv([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_bim_act/subjects/total' '/' 'pap_bim_act' '_' 'fst' ],bim_fst)

% Rep
aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_rep_nob/subjects/total' '/' 'pap_aud_rep_nob' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_rep_nob/subjects/total' '/' 'pap_vis_rep_nob' '_' 'fst' ]);

bim_fst(:,1:2) = aud_fst(:,1:2);
bim_fst(:,3)   = {'bim'};
for iR = 2:size(bim_fst,1)
    if ~isnan(aud_fst{iR,3}) && ~isnan(vis_fst{iR,3})
        bim_fst{iR,3} = vis_fst{iR,3}-aud_fst{iR,3};
    else
        bim_fst{iR,3} = NaN;
    end
end

cell2csv([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_bim_rep_nob/subjects/total' '/' 'pap_bim_rep_nob' '_' 'fst' ],bim_fst)

% Amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%
aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/amp' '/' 'pap_aud_act' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/amp' '/' 'pap_vis_act' ]);

bim_fst(:,1:2) = aud_fst(:,1:2);
bim_fst(:,3)   = {'bim'};
for iR = 2:size(bim_fst,1)
    if ~isnan(aud_fst{iR,3}) && ~isnan(vis_fst{iR,3})
        bim_fst{iR,3} = aud_fst{iR,3}- vis_fst{iR,3};
    else
        bim_fst{iR,3} = NaN;
    end
end

cell2csv([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/amp' '/' 'pap_bim_act' ],bim_fst)

%% Make Bimodal comparison
% ACTIVATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportion Responsive - Activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_act/subjects/total' '/' 'pap_aud_act' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_act/subjects/total' '/' 'pap_vis_act' '_' 'plt' ]);

%
aud_bim_plt(:,1:2) = aud_plt(:,1:2);
aud_bim_plt(1,3)   = {'UniModal-Auditory'};
aud_bim_plt(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_plt{iR,3} = 1;
        aud_bim_plt{iR,4} = 0;
    else
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 0;
    end
end

%
vis_bim_plt(:,1:2) = aud_plt(:,1:2);
vis_bim_plt(1,3)   = {'UniModal-Visual'};
vis_bim_plt(1,4)   = {'BiModal-Visual'};
vis_bim_plt(1,5)   = {'UniModal-Overall'};
for iR = 2:size(vis_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 1;
        vis_bim_plt{iR,4} = 0;
    else
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 0;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) || (aud_plt{iR,3}==1 && vis_plt{iR,3}==0)
        vis_bim_plt{iR,5} = 1;
    else
        vis_bim_plt{iR,5} = 0;
    end
end

cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/subjects/total' '/' 'pap_aud_bim_act' '_' 'plt' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/subjects/total' '/' 'pap_vis_bim_act' '_' 'plt' ],vis_bim_plt)

%
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/subjects/total' '/' 'pap_aud_bim_act' '_' 'tbl' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/subjects/total' '/' 'pap_vis_bim_act' '_' 'tbl' ],vis_bim_plt)

tcfg = [];

tcfg.sbj_clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
tcfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ//epoch_data';
tcfg.sbj_nme = 'total';
tcfg.typ     = {'hgp'};
tcfg.ele_typ = {'ecog'};
tcfg.cmb_nme = {'pap_vis_bim_act'};
tcfg.loc_typ = {'split'};
mmil_sig_cmb_tbl(tcfg)

tcfg.cmb_nme = {'pap_aud_bim_act'};
mmil_sig_cmb_tbl(tcfg)

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/table/lhs/total_pap_vis_bim_act_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/total/total_pap_vis_bim_act_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/table/rhs/total_pap_vis_bim_act_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/total/total_pap_vis_bim_act_rhs_table_plot')
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/table/lhs/total_pap_aud_bim_act_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/total/total_pap_aud_bim_act_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/table/rhs/total_pap_aud_bim_act_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/total/total_pap_aud_bim_act_rhs_table_plot')

% Timing - Activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_act/subjects/total' '/' 'pap_aud_act' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_act/subjects/total' '/' 'pap_vis_act' '_' 'plt' ]);

aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_act/subjects/total' '/' 'pap_aud_act' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_act/subjects/total' '/' 'pap_vis_act' '_' 'fst' ]);

%
aud_bim_fst(:,1:2) = aud_fst(:,1:2);
aud_bim_fst(1,3)   = {'UniModal-Auditory'};
aud_bim_fst(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = aud_fst{iR,3};
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_fst{iR,3} = aud_fst{iR,3};
        aud_bim_fst{iR,4} = NaN;
    else
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = NaN;
    end
end

%
vis_bim_fst(:,1:2) = aud_fst(:,1:2);
vis_bim_fst(1,3)   = {'UniModal-Visual'};
vis_bim_fst(1,4)   = {'BiModal-Visual'};
for iR = 2:size(vis_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = vis_fst{iR,3};
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = vis_fst{iR,3};
        vis_bim_fst{iR,4} = NaN;
    else
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = NaN;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = vis_fst{iR,3};
    elseif (aud_plt{iR,3}==1 && vis_plt{iR,3}==0) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = aud_fst{iR,3};
    else
        vis_bim_fst{iR,5} = NaN;
    end
end

mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_act/subjects/total' '/' 'pap_aud_bim_act' '_' 'fst' ],aud_bim_fst)
mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_act/subjects/total' '/' 'pap_vis_bim_act' '_' 'fst' ],vis_bim_fst)
    
% REPETITION HGP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportion Repetition HGP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_rep/subjects/total' '/' 'pap_aud_rep' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_rep/subjects/total' '/' 'pap_vis_rep' '_' 'plt' ]);

%
aud_bim_plt(:,1:2) = aud_plt(:,1:2);
aud_bim_plt(1,3)   = {'UniModal-Auditory'};
aud_bim_plt(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_plt{iR,3} = 1;
        aud_bim_plt{iR,4} = 0;
    else
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 0;
    end
end

%
vis_bim_plt(:,1:2) = aud_plt(:,1:2);
vis_bim_plt(1,3)   = {'UniModal-Visual'};
vis_bim_plt(1,4)   = {'BiModal-Visual'};
vis_bim_plt(1,5)   = {'UniModal-Overall'};
for iR = 2:size(vis_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 1;
        vis_bim_plt{iR,4} = 0;
    else
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 0;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) || (aud_plt{iR,3}==1 && vis_plt{iR,3}==0)
        vis_bim_plt{iR,5} = 1;
    else
        vis_bim_plt{iR,5} = 0;
    end
end

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/subjects/total')
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/subjects/total')

cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/subjects/total' '/' 'pap_aud_bim_rep' '_' 'plt' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/subjects/total' '/' 'pap_vis_bim_rep' '_' 'plt' ],vis_bim_plt)

%
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/subjects/total' '/' 'pap_aud_bim_rep' '_' 'tbl' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/subjects/total' '/' 'pap_vis_bim_rep' '_' 'tbl' ],vis_bim_plt)

tcfg = [];

tcfg.sbj_clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
tcfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ//epoch_data';
tcfg.sbj_nme = 'total';
tcfg.typ     = {'hgp'};
tcfg.ele_typ = {'ecog'};
tcfg.cmb_nme = {'pap_vis_bim_rep'};
tcfg.loc_typ = {'split'};
mmil_sig_cmb_tbl(tcfg)

tcfg.cmb_nme = {'pap_aud_bim_rep'};
mmil_sig_cmb_tbl(tcfg)

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/table/lhs/total_pap_vis_bim_rep_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/total/total_pap_vis_bim_rep_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/table/rhs/total_pap_vis_bim_rep_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/total/total_pap_vis_bim_rep_rhs_table_plot')
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/table/lhs/total_pap_aud_bim_rep_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/total/total_pap_aud_bim_rep_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/table/rhs/total_pap_aud_bim_rep_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/total/total_pap_aud_bim_rep_rhs_table_plot')

% Timing - Repetition HGP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_rep/subjects/total' '/' 'pap_aud_rep' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_rep/subjects/total' '/' 'pap_vis_rep' '_' 'plt' ]);

aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_rep/subjects/total' '/' 'pap_aud_rep' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_rep/subjects/total' '/' 'pap_vis_rep' '_' 'fst' ]);

%
aud_bim_fst(:,1:2) = aud_fst(:,1:2);
aud_bim_fst(1,3)   = {'UniModal-Auditory'};
aud_bim_fst(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = aud_fst{iR,3};
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_fst{iR,3} = aud_fst{iR,3};
        aud_bim_fst{iR,4} = NaN;
    else
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = NaN;
    end
end

%
vis_bim_fst(:,1:2) = aud_fst(:,1:2);
vis_bim_fst(1,3)   = {'UniModal-Visual'};
vis_bim_fst(1,4)   = {'BiModal-Visual'};
for iR = 2:size(vis_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = vis_fst{iR,3};
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = vis_fst{iR,3};
        vis_bim_fst{iR,4} = NaN;
    else
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = NaN;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = vis_fst{iR,3};
    elseif (aud_plt{iR,3}==1 && vis_plt{iR,3}==0) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = aud_fst{iR,3};
    else
        vis_bim_fst{iR,5} = NaN;
    end
end

mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_bim_rep/subjects/total' '/' 'pap_aud_bim_rep' '_' 'fst' ],aud_bim_fst)
mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_bim_rep/subjects/total' '/' 'pap_vis_bim_rep' '_' 'fst' ],vis_bim_fst)
    
% REPETITION LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportion Repetition lfp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_rep_nob/subjects/total' '/' 'pap_aud_rep_nob' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_rep_nob/subjects/total' '/' 'pap_vis_rep_nob' '_' 'plt' ]);

%
aud_bim_plt(:,1:2) = aud_plt(:,1:2);
aud_bim_plt(1,3)   = {'UniModal-Auditory'};
aud_bim_plt(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_plt{iR,3} = 1;
        aud_bim_plt{iR,4} = 0;
    else
        aud_bim_plt{iR,3} = 0;
        aud_bim_plt{iR,4} = 0;
    end
end

%
vis_bim_plt(:,1:2) = aud_plt(:,1:2);
vis_bim_plt(1,3)   = {'UniModal-Visual'};
vis_bim_plt(1,4)   = {'BiModal-Visual'};
vis_bim_plt(1,5)   = {'UniModal-Overall'};
for iR = 2:size(vis_bim_plt,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 1;
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_plt{iR,3} = 1;
        vis_bim_plt{iR,4} = 0;
    else
        vis_bim_plt{iR,3} = 0;
        vis_bim_plt{iR,4} = 0;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) || (aud_plt{iR,3}==1 && vis_plt{iR,3}==0)
        vis_bim_plt{iR,5} = 1;
    else
        vis_bim_plt{iR,5} = 0;
    end
end

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/subjects/total')
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/subjects/total')

cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/subjects/total' '/' 'pap_aud_bim_rep_nob' '_' 'plt' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/subjects/total' '/' 'pap_vis_bim_rep_nob' '_' 'plt' ],vis_bim_plt)

%
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/subjects/total' '/' 'pap_aud_bim_rep_nob' '_' 'tbl' ],aud_bim_plt)
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/subjects/total' '/' 'pap_vis_bim_rep_nob' '_' 'tbl' ],vis_bim_plt)

tcfg = [];

tcfg.sbj_clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
tcfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ//epoch_data';
tcfg.sbj_nme = 'total';
tcfg.typ     = {'lfp'};
tcfg.ele_typ = {'ecog'};
tcfg.cmb_nme = {'pap_vis_bim_rep_nob'};
tcfg.loc_typ = {'split'};
mmil_sig_cmb_tbl(tcfg)

tcfg.cmb_nme = {'pap_aud_bim_rep_nob'};
mmil_sig_cmb_tbl(tcfg)

mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/table/lhs/total_pap_vis_bim_rep_nob_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/total/total_pap_vis_bim_rep_nob_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/table/rhs/total_pap_vis_bim_rep_nob_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/total/total_pap_vis_bim_rep_nob_rhs_table_plot')
mmil_chk_dir('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/total')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/table/lhs/total_pap_aud_bim_rep_nob_lhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/total/total_pap_aud_bim_rep_nob_lhs_table_plot')
    copyfile('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/table/rhs/total_pap_aud_bim_rep_nob_rhs_table_plot','/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/total/total_pap_aud_bim_rep_nob_rhs_table_plot')

% Timing - Repetition lfp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aud_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_rep_nob/subjects/total' '/' 'pap_aud_rep_nob' '_' 'plt' ]);
vis_plt = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_rep_nob/subjects/total' '/' 'pap_vis_rep_nob' '_' 'plt' ]);

aud_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_rep_nob/subjects/total' '/' 'pap_aud_rep_nob' '_' 'fst' ]);
vis_fst = mmil_readtext([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_rep_nob/subjects/total' '/' 'pap_vis_rep_nob' '_' 'fst' ]);

%
aud_bim_fst(:,1:2) = aud_fst(:,1:2);
aud_bim_fst(1,3)   = {'UniModal-Auditory'};
aud_bim_fst(1,4)   = {'BiModal-Auditory'};
for iR = 2:size(aud_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = aud_fst{iR,3};
    elseif aud_plt{iR,3}==1 && vis_plt{iR,3}==0
        aud_bim_fst{iR,3} = aud_fst{iR,3};
        aud_bim_fst{iR,4} = NaN;
    else
        aud_bim_fst{iR,3} = NaN;
        aud_bim_fst{iR,4} = NaN;
    end
end

%
vis_bim_fst(:,1:2) = aud_fst(:,1:2);
vis_bim_fst(1,3)   = {'UniModal-Visual'};
vis_bim_fst(1,4)   = {'BiModal-Visual'};
for iR = 2:size(vis_bim_fst,1)
    if aud_plt{iR,3}==1 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = vis_fst{iR,3};
    elseif aud_plt{iR,3}==0 && vis_plt{iR,3}==1
        vis_bim_fst{iR,3} = vis_fst{iR,3};
        vis_bim_fst{iR,4} = NaN;
    else
        vis_bim_fst{iR,3} = NaN;
        vis_bim_fst{iR,4} = NaN;
    end
    if (aud_plt{iR,3}==0 && vis_plt{iR,3}==1) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = vis_fst{iR,3};
    elseif (aud_plt{iR,3}==1 && vis_plt{iR,3}==0) && ~(aud_plt{iR,3}==1 && vis_plt{iR,3}==1)
        vis_bim_fst{iR,5} = aud_fst{iR,3};
    else
        vis_bim_fst{iR,5} = NaN;
    end
end

mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_aud_bim_rep_nob/subjects/total' '/' 'pap_aud_bim_rep_nob' '_' 'fst' ],aud_bim_fst)
mmil_chk_dir([ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/subjects/total']);
cell2csv(    [ '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/lfp/ecog/split/pap_vis_bim_rep_nob/subjects/total' '/' 'pap_vis_bim_rep_nob' '_' 'fst' ],vis_bim_fst)

























