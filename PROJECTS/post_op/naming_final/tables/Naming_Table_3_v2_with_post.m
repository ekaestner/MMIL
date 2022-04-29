out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/tables/Table3';

cog_cut = -1.5;

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_with_post_raw_QC.csv']);
dta_inp{2} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Cognitive/ANOVA/Pre_Omnibus_anova/TLE.Controls.pre.anova/output_table.csv');
dta_inp{3} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Cognitive/ttest/Post_TLE_ttest/TLE.post.ttest/output_table.csv');
dta_inp{4} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Cognitive/Fisher/output_table.csv');

%%
dcl_dta = cog_dta;
for iR = 1:size(dcl_dta,1)
    for iC = 1:size(dcl_dta,2)
        if dcl_dta{iR,iC} <= cog_cut
            dcl_dta{iR,iC} = 'Impaired';
        elseif  dcl_dta{iR,iC} > cog_cut
            dcl_dta{iR,iC} = 'NotImpaired';
        elseif  isnan(dcl_dta{iR,iC})
            dcl_dta{iR,iC} = '';
        end
    end
end
dta_inp{5} = [ {''} cog_dta_col ; cog_dta_sbj dcl_dta];

%%
fcfg = [];

fcfg.tbl(1,:) = {'mean/std,1,controls_pre_3T_allSurg_all,bnt_raw_scr' ... 
                 'mean/std,1,tle_pre_3T_allSurg_left,bnt_raw_scr' ...
                 'mean/std,1,tle_pre_3T_allSurg_right,bnt_raw_scr' ...
                 'copy,2,bnt_raw_scr,report'};

fcfg.tbl(2,:) = {'empty' ... 
                 'mean/std,1,tle_post_3T_ATLonly_left,bnt_raw_scr_pst_raw' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,bnt_raw_scr_pst_raw' ...
                 'empty'};
             
fcfg.tbl(3,:) = {'empty' ... 
                 'mean/std,1,tle_post_3T_ATLonly_left,bnt_raw_scr_pst' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,bnt_raw_scr_pst' ...
                 'copy,3,bnt_raw_scr_pst,report'};
             
fcfg.tbl(4,:) = {'empty' ... 
                 'count,5,tle_post_3T_ATLonly_left,bnt_raw_scr_pst,Impaired/NotImpaired' ...
                 'count,5,tle_post_3T_ATLonly_right,bnt_raw_scr_pst,Impaired/NotImpaired' ...
                 'copy,4,bnt_raw_scr_pst,report'};

fcfg.tbl(5,:) = {'empty' ... 
                 'empty' ...
                 'empty' ...
                 'empty' };

             
fcfg.tbl(6,:) = {'mean/std,1,controls_pre_3T_allSurg_all,ant_mem_raw_scr' ... 
                 'mean/std,1,tle_pre_3T_allSurg_left,ant_mem_raw_scr' ...
                 'mean/std,1,tle_pre_3T_allSurg_right,ant_mem_raw_scr' ...
                 'copy,2,ant_mem_raw_scr,report'};         
             
fcfg.tbl(7,:) = {'empty' ... 
                 'mean/std,1,tle_post_3T_ATLonly_left,ant_mem_raw_scr_pst_raw' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,ant_mem_raw_scr_pst_raw' ...
                 'empty'};             
             
fcfg.tbl(8,:) = {'empty' ... 
                 'mean/std,1,tle_post_3T_ATLonly_left,ant_mem_raw_scr_pst' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,ant_mem_raw_scr_pst' ...
                 'copy,3,ant_mem_raw_scr_pst,report'};             
  
fcfg.tbl(9,:) = {'empty' ... 
                 'count,5,tle_post_3T_ATLonly_left,ant_mem_raw_scr_pst,Impaired/NotImpaired' ...
                 'count,5,tle_post_3T_ATLonly_right,ant_mem_raw_scr_pst,Impaired/NotImpaired' ...
                 'copy,4,ant_mem_raw_scr_pst,report'};
             
fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );
             
%%
col_lbl = { '' 'HC' 'L-TLE' 'R-TLE' 'Test' };
row_lbl = { 'BNT Pre-operative correct (total sample)' 'BNT Pre-operative correct (post-operative sample)' 'BNT Post-operative RCI (total sample)' 'BNT Decliners' ...
            '' ...
            'ANT Pre-operative correct (total sample)' 'ANT Pre-operative correct (post-operative sample)' 'ANT Post-operative RCI (total sample)' 'ANT Decliners'}';


out_tbl = [ col_lbl ; row_lbl tbl_out ];

%%
cell2csv( [ out_dir '/' 'Table3_with_post_raw.csv' ], out_tbl)

%% Test
% Overall
( ( 45.7 - 46.3 ) - (49.42-49.08) ) / 1.09
( ( 48.3 - 46.8 ) - (49.42-49.08) ) / 1.09

( ( 45.5 - 45.2 ) - (50.06 - 48.78) ) / 2.67
( ( 48.0 - 47.9 ) - (50.06 - 48.78) ) / 2.67

% Subset
% ANT
( ( 45.7 - 46.5 ) - (49.42-49.08) ) / 1.09
( ( 48.3 - 47.1 ) - (49.42-49.08) ) / 1.09
% BNT
( ( 45.5 - 44.8 ) - (50.06 - 48.78) ) / 2.67
( ( 49.8 - 48.3 ) - (50.06 - 48.78) ) / 2.67

% BNT - LTLE
clear dta_use
iR = 2; iC = 2;
cut_hld = regexp( cfg.tbl{iR,iC},',','split');
dta_use = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})));
dta_use(:,2) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'bnt_raw_scr')));
dta_use(:,3) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'bnt_raw_scr_pst')));
dta_use(:,4) = ( ( dta_use(:,1) - dta_use(:,2) ) - (50.06 - 48.78) ) / 2.67;
roundsd(nanmean(dta_use(~isnan(dta_use(:,3)),:),1),3)
roundsd(nanstd(dta_use(~isnan(dta_use(:,3)),:),1),3)
( ( 45.3 - 45.3 ) - (50.06 - 48.78) ) / 2.67

% BNT - RTLE
clear dta_use
iR = 2; iC = 3;
cut_hld = regexp( cfg.tbl{iR,iC},',','split');
dta_use = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})));
dta_use(:,2) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'bnt_raw_scr')));
dta_use(:,3) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'bnt_raw_scr_pst')));
dta_use(:,4) = ( ( dta_use(:,1) - dta_use(:,2) ) - (50.06 - 48.78) ) / 2.67;
roundsd(nanmean(dta_use(~isnan(dta_use(:,4)),:),1),3)
roundsd(nanstd(dta_use(~isnan(dta_use(:,4)),:),1),3)
( ( 49.8 - 48.3 ) - (50.06 - 48.78) ) / 2.67

% ANT - LTLE
clear dta_use
iR = 7; iC = 2;
cut_hld = regexp( cfg.tbl{iR,iC},',','split');
dta_use = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})));
dta_use(:,2) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'ant_mem_raw_scr')));
dta_use(:,3) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'ant_mem_raw_scr_pst')));
dta_use(:,4) = ( ( dta_use(:,1) - dta_use(:,2) )  - (49.42-49.08) ) / 1.09;
roundsd(nanmean(dta_use(~isnan(dta_use(:,3)),:),1),3)
roundsd(nanstd(dta_use(~isnan(dta_use(:,3)),:),1),3)
( ( 45.7 - 46.6 ) - (49.42-49.08) ) / 1.09

% ANT - RTLE
clear dta_use
iR = 7; iC = 3;
cut_hld = regexp( cfg.tbl{iR,iC},',','split');
dta_use = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})));
dta_use(:,2) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'ant_mem_raw_scr')));
dta_use(:,3) = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),'ant_mem_raw_scr_pst')));
dta_use(:,4) = ( ( dta_use(:,1) - dta_use(:,2) )  - (49.42-49.08) ) / 1.09;
roundsd(nanmean(dta_use(~isnan(dta_use(:,4)),:),1),3)
roundsd(nanstd(dta_use(~isnan(dta_use(:,4)),:),1),3)
( ( 48.3 - 47.9 ) - (49.42-49.08) ) / 1.09
