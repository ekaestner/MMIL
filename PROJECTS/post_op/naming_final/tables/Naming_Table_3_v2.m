out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/tables/Table3';

cog_cut = -1.5;

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
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
                 'mean/std,1,tle_post_3T_ATLonly_left,bnt_raw_scr' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,bnt_raw_scr' ...
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
                 'mean/std,1,tle_post_3T_ATLonly_left,ant_mem_raw_scr' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,ant_mem_raw_scr' ...
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
cell2csv( [ out_dir '/' 'Table3.csv' ], out_tbl)