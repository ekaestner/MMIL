out_dir =  [ prj_dir '/' prj_nme '/' 'Tables' ];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

clear out_tbl

%% Out Table - OUTPUT (r-value)
clear dta_inp
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Clinical/Correlation/RTLE_post_cln/cross_correlation_rvalues.csv']);
dta_inp{2} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/MRI/subcort_vol_ICV_cor/Raw/tle_post_3T_ATLonly_right/cross_correlation_rvalues.csv']);
dta_inp{3} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_right/cross_correlation_rvalues.csv']);
dta_inp{4} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_right/cross_correlation_rvalues.csv']);
dta_inp{5} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Cognitive/Correlation/RTLE_pre_post/cross_correlation_rvalues.csv']);

% %%%
fcfg = [];

fcfg.tbl(1,:) = { 'copy,5,bnt.raw.scr,xbnt.raw.scr.pst' ... 
                  'copy,5,ant.mem.raw.scr,xant.mem.raw.scr.pst'  };

fcfg.tbl(2,:) = { 'copy,1,Educ,bnt.raw.scr.pst' ... 
                  'copy,1,Educ,ant.mem.raw.scr.pst' };
              
fcfg.tbl(3,:) = { 'copy,1,AgeOfSeizureOnset,bnt.raw.scr.pst' ... 
                  'copy,1,AgeOfSeizureOnset,ant.mem.raw.scr.pst' };

fcfg.tbl(4,:) = { 'copy,1,NumAEDs,bnt.raw.scr.pst' ... 
                  'copy,1,NumAEDs,ant.mem.raw.scr.pst' };
                            
fcfg.tbl(5,:) = {'copy,2,xLeft.Hippocampus,bnt.raw.scr.pst' ... 
                 'copy,2,xLeft.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(6,:) = {'copy,2,xRight.Hippocampus,bnt.raw.scr.pst' ...
                 'copy,2,xRight.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(7,:) = {'copy,3,xL.ILF,bnt.raw.scr.pst' ...
                 'copy,3,xL.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(8,:) = {'copy,3,xR.ILF,bnt.raw.scr.pst' ...
                 'copy,3,xR.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(9,:) = {'copy,3,xL.IFO,bnt.raw.scr.pst' ...
                 'copy,3,xL.IFO,ant.mem.raw.scr.pst' };
             
fcfg.tbl(10,:) = {'copy,3,xR.IFO,bnt.raw.scr.pst' ...
                 'copy,3,xR.IFO,ant.mem.raw.scr.pst' };

fcfg.tbl(11,:) = {'copy,4,xlh.fusiform,bnt.raw.scr.pst' ...
                 'copy,4,xlh.fusiform,ant.mem.raw.scr.pst' };
             
fcfg.tbl(12,:) ={'copy,4,xrh.fusiform,bnt.raw.scr.pst' ...
                 'copy,4,xrh.fusiform,ant.mem.raw.scr.pst' };

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%% Supp TABLE 1 - Significance (p-value)
clear dta_inp
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Clinical/Correlation/RTLE_post_cln/cross_correlation_pvalues.csv']);
dta_inp{2} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/MRI/subcort_vol_ICV_cor/Raw/tle_post_3T_ATLonly_right/cross_correlation_pvalues.csv']);
dta_inp{3} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_right/cross_correlation_pvalues.csv']);
dta_inp{4} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_right/cross_correlation_pvalues.csv']);
dta_inp{5} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Cognitive/Correlation/RTLE_pre_post/cross_correlation_pvalues.csv']);

% %%%
fcfg.tbl(1,:) = { 'copy,5,bnt.raw.scr,xbnt.raw.scr.pst' ... 
                  'copy,5,ant.mem.raw.scr,xant.mem.raw.scr.pst'  };

fcfg.tbl(2,:) = { 'copy,1,Educ,bnt.raw.scr.pst' ... 
                  'copy,1,Educ,ant.mem.raw.scr.pst' };
              
fcfg.tbl(3,:) = { 'copy,1,AgeOfSeizureOnset,bnt.raw.scr.pst' ... 
                  'copy,1,AgeOfSeizureOnset,ant.mem.raw.scr.pst' };

fcfg.tbl(4,:) = { 'copy,1,NumAEDs,bnt.raw.scr.pst' ... 
                  'copy,1,NumAEDs,ant.mem.raw.scr.pst' };
                            
fcfg.tbl(5,:) = {'copy,2,xLeft.Hippocampus,bnt.raw.scr.pst' ... 
                 'copy,2,xLeft.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(6,:) = {'copy,2,xRight.Hippocampus,bnt.raw.scr.pst' ...
                 'copy,2,xRight.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(7,:) = {'copy,3,xL.ILF,bnt.raw.scr.pst' ...
                 'copy,3,xL.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(8,:) = {'copy,3,xR.ILF,bnt.raw.scr.pst' ...
                 'copy,3,xR.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(9,:) = {'copy,3,xL.IFO,bnt.raw.scr.pst' ...
                 'copy,3,xL.IFO,ant.mem.raw.scr.pst' };
             
fcfg.tbl(10,:) = {'copy,3,xR.IFO,bnt.raw.scr.pst' ...
                 'copy,3,xR.IFO,ant.mem.raw.scr.pst' };

fcfg.tbl(11,:) = {'copy,4,xlh.fusiform,bnt.raw.scr.pst' ...
                 'copy,4,xlh.fusiform,ant.mem.raw.scr.pst' };
             
fcfg.tbl(12,:) ={'copy,4,xrh.fusiform,bnt.raw.scr.pst' ...
                 'copy,4,xrh.fusiform,ant.mem.raw.scr.pst' };

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

%% Supp TABLE 2 - FDR (summary)
clear dta_inp
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Summary/apriori/include_total_raw_tle_post_3T_ATLonly_right_FDR.csv']);
    dta_inp{1} = dta_inp{1}(:,4:end);

    % %%%
fcfg = [];

fcfg.tbl(1,:) = { 'copy,1,bnt.raw.scr,bnt.raw.scr.pst' ... 
                  'copy,1,ant.mem.raw.scr,ant.mem.raw.scr.pst'  };

fcfg.tbl(2,:) = { 'copy,1,Educ,bnt.raw.scr.pst' ... 
                  'copy,1,Educ,ant.mem.raw.scr.pst' };
              
fcfg.tbl(3,:) = { 'copy,1,AgeOfSeizureOnset,bnt.raw.scr.pst' ... 
                  'copy,1,AgeOfSeizureOnset,ant.mem.raw.scr.pst' };

fcfg.tbl(4,:) = { 'copy,1,NumAEDs,bnt.raw.scr.pst' ... 
                  'copy,1,NumAEDs,ant.mem.raw.scr.pst' };
                            
fcfg.tbl(5,:) = {'copy,1,xLeft.Hippocampus,bnt.raw.scr.pst' ... 
                 'copy,1,xLeft.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(6,:) = {'copy,1,xRight.Hippocampus,bnt.raw.scr.pst' ...
                 'copy,1,xRight.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(7,:) = {'copy,1,xL.ILF,bnt.raw.scr.pst' ...
                 'copy,1,xL.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(8,:) = {'copy,1,xR.ILF,bnt.raw.scr.pst' ...
                 'copy,1,xR.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(9,:) = {'copy,1,xL.IFO,bnt.raw.scr.pst' ...
                 'copy,1,xL.IFO,ant.mem.raw.scr.pst' };
             
fcfg.tbl(10,:) = {'copy,1,xR.IFO,bnt.raw.scr.pst' ...
                 'copy,1,xR.IFO,ant.mem.raw.scr.pst' };

fcfg.tbl(11,:) = {'copy,1,xlh.fusiform,bnt.raw.scr.pst' ...
                 'copy,1,xlh.fusiform,ant.mem.raw.scr.pst' };
             
fcfg.tbl(12,:) ={'copy,1,xrh.fusiform,bnt.raw.scr.pst' ...
                 'copy,1,xrh.fusiform,ant.mem.raw.scr.pst' };

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_fdr = ejk_create_table( fcfg );

%% Supp TABLE 2 - Out stat (summary)
clear dta_inp
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' '/SpecificCor/Summary/apriori/include_total_raw_tle_post_3T_ATLonly_right.csv']);
    dta_inp{1} = dta_inp{1}(:,4:end);

    % %%%
fcfg = [];

fcfg.tbl(1,:) = { 'copy,1,bnt.raw.scr,bnt.raw.scr.pst' ... 
                  'copy,1,ant.mem.raw.scr,ant.mem.raw.scr.pst'  };

fcfg.tbl(2,:) = { 'copy,1,Educ,bnt.raw.scr.pst' ... 
                  'copy,1,Educ,ant.mem.raw.scr.pst' };
              
fcfg.tbl(3,:) = { 'copy,1,AgeOfSeizureOnset,bnt.raw.scr.pst' ... 
                  'copy,1,AgeOfSeizureOnset,ant.mem.raw.scr.pst' };

fcfg.tbl(4,:) = { 'copy,1,NumAEDs,bnt.raw.scr.pst' ... 
                  'copy,1,NumAEDs,ant.mem.raw.scr.pst' };
                            
fcfg.tbl(5,:) = {'copy,1,xLeft.Hippocampus,bnt.raw.scr.pst' ... 
                 'copy,1,xLeft.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(6,:) = {'copy,1,xRight.Hippocampus,bnt.raw.scr.pst' ...
                 'copy,1,xRight.Hippocampus,ant.mem.raw.scr.pst' };
             
fcfg.tbl(7,:) = {'copy,1,xL.ILF,bnt.raw.scr.pst' ...
                 'copy,1,xL.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(8,:) = {'copy,1,xR.ILF,bnt.raw.scr.pst' ...
                 'copy,1,xR.ILF,ant.mem.raw.scr.pst' };
             
fcfg.tbl(9,:) = {'copy,1,xL.IFO,bnt.raw.scr.pst' ...
                 'copy,1,xL.IFO,ant.mem.raw.scr.pst' };
             
fcfg.tbl(10,:) = {'copy,1,xR.IFO,bnt.raw.scr.pst' ...
                 'copy,1,xR.IFO,ant.mem.raw.scr.pst' };

fcfg.tbl(11,:) = {'copy,1,xlh.fusiform,bnt.raw.scr.pst' ...
                 'copy,1,xlh.fusiform,ant.mem.raw.scr.pst' };
             
fcfg.tbl(12,:) ={'copy,1,xrh.fusiform,bnt.raw.scr.pst' ...
                 'copy,1,xrh.fusiform,ant.mem.raw.scr.pst' };

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_stt = ejk_create_table( fcfg );

%% Format table
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end
        
        if isempty(strfind(tbl_fdr{iR,iC},'p = NaN'))
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '#'];
        end
        
    end
end

col_lbl = { '' 'Pre-operative BNT' 'Pre-operative ANT' };
row_lbl = { 'Pre-operative Score' ...
            'Education' ...
            'Age of Seizure Onset' ...
            '# AEDs' ...
            'L-Hippocampus' ...
            'R-Hippocampus' ...
            'L-ILF' ...
            'R-ILF' ...
            'L-IFOF' ...
            'R-IFOF' ...
            'L-Fusiform' ...
            'R-Fusiform' }';

out_tbl = [ col_lbl ; row_lbl out_tbl ];
out_stt = [ col_lbl ; row_lbl tbl_stt ];

%% Save Table
cell2csv( [ out_dir '/' 'Table5_2.csv' ], out_tbl)
cell2csv( [ out_dir '/' 'Table5_2_stats.csv' ], out_stt)