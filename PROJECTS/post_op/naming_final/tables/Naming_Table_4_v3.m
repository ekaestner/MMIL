out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/tables/Table4';

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

%% Setup
dta_loc = { [ prj_dir '/' prj_nme '/' 'SpecificCor/Clinical/Correlation/TLE_Controls_pre_cln'  '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/MRI/subcort_vol_ICV_cor/Raw/tle_controls_pre_3T_allSurg_all' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/fiber_FA/Raw/tle_controls_pre_3T_allSurg_all' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/wmparc_FA_wm/Raw/tle_controls_pre_3T_allSurg_all' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/Cognitive/Correlation/TLE_pre_post' '/' ] };

tbl_nme = { 'Educ' 1 ; ...
            'AgeOfSeizureOnset' 1 ; ...
            'AgeAtImaging' 1 ; ...
            'NumAEDs' 1 ; ...
            'xLeft.Hippocampus' 2 ; ...
            'xRight.Hippocampus' 2 ; ...
            'xL.ILF' 3 ; ...
            'xR.ILF' 3 ; ...
            'xL.IFO' 3 ; ...
            'xR.IFO' 3 ; ...
            'xlh.fusiform' 4 ; ...
            'xrh.fusiform' 4 };
        
col_lbl = { '' 'Pre-operative BNT' 'Pre-operative ANT' };
row_lbl = { 'Education' 1 ;...
            'AgeofSeizureOnset' 1 ;...
            'Age' 1 ; ...
            'AED #' 1 ;...
            'L-Hippocampus' 2 ;...
            'R-Hippocampus' 2 ;...
            'L-ILF' 2 ;...
            'R-ILF' 2 ;...
            'L-IFOF' 2 ;...
            'R-IFOF' 2 ;...
            'L-Fusiform' 2 ;...
            'R-Fusiform' 2 };

%% Out Table - OUTPUT (r-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%% Supp TABLE 1 - Significance (p-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

%% Supp TABLE 2 - Stat Out (summary)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_stats.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_stt = ejk_create_table( fcfg );

%% Supp TABLE 3 - FDR (summary)
fdr_grp_hld = cell2mat(row_lbl(:,2));
fdr_typ_hld = unique(fdr_grp_hld);

for iF = 1:numel(fdr_typ_hld)
    fdr_grp_mbm = fdr_grp_hld==fdr_typ_hld(iF);
    fdr_hld = cell2mat(tbl_pvl(fdr_grp_mbm,:));
    fdr_use{iF} = FDR(fdr_hld(:),.05);
    if isempty(fdr_use{iF}); fdr_use{iF} = 0; end
end

tbl_fdr=tbl_pvl;
for iR = 1:size(tbl_out,1)
    fdr_ind = fdr_typ_hld==row_lbl{iR,2};
    for iC = 1:size(tbl_out,2)
       if tbl_fdr{iR,iC}>fdr_use{fdr_ind}
           tbl_fdr{iR,iC}=NaN;
       end       
    end
end

%% Format table
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end
        
        if ~isnan(tbl_fdr{iR,iC})
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '#'];
        end
        
    end
end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];
out_stt = [ col_lbl ; row_lbl(:,1) tbl_stt ];

%% Save Table
cell2csv( [ out_dir '/' 'Table4.csv' ], out_tbl)
cell2csv( [ out_dir '/' 'Table4_stats.csv' ], out_stt)