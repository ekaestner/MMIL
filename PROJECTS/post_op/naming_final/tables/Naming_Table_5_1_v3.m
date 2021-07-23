out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/tables/Table5';

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

%% Setup
dta_loc = { [ prj_dir '/' prj_nme '/' 'SpecificCor/Clinical/Correlation/LTLE_post_cln'  '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/MRI/subcort_vol_ICV_cor/Raw/tle_post_3T_ATLonly_left' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_left' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_left' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/Cognitive/Correlation/LTLE_pre_post' '/' ] };

tbl_nme = [ {'bnt.raw.scr'} {'ant.mem.raw.scr'} {5} ; ...
            repmat({'Educ'},1,2)  {1} ; ...
            repmat({'AgeOfSeizureOnset'},1,2) {1} ; ...
            repmat({'AgeAtImaging'},1,2) {1} ; ...
            repmat({'NumAEDs'},1,2) {1} ; ...
            repmat({'xLeft.Hippocampus'},1,2) {2} ; ...
            repmat({'xRight.Hippocampus'},1,2) {2} ; ...
            repmat({'xL.ILF'},1,2) {3} ; ...
            repmat({'xR.ILF'},1,2) {3} ; ...
            repmat({'xL.IFO'},1,2) {3} ; ...
            repmat({'xR.IFO'},1,2) {3} ; ...
            repmat({'xlh.fusiform'},1,2) {4} ; ...
            repmat({'xrh.fusiform'},1,2) {4} ];
        
col_lbl = { '' 'Post-operative BNT' 'Post-operative ANT' };
row_lbl = { 'Pre-operative' 1
            'Education' 1 ;...
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
      if iR==1
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',xbnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',xant.mem.raw.scr.pst'] };
    else
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',ant.mem.raw.scr.pst'] };
    end
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%% Supp TABLE 1 - Significance (p-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
     if iR==1
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',xbnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',xant.mem.raw.scr.pst'] };
    else
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',ant.mem.raw.scr.pst'] };
    end
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

%% Supp TABLE 2 - Stat Out (summary)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_stats.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
    if iR==1
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',xbnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',xant.mem.raw.scr.pst'] };
    else
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ...
                           ['copy,' num2str(tbl_nme{iR,3}) ',' num2str(tbl_nme{iR,2}) ',ant.mem.raw.scr.pst'] };
    end
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
cell2csv( [ out_dir '/' 'Table5_1.csv' ], out_tbl)
cell2csv( [ out_dir '/' 'Table5_1_stats.csv' ], out_stt)