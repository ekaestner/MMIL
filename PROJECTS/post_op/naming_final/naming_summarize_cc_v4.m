dta_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor']; % CrossCorrelation_3T_ATLonly % CrossCorrelation_3T_ATLonly_noQC

run_grp = { 'tle_controls_pre_3T_allSurg_all'   ...
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_left' ...
            'tle_post_3T_ATLonly_right' };

ejk_chk_dir([ dta_dir '/' 'Summary' '/' 'apriori']);

%% Include
mes_dir{1} = 'Clinical'; mes_typ_dir{1} = {'Correlation/TLE_Controls_pre_cln' 'Correlation/LTLE_Controls_pre_cln' 'Correlation/RTLE_Controls_pre_cln' 'Correlation/LTLE_post_cln' 'Correlation/RTLE_post_cln'}; mes_sub_dir{1} = '';
roi_int{1} = { 'AgeAtSurgery' 'Educ' 'AgeOfSeizureOnset' 'NumAEDs' 'SeizureFreq' };

mes_dir{2} = 'Cognitive'; mes_typ_dir{2} = {'Correlation/TLE_Controls_pre_pre' 'Correlation/LTLE_Controls_pre_pre' 'Correlation/RTLE_Controls_pre_pre' 'Correlation/LTLE_pre_post' 'Correlation/RTLE_pre_post'}; mes_sub_dir{2} = '';
roi_int{2} = { 'bnt.raw.scr' 'ant.mem.raw.scr' };

mes_dir{3} = 'MRI'; mes_typ_dir{3} = 'subcort_vol_ICV_cor'; mes_sub_dir{3} = 'Raw';
roi_int{3} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

mes_dir{4} = 'DTI'; mes_typ_dir{4} = 'fiber_FA'; mes_sub_dir{4} = 'Raw';
roi_int{4} = { 'xL.ILF' 'xL.IFO' 'xR.ILF' 'xR.IFO' };

mes_dir{5} = 'DTI'; mes_typ_dir{5} = 'wmparc_FA_wm'; mes_sub_dir{5} = 'Raw';
roi_int{5} = { 'xlh.fusiform' ...
               'xrh.fusiform' };         

grp_mes_cut = { [1 2] [3 4 5] };
           
clear pvl_hld_fdr
           
% Concatenate
for iG = 1:numel(run_grp)
    
    for iM = 1:numel(mes_dir)
        
        if ~iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_pvalues.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_n.csv' ]);
        elseif iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_pvalues.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_n.csv' ]);
        end
        
         % Get p-values
        for iC = 1:size(rvl_hld,2)-1
            cog_col_hld = iC+1;
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                pvl_hld_fdr{iG,iM}(iR,iC) =  pvl_hld{neu_col, cog_col_hld};
            end
        end
        
        for iC = 1:size(rvl_hld,2)-1
            
            cog_col_hld = iC+1;
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                         
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col_hld},2));
                if strcmpi(rvl_str(1),'-')
                    rvl_str = rvl_str([1 3:end]);
                else
                    rvl_str = rvl_str(2:end);
                end
                num_str = num2str(num_hld{neu_col, cog_col_hld});
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col_hld},2));
                pvl_str = pvl_str(2:end);
                tot_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                if ~iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                elseif iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = '';
                end
                tot_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                tot_out_hld{iG}{iM}{row_cnt, iC+4} = [ 'r(' num_str ') = ' rvl_str '; p = ' pvl_str];
                
                row_cnt = row_cnt + 1;
            end
        end
        
    end
        
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'include_total_raw' '_' run_grp{iG} '.csv' ], [ {''} {''} {''} {''}  rvl_hld(1,2:end) ; tot_out_hld{iG} ])
    
    clear tot_out_hld
    
end

%% FDR Correction
% Put together for FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(run_grp)
    
    for iTM = 1:numel(grp_mes_cut)
        hld_pvl              = cat(1,pvl_hld_fdr{ iG, grp_mes_cut{iTM} });
        fdr_cut_off_grp{iTM} = FDR( hld_pvl(:), .05);        
    end

    for iM = 1:numel(mes_dir)        
        % Individual Measures
        fdr_cut_off_ind = FDR( pvl_hld_fdr{iG,iM}, .05 );
        if ~isempty(fdr_cut_off_ind)
            fdr_cut_ind{iG,iM} = pvl_hld_fdr{iG,iM}<=fdr_cut_off_ind;
        else
            fdr_cut_ind{iG,iM} = zeros(size(pvl_hld_fdr{iG,iM}));
        end
        % Grouped measures
        if ~isempty(fdr_cut_off_grp{~cellfun(@isempty,cellfun(@(x) intersect(x,iM),grp_mes_cut,'UniformOutput',false))})
            fdr_cut_grp{iG,iM} = pvl_hld_fdr{iG,iM}<=fdr_cut_off_grp{~cellfun(@isempty,cellfun(@(x) intersect(x,iM),grp_mes_cut,'UniformOutput',false))};
        else
            fdr_cut_grp{iG,iM} = zeros(size(pvl_hld_fdr{iG,iM}));
        end
    end    
    
end

% Concatenate again %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(run_grp)
    
    for iM = 1:numel(mes_dir)
        
        if ~iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_pvalues.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_n.csv' ]);
        elseif iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_pvalues.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_n.csv' ]);
        end
                  
        % Put together
        for iC = 1:size(rvl_hld,2)-1
            
            cog_col_hld = iC+1;
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                   
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col_hld},2));
                if strcmpi(rvl_str(1),'-')
                    rvl_str = rvl_str([1 3:end]);
                else
                    rvl_str = rvl_str(2:end);
                end
                num_str = num2str(num_hld{neu_col, cog_col_hld});                
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col_hld},2));
                pvl_str = pvl_str(2:end);                
                tot_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                if ~iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                elseif iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = '';
                end
                tot_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                if fdr_cut_grp{iG,iM}(row_cnt,iC)
                tot_out_hld{iG}{iM}{row_cnt, iC+4} = [ 'r(' num_str ') = ' rvl_str '; p = ' pvl_str];
                else
                tot_out_hld{iG}{iM}{row_cnt, iC+4} = [ 'r(' num_str ') = ' rvl_str '; p = ' 'NaN'];
                end
                row_cnt = row_cnt + 1;
            end
        end
        
    end
        
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'include_total_raw' '_' run_grp{iG} '_FDR' '.csv' ], [ {''} {''} {''} {''}  rvl_hld(1,2:end) ; tot_out_hld{iG} ])
    
    clear tot_out_hld
    
end

%% Shuffle Correction
% Concatenate
for iG = 1:numel(run_grp)
    
    for iM = 1:numel(mes_dir)
        
        if ~iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_pvalues_shuffle.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_n.csv' ]);
        elseif iscell(mes_typ_dir{iM})
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_rvalues.csv' ]);
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_pvalues_shuffle.csv' ]);
            num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM}{iG} '/' 'cross_correlation_n.csv' ]);
        end
        
         % Get p-values
        for iC = 1:size(rvl_hld,2)-1
            cog_col_hld = iC+1;
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                pvl_hld_fdr{iG,iM}(iR,iC) =  pvl_hld{neu_col, cog_col_hld};
            end
        end
        
        for iC = 1:size(rvl_hld,2)-1
            
            cog_col_hld = iC+1;
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                         
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col_hld},2));
                if strcmpi(rvl_str(1),'-')
                    rvl_str = rvl_str([1 3:end]);
                else
                    rvl_str = rvl_str(2:end);
                end
                num_str = num2str(num_hld{neu_col, cog_col_hld});
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col_hld},2));
                if ~strcmpi(pvl_str,'NaN')
                    pvl_str = pvl_str(2:end);
                end
                tot_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                if ~iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                elseif iscell(mes_typ_dir{iM})
                    tot_out_hld{iG}{iM}{row_cnt,2} = '';
                end
                tot_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                tot_out_hld{iG}{iM}{row_cnt, iC+4} = [ 'r(' num_str ') = ' rvl_str '; p = ' pvl_str];
                
                row_cnt = row_cnt + 1;
            end
        end
        
    end
        
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'include_total_raw' '_' run_grp{iG} '_' 'shuffled' '.csv' ], [ {''} {''} {''} {''}  rvl_hld(1,2:end) ; tot_out_hld{iG} ])
    
    clear tot_out_hld
    
end

%% FDR Play
grp_run = 4; % [1 4 5];
mes_grp = [4 5 6];

%
clc

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(1) }(:), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(1) } pvl_hld_fdr{ grp_run, mes_grp(1) }<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(1) } zeros(size(pvl_hld_fdr{ grp_run, mes_grp(1) })) ]
end

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(2) }(:), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(2) } pvl_hld_fdr{ grp_run, mes_grp(2) }<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(2) } zeros(size(pvl_hld_fdr{ grp_run, mes_grp(2) })) ]
end

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(3) }(:), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(3) } pvl_hld_fdr{ grp_run, mes_grp(3) }<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(3) } zeros(size(pvl_hld_fdr{ grp_run, mes_grp(3) })) ]
end

tot_pvl_hld = cat(1, pvl_hld_fdr{ grp_run, mes_grp });
fdr_cut_off = FDR(tot_pvl_hld, .05);
if ~isempty(fdr_cut_off)
    [tot_pvl_hld tot_pvl_hld<=fdr_cut_off]
else
  [tot_pvl_hld zeros(size(tot_pvl_hld)) ]  
end

%
clc

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(1) }(:,1), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(1) }(:,1) pvl_hld_fdr{ grp_run, mes_grp(1) }(:,1)<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(1) }(:,1) zeros(size(pvl_hld_fdr{ grp_run, mes_grp(1) }(:,1))) ]
end

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(2) }(:,1), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(2) }(:,1) pvl_hld_fdr{ grp_run, mes_grp(2) }(:,1)<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(2) }(:,1) zeros(size(pvl_hld_fdr{ grp_run, mes_grp(2) }(:,1))) ]
end

fdr_cut_off = FDR(pvl_hld_fdr{ grp_run, mes_grp(3) }(:,1), .05)
if ~isempty(fdr_cut_off);
[ pvl_hld_fdr{ grp_run, mes_grp(3) }(:,1) pvl_hld_fdr{ grp_run, mes_grp(3) }(:,1)<=fdr_cut_off]
else
[ pvl_hld_fdr{ grp_run, mes_grp(3) }(:,1) zeros(size(pvl_hld_fdr{ grp_run, mes_grp(3) }(:,1))) ]
end

tot_pvl_hld = cat(1, pvl_hld_fdr{ grp_run, mes_grp });
tot_pvl_hld = tot_pvl_hld(:,1);
fdr_cut_off = FDR(tot_pvl_hld, .05);
if ~isempty(fdr_cut_off)
    [tot_pvl_hld tot_pvl_hld<=fdr_cut_off]
else
  [tot_pvl_hld zeros(size(tot_pvl_hld)) ]  
end







