dta_dir = [ prj_dir '/' prj_nme '/' 'TotalCor']; % CrossCorrelation_3T_ATLonly % CrossCorrelation_3T_ATLonly_noQC

run_grp = fieldnames(grp);

ejk_chk_dir([ dta_dir '/' 'Summary' '/' 'apriori']);

cog_scr_nme = { { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr'     } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr'     } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr' 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' 'cat.flu.nor.scr.pst'  } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr' 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' 'cat.flu.nor.scr.pst'  } };

%% Raw
mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA'; mes_sub_dir{1} = 'Raw';
roi_int{1} = { 'xL.ILF' 'xL.tSLF' 'xL.Unc' 'xR.ILF' 'xR.tSLF' 'xR.Unc' 'xL.IFO' 'xR.IFO' 'xL.CgH' 'xR.CgH' };

mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm'; mes_sub_dir{2} = 'Raw';
roi_int{2} = { 'xlh.supramarginal' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.fusiform' 'xlh.inferiortemporal' 'xlh.precuneus' ...
               'xrh.supramarginal' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.fusiform' 'xrh.inferiortemporal' 'xrh.precuneus' };

mes_dir{3} = 'MRI'; mes_typ_dir{3} = 'subcort_vol_ICV_cor'; mes_sub_dir{3} = 'Raw';
roi_int{3} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };
           
mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'cort_thick_ctx'; mes_sub_dir{4} = 'Raw';
roi_int{4} = { 'xlh.supramarginal' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.fusiform' 'xlh.inferiortemporal' 'xlh.precuneus' ...
               'xrh.supramarginal' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.fusiform' 'xrh.inferiortemporal' 'xrh.precuneus' };

mes_dir{5} = 'rsfMRI'; mes_typ_dir{5} = 'var_ctx'; mes_sub_dir{5} = 'Raw';
roi_int{5} = {'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus'...
              'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_vol'; mes_sub_dir{6} = 'Raw';
roi_int{6} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

% Concatenate
for iG = 1:numel(run_grp)
    
    for iM = 1:numel(mes_dir)
                
        rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_n.csv' ]);
        
        for iC = 1:size(rvl_hld,2)-1
            
            cog_col_hld = iC+1;
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                
                rvl_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                
                rvl_out_hld{iG}{iM}{row_cnt, iC+4} = num2str(roundsd(rvl_hld{neu_col, cog_col_hld},2));
                if rvl_hld{neu_col, cog_col_hld} > 0
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}(2:end);
                else
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}([1 3:end]);
                end
                if pvl_hld{neu_col, cog_col_hld}<=.05;                                          rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' *']; end
                if pvl_hld{neu_col, cog_col_hld}<=.10 && pvl_hld{neu_col, cog_col_hld}>.05; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' #']; end
                if pvl_hld{neu_col, cog_col_hld}>.10; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' ']; end
                
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col_hld},2));
                rvl_str = rvl_str(2:end);
                num_str = num2str(num_hld{neu_col, cog_col_hld});
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col_hld},2));
                pvl_str = pvl_str(2:end);
                tot_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                tot_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                tot_out_hld{iG}{iM}{row_cnt, iC+4} = [ 'r(' num_str ') = ' rvl_str '; p = ' pvl_str];
                
                row_cnt = row_cnt + 1;
            end
        end
        
    end
    
    rvl_out_hld{iG} = cat(1, rvl_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'rvalues_raw' '_' run_grp{iG} '.csv' ], [ {''} {''} {''} {''} rvl_hld(1,2:end) ; rvl_out_hld{iG} ])
    
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'total_raw' '_' run_grp{iG} '.csv' ], [ {''} {''} {''} {''}  rvl_hld(1,2:end) ; tot_out_hld{iG} ])
    
    clear rvl_out_hld tot_out_hld
    
end

%% Brain ROI Summarize
mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA';            mes_sub_dir{1} = 'Raw'; roi_use{1} = [0];
mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm';        mes_sub_dir{2} = 'Raw'; roi_use{2} = [2];
mes_dir{3} = 'MRI'; mes_typ_dir{3} = 'cort_thick_ctx';      mes_sub_dir{3} = 'Raw'; roi_use{3} = [2];
mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'subcort_vol_ICV_cor'; mes_sub_dir{4} = 'Raw'; roi_use{4} = [0];

pvl_cut = .05;

run_grp = run_grp([ 7 8 9 10 ]);

for iM = 1:numel(mes_dir)
    if roi_use{iM} > 0
        for iG = 1:numel(run_grp)
            
            out_dir_mes = [dta_dir '/' 'Summary' '/' 'ROI_plots' '/' mes_dir{iM} '/'];
            ejk_chk_dir(out_dir_mes);
            
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_rvalues.csv' ]);
            rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
            
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' run_grp{iG} '/' 'cross_correlation_pvalues.csv' ]);
            pvl_hld(strcmpi(pvl_hld,'NA')) = {NaN};
            
            for iC = 2:size(rvl_hld,2)
                
                % lhs
                plt_tbl{1} = rvl_hld( string_find(rvl_hld(:,1),'lh.'), [1 iC]);
                pvl_tbl{1} = pvl_hld( string_find(pvl_hld(:,1),'lh.'), [1 iC]);
                plt_tbl{1}(cell2mat(pvl_tbl{1}(:,2))>pvl_cut,2) = {0};
                if strcmpi(plt_tbl{1}{1}(1),'x'); plt_tbl{1}(:,1) = cellfun( @(x) x(2:end),plt_tbl{1}(:,1),'uni',0); end
                % rhs
                plt_tbl{2} = rvl_hld( string_find(rvl_hld(:,1),'rh.'), [1 iC]);
                plt_tbl{2}(1,:) = [];
                pvl_tbl{2} = pvl_hld( string_find(pvl_hld(:,1),'rh.'), [1 iC]);
                pvl_tbl{2}(1,:) = [];
                plt_tbl{2}(cell2mat(pvl_tbl{2}(:,2))>pvl_cut,2) = {0};
                if strcmpi(plt_tbl{2}{1}(1),'x'); plt_tbl{2}(:,1) = cellfun( @(x) x(2:end),plt_tbl{2}(:,1),'uni',0); end
                
                %% Plot
                cfg = [];
                
                cfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc';
                cfg.fsr_nme = 'fsaverage';
                
                cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
                cfg.prc_nme = '.aparc.annot';
                
                cfg.col_map = { rgb('purplish blue') rgb('blue') rgb('grayish blue') rgb('light grey') rgb('reddish grey') rgb('red') rgb('orangish red') };
                cfg.low_rng_num = [ -0.20 0.20 ];
                cfg.hgh_rng_num = [ -0.50 0.50 ];
                
                cfg.roi_tbl_loc = [ {plt_tbl{1}(:,1)}           {plt_tbl{2}(:,1)}    ];
                cfg.roi_tbl = [ {cell2mat(plt_tbl{1}(:,2))} {cell2mat(plt_tbl{2}(:,2))}    ];
                
                cfg.sph = { 'lh' 'rh' };
                cfg.sph_vew = { 'lat' 'ven' 'med' };
                
                cfg.out_dir = [ dta_dir '/' 'Summary' '/' 'ROI_plots' '/' mes_dir{iM} ];
                cfg.out_nme = [ mes_typ_dir{iM} '_' mes_sub_dir{iM} '_' run_grp{iG} '_' mmil_spec_char(rvl_hld{1,iC},{'.'})  ];
                
                ejk_roi_display_plot(cfg);
                
            end
            
        end
    end
end