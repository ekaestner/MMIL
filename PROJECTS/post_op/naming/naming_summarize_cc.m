dta_dir = [ prj_dir '/' prj_nme '/' 'CrossCorrelation_3T_ATLonly' '/' ]; % CrossCorrelation_3T_ATLonly % CrossCorrelation_3T_ATLonly_noQC

ejk_chk_dir([ dta_dir '/' 'Summary' '/' 'apriori']);

nme_typ = { 'left_pre' 'right_pre' 'left_post' 'right_post'};

cog_scr_nme = { { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr'     } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr'     } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr' 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' 'cat.flu.nor.scr.pst'  } ...
                { 'bnt.raw.scr' 'ant.mem.raw.scr' 'cat.flu.nor.scr' 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' 'cat.flu.nor.scr.pst'  } };

%% Raw
% mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'aseg_FA'; mes_sub_dir{1} = 'Raw';
% roi_int{1} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };
% 
% mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'aseg_MD'; mes_sub_dir{2} = 'Raw';
% roi_int{2} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA'; mes_sub_dir{1} = 'Raw';
roi_int{1} = { 'xL.ILF' 'xL.tSLF' 'xL.Unc' 'xR.ILF' 'xR.tSLF' 'xR.Unc' };

% mes_dir{4} = 'DTI'; mes_typ_dir{4} = 'fiber_MD'; mes_sub_dir{4} = 'Raw';
% roi_int{4} = { 'xL.ILF' 'xL.tSLF' 'xL.Unc' 'xR.ILF' 'xR.tSLF' 'xR.Unc' };

mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm'; mes_sub_dir{2} = 'Raw';
roi_int{2} = { 'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus' ...
    'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

% mes_dir{6} = 'DTI'; mes_typ_dir{6} = 'wmparc_MD_wm'; mes_sub_dir{6} = 'Raw';
% roi_int{6} = { 'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus' ...
%     'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

mes_dir{3} = 'fMRI'; mes_typ_dir{3} = 'alicia'; mes_sub_dir{3} = 'Raw';
roi_int{3} = { 'xlh.LateralTemporal02.cortical' 'xlh.InferiorFrontal.cortical' 'xlh.VentralTemporal.cortical'...
               'xrh.LateralTemporal02.cortical' 'xrh.InferiorFrontal.cortical' 'xrh.VentralTemporal.cortical' };
% 
% mes_dir{8} = 'fMRI'; mes_typ_dir{8} = 'destr'; mes_sub_dir{8} = 'Raw';
% roi_int{8} = { 'lh.Pole.temporal' 'lh.G.and.S.cingul.Mid.Post' 'lh.G.cuneus' 'lh.G.front.inf.Triangul' 'lh.G.front.inf.Opercular' 'lh.S.intrapariet.and.P.trans' ...
%     'rh.Pole.temporal' 'rh.G.and.S.cingul.Mid.Post' 'rh.G.cuneus' 'rh.G.front.inf.Triangul' 'rh.G.front.inf.Opercular' 'rh.S.intrapariet.and.P.trans'};

mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'cort_thick_ctx'; mes_sub_dir{4} = 'Raw';
roi_int{4} = { 'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus'...
    'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

mes_dir{5} = 'MRI'; mes_typ_dir{5} = 'subcort_vol_ICV_cor'; mes_sub_dir{5} = 'Raw';
roi_int{5} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_ctx'; mes_sub_dir{6} = 'Raw';
roi_int{6} = {'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus'...
    'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

mes_dir{7} = 'rsfMRI'; mes_typ_dir{7} = 'var_vol'; mes_sub_dir{7} = 'Raw';
roi_int{7} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

% Concatenate
for iG = 1:numel(nme_typ)
    
    for iM = 1:numel(mes_dir)
                
        rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_n.csv' ]);
        
        for iC = 1:numel(cog_scr_nme{iG})
            
            cog_col = find(ismember( rvl_hld(1,:), cog_scr_nme{iG}{iC}));
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                
                rvl_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                
                rvl_out_hld{iG}{iM}{row_cnt, iC+4} = num2str(roundsd(rvl_hld{neu_col, cog_col},2));
                if rvl_hld{neu_col, cog_col} > 0
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}(2:end);
                else
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}([1 3:end]);
                end
                if pvl_hld{neu_col, cog_col}<=.05;                                          rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' *']; end
                if pvl_hld{neu_col, cog_col}<=.10 && pvl_hld{neu_col, cog_col}>.05; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' #']; end
                if pvl_hld{neu_col, cog_col}>.10; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' ']; end
                
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col},2));
                rvl_str = rvl_str(2:end);
                num_str = num2str(num_hld{neu_col, cog_col});
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col},2));
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
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'rvalues_raw' '_' nme_typ{iG} '.csv' ], [ {''} {''} {''} {''}  cog_scr_nme{iG} ; rvl_out_hld{iG} ])
    
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'total_raw' '_' nme_typ{iG} '.csv' ], [ {''} {''} {''} {''}  cog_scr_nme{iG} ; tot_out_hld{iG} ])
    
    clear rvl_out_hld tot_out_hld
    
end

%% LI
% mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'aseg_FA'; mes_sub_dir{1} = 'Raw';
% roi_int{1} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };
% 
% mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'aseg_MD'; mes_sub_dir{2} = 'Raw';
% roi_int{2} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };

mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA'; mes_sub_dir{1} = 'LI';
roi_int{1} = { 'ILF' 'tSLF' 'Unc' };

% mes_dir{4} = 'DTI'; mes_typ_dir{4} = 'fiber_MD'; mes_sub_dir{4} = 'Raw';
% roi_int{4} = { 'xL.ILF' 'xL.tSLF' 'xL.Unc' 'xR.ILF' 'xR.tSLF' 'xR.Unc' };

mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm'; mes_sub_dir{2} = 'LI';
roi_int{2} = { 'temporalpole' 'superiortemporal' 'parstriangularis' 'parsopercularis' 'entorhinal' 'fusiform' 'cuneus' 'precuneus' };

% mes_dir{6} = 'DTI'; mes_typ_dir{6} = 'wmparc_MD_wm'; mes_sub_dir{6} = 'Raw';
% roi_int{6} = { 'xlh.temporalpole' 'xlh.superiortemporal' 'xlh.parstriangularis' 'xlh.parsopercularis' 'xlh.entorhinal' 'xlh.fusiform' 'xlh.cuneus' 'xlh.precuneus' ...
%     'xrh.temporalpole' 'xrh.superiortemporal' 'xrh.parstriangularis' 'xrh.parsopercularis' 'xrh.entorhinal' 'xrh.fusiform' 'xrh.cuneus' 'xrh.precuneus' };

mes_dir{3} = 'fMRI'; mes_typ_dir{3} = 'alicia'; mes_sub_dir{3} = 'LI';
roi_int{3} = { 'LateralTemporal02.cortical' 'InferiorFrontal.cortical' 'VentralTemporal.cortical' };
% 
% mes_dir{8} = 'fMRI'; mes_typ_dir{8} = 'destr'; mes_sub_dir{8} = 'Raw';
% roi_int{8} = { 'lh.Pole.temporal' 'lh.G.and.S.cingul.Mid.Post' 'lh.G.cuneus' 'lh.G.front.inf.Triangul' 'lh.G.front.inf.Opercular' 'lh.S.intrapariet.and.P.trans' ...
%     'rh.Pole.temporal' 'rh.G.and.S.cingul.Mid.Post' 'rh.G.cuneus' 'rh.G.front.inf.Triangul' 'rh.G.front.inf.Opercular' 'rh.S.intrapariet.and.P.trans'};

mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'cort_thick_ctx'; mes_sub_dir{4} = 'LI';
roi_int{4} = { 'temporalpole' 'superiortemporal' 'parstriangularis' 'parsopercularis' 'entorhinal' 'fusiform' 'cuneus' 'precuneus' };

mes_dir{5} = 'MRI'; mes_typ_dir{5} = 'subcort_vol_ICV_cor'; mes_sub_dir{5} = 'LI';
roi_int{5} = { 'Hippocampus' };

mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_ctx'; mes_sub_dir{6} = 'LI';
roi_int{6} = {'temporalpole' 'superiortemporal' 'parstriangularis' 'parsopercularis' 'entorhinal' 'fusiform' 'cuneus' 'precuneus' };

mes_dir{7} = 'rsfMRI'; mes_typ_dir{7} = 'var_vol'; mes_sub_dir{7} = 'LI';
roi_int{7} = { 'Hippocampus' };

% Concatenate
for iG = 1:numel(nme_typ)
    
    for iM = 1:numel(mes_dir)
                
        rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_n.csv' ]);
        
        for iC = 1:numel(cog_scr_nme{iG})
            
            cog_col = find(ismember( rvl_hld(1,:), cog_scr_nme{iG}{iC}));
            
            row_cnt = 1;
            for iR = 1:numel(roi_int{iM})
                
                neu_col = find(ismember( rvl_hld(:,1), roi_int{iM}{iR}));
                
                rvl_out_hld{iG}{iM}{row_cnt,1} = mes_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,2} = mes_typ_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,3} = mes_sub_dir{iM};
                rvl_out_hld{iG}{iM}{row_cnt,4} = roi_int{iM}{iR};
                
                rvl_out_hld{iG}{iM}{row_cnt, iC+4} = num2str(roundsd(rvl_hld{neu_col, cog_col},2));
                if rvl_hld{neu_col, cog_col} > 0
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}(2:end);
                else
                    rvl_out_hld{iG}{iM}{row_cnt, iC+4} = rvl_out_hld{iG}{iM}{row_cnt, iC+4}([1 3:end]);
                end
                if pvl_hld{neu_col, cog_col}<=.05;                                          rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' *']; end
                if pvl_hld{neu_col, cog_col}<=.10 && pvl_hld{neu_col, cog_col}>.05; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' #']; end
                if pvl_hld{neu_col, cog_col}>.10; rvl_out_hld{iG}{iM}{row_cnt, iC+4} = [rvl_out_hld{iG}{iM}{row_cnt, iC+4} ' ']; end
                
                rvl_str =  num2str(roundsd(rvl_hld{neu_col, cog_col},2));
                rvl_str = rvl_str(2:end);
                num_str = num2str(num_hld{neu_col, cog_col});
                pvl_str =  num2str(roundsd(pvl_hld{neu_col, cog_col},2));
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
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'rvalues_raw' '_' nme_typ{iG} '.csv' ], [ {''} {''} {''} {''}  cog_scr_nme{iG} ; rvl_out_hld{iG} ])
    
    tot_out_hld{iG} = cat(1, tot_out_hld{iG}{:});
    cell2csv([ dta_dir '/' 'Summary' '/' 'apriori' '/' 'total_raw' '_' nme_typ{iG} '.csv' ], [ {''} {''} {''} {''}  cog_scr_nme{iG} ; tot_out_hld{iG} ])
    
    clear rvl_out_hld tot_out_hld
    
end

%% Look for any missed good ROIs
rvl_dir = dta_dir;
out_dir = [dta_dir '/' 'Summary' '/' 'check'];

ejk_chk_dir(out_dir)

mes_dir = { 'DTI' ...
            'fMRI' ...
            'MRI' ...
            'rsfMRI'};
mes_typ = { { 'fiber_FA' 'wmparc_FA_wm' } ...
            { 'alicia'  } ...
            { 'cort_thick_ctx' 'subcort_vol_ICV_cor' } ...
            { 'var_ctx' 'var_vol' } };
mes_sub = { 'Raw' 'LI' };

nme_typ = { 'left_pre' 'right_pre' 'left_post' 'right_post'};
tst_nme = { { 'bnt.raw.scr' 'ant.mem.raw.scr' } ...
            { 'bnt.raw.scr' 'ant.mem.raw.scr' } ...
            { 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' } ...
            { 'bnt.raw.scr.pst' 'ant.mem.raw.scr.pst' } };


for iG = 1:numel(nme_typ)
    for iMD = 1:numel(mes_dir)
        for iMT = 1:numel(mes_typ{iMD})
            for iMS = 1:numel(mes_sub)
                
                lod_rvl = mmil_readtext( [ rvl_dir '/' mes_dir{iMD} '/' mes_typ{iMD}{iMT} '/' mes_sub{iMS} '/' nme_typ{iG} '/' 'cross_correlation_rvalues.csv'  ] );
                lod_rvl(strcmpi(lod_rvl,'NA')) = {NaN};
                
                [ ~, col_ind] = intersect( lod_rvl(1,:), tst_nme{iG});
                
                [~, avg_ind_hld] = sort(nanmean(cell2mat(lod_rvl(2:end,col_ind)),2));
                
                if numel(avg_ind_hld) > 10
                    out_rvl = lod_rvl( [1 avg_ind_hld(1:5)'+1 avg_ind_hld(end-4:end)'+1], :);
                    cell2csv( [ out_dir '/' mes_dir{iMD} '_' mes_typ{iMD}{iMT}  '_' mes_sub{iMS} '_' nme_typ{iG} '.csv' ], out_rvl)
                    clear out_rvl
                else
                    out_rvl = lod_rvl( [1 avg_ind_hld'+1], :);
                    cell2csv( [ out_dir '/' mes_dir{iMD} '_' mes_typ{iMD}{iMT}  '_' mes_sub{iMS} '_' nme_typ{iG} '.csv' ], out_rvl)
                    clear out_rvl
                end
                
            end
        end
    end
end

%% Initial Summarize
clear mes_dir

out_dir = [dta_dir '/' 'Summary' '/' ];

mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA'; mes_sub_dir{1} = 'LI';
mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm'; mes_sub_dir{2} = 'LI';
mes_dir{3} = 'fMRI'; mes_typ_dir{3} = 'alicia'; mes_sub_dir{3} = 'LI';
mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'cort_thick_ctx'; mes_sub_dir{4} = 'LI';
mes_dir{5} = 'MRI'; mes_typ_dir{5} = 'subcort_vol_ICV_cor'; mes_sub_dir{5} = 'LI';
roi_int{5} = { 'xLeft.Hippocampus' 'xRight.Hippocampus' };
mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_ctx'; mes_sub_dir{6} = 'LI';
mes_dir{7} = 'rsfMRI'; mes_typ_dir{7} = 'var_vol'; mes_sub_dir{7} = 'LI';

nme_typ = { 'left_pre' 'right_pre' 'left_post' 'right_post'};

dta_col_use = { { 5:6                           5:6                           7          } ...
                { 5:6                           5:6                           7          } ...
                { 5:6                           5:6                           7          } ...
                { 8:9                           8:9                           10          } ...
                { 5:6                           5:6                           7          } ...
                { 8:9                           8:9                           10          } };
dta_col_typ = { { 'average'                     'subtract'                    'average'   } ...
                { 'average'                     'subtract'                    'average'   } ...
                { 'average'                     'subtract'                    'average'   } ...
                { 'average'                     'subtract'                    'average'   } ...
                { 'average'                     'subtract'                    'average'   } ...
                { 'average'                     'subtract'                    'average'   } };
dta_out_nme = { { 'Pre_Total_NamingTogether'    'Pre_Total_NamingSeparate'    'Pre_Total_CategoryFluency'} ...
                { 'Pre_Total_NamingTogether'    'Pre_Total_NamingSeparate'    'Pre_Total_CategoryFluency'} ...
                { 'Pre_Surgery_NamingTogether'  'Pre_Surgery_NamingSeparate'  'Pre_Surgery_CategoryFluency'} ...
                { 'Post_Surgery_NamingTogether' 'Post_Surgery_NamingSeparate' 'Post_Surgery_CategoryFluency'} ...
                { 'Pre_Surgery_NamingTogether'  'Pre_Surgery_NamingSeparate'  'Pre_Surgery_CategoryFluency'} ...
                { 'Post_Surgery_NamingTogether' 'Post_Surgery_NamingSeparate' 'Post_Surgery_CategoryFluency'} };

for iG = 1:numel(nme_typ)
    
    out_rvl = cell(0);
    for iM = 1:numel(mes_dir)
        
        rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_rvalues.csv' ]);
        rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
        
        out_rvl_add      = repmat(mes_dir(iM),size(rvl_hld,1)-1,1);
        out_rvl_add(:,2) = repmat(mes_typ_dir(iM),size(rvl_hld,1)-1,1);
        out_rvl_add(:,3) = repmat(mes_sub_dir(iM),size(rvl_hld,1)-1,1);
        out_rvl_add      = [out_rvl_add rvl_hld(2:end,:)];
        
        out_rvl = [ out_rvl ; out_rvl_add];
        
    end    
    
    for iT = 1:numel(dta_col_use{iG})
        fcfg         = [];
        fcfg.dta     = out_rvl;
        fcfg.dta_col = dta_col_use{iG}{iT};
        fcfg.typ     = dta_col_typ{iG}{iT};
        dta_out = ejk_roi_summarize(fcfg);
        
        cell2csv( [ out_dir '/' dta_out_nme{iG}{iT} '_' nme_typ{iG} '.csv'], [ {''} {''} {''} {''}  rvl_hld(1,2:end)  dta_col_typ{iG}{iT} ; dta_out]);
        clear dta_out
    end
    
    clear rvl_hld
    
end

%% Brain ROI Summarize
nme_typ = { 'left_pre' 'right_pre' 'left_post' 'right_post'};

mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA';            mes_sub_dir{1} = 'Raw'; roi_use{1} = [0];
mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm';        mes_sub_dir{2} = 'Raw'; roi_use{2} = [2];
mes_dir{3} = 'MRI'; mes_typ_dir{3} = 'cort_thick_ctx';      mes_sub_dir{3} = 'Raw'; roi_use{3} = [2];
mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'subcort_vol_ICV_cor'; mes_sub_dir{4} = 'Raw'; roi_use{4} = [0];
mes_dir{5} = 'rsfMRI'; mes_typ_dir{5} = 'var_ctx';          mes_sub_dir{5} = 'Raw'; roi_use{5} = [2];
mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_vol';          mes_sub_dir{6} = 'Raw'; roi_use{6} = [0];

pvl_cut = .10;

for iM = 5:numel(mes_dir)
    if roi_use{iM} > 0
        for iG = 1:numel(nme_typ)
            
            out_dir_mes = [dta_dir '/' 'Summary' '/' 'ROI_plots' '/' mes_dir{iM} '/'];
            ejk_chk_dir(out_dir_mes);
            
            rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_rvalues.csv' ]);
            rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
            
            pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' nme_typ{iG} '/' 'cross_correlation_pvalues.csv' ]);
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
                cfg.out_nme = [ mes_typ_dir{iM} '_' mes_sub_dir{iM} '_' nme_typ{iG} '_' mmil_spec_char(rvl_hld{1,iC},{'.'})  ];
                
                ejk_roi_display_plot(cfg);
                
            end
            
        end
    end
end

    
%% Brain ROI Summarize
% pvl_cut = 0.10;
% 
% %
% sph = { 'lh' 'rh' };
% sph_vew = { 'lat' 'med' 'ven' };
% plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};
% 
% srf_brn{1} = fs_read_surf(['/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'fsaverage' '/' 'surf' '/' 'lh.pial']);
% srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
% srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
% srf_brn{2} = fs_read_surf(['/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'fsaverage' '/' 'surf' '/' 'rh.pial']);
% srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
% srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;
% 
% col{1} = rgb('purplish blue');
% col{2} = rgb('blue');
% col{3} = rgb('grayish blue');
% col{4} = rgb('light grey');
% col{5} = rgb('light grey');
% col{6} = rgb('light grey');
% col{7} = rgb('light grey');
% col{8} = rgb('light grey');
% col{9} = rgb('reddish grey');
% col{10} = rgb('red');
% col{11} = rgb('orangish red');
% col_map = [];
% top_pct = 1;
% for iC = 1:numel(col)-1
%     col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
% end
% top_pct = 0.60;
% 
% %
% clear mes_dir
% 
% out_dir = [dta_dir '/' 'Summary' '/' 'ROI_plots' '/' ];
% ejk_chk_dir(out_dir);
% 
% mes_dir{1} = 'DTI'; mes_typ_dir{1} = 'fiber_FA';            mes_sub_dir{1} = 'Raw'; roi_use{1} = [0];
% mes_dir{2} = 'DTI'; mes_typ_dir{2} = 'wmparc_FA_wm';        mes_sub_dir{2} = 'Raw'; roi_use{2} = [2];
% mes_dir{3} = 'MRI'; mes_typ_dir{3} = 'cort_thick_ctx';      mes_sub_dir{3} = 'Raw'; roi_use{3} = [2];
% mes_dir{4} = 'MRI'; mes_typ_dir{4} = 'subcort_vol_ICV_cor'; mes_sub_dir{4} = 'Raw'; roi_use{4} = [0];
% mes_dir{5} = 'rsfMRI'; mes_typ_dir{5} = 'var_ctx';          mes_sub_dir{5} = 'Raw'; roi_use{5} = [2];
% mes_dir{6} = 'rsfMRI'; mes_typ_dir{6} = 'var_vol';          mes_sub_dir{6} = 'Raw'; roi_use{6} = [0];
% 
% for iM = 1:numel(mes_dir)
%     if roi_use{iM} ~= 0
%         
%         out_dir_mes = [dta_dir '/' 'Summary' '/' 'ROI_plots' '/' mes_dir{iM} '/'];
%         ejk_chk_dir(out_dir_mes);
%         
%         rvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' 'cross_correlation_rvalues.csv' ]);
%         rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
%         
%         pvl_hld = mmil_readtext([ dta_dir '/' mes_dir{iM} '/' mes_typ_dir{iM} '/' mes_sub_dir{iM} '/' 'cross_correlation_pvalues.csv' ]);
%         pvl_hld(strcmpi(pvl_hld,'NA')) = {NaN};
%         
%         for iC = 2:size(rvl_hld,2)
%             
%             fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
%             
%             for iH = 1:2
%                 for iSP = 1:numel(sph_vew)
%                     % Setup Data
%                     
%                     if iH == 1
%                         plt_tbl = rvl_hld( string_find(rvl_hld(:,1),'lh.'), [1 iC]);
%                         pvl_tbl = pvl_hld( string_find(pvl_hld(:,1),'lh.'), [1 iC]);
%                         plt_tbl(cell2mat(pvl_tbl(:,2))>pvl_cut,2) = {0};
%                         plt_tbl(:,2) = num2cell((cell2mat(plt_tbl(:,2))+1)/2);
%                     elseif iH == 2
%                         plt_tbl = rvl_hld( string_find(rvl_hld(:,1),'rh.'), [1 iC]);
%                         pvl_tbl = pvl_hld( string_find(pvl_hld(:,1),'rh.'), [1 iC]);
%                         plt_tbl(cell2mat(pvl_tbl(:,2))>pvl_cut,2) = {0};
%                         plt_tbl(:,2) = num2cell((cell2mat(plt_tbl(:,2))+1)/2);
%                     end
%                     
%                     %                     if strcmpi(plt_tbl{1,1}(1), 'x')
%                     plt_tbl(:,1) = cellfun( @(x) x(2:end), plt_tbl(:,1), 'uni', 0);
%                     %                     else
%                     %                         plt_tbl(:,1) = cellfun( @(x) x(4:end), plt_tbl(:,1), 'uni', 0);
%                     %                     end
%                     
%                     % Plot
%                     pcfg = [];
%                     
%                     pcfg.surf_brain  = srf_brn{iH};
%                     pcfg.aparc       = ['/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' sph{iH} '.aparc.annot'];
%                     
%                     pcfg.sph         = sph{iH};
%                     pcfg.sph_vew     = sph_vew{iSP};
%                     
%                     pcfg.label       = 0;
%                     pcfg.radius      = [];
%                     pcfg.alpha       = 1;
%                     
%                     pcfg.non_ele     = [];
%                     pcfg.sve_img     = 0; % ###
%                     
%                     pcfg.fig_hdl = fig_hld(1);
%                     pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
%                     
%                     pcfg.col_map = col_map;
%                     
%                     pcfg.tbl_pct = cell2mat(plt_tbl(:,2)) ./ top_pct;
%                     pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
%                     
%                     pcfg.tbl_loc = plt_tbl(:,1);
%                     
%                     pcfg.top_pct = 0.5;
%                     
%                     nyu_plot2(pcfg);
%                     
%                 end
%             end
%             
%             ax1 = axes('OuterPosition',[0.85 0.20 0.04 0.60],'visible','off','Parent',fig_hld(1));
%             
%             colormap(ax1,col_map)
%             clb = colorbar('west','Position',[0.92 0.20 0.02 0.60]);
%             clb.TickLength = 0;
%             col_num = [ num2cell(roundsd(linspace(-top_pct,-top_pct,5),2)) ...
%                 {0} ...
%                 num2cell(roundsd(linspace(top_pct,cfg.hgh_rng_num(2),5),2)) ];
%             clb.TickLabels = cellfun(@num2str,col_num,'uni',0);
%             
%             print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/ROI/' '/' 'ROI_guide' '_' sph{iH} '_' sph_vew{iSP} '.png'],'-dpng','-r200')
%             close all
%             
%         end
%         
%     end
% end





