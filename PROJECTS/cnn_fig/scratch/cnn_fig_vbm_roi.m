
load( [ prj_dir '/' 'vbm' '/' dte_str '/' 'vbm_all.mat'], 'vbm_out_tvl','vbm_out_pvl');

%% Extract
for iRT = 1:numel(run_typ)
    use_nme = fieldnames(grp.(run_typ{iRT}));
    for iUN = 1:numel(use_nme)
        for iAT = 1:numel(ana_typ)
            for iCT = 1:numel(cmp_typ{iAT})


                for iA = 1:numel(atl_ovr_nme)

                    atl_nii_dta = niftiread([ atl_dir '/' '113' '/' atl_ovr_nme{iA} '.nii' ]);
                    if size(atl_nii_dta,2) == 138; atl_nii_dta = atl_nii_dta(:,1:137,:); elseif size(atl_nii_dta,2) ~= 138 || size(atl_nii_dta,2) ~= 137; error("Size mismatch"); end

                    fcfg = [];
                    fcfg.dta_loc = [ atl_dir '/' 'labels' '/' atl_ovr_nme{iA} '.csv'];
                    fcfg.delim = ';';
                    fcfg.dta_col = 2;
                    [ atl_dta, atl_num, atl_col] = ejk_dta_frm( fcfg );
                    atl_num = cell2mat(atl_num);

                    atl_dta_out = cell(numel(atl_num),4);
                    for iT = 1:numel(atl_num)
                        
                        atl_dta_out{iT,1} = atl_dta{iT,strcmpi(atl_col,'ROIname')};
                        atl_dta_out{iT,2} = atl_dta{iT,strcmpi(atl_col,'ROIabbr')};

                        mtc_ary = atl_nii_dta==atl_num(iT);
                        mtc_ary_num = sum(mtc_ary(:));

                        atl_dta_out{iT,3} = nansum(vbm_out_tvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN})(mtc_ary)) / mtc_ary_num;
                        atl_dta_out{iT,4} = nansum(vbm_out_pvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN})(mtc_ary)) / mtc_ary_num;

                    end

                    cell2csv( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/' cmp_typ{iAT}{iCT} '/' atl_ovr_nme{iA} '_output.csv'], atl_dta_out );

                end
            end
        end
    end
end

%%