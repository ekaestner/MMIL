
slc_ind = 71:2:79;

run_typ = { 'total' 'dataset' 'site' };
ana_typ = { 'diagnosis' 'lateralization' };
cmp_typ = { { 'TLE_vs_HC' } { 'LeftTLE_vs_HC' 'RightTLE_vs_HC' 'LeftTLE_vs_RightTLE' } };
grp_cmp = { { {'EPD' 'HC'} } { {'left' 'HC'} {'right' 'HC'} {'left' 'right'} } };

%%
ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str] )
for iRT = 1:numel(run_typ); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/'] )
    for iAT = 1:numel(ana_typ); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/'] )
        for iCT = 1:numel(cmp_typ{iAT}); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/' cmp_typ{iAT}{iCT} '/'] ); end; end; end

%%
fcfg = [];
fcfg.dta_loc = [ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'];
fcfg.dta_col = 2;
[ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

fle_nme = mmil_readtext([ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'files' '_' dte_str '.csv']);

%%
mtl_hld = NaN(113,numel(slc_ind),113,numel(cov_sbj));

fprintf("Load nii\n\n")
for iS = 1:numel(cov_sbj)
    
    if rem(iS,50)==0
        fprintf("Load nii | loading | subject #%i\n", iS)
    end

    nii_hld = load_nifti(fle_nme{iS});
    mtl_hld(:,:,:,iS) = nii_hld.vol(:,slc_ind,:);

end

%%
for iRT = 1:numel(run_typ)
    use_nme = fieldnames(grp.(run_typ{iRT}));
    for iUN = 1:numel(use_nme)
        for iAT = 1:numel(ana_typ)
            for iCT = 1:numel(cmp_typ{iAT})

                if size(mtl_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1})),4) > 1 && ...
                        size(mtl_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2})),4) > 1 && ...
                        numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).HC) > 0

                %
                fcfg = [];
                fcfg.dta_one = mtl_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1}));
                fcfg.dta_two = mtl_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2}));
                vbm_out = ejk_vbm(fcfg);

                %
                avg_brn = nanmean(mtl_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).HC),4);
                avg_brn_cut = 0.15;

                slc = 1:size(mtl_hld,2);

                ylm_low = -12;
                ylm_hgh = 12;

                if (numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1})) + numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2}))) > 1000
                    pvl_msk = .0001;
                elseif (numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1})) + numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2}))) > 100
                    pvl_msk = .001;
                else
                    pvl_msk = .01;
                end

                fcfg = [];

                fcfg.plt_dta = vbm_out.tvl;

                fcfg.brn_dta = avg_brn;
                fcfg.brn_cut_off = avg_brn_cut;

                fcfg.pvl_dta = vbm_out.pvl;
                fcfg.pvl_msk = pvl_msk;

                fcfg.slc      = slc;

                fcfg.ylm_low = ylm_low;
                fcfg.ylm_hgh = ylm_hgh;

                fcfg.plt_dir = [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/' cmp_typ{iAT}{iCT} '/'];
                fcfg.plt_nme = [ 'montage' '_'  run_typ{iRT} '_' ana_typ{iAT} '_' cmp_typ{iAT}{iCT} '_' use_nme{iUN} '_' grp_cmp{iAT}{iCT}{1} '_VS_' grp_cmp{iAT}{iCT}{2} '_' num2str(pvl_msk) '_' dte_str];

                ejk_montage(fcfg)

                end
            end
        end
    end
end

