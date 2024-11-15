slc_btc = 10;

run_typ = { 'total' 'dataset' 'site' };
ana_typ = { 'diagnosis' 'lateralization' };
cmp_typ = { { 'TLE_vs_HC' } { 'LeftTLE_vs_HC' 'RightTLE_vs_HC' 'LeftTLE_vs_RightTLE' } };
grp_cmp = { { {'EPD' 'HC'} } { {'left' 'HC'} {'right' 'HC'} {'left' 'right'} } };

%% Get info
fcfg = [];
fcfg.dta_loc = [ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'];
fcfg.dta_col = 2;
[ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

fle_nme = mmil_readtext([ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'files' '_' dte_str '.csv']);
exp_sbj = load_nifti(fle_nme{1});

%% Folder & VBM out setup
ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str] )
for iRT = 1:numel(run_typ); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/'] );
    for iAT = 1:numel(ana_typ); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/'] ); 
        for iCT = 1:numel(cmp_typ{iAT}); ejk_chk_dir( [ prj_dir '/' 'vbm' '/' dte_str '/' run_typ{iRT} '/' ana_typ{iAT} '/' cmp_typ{iAT}{iCT} '/'] ); 
        
            use_nme = fieldnames(grp.(run_typ{iRT}));
            for iUN = 1:numel(use_nme)
                vbm_out_tvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN}) = nan(size(exp_sbj.vol));
                vbm_out_pvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN}) = nan(size(exp_sbj.vol));
%                 vbm_out_deg.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}) = nan(size(exp_sbj.vol));
            end
        end; end; end

%% Batch Setup
slc_btc_cut = [ 1:slc_btc:size(exp_sbj.vol,2) size(exp_sbj.vol,2)+1 ];

for iSB = 1:numel(slc_btc_cut)-1

    dta_hld = nan( size(exp_sbj.vol,1), slc_btc_cut(iSB+1) - slc_btc_cut(iSB), size(exp_sbj.vol,3), numel(cov_sbj));

    %% Load Data in batches
    fprintf("Load nii | Slice %i to %i\n\n",slc_btc_cut(iSB), slc_btc_cut(iSB+1)-1)
    for iS = 1:numel(cov_sbj)

        if rem(iS,200)==0
            fprintf("Load nii | Slice %i to %i | subject #%i\n", slc_btc_cut(iSB), slc_btc_cut(iSB+1)-1, iS)
        end

        nii_hld = load_nifti(fle_nme{iS});
        dta_hld(:,:,:,iS) = nii_hld.vol(:,slc_btc_cut(iSB):slc_btc_cut(iSB+1)-1,:);

    end

    %% Run VBM in batches
    for iRT = 1:numel(run_typ)
        use_nme = fieldnames(grp.(run_typ{iRT}));
        for iUN = 1:numel(use_nme)
            for iAT = 1:numel(ana_typ)
                for iCT = 1:numel(cmp_typ{iAT})

                    fprintf("\n\nVBM: %i of %i | %i of %i | %i of %i | %i of %i | Slice %i to %i\n\n",iRT,numel(run_typ),iUN,numel(use_nme),iAT,numel(ana_typ),iCT,numel(cmp_typ{iAT}),slc_btc_cut(iSB), slc_btc_cut(iSB+1)-1)

                    if size(dta_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1})),4) > 1 && ...
                            size(dta_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2})),4) > 1 && ...
                            numel(grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).HC) > 0

                        %
                        fcfg = [];
                        fcfg.dta_one = dta_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{1}));
                        fcfg.dta_two = dta_hld(:,:,:,grp.(run_typ{iRT}).(use_nme{iUN}).(ana_typ{iAT}).(grp_cmp{iAT}{iCT}{2}));
                        vbm_out = ejk_vbm(fcfg);

                        vbm_out_tvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN})(:,slc_btc_cut(iSB):slc_btc_cut(iSB+1)-1,:) = vbm_out.tvl;
                        vbm_out_pvl.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN})(:,slc_btc_cut(iSB):slc_btc_cut(iSB+1)-1,:) = vbm_out.pvl;
%                         vbm_out_deg.(run_typ{iRT}).(ana_typ{iAT}).(cmp_typ{iAT}{iCT}).(use_nme{iUN})(:,slc_btc_cut(iSB):slc_btc_cut(iSB+1)-1,:) = vbm_out.df;
                    end
                end
            end
        end
    end

end

save( [ prj_dir '/' 'vbm' '/' dte_str '/' 'vbm_all.mat'], 'vbm_out_tvl','vbm_out_pvl', '-v7.3'); %'vbm_out_deg'



