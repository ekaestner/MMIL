function vbm_out = ejk_vbm(cfg)

%% 
vbm_sze = size(cfg.dta_one);

tvl_hld = nan(vbm_sze(1),vbm_sze(2),vbm_sze(3));
pvl_hld = nan(vbm_sze(1),vbm_sze(2),vbm_sze(3));
df_hld  = nan(vbm_sze(1),vbm_sze(2),vbm_sze(3));

%% 
tic
for i1 = 1:vbm_sze(1) % parfor?
    for  i2 = 1:vbm_sze(2)
        tic
        %         % Pre-allocate arrays for efficiency within the worker (optional)
        %         tvl_hld_local = nan(vbm_sze(2),vbm_sze(3));
        %         pvl_hld_local = nan(vbm_sze(2),vbm_sze(3));
        %         df_hld_local = nan(vbm_sze(2),vbm_sze(3));

        for i3 = 1:vbm_sze(3)
            [~, pvl, ~, stt] = ttest2(cfg.dta_one(i1, i2, i3, :), cfg.dta_two(i1, i2, i3, :));
            tvl_hld(i1,i2,i3) = stt.tstat;
            pvl_hld(i1,i2,i3) = pvl;
            df_hld(i1,i2,i3) = stt.df;
        end

    end

    % Assign local arrays to output arrays after loop (avoid race conditions)
    %     tvl_hld(i1, :, :) = tvl_hld_local;
    %     pvl_hld(i1, :, :) = pvl_hld_local;
    %     df_hld(i1, :, :) = df_hld_local;

    if rem(i1,25)==0
        dsp_tme = toc;
        fprintf('%i out of %i, Time: %f\n',i1,vbm_sze(1),dsp_tme)
        tic
    end
end

vbm_out.tvl = tvl_hld;
vbm_out.pvl = pvl_hld;
vbm_out.df  = df_hld;

end