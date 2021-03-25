function mmil_app_all_nse(cfg,dat)

for iC = 1:size(cfg.app_all,2)
    
    inp_ind = cfg.app_all(1,iC);
    out_ind = cfg.app_all(2,iC);
    
    for iPR = 1:numel(cfg.pre_fix{inp_ind})
        if ~exist([cfg.out_dir '/' cfg.pre_fix{out_ind}{iPR}]); mkdir([cfg.out_dir '/' cfg.pre_fix{out_ind}{iPR}]); end
        copyfile([cfg.out_dir '/' cfg.pre_fix{inp_ind}{iPR} '/' 'rmv_num.mat'],[cfg.out_dir '/' cfg.pre_fix{out_ind}{iPR} '/' 'rmv_num.mat']);
    end

end