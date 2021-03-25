

function [ cnn_dta , row_nme , col_nme ] = ejk_extract_connectome(cfg)

%%

if exist([cfg.ovr_dir '/' cfg.sbj_nme '_norm.csv']) || exist([cfg.ovr_dir '/' cfg.sbj_nme([1:3 5:end]) '_norm.csv'])
    
    % Load Connectome
    try
        cnn_dat_hld = mmil_readtext([cfg.ovr_dir '/' cfg.sbj_nme '.csv']); % cnn_dat_hld = mmil_readtext([cfg.ovr_dir '/' cfg.sbj_nme '_norm.csv']);
    catch
        cnn_dat_hld = mmil_readtext([cfg.ovr_dir '/' cfg.sbj_nme([1:3 5:end]) '.csv']); % cnn_dat_hld = mmil_readtext([cfg.ovr_dir '/' cfg.sbj_nme([1:3 5:end]) '_norm.csv']);
    end
    
    % Extract Data
    cnn_dta = cell2mat(cnn_dat_hld(2:end,2:end));
    
    % Extract Label Name
    row_nme = cnn_dat_hld(1,2:end);
    col_nme = cnn_dat_hld(2:end,1);
    
else
    cnn_dta = [];
    row_nme = [];
    col_nme = [];
end

end