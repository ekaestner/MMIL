function [ dta_out, wrn_out ] = ejk_extract_mmps_roi( cfg )

if ~isfield(cfg,'roi_nme'); cfg.roi_nme = []; end

wrn_out = cell(0); wrn_cnt = 1;

%% Load Files
roi_hld = mmil_readtext( cfg.fle_nme );
    roi_hld(1,:) = cellfun(@(x) mmil_spec_char(x,{'-'}),roi_hld(1,:),'uni',0);
    rcn_col     = find(strcmpi( roi_hld(1,:), 'VisitID' ));
    fld_str_col = find(strcmpi( roi_hld(1,:), 'MagneticFieldStrength' ));
    mps_col     = find(strcmpi( roi_hld(1,:), 'MMPS_version' ) | strcmpi( roi_hld(1,:), 'AFNI_version' ));
    
rcn_hld = mmil_readtext( cfg.rcn_nme );

for iR = 1:numel(cfg.roi_nme); [ ~, roi_nme_hld{iR} ] = fs_read_annotation( cfg.roi_nme{iR} ); end

%% Find Column
if ~isempty(cfg.roi_nme)
    col_use = [];
    col_plc_ind = 4;
    for iR = 1:numel(cfg.roi_nme)        
        col_use = [col_use find(ismember( roi_hld(1, :), cellfun(@(x) mmil_spec_char(x,{'-'}),strcat(cfg.mes_nme, '_', cfg.roi_hms{iR}, '_', roi_nme_hld{iR}),'uni',0)))];  
        col_plc{iR} = col_plc_ind(iR) + (1:numel(find(ismember( roi_hld(1, :), strcat(cfg.mes_nme, '_', cfg.roi_hms{iR}, '_', roi_nme_hld{iR})))));
        col_plc_ind(iR+1) = col_plc_ind(iR) + numel(find(ismember( roi_hld(1, :), strcat(cfg.mes_nme, '_', cfg.roi_hms{iR}, '_', roi_nme_hld{iR}))));
    end
    col_plc = cat(2,col_plc{:});
else
    col_plc_ind = 4;
    col_use = string_find( roi_hld(1, :), cfg.mes_nme);
    col_plc = col_plc_ind + (1:numel(col_use));
end

%% Extract Data
dta_out = cell( size(cfg.sbj_nme,1)+1, numel(col_use)+4);

dta_out(1,1) = {'sbj_nme'};
dta_out(1,2) = {'recon'};
dta_out(1,3) = {'field_strength'};
dta_out(1,4) = {'mmps'};
dta_out(1,col_plc) = roi_hld(1, col_use);

for iS = 1:numel(cfg.sbj_nme)
    
    sbj_row = find( strcmpi( roi_hld(:,1), cfg.sbj_nme{iS}) ) ;
    if isempty(sbj_row); sbj_row = string_find( roi_hld(:,1), cfg.sbj_nme(iS)); end
    
    sbj_rcn     = rcn_hld( strcmpi( rcn_hld(:,1), cfg.sbj_nme{iS}), 3 ) ;
    if ~isempty(sbj_rcn) && ~isempty(sbj_rcn{1}) && ~isempty(sbj_row)
        sbj_rcn_num = strnearest( sbj_rcn, roi_hld( sbj_row, rcn_col));
        sbj_row = sbj_row( sbj_rcn_num{1});
        
        if ~strcmpi( sbj_rcn{1}, roi_hld{sbj_row,rcn_col})
            
            exp_cut = strfind( sbj_rcn{1}, '_');
            obt_cut = strfind( roi_hld{sbj_row,rcn_col}, '_');
            if      numel(exp_cut) > numel(obt_cut)
                exp_hld = sbj_rcn{1}(1:numel(roi_hld{sbj_row,rcn_col}));
                obt_hld = roi_hld{sbj_row,rcn_col};
            elseif  numel(exp_cut) < numel(obt_cut)
                exp_hld = sbj_rcn{1};
                obt_hld = roi_hld{sbj_row,rcn_col}(1:numel(sbj_rcn{1}));
            else
                exp_hld = sbj_rcn{1};
                obt_hld = roi_hld{sbj_row,rcn_col};
            end
            
            if ~strcmpi(exp_hld,obt_hld)
            wrn_out{wrn_cnt,1} = ['EXPECTED: ' sbj_rcn{1} ' ; ' 'OBTAINED: ' roi_hld{sbj_row,rcn_col} ];
            wrn_cnt = wrn_cnt + 1;
            end
        end
            
        dta_out(iS+1,1) = roi_hld(sbj_row,1);
        dta_out(iS+1,2) = roi_hld(sbj_row,rcn_col);
        if ~isempty(fld_str_col); dta_out(iS+1,3) = roi_hld(sbj_row,fld_str_col); end;
        if ~isempty(mps_col); dta_out(iS+1,4) = roi_hld(sbj_row,mps_col); end;
        dta_out(iS+1,col_plc) = roi_hld(sbj_row, col_use);
    else
        dta_out(iS+1,1) = cfg.sbj_nme(iS);
        dta_out(iS+1,2) = {''};
        dta_out(iS+1,3) = {''};
        dta_out(iS+1,4) = {''};
        dta_out(iS+1,col_plc) = num2cell(nan(1, numel(col_use)));
    end    
end

end