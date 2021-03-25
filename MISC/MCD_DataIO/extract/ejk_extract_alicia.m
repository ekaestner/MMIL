function ejk_extract_alicia(cfg)

%% Alicia 
tot_roi = dir([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' '*.csv']); tot_roi = {tot_roi(:).name};

rep_lbl = 'a2009s';
tot_roi = tot_roi(string_find(tot_roi,{rep_lbl}));

for iTR = 1:numel(tot_roi)
        
    bse_lbl = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' tot_roi{iTR}]);
    
    msk_crt = dir([cfg.prj_dir '/' 'PARCELLATIONS' '/' 'Alicia_ROI_Holder']); msk_crt = {msk_crt(:).name}; msk_crt = msk_crt(string_find(msk_crt,'.csv')); msk_crt(string_find(msk_crt,'.~')) = []; msk_crt(strcmpi(msk_crt,'AliciaKey.csv')) = [];
    
    for iF = 1:numel(msk_crt)
        
        out_msk{1,iF+1} = mmil_spec_char(msk_crt{iF}(1:end-4),{'.'});
        
        for iS = 2:size(bse_lbl,1)
            
            out_msk{iS,1} = bse_lbl{iS,1};
            
            msk_hld = mmil_readtext([cfg.prj_dir '/' 'PARCELLATIONS' '/' 'Alicia_ROI_Holder' '/' msk_crt{iF}]);
            
            [~,dat_loc,~] = intersect( cellfun(@(x) mmil_spec_char(x,{'.' '-'}),bse_lbl(1,:),'uni',0) , cellfun(@(x) mmil_spec_char(x,{'.' '-'}),msk_hld(:,2),'uni',0));
            
            if isempty(strfind(tot_roi{iTR},'nzvoxels'))
                out_msk{iS,iF+1} = nanmean(cell2mat(bse_lbl(iS,dat_loc)));
            else
                out_msk{iS,iF+1} = sum(cell2mat(bse_lbl(iS,dat_loc)));
            end
            
        end
    end
    
    cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' regexprep(tot_roi{iTR},rep_lbl,'Alicia')],out_msk);
    
end

end