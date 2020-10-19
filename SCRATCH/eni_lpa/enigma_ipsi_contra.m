function [ lat_dta_sbj_nme, lat_dta_dta, lat_dta_roi_nme, lat_cov_sbj_nme, lat_cov_dta, lat_cov_dta_lbl] = enigma_ipsi_contra( dta_sbj_nme, dta_dta, dta_roi_nme, cov_sbj_nme, cov_dta, cov_dta_lbl)

%% Determine relevant subjects
% Subject Clean
use_ind = strcmpi(cov_dta_lbl(:,2),'L-TLE') | strcmpi(cov_dta_lbl(:,2),'R-TLE');

lat_dta_sbj_nme = dta_sbj_nme(use_ind,:);
lat_cov_sbj_nme = cov_sbj_nme(use_ind,:);
lat_cov_dta     = cov_dta(use_ind,:);
lat_cov_dta_lbl = cov_dta_lbl(use_ind,:);

dta_hld = dta_dta(use_ind, :);

% Column Clean
roi_nme_hld = dta_roi_nme([ string_find(dta_roi_nme(1,:),{'-L'}) string_find(dta_roi_nme(1,:),{'-R'}) ]);
roi_nme_hld = cellfun(@(x) x(1:end-2), roi_nme_hld, 'uni', 0);
roi_nme_hld = unique( roi_nme_hld );

%% identify ROIs
lat_dta_dta = nan( sum(use_ind), numel(roi_nme_hld)*2 );

for iRN = 1:numel(roi_nme_hld)
    
    lat_dta_roi_nme{((iRN-1)*2) + 1} = [ roi_nme_hld{iRN} '_ipsi'  ];
    lat_dta_roi_nme{((iRN-1)*2) + 2} = [ roi_nme_hld{iRN} '_contra'];
    
    lft_ind = strcmpi( dta_roi_nme, [ roi_nme_hld{iRN} '-L' ]);
    rgh_ind = strcmpi( dta_roi_nme, [ roi_nme_hld{iRN} '-R' ]);
    
    for iS = 1:size(lat_dta_dta,1)
        if strcmpi(lat_cov_dta_lbl(iS,2),'L-TLE')
            lat_dta_dta(iS, ((iRN-1)*2) + 1 ) = dta_hld(iS,lft_ind);
            lat_dta_dta(iS, ((iRN-1)*2) + 2 ) = dta_hld(iS,rgh_ind);
        elseif strcmpi(lat_cov_dta_lbl(iS,2),'R-TLE')
            lat_dta_dta(iS, ((iRN-1)*2) + 1 ) = dta_hld(iS,rgh_ind);
            lat_dta_dta(iS, ((iRN-1)*2) + 2 ) = dta_hld(iS,lft_ind);
        end        
    end
    
end

end