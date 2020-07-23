

function dta_zsc = ejk_create_zscore( cfg , dta )

dta_zsc = zeros( size(dta,1) , size(dta,2) );

%% 
for iC = 1:size(dta,2)

    %
    dta_men = nanmean( dta( cfg.con_ind , iC) );
    dta_std = nanstd( dta( cfg.con_ind , iC) );
    
    %
    if dta_std==0
        dta_zsc(:,iC) = nan( numel(dta_zsc(:,iC)) , 1 );
    else
        dta_zsc(:,iC) = ( dta( : , iC ) - dta_men ) ./ dta_std;
    end
    
end



