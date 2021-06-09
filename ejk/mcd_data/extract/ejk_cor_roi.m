function [ dta_out ] = ejk_cor_roi( cfg )

cor_col_ind = find(strcmpi( cfg.dta(1,:) , cfg.cor_col ));

%%
dta_out = cell( size(cfg.dta) );

dta_out(:,1:4) = cfg.dta(:,1:4);
dta_out(1,:) = cfg.dta(1,:);

for iC = 5:size(dta_out,2)
    for iR = 2:size(dta_out,1)
        
        dta_out{iR,iC} = (cfg.dta{iR,iC} / cfg.dta{iR,cor_col_ind}) * 100;
        
    end
end


end