
cor_hld_kep = nan(1,size(sal_dta_nor.zsc_dta.keep.x3D.ep,2));
cor_hld_rmv = nan(1,size(sal_dta_nor.zsc_dta.keep.x3D.ep,2));
for iCO = 1:size(sal_dta_nor.zsc_dta.keep.x3D.ep,2)    
    dta_hld = mean(squeeze(sal_dta_nor.zsc_dta.keep.x3D.ep(:,iCO,:,:)),3);
    cor_hld_kep(iCO) = nanmean(dta_hld(:));
    dta_hld = mean(squeeze(sal_dta_nor.zsc_dta.remove.x3D.ep(:,iCO,:,:)),3);
    cor_hld_rmv(iCO) = nanmean(dta_hld(:));
end

plot(1:numel(cor_hld),cor_hld_kep,'k'); hold on;
plot(1:numel(cor_hld),cor_hld_rmv,'g');