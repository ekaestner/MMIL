%% ft_plot_spectrum
% 
% cfg.outdir = string defining directory to output data
% cfg.prefix = name of file to save
% 
% Created 9/29/14 by Erik Kaestner (ekaestne@ucsd.edu)

function ft_plot_spectrum(cfg,ft_dat)

if ~isdir(cfg.outdir); mkdir(cfg.outdir); end

num_chn = numel(ft_dat.label);
tic
for iCH = 1:ceil(num_chn/10)
    iPL = 1; figure('Visible','off');
    while iPL <= 10 && (iPL + (iCH-1)*10 < num_chn)
        
        spc_hld = mmil_welch(ft_dat.trial{1}(iPL + (iCH-1)*10,:),[],[],[],ft_dat.fsample);
        spc_num = zeros(numel(ft_dat.trial),size(spc_hld.Data,1));
        
        trl_ind = randperm(numel(ft_dat.trial),min(40,numel(ft_dat.trial)));
        
        for iT = 1:numel(trl_ind); 
            spc_hld = mmil_welch(ft_dat.trial{trl_ind(iT)}(iPL + (iCH-1)*10,:),[],[],[],ft_dat.fsample);
            spc_num(iT,:) = spc_hld.Data;
        end  
        wlc = dspdata.psd(mean(spc_num),'Fs',spc_hld.Fs,'SpectrumType',spc_hld.SpectrumType);
        
        subplot(5,2,iPL);
        plot(wlc);
        ylim([-75 75]);
        xlabel('');ylabel('');title([num2str(iPL + (iCH-1)*10) ' ' ft_dat.label{iPL + (iCH-1)*10}]);
        iPL = iPL + 1;
    end
    print(gcf,[cfg.outdir '/' cfg.prefix '_' num2str(iCH)],'-dpng'); close all;
end
toc

end
