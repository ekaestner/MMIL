function mmil_save_plot(cfg)

if isfield(cfg,'sve_plt')
    savefig(cfg.pth);
    fid = fopen('/home/ekaestne/matlab/bespoke/misc/saved_figures.csv','a');
    fprintf(fid,['\n' cfg.pth]);
    fclose(fid);
elseif isfield(cfg,'prt_plt')
    plt = mmil_readtext('/home/ekaestne/matlab/bespoke/misc/saved_figures.csv');
    for iS = 1:numel(plt)
        if exist(plt{iS},'file')
        hgload(plt{iS}); 
        fprintf('Saving plot %i out of %i : %s \n',iS,numel(plt),plt{iS})
        print(gcf,[plt{iS}(1:end-4) '.png'],'-dpng','-r500');
        delete(plt{iS})
        close all
        end         
    end
    plt = {};
    cell2csv('/home/ekaestne/matlab/bespoke/misc/saved_figures.csv',plt)
end