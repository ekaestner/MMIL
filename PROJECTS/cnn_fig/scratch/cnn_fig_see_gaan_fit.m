

for iM = 1:numel(mdl_nme)

    ejk_chk_dir( [ mdl_fld '/' mdl_nme{iM} '/' 'fit' '/' dte_str ] );

    % Get Data
    fcfg = [];
    fcfg.dta_loc = [ mdl_fld '/' mdl_nme{iM} '/' 'results_figaan.csv' ];
    fcfg.dta_col = 2;
    [ fit_dta, fit_sbj, fit_col] = ejk_dta_frm( fcfg );
    fit_sbj = cellfun(@(x) (x+1)*100,fit_sbj);

    % Plot Omnibus Figure
    figure('visible','off')
    for iP = 1:numel(fit_col)
        subplot(4,4,iP)
        plot(fit_sbj',cell2mat(fit_dta(:,iP))')
        subtitle(mmil_spec_char(fit_col{iP},{'_'},{' '}))
    end
    tightfig()
    print([ mdl_fld '/' mdl_nme{iM} '/' 'fit' '/' dte_str '/' '00_omnibus.png' ],'-dpng');
    close all

    % Plot individual figures    
    for iP = 1:numel(fit_col)
        figure('visible','off')
        plot(fit_sbj',cell2mat(fit_dta(:,iP))')
        title(mmil_spec_char(fit_col{iP},{'_'},{' '}))
        tightfig()
        print([ mdl_fld '/' mdl_nme{iM} '/' 'fit' '/' dte_str '/' fit_col{iP} '.png' ],'-dpng');
        close all
    end

end