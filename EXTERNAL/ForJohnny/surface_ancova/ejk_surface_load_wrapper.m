%% Load Surface
smt_stp = [ 176 313 705 2819 ];

prj_dir = '/home/ekaestne/PROJECTS/';

hms     = {'lhs' 'rhs'};

for iSM = [1 3 4]
    tic;
    fcfg = [];
    
    fcfg.prj_dir = prj_dir;
    
    fcfg.mes_typ = 'aMRI_thickness';
    fcfg.smt_stp  = smt_stp(iSM);
    
    fcfg.anl_dir = 'analysis';
    
    for iH = 1:numel(hms)
        
        srf_dta  = nan(size(sbj_nme,1),163842);
        fcfg.hms = hms{iH}(1:2);
        
        for iS = 1:size(sbj_nme,1)
            
            fcfg.prc_dir     = dta_loc{iS,1};
            fcfg.sbj_fsr_dir = sbj_nme{iS,2};
            
            srf_dta_hld = ejk_extract_vertices(fcfg);
            if ~isempty(srf_dta_hld); srf_dta(iS,:) = srf_dta_hld; end
            
        end
        srf_dta_sbj = sbj_nme(1:end,1);
        save([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'surf' '_' fcfg.mes_typ '_' fcfg.hms 's' '_' 'sm' num2str(fcfg.smt_stp) '.mat'],'srf_dta_sbj','srf_dta');
        clear srf_dta
        
    end
    toc
end