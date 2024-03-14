clear; clc;

dta_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/05_Developement/wda';

%%
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Wada_Raw_Data.csv' ];
[ wda_dta, wda_sbj, wda_col ]  = ejk_dta_frm(fcfg);

%%
sde_col = strcmpi(wda_col,'Side of Injection');
mem_col = strcmpi(wda_col,'Memory Type');
dom_col = strcmpi(wda_col,'Domain');
prf_col = strcmpi(wda_col,'Response Type');
foi_col = strcmpi(wda_col,'Foil');

wda_dta(:,mem_col) = mmil_spec_char(wda_dta(:,mem_col),{' '},{'_'});
wda_dta(:,dom_col) = mmil_spec_char(wda_dta(:,dom_col),{' '},{'_'});

tot_sbj = unique(wda_sbj);
tot_sde = unique(wda_dta(:,sde_col));
tot_mem = unique(wda_dta(:,mem_col));
tot_dom = unique(wda_dta(:,dom_col));

wda_prf_col = [ strcat(tot_mem,'_',tot_sde{1})' strcat(tot_dom,'_',tot_sde{1})' strcat('total','_',tot_sde{1}) strcat(tot_mem,'_',tot_sde{2})' strcat(tot_dom,'_',tot_sde{2})' strcat('total','_',tot_sde{2})];
wda_prf.hit = nan(numel(tot_sbj),numel(wda_prf_col));
wda_prf.fls = nan(numel(tot_sbj),numel(wda_prf_col));
wda_prf.dpr = nan(numel(tot_sbj),numel(wda_prf_col));

%%
foi_ind = cell2mat(wda_dta(:,foi_col))==1;
tgt_ind = cell2mat(wda_dta(:,foi_col))==0;

for iS = 1:numel(tot_sbj)
    sbj_ind = strcmpi(wda_sbj,tot_sbj{iS});
    
    for iH = 1:numel(tot_sde)
        hms_ind = strcmpi(wda_dta(:,sde_col), tot_sde{iH});
        
        % Memory Types
        for iM = 1:numel(tot_mem)           
            mem_ind = strcmpi(wda_dta(:,mem_col), tot_mem{iM});
            
            col_ind = strcmpi(wda_prf_col,[tot_mem{iM} '_' tot_sde{iH}]);
            
            wda_prf.hit(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & mem_ind & tgt_ind,prf_col))==1) / numel(wda_dta( sbj_ind & hms_ind & mem_ind & tgt_ind,prf_col))) * 100;
            wda_prf.fls(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & mem_ind & foi_ind,prf_col))==3) / numel(wda_dta( sbj_ind & hms_ind & mem_ind & foi_ind,prf_col))) * 100;
            
            if wda_prf.hit(iS,col_ind)==100; hit_cor = 99; elseif wda_prf.hit(iS,col_ind)==0; hit_cor = 1;  else hit_cor = wda_prf.hit(iS,col_ind);end
            if wda_prf.fls(iS,col_ind)==100; fls_cor = 99; elseif wda_prf.fls(iS,col_ind)==0; fls_cor = 1;  else fls_cor = wda_prf.fls(iS,col_ind);end
            
            wda_prf.dpr(iS,col_ind) = norminv(hit_cor/100)-norminv(fls_cor/100);
        end
        
        % Domains
        for iD = 1:numel(tot_dom)           
            dom_ind = strcmpi(wda_dta(:,dom_col), tot_dom{iD});
            
            col_ind = strcmpi(wda_prf_col,[ tot_dom{iD} '_' tot_sde{iH}]);
            
            wda_prf.hit(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & dom_ind & tgt_ind,prf_col))==1) / numel(wda_dta( sbj_ind & hms_ind & dom_ind & tgt_ind,prf_col))) * 100;
            wda_prf.fls(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & dom_ind & foi_ind,prf_col))==3) / numel(wda_dta( sbj_ind & hms_ind & dom_ind & foi_ind,prf_col))) * 100;
            
            if wda_prf.hit(iS,col_ind)==100; hit_cor = 99; elseif wda_prf.hit(iS,col_ind)==0; hit_cor = 1;  else hit_cor = wda_prf.hit(iS,col_ind);end
            if wda_prf.fls(iS,col_ind)==100; fls_cor = 99; elseif wda_prf.fls(iS,col_ind)==0; fls_cor = 1;  else fls_cor = wda_prf.fls(iS,col_ind);end
            
            wda_prf.dpr(iS,col_ind) = norminv(hit_cor/100)-norminv(fls_cor/100);
        end
        
        % Total        
        col_ind = strcmpi(wda_prf_col,['total' '_' tot_sde{iH}]);
        
        wda_prf.hit(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & tgt_ind,prf_col))==1) / numel(wda_dta( sbj_ind & hms_ind & tgt_ind,prf_col))) * 100;
        wda_prf.fls(iS,col_ind) = (sum(cell2mat(wda_dta( sbj_ind & hms_ind & foi_ind,prf_col))==3) / numel(wda_dta( sbj_ind & hms_ind & foi_ind,prf_col))) * 100;
        
        if wda_prf.hit(iS,col_ind)==100; hit_cor = 99; elseif wda_prf.hit(iS,col_ind)==0; hit_cor = 1;  else hit_cor = wda_prf.hit(iS,col_ind);end
        if wda_prf.fls(iS,col_ind)==100; fls_cor = 99; elseif wda_prf.fls(iS,col_ind)==0; fls_cor = 1;  else fls_cor = wda_prf.fls(iS,col_ind);end
        
        wda_prf.dpr(iS,col_ind) = norminv(hit_cor/100)-norminv(fls_cor/100);
        
        
    end
end

%%
plt_col = distinguishable_colors(numel(tot_sbj));

plt_tst     = { wda_prf_col([9 8 7 18 17 16]) wda_prf_col([1 10 2 11 3 12 4 13 5 14 6 15]) };
plt_tst_nme = { 'summary' 'individual'};
plt_stt     = { 'hit' 'fls' 'dpr' };

% iTS = 1;
% iST = 3;

for iTS = 1:numel(plt_tst)
    for iST = 1:numel(plt_stt)
        
        % Side by Side
        for iM = 1:numel(plt_tst{iTS})
            use_col(iM) = find(strcmpi(wda_prf_col,plt_tst{iTS}{iM}));
        end
        
        fcfg = [];
        
        for iS = 1:numel(tot_sbj)
            fcfg.ydt{iS}     = wda_prf.(plt_stt{iST})(iS,use_col);
            fcfg.xdt{iS}     = 1:numel(fcfg.ydt{iS});
            fcfg.fce_col{iS} = plt_col(iS,:);
        end
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
        
        fcfg.xlb = plt_tst{iTS};
        fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
        fcfg.ylb = plt_stt(iST);
        
        fcfg.jtr = 1;
        fcfg.jtr_wdt = 0.15;
        
        fcfg.out_dir = [ dta_dir '/' 'plots' '/' 'side_by_side' '/' ];
        fcfg.out_nme = [ plt_tst_nme{iTS} '_' plt_stt{iST} ];

        ejk_scatter(fcfg)
        
    end
end

% Subtraction Plots
for iTS = 1:numel(plt_tst)
    
    unq_tst = cellfun(@(x) strsplit(x,'_'),plt_tst{iTS},'uni',0);
    unq_tst = cellfun(@(x) strcat(x(1:end-1),'_'),unq_tst,'uni',0);
    unq_tst = unique(cellfun(@(x) x{1}(1:end-1),unq_tst,'uni',0));
    
    for iST = 1:numel(plt_stt)
        
        for iM = 1:numel(plt_tst{iTS})
            use_col(iM) = find(strcmpi(wda_prf_col,plt_tst{iTS}{iM}));
        end
        
        fcfg = [];
        
        for iS = 1:numel(tot_sbj)            
            ydt_hld = wda_prf.(plt_stt{iST})(iS,use_col);
            for iD = 1:numel(unq_tst)
                ind = string_find(plt_tst{iTS},unq_tst{iD});
                fcfg.ydt{iS}(iD) = ydt_hld(ind(1)) - ydt_hld(ind(2));
            end
            fcfg.xdt{iS}     = 1:numel(fcfg.ydt{iS});
            fcfg.fce_col{iS} = plt_col(iS,:);
        end
                
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
        
        fcfg.xlb = unq_tst;
        fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
        fcfg.ylb = {[ plt_stt{iST} ' ' '(left-right)']};
        fcfg.ylm = [-ceil(max(abs(cat(2,fcfg.ydt{:})))) ceil(max(abs(cat(2,fcfg.ydt{:}))))];
        
        fcfg.hln = 0;
        
        fcfg.jtr = 1;
        fcfg.jtr_wdt = 0.15;
        
        fcfg.out_dir = [ dta_dir '/' 'plots' '/' 'difference' '/' ];
        fcfg.out_nme = [ plt_tst_nme{iTS} '_' plt_stt{iST} ];
                
        ejk_scatter(fcfg)
        
    end
end



