clear; clc;

%% Constants
ovr_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Developement/cse_fnx_cog/';

dta_dir = [ ovr_dir '/' 'Data' ];

plt_dir = [ ovr_dir '/' 'Figures' ];

sbj_typ = { 'HC' 'LTLE' 'RTLE' };
dem_col = { 'sbj_age' 'sbj_sex' 'sbj_edu' };

scr_typ         = { 'raw'                  'nor' };
fnx_scr_typ_col = { { 'T1_raw' 'T2_raw' }  { 'T1_tscore' 'T2_tscore' } };

plt_ord = { 'fornix_one' 'fornix_two' 'LTLE' 'RTLE' 'HC' };
plt_col = { rgb('dark teal') rgb('light teal') rgb('blue') rgb('red') rgb('dark grey') };

%% Load data & Collate
% Fornix NP data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'NP_Scores_2022_11_10_updated.csv' ];
fcfg.dta_col = 2;
[ fnx_dta, fnx_dta_sbj, fnx_dta_col] = ejk_dta_frm( fcfg );
    fnx_dta = [ fnx_dta_sbj fnx_dta ];
    fnx_dta_col = [ 'test' fnx_dta_col ];

fnx_dta(cellfun(@isempty,fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1'))),strcmpi(fnx_dta_col,'Subtest_T1')) = fnx_dta(cellfun(@isempty,fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1'))),strcmpi(fnx_dta_col,'Subtest_T2'));

ovr_col = strcat(fnx_dta( :, strcmpi(fnx_dta_col,'Test')),'_',fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1')));

% Redcap
fcfg = [];
fcfg.red_fle = [ dta_dir '/' 'NP_Redcap_2023_03_17.csv' ];
fcfg.sep     = '|';
[ sbj_dem , sbj_sze , ~ , sbj_cog, ~] = ejk_load_redcap(fcfg);

%% Cipher
% Save out fornix half of cipher
% cell2csv([ dta_dir '/' 'cipher_fnx.csv' ],ovr_col)

% Load modified cipher
neu_psy_cph = mmil_readtext([ dta_dir '/' 'cipher_fnx_ejk.csv' ]);

%% Collate Demographic Data
dem_dta = cell( numel(sbj_dem.sbj_nme), numel(dem_col) );
dem_sbj = sbj_dem.sbj_nme;
for iC = 1:numel(dem_col)
    if isnumeric(sbj_dem.(dem_col{iC}))
        dem_dta(:,iC) = num2cell(sbj_dem.(dem_col{iC})); 
    else
        dem_dta(:,iC) = sbj_dem.(dem_col{iC});
    end
end

% Add Group type
sde_col = numel(dem_col)+1;
dem_col = [ dem_col 'sbj_sde_ons' ];
dem_dta(:,sde_col) = sbj_sze.sbj_sde_ons;

% Add Group type
typ_col = numel(dem_col)+1;
dem_col = [ dem_col 'sbj_typ' ];
dem_dta(:,typ_col) = repmat({''},numel(sbj_dem.sbj_nme),1);
dem_dta(intersect(string_find(dem_sbj,{'epd'}),string_find(dem_dta(:,sde_col),'L')),typ_col) = {'LTLE'};
dem_dta(intersect(string_find(dem_sbj,{'epd'}),string_find(dem_dta(:,sde_col),'R')),typ_col) = {'RTLE'};
dem_dta(string_find(dem_sbj,{'fc'}),typ_col) = {'HC'};

%% Check Education/Age/Sex Matching
% Identify candidates
clear mle_ind_mss fml_ind_mss

sbj_edu = 20;
sbj_age = 31;

for iT = 1:numel(sbj_typ)
       
    mle_ind = find(strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_sex')),'M') & strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_typ')),sbj_typ{iT}));
    fml_ind = find(strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_sex')),'F') & strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_typ')),sbj_typ{iT}));
    
    mle_edu_mss = (abs(sbj_edu - cell2mat(dem_dta(mle_ind,strcmpi(dem_col,'sbj_edu')))) ./ max(cell2mat(dem_dta(:,strcmpi(dem_col,'sbj_edu')))));
    mle_age_mss =  abs(sbj_age - cell2mat(dem_dta(mle_ind,strcmpi(dem_col,'sbj_age')))) ./ max(cell2mat(dem_dta(:,strcmpi(dem_col,'sbj_age')))); % * .66;
    [ ~, mle_ind_mss{iT} ] = sort(mle_age_mss + mle_edu_mss);
    mle_ind_mss{iT} = mle_ind(mle_ind_mss{iT}(1:10));
    
    fml_edu_mss = (abs(sbj_edu - cell2mat(dem_dta(fml_ind,strcmpi(dem_col,'sbj_edu')))) ./ max(cell2mat(dem_dta(:,strcmpi(dem_col,'sbj_edu')))));
    fml_age_mss =  abs(sbj_age - cell2mat(dem_dta(fml_ind,strcmpi(dem_col,'sbj_age')))) ./ max(cell2mat(dem_dta(:,strcmpi(dem_col,'sbj_age')))); % * .66;
    [ ~, fml_ind_mss{iT} ] = sort(fml_age_mss + fml_edu_mss);
    fml_ind_mss{iT} = fml_ind(fml_ind_mss{iT}(1:10));
    
end


% Plot to examine
for iT = 1:numel(sbj_typ)
       
    mle_ind = strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_sex')),'M') & strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_typ')),sbj_typ{iT});
    fml_ind = strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_sex')),'F') & strcmpi(dem_dta(:,strcmpi(dem_col,'sbj_typ')),sbj_typ{iT});
    
    fcfg = [];
    
    fcfg.xdt = { cell2mat(dem_dta(mle_ind,strcmpi(dem_col,'sbj_age'))) cell2mat(dem_dta(fml_ind,strcmpi(dem_col,'sbj_age'))) cell2mat(dem_dta(mle_ind_mss{iT},strcmpi(dem_col,'sbj_age'))) cell2mat(dem_dta(fml_ind_mss{iT},strcmpi(dem_col,'sbj_age')))};
    fcfg.ydt = { cell2mat(dem_dta(mle_ind,strcmpi(dem_col,'sbj_edu'))) cell2mat(dem_dta(fml_ind,strcmpi(dem_col,'sbj_edu'))) cell2mat(dem_dta(mle_ind_mss{iT},strcmpi(dem_col,'sbj_edu'))) cell2mat(dem_dta(fml_ind_mss{iT},strcmpi(dem_col,'sbj_edu')))};
    
    fcfg.fce_col     = { rgb('dark blue')   rgb('dark purple')  rgb('neon blue')   rgb('neon purple')};
    fcfg.edg_col     = { [0 0 0]             [0 0 0]            [0 0 0]            [0 0 0]};
    fcfg.box_plt_col = { rgb('light blue')  rgb('light purple') rgb('light blue')  rgb('light purple')};
    
    fcfg.xlb = {'Age'};
    fcfg.ylb = {'Education'};
    fcfg.xlm = [ 0 80 ];
    fcfg.ylm = [ 0 22 ];
    
    fcfg.hln = 20;
    fcfg.vln = 31;
    
    fcfg.ttl = [ sbj_typ{iT} ];
    
    fcfg.out_dir = plt_dir;
    fcfg.out_nme = [ 'Matching' '_' sbj_typ{iT} ];
    
    ejk_scatter(fcfg)
end

%% Collate Neuropsychological Data
% Choose tests of interest
tst_ind = find(~cellfun(@isempty,neu_psy_cph(:,2)));
cvt_ind = cell2mat(neu_psy_cph(tst_ind,4));

% Put together
for iST = 1:numel(scr_typ)
    
    cog_sbj.(scr_typ{iST}) = sbj_cog.sbj_nme;
    cog_sbj.(scr_typ{iST}) = [ 'fornix_one' ; 'fornix_two' ; cog_sbj.(scr_typ{iST}) ];
    
    cog_dta.(scr_typ{iST}) = cell( numel(cog_sbj.(scr_typ{iST})), numel(tst_ind) );
    cog_col.(scr_typ{iST}) = strrep(neu_psy_cph(tst_ind,2),'_raw_scr','');
    
    for iC = 1:numel(cog_col.(scr_typ{iST}))
        
        fnx_one_dta = fnx_dta( tst_ind(iC) , strcmpi( fnx_dta_col, fnx_scr_typ_col{iST}{1}) );
        fnx_two_dta = fnx_dta( tst_ind(iC) , strcmpi( fnx_dta_col, fnx_scr_typ_col{iST}{2}) );
        if cvt_ind(iC)==1
            fcfg = [];
            fcfg.typ = 'ss_to_t';
            red_cap_dta = num2cell(ejk_convert_neuropsych( fcfg, sbj_cog.(neu_psy_cph{tst_ind(iC),iST+1}) ));
        elseif cvt_ind(iC)==2
            fcfg = [];
            fcfg.typ = 'z_to_t';
            red_cap_dta = num2cell(ejk_convert_neuropsych( fcfg, sbj_cog.(neu_psy_cph{tst_ind(iC),iST+1}) ));
        else
            red_cap_dta = num2cell(sbj_cog.(neu_psy_cph{tst_ind(iC),iST+1}));
        end
        cog_dta.(scr_typ{iST})(:,iC) = [ fnx_one_dta ; fnx_two_dta ; red_cap_dta];
        
    end
end

% Modify demographic data to include patient
dem_sbj = [ 'fornix_one' ; 'fornix_two' ; dem_sbj];
dem_dta = [ { 31 'F' 20 '' 'fornix_one' } ; { 31 'F' 20 '' 'fornix_two'}  ; dem_dta ];

%% Make basic plots
out_plt = [ plt_dir '/' 'Cognitive' '/' ];

% Plot
for iST = 1:numel(scr_typ)
     
    out_plt_use = [ out_plt '/' scr_typ{iST} ];
    out_plt_mtc = [ out_plt '/' scr_typ{iST} '_' 'matched' ];
    ejk_chk_dir(out_plt_use);
    ejk_chk_dir(out_plt_mtc);
           
    for iT = 1:numel(cog_col.(scr_typ{iST}))
     
        % Plot all
        fcfg = [];
        
        for iPO = 1:numel(plt_ord)
            fcfg.ydt{iPO} = cell2mat(cog_dta.(scr_typ{iST})(strcmpi(dem_dta(:,typ_col),plt_ord{iPO}),iT));         
        end
        fcfg.ydt( cellfun(@isempty,fcfg.ydt) ) = {NaN};
        
        fcfg.xdt = num2cell(1:numel(fcfg.ydt));
        
        fcfg.fce_col     = plt_col;
        fcfg.edg_col     = repmat({[0 0 0]}, 1, numel(fcfg.xdt));
       
        fcfg.box_plt = [ 0 0 ones(1,numel(fcfg.xdt(3:end))) ];
        fcfg.box_plt_col = plt_col;        
        
        fcfg.xlb = plt_ord;
        fcfg.ylb = cog_col.(scr_typ{iST})(iT);
        fcfg.xlm = [ 0.5 numel(fcfg.ydt)+0.5 ];
        
        fcfg.ttl = cog_col.(scr_typ{iST}){iT};
        
        fcfg.out_dir = out_plt_use;
        fcfg.out_nme = [ scr_typ{iST} '_' cog_col.(scr_typ{iST}){iT} ];
        
        ejk_scatter(fcfg)
        
        %
        fcfg = [];
        
        cnt = 1;
        for iPO = 1:numel(plt_ord)
            if ~any(ismember(sbj_typ,plt_ord{iPO}))
                fcfg.ydt{cnt} = cell2mat(cog_dta.(scr_typ{iST})(strcmpi(dem_dta(:,typ_col),plt_ord{iPO}),iT));
                fcfg.xdt{cnt} = iPO;
                fcfg.box_plt(cnt) = 0;
                fcfg.fce_col{cnt} = plt_col{iPO};
                fcfg.box_plt_col{cnt} = plt_col{iPO};
                fcfg.xlb{cnt} = plt_ord{iPO};
                cnt = cnt + 1;
            else
                fcfg.ydt{cnt} = cell2mat(cog_dta.(scr_typ{iST})(fml_ind_mss{ismember(sbj_typ,plt_ord{iPO})}+2,iT));
                fcfg.xdt{cnt} = iPO-.25;
                fcfg.xlb{cnt} = [ plt_ord{iPO} ' ' 'F'];
                cnt = cnt + 1;
                
                fcfg.ydt{cnt} = cell2mat(cog_dta.(scr_typ{iST})(mle_ind_mss{ismember(sbj_typ,plt_ord{iPO})}+2,iT));
                fcfg.xdt{cnt} = iPO+.25;
                fcfg.xlb{cnt} = [ plt_ord{iPO} ' ' 'M'];
                
                fcfg.box_plt(cnt-1:cnt) = 1;
                fcfg.fce_col(cnt-1:cnt) = [ plt_col(iPO) plt_col(iPO) ];
                fcfg.box_plt_col(cnt-1:cnt) = [ plt_col(iPO) plt_col(iPO) ];
                
                cnt = cnt + 1;
                
                
            end
        end
        fcfg.ydt( cellfun(@isempty,fcfg.ydt) ) = {NaN};
        
        fcfg.edg_col     = repmat({[0 0 0]}, 1, numel(fcfg.xdt));
        
        fcfg.jtr_wdt = 0.10;
        fcfg.box_wdt = 0.20;
        
        fcfg.ylb = cog_col.(scr_typ{iST})(iT);
        fcfg.xlm = [ 0.5 max(cell2mat(fcfg.xdt))+0.5 ];
        
        fcfg.ttl = cog_col.(scr_typ{iST}){iT};
        
        fcfg.out_dir = out_plt_mtc;
        fcfg.out_nme = [ scr_typ{iST} '_' cog_col.(scr_typ{iST}){iT} '_matched' ];
        
        ejk_scatter(fcfg)
        
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normative scores plotting - Original
% nor_plt = { 'Age_Sex' 'Age_Sex_Education' };
% nor_nme = { { 'T1 (Age,Sex)' 'T2 (Age,Sex)' } {} };
% nor_col = { { 'T1_tscore' 'T2_tscore' }       {} };
% 
% typ_nme = { 'List_Encoding' 'List_Retrieval' };
% typ_scr     = { { 'CVLT-II_List A Trial 1' 'CVLT-II_List A Trial 2' 'CVLT-II_List A Trial 3' 'CVLT-II_List A Trial 4' 'CVLT-II_List A Trial 5' 'CVLT-II_List A 1-5 Total' 'CVLT-II_List B' } ...
%                 { 'CVLT-II_Short Delay Free Recall' 'CVLT-II_Short Delay Cued Recall' 'CVLT-II_Long Delay Free Recal' 'CVLT-II_Short Delay Free Recall'  } };
% typ_scr_nme = { { 'Trial 1'                'Trial 2'                'Trial 3'                'Trial 4'                'Trial 5'                'Trial Total'              'List B'         } ...
%                 { 'CVLT_SDFR'                       'CVLT_SDCR'                       'CVLT_LDFR'                     'CVLT_LDCR'                        } };
% 
% % Norm Plotting %%%%%%%%%%%%%%%%%%%%%%%%            
% iNO = 1;
% iTY = 2;
% 
% xdt = num2cell([[1:numel(typ_scr_nme{iTY})] [1:numel(typ_scr_nme{iTY})]]);
% 
% 
% [~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
% ydt = [[(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{1}))-50)./10]' [(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{2}))-50)./10]'];
% 
% fcfg = [];
% 
% fcfg.xdt     = xdt;
% fcfg.ydt     = ydt;
% 
% fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY}))];
% fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
% 
% fcfg.xlb = typ_scr_nme{iTY};
% fcfg.ylb = { 'z-score' };
% fcfg.ttl = typ_nme{iTY};
% 
% fcfg.hln     = -1;
% fcfg.hln_col = rgb('red');
% 
% fcfg.jtr = 1;
% fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+0.5 ];
% fcfg.ylm = [ -4 4 ];
% 
% fcfg.out_dir = [ ovr_dir '/' 'plots' '/' nor_plt{iNO} '/' ];
% fcfg.out_nme = typ_nme{iTY};
% 
% ejk_scatter(fcfg)

%% Normative Data Group Plotting
% typ_nme = { 'CVLT' };
% typ_scr     = { { 'CVLT-II_List A 1-5 Total' 'CVLT-II_Long Delay Free Recal'   } };
% typ_scr_nme = { { 'CVLT_Total'               'CVLT_LDFR'           } };
% 
% raw_col = { 'T1_raw' 'T2_raw' };
% nor_plt = { 'Age_Sex' 'Age_Sex_Education' };
% nor_col = { { 'T1_tscore' 'T2_tscore' }       {} };
% 
% cph_typ = mmil_readtext([ ovr_dir '/' 'data' '/' 'cipher_bth.csv' ]);
% 
% % Norm Group Plotting %%%%%%%%%%%%%%%%%%%%%%%%
% iNO = 1;
% iTY = 1;
% 
% xdt = num2cell([[1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]+.1 ]);
% 
% [~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
% ydt = [[(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{1}))-50)./10]' [(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{2}))-50)./10]'];
% xlb = [];
% for iRA = 1:numel(typ_scr{iTY})
%     ydt = [ ydt (sbj_cog.(cph_typ{strcmpi(cph_typ(:,1),typ_scr{iTY}{iRA}),3})(string_find(sbj_cog.sbj_nme,'fc'))-50)./10 ];
%     xlb = [ xlb typ_scr_nme{iTY}(iRA) {'HC'} ];
% end
% 
% fcfg = [];
% 
% fcfg.xdt     = xdt;
% fcfg.ydt     = ydt;
% 
% fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY})) repmat({rgb('grey')},1,numel(typ_scr_nme{iTY}))];
% fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
% 
% fcfg.xlb = xlb;
% fcfg.ylb = { 'Z-Score' };
% fcfg.ttl = typ_nme{iTY};
% 
% fcfg.hln     = -1;
% fcfg.hln_col = rgb('red');
% 
% fcfg.jtr = 1;
% fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+0.5 ];
% fcfg.ylm = [ -4 4 ];
% 
% fcfg.out_dir = [ ovr_dir '/' 'plots' '/' nor_plt{iNO} '/' ];
% fcfg.out_nme = [ typ_nme{iTY} '_' 'group' ];
% 
% ejk_scatter(fcfg)
% 
% % Raw Group Plotting %%%%%%%%%%%%%%%%%%%%%%%%
% raw_col = { 'T1_raw' 'T2_raw' };
% 
% cph_typ = mmil_readtext([ ovr_dir '/' 'data' '/' 'cipher_bth.csv' ]);
% 
% 
% iTY = 1;
% 
% xdt = num2cell([[1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]+.1 ]);
% 
% [~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
% ydt = [[fnx_dta(ydt_row,strcmpi(fnx_dta_col,raw_col{1}))]' [fnx_dta(ydt_row,strcmpi(fnx_dta_col,raw_col{2}))]'];
% xlb = [];
% for iRA = 1:numel(typ_scr{iTY})
%     ydt = [ ydt sbj_cog.(cph_typ{strcmpi(cph_typ(:,1),typ_scr{iTY}{iRA}),2})(string_find(sbj_cog.sbj_nme,'fc')) ];
%     xlb = [ xlb typ_scr_nme{iTY}(iRA) {'HC'} ];
% end
% 
% fcfg = [];
% 
% fcfg.xdt     = xdt;
% fcfg.ydt     = ydt;
% 
% fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY})) repmat({rgb('grey')},1,numel(typ_scr_nme{iTY}))];
% fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));
% 
% fcfg.xlb = xlb;
% fcfg.ylb = { 'Raw Score' };
% fcfg.ttl = typ_nme{iTY};
% 
% fcfg.jtr = 1;
% fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+1 ];
% 
% fcfg.out_dir = [ ovr_dir '/' 'plots' '/' 'Raw' '/' ];
% fcfg.out_nme = typ_nme{iTY};
% 
% ejk_scatter(fcfg)








