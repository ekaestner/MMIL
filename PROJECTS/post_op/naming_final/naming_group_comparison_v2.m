out_put = [ prj_dir '/' prj_nme];

%% Make groups
sbj_hld = { { 'epd' } { 'fc' 'epd' } {'fc'} };
sbj_nme = { 'tle'     'tle_controls' 'controls' };

tst_hld = { [ 1 2 ] [ 3 4 ] };
tst_nme = { 'pre'   'post'};

str_hld = { [1.5 3] [3] };
str_nme = { 'allT' '3T' };

srg_hld = { ''        'ATL' };
srg_nme = { 'allSurg' 'ATLonly' };

% sbj_typ x tst_typ x str_typ x srg_typ x grp_typ
run_typ = { [ 1 1 2 1 1 ] ... % epd,          pre,  3T, Non-Surgical, all patients 
            [ 1 1 2 1 2 ] ... % epd,          pre,  3T, Non-Surgical, L/R split 
            [ 1 1 2 2 1 ] ... % epd,          pre,  3T, ATL,          all patients 
            [ 1 1 2 2 2 ] ... % epd,          pre,  3T, ATL,          L/R split
            [ 1 2 2 2 1 ] ... % epd,          post, 3T, ATL,          all patients 
            [ 1 2 2 2 2 ] ... % epd,          post, 3T, ATL,          L/R split
            [ 2 1 2 1 1 ] ... % controls/epd, pre,  3T, Non-Surgical, all patients 
            [ 3 1 2 1 1 ] ... % controls,     pre,  3T, Non-Surgical, all  
            [ 2 1 2 1 3 ] ... % controls/epd, pre,  3T, Non-Surgical, L/R split including controls
            }; 

cln_dta = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ]);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
cog_dta_col = cog_dta(1,2:end);
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);
        
%% Setup Groups
clear grp

for iG = 1:numel(run_typ)
    
    sbj_typ = run_typ{iG}(1);
    tst_typ = run_typ{iG}(2);
    str_typ = run_typ{iG}(3);
    srg_typ = run_typ{iG}(4);
    grp_typ = run_typ{iG}(5);
    
    % Choose subjects to use %%%%%%%%
    % Find subject type
    tot_sbj = zeros(0);
    for iT = 1:numel(sbj_hld{sbj_typ})
        tot_sbj = [ tot_sbj ; string_find(cln_dta_sbj(:,1),sbj_hld{sbj_typ}{iT}) ];
    end
    tot_sbj = unique(tot_sbj);
    
    % Find patients with cognitive tests
    cog_sbj = zeros(0);
    for iT = 1:numel(tst_hld{tst_typ})
        cog_sbj = [ cog_sbj ; find(~isnan(cell2mat(cog_dta(:,tst_hld{tst_typ}(iT))))) ];
    end
    cog_sbj = unique(cog_sbj);
    
    % Find Strength
    str_sbj = zeros(0);
    for iT = 1:numel(str_hld{str_typ})
        str_sbj = [ str_sbj ; find(strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), num2str(str_hld{str_typ}(iT)))) ];
    end
    str_sbj = unique(str_sbj);
    
    % Find Surgery
    srg_sbj = zeros(0);
    if ~strcmpi(srg_hld{srg_typ},'')
        for iT = 1:numel(srg_hld{srg_typ})
            srg_sbj = [ srg_sbj ; find(strcmpi( cln_dta(:,10), srg_hld{srg_typ})) ];
        end
        srg_sbj = unique(srg_sbj);
    else
        srg_sbj = 1:size(cln_dta,1);
    end
    
    % Type of Group
    % Grab
    lft_tle = find( strcmpi(cln_dta(:,2),'L') );
    rgh_tle = find( strcmpi(cln_dta(:,2),'R') );
    
    if grp_typ==1
        grp.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'all'] ) = intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj);
        cog_col.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'all'] ) = tst_hld{tst_typ};
    elseif grp_typ==2
        grp.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'left'] )  = intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), lft_tle);
        grp.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'right'] ) = intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), rgh_tle);
        cog_col.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'left'] ) = tst_hld{tst_typ};
        cog_col.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'right'] ) = tst_hld{tst_typ};
    elseif grp_typ==3
        grp.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'left'] )  = ...
            [ intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), lft_tle) ; ...
              grp.controls_pre_3T_allSurg_all ] ;
        grp.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'right'] ) = ...
            [ intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), rgh_tle) ; ...
              grp.controls_pre_3T_allSurg_all ];
        cog_col.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'left'] ) = tst_hld{tst_typ};
        cog_col.( [ sbj_nme{sbj_typ} '_' tst_nme{tst_typ} '_' str_nme{str_typ} '_' srg_nme{srg_typ} '_' 'right'] ) = tst_hld{tst_typ};
    end
        
end

%% Save
save([ out_put '/' 'groups.mat' ], 'grp', 'cog_col');






