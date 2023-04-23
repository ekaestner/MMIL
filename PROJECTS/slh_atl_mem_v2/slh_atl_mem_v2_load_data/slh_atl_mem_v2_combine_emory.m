
%% Combine Emory
% Create combined data structures
cmb_emy_sbj = unique( [ emy_cog_sbj ; emy_dem_sbj]);
cmb_emy_col = dem_nme_row';
cmb_emy_dta = cell(numel(cmb_emy_sbj),numel(cmb_emy_col));

% Get Clinical Data subject ind
dem_sbj_ind = nan(numel(cmb_emy_sbj),2);
for iS = 1:size(dem_sbj_ind,1)
    dem_sbj_ind(iS,1) = iS;
    try dem_sbj_ind(iS,2) = find(strcmpi(emy_dem_sbj,cmb_emy_sbj{iS})); catch; end;
end
dem_sbj_ind(isnan(dem_sbj_ind(:,2)),:) = [];

% Get Cognitive Data subject in
cog_sbj_ind = nan(numel(cmb_emy_sbj),2);
for iS = 1:size(cog_sbj_ind,1)
    cog_sbj_ind(iS,1) = iS;
    try cog_sbj_ind(iS,2) = find(strcmpi(emy_cog_sbj,cmb_emy_sbj{iS})); catch; end;
end
cog_sbj_ind(isnan(cog_sbj_ind(:,2)),:) = [];

% Go through
for iR = 1:numel(cmb_emy_col)
    
    emy_cog_nme = dem_nme_dta{ strcmpi(dem_nme_row,cmb_emy_col{iR}),strcmpi(dem_nme_col,'emy_cog_col_nme') };
    emy_dem_nme = dem_nme_dta{ strcmpi(dem_nme_row,cmb_emy_col{iR}),strcmpi(dem_nme_col,'emy_cln_col_nme') };
    
    if ~isempty(emy_dem_nme)
        cmb_emy_dta(dem_sbj_ind(:,1),iR) = emy_dem_dta(dem_sbj_ind(:,2),strcmpi(emy_dem_col,emy_dem_nme));
    elseif ~isempty(emy_cog_nme)
        cmb_emy_dta(cog_sbj_ind(:,1),iR) = emy_cog_dta(cog_sbj_ind(:,2),strcmpi(emy_cog_col,emy_cog_nme));
    end    
end
cmb_emy_dta(:,strcmpi(cmb_emy_col,'sbj_nme')) = cmb_emy_sbj;

%% Combine Emory w/ Redcap
for iR = 1:numel(dem_nme_row)
    if ~isempty(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}) && ~(strcmpi(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')},'pst_cog_raw') || strcmpi(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')},'pst_cog_rci')) 
        if isnumeric(red_cap_dta.(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld_nme')})(1))
            cmb_emy_dta(cellfun(@isempty,cmb_emy_dta(:,strcmpi(cmb_emy_col,dem_nme_row{iR}))),strcmpi(cmb_emy_col,dem_nme_row{iR})) = {NaN};
            red_cap_dta.(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld_nme')}) = ...
                [ red_cap_dta.(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld_nme')}) ; 
                  cell2mat(cmb_emy_dta(:,strcmpi(cmb_emy_col,dem_nme_row{iR})))  ];  
        else
            red_cap_dta.(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld_nme')}) = ...
                [ red_cap_dta.(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iR,strcmpi(dem_nme_col,'red_cap_fld_nme')}) ; 
                  cmb_emy_dta(:,strcmpi(cmb_emy_col,dem_nme_row{iR}))  ];  
        end
    end
end

%% Recode


%% Add post-cognitive
fcfg = [];
pst_cog_dta = ejk_post_cognitive(fcfg,red_cap_dta.sbj_cog);

red_cap_dta.pst_cog_raw = pst_cog_dta.raw;
red_cap_dta.pst_cog_rci = pst_cog_dta.rci;

%% Extract
% Create final data structures
cmb_dta_sbj = red_cap_dta.(dem_nme_dta{strcmpi(dem_nme_row,'sbj_nme'),strcmpi(dem_nme_col,'red_cap_fld')}).sbj_nme;
cmb_dta_col = dem_nme_row';
cmb_dta_dta = cell(numel(cmb_dta_sbj),numel(cmb_dta_col));

%
for iC = 1:numel(cmb_dta_col)
    if ~isempty(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld')})
        if isnumeric(red_cap_dta.(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld_nme')})(1))
            cmb_dta_dta(:,iC) = num2cell(red_cap_dta.(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld_nme')}));
        else
            cmb_dta_dta(:,iC) = red_cap_dta.(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld')}).(dem_nme_dta{iC,strcmpi(dem_nme_col,'red_cap_fld_nme')});
        end        
    else
        cmb_dta_dta(:,iC) = cell(numel(cmb_dta_sbj),1);
    end
end

%% Save
cell2csv([ dta_dir '/' 'Total_Demographic_Clinical_Data.csv'], [ {'sbj_nme'} cmb_dta_col ; cmb_dta_sbj cmb_dta_dta]);

