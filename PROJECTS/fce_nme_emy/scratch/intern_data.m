clear; clc;

dta_dir = '/home/ekaestner/Dropbox/McDonald Lab/RA_Stuff/HighTechHigh_Interns_2022/database_check';

%% Load data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Database_Solei.csv'];
[ sol_dta, sol_dta_col ] = ejk_dta_frm(fcfg);
    sol_dta_sbj = sol_dta(1,:)';
    sol_dta = sol_dta(2:end,:)';
    sol_dta_col = sol_dta_col(2:end);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Database_Bianca.csv'];
[ bnc_dta, bnc_dta_col ] = ejk_dta_frm(fcfg);
    bnc_dta_sbj = bnc_dta(1,:)';
    bnc_dta = bnc_dta(2:end,:)';
    bnc_dta_col = bnc_dta_col(2:end);


%% IRR
ovr_sbj = intersect(sol_dta_sbj,bnc_dta_sbj);    

out_mtc = { ['Bianca' '_' 'Subject']  ['Solei' '_' 'Subject']  ...
            ['Bianca' '_' 'Variable'] ['Solei' '_' 'Variable'] ...
            ['Bianca' '_' 'Value']    ['Solei' '_' 'Value'] };
        
for iS = 1:numel(ovr_sbj)
    
    bia_sbj_num = find(strcmpi(bnc_dta_sbj,ovr_sbj{iS}));
    sol_sbj_num = find(strcmpi(sol_dta_sbj,ovr_sbj{iS}));
    
    for iV = 1:numel(sol_dta_col)
                
        bia_num = isnumeric(bnc_dta{bia_sbj_num,iV});
        sol_num = isnumeric(sol_dta{sol_sbj_num,iV});
        
        if bia_num && sol_num
            bia_nan = isnan(bnc_dta{bia_sbj_num,iV});
            sol_nan = isnan(sol_dta{sol_sbj_num,iV}); 
            
            if ~(bia_nan && sol_nan) && (bnc_dta{bia_sbj_num,iV} ~= sol_dta{sol_sbj_num,iV})
                out_mtc = [ out_mtc ; bnc_dta_sbj(bia_sbj_num) sol_dta_sbj(sol_sbj_num) bnc_dta_col(iV) sol_dta_col(iV) bnc_dta(bia_sbj_num,iV) sol_dta(sol_sbj_num,iV) ];
            end
        elseif  ~bia_num && ~sol_num
            if ~strcmpi(bnc_dta{bia_sbj_num,iV},sol_dta{sol_sbj_num,iV})
               out_mtc = [ out_mtc ; bnc_dta_sbj(bia_sbj_num) sol_dta_sbj(sol_sbj_num) bnc_dta_col(iV) sol_dta_col(iV) bnc_dta(bia_sbj_num,iV) sol_dta(sol_sbj_num,iV) ]; 
            end
        elseif  bia_num ~= sol_num
            out_mtc = [ out_mtc ; bnc_dta_sbj(bia_sbj_num) sol_dta_sbj(sol_sbj_num) bnc_dta_col(iV) sol_dta_col(iV) bnc_dta(bia_sbj_num,iV) sol_dta(sol_sbj_num,iV) ];
        end        
    end    
end

cell2csv([ dta_dir '/' 'Mismatch.csv' ],out_mtc)

%% Integrate Data
tot_sbj = unique([ sol_dta_sbj ; bnc_dta_sbj ]);
tot_col = bnc_dta_col;
tot_dta = cell( numel(tot_sbj), numel(tot_col) );

for iS = 1:numel(tot_sbj)
    
    sol_sbj = ismember(tot_sbj{iS},sol_dta_sbj);
    bnc_sbj = ismember(tot_sbj{iS},bnc_dta_sbj);
    
    if bnc_sbj
        sbj_ind = strcmpi(bnc_dta_sbj,tot_sbj{iS});
        tot_dta(iS,:) = bnc_dta(sbj_ind,:);
    elseif sol_sbj
        sbj_ind = strcmpi(sol_dta_sbj,tot_sbj{iS});
        tot_dta(iS,:) = sol_dta(sbj_ind,:);
    end
end

cell2csv([dta_dir '/' 'Total_Database.csv'],[ {'sbj_nme'} tot_col' ; tot_sbj tot_dta])

%% Create & check data entry type
cat_typ_out = cell(numel(tot_col),7);
wrg_cat_out = cell(0);
for iC = 1:numel(tot_col)
    
    cat_typ_out{iC,1} = tot_col{iC};
    
    % Identify types of data present
    tot_ind = 1:size(tot_dta,1);
    emp_ind = find(cellfun(@isempty,tot_dta(:,iC)));
    mss_ind = setxor(find(cellfun(@all,cellfun(@(x) x==-9999,tot_dta(:,iC),'uni',0))),emp_ind);
    nsc_ind = setxor(find(cellfun(@all,cellfun(@isnan,tot_dta(:,iC),'uni',0))),emp_ind);
    
    tot_ind = setxor( tot_ind, unique([ emp_ind ; mss_ind ; nsc_ind ]));    
   
    % Identify types of data
    if ~isempty(tot_ind)
        cat_tbl = tabulate( cellfun(@class,tot_dta(tot_ind,iC),'uni',0) );
        [~,cat_typ_ind] = max(cell2mat(cat_tbl(:,2)));
        cat_typ_out{iC,2} = cat_tbl{cat_typ_ind,1};
    else
        cat_typ_out{iC,2} = '';
    end

    cat_typ_out{iC,3} = size(tot_dta,1);
    cat_typ_out{iC,4} = numel(tot_ind);
    cat_typ_out{iC,5} = numel(emp_ind);
    cat_typ_out{iC,6} = numel(mss_ind);
    cat_typ_out{iC,7} = numel(nsc_ind);
        
    % Identify conflicting data
    if ~isempty(tot_ind)
        for iS = 1:numel(tot_ind)
            if ~strcmpi(class(tot_dta{tot_ind(iS),iC}),cat_typ_out{iC,2})
                wrg_cat_out(end+1,1) = tot_sbj(tot_ind(iS));
                wrg_cat_out(end,2)   = tot_col(iC);
                wrg_cat_out(end,3)   = cat_typ_out(iC,2);
                wrg_cat_out(end,4)   = tot_dta(tot_ind(iS),iC);
            end
        end
    end
    
end

cell2csv([ dta_dir '/' 'Data_type.csv'],[ {'' 'Type' 'PossibleN' 'Present' 'Empty' 'NotScored' 'Missing'} ; cat_typ_out ]);
cell2csv([ dta_dir '/' 'Conflicting_type.csv'],wrg_cat_out);

%% Intern Figures
fcfg = [];
fcfg.dta_loc    = [dta_dir '/' 'Total_Database.csv'];
[tot_dta, tot_sbj, tot_col] = ejk_dta_frm(fcfg);

sde_col = string_find(tot_col,'side_seizure_onset');
edu_col = string_find(tot_col,'education');
age_ons_col = string_find(tot_col,'age_of_seizure_onset');
tme_one = string_find(tot_sbj,'-T1');

bnt_col = string_find(tot_col,'BNT_with_semantic_cues_raw');
bnt_dta = tot_dta(:,bnt_col);
bnt_dta(10) = {55};
bnt_dta(25) = {53};
bnt_dta(27) = {39};
bnt_dta(28) = {51};
bnt_dta(29) = {47};

rvl_col = string_find(tot_col,'RAVLT_Long_Delay_Recall_raw');
rvl_dta = tot_dta(:,rvl_col);
rvl_dta(26) = {NaN};

rvl_lrn_col = string_find(tot_col,'RAVLT_Total_Trials_1_to_5_raw');
rvl_lrn_dta = tot_dta(:,rvl_lrn_col);
rvl_lrn_dta(14) = {NaN};
rvl_lrn_dta(end-7) = {NaN};
rvl_lrn_dta(end-13) = {NaN};

sde_src = tot_dta(:,sde_col); sde_src(cell2mat(cellfun(@isnumeric,sde_src,'uni',0))) = {''};
lft_tle = intersect(unique( [string_find(sde_src,'L '); string_find(sde_src,'left') ; string_find(sde_src,'l temp')]),tme_one);
rgh_tle = intersect(unique( [string_find(sde_src,'R '); string_find(sde_src,'right'); string_find(sde_src,'r temp')]),tme_one);

% BNT
fcfg = [];

fcfg.xdt = { 1                            2 };
ccc

fcfg.fce_col     = { rgb('light blue')   rgb('light red')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark blue')    rgb('dark red') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L-TLE' 'R-TLE' };
fcfg.xlm = [ 0.5 2.5 ];
fcfg.ylb = {'BNT'};
fcfg.ylm = [ 0 65 ];

fcfg.hln = 60;
fcfg.hln_col = rgb('dark green');

fcfg.out_dir = dta_dir;
fcfg.out_nme = 'BNT';

ejk_scatter(fcfg)

% RVLT
fcfg = [];

fcfg.xdt = { 1                                  2 };
fcfg.ydt = { cell2mat(rvl_dta(lft_tle,1)) cell2mat(rvl_dta(rgh_tle,1))  };

fcfg.fce_col     = { rgb('light blue')   rgb('light red')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark blue')    rgb('dark red') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L-TLE' 'R-TLE' };
fcfg.xlm = [ 0.5 2.5 ];
fcfg.ylb = {'RVLT'};
fcfg.ylm = [ 0 15 ];

fcfg.hln = 15;
fcfg.hln_col = rgb('dark green');

% fcfg.out_dir = dta_dir;
% fcfg.out_nme = 'RVLT';

ejk_scatter(fcfg)

% RVLT Learn
fcfg = [];

fcfg.xdt = { 1                                  2 };
fcfg.ydt = { cell2mat(rvl_lrn_dta(lft_tle,1)) cell2mat(rvl_lrn_dta(rgh_tle,1))  };

fcfg.fce_col     = { rgb('light blue')   rgb('light red')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark blue')    rgb('dark red') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L-TLE' 'R-TLE' };
fcfg.xlm = [ 0.5 2.5 ];
fcfg.ylb = {'RVLT'};
fcfg.ylm = [ 0 70 ];

fcfg.hln = 15;
fcfg.hln_col = rgb('dark green');

% fcfg.out_dir = dta_dir;
% fcfg.out_nme = 'RVLT';

ejk_scatter(fcfg)

%% Get Change scores
sde_col = string_find(tot_col,'side_seizure_onset');
tme_one = string_find(tot_sbj,'-T1');
tme_two = string_find(tot_sbj,'-T2');

rvl_col = string_find(tot_col,'RAVLT_Long_Delay_Recall_raw');
rvl_dta = tot_dta(:,rvl_col);
rvl_dta(26) = {NaN};
rvl_dta(35) = {NaN};

tot_sbj_nme = regexpi(tot_sbj,'-','split');
tot_sbj_nme = cat(1,tot_sbj_nme{:});
for iS = 1:numel(tme_two)
    sbj_ind(iS,1) = tme_two(iS);
    sbj_ind(iS,2) = setxor(find(strcmpi(tot_sbj_nme,tot_sbj_nme(tme_two(iS)))),tme_two(iS));
    pst_val(iS) = rvl_dta{sbj_ind(iS,1)} - rvl_dta{sbj_ind(iS,2)};
end

% RVLT
fcfg = [];

fcfg.xdt = { 1       };
fcfg.ydt = { pst_val };

fcfg.fce_col     = { rgb('light orange') };
fcfg.edg_col     = { [0 0 0]             };
fcfg.box_plt_col = { rgb('dark orange')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'TLE' };
fcfg.xlm = [ 0.5 1.5 ];
fcfg.ylb = {'RVLT Post-operative change'};
fcfg.ylm = [ -7 7 ];

fcfg.out_dir = dta_dir;
fcfg.out_nme = 'RVLT_post_operative';

fcfg.hln = 0;
fcfg.hln_col = rgb('dark green');

ejk_scatter(fcfg)
