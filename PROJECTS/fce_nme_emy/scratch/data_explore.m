clear; clc

dta_loc = '/home/ekaestner/';

dta_loc_sub = 'Dropbox/McDonald Lab/Erik/Projects/Imaging/fce_nme_emy/';

out_plt_dir = [ dta_loc '/' dta_loc_sub '/' 'Behavior' '/' ];    

out_sht_col = { 'sbj_nme'     'tme_pnt_one'     'tme_pnt_two' 'scn_one' 'scn_two' ...
                'run_one_rcl' 'run_one_frc_chc' 'run_two_rcl' 'run_two_frc_chc' 'chg_rcl' 'chg_chc' };
            
%% Load Data
sbj_nme = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'ScanData' '/' 'sub_list.txt']);
    sbj_nme(cellfun(@isnumeric,sbj_nme(:,1)),1) = cellfun(@num2str,sbj_nme(cellfun(@isnumeric,sbj_nme(:,1)),1),'uni',0);
    cell2csv([ out_plt_dir '/' 'subject_names.csv'],sbj_nme)

dta_cph = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'ScanData' '/' 'cipher.tsv'],['\t']);
    dta_cph(cellfun(@isnumeric,dta_cph(:,1)),1) = cellfun(@num2str,dta_cph(cellfun(@isnumeric,dta_cph(:,1)),1),'uni',0);
    
red_cap_dta = mmil_readtext([ dta_loc '/' dta_loc_sub '/' 'RedcapData' '/' 'RedcapData_pipe.csv'],'|');
    red_cap_dta(cellfun(@isnumeric,red_cap_dta(:,1)),1) = cellfun(@num2str,red_cap_dta(cellfun(@isnumeric,red_cap_dta(:,1)),1),'uni',0);

out_csv = cell(numel(sbj_nme),numel(out_sht_col));

%% Populate sheet

% Redcap
sbj_nme_col = strcmpi(out_sht_col,'sbj_nme');
tme_pnt_one_col = strcmpi(out_sht_col,'tme_pnt_one');
tme_pnt_two_col = strcmpi(out_sht_col,'tme_pnt_two');

for iS = 1:numel(sbj_nme)
    
    out_csv{iS,sbj_nme_col} = sbj_nme{iS};
    
    ind_one = find(strcmpi(red_cap_dta(:,1),sbj_nme{iS}));
    ind_two = find(strcmpi(red_cap_dta(:,1),[sbj_nme{iS} '_2']));
    
    if ~isempty(ind_one); out_csv{iS,tme_pnt_one_col} = 1; else out_csv{iS,tme_pnt_one_col} = 0; end
    if ~isempty(ind_two); out_csv{iS,tme_pnt_two_col} = 1; else out_csv{iS,tme_pnt_two_col} = 0; end
        
end

% Number of complete records
scn_one_col = strcmpi(out_sht_col,'scn_one');
scn_two_col = strcmpi(out_sht_col,'scn_two');

for iS = 1:numel(sbj_nme)
    sbj_ind = strcmpi(dta_cph(:,1),sbj_nme{iS});
    
    if any(cell2mat(dta_cph(sbj_ind,2)==1)); out_csv{iS,scn_one_col} = 1; else out_csv{iS,scn_one_col} = 0; end
    if any(cell2mat(dta_cph(sbj_ind,2)==2)); out_csv{iS,scn_two_col} = 1; else out_csv{iS,scn_two_col} = 0; end
end
   
% Subject behavioral numbers
sum_tme_pnt_one = sum(cell2mat(out_csv(:,tme_pnt_one_col)));
sum_tme_pnt_two = sum(cell2mat(out_csv(:,tme_pnt_two_col)));

sum_scn_one_and_two = sum(logical(cell2mat(out_csv(:,tme_pnt_one_col))) & logical(cell2mat(out_csv(:,tme_pnt_two_col))));
sum_scn_one_oly     = sum(logical(cell2mat(out_csv(:,tme_pnt_one_col))) & ~logical(cell2mat(out_csv(:,tme_pnt_two_col))));
sum_tme_pnt_two_oly     = sum(~logical(cell2mat(out_csv(:,tme_pnt_one_col))) & logical(cell2mat(out_csv(:,tme_pnt_two_col))));

% Scan Subject numbers
sum_scn_one = sum(cell2mat(out_csv(:,scn_one_col)));
sum_scn_two = sum(cell2mat(out_csv(:,scn_two_col)));

sum_scn_one_and_two = sum(logical(cell2mat(out_csv(:,scn_one_col))) & logical(cell2mat(out_csv(:,scn_two_col))));
sum_scn_one_oly     = sum(logical(cell2mat(out_csv(:,scn_one_col))) & ~logical(cell2mat(out_csv(:,scn_two_col))));
sum_scn_two_oly     = sum(~logical(cell2mat(out_csv(:,scn_one_col))) & logical(cell2mat(out_csv(:,scn_two_col))));

% Scan & Behavioral Subject numbers
scn_one_and_beh_one = logical(cell2mat(out_csv(:,scn_one_col))) & logical(cell2mat(out_csv(:,tme_pnt_one_col)));
scn_two_and_beh_two = logical(cell2mat(out_csv(:,scn_two_col))) & logical(cell2mat(out_csv(:,tme_pnt_two_col)));

sum_bth_one = sum(scn_one_and_beh_one);
sum_bth_two = sum(scn_two_and_beh_two);

sum_bth_one_and_two = sum(scn_one_and_beh_one & scn_two_and_beh_two);
sum_bth_one_oly     = sum(scn_one_and_beh_one & ~scn_two_and_beh_two);
sum_bth_two_oly     = sum(~scn_one_and_beh_one & scn_two_and_beh_two);

% Get missing patients
mss_scn_one = out_csv(~logical(cell2mat(out_csv(:,scn_one_col))) & logical(cell2mat(out_csv(:,tme_pnt_one_col))),1);
mss_scn_two = out_csv(~logical(cell2mat(out_csv(:,scn_two_col))) & logical(cell2mat(out_csv(:,tme_pnt_two_col))),1);

mss_beh_one = out_csv(logical(cell2mat(out_csv(:,scn_one_col))) & ~logical(cell2mat(out_csv(:,tme_pnt_one_col))),1);
mss_beh_two = out_csv(logical(cell2mat(out_csv(:,scn_two_col))) & ~logical(cell2mat(out_csv(:,tme_pnt_two_col))),1);

%% Populate behavior
cor_ans = { 'Shawn'  3 ; ...'Shawn'  3 ;
            'Sandra' 2 ; ...'Sandra' 2 ;
            'Eugene' 3 ; ...'Eugene' 3 ;
            'Andrea' 1 ; ...'Connie' 2 ;
            'Daniel' 1 ; ...'Daniel' 1 ;
            'Erica'  2 ; ...'Erica'  2 ;
            'Laura'  3 ; ...'Laura'  3 ;
            'Gerald' 1 ; ...'Joseph' 3 ;
            'Roger'  3 ; ...'Roger'  3 ;
            'Larry'  2 ; ...'Shawn'  1 ;
            'Joseph' 1 ; ...'Daniel' 2 ;
            'Willie' 3 ; ...'Willie' 3 ;
            'Irene'  2 ; ...'Irene'  2 ;
            'Tracy'  3 ; ...'Tracy'  3 ;
            'Connie' 3 ; ...'Joyce'  2 ;
            'Elaine' 2 ; ...'Elaine' 2 ;
            'Monica' 1 ; ...'Joyce'  2 ;
            'Nancy'  2 ; ...'Erica'  1 ;
            'Brenda' 3 ; ...'Carol'  1 ;
            'Dennis' 1 ; ...'Larry'  3 ;
            'Jacob'  2 ; ...'Daniel' 1 ;
            'Harold' 1 ; ...'Harold' 1 ;
            'David'  2 ; ...'David'  2 ;
            'Randy'  1 ; ...'Shawn'  3 ;
            'Scott'  2 ; ...'David'  3 ;
            'Lauren' 1 ; ...'Lauren' 1 ;
            'Wanda'  3 ; ...'Wanda'  3 ;
            'Sheila' 2 ; ...'Elaine' 1 ;
            'Joyce'  3 ; ...'Irene'  2 ;
            'Carol'  1 }; %'Carol'  1 };

run_one_rcl_col     = strcmpi(out_sht_col,'run_one_rcl');
run_one_frc_chc_col = strcmpi(out_sht_col,'run_one_frc_chc');

run_two_rcl_col     = strcmpi(out_sht_col,'run_two_rcl');
run_two_frc_chc_col = strcmpi(out_sht_col,'run_two_frc_chc');

chg_rcl_col     = strcmpi(out_sht_col,'chg_rcl');
chg_chc_col = strcmpi(out_sht_col,'chg_chc');

for iS = 1:numel(sbj_nme)
        
    ind_one = find(strcmpi(red_cap_dta(:,1),sbj_nme{iS}));
    ind_two = find(strcmpi(red_cap_dta(:,1),[sbj_nme{iS} '_2']));
        
    if ~isempty(ind_one)
        
        rcl_ind = string_find(red_cap_dta(1,:),{'free_recall'});
        chc_ind = string_find(red_cap_dta(1,:),{'mult_choice'});
        
        rcl_cor = 0;
        chc_cor = 0;
        for iN = 1:numel(rcl_ind)
            if strcmpi(red_cap_dta{ind_one,rcl_ind(iN)},cor_ans{iN,1}); rcl_cor=rcl_cor+1; end
            if red_cap_dta{ind_one,chc_ind(iN)}==cor_ans{iN,2};         chc_cor=chc_cor+1; end
        end
        
        out_csv{iS,run_one_rcl_col}     = round((rcl_cor/numel(rcl_ind)) * 100);
        out_csv{iS,run_one_frc_chc_col} = round((chc_cor/numel(chc_ind)) * 100);  
        
    else
        out_csv{iS,run_one_rcl_col}     = NaN;
        out_csv{iS,run_one_frc_chc_col} = NaN;        
    end
    
    
    if ~isempty(ind_two)
        
        rcl_ind = string_find(red_cap_dta(1,:),{'free_recall'});
        chc_ind = string_find(red_cap_dta(1,:),{'mult_choice'});
        
        rcl_cor = 0;
        chc_cor = 0;
        for iN = 1:numel(rcl_ind)
            if strcmpi(red_cap_dta{ind_two,rcl_ind(iN)},cor_ans{iN,1}); rcl_cor=rcl_cor+1; end
            if red_cap_dta{ind_two,chc_ind(iN)}==cor_ans{iN,2};         chc_cor=chc_cor+1; end
        end
        
        out_csv{iS,run_two_rcl_col}     = round((rcl_cor/numel(rcl_ind)) * 100);
        out_csv{iS,run_two_frc_chc_col} = round((chc_cor/numel(chc_ind)) * 100);  
        
    else
        out_csv{iS,run_two_rcl_col}     = NaN;
        out_csv{iS,run_two_frc_chc_col} = NaN;        
    end

    out_csv{iS,chg_rcl_col} = out_csv{iS,run_two_rcl_col} - out_csv{iS,run_one_rcl_col};
    out_csv{iS,chg_chc_col} = out_csv{iS,run_two_frc_chc_col} - out_csv{iS,run_one_frc_chc_col};
            
end

% Average performance
[   round(nanmean(cell2mat(out_csv(:,run_one_rcl_col)))),round(nanstd(cell2mat(out_csv(:,run_one_rcl_col)))) ...
    round(nanmean(cell2mat(out_csv(:,run_one_frc_chc_col)))),round(nanstd(cell2mat(out_csv(:,run_one_frc_chc_col)))); ...
    round(nanmean(cell2mat(out_csv(:,run_two_rcl_col)))),round(nanstd(cell2mat(out_csv(:,run_two_rcl_col)))) ...
    round(nanmean(cell2mat(out_csv(:,run_two_frc_chc_col)))),round(nanstd(cell2mat(out_csv(:,run_two_frc_chc_col)))); ...
    round(nanmean(cell2mat(out_csv(:,chg_rcl_col)))),round(nanstd(cell2mat(out_csv(:,chg_rcl_col)))) ...
    round(nanmean(cell2mat(out_csv(:,chg_chc_col)))),round(nanstd(cell2mat(out_csv(:,chg_chc_col))))]


% First Run Plot
fcfg = [];

fcfg.xdt     = { 1 2 3 4 5 };
fcfg.ydt     = { cell2mat(out_csv(:,run_one_rcl_col)) cell2mat(out_csv(:,run_one_frc_chc_col)) [] cell2mat(out_csv(:,run_two_rcl_col)) cell2mat(out_csv(:,run_two_frc_chc_col)) };

fcfg.fce_col = { rgb('light mauve') rgb('light tan') rgb('grey') rgb('dark mauve') rgb('dark tan')  };
fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'Recall Change' 'Forced Choice Change' };
fcfg.ylb = { '% Correct Run2 - Run1'  };

fcfg.hln     = 0;
fcfg.hln_col = rgb('black');

fcfg.xlm = [ 0 6 ];
fcfg.ylm = [ -10 100 ];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Run1_Run2';

ejk_scatter(fcfg)

% Change plot
fcfg = [];

fcfg.xdt     = { 1 2 };
fcfg.ydt     = { cell2mat(out_csv(:,chg_rcl_col)) cell2mat(out_csv(:,chg_chc_col)) };

fcfg.fce_col = { rgb('mauve') rgb('tan') };
fcfg.edg_col = repmat({rgb('black')},1,numel(fcfg.xdt));

fcfg.xlb = { 'Run1 Recall' 'Run1 Forced Choice' '' 'Run2 Recall' 'Run2 Forced Choice' };
fcfg.ylb = { '% Correct'  };

fcfg.hln     = 0;
fcfg.hln_col = rgb('black');

fcfg.xlm = [ 0 3 ];
fcfg.ylm = [ -40 40 ];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Change';

ejk_scatter(fcfg)

%% Missing subjects

