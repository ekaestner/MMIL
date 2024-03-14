clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Adam/epd_t2h/';
dta_dir = [ prj_dir '/' 'Data' '/'];
plt_dir = [ prj_dir '/' 'Output' '/' 'v2'];
plt_dir_new = [ prj_dir '/' 'Output' '/' 'v2_subgroup'];

dta_fle = 'ashs_combined_v2.csv';

cpe_fle = 'T2AgeReport_forErik.csv';

%%
grp_loc = 'Group_TLE';
grp_nme = {  'Y-HC'          'Y-LTLE'           'Y-RTLE'         'O-Epilepsy'       }; % {  'Y-HC' 'Y-TLE' 'O-HC' 'O-Epilepsy' 'MCI-C'};
grp_pos = [  1                2                 3                5                  ]; % [  1      2       4      5            6];
grp_col = { rgb('light grey') rgb('light blue') rgb('light red') rgb('dark orange') }; % { rgb('light grey') rgb('light orange') rgb('dark grey') rgb('dark orange') rgb('dark purple') };

sub_grp_loc = 'Group_Old';
sub_grp_nme = {  'O-HC'          'Refractory-LTLE-Early'  'Refractory-RTLE-Early' 'Benign-Early'     'Benign-Late'    };
sub_grp_pos = [  1                3                 4                 5                  6                  ]; 
sub_grp_col = { rgb('light grey') rgb('light blue') rgb('light red')  rgb('dark orange') rgb('dark yellow') };

%%
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_fle ];
[ hip_sub_dta, hip_sub_sbj, hip_sub_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' cpe_fle ];
[ cpe_sub_dta, cpe_sub_sbj, cpe_sub_col ] = ejk_dta_frm(fcfg);

% unique(hip_sub_dta(:,strcmpi(hip_sub_col,'Group_TLE')));

%%
grp_var = hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc));

grp_var( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),'O-Epilepsy') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age of Disease Onset')))>=55, 1) = {'Benign-Late'};
grp_var( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),'O-Epilepsy') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age of Disease Onset')))<55, 1) = {'Benign-Early'};
grp_var( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),'Y-LTLE') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age')))>=55, 1) = {'Refractory-LTLE-Early'};
grp_var( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),'Y-RTLE') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age')))>=55, 1) = {'Refractory-RTLE-Early'};
grp_var( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),'Y-HC') & cell2mat(hip_sub_dta(:,strcmpi(hip_sub_col,'Age')))>=55, 1) = {'O-HC'};


hip_sub_col = [ hip_sub_col 'Group_Old' ];
hip_sub_dta = [ hip_sub_dta grp_var ];

plt_fld = hip_sub_col;

%%
for iF = 1:numel(plt_fld)
    
    if isnumeric(hip_sub_dta{1,iF}) && ~isnan(hip_sub_dta{1,iF})
        
        fcfg = [];
        
        for iG = 1:numel(grp_nme)
            fcfg.xdt{iG}         = grp_pos(iG);
            fcfg.ydt{iG}         = cell2mat(hip_sub_dta( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),grp_nme{iG}), strcmpi(hip_sub_col,plt_fld{iF}) ));
            fcfg.fce_col{iG}     = grp_col{iG};
            fcfg.box_plt_col{iG} = grp_col{iG};
        end
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        
        fcfg.xlm = [ 0.5 max(cell2mat(fcfg.xdt))+0.5 ];
        fcfg.xlb = grp_nme;
        
        fcfg.ylb = plt_fld(iF);
        
        fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
        fcfg.aph_val = 0.45;
        
        fcfg.out_dir = plt_dir;
        fcfg.out_nme = plt_fld{iF};
        
        ejk_scatter(fcfg)
        
    end
end

%%
for iF = 1:numel(plt_fld)
    
    if isnumeric(hip_sub_dta{1,iF}) && ~isnan(hip_sub_dta{1,iF})
        
        fcfg = [];
        
        for iG = 1:numel(sub_grp_nme)
            fcfg.xdt{iG}         = sub_grp_pos(iG);
            fcfg.ydt{iG}         = cell2mat(hip_sub_dta( strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,sub_grp_loc)),sub_grp_nme{iG}), strcmpi(hip_sub_col,plt_fld{iF}) ));
            fcfg.fce_col{iG}     = sub_grp_col{iG};
            fcfg.box_plt_col{iG} = sub_grp_col{iG};
        end
        
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        
        fcfg.xlm = [ 0.5 max(cell2mat(fcfg.xdt))+0.5 ];
        fcfg.xlb = sub_grp_nme;
        
        fcfg.ylb = plt_fld(iF);
        
        fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
        fcfg.aph_val = 0.45;
        
        fcfg.out_dir = plt_dir_new;
        fcfg.out_nme = plt_fld{iF};
        
        ejk_scatter(fcfg)
        
    end
end

%% Age & Age of Onset
cpe_has_aon = ~cellfun(@isempty,cpe_sub_dta(:,strcmpi(cpe_sub_col,'Age at Onset')));

% CAPES age of onset & Age
fcfg = [];

fcfg.xdt{1}          = cell2mat(cpe_sub_dta(cpe_has_aon,strcmpi(cpe_sub_col,'Age at Onset')));
fcfg.ydt{1}          = cell2mat(cpe_sub_dta(cpe_has_aon,strcmpi(cpe_sub_col,'Age')));
fcfg.fce_col{1}     = rgb('light orange');

fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];

fcfg.xlm = [ 0 75 ];
fcfg.xlb = {'Age at Onset'};
fcfg.ylb = {'Age'};

fcfg.hln = 55;
fcfg.hln_col = rgb('dark grey');

fcfg.vln = 55;
fcfg.vln_col = rgb('dark grey');

fcfg.out_dir = plt_dir;
fcfg.out_nme = 'Age_CAPES';

ejk_scatter(fcfg)

% CAPES missing age of onset
fcfg = [];

fcfg.xdt{1}          = 1;
fcfg.ydt{1}          = cell2mat(cpe_sub_dta(~cpe_has_aon,strcmpi(cpe_sub_col,'Age')));
fcfg.fce_col{1}     = rgb('light orange');

fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];

fcfg.xlm = [ 0.5 1.5 ];

fcfg.xlb = {'CAPES EPD'};
fcfg.ylb = {'Age'};

fcfg.hln = 55;
fcfg.hln_col = rgb('dark grey');

fcfg.out_dir = plt_dir;
fcfg.out_nme = 'Age_CAPES_missing_onset';

ejk_scatter(fcfg)

% UCSD Data
fcfg = [];

for iG = 1:numel(grp_nme)
    grp_row              = strcmpi(hip_sub_dta(:,strcmpi(hip_sub_col,grp_loc)),grp_nme{iG});
    fcfg.xdt{iG}         = cell2mat(hip_sub_dta(grp_row,strcmpi(hip_sub_col,'Age of Disease Onset')));
    fcfg.ydt{iG}         = cell2mat(hip_sub_dta(grp_row,strcmpi(hip_sub_col,'Age')));
    fcfg.fce_col{iG}     = grp_col{iG};
end

fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];

fcfg.xlm = [ 0 75 ];
fcfg.xlb = {'Age at Onset'};
fcfg.ylb = {'Age'};

fcfg.hln = 55;
fcfg.hln_col = rgb('dark grey');

fcfg.vln = 55;
fcfg.vln_col = rgb('dark grey');

fcfg.out_dir = plt_dir;
fcfg.out_nme = 'Age_UCSD';

ejk_scatter(fcfg)












