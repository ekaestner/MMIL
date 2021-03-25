clear; clc;

%% Alt epd006 - 313
ejk_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313.mat');
ejk_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313.mat');

rep_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd006_meg.mat');
rep_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd006_meg.mat');

row_nme_lhs = find(strcmpi( ejk_lhs.srf_dta_sbj, 'epd006'));
row_nme_rhs = find(strcmpi( ejk_rhs.srf_dta_sbj, 'epd006'));

% LHS
srf_dta     = ejk_lhs.srf_dta;
srf_dta_sbj = ejk_lhs.srf_dta_sbj;

srf_dta(row_nme_lhs,:) = rep_lhs.srf_dta(1,:);

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd006switch.mat', 'srf_dta', 'srf_dta_sbj' );

% RHS
srf_dta     = ejk_rhs.srf_dta;
srf_dta_sbj = ejk_rhs.srf_dta_sbj;

srf_dta(row_nme_rhs,:) = rep_rhs.srf_dta(1,:);

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd006switch.mat', 'srf_dta', 'srf_dta_sbj' );

%% Alt epd006 - 176
ejk_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm176.mat');
ejk_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm176.mat');

rep_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm176_epd006_meg.mat');
rep_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm176_epd006_meg.mat');

row_nme_lhs = find(strcmpi( ejk_lhs.srf_dta_sbj, 'epd006'));
row_nme_rhs = find(strcmpi( ejk_rhs.srf_dta_sbj, 'epd006'));

% LHS
srf_dta     = ejk_lhs.srf_dta;
srf_dta_sbj = ejk_lhs.srf_dta_sbj;

srf_dta(row_nme_lhs,:) = rep_lhs.srf_dta(1,:);

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm176_epd006switch.mat', 'srf_dta', 'srf_dta_sbj' );

% RHS
srf_dta     = ejk_rhs.srf_dta;
srf_dta_sbj = ejk_rhs.srf_dta_sbj;

srf_dta(row_nme_rhs,:) = rep_rhs.srf_dta(1,:);

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm176_epd006switch.mat', 'srf_dta', 'srf_dta_sbj' );

%% Add epd096 - 313
ejk_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd006switch.mat');
ejk_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd006switch.mat');

add_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd096.mat');
add_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd096.mat');

% LHS
srf_dta     = [ ejk_lhs.srf_dta ;   add_lhs.srf_dta ];
srf_dta_sbj = [ ejk_lhs.srf_dta_sbj 'epd096' ];

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd006switch_epd096add.mat', 'srf_dta', 'srf_dta_sbj' );

% RHS
srf_dta     = [ ejk_rhs.srf_dta ;   add_rhs.srf_dta ];
srf_dta_sbj = [ ejk_rhs.srf_dta_sbj 'epd096' ];

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd006switch_epd096add.mat', 'srf_dta', 'srf_dta_sbj' );

% Check Sbj
[sbj_inc([2:27 29:end]) srf_dta_sbj(1:end)']

%% Add epd096 to original files - 313
ejk_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313.mat');
ejk_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313.mat');

add_lhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd096.mat');
add_rhs = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd096.mat');

% LHS
srf_dta     = [ ejk_lhs.srf_dta ;   add_lhs.srf_dta ];
srf_dta_sbj = [ ejk_lhs.srf_dta_sbj 'epd096' ];

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd096add.mat', 'srf_dta', 'srf_dta_sbj' );

% RHS
srf_dta     = [ ejk_rhs.srf_dta ;   add_rhs.srf_dta ];
srf_dta_sbj = [ ejk_rhs.srf_dta_sbj 'epd096' ];

save( '/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd096add.mat', 'srf_dta', 'srf_dta_sbj' );

% Check Sbj
[sbj_inc([2:27 29:end]) srf_dta_sbj(1:end)']

