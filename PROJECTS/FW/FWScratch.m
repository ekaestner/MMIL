% GRAINSIZE COMPARING BETWEEN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
fcfg.col     = [1 2 3];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
bet_sig_lhs_grn = mmil_ranksum_region_test(fcfg);

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  num2cell(cat(1,cellfun(@median,tme_plt_col.cmb_reg_lhs(:,2)))) ...
  num2cell(cat(1,cellfun(@median,tme_plt_col.cmb_reg_lhs(:,3)))) ...
  num2cell(cat(1,cellfun(@median,tme_plt_col.cmb_reg_lhs(:,4)))) ]

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  num2cell(cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(:,3)))) ...
  num2cell(cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(:,4)))) ]

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  num2cell(cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(:,3)))) ...
  num2cell(cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(:,3)))) ...
  num2cell(cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(:,4)))) ]

% REPETITION/FREQUENCY COMPARING BETWEEN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'rep'                'frq'            };

fcfg.col     = [1 3];

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 0;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
bet_sig_lhs_lex = mmil_ranksum_region_test(fcfg);

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  num2cell(cat(1,cellfun(@median,tme_plt_col.cmb_reg_lhs(:,2)))) ...
  num2cell(cat(1,cellfun(@median,tme_plt_col.cmb_reg_lhs(:,3)))) ]

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(:,2),'uni',0)) ...
  cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(:,3),'uni',0)) ]

[ tme_plt_col.cmb_reg_lhs(:,1) ...
  cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(:,2),'uni',0)) ...
  cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(:,3),'uni',0)) ]

%% New Table
ord_hld = [2 1 3 4 5 8 10 9];

tme_plt_col.cmb_reg_lhs(ord_hld,1)

[ tme_plt_col.cmb_reg_lhs(ord_hld,1) ...
  cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(ord_hld,3),'uni',0)) ...
  cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(ord_hld,4),'uni',0)) ...
  ...cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(ord_hld,3),'uni',0)) ...
  ...cat(1,cellfun(@max,tme_plt_col.cmb_reg_lhs(ord_hld,2),'uni',0)) ]...
  cat(1,cellfun(@min,tme_plt_col.cmb_reg_lhs(ord_hld,5),'uni',0)) ]

openfig('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure4/grainsize/scatter/combined_hgp_lhs.fig')
set(gcf,'Visible','on')

scr_sze = get(0,'Screensize');
set(gcf,'Position',[1 1 scr_sze(3) scr_sze(4)]);
colormap jet;
exportfig(gcf,'/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure4/grainsize/scatter/updatedlhs.eps','Color','rgb','SeparateText',1,'output','eps')
close all



