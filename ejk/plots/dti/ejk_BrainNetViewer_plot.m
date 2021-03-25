%%


%%
function connectome_plot(cfg)

% Hardcode
if isfield(cfg,'edg_cut_off'); cfg.edg_cut_off = 0; end

bse_dir = '/home/ekaestne/PROJECTS/SCRIPTS/BCT/BrainNetViewer_20181219/Data/ModDesikan';

% Color Map
col_map = [ cat(1,cfg.edg_col{:}) ; [1 1 1] ];

% Edge Names
edg_clr = [];
for iEN = 1:numel(cfg.edg_nme)
    edg_nme{iEN} = cfg.edg_nme{iEN}(:,1);
    edg_nme{iEN} = regexp(edg_nme{iEN},'==','split');
    edg_nme{iEN} = cat(1,edg_nme{iEN}{:});
    edg_nme{iEN} = edg_nme{iEN}(cfg.edg_wgh{iEN}>=cfg.edg_cut_off,:);
    edg_clr = [ edg_clr ; ones(size(cfg.edg_nme{iEN},1),1) * iEN ]; 
end
edg_nme = cat(1,edg_nme{:});
edg_wgh = cat(1,cfg.edg_wgh{:});

% Brain
brn_fle = [bse_dir '/' cfg.brn_nme '.nv' ];
nde_fle = [bse_dir '/' cfg.nde_nme '.node' ];

% Nodes
nde_lst = mmil_readtext([ bse_dir '/' cfg.nde_nme '.node'],[' ']);
nde_nme = unique([edg_nme(:,1) ; edg_nme(:,2)]);
nde_fcs = nde_lst(ismember(nde_lst(:,end),cellfun(@(x) x([1 3:end]),nde_nme,'uni',0)),:);
cell2csv([ bse_dir '/' cfg.nde_nme '_' 'tmp' '.node'],[repmat({'#'},1,size(nde_fcs,2)) ; nde_fcs ],' ')
nde_fle_tmp = [ bse_dir '/' cfg.nde_nme '_' 'tmp' '.node'];

% Edge Weight & Color Matrix
plt_mtx     = zeros(size(nde_fcs,1),size(nde_fcs,1));
clr_cus_mtx = zeros(size(nde_fcs,1),size(nde_fcs,1));
edg_sze     = zeros(size(edg_nme,1),1);
for iCN = 1:size(edg_nme,1)
    
    plt_mtx(strcmpi(nde_fcs(:,end),edg_nme{iCN,1}([1 3:end])),strcmpi(nde_fcs(:,end),edg_nme{iCN,2}([1 3:end]))) = 1;
    plt_mtx(strcmpi(nde_fcs(:,end),edg_nme{iCN,2}([1 3:end])),strcmpi(nde_fcs(:,end),edg_nme{iCN,1}([1 3:end]))) = 1;
    
    clr_cus_mtx(strcmpi(nde_fcs(:,end),edg_nme{iCN,1}([1 3:end])),strcmpi(nde_fcs(:,end),edg_nme{iCN,2}([1 3:end]))) = edg_clr(iCN,1);
    clr_cus_mtx(strcmpi(nde_fcs(:,end),edg_nme{iCN,2}([1 3:end])),strcmpi(nde_fcs(:,end),edg_nme{iCN,1}([1 3:end]))) = edg_clr(iCN,1);
    
    edg_sze(iCN,1) = edg_wgh(iCN,1);
    
end

% Put together Plotting CFG 
cnn_lbl = cell(1,size(plt_mtx,2));
cnn_lbl{1} = '#';
cell2csv([bse_dir  '/' 'edge' '.edge'],[cnn_lbl ; num2cell(plt_mtx)],'\t')
edg_fle = [bse_dir '/' 'edge' '.edge'];

load([bse_dir '/' 'EC.mat' ])
EC.msh.alpha = 0.4;
EC.msh.color = [0.8,0.8,0.8];
EC.edg.interhemiedges = 0;
EC.edg.size           = 1;
EC.edg.size_size      = edg_sze(edg_sze~=0);
EC.edg.CM             = col_map; % [rgb('bright blue') ; repmat(rgb('bright red'),62,1) ; rgb('bright green')];
EC.edg.color          = 6;
EC.edg.color_custom_matrix = clr_cus_mtx;
EC.edg.color_custom_index  = 1:size(col_map,1);
EC.nod.CM = repmat(rgb('black') / 2 ,64,1);
if all(EC.nod.CM(1,:)==0); EC.nod.CM = repmat(rgb('black'),64,1); end
EC.lot.view  = 2;
EC.img.width  = 2000;
EC.img.height = 1500;
EC.img.dpi = 100;close all

save([bse_dir '/' 'EC_tmp.mat' ],'EC')
if ~exist(cfg.out_dir); mkdir(cfg.out_dir); end
mmil_BrainNet_MapCfg(brn_fle,nde_fle_tmp,edg_fle,[bse_dir '/' 'EC_tmp.mat' ],[cfg.out_dir '/' cfg.out_nme '_' 'brain' '.png'])
close all

end