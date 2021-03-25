%% nyu_plot
function nyu_plot2(pcfg)

if ~exist('pcfg.label','var'); pcfg.label = 0; end
if ~exist('pcfg.alpha','var') || isempty(pcfg.alpha); pcfg.alpha = 1; end
if isfield(pcfg,'rsp_ele') && isempty(pcfg.radius); pcfg.radius = repmat(2,1,numel(pcfg.rsp_ele)); end
if isfield(pcfg,'rsp_ele') && (numel(pcfg.radius) ~= numel(pcfg.rsp_ele)); pcfg.radius = repmat(pcfg.radius,1,numel(pcfg.rsp_ele)); end
if ~isfield(pcfg,'aparc'); pcfg.aparc = []; end

if ~isfield(pcfg,'par_ele_lne') && isfield(pcfg,'par_ele'); pcfg.par_ele_lne = repmat(1.25,1,numel(pcfg.par_ele)); end
if ~isfield(pcfg,'par_ele_col') && isfield(pcfg,'par_ele'); pcfg.par_ele_col = repmat({rgb('black')},1,numel(pcfg.par_ele)); end

if isfield(pcfg.surf_brain,'coords')==0
    sub.vert = pcfg.surf_brain.surf_brain.coords;
    sub.tri = pcfg.surf_brain.surf_brain.faces;
else
    sub.vert = pcfg.surf_brain.coords;
    sub.tri = pcfg.surf_brain.faces;
end

if ~isfield(pcfg,'aparc') || isempty(pcfg.aparc) && ~isfield(pcfg,'tbl_pct')
    col = [.7 .7 .7];
    col = repmat(col(:)', [size(sub.vert, 1) 1]);
elseif ~isfield(pcfg,'tbl_pct')
    col = repmat(rgb('white'), [size(sub.vert, 1) 1]);
    [col_loc,albl,actbl]=fs_read_annotation(pcfg.aparc);
    if size(albl,1) == size(col_loc,1)
        for iCL = 1:size(col_loc,1); if albl(iCL);    col(iCL,:) = (actbl.table(albl(iCL)==actbl.table(:,5),1:3))/255; end; end
    else
        for iCL = 1:size(col_loc,1); if col_loc(iCL); col(iCL,:) = (actbl.table(col_loc(iCL),1:3))/255; end; end
    end
elseif isfield(pcfg,'tbl_pct') && isnumeric(pcfg.tbl_loc)   
    
    col = repmat(rgb('white'), [size(sub.vert, 1) 1]);
    for iR = 1:numel(pcfg.tbl_loc)
        if pcfg.tbl_pct(iR)==0; pcfg.tbl_pct(iR)= .0001 ; end
        col(pcfg.tbl_loc(iR),:) = pcfg.col_map(ceil(pcfg.tbl_pct(iR)*1000),:);        
    end
    
else
    [col_loc,albl,actbl]=fs_read_annotation(pcfg.aparc);
        
    for iLC = 1:numel(albl)
        tbl_ind = string_find(cellfun(@(x) x(4:end),pcfg.tbl_loc,'uni',0),albl{iLC});
        if ~isempty(tbl_ind) && sum(col_loc==iLC)>0
            if round(pcfg.tbl_pct(tbl_ind)*(round(pcfg.top_pct*1000)))==0; pcfg.tbl_pct(tbl_ind)=0.001; end
            try
                col_hld(iLC,:) = pcfg.col_map(round(pcfg.tbl_pct(tbl_ind)*1000),:);
            catch
                col_hld(iLC,:) = pcfg.col_map(size(pcfg.col_map,1),:);
            end
        else
            col_hld(iLC,:) = rgb('white');%col_hld(iLC,:) = rgb('dark grey');
        end
    end
    
    col = repmat(rgb('white'), [size(sub.vert, 1) 1]);
    col_loc_typ = unique(col_loc);
    
    if size(albl,1) == size(col_loc,1)
        error('Fix Me')
    else
        for iCL = 1:size(col_hld,1); 
            if any(iCL==col_loc)
                col(iCL==col_loc,:) = repmat(col_hld(iCL,:),sum(iCL==col_loc),1); 
            end
        end
    end
    
end

%% Plot 3d object
trisurf(sub.tri, sub.vert(:, 1), sub.vert(:, 2), sub.vert(:, 3),...
    'FaceVertexCData', col,'FaceColor', 'interp','FaceAlpha',pcfg.alpha,'parent',pcfg.axe_hnd);

shading(pcfg.axe_hnd,'interp');
lighting(pcfg.axe_hnd,'gouraud');
material(pcfg.axe_hnd,'dull');
axis(pcfg.axe_hnd,'off')
% axes(pcfg.axe_hnd,'Visible','off'); %light; 
if ~isempty(pcfg.sve_img) || pcfg.plv; set(pcfg.fig_hdl,'Visible','off'); end

%% Tighten Plot
x_limit = [min(sub.vert(:, 1)) max(sub.vert(:, 1))];
x_limit(3) = (x_limit(2)-x_limit(1))*0.05;
xlim([x_limit(1)-x_limit(3) x_limit(2)+x_limit(3)])

y_limit = [min(sub.vert(:, 2)) max(sub.vert(:, 2))];
y_limit(3) = (y_limit(2)-y_limit(1))*0.05;
ylim([y_limit(1)-y_limit(3) y_limit(2)+y_limit(3)])

z_limit = [min(sub.vert(:, 3)) max(sub.vert(:, 3))];
z_limit(3) = (z_limit(2)-z_limit(1))*0.05;
zlim([z_limit(1)-z_limit(3) z_limit(2)+z_limit(3)])

set(pcfg.axe_hnd ,'LooseInset',get(pcfg.axe_hnd ,'TightInset')) % Shrink everything to fit);

%% Add Responsive electrodes
if isfield(pcfg,'rsp_ele') && isfield(pcfg,'color')
    hold on;
    
    if iscell(pcfg.color) && ~isempty(pcfg.nsl_col)
        colormap(cell2mat(cat(1,pcfg.color{:},pcfg.nsl_col))); freezeColors
    elseif iscell(pcfg.color)
        colormap(cell2mat(cat(1,pcfg.color{:},pcfg.nsl_col(1)))); freezeColors
    else
        colormap(cell2mat(cat(1,pcfg.color,pcfg.nsl_col(1)))); freezeColors
    end
    
    for iGP = 1:numel(pcfg.rsp_ele)
        dbl_skp{iGP} = zeros(1,numel(pcfg.ele_lbl{iGP}));
        oth_num = logical(ones(1,numel(pcfg.rsp_ele))); oth_num(iGP) = 0;
        for i=1:size(pcfg.rsp_ele{iGP},1)
            dbl_skp{iGP}(i) = any(cellfun(@any,cellfun(@(x) ismember(x,pcfg.ele_lbl{iGP}(i)),pcfg.ele_lbl(oth_num),'uni',0)));
            if iGP>1 && any(cellfun(@any,cellfun(@(x) ismember(x,pcfg.ele_lbl{iGP}(i)),pcfg.ele_lbl(1:iGP-1),'uni',0)))
                dbl_skp{iGP}(i) = dbl_skp{iGP}(i)+1;
            end
        end
    end
    
    for iGP = 1:numel(pcfg.rsp_ele)
%         dbl_skp = zeros(1,numel(pcfg.rsp_ele{iGP}));
        for i=1:size(pcfg.rsp_ele{iGP},1)
            
            %if dbl_skp(i) < 1 && sum(~cellfun(@isempty,cellfun(@(x) intersect(pcfg.ele_lbl{iGP}{i},x),pcfg.ele_lbl,'uni',0)))>1
            if dbl_skp{iGP}(i)==1
                tmp_col = ~cellfun(@isempty,cellfun(@(x) intersect(pcfg.rsp_ele{iGP}(i,1),x),pcfg.rsp_ele,'uni',0));
                col = repmat(repmat(sort(repmat([find(tmp_col)],1,2)),1,ceil(15/numel(find(tmp_col)))),30,1);
                col = col(1:30,1:30)/size(colormap,1)-(1/(size(colormap,1)*2));
            %elseif dbl_skp(i) < 1
            elseif dbl_skp{iGP}(i)==0
                col = repmat(iGP,30,30)/size(colormap,1)-(1/(size(colormap,1)*2));
            end
            
            if dbl_skp{iGP}(i) < 2
                if ~isfield(pcfg,'par_ele')
                    plotSpheres(pcfg.rsp_ele{iGP}(i,1),pcfg.rsp_ele{iGP}(i,2),pcfg.rsp_ele{iGP}(i,3),pcfg.radius(iGP),col,pcfg.sph_vew);
                elseif isfield(pcfg,'par_ele')
                    if mean(sub.vert(:, 1))<0
                       plotSpheres(min(sub.vert(:, 1))+(min(sub.vert(:, 1))*.025),pcfg.rsp_ele{iGP}(i,2),pcfg.rsp_ele{iGP}(i,3),pcfg.radius(iGP),col,pcfg.sph_vew);
                    elseif mean(sub.vert(:, 1))>0
                        plotSpheres(max(sub.vert(:, 1))+(max(sub.vert(:, 1))*.025),pcfg.rsp_ele{iGP}(i,2),pcfg.rsp_ele{iGP}(i,3),pcfg.radius(iGP),col,pcfg.sph_vew);
                    end
                end
                if pcfg.label==1
                    [x, y, z] = adjust_elec_label(pcfg.rsp_ele{iGP}(i,:),pcfg.radius(iGP));
                    text('Position',[x y z],'String',ele_lbl{iGP}(i,:),'color','w','VerticalAlignment','top','FontSize',3);
                    text('Position',[x y z],'String',ele_lbl{iGP}(i,:),'color','k','VerticalAlignment','top','FontSize',4,'FontWeight','bold');
                elseif isfield(pcfg,'ele_txt')
                    [x, y, z] = adjust_elec_label(pcfg.rsp_ele{iGP}(i,:),pcfg.radius(iGP));
                    text('Position',[x y z]+[0 0 0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]-[0 0 0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]+[0 0.2 0],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]-[0 0.2 0],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]+[0 0.2 0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]-[0 0.2 0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]+[0 -0.2 0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z]-[0 0.2 -0.2],'String',pcfg.ele_txt{iGP}(i,:),'color',[0 0 0],'VerticalAlignment','baseline','FontSize',12,'FontWeight','bold');
                    text('Position',[x y z],'String',pcfg.ele_txt{iGP}(i,:),'color',[1 1 1],'VerticalAlignment','baseline','FontSize',7); % pcfg.txt_col{iGP}
                end
            end
        end
    end
end

%% Add in lines if desired
if isfield(pcfg,'par_ele')
    
    %
    rsp_ele_nme = cat(1,pcfg.ele_lbl{:});
    rsp_ele_loc = cat(1,pcfg.rsp_ele{:});
    
    [~,iX,~] = unique(rsp_ele_nme);
    
    rsp_ele_nme = rsp_ele_nme(iX,:);
    rsp_ele_loc = rsp_ele_loc(iX,:);
    
    %
    for iGP = 1:numel(pcfg.par_ele)
        
        for iPR = 1:numel(pcfg.par_ele{iGP})
            if ~isempty(pcfg.rsp_ele{iGP})
                if mean(sub.vert(:, 1))<0
                     ele_1st_loc = find(strcmpi(rsp_ele_nme,pcfg.par_ele{iGP}{iPR}{1}));
                     ele_2nd_loc = find(strcmpi(rsp_ele_nme,pcfg.par_ele{iGP}{iPR}{2}));
            
                     line( [min(sub.vert(:, 1))+(min(sub.vert(:, 1))*.025) min(sub.vert(:, 1))+(min(sub.vert(:, 1))*.025)], ...
                           [rsp_ele_loc(ele_1st_loc,2) rsp_ele_loc(ele_2nd_loc,2)], ...
                           [rsp_ele_loc(ele_1st_loc,3) rsp_ele_loc(ele_2nd_loc,3)], ...
                           'Color',pcfg.par_ele_col{iGP},'LineWidth',pcfg.par_ele_lne(iGP))
                elseif mean(sub.vert(:, 1))>0
                     ele_1st_loc = find(strcmpi(rsp_ele_nme,pcfg.par_ele{iGP}{iPR}{1}));
                     ele_2nd_loc = find(strcmpi(rsp_ele_nme,pcfg.par_ele{iGP}{iPR}{2}));
            
                     line( [max(sub.vert(:, 1))+(max(sub.vert(:, 1))*.025) max(sub.vert(:, 1))+(max(sub.vert(:, 1))*.025)], ...
                     [rsp_ele_loc(ele_1st_loc,2) rsp_ele_loc(ele_2nd_loc,2)], ...
                     [rsp_ele_loc(ele_1st_loc,3) rsp_ele_loc(ele_2nd_loc,3)], ...
                     'Color',pcfg.par_ele_col{iGP},'LineWidth',pcfg.par_ele_lne(iGP))
                end
            end
        end
    end
end

%% Add Non-responsive electrodes
if isfield(pcfg,'rsp_ele') && isfield(pcfg,'nsl_ele')
if pcfg.non_ele
    c = repmat(size(colormap,1),30,30)/size(colormap,1)-(1/(size(colormap,1)*2));
    for i=1:size(pcfg.nsl_ele,1)
        plotSpheres(pcfg.nsl_ele(i,1),pcfg.nsl_ele(i,2),pcfg.nsl_ele(i,3),1,c,pcfg.sph_vew);
    end
    hold off;
end
end

%% Rotate View
if strcmpi(pcfg.sph,'lh')
    if strcmpi(pcfg.sph_vew,'lat')
        view(pcfg.axe_hnd,270, 0);
%         shadem('dull',[270 0])
    elseif strcmpi(pcfg.sph_vew,'ven')
        if ~isfield(pcfg,'sml_vew') || ~pcfg.sml_vew;
        view(pcfg.axe_hnd,90, 270);
        else
            view(pcfg.axe_hnd,270, 270);
        end
%         shadem('dull',[270 270])
    elseif strcmpi(pcfg.sph_vew,'med')
        view(pcfg.axe_hnd,90,0);
%         shadem('dull',[90 0])
    end
elseif strcmp(pcfg.sph,'rh')
    if strcmpi(pcfg.sph_vew,'lat')
        view(pcfg.axe_hnd,90, 0);
%         shadem('dull',[90 0])
    elseif strcmpi(pcfg.sph_vew,'ven')
        if ~isfield(pcfg,'sml_vew') || ~pcfg.sml_vew;
        view(pcfg.axe_hnd,270, 270);
        else
            view(pcfg.axe_hnd,90, 270);
        end
%         shadem('dull',[90 270])
    elseif strcmpi(pcfg.sph_vew,'med')
        view(pcfg.axe_hnd,270,0);
%         shadem('dull',[270 0])
    end
end
light('Position',campos);

set(pcfg.fig_hdl, 'color','white','InvertHardCopy', 'off');
axis tight;
% axis equal;
end
