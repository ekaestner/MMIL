% cfg.measure = string of which rejections to use (default is all)
% cfg.apply   = string of which method to use (either ieeg or meg)

function ft_dat = ft_apply_rejection(cfg,ft_dat)

if ~isfield(cfg,'measure'); cfg.measure = 'all'; end

mea_num = find(strcmpi(cfg.measure,ft_dat.cfg.badtrl_type));

if isfield(ft_dat,'trial')
    dat     = permute(cat(3,ft_dat.trial{:}),[3 1 2]);
end

if strcmpi(cfg.apply,'ieeg') && isfield(ft_dat,'trial')
    
    [ib,ic] = find(squeeze(ft_dat.cfg.badtrl(mea_num,:,:)));
    rmv_ind = sub2ind(size(dat),cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ic,repmat(numel(ft_dat.time{1}),numel(ic),1),'UniformOutput',0)')',cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ib,repmat(numel(ft_dat.time{1}),numel(ib),1),'UniformOutput',0)')',repmat(1:numel(ft_dat.time{1}),1,numel(ib))');
    dat(rmv_ind) = nan;
    
elseif strcmpi(cfg.apply,'ieeg') && isfield(ft_dat,'power') || isfield(ft_dat,'phase') || isfield(ft_dat,'fourierspctrm')
    
    if isfield(ft_dat,'power'); dat_fld = 'power'; elseif isfield(ft_dat,'phase'); dat_fld = 'phase'; elseif isfield(ft_dat,'fourierspctrm'); dat_fld = 'fourierspctrm'; end
        
    for iCH = 1:numel(ft_dat.label)
        fprintf(['Channel' num2str(iCH) '\n'])
        
        [ib,ic] = find(squeeze(ft_dat.cfg.badtrl(mea_num,iCH,:)));
        
        tme = repmat(1:numel(ft_dat.time), size(ft_dat.(dat_fld),3), 1); tme = tme(:);
        rmv_ind = sub2ind(size(ft_dat.(dat_fld)), ...
            cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ib,repmat(numel(ft_dat.time)*size(ft_dat.(dat_fld),3),numel(ib),1),'UniformOutput',0)')', ...
            cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ic,repmat(numel(ft_dat.time)*size(ft_dat.(dat_fld),3),numel(ic),1),'UniformOutput',0)')', ...
            repmat(1:size(ft_dat.(dat_fld),3),1,numel(ib)*numel(ft_dat.time))', ...cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),[1:size(ft_dat.(dat_fld),3)]',repmat(numel(ft_dat.time),size(ft_dat.(dat_fld),3),1)*numel(ib),'UniformOutput',0)')', ...
            repmat(tme',1,numel(ib))');
        
        ft_dat.(dat_fld)(rmv_ind) = nan;
    end
    
elseif strcmpi(cfg.apply,'meg')
    
    if isfield(ft_dat,'trial')
        [~,ic] = find(squeeze(ft_dat.cfg.badtrl(mea_num,:,:)));
        dat(unique(ic),:,:) = nan;
    elseif isfield(ft_dat,'powspctrm')
        [~,ic] = find(squeeze(ft_dat.cfg.badtrl(mea_num,:,:)));
        ft_dat.powspctrm(unique(ic),:,:,:) = nan;
    elseif isfield(ft_dat,'fourierspctrm')
        [~,ic] = find(squeeze(ft_dat.cfg.badtrl(mea_num,:,:)));
        ft_dat.fourierspctrm(unique(ic),:,:,:) = nan;
    end
end

if isfield(ft_dat,'trial')
    ft_dat.trial  = cellfun(@squeeze,num2cell(dat,[2 3]),'uni',0);
end

end


% squeeze(tmp_frq_dat(64,64,10,:))
% phs = abs(squeeze(tmp_frq_dat(64,64,:,:)));
% pow = angle(squeeze(tmp_frq_dat(64,64,:,:)));




