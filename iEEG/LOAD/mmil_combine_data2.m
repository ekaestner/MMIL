function dat = mmil_combine_data2(fcfg,dat)

num = unique(fcfg.cmb);
frs_num = [];
scd_num = [];
for iN = 1:numel(num)
    frs_num = [frs_num repmat(find(num(iN)==fcfg.cmb,1),1,sum(fcfg.cmb==num(iN))-1)];
    loc = find(fcfg.cmb==num(iN));
    scd_num = [scd_num loc(2:end)];
end

if ~isempty(frs_num)
if numel(intersect(dat.(dat.data_name{frs_num(1)}).label,dat.(dat.data_name{scd_num(1)}).label))==numel(dat.(dat.data_name{frs_num(1)}).label)
    cfg           = [];
    cfg.data_name = [frs_num ; scd_num];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'trials';
    dat = ft_func(@ft_appenddata,cfg,dat);
else
    cfg           = [];
    cfg.data_name = [frs_num ; scd_num];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'channels';
    dat = ft_func(@ft_appenddata,cfg,dat);
end

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name =  scd_num;
dat = ft_func([],cfg,dat);
end

mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd',['[' strrep(num2str(ones(1,numel(num))),'  ',' ') ']']);


end