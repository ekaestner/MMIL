%% SL Analysis Script
clear; clc;

sbj_num = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting on Subject %s \n\n\n',subj)

indir       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/indir/' subj]);  
cln_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cln_dir/' subj]); 
end_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/end_dir/' subj]); 
tsk         = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/tsk/' subj]);     
rmv_chn     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/rmv_chn/' subj]),'uni',0)';
bse_frq     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/bse_frq/' subj]),'uni',0)';
trl_fun = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/trialfun/' subj]);

chn_nse_grp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' subj]);
tmp = strfind(chn_nse_grp,';');
if any(~isempty(tmp{:}))
    for i = 1:numel(chn_nse_grp)
        chn_nse_grp_dat_nme(i) = str2num(chn_nse_grp{i}(tmp{1}+1:end))';
        chn_nse_grp{i} = str2num(chn_nse_grp{i}(1:tmp{1}-1))'; 
    end
else
    chn_nse_grp = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' subj]),'uni',0)';
end

% In Script Switches 
plt_spc = 0;

%% Data Paths
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_freq/';

for iIN = 1:numel(indir)
    inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
end

inpath = [];
for iIN = 1:numel(indir)
    inpath = [inpath inpath_holder.(tsk{iIN})];
end
inpath = sort(inpath);

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sll_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0);

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
if any(~cellfun(@isempty,strfind(end_dir,'set')) | ~cellfun(@isempty,strfind(end_dir,'edf'))); cfg.continuous = 'yes'; else cfg.continuous = 'no'; end
cfg.dataset    = inpath;
sll_dat        = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Remove Unimportant Channels - EJK
rmv_chn_uni = find(~cellfun(@isempty,strfind(sll_dat.(sll_dat.data_name{1}).label,'DC')))';
rmv_chn_uni = [rmv_chn_uni find(~cellfun(@isempty,strfind(sll_dat.(sll_dat.data_name{1}).label,'EKG')))'];

if ~isempty(rmv_chn_uni)
    cfg = []; 
    cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{1}).label(rmv_chn_uni))'];
    sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'before' '/' subj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Remove identified Bad Channels
if ~isempty(rmv_chn)
    if numel(rmv_chn) == 1
    cfg = []; 
    cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{1}).label(rmv_chn))'];
    sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
    else
        for iRM = 1:numel(rmv_chn)
            cfg = [];
            cfg.data_name = iRM;
            cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{iRM}).label(rmv_chn{iRM}))'];
            sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
        end
    end
end

%% Remove noise through removing common noise
% if ~isvar('chn_nse_grp_dat_nme'); chn_nse_grp_dat_nme = repmat({1:numel(tsk)},1,numel(chn_nse_grp)); end
for i = 1:numel(chn_nse_grp)
    if ~isempty(chn_nse_grp{i})
        cfg             = [];
        cfg.data_name   = chn_nse_grp_dat_nme{1};
        cfg.chn_nse_grp = chn_nse_grp(i);
        sll_dat = ft_func(@ft_remove_common_noise,cfg,sll_dat);
    end
end

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sll_dat      = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'after' '/' subj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Epoch Data
if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'])
    cfg = [];
    cfg.trialfun = trl_fun{1};
    cfg.dataset  = inpath;
    cfg.specific = {'dataset' ; 1:numel(inpath)};
    cfg.minduration = 0.500;
    cfg.pre         = 1.5;
    cfg.pos         = 2.5;
    cfg.ignore      = 1;
    cfg.trg_chn     = 130:136;
    trl = ft_func(@ft_definetrial,cfg,trl);  
    
    trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end
    
    % - EJK Grab timings & events here
    
    for iTS = 1:numel(trl_hld)
       
        blk_ind = regexpi(inpath{iTS},'\d_\d');
        
        blk = num2str(inpath{iTS}(blk_ind));
        run = num2str(inpath{iTS}(blk_ind+2)); 
        
        % - EJK Check timings & add events here
        
        lst = mmil_readtext(['/home/ekaestne/Desktop/N439_SL_Behav/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv']);
        trl_hld{iTS}(:,4) = cell2mat(lst(:,3));
        
        % - EJK add auditory specific timings here
        tmp_trl = [trl_hld{iTS}(:,1:2)+(round(0.450*sll_dat.(sll_dat.data_name{iTS}).fsample)) trl_hld{iTS}(:,3) trl_hld{iTS}(:,4)+10];
        trl_hld{iTS} = [trl_hld{iTS} ; tmp_trl];
    end
    
    tmp = input('Satisfied with defined trials?: '); clear tmp
    
    save(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
else
    load(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
end

cfg = [];
cfg.specific  = {'trl';[1:numel(inpath) 1:numel(inpath) 1:numel(inpath)]};
cfg.trl       = trl_hld;
sll_dat = ft_func(@ft_redefinetrial,cfg,sll_dat);

%% Combine Data
% combine trials
num_trl       = numel(inpath);
dat_cmb_ind = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2; ...
        3:2:num_trl 4:2:num_trl];
    
for iCM = 1:numel(dat_cmb_ind(2,:)) % 
    cfg           = [];
    cfg.data_name = [dat_cmb_ind(1,iCM) ; dat_cmb_ind(2,iCM)];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'trials';
    sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);
        
end

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

% combine Clin1 & Clin2
dat_cmb_ind = [1:2:numel(sll_dat.data_name) ; ...
    2:2:numel(sll_dat.data_name)];
    
for iCM = 1:numel(dat_cmb_ind(2,:)) % 
    cfg           = [];
    cfg.data_name = [dat_cmb_ind(1,iCM) ; dat_cmb_ind(2,iCM)];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'channels';
    sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);
       
end

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

%% Grab previous work
load([outpath '/' 'data_information' '/' subj '/' 'badtrl_hld.mat']);  % badtrl_hld, badtrl_type_hld
load([outpath '/' 'data_information' '/' subj '/' 'alt_eve_hld.mat']); % alt_eve_hld
load([outpath '/' 'data_information' '/' subj '/' 'alt_lab_hld.mat']); % alt_lab_hld

sll_dat.(sll_dat.data_name{1}).cfg.badtrl      = badtrl_hld;
sll_dat.(sll_dat.data_name{1}).cfg.badtrl_type = badtrl_type_hld;
sll_dat.(sll_dat.data_name{1}).cfg.alt_eve     = alt_eve_hld;
sll_dat.(sll_dat.data_name{1}).cfg.alt_lab     = alt_lab_hld;

%% iterate through wavelets
foi =  [2:12,14:2:24,25,30:5:55,70:10:170];    %frequency of interest
sf  =   [1,1,repmat(2,1,9),repmat(3,1,6),repmat(5,1,7),repmat(10,1,11)]; %specific frequency
gwidth = ones(size(foi))*pi; %wavelet

eve_cmp     = {[1 2 3 4] 1 2 3 4};
eve_lab     = {'over all' 'match' 'mismatch' 'visual noise' 'auditory noise'};
alt_eve_nme = 'trialinfo';
sve_lab     = 'repair_macro_ejk1st_meso_label';

for iCH = 1:numel(sll_dat.(sll_dat.data_name{1}).label)

    cfg = [];
    cfg.channel   = sll_dat.(sll_dat.data_name{1}).label(iCH); % 182:198
    itc_dat       = ft_func(@ft_preprocessing,cfg,sll_dat);
    
    cfg = [];
    cfg.data_name  = 1;
    cfg.method     = 'wavelet';
    cfg.output     = 'fourier';
    cfg.foi        = foi;	      % frequency of interest
    cfg.sf         = sf;        % specific frequency
    cfg.width      = foi./sf;
    cfg.gwidth     = gwidth;
    cfg.toi        = itc_dat.(itc_dat.data_name{1}).time{1}; %cellfun(@(x) sll_dat.(x).time{1},sll_dat.data_name,'uni',0);
    cfg.keeptrials = 'yes';
    itc_dat        = ft_func(@ft_freqanalysis,cfg,itc_dat);
    
    % Trim Timing
    cfg = [];
    cfg.data_name = 1;
    cfg.latency   = [-0.5 1.5];
    itc_dat       = ft_func(@ft_selectdata,cfg,itc_dat);
        
    % grab events %%%%%%%
    % remove bad trials & bad events   
    eve_hld = cell(1,numel(eve_cmp));
    for iEH = 1:numel(eve_cmp)
        for iEN = 1:numel(eve_cmp{iEH})
            eve_hld{iEH} = [eve_hld{iEH} ; find(sll_dat.(sll_dat.data_name{1}).cfg.alt_eve.(alt_eve_nme)==eve_cmp{iEH}(iEN))];
        end
    end
    
    figure('Visible','off')
    for iEV = 1:numel(eve_hld)
        
        % remove bad trials %%%%%%%
        bad_trl = find(squeeze(sll_dat.(sll_dat.data_name{1}).cfg.badtrl(end,iCH,:)));
        
        eve_inc = setdiff(eve_hld{iEV},bad_trl);
        
        % calculate %%%%%%%
        % Calculate inter-trial phase coherence (itpc)
        F = itc_dat.(itc_dat.data_name{1}).fourierspctrm(eve_inc,:,:,:);   % copy the Fourier spectrum
        N = size(F,1);           % number of trials
        
        itc_eeg = sum(F,1) ./ sum(abs(F),1);
        itc_eeg = abs(itc_eeg);
        itc_eeg = squeeze(itc_eeg);
        
        % Calculate Rayleigh
        F_ang = angle(F);
        pval = zeros(size(F,3),size(F,4));
        z    = zeros(size(F,3),size(F,4));
        for iF = 1:size(F,3)
            for iT = 1:size(F,4)
                [pval(iF,iT),z(iF,iT)] = circ_rtest(squeeze(F_ang(:,1,iF,iT)));
            end
        end
        
        % Add FDR correction
        [~,fdrmask] = fdr_cor([pval(:)]',0.01);
        fdrmask     = vec2mat(fdrmask,size(pval,1))';
        pval_fix  = zeros(size(pval,1),size(pval,2));
        pval_fix(logical(fdrmask)) = 1-pval(logical(fdrmask));
        
        % plot %%%%%%%
        sub_plt(:,1) = repmat(0:1/3:1-1/3,1,numel(eve_hld));
        sub_plt(:,2) = rude(repmat(3,1,numel(1-1/numel(eve_hld):-1/numel(eve_hld):0)),1-1/numel(eve_hld):-1/numel(eve_hld):0);
        sub_plt(:,3) = repmat(1/3,1,3*numel(eve_hld));
        sub_plt(:,4) = repmat(1/numel(eve_hld),1,3*numel(eve_hld));
        
        axes('OuterPosition',sub_plt(1+(3*(iEV-1)),:))
        surf(itc_dat.(itc_dat.data_name{1}).time, itc_dat.(itc_dat.data_name{1}).freq(1:20),itc_eeg(1:20,:),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[0 1],'xlim',[itc_dat.(itc_dat.data_name{1}).time(1) itc_dat.(itc_dat.data_name{1}).time(end)]);
        ylabel(eve_lab{iEV});
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        axes('OuterPosition',sub_plt(2+(3*(iEV-1)),:))
        surf(itc_dat.(itc_dat.data_name{1}).time, itc_dat.(itc_dat.data_name{1}).freq(1:20),z(1:20,:),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[0 50],'xlim',[itc_dat.(itc_dat.data_name{1}).time(1) itc_dat.(itc_dat.data_name{1}).time(end)]);
        line('XData',[0 0],'YData',ylim,'ZData',[50 50],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[50 50],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        axes('OuterPosition',sub_plt(3+(3*(iEV-1)),:))
        surf(itc_dat.(itc_dat.data_name{1}).time, itc_dat.(itc_dat.data_name{1}).freq(1:20),pval_fix(1:20,:),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[0.98 1],'xlim',[itc_dat.(itc_dat.data_name{1}).time(1) itc_dat.(itc_dat.data_name{1}).time(end)]);
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
    end
    
    ppos = get(gcf,'PaperPosition');
    su = get(gcf,'Units');
    pu = get(gcf,'PaperUnits');
    set(gcf,'Units',pu);
    spos = get(gcf,'Position');
    set(gcf,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
    set(gcf,'Units',su)
    print(gcf,sprintf('%s/initial_plot/%s/itc/%s.png',outpath,subj,sll_dat.(sll_dat.data_name{1}).cfg.alt_lab.(sve_lab){iCH}),'-dpng','-r300') % - EJK implement chan-by-chan
    
end




















