clear; clc;

%% Alignment Data
% sEEG file containing HH task data 
fle_nme = '/space/seh8/3/halgdev/projects/bqrosen/MULTI_initial/seeg_UCSD/SD018/SD018_Day08.edf';

ntc_flt = [60,120,180,240]; % [60,120,180,240,300,360,420,480];

% optional: load only the header to inspect metadata
hdr_hld = edfread(fle_nme); 
smp_frq = hdr_hld.frequency(1);
trg_frq = 500;

% load referential sEEG data recorded from the macro contacts of the
% micromacroelectrode (Macro1-7) as well as the trigger data (TRIG)
targetSignals = {'Macro1', 'Macro2', 'Macro7', 'TRIG' };
tic
[header, dta_cln] = edfread(fle_nme, 'targetSignals', targetSignals);
toc

dta_cln = resample(dta_cln', trg_frq, smp_frq)';

% isolate trigger and sEEG data
dta_trg = dta_cln(3,:);
dta_cln = dta_cln(1:2,:);

% apply 60 Hz notch filter with harmonics up to Nyquist
for nf = ntc_flt
    fprintf('Performing Notch: %i Hz\n',nf)
    Wo = nf/(trg_frq/2);
    BW = Wo/35;
    [b,a] = iirnotch(Wo, BW);
    dta_cln = filtfilt(b,a,dta_cln')';
end
dta_cln_wpp = dta_cln(1,:)-dta_cln(2,:);

clear dta_cln

% 
itn_dir = '/space/seh9/1/halgdev/projects/micromacro/3_25_19/';

% Intan files (each 1 min long) containing HH task data -> From 8:45 to 9:30
itn_fle_nme = { 'SD018_Mon_190325_084425.rhd','SD018_Mon_190325_084525.rhd','SD018_Mon_190325_084625.rhd','SD018_Mon_190325_084725.rhd','SD018_Mon_190325_084825.rhd','SD018_Mon_190325_084925.rhd','SD018_Mon_190325_085025.rhd' ...
                'SD018_Mon_190325_085125.rhd','SD018_Mon_190325_085225.rhd','SD018_Mon_190325_085326.rhd','SD018_Mon_190325_085426.rhd','SD018_Mon_190325_085526.rhd','SD018_Mon_190325_085626.rhd','SD018_Mon_190325_085726.rhd','SD018_Mon_190325_085826.rhd','SD018_Mon_190325_085926.rhd','SD018_Mon_190325_090026.rhd' ...
                'SD018_Mon_190325_090126.rhd','SD018_Mon_190325_090226.rhd','SD018_Mon_190325_090326.rhd','SD018_Mon_190325_090426.rhd','SD018_Mon_190325_090526.rhd','SD018_Mon_190325_090626.rhd','SD018_Mon_190325_090726.rhd','SD018_Mon_190325_090826.rhd','SD018_Mon_190325_090926.rhd','SD018_Mon_190325_091026.rhd' ...
                'SD018_Mon_190325_091126.rhd','SD018_Mon_190325_091226.rhd','SD018_Mon_190325_091326.rhd','SD018_Mon_190325_091426.rhd','SD018_Mon_190325_091526.rhd','SD018_Mon_190325_091626.rhd','SD018_Mon_190325_091726.rhd','SD018_Mon_190325_091826.rhd' ...,'SD018_Mon_190325_091926.rhd','SD018_Mon_190325_092026.rhd' ...
                ...'SD018_Mon_190325_092126.rhd','SD018_Mon_190325_092226.rhd','SD018_Mon_190325_092326.rhd','SD018_Mon_190325_092427.rhd','SD018_Mon_190325_092527.rhd','SD018_Mon_190325_092627.rhd','SD018_Mon_190325_092727.rhd','SD018_Mon_190325_092827.rhd','SD018_Mon_190325_092927.rhd','SD018_Mon_190325_093027.rhd' 
                };

% itn_dir = '/space/seh9/1/halgdev/projects/micromacro/3_24_19/';
% 
% % Intan files (each 1 min long) containing HH task data -> From 8:45 to 9:30
% itn_fle_nme = { 'SD018_Sun1_uM__190324_132009.rhd','SD018_Sun1_uM__190324_132109.rhd','SD018_Sun1_uM__190324_132209.rhd','SD018_Sun1_uM__190324_132309.rhd','SD018_Sun1_uM__190324_132409.rhd','SD018_Sun1_uM__190324_132509.rhd','SD018_Sun1_uM__190324_132609.rhd','SD018_Sun1_uM__190324_132709.rhd','SD018_Sun1_uM__190324_132809.rhd','SD018_Sun1_uM__190324_132909.rhd' };

            
% load and concatenate Intan files from list
dta_itn = [];
fss_itn = [];
for file = 1:numel(itn_fle_nme)
    fprintf('%s %i complete\n',itn_fle_nme{file},file)
    dta_itn_tmp = read_Intan_RHD2000_file_nogui([itn_dir '/' itn_fle_nme{file}]);
    dta_itn = [dta_itn dta_itn_tmp.amplifier_data];
    fss_itn = [fss_itn dta_itn_tmp.frequency_parameters.amplifier_sample_rate];
end

clear dta_itn_tmp

% check whether sample rate for all Intan files is the same
if range(fss_itn) == 0
    fss_itn = fss_itn(1);
else
    disp('ERROR: Sample rates differ across Intan files.')
end

% downsample Intan data to clinical sampling rate (includes auto anti-aliasing)
dta_itn_dsm = resample(dta_itn', trg_frq, fss_itn)';

clear dta_itn

% 60 Hz notch filter with harmonics up to 480
for nf = ntc_flt
    fprintf('Performing Notch: %i Hz\n',nf)
    Wo = nf/(trg_frq/2);
    BW = Wo/35;
    [b,a] = iirnotch(Wo, BW);
    dta_itn_dsm = filtfilt(b,a,dta_itn_dsm')';
end

% align data using whole probe bipolar
% channel 16 = whole probe bipolar (Macro1 - Macro7)
itn_chn = 16;

% compute cross-correlation and find alignment time lag based on max value
xc = xcorr(dta_cln_wpp, dta_itn_dsm(itn_chn,:));
[~,lag] = max(xc);
offset = lag-size(dta_cln_wpp,2)+1; % use this value for the alignment

ttt = sort(xc);
for iST = 1:10
   ind_hld(iST) = find(ttt(iST)==xc); 
   off_set_ind(iST) = ind_hld(iST)-size(dta_cln_wpp,2)+1;
end

% plot total alignment
start = 1;
stop = size(dta_itn_dsm,2);
times = (1:stop-start+1)./trg_frq;

[start+offset stop+offset]% Alignment times

subplot(2,1,1); hold on
plot(times, dta_itn_dsm(itn_chn, start:stop), 'r')
plot(times, dta_cln_wpp(start+offset:stop+offset), '--b')
xlim([0 times(end)])
title('alignment of Intan and clinical recordings')
legend('intan whole probe bipolar', 'clinical whole probe bipolar', 'Location', 'southwest')
xlabel('time (s)')
ylabel('amplitude (\muV)')
subplot(2,1,2); hold on
plot(times, dta_trg(start+offset:stop+offset), 'k')

tightfig();
print(['/home/ekaestne/Downloads' '/' 'total_alignment.png'],'-dpng')
close all

% plot random alignment checks
start = randperm(size(dta_itn_dsm,2),10);
stop = start+(5*trg_frq);
times = (1:stop-start+1)./trg_frq;

for iP = 1:numel(start)
    subplot(2,5,iP); hold on
    plot(times, dta_itn_dsm(itn_chn, start(iP):stop(iP)), 'r')
    plot(times, dta_cln_wpp(start(iP)+offset:stop(iP)+offset), '--b')
    xlim([0 times(end)])
    title('alignment of Intan and clinical recordings')
%     legend('intan whole probe bipolar', 'clinical whole probe bipolar', 'Location', 'southwest')
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
end
tightfig();
print(['/home/ekaestne/Downloads' '/' 'random_alignment_3.png'],'-dpng')
close all

%% Get Burke Data
brk_dta = load('/space/seh8/6/halgdev/projects/bqrosen/MULTI/out/UCSD/SD018/SD018_SEEG_007.mat');

% Check overlap
tst_dta

%% Behavior 
log_fle_dir = '/space/seh9/1/halgdev/projects/micromacro/task_data/MST (Erik)';
log_fle_nme = { 'SD0018-MST_Rand1_PRACTICE.log' 'SD0018-MST_Rand2.log' 'SD0018_rand3-MST_Rand3.log'  }; % 'Austin_Test-MST_Rand1.log' 'Austin_fMRI_test-MST_Rand2.log' 'Austin_fMRI_test-MST_Rand3.log' 'Austin_fMRI_test-MST_Rand4.log'

%
clear out_dta

trl_col = 2;
stm_col = 3;
rsp_col = 4;
tme_col = 5;

nov_btn = 1001; %1001;
rep_btn = 1002; %1002;
sim_btn = 1003; %1003;

for iF = 1:numel(log_fle_nme)
    
    out_dta{iF} = zeros(3,4);
    out_tot_dta{iF} = zeros(3,4);
    out_prf_dta{iF} = zeros(3,4);
    
    log_fle = mmil_readtext([log_fle_dir '/' log_fle_nme{iF}],['\t']);
    fnd_beg = 0;
    beg_row = 1;
    while fnd_beg == 0
        if isnumeric(log_fle{beg_row,trl_col}) & log_fle{beg_row,trl_col}==3
            fnd_beg = 1;
        else
            beg_row = beg_row + 1;
        end
    end
    log_fle = log_fle(beg_row:end,:);
    
    nov_trl = 0;
    nov_trl_crr = 0;
    nov_trl_rep = 0;
    nov_trl_sim = 0;
    nov_trl_mss = 0;
    
    rep_trl = 0;
    rep_trl_crr = 0;
    rep_trl_sim = 0;
    rep_trl_nov = 0;
    rep_trl_mss = 0;
    
    sim_trl = 0;
    sim_trl_crr = 0;
    sim_trl_rep = 0;
    sim_trl_nov = 0;
    sim_trl_mss = 0;    
    
    for iR = 1:size(log_fle,1)
        if strcmpi(log_fle{iR,stm_col},'Picture') && log_fle{iR,rsp_col}>=1 && log_fle{iR,rsp_col}<100
            
            nov_trl = nov_trl + 1;
            
            nxt_trl = find(diff(cell2mat(log_fle(iR:end,trl_col)))==1);
            try nxt_trl = nxt_trl(1); rsp_trl = find(strcmpi(log_fle(iR:iR+nxt_trl,stm_col),'Response'));
                rsp_trl = rsp_trl(end); rsp_trl = log_fle{iR+rsp_trl-1,rsp_col}; catch rsp_trl = []; end
            
            if ~isempty(rsp_trl) && rsp_trl == nov_btn
                nov_trl_crr = nov_trl_crr + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == rep_btn
                nov_trl_rep = nov_trl_rep + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == sim_btn
                nov_trl_sim = nov_trl_sim + 1;
            elseif isempty(rsp_trl)
                nov_trl_mss = nov_trl_mss + 1; 
            end
            
        elseif strcmpi(log_fle{iR,stm_col},'Picture') && log_fle{iR,rsp_col}>=101 && log_fle{iR,rsp_col}<200
            
            rep_trl = rep_trl + 1;
            
            nxt_trl = find(diff(cell2mat(log_fle(iR:end,trl_col)))==1);
            try nxt_trl = nxt_trl(1); rsp_trl = find(strcmpi(log_fle(iR:iR+nxt_trl,stm_col),'Response'));
                rsp_trl = rsp_trl(end); rsp_trl = log_fle{iR+rsp_trl-1,rsp_col}; catch rsp_trl = []; end
            
            if ~isempty(rsp_trl) && rsp_trl == nov_btn
                rep_trl_nov = rep_trl_nov + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == rep_btn
                rep_trl_crr = rep_trl_crr + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == sim_btn
                rep_trl_sim = rep_trl_sim + 1;
            elseif isempty(rsp_trl)
                rep_trl_mss = rep_trl_mss + 1; 
            end
            
        elseif strcmpi(log_fle{iR,stm_col},'Picture') && log_fle{iR,rsp_col}>=201 && log_fle{iR,rsp_col}<300 && log_fle{iR,rsp_col}~=255
            
            nov_trl = nov_trl + 1;
            
            nxt_trl = find(diff(cell2mat(log_fle(iR:end,trl_col)))==1);
            try nxt_trl = nxt_trl(1); rsp_trl = find(strcmpi(log_fle(iR:iR+nxt_trl,stm_col),'Response'));
                rsp_trl = rsp_trl(end); rsp_trl = log_fle{iR+rsp_trl-1,rsp_col}; catch rsp_trl = []; end
            
            if ~isempty(rsp_trl) && rsp_trl == nov_btn
                nov_trl_crr = nov_trl_crr + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == rep_btn
                nov_trl_rep = nov_trl_rep + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == sim_btn
                nov_trl_sim = nov_trl_sim + 1;
            elseif isempty(rsp_trl)
                nov_trl_mss = nov_trl_mss + 1;  
            end
            
        elseif strcmpi(log_fle{iR,stm_col},'Picture') && log_fle{iR,rsp_col}>=301 && log_fle{iR,rsp_col}<400
            
            sim_trl = sim_trl + 1;
            
            nxt_trl = find(diff(cell2mat(log_fle(iR:end,trl_col)))==1);
            try nxt_trl = nxt_trl(1); rsp_trl = find(strcmpi(log_fle(iR:iR+nxt_trl,stm_col),'Response'));
                rsp_trl = rsp_trl(end); rsp_trl = log_fle{iR+rsp_trl-1,rsp_col}; catch rsp_trl = []; end
            
            if ~isempty(rsp_trl) && rsp_trl == nov_btn
                sim_trl_nov = sim_trl_nov + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == rep_btn
                sim_trl_rep = sim_trl_rep + 1;
            elseif ~isempty(rsp_trl) && rsp_trl == sim_btn
                sim_trl_crr = sim_trl_crr + 1;
            elseif isempty(rsp_trl)
                sim_trl_mss = sim_trl_mss + 1;    
            end
            
        elseif strcmpi(log_fle{iR,stm_col},'Picture') && log_fle{iR,rsp_col}>=401
            
        end
    end
    
    % Out Data
    out_dta{iF}(1,1) = nov_trl_crr / nov_trl; out_dta{iF}(1,2) = nov_trl_rep / nov_trl; out_dta{iF}(1,3) = nov_trl_sim / nov_trl; out_dta{iF}(1,4) = nov_trl_mss / nov_trl;
    out_dta{iF}(2,1) = rep_trl_nov / rep_trl; out_dta{iF}(2,2) = rep_trl_crr / rep_trl; out_dta{iF}(2,3) = rep_trl_sim / rep_trl; out_dta{iF}(2,4) = rep_trl_mss / rep_trl;
    out_dta{iF}(3,1) = sim_trl_nov / sim_trl; out_dta{iF}(3,2) = sim_trl_rep / sim_trl; out_dta{iF}(3,3) = sim_trl_crr / sim_trl; out_dta{iF}(3,4) = sim_trl_mss / sim_trl;
    
    % Total #'s
    out_tot_dta{iF}(1,1) = nov_trl; out_tot_dta{iF}(1,2) = nov_trl; out_tot_dta{iF}(1,3) = nov_trl; out_tot_dta{iF}(1,4) = nov_trl;
    out_tot_dta{iF}(2,1) = rep_trl; out_tot_dta{iF}(2,2) = rep_trl; out_tot_dta{iF}(2,3) = rep_trl; out_tot_dta{iF}(2,4) = rep_trl;
    out_tot_dta{iF}(3,1) = sim_trl; out_tot_dta{iF}(3,2) = sim_trl; out_tot_dta{iF}(3,3) = sim_trl; out_tot_dta{iF}(3,4) = sim_trl;
    
    % Correct #'s
    out_prf_dta{iF}(1,1) = nov_trl_crr ; out_prf_dta{iF}(1,2) = nov_trl_rep ; out_prf_dta{iF}(1,3) = nov_trl_sim ; out_prf_dta{iF}(1,4) = nov_trl_mss ;
    out_prf_dta{iF}(2,1) = rep_trl_nov ; out_prf_dta{iF}(2,2) = rep_trl_crr ; out_prf_dta{iF}(2,3) = rep_trl_sim ; out_prf_dta{iF}(2,4) = rep_trl_mss ;
    out_prf_dta{iF}(3,1) = sim_trl_nov ; out_prf_dta{iF}(3,2) = sim_trl_rep ; out_prf_dta{iF}(3,3) = sim_trl_crr ; out_prf_dta{iF}(3,4) = sim_trl_mss ;
    
end

tot_dta = zeros(size(out_dta{1}));
prf_dta = zeros(size(out_dta{1}));
for iF = 1:numel(out_dta)
    for iR = 1:size(out_dta{iF},1)
        for iC = 1:size(out_dta{iF},2)
            tot_dta(iR,iC) = tot_dta(iR,iC) + out_tot_dta{iF}(iR,iC);
            prf_dta(iR,iC) = prf_dta(iR,iC) + out_prf_dta{iF}(iR,iC);
        end
    end
end

prf_dta ./ tot_dta

%% Timing
log_fle_dir = '/space/seh9/1/halgdev/projects/micromacro/task_data/MST (Erik)';
log_fle_nme = { 'SD0018-MST_Rand1_PRACTICE.log' 'SD0018-MST_Rand2.log' 'SD0018_rand3-MST_Rand3.log'  }; % 'Austin_Test-MST_Rand1.log' 'Austin_fMRI_test-MST_Rand2.log' 'Austin_fMRI_test-MST_Rand3.log' 'Austin_fMRI_test-MST_Rand4.log'

for iF = 1:numel(log_fle_nme)
        
    log_fle = mmil_readtext([log_fle_dir '/' log_fle_nme{iF}],['\t']);
    fnd_beg = 0;
    beg_row = 1;
    while fnd_beg == 0
        if isnumeric(log_fle{beg_row,trl_col}) & log_fle{beg_row,trl_col}==3
            fnd_beg = 1;
        else
            beg_row = beg_row + 1;
        end
    end
    log_fle = log_fle(beg_row:end-1,:);
    
    log_fle(strcmpi(log_fle(:,stm_col),'Picture') & cell2mat(log_fle(:,rsp_col))==0,:) = [];
    log_fle(strcmpi(log_fle(:,stm_col),'Response'),:) = [];
    
    tme_isi{iF} = [ cell2mat(log_fle(:,rsp_col))  [0 ; diff(cell2mat(log_fle(:,tme_col)))] ];
    
end
