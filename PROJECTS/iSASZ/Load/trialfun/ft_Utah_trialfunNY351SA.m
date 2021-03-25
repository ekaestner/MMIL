function [trl, evt] = ft_Utah_trialfunNY351SA(data)

addpath(genpath('/home/qdeng/qdeng'));
addpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

%% Load the data
data = load('/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel/NY351_SA_edf.mat');
data_orig_name = cell2mat(fieldnames(data));
dat = getfield(data,data_orig_name);   % channel 5,6,8 should be removed

%% Define parameters 
Fs = 512;
pre   = round(1.5*Fs);
pos   = round(2.5*Fs);
t = dat.time{1};

thresh = 2e6;
nsamp = numel(t);

%% Find the timings
ichan=[1 2 3 4 7];  % valid channels for this data
ntrigchan = numel(ichan);
[crsind{1:ntrigchan}] = deal([]);

for k = 1:ntrigchan
    rmdat = dat.trial{1}(ichan(k),:);
    rmdat([1:1.12e5])=0;    % remove unuseful trigger
    rmdat([3.6e5:4e5])=0;
    dat.trial{1}(ichan(k),:) = rmdat;
    clear rmdat
end
        
for k = 1:ntrigchan
    tmpdat  = dat.trial{1}(ichan(k),:);
    ind     = crossing(tmpdat,[],thresh);
    % only keep left edge of pulse
    sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
    ind     = ind(sel);
    crsind{k} = ind;
    clear ind sel tmpdat
end

% combine all the triggers 
sampid = unique([crsind{1} crsind{2} crsind{3} crsind{4} crsind{5}]);

% dif_smp = diff(sampid);
% dif_smp = find(dif_smp<20);
% sampid(dif_smp) = [];

tri_end = numel(sampid);
tri_blc = find(diff(sampid)>3000);

tri_blc_1 = [1:1:tri_blc(1)];
tri_blc_2 = [(tri_blc(1)+1):1:tri_blc(2)];
tri_blc_3 = [(tri_blc(2)+1):1:tri_end];

%  for block 1
smp_1st = sampid(tri_blc_1);

gss_1st = sampid(tri_blc_1(2));
for iGS = 2:160
    gss_1st(iGS) = gss_1st(iGS-1) + round(512*2.2);
    [~,ejk] = min(abs(sampid-gss_1st(iGS)));
    gss_1st(iGS) = sampid(ejk);
end

% for block 2
smp_2nd = sampid(tri_blc_2);

gss_2nd = sampid(tri_blc_2(3));
for iGS = 2:160
    gss_2nd(iGS) = gss_2nd(iGS-1) + round(512*2.2);
    [~,ejk] = min(abs(sampid-gss_2nd(iGS)));
    gss_2nd(iGS) = sampid(ejk);
end

% for block 3
smp_3rd = sampid(tri_blc_3);

gss_3rd = sampid(tri_blc_3(3));
for iGS = 2:160
    gss_3rd(iGS) = gss_3rd(iGS-1) + round(512*2.2);
    [~,ejk] = min(abs(sampid-gss_3rd(iGS)));
    gss_3rd(iGS) = sampid(ejk);
end

%% make trigger matrix
trl(:,1) = [gss_1st' ; gss_2nd' ; gss_3rd'] - pre;
trl(:,2) = [gss_1st' ; gss_2nd' ; gss_3rd'] + pos;
trl(:,3) = repmat(-pre,size(trl,1),1);

trigdata = mmil_readtext('/home/qdeng/qdeng/dataset/SAi-actual.csv');
trig =trigdata(:,3);
evt_val = cell2num(trig);
trl(:,4) = evt_val(1:size(trl,1));
evt      = trl(:,4);

rmpath(genpath('/home/qdeng/qdeng'));
rmpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

end