%% Function to read events & create file to epoch NYU iEEG data from a edf.mat format
% Function to epoch NYU iEEG data. 
% dat has 8 channels:each channel has the same trigger patten as the others
% Created by Jane Deng (4-13-15) q5deng@ucsd.edu

function [trl, evt] = ft_Utah_trialfunNY349SZ(data)

addpath(genpath('/home/qdeng/qdeng'));
addpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

data = load('/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel/NY349_SZ_edf.mat');
data_orig_name = cell2mat(fieldnames(data));
dat = getfield(data,data_orig_name);

thresh = 1 ;

Fs = 512;
pre = 1.5*Fs;
pos = 2.5*Fs;
%% Find the timings

tmpdat = dat.trial{1}(8,:);  
   
ind     = crossing(tmpdat,[],thresh);  
ind     = ind(1:2:length(ind));         % only keep left edge of pulse only the odd number
tmpdat2 = tmpdat(ind+2);                % select the peak of each pulses

trig_pos   =  [1 find(diff(ind)*1000/Fs > 100)+1];
samp   = ind(trig_pos);
% if sum(samp-pre < 1) || sum(samp+pos > numel(dat.time{1}));%delete the remote trial, and
%    trig_pos(samp+pos > numel(dat.time{1})) =[];
%    trig_pos(samp-pre < 1) =[];
%    samp(samp+pos > numel(dat.time{1})) =[];
%    samp(samp-pre < 1) =[];
% end


Nevents = numel(samp);
sel = zeros(1,numel(ind));
sel(trig_pos) = 1;
r = round(tmpdat2);
trigs = cell(1,Nevents);

%ft
for i = 1:Nevents;
    trigs{i} = r([trig_pos(i)]);
    trig_num(i) = numel(trigs{i});
end
% trigs{Nevents} = r(trig_pos(Nevents)+1:end);
% trig_num(Nevents) = numel(trigs{Nevents});

% if sum(ismember(trig_num,0)) >= 1 
%     ismember(trig_num,0);
% end
% 
%  z_ind = ismember(trig_num,0);
% % trigs(z_ind) = 1;              % throw away 0 member in trigger information
% % samp(z_ind) = [];               % throw away 0 member in sample
%  trig_num(z_ind) = 1;           % throw away 0 member in trigger number information 

%% Define Events
evt = zeros(1,numel(trig_num));

trigdata = mmil_readtext('/home/qdeng/qdeng/dataset/SZi-actual.csv');
trig =trigdata(:,4);
evt_val = cell2num(trig);

%% Adjust time to make it in the correct sampling frequency

% this one has been downsampled already.no need to run this step
%     dsfact = 3e4/500;
%     samp  = round(samp/dsfact);  
%     pre = round(pre/dsfact);
%     pos = round(pos/dsfact);


%% Make new trl file
mtx = zeros(length(evt),4);
for i=1:length(evt)
    mtx(i,1) = samp(i)-pre;
    mtx(i,2) = samp(i)+pos;
    mtx(i,3) = -pre;
    mtx(i,4) = evt(i);
end
tri_end = size(mtx,1);
for i = 1:tri_end 
    sen(i) = abs(mtx(2,1)-mtx(i,1));
end
sen2 = zeros(1,length(sen)-1);

for i = 2:(length(sen)-2)
    if sen(i) == 0
       sen2(i) = 1;
    end
    if sen(i+1)-sen(i) > 1000
       sen2(i) =sen(i);
       sen2(i+1)= sen(i+1);
    elseif (sen(i+1)-sen(i) < 1000) && sen2(i) ~= 0 && (sen(i+2)-sen(i) > 1000)
            sen2(i+2) = sen(i+2);      
    end  
    if sen(end)-sen(end-1) > 1100
       sen2(end) = sen(end);
    end
end
idd = find(~ismember(sen2,0));
if mtx(end,1)-mtx((end-1),1) >1000
   idd(end+1) = tri_end;
end
seq = sort(idd);

for i = 1:length(seq)
    vline(samp(seq(i)),'r')
end

%% make trigger matrix
trl=mtx(seq,:);
trl(:,4) = evt_val(1:size(trl,1));

rmpath(genpath('/home/qdeng/qdeng'));
rmpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

end


