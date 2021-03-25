
%% Function to read events & create file to epoch Utah array data

% Modified from ft_NYU_iEEG_trialfun created by Erik Kaestner (4-4-14) ekaestne@ucsd.edu
% by Jane Deng q5deng@ucsd.edu (03/18/2015)
% use sampling frequency/rate = 500; Fs = 500; 
function ft_Utah_trialfunMG44SZ(data)

addpath(genpath('/home/qdeng/qdeng'));
addpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

data = load('/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel/MG44_SZ__Utah.mat');
data_orig_name = cell2mat(fieldnames(data));
dat = getfield(data,data_orig_name);

thresh = 1.5e4 ;%15000

Fs = 500;
pre = 0.5*Fs;
pos = 1.5*Fs;
%% Find the timings
tmpdat = dat;    tmpdat([1:9.5e5])=0;   tmpdat([1.5e6:end])=0;
   
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
for i = 1:Nevents-1;
    trigs{i} = r([trig_pos(i)+1:trig_pos(i+1)-1]);
    trig_num(i) = numel(trigs{i});
end
trigs{Nevents} = r(trig_pos(Nevents)+1:end);
trig_num(Nevents) = numel(trigs{Nevents});
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
evt = cell2num(trig);

%% Adjust time to make it in the correct sampling frequency

% this one has been downsampled already.no need to run this step
%     dsfact = 3e4/500;
%     samp  = round(samp/dsfact);  
%     pre = round(pre/dsfact);
%     pos = round(pos/dsfact);


%% Make new trl file
trl = zeros(length(evt),4);
for i=1:length(evt)
    trl(i,1) = samp(i)-pre;
    trl(i,2) = samp(i)+pos;
    trl(i,3) = -pre;
    trl(i,4) = evt(i);
end



end

