%% Function to read events & create file to epoch Utah array data

% Modified from ft_NYU_iEEG_trialfun created by Erik Kaestner (4-4-14) ekaestne@ucsd.edu
% by Jane Deng q5deng@ucsd.edu (03/17/2015)
% use sampling frequency/rate = 500; Fs = 500; 
% has missing triggers. add the pre&post timing at the end part

function ft_Utah_trialfunMG63SA(data,trl)
addpath(genpath('/home/qdeng/qdeng'));
addpath(genpath('/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/MMPS/MMPS_236'));

data = load('/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel/MG63_SA__Utah.mat');
data_orig_name = cell2mat(fieldnames(data));
dat = getfield(data,data_orig_name);

thresh = 1.5e4 ;%15000

Fs = 500;
pre = 0.5*Fs;
pos = 1.5*Fs;
%% Find the timings
tmpdat = dat;    tmpdat([1:1.5e5])=0;   tmpdat([6.4e5:end])=0;
   
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

trigdata = mmil_readtext('/home/qdeng/qdeng/dataset/SAi-actual.csv');
trig =trigdata(:,3);
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

idx = [2 4 5 7 8 10 12 13 15 16 18 20 21 23 25 26 28 30 32:1:38 39 41 43 45 46 47 49 50 52 54 55 56 58 60 62 64 66 67 68 69 71 73 74 76 77 79 81 83 85 87 89 90 92 93 94 95 97 99];
idx2=[101 102 103 104 106 107 109 111 112 114 115 116 117 119 120 122 124 126 128 130 131 133 134 135 136 137 139 141 143 145 146 148 150 151 153 155 156 157 158 159 161 163 165 166 168 170 171:174 176 178 180 182 183 184 186 188 190 192 194 196 197 199 201];    
idx3 =[203 205 207 208 209 211 212 214 215 217 218 220 222 223 225 227 228 230 231 233 235 236 238 240 242 243 245 246 247 249 251 252];
idx4 =[254 256 257 259 261 263 265 266 268 269 271 272 274 276 277 279 281 282 283 285 287 288 290 291 292];
idx5 =[294 295 296 298 300 301 302 303 304  306 307 309 311 312 313 315 317 318 319 320 322 324 325 326 328 330 332 333 334 335 337 338 340 341 343 345 347 349 350];
idx6 = [352 353 355 356 358 360 361 363 365 366 368 369 371 372 373 375 377 378 380 382 384 386 387 389 391 392 393 395 397 399 400];
idx7 = [401 403 404 406 408 409 410 412 413 415 416 417 418 419 421 422 424 426 427 428 430 431 433 435 437 438 440 442 444 445 446 448 450];
idx8 = [452 453 455 456 457 458 459 460 461 463 465 466 468 470 471 473 474 476 478 480 482 484 485 487 489 491 492 494 495 496 498 499];

id = [idx idx2 idx3 idx4 idx5 idx6 idx7 idx8];

trl=mtx(id,:);
trl(:,4) = evt_val(1:320);

end
