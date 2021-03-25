function abcd_pqa_fbirn(vol4D, meta, output, fwhm)

%   This code is adapted from Gary Glover fBIRN phantom qa.m matlab code
%   http://www.birncommunity.org/tools-catalog/function-birn-stability-phantom-qa-procedures/
%   It also contains contributions from code provided by Thomas Liu (Center
%   for Functional MRI, UCSD)
%
%   Contact author:
%   Jose Teruel. jteruelantolin@ucsd.edu
%   =======================================================================
%
%   vol4D = 4D matrix containing all phantom volumes over time
%   meta = struct containing TR, imageFreq, transmitGain and ReceiverGain
%   from dicom header.
%   output = folder to export results
%
%   =======================================================================
version = '0.0.1'; %Script version to be included as output of json file 

tr = meta.TR;
imageFreq = meta.imageFreq;
gains = [meta.transmitGain, meta.aRecGain];

path = output;

TR = tr/1000;
vol_dims = size(vol4D);
NPIX = vol_dims(1);


vol4D = vol4D(:,:,:,3:end); % Skip two first volumes

    
selectedSlice = ceil(size(vol4D,3)/2);

data = zeros(vol_dims(1),vol_dims(2),vol_dims(4));
data = permute(vol4D(:,:,selectedSlice,:),[1 2 4 3]);

[spatial_drift] = motion_estimates(vol4D,meta,path);


if(NPIX == 128)
  R = 30;                   % ROI width
elseif(NPIX == 104)
  R = 24;
else
  R = 15;                   % probably 64x64
end                         %ROI size fixed

npo2 = NPIX/2;              %half FOV 
ro2 = fix(R/2);             %half ROI size
X1 = npo2 - ro2;            %Beginging of masked image in X
X2 = X1 + R - 1;            %End of masked image X
Y1 = X1;                    %Beginning of masked image in Y
Y2 = X2;                    %Beginning of masked image in Y
r1 = 1;                     %First ROI voxel
r2 = R;                     %Last ROI voxel

%  set up input and ROI mask

mask = ones(R);                 %ROI mask
npx_roi = sum(sum(mask));       %Total number of voxels in ROI


i1=1;                           % initial time frame
i2=size(data,3);                % final time frame

N = i2 - i1 + 1;                % num time frames
M = r2 - r1 + 1;                % num ROI's
weisskoff = zeros(N,M);         % matrix with frames by ROIs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  begin loop through images

Iodd = zeros(NPIX^2,1);     %Images in odd positions as vector
Ieven = zeros(NPIX^2,1);    %Images in even positions as vector
Sy = zeros(NPIX^2,1);
Syy = zeros(NPIX^2,1);      % Sum t=1,nFrame; It*It
Syt = zeros(NPIX^2,1);      % Sum t=1,nFrame; t*It
St = 0;                     % Sum t=1,nFrame; t
Stt = 0;                    % Sum t=1,nFrame; t^2 
S0 = 0;                     % Time counter 1,2,3... nFrame
img = zeros(NPIX);


if(mod(N,2)==1) 
    even_tf_flag=0;
else
    even_tf_flag=1;
end

for j = i1:i2                       %For each tFrame
    workingSlice = data(:,:,j);     %Get appropiate slice
    I = workingSlice(:);            %Image slice as vector
    clear workingSlice;
    if(mod(j,2)==1)
        if even_tf_flag
            Iodd = Iodd + I; %Add odd images together
        else
            even_tf_flag=1;
        end
    else
      Ieven = Ieven + I;        %Add even images together
    end
    Sy = Sy + I;             % Add all time frimes
    Syt = Syt + I*j;         % Calc Sum t=1,nFrame; t*It
    Syy = Syy + I.*I;        % Calc Sum t=1,nFrame; It*It
    S0 = S0 + 1;             % Update counter
    St = St + j;             % Calc Sum t=1,nFrame; t
    Stt = Stt + j*j;         % Calc Sum t=1,nFrame; t^2 
    img(:) = I;
    sub = img(X1:X2,Y1:Y2);         % masked part of image phantom
    avg_signal_roi(S0) = sum(sum(sub))/npx_roi;    % average signal intensity trough ROI
    for r = r1:r2                   % each roi size
      ro2 = fix(r/2);
      x1 = npo2 - ro2;
      x2 = x1 + r - 1;
      sub = img(x1:x2,x1:x2);
      weisskoff(j-i1+1, r) = mean(sub(:));   %mean of masked image through time and roi size
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  write out diff image

Isub = Iodd - Ieven;        %Isub in vector form
img(:) = Isub;              %Isub in 2D form
sub = img(X1:X2,Y1:Y2);     %Isub within ROI
varI = var(sub(:));         %Variance of Isub within ROI only

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
imagesc(img/(fix(N/2))); 
colormap(gray);
colorbar;
file = 'Idiff_noise.tif';
filename = fullfile(path, file);
saveas(gcf, filename); % save the figure in the current folder as a .tif file
close; % closes current figure

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  write out ave image

Iave = Sy/N;                %Iave in vector form
img(:) = Iave;              %Iave in 2D form
sub = img(X1:X2,Y1:Y2);     %Iave with ROI masked
meanI = mean(sub(:));       %Mean of Iave within ROI

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
imagesc(img); % draw a simple figure containing a sine wave, title, etc.
colormap(gray);
colorbar;
file = 'Iave.tif';
filename = fullfile(path, file);
saveas(gcf, filename) % save the figure in the current folder as a .tif file
close; % closes current figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find trend line at + b

D = (Stt*S0 - St*St);
a = (Syt*S0 - St*Sy)/D;
b = (Stt*Sy - St*Syt)/D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sd image

Var = Syy + a.*a*Stt +b.*b*S0 + 2*a.*b*St - 2*a.*Syt - 2*b.*Sy;
Isd = sqrt(Var/(N-1));
img(:) = Isd;   

% print as 10*Isd if you want. Not now

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
imagesc(img); % draw a simple figure containing a sine wave, title, etc.
colormap(gray);
colorbar;
file = 'Isd.tif';
filename = fullfile(path, file);
saveas(gcf, filename) % save the figure in the current folder as a .tif file
close; % closes current figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sfnr image

sfnr = Iave./(Isd + eps);     %sfnr image in vector
img(:) = sfnr;                %sfnr image in 2D
sub = img(X1:X2,Y1:Y2);       %sfnr image within ROI
sfnrI = mean(sub(:));         %mean sfnr value within ROI

% print as 10*sfnr if you want. Not now
figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
imagesc(img); % draw a simple figure containing a sine wave, title, etc.
colormap(gray);
colorbar;
file = 'sfnr.tif';
filename = fullfile(path, file);
saveas(gcf, filename) % save the figure in the current folder as a .tif file
close; % closes current figure



snr = meanI/sqrt(varI/N); %SNR value
%fprintf('\nmean, SNR, SFNR = %5.1f  %5.1f  %5.1f\n', meanI, snr, sfnrI);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Do fluctation analysis

% First plot: Percent fluctuation and percent drift

x=(1:N);
p=polyfit(x,avg_signal_roi,2);
yfit = polyval(p, x);
y = avg_signal_roi - yfit;
time_minutes = x*TR/60;

m=mean(avg_signal_roi);
sd=std(y);
drift = ((yfit(N)-yfit(1))/m);
drift_per_minute = ((yfit(N)-yfit(1))/m)/time_minutes(end); %Drift using trend line. Last intensity-initial divided by mean of trend line (averaged by minute)
maxabsdrift = (max(avg_signal_roi) - min(avg_signal_roi))/m;

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
subplot(2,1,1);
plot(time_minutes,avg_signal_roi,time_minutes,yfit);
xlabel('Time (minutes)');
ylabel('Raw signal');
grid
title(sprintf('mean = %5.1f // SNR = %5.1f // SFNR = %5.1f \n drift (%% of mean) = %5.2f // drift per minute (%% of mean per minute) = %5.2f \n max absolute drift (%% of mean) = %5.2f',  meanI, snr, sfnrI, 100*drift, 100*drift_per_minute, 100*maxabsdrift));

subplot(2,1,2);
plot(time_minutes,y);
xlabel('Time (minutes)');
ylabel('Raw signal (detrended)');
grid
title(sprintf('sd = %5.1f // detrended fluctuation (%% of mean) = %5.2f', sd, 100*sd/m));


%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''off'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
resultFile = 'SignalTime_Graph.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .tif file
close;

close;

%fprintf('std, percent fluc, drift = %5.2f  %6.2f %6.2f \n', sd, 100*sd/m, 100*drift);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second plot frequency analysis of fluctuations (detrended)

z = fft(y);
fs = 1/TR;
nf = N/2+1;
f = 0.5*(1:nf)*fs/nf;
figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
plot(f, abs(z(1:fix(nf))));
grid
ylabel('spectrum');
xlabel('frequency, Hz');


resultFile = 'Freq_Graph.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .tiff file
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  now do analysis for each roi size

t = (1:N);
for r = r1:r2
  y = weisskoff(:, r)';
  yfit = polyval(polyfit(t, y, 2), t);  % 2nd order trend
  F(r) = std(y - yfit)/mean(yfit);
end
rr = (r1:r2);
F = 100*F;              % percent
fcalc = F(1)./rr;
rdc = F(1)/F(r2);	% decorrelation distance

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
loglog(rr, F, '-x', rr, fcalc, '--');
grid
xlabel('ROI full width, voxels');
ylabel('Relative std, %');
xlim([r1 r2]);
legend('measured', 'theoretical')

title(sprintf('rdc = %3.1f voxels // ROI size = %d voxels',rdc,R));
resultFile = 'Weisskoff_Graph.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .fig file
close;


[mean_ghost, top_ghost] = ghosting_metrics(data, NPIX, TR, path);

resultFile = 'QA_metrics.json';
file = fullfile(path, resultFile);

%data2json=struct('mean',meanI,'sd',sd,'snr',snr,'sfnr',sfnrI,'rms',100*sd/m,'temp_drift',100*drift,...
 %             'max_temp_drift',100*maxabsdrift,'rdc',rdc,'mean_ghost',mean_ghost,'top_ghost',top_ghost,...
  %            'spatial_drift_relx', spatial_drift.max_sdx, 'spatial_drift_rely', spatial_drift.max_sdy,...
   %           'spatial_drift_relz', spatial_drift.max_sdz,'image_frequency', imageFreq,'transmit_gain',...
    %          gains(1),'receiver_gain', gains(2));


qc_metrics=struct('mean',meanI,'sd',sd,'snr',snr,'sfnr',sfnrI,'rms',100*sd/m,'temp_drift',100*drift, 'temp_drift_per_minute', 100*drift_per_minute,...
          'max_temp_drift',100*maxabsdrift,'rdc',rdc, 'FWHM_x', fwhm(1), 'FWHM_y', fwhm(2), 'FWHM_z', fwhm(3), 'mean_ghost',mean_ghost,'top_ghost',top_ghost,...
          'spatial_drift_pe', spatial_drift.max_sdy,'image_frequency', imageFreq,'transmit_gain',...
          gains(1),'receiver_gain', gains(2));
series_info=struct('Multiband', 0, 'RepetitionTime', meta.TR, 'EchoTime', meta.TE, 'FlipAngle', meta.FA, 'StudyTime', meta.s_time, 'StudyDate', meta.s_date, 'StudyInstanceUID', meta.si_UID,...,
          'Manufacturer', meta.manufact, 'ManufacturerModelName', meta.model, 'Coil', meta.coil, 'SeriesDescription', meta.sDes, 'SeriesNumber', meta.se_number, 'SeriesTime', meta.se_time,...,
          'ScannerSerialNumber', meta.serialNumber, 'OperatingSystemVersion', meta.softVersion, 'OSLevel', meta.osLevel, 'ImageComments', meta.imComments); 
      
data2json=struct('QA_metrics', qc_metrics, 'SeriesInfo', series_info, 'version', version);       
          
jsonfile = savejson('fBIRN_Phantom_QA',data2json,struct('FloatFormat','%.3f'));

fileID = fopen(file,'w');
fprintf(fileID, jsonfile);
fclose(fileID);

end


function [spatial_drift] = motion_estimates(vol4D,meta,path)

TR=meta.TR/1000;
step = ceil(30/TR);
nFrames = size(vol4D,4);

if ~mod(nFrames,step)
    timepoints = 1 + fix(nFrames/step);
    lastpoint = 0;
else
    timepoints = 2 + fix(nFrames/step);
    lastpoint = 1;
end

[optimizer,metric] = imregconfig('monomodal');
%optimizer.MinimumStepLength = 1e-4;
%optimizer.MaximumIterations = 50;

fixedImage = vol4D(:,:,:,1);
Rfixed  = imref3d(size(fixedImage),meta.sx,meta.sy,meta.sz);
Rmoving = imref3d(size(fixedImage),meta.sx,meta.sy,meta.sz);

motionMatrix=zeros(4,4,timepoints);
spatial_drift_relx = zeros(timepoints,1);
spatial_drift_rely = zeros(timepoints,1);
spatial_drift_relz = zeros(timepoints,1);
for i=step:step:nFrames
   movingImage = vol4D(:,:,:,i);
   geomtform = imregtform(movingImage,Rmoving, fixedImage,Rfixed, 'rigid', optimizer, metric);
   motionMatrix(:,:,(i/step)+1) = geomtform.T;
   spatial_drift_relx((i/step)+1) = geomtform.T(4,1);
   spatial_drift_rely((i/step)+1) = geomtform.T(4,2);
   spatial_drift_relz((i/step)+1) = geomtform.T(4,3);
end

if (lastpoint)
   movingImage = vol4D(:,:,:,end);
   geomtform = imregtform(movingImage,Rmoving, fixedImage,Rfixed, 'rigid', optimizer, metric);
   motionMatrix(:,:,end) = geomtform.T;
   spatial_drift_relx(end) = geomtform.T(4,1);
   spatial_drift_rely(end) = geomtform.T(4,2);
   spatial_drift_relz(end) = geomtform.T(4,3);
end

spatial_drift = struct ('max_sdx', max(abs(spatial_drift_relx)), 'max_sdy', max(abs(spatial_drift_rely)), 'max_sdz', max(abs(spatial_drift_relz)));


x=(0:timepoints-1);
x=x*step*(TR/60);
x(end)=nFrames*TR/60;

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously

%plot(x,spatial_drift_relx,'-bo',x,spatial_drift_rely,'-ro',x,spatial_drift_relz,'-go');
plot(x,spatial_drift_rely,'-ro');
xlabel('Time (minutes)');
ylabel('Spatial Drift (mm)');
%legend('relx(LR)','rely(PA)', 'relz(IS)','Location','northwest')
legend('rely(PA)','Location','northwest')
grid

%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''off'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
resultFile = 'SpatialDrift_Graph.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .tif file
close;

end

function [ghostperc, topghostperc] = ghosting_metrics(data, NPIX, TR, path)

datav= mean(data,3);
maskshift=zeros(size(datav));

nhood = strel('disk',1,0);
phantom_edges = edge(datav,'canny');
phantom_edges = imdilate(phantom_edges,nhood);
phantom_mask = imfill(phantom_edges, 'holes');
phantom_mask_dil = imdilate(phantom_mask,nhood);

maskshift = [phantom_mask(fix(NPIX/2)+1:end,:); phantom_mask(1:fix(NPIX/2),:)];

ghostarea = maskshift & (~phantom_mask_dil);
nhood = strel('disk',3,0);
phantom_mask_erod = imerode(phantom_mask,nhood);

ghostsig = datav(find(ghostarea));
meanghost = mean(ghostsig);
meansig =   mean(datav(find(phantom_mask_erod)));
ghostperc = meanghost/meansig*100;

[topghost,s0] = (sort(ghostsig(:)));
topghost = flipud(topghost);
s0 = flipud(s0);
nvox = length(topghost);
topghost10 = topghost(1:round(nvox/10));
topghostperc =  mean(topghost10)/meansig*100;

figure; % plot ghost as a function of time.
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
sz = size(data);
usedat = reshape(data,sz(1)*sz(2),sz(3));
usedat = usedat(find(ghostarea),:);
mghost = squeeze(mean(usedat,1));
mghost = 100*mghost/meansig;

topghost_t = squeeze(mean(usedat(s0(1:round(nvox/10)),:),1));
topghost_t = 100*topghost_t/meansig;

tax = 1:length(mghost);
tax = tax*TR/60;

subplot(2,1,1);

plot(tax,mghost);grid;
title(sprintf('mean ghost: %4.2f // 1st: %4.2f // last %4.2f',...
ghostperc,mghost(1),mghost(end)));

xlabel('Time (minutes)');
ylabel('Ghost Percentage');

subplot(2,1,2);
plot(tax,topghost_t);grid;
title(sprintf('mean top 10%% ghost: %4.2f // 1st: %4.2f // last %4.2f',...
topghostperc,topghost_t(1),topghost_t(end)));

xlabel('Time (minutes)');
ylabel('Ghost Percentage');

%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''off'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
resultFile = 'Ghost_Graph.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .tiff file
close;

%%%%% DISPLAY MASKS

b_phantom=bwperim(phantom_mask_erod);
b_ghost=bwperim(ghostarea);

green=zeros(size(ghostarea,1),size(ghostarea,2),3);
green(:,:,2)=0.6;
red=zeros(size(green));
red(:,:,1)=0.8;

figure;
axes('Visible','off');
set(gcf,'Visible','off'); % this disables the figure and set the 'CreateFcn' property simultaneously
%set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
imagesc(datav);
colormap(gray);
hold all;
%h=imshow(green);
h=image(green);
set(h,'AlphaData',b_phantom);
%h=imshow(red);
h=image(red);
set(h,'AlphaData',b_ghost)

%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''off'')'); % this disables the figure and set the 'CreateFcn' property simultaneously
resultFile = 'Ghost_Region.tif';
file = fullfile(path, resultFile);
saveas(gcf, file) % save the figure in the current folder as a .tiff file
close;

end
