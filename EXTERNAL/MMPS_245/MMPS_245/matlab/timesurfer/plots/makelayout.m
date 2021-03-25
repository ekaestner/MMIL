% % collect indices of significant channels
% clear sig
% n=numel(stat_data)-1; 
% sig.pos=zeros(1,n); 
% sig.neg=zeros(1,n);
% for i=1:n,
%     try if stat_data{i}.posclusters(1).prob<=0.05, sig.pos(i)=1; end; end
%     try if stat_data{i}.negclusters(1).prob<=0.05, sig.neg(i)=1; end; end
% end
% 
% idx=find(sig.pos); for i=1:numel(idx), pvals(i)=stat_data{idx(i)}.posclusters(1).prob; end;
    
%% MAKE LAYOUT
% make a grid layout w/ columns: chan# x-pos y-pos subplot_width subplot_height channel_label

filename='/home/jsherfey/matlab/ts_iEEG_stream/development/output4/matfiles/TFstats/128layout.lay';

pad=0; cha=40; chb=119;
ncols=10; wlim=[-20 20]; 
nrows=8;  hlim=[-20 20];  

% for subplot size
wt=wlim(2)-wlim(1); pwt=(wt-pad)/ncols;
ht=hlim(2)-hlim(1); pht=(ht-pad)/nrows; 

% for plot positioning
gwt=(wt+2*pad)/ncols; lft=wlim(1)-pad+gwt/2;
ght=(ht+2*pad)/ncols; top=hlim(2)+pad-ght/2; 

% column 1: channum
channum=[1:nrows*ncols];

% columns 2 (xpos) and 3 (ypos)
cols=0:ncols-1; gridcols=lft+cols*gwt;
rows=0:nrows-1; gridrows=top-rows*ght;

xpos=gridcols; xpos=repmat(xpos,1,nrows);
ypos=gridrows'; ypos=repmat(ypos,1,ncols)'; ypos=ypos(:)';

% columns 4 (pwt) and 5 (pht)
pwt=pwt*ones(1,nrows*ncols);
pht=pht*ones(1,nrows*ncols);

% column 6 (col6.lbl)
for i=1:nrows*ncols,col6(i).lbl=stat_data{end}.label{cha+i-1}; end

out.col1=channum';
out.col2=xpos';
out.col3=ypos';
out.col4=pwt';
out.col5=pht';
out.col6=col6;

fid=fopen(filename,'w+');
for i=1:nrows*ncols,
    fprintf(fid,'%d %f %f %f %f %s\n',out.col1(i),out.col2(i),out.col3(i),out.col4(i),out.col5(i),out.col6(i).lbl);
end
fclose(fid);

%% Plot stats with multiplotTFR
cfg=[];
clear cond1vs2
samples=161:281;
cond1vs2=stat_data{end}; 
cond1vs2.label=cond1vs2.label(cha:chb);
cond1vs2.powspctrm=cond1vs2.powspctrm(cha:chb,:,samples);
cond1vs2.time=cond1vs2.time(samples);
cond1vs2.dof=cond1vs2.dof(cha:chb,:,samples);
cond1vs2.cfg.channel=cond1vs2.cfg.channel(cha:chb);

% calculate z-scores
avgs=nanmean(cond1vs2.powspctrm(:,:,:),3);
avgs=repmat(avgs,[1 1 size(cond1vs2.powspctrm,3)]);
stds=nanstd(cond1vs2.powspctrm(:,:,:),[],3);
stds=repmat(stds,[1 1 size(cond1vs2.powspctrm,3)]);
power=(cond1vs2.powspctrm-avgs)./stds;

mask=[];
for i=cha:chb
    mask=cat(1,mask,stat_data{i}.mask);
end

cond1vs2.mask=mask;
cfg.maskparameter='mask';

cfg.xparam='time';
cfg.yparam='freq';
cfg.zparam='powspctrm';
cfg.showlabels='yes';
cfg.xlim=[.3 .6];
cfg.ylim=[1 50];
cfg.zlim=[-50 50];
cfg.rotate=0;
cfg.layout=filename;
cfg.colorbar='yes';
figure; multiplotTFR(cfg,cond1vs2);

% 
% cfg.zlim=[-.5 .5];
% figure; multiplotTFR(cfg,cond1vs2);
% cfg.zlim=[-5 5];
% figure; multiplotTFR(cfg,cond1vs2);
% cfg.zlim=[-50 50];
% figure; multiplotTFR(cfg,cond1vs2);



