function iplsplg(datafile, event2file, savefile, factor, freqVec, Fs, xmin, xmax, wfactor, CompCh);
%
% Induced phase lock statistics (and phase lag 7/6/01)
%
% Chunmao Wang, 6/7/2001

%datafile='BI4WME1H.ASC';
%event2file='bi4wme240.ev2';
%savefile=['BI4WME1H-all-1-1-50'];
%factor=[1, 2];
%freqVec=1:1:50;
%Fs=200;
%width=12;

%xmin=-500;
%xmax=1500;
PointNumber=Fs*(xmax-xmin)/1000;
dx = (xmax-xmin)/(PointNumber-1);
x=xmin*ones(1,PointNumber)+dx*(0:PointNumber-1); % compute x-values
xmax = xmax*PointNumber/PointNumber;
FigureTitle=datafile;

data=load(datafile);
E = load(event2file);

[EventNumber ChnNumber] = size(data);
EventNumber = EventNumber/PointNumber;

data = reshape(data,PointNumber,EventNumber,ChnNumber);

for f=1:size(factor,1),
   n=0;
   PLG=[];B=[];Bstat=[];Bplf=[];
   for i=1:(size(CompCh,2)-1),
      for l=(i+1):size(CompCh,2),
         n=n+1;
         k=0;
         for j=1:EventNumber
            if (findstr(E(j,3),factor(f,:))>0)
               k=k+1;
               d1(:,k)=squeeze(data(:,j,CompCh(i)));
               d2(:,k)=squeeze(data(:,j,CompCh(l)));
            end   
         end
         w=0;
         for fq=1:(size(freqVec,2))
            [plg,b,bstat,bplf,fv,timeVec] = traces2PLSPLG(d1,d2,freqVec(fq),Fs,wfactor(fq));
            PLG(n,fq,:)=plg;
            B(n,fq,:)=abs(b);
            Bstat(n,fq,:)=bstat;
            Bplf(n,fq,:)=bplf;
         end
         fprintf(1,'%d ',i);
      end   
      fprintf('\n'); 
   end
   fprintf('\n'); 
   PLG=reshape(PLG,(size(CompCh,2)-1)*size(CompCh,2)/2,size(freqVec,2)*PointNumber);
   eval (['save '  savefile(f,:) '.plg PLG -ascii'])
   B=reshape(abs(B),(size(CompCh,2)-1)*size(CompCh,2)/2,size(freqVec,2)*PointNumber);
   eval (['save '  savefile(f,:) '.B B -ascii'])
   Bstat=reshape(Bstat,(size(CompCh,2)-1)*size(CompCh,2)/2,size(freqVec,2)*PointNumber);
   eval (['save '  savefile(f,:) '.Bstat Bstat -ascii'])
   Bplf=reshape(abs(Bplf),(size(CompCh,2)-1)*size(CompCh,2)/2,size(freqVec,2)*PointNumber);
   eval (['save '  savefile(f,:) '.Bplf Bplf -ascii'])
end