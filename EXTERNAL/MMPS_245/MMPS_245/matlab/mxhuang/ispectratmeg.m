%function ispectratmeg(datafile, event2file, freqVec, savefile, factor, Fs, xmin, xmax, wfactor, ChnVec,ChannelNumber);
% modified ispectrat.m to adapt MEG data
% Chunmao Wang 10/4/2001
% 
% Induced spectrat - also output SE
% ispectratmeg(datafile, event2file, freqVec, savefile, factor, Fs, xmin, xmax, wfactor);
%
% Chunmao Wang, 12/01/2000

if size(savefile,1)~=size(factor,1),
   help ispectra
   return
end
spline=3;

%ChannelNumber=379;	% all channels
PointNumber=ceil(Fs*(xmax-xmin)/1000);

[fid message] = fopen(datafile,'r','ieee-be');
E = load(event2file);

EventNumber = size(E,1);
ChnNumber = size(ChnVec,2);

PointNumber2=floor(PointNumber/spline);

for f=1:size(factor,1),
   TFR=zeros(ChnNumber,size(freqVec,2),PointNumber2);
   TFRse=zeros(ChnNumber,size(freqVec,2),PointNumber2);
   for i=1:ChnNumber,
      k=0;
      %tf=zeros(EventNumber,size(freqVec,2),PointNumber2);
      tfr=zeros(size(freqVec,2),PointNumber2);
      for j=1:EventNumber
         if (findstr(E(j,1),factor(f,:))>0)
%            fprintf(1,'%d',j);
            k=k+1;
            STATUS=fseek(fid,(ChannelNumber+2)*2*(E(j,2)-ceil(Fs*(0-xmin)/1000)-1)+(ChnVec(i)+1)*2,'bof');
            if STATUS==-1,
               fprintf(1,'\n%s\n','fseek error!');
               pause;
            end
            for n=1:PointNumber
               [d(n) count]=fread(fid,1,'short');
               STATUS=fseek(fid,(ChannelNumber+2)*2-2,'cof');
            end
            %figure;plot(d);
            
            % spline fit d
            d=d(1:spline*floor(PointNumber/spline));
            d=reshape(d,spline,floor(PointNumber/spline));
            d=squeeze(mean(d));
%            for ni=1:PointNumber2,
%               fprintf(1,'ch%d event %d pnt%d d %f\n', i,j,ni,d(ni));
%            end
%            pause;
            
            Fs2=Fs/spline;
            
            w=0;
            d=d';
            for fq=1:size(freqVec,2)
               [tf(k,fq,:),timeVec,fq] = traces2TFR(d,freqVec(fq),Fs2,wfactor(fq));
            end
            clear d;
            tfr=tfr+squeeze(tf(k,:,:));
%            fprintf(1,'\b\b\b\b');
         end
      end
      TFR(i,:,:)=tfr/k;
%      for fq=1:size(freqVec,2),
%         for ni=1:PointNumber2,
%            fprintf(1,'ch%d fq%f pnt%d tfr %f\n', i,freqVec(fq),ni,TFR(i,fq,ni));
%         end
%         pause;
%      end

      for j=1:size(freqVec,2)
         TFRse(i,j,:)=squeeze(std(tf(:,j,:),0,1))/(k^0.5);
%         for en=1:k,
%            fprintf(1,'%f\n', tf(en,1,1));
%         end
%         fprintf(1,'std=%f\n', TFRse(1,j,1));
%         pause;
      end
%      for fq=1:size(freqVec,2),
%         for ni=1:PointNumber2,
%            fprintf(1,'ch%d fq%f pnt%d tfr %f tfrse %f\n', i,freqVec(fq),ni,TFR(i,fq,ni),TFRse(i,fq,ni));
%         end
%         pause;
%      end
      fprintf(1,'%d',i);
   end
   fprintf(1,'\n');
   clear tfr;
   clear tf;
   
   fprintf('\n');  
   sTFR=reshape(TFR,ChnNumber,size(freqVec,2)*PointNumber2)';
   clear TFR;
   eval (['save '  savefile(f,:) '.itf.mat k sTFR'])
   clear sTFR;
   sTFRse=reshape(TFRse,ChnNumber,size(freqVec,2)*PointNumber2)';
   clear TFRse;
   eval (['save '  savefile(f,:) 'se.itf.mat k sTFRse'])
   clear sTFRse;
end
fclose(fid);