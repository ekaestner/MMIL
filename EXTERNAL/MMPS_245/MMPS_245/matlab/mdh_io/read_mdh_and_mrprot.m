function [kdat,krev,ktim,pcdat,pcrev,pctim,h] = read_mdh_and_mrprot(fdir,chan_id,replist,echolist,klist,idim,jdim,npclines,flag3d)

VersionNum = 15;
fname = sprintf('%s/mrprot.asc',fdir);
fid = fopen(fname, 'r', 'l');
if fid < 0
  fname = sprintf('%s/MrProt.asc',fdir);
  fid = fopen(fname, 'r', 'l');
  if fid < 0
    fname = sprintf('%s/meas.asc',fdir);
    fid = fopen(fname, 'r', 'l');
    if (fid >= 0)
      VersionNum = 21;
    else
      disp(['can''t open file mrprot.asc, MrProt.asc or meas.asc']);
      return;
    end
  end
end

tic

if 1
  h.SequenceFileName     = get_asch_value(fid,'tSequenceFileName');
  h.FlipAngleDegree      = get_asch_value(fid,'dFlipAngleDegree');
  h.TxFrequency          = get_asch_value(fid,'sTXSPEC.lFrequency');
  h.BaseResolution       = get_asch_value(fid,'sKSpace.lBaseResolution');
  h.PhaseEncodingLines   = get_asch_value(fid,'sKSpace.lPhaseEncodingLines');
  h.PhaseResolution      = get_asch_value(fid,'sKSpace.dPhaseResolution');
  h.Partitions           = get_asch_value(fid,'sKSpace.lPartitions');
  h.TR                   = get_asch_value(fid,'alTR[0]');
  h.TI                   = get_asch_value(fid,'alTI[0]');
  echo_num = 1;
  f                      = get_asch_value(fid,sprintf('alTE[%d]',echo_num-1));
  h.TE(echo_num)         = f;
  while (f~=0)
    echo_num             = echo_num+1;
    f                    = get_asch_value(fid,sprintf('alTE[%d]',echo_num-1));
    h.TE(echo_num)       = f;
  end
  h.SliceArraySize       = get_asch_value(fid,'sSliceArray.lSize');
  for slice_num=[1,h.SliceArraySize] % read info for first and last slices only
    %for slice_num=[1:h.SliceArraySize] % read info for all slices
    h.SlicePositionSag(slice_num) = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sPosition.dSag',slice_num-1));
    h.SlicePositionCor(slice_num) = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sPosition.dCor',slice_num-1));
    h.SlicePositionTra(slice_num) = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sPosition.dTra',slice_num-1));
    h.SliceNormalSag(slice_num)   = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sNormal.dSag',slice_num-1));
    h.SliceNormalCor(slice_num)   = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sNormal.dCor',slice_num-1));
    h.SliceNormalTra(slice_num)   = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].sNormal.dTra',slice_num-1));
    h.Thickness(slice_num)        = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].dThickness',slice_num-1));
    h.PhaseFOV(slice_num)         = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].dPhaseFOV',slice_num-1));
    h.ReadoutFOV(slice_num)       = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].dReadoutFOV',slice_num-1));
    h.InPlaneRot(slice_num)       = get_asch_value(fid,sprintf('sSliceArray.asSlice[%d].dInPlaneRot',slice_num-1));
  end

  h.Center = [mean(h.SlicePositionSag([1,h.SliceArraySize])); ...
    mean(h.SlicePositionCor([1,h.SliceArraySize])); ...
    mean(h.SlicePositionTra([1,h.SliceArraySize]))];

end
fclose(fid);

%keyboard

% NEED TO CONVERT HEADER INFO TO vox2ras matrix!


fname = sprintf('%s/meas.out',fdir);
fid = fopen(fname, 'r', 'l');
if fid < 0
  disp(['can''t open file ' fname]);
  return;
end

kdim = length(klist);
repdim = length(replist);
echodim = length(echolist);

kdat = complex(zeros(kdim,idim,jdim,repdim,echodim,'single'),zeros(kdim,idim,jdim,repdim,echodim,'single'));
krev = zeros(kdim,idim,repdim,echodim,'single');
ktim = zeros(kdim,idim,repdim,echodim,'single');
echocntr = zeros(kdim,idim,repdim,'single');
if (npclines>0)
  pcdat = complex(zeros(kdim,npclines,jdim,repdim,echodim,'single'),zeros(kdim,npclines,jdim,repdim,echodim,'single'));
  pcrev = zeros(kdim,npclines,repdim,echodim,'single');
  pctim = zeros(kdim,npclines,repdim,echodim,'single');
  pcnum = zeros(kdim,repdim,echodim,'single');
else
  pcdat = [];
  pcrev = [];
  pctim = [];
  pcnum = [];
end

Repetition = 0;
meas_out_start_offset = fread(fid, 1, 'long');
fseek(fid, meas_out_start_offset, 'bof');
[adc_data, mdh] = read_mdh_adc(fid,VersionNum);
firstTS = mdh.TimeStamp;
while (~feof(fid) & ~isempty(adc_data) & (Repetition<=max(replist)))
  i = mdh.LoopCounter.Line+1;
  if (flag3d)
    k = mdh.LoopCounter.Partition+1;
  else
    k = mdh.LoopCounter.Slice+1;
  end
  Repetition = mdh.LoopCounter.Repetition;
  Echo = mdh.LoopCounter.Echo;
  online = mdh.EvalInfoMask(4);
  phasecor = mdh.EvalInfoMask(22);
  reflect = mdh.EvalInfoMask(25);
  repindx = find(Repetition==replist);
  kindx = find(k==klist);

  %% hack to work around incorrect Echo numbering in multi-echo EPI seq.
  if (~phasecor & ~isempty(kindx) & ~isempty(repindx) & (chan_id == mdh.ChannelId))
    Echo = echocntr(kindx,i,repindx);
    echocntr(kindx,i,repindx) = echocntr(kindx,i,repindx)+1;
  end

  echoindx = find(Echo==echolist);
  if (online & ~isempty(repindx) & ~isempty(kindx) & (~isempty(echoindx) | phasecor) & (chan_id == mdh.ChannelId))
    if (reflect)
      adc_data = adc_data(length(adc_data):-1:1) ;
    end
    if (phasecor)
      echoindx = 1; % Make all phase correction navigators belong to echo 1
      pcnum(kindx,repindx,echoindx) = pcnum(kindx,repindx,echoindx)+1;
      pcdat(kindx,pcnum(kindx,repindx,echoindx),:,repindx,echoindx) = adc_data;
      pcrev(kindx,pcnum(kindx,repindx,echoindx),repindx,echoindx) = reflect;
      pctim(kindx,pcnum(kindx,repindx,echoindx),repindx,echoindx) = mdh.TimeStamp;
      %      fprintf(1,'### phasecor=%d, Chan=%d, Acq=%d, Rep=%d, Echo=%d (%d, %d), k=%d, i=%d\n', ...
      %              phasecor,mdh.ChannelId,mdh.LoopCounter.Acquisition,mdh.LoopCounter.Repetition,mdh.LoopCounter.Echo,Echo,echoindx,k,i);
    else
        kdat(kindx,i,:,repindx,echoindx) = single(adc_data);
        krev(kindx,i,repindx,echoindx) = reflect;
        ktim(kindx,i,repindx,echoindx) = mdh.TimeStamp;
    end
    if (toc>1) %(k==1) % Echo==1 %(i==1)
      fprintf(1,'*** phasecor=%d, Chan=%d, Acq=%d, Rep=%d, Echo=%d (%d), k=%d, i=%d\n', ...
        phasecor,mdh.ChannelId,mdh.LoopCounter.Acquisition,mdh.LoopCounter.Repetition,mdh.LoopCounter.Echo,Echo,k,i);
      tic;
    end
  else
    if (toc>1) %(k==1) % Echo==1 %(i==1)
      fprintf(1,'phasecor=%d, Chan=%d, Acq=%d, Rep=%d, Echo=%d (%d), k=%d, i=%d\n', ...
        phasecor,mdh.ChannelId,mdh.LoopCounter.Acquisition,mdh.LoopCounter.Repetition,mdh.LoopCounter.Echo,Echo,k,i);
      tic;
    end
  end
  [adc_data, mdh] = read_mdh_adc(fid,VersionNum);
end
fclose(fid);
