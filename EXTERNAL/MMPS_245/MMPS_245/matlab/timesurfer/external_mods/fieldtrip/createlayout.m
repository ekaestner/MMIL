function lay = createlayout(input);

% CREATELAYOUT will make a nice 2D layout of your channel positions
% that is suitable for topoplotting or multiplotting
% 
% Use as
%   lay = createlayout(input)
% where the input can be any of the following
%   string, indicating an ascii layout file (*.lay)
%   string, indicating a 3D electrode positions file
%   elec-structure, containing 3D electrode positions
%   grad-structure, containing 3D magnetometer/gradiometer positions
%
% See also READ_LAY, WRITE_LAY, READ_FCDC_ELEC

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: createlayout.m,v $
% Revision 1.7  2005/06/15 08:21:47  roboos
% added boxes to debug plotting of layout
%
% Revision 1.6  2005/05/17 17:50:49  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.5  2005/04/05 19:33:27  roboos
% made the addition of the COMNT and SCALE location optional, i.e. only if not already present
%
% Revision 1.4  2005/04/05 16:02:41  roboos
% added placeholder for COMNT and SCALE location
%
% Revision 1.3  2005/02/08 13:45:40  roboos
% changed from %s to %q for reading channel names with a space in them
%
% Revision 1.2  2005/02/07 17:37:16  roboos
% swapped order of channels in the triplet according to Dasha's specification
%
% Revision 1.1  2005/02/07 17:15:27  roboos
% new implementations, are used in multiplotting and topoplotting
%

% the different options for the input are handled below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isstr(input) && filetype(input, 'layout')
  fprintf('reading layout from file %s\n', input);
  if ~exist(input, 'file')
    error(sprintf('could not open layout file: %s', input));
  end
  % read the layout information from the ascii file
  [chNum,X,Y,Width,Height,Lbl] = textread(input,'%f %f %f %f %f %q');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isstr(input) && ~filetype(input, 'layout')
  % assume that it is an electrode file
  fprintf('creating layout from electrode file %s\n', input);
  if ~exist(input, 'file')
    error(sprintf('could not open electrode file: %s', input));
  end
  elec = read_fcdc_elec(input);
  % convert electrode positions into layout
  prj = elproj(elec.pnt, 'stereographic'); % * [0 1; -1 0]
  d = dist(prj');
  d(find(eye(size(d)))) = inf;
  mindist = min(d(:));
  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;
  Lbl = elec(1).label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isstruct(input) && isfield(input, 'pnt') && ~isfield(input, 'ori')
  % rename the input for convenience
  elec = input;
  % convert electrode positions into layout
  Lbl = elec.label;
  prj = elproj(elec.pnt, 'stereographic'); % * [0 1; -1 0];
  d = dist(prj');
  d(find(eye(size(d)))) = inf;
  mindist = min(d(:));
  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isstruct(input) && isfield(input, 'pnt') && isfield(input, 'ori')
  % rename the input for convenience
  grad = input;
  % convert magnetometer/gradiometer positions into layout
  fprintf('creating layout for %s MEG system\n', megsystem(grad));

  switch lower(megsystem(grad))

  case {'ctf151', 'ctf275'}
    Lbl = channelselection('MEG', grad.label);
    ind = match_str(grad.label, Lbl);
    prj = elproj(grad.pnt(ind,:), 'stereographic') * [0 1; -1 0];
    d = dist(prj');
    d(find(eye(size(d)))) = inf;
    mindist = min(d(:));
    X = prj(:,1);
    Y = prj(:,2);
    Width  = ones(size(X)) * mindist * 0.8;
    Height = ones(size(X)) * mindist * 0.6;

  case {'ctf151_planar', 'ctf275_planar'}
    % create a list with planar channel names
    chan = {};
    for i=1:length(grad.label)
      if findstr(grad.label{i}, '_dH') || findstr(grad.label{i}, '_dV')
        chan{i} = grad.label{i}(1:(end-3));
      end
    end
    chan = unique(chan);
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:length(chan)
      ch1 =  [chan{i} '_dH'];
      ch2 =  [chan{i} '_dV'];
      sel = match_str(grad.label, {ch1, ch2});
      if length(sel)==2
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(grad.pnt(find(grad.tra(sel(1),:)),:), 1);
        meanpnt2 = mean(grad.pnt(find(grad.tra(sel(2),:)),:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    prj = elproj(pnt, 'stereographic') * [0 1; -1 0];
    X = prj(:,1);
    Y = prj(:,2);
    d = dist(prj');
    d(find(eye(size(d)))) = inf;
    mindist = min(d(:));
    X1 = X; Y1 = Y + 0.21 * mindist;
    X2 = X; Y2 = Y - 0.21 * mindist;
    X = [X1; X2];
    Y = [Y1; Y2];
    Lbl = [lab(:,1); lab(:,2)];
    Width  = ones(size(X))*mindist/2;
    Height = ones(size(X))*mindist/3.0;

  case 'neuromag122'
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:2:140
      ch1 = sprintf('MEG %03d', i);
      ch2 = sprintf('MEG %03d', i+1);
      sel = match_str(grad.label, {ch1, ch2});
      if length(sel)==2
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(grad.pnt(find(grad.tra(sel(1),:)),:), 1);
        meanpnt2 = mean(grad.pnt(find(grad.tra(sel(2),:)),:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    prj = elproj(pnt, 'stereographic'); % * [0 1; -1 0];
    X = prj(:,1);
    Y = prj(:,2);
    d = dist(prj');
    d(find(eye(size(d)))) = inf;
    mindist = min(d(:));
    X1 = X; Y1 = Y + 0.21 * mindist;
    X2 = X; Y2 = Y - 0.21 * mindist;
    X = [X1; X2];
    Y = [Y1; Y2];
    Lbl = [lab(:,1); lab(:,2)];
    Width  = ones(size(X))*mindist/2;
    Height = ones(size(X))*mindist/3.0;

  case 'neuromag306'
    % find the matching channel-triplets
    ind = [];
    lab = {};
    for i=1:300
      ch1 = sprintf('MEG %03d1', i);
      ch2 = sprintf('MEG %03d2', i);
      ch3 = sprintf('MEG %03d3', i);
      sel = match_str(grad.label, {ch1, ch2, ch3});
      if length(sel)==3
        ind = [ind; i];
        lab(i,:) = {ch1, ch2, ch3};
        meanpnt1 = mean(grad.pnt(find(grad.tra(sel(1),:)),:), 1);
        meanpnt2 = mean(grad.pnt(find(grad.tra(sel(2),:)),:), 1);
        meanpnt3 = mean(grad.pnt(find(grad.tra(sel(3),:)),:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2; meanpnt3], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    prj = elproj(pnt, 'stereographic'); % * [0 1; -1 0];
    X = prj(:,1);
    Y = prj(:,2);
    d = dist(prj');
    d(find(eye(size(d)))) = inf;
    mindist = min(d(:));
    X1 = X - 0.2 * mindist; Y1 = Y + 0.2 * mindist;
    X2 = X - 0.2 * mindist; Y2 = Y - 0.2 * mindist;
    X3 = X + 0.2 * mindist; Y3 = Y + 0.0 * mindist;
    X = [X1; X2; X3];
    Y = [Y1; Y2; Y3];
    Lbl = [lab(:,1); lab(:,2); lab(:,3)];
    Width  = ones(size(X))*mindist/3;
    Height = ones(size(X))*mindist/3;

  otherwise
    error('unknown type of magnetometer/gradiometer system');
  end
end

if ~any(strcmp('COMNT', Lbl))
  % add a placeholder for the comment in the upper left corner
  Lbl{end+1}    = 'COMNT';
  Width(end+1)  = Width(1);
  Height(end+1) = Height(1);
  X(end+1)      = min(X);
  Y(end+1)      = max(Y);
end

if ~any(strcmp('COMNT', Lbl))
  % add a placeholder for the scale in the upper right corner
  Lbl{end+1}    = 'SCALE';
  Width(end+1)  = Width(1);
  Height(end+1) = Height(1);
  X(end+1)      = max(X);
  Y(end+1)      = max(Y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine the output info in a convenient structure that is similar to a elec/grad structure
lay.prj   = [X Y];
lay.box   = [Width Height];
lay.label = Lbl;

% to plot the layout for debugging, you can use this code snippet
if 0
  X      = lay.prj(:,1);
  Y      = lay.prj(:,2);
  Width  = lay.box(:,1);
  Height = lay.box(:,2);
  Lbl    = lay.label;
  figure
  plot(X, Y, '.');
  text(X, Y, Lbl);
  line([X X+Width X+Width X X]',[Y Y Y+Height Y+Height Y]');
end

