function syncview
% Idea:
% 1. Elements are arranged in space according to meaningful layouts
% 2. Element selection activates corresponding object
% Specifics:
% Element = subplot in hControl
% Element ID (tag) = sensor label
% Object = image matrix in file on disk: files{ID}
% Action: subplot ButtonDownFcn displays image in hDisplay

% images(page#).category(condition#).matrix{RefID} => 2D real
clear global
global files images

% Draw UI Controls
figure('name','Select files','numbertitle','off');
geometry = {[1 2] [1 2] [1 2]};
geomvert = [1 1 1];
listui   = { ...
           {'style','pushbutton','string','layouts' ,'callback',{@getfiles,'lay'}} ...  
           {'style','edit'      ,'string','layout file names','tag','lay'} ...
           {'style','pushbutton','string','images'  ,'callback',{@getfiles,'jpg'}} ...           
           {'style','edit'      ,'string','image file names' ,'tag','jpg'} ...
           {'style','pushbutton','string','load'    ,'callback',{@setup}} ...
           {'style','text'      ,'string',''} ...           
           };
supergui(gcf,geometry,geomvert,listui{:});
fac     = get(0,'screensize');
pos     = get(gcf,'position');
pos(3)  = 300;
pos(4)  = 100;
pos(1)  = (fac(3)-pos(3))/2;
pos(2)  = fac(4) - 1.5*pos(4);%(fac(4)-pos(4));
set(gcf,'position',pos);

% CALLBACKS
function getfiles(src,evnt,format)
global files
% layout (.lay): ID | xpos | ypos | width | height | label
% images (.jpg): jpegs
[f,p]=uigetfile(['*.' format],'multiselect','on');
if isequal(f,0) || isequal(p,0),return; end; 
if iscell(f)
  for n = 1:length(f)
    files.(format){n} = fullfile(p,f{n});
  end; 
else
  files.(format){1}=fullfile(p,f);
end;
set(findobj('tag',format),'string',files.(format),'max',length(files.(format)));

function setup(src,evnt)
global files images % files.lay,files.jpg
close(gcf)
Nlay = length(files.lay);
Njpg = length(files.jpg);
files.CurrID = 1;

% combine layouts
for n = 1:Nlay
  cfg.layout = files.lay{n};
  this  = prepare_layout(cfg); % this layout:pos [x y],width,height,label
  % shift so (0,0) is lower left if any positions are negative
  if min(this.pos(:,1))<0
    this.pos(:,1) = this.pos(:,1) - min(this.pos(:,1));
  end
  if min(this.pos(:,2))<0
    this.pos(:,2) = this.pos(:,2) - min(this.pos(:,2));
  end  
  if n == 1
    lay = rmfield(this,'pos');
    lay.xpos = this.pos(:,1);
    lay.ypos = this.pos(:,2);
  else
    ind = length(lay.label) + [1:length(this.label)];
    lay.label(ind)  = this.label;
    lay.width(ind)  = this.width;
    lay.height(ind) = this.height;
    lay.xpos(ind)   = this.pos(:,1);
    lay.ypos(ind)   = this.pos(:,2) + max(lay.ypos);
  end
end

% 
% select files with refchans in lay
Nelm   = length(lay.label);
keep   = [];
sel    = [];
cnt    = 0;
for n  = 1:Nelm
  this = lay.label{n};
  % find indices to files containing this label
  tmp  = find(~cellfun(@isempty,strfind(files.jpg,['-' this '_'])));
  if length(tmp)==0, continue; end
  keep = [keep tmp];
  sel  = [sel n];
  cnt  = cnt + 1;
  jpg2lay{cnt} = tmp;
end
if isempty(keep)
  warning('no file names contained the sensor labels in the layout file');
  return;
else
  keep = sort(keep);
end

x = lay.xpos(sel);
y = lay.ypos(sel);
w = lay.width(sel);  a = max(x)+1.5*mean(w(x==max(x)));
h = lay.height(sel); b = max(y)+1.5*mean(h(y==max(y)));
w = w / a;
h = h / b;
x = x / a;
y = y / b;
labels = lay.label(sel);
Nelm   = length(labels);

% get indices for grouping files by category & page
fprintf('Loading %g images...',Njpg); tic
res = regexp(files.jpg,'event(\w*)(\d)+','match');
res = res(~cellfun(@isempty,res));
res = cellfun(@(x)(x{1}),res,'UniformOutput',false);
sak = [];
images     = [];
Categories = unique({res{:}});
Ncategory  = length(Categories) + 1;
for n = 1:Ncategory
  if n < Ncategory
    CatString = Categories{n};
    % get indices to files for this category
    ind   = regexp(files.jpg,CatString);
    ind   = find(~cellfun(@isempty,ind));
  else
    CatString = 'MiscCategory';
    % all remaining indices
    ind   = setdiff(1:Njpg,sak);
  end
  % select indices to files with reference channels in lay
  ind = ind(ismember(ind,keep));
  if isempty(ind), continue; end % n==Ncategory && ...
  images(n).category = CatString;
  images(n).index    = ind;
  % track indices
  sak = [sak ind];
  % divide this category by page
  res = regexp(files.jpg(ind),'page(\d)+','match'); 
  res = res(~cellfun(@isempty,res));
  res = cellfun(@(x)(x{1}),res,'UniformOutput',false); 
  Pages  = unique({res{:}});                        % page IDs
  Npage  = length(Pages) + 1;                         % need one group per 'page' plus others
  catsak = [];
  for m = 1:Npage
    if m < Npage
      PageString = Pages{m};
      % get indices to files for this page in this category
      subind = regexp(files.jpg(ind),PageString);
      subind = find(~cellfun(@isempty,subind));
    else
      PageString = 'MiscPage';
      % all remaining indices in this category
      subind = setdiff(1:length(ind),catsak);
    end
    catsak   = [catsak subind];
    if m==Npage && isempty(subind), continue; end
    % record index info
    index = ind(subind);
    images(n).pages(m).page     = PageString;
    images(n).pages(m).index    = index;
    images(n).pages(m).relindex = subind;
    images(n).pages(m).FigureTag= PageString;
    images(n).pages(m).image    = {};
    images(n).pages(m).label    = {};
    
    % loop over jpegs
    for k = 1:length(index)
      % get lay index for each jpeg
      refindex = find(cellfun(@(x)(ismember(index(k),x)),jpg2lay));
      if isempty(refindex), continue; end
      % load image
      images(n).pages(m).image{refindex} = imread(files.jpg{index(k)});
      images(n).pages(m).label{refindex} = labels{refindex};
    end
  end
end
fprintf('done.\n'); toc

% draw hControl
hControl = figure('color','w','tag','hControl','Name','Reference Channel Controller');
fac = get(0,'screensize');
pos = get(hControl,'position');
pos = [.05 .05 .3*fac(3) min(Nlay*.3*fac(3),fac(4))];%.9*fac(4)];
set(hControl,'position',pos);
% draw control subplots
for n = 1:Nelm
  ax(n) = subplot('Position',[x(n),y(n),w(n),h(n)]);
  axis on; box on; title(labels{n},'verticalalignment','middle');
  set(ax(n),'color','w','tag','elm','ButtonDownFcn',{@action,n},...
    'XTickLabel','','YTickLabel','');
end
files.Nelm = Nelm;
files.hax  = ax;
% add listbox to hControl
uicontrol('Parent',hControl,'style','popupmenu','tag','CategoryList','value',1,'string',{images.category},'Callback',{@action,0},'Units','normalized','Position',[.3 0 .15 .05]);
uicontrol('Parent',hControl,'style','togglebutton','tag','ToggleRun','string','Run','Callback',@run,'Units','normalized','Position',[.55 0 .15 .05]);

% initialize
InitID   = files.CurrID; 
action(ax(InitID),[],InitID);

function action(src,evnt,ID)
global images files
if ID==0
  ID = files.CurrID; 
else
  files.CurrID = ID; 
end
% get category from uicontrol
imageset = images(get(findobj('tag','CategoryList'),'value'));
fac = get(0,'screensize');
step= .1;
w   = .6*fac(3);
N   = length(imageset.pages);
xoff= step*(N-1)*w;
fac = [.35*fac(3) .05 w-xoff .9*fac(4)];%[1 1]*min(.6*fac(3),.95*fac(4))];
% loop over pages
for p  = 1:length(imageset.pages)
  hFig = findobj('tag',imageset.pages(p).FigureTag);
  if isempty(hFig)
    hDisplay = figure('color','w','tag',imageset.pages(p).FigureTag); 
    pos = fac + [step*(p-1)*w 0 0 0];
    set(hDisplay,'position',pos,'Name',sprintf('%s - %s',imageset.category,imageset.pages(p).page));
  else
    figure(hFig);
  end
  if ID <= numel(imageset.pages(p).image)
    image(imageset.pages(p).image{ID});
    title(imageset.pages(p).label{ID}); axis off;
  end
end
set(findobj('tag','elm'),'color','w');
if strcmp(get(src,'tag'),'elm')
  set(src,'color','k');
end
% set(gcf,'name',files.jpg{ID});

function run(src,evnt)
global files
N  = files.Nelm;
ID = files.CurrID;
h  = files.hax;
while get(src,'value')  
  action(h(ID),[],ID); pause(1);
  if ID == N, ID = 1; else ID = ID + 1; end
end
