function txtHandle = TextZoomable(x,y,varargin)
% txtHandle = FixedSizeText(x,y,varargin)
% 
% Adds text to a figure in the normal manner, except that this text
% grows/shrinks with figure scaling and zooming, unlike normal text that
% stays at a fixed font size during figure operations. Note it scales with
% figure height - for best scaling use 'axis equal' before setting up the
% text.
%
% All varargin{:} arguments will be passed directly on to the text
% function (text properties, etc.)
%
% (doesn't behave well with FontUnits = 'normalized')
%
% example:
% 
% figure(1); clf;
% rectangle('Position', [0 0 1 1]);
% rectangle('Position', [.25 .25 .5 .5]);
% 
% th = TextZoomable(.5, .5, 'red', 'color', [1 0 0], 'Clipping', 'on');
% th2 = TextZoomable(.5, .1, 'blue', 'color', [0 0 1]);
%
%
% Ken Purchase, 4-25-2013
%
%

    % create the text
    txtHandle = text(x,y,varargin{:});

    % detect its size relative to the figure, and set up listeners to resize
    % it as the figure resizes, or axis limits are changed.
    hAx = gca;
    hFig = get(hAx,'Parent');
    
    fs        = get(txtHandle, 'FontSize');
    ratios    = fs * diff(get(hAx,'YLim')) / max(get(hFig,'Position') .* [0 0 0 1]);
    ax_ratios = fs * max(get(hAx,'Position') .* [0 0 0 1]);
    
    % append the handles and ratios to the user data - repeated calls will
    % add each block of text to the list
    ud = get(hAx, 'UserData');
    if isfield(ud, 'fs')
        ud.fs = [ud.fs(:) fs];
    else
        ud.fs = fs;
    end
    if isfield(ud, 'ratios')
        ud.ratios = [ud.ratios(:); ratios];
    else
        ud.ratios = ratios;
    end
    if isfield(ud, 'ax_ratios')
        ud.ax_ratios = [ud.ratios(:); ratios];
    else
        ud.ax_ratios = ax_ratios;
    end
    if isfield(ud, 'handles')
        ud.handles = [ud.handles(:); txtHandle];
    else
        ud.handles = txtHandle;
    end
    
    set(hAx,'UserData', ud);
    localSetupPositionListener(hFig,hAx);
    localSetupLimitListener(hAx);
    localSetupAxesListener(hAx);
    
end


%% Helper Functions

function localSetupPositionListener(hFig,imAxes)
    % helper function to sets up listeners for resizing, so we can detect if
    % we would need to change the fontsize
    PostPositionListener = handle.listener(hFig,'ResizeEvent',...
        {@localPostPositionListener,imAxes});
    setappdata(hFig,'KenFigResizeListeners',PostPositionListener);
end

function localPostPositionListener(~,~,imAxes) 
    % when called, rescale all fonts in image
    ud = get(imAxes,'UserData');
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
        set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end   
end

function localSetupLimitListener(imAxes)
    % helper function to sets up listeners for zooming, so we can detect if
    % we would need to change the fontsize
    hgp     = findpackage('hg');
    axesC   = findclass(hgp,'axes');
    LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
        'PropertyPostSet',@localLimitListener);
    hFig = get(imAxes,'Parent');
    setappdata(hFig,'KenAxeResizeListeners',LimListener);
end

function localLimitListener(~,event)
    % when called, rescale all fonts in image
    imAxes = event.AffectedObject;
    ud = get(imAxes,'UserData');
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
        set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end
end

function localSetupAxesListener(imAxes)
hgp     = findpackage('hg');
axesC   = findclass(hgp,'axes');
AxesListener = handle.listener(imAxes,axesC.findprop('Position'), ...
    'PropertyPostSet',{@localAxesListener,imAxes});
hFig = get(imAxes,'Parent');
    setappdata(hFig,'KenAxeResizeListeners',AxesListener);
end

function localAxesListener(~,~,imAxes) 
%     fprintf('Detected!!\n')
    % when called, rescale all fonts in image
    ud = get(imAxes,'UserData');
    fs = getBestFontSize(imAxes);
%     fprintf('%i\n',fs)
    for ii = 1:length(ud.handles)
        set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end   
end

function fs = getBestFontSize(imAxes)
% Try to keep font size reasonable for text
hFig        = get(imAxes,'Parent');
hFigFactor  = max(get(hFig,'Position') .* [0 0 0 1]);
hAxesFactor = max(get(imAxes,'Position') .* [0 0 0 1]);
axHeight    = diff(get(imAxes,'YLim'));
ud          = get(imAxes,'UserData');  % stored in teh first user data.
if ud.fs ~= round(ud.ratios * hFigFactor / axHeight)
%     fprintf('%f ~= %f',ud.fs,round(ud.ratios * hFigFactor / axHeight))
    fs = round(ud.ratios * hFigFactor / axHeight);
    fs = max(fs, 3);
elseif ud.fs ~= ud.ax_ratios / hAxesFactor;
%     fprintf('%f ~= %f',ud.fs,ud.ax_ratios / max(get(imAxes,'Position') .* [0 0 0 1]))
    fs = ud.fs * max(get(imAxes,'Position') .* [0 0 0 1]);
%     fprintf('fs = %f\n',fs)
    fs = max(fs, 3);
    ud.ratios = fs * axHeight / hFigFactor;
    set(imAxes,'UserData', ud);
else
    fs = ud.fs;
end
end