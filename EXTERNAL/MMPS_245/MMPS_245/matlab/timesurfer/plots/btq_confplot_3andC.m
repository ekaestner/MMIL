function varargout = btq_confplot_3andC(x,y,err,c)
%BTQ_CONFPLOT_3ANDC Linear plot with continuous confidence/error boundaries.
%
%   BTQ_CONFPLOT_3ANDC(X,Y,ERR,C) plots the graph of vector X vs. vector Y with
%   'continuous' confidence/error boundaries specified by the vector
%   ERR in the color specified by C.
%
%   H = BTQ_CONFPLOT_3ANDC(...) returns a vector of line handles.
%
%   For example,
%      x = 1:0.1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      btq_confplot_3andC(x,y,e,[1 0 0])
%   draws symmetric continuous confidence/error boundaries of unit standard deviation.
%

if (nargin<4)
    disp('ERROR: not enough input arguments!');
    return;
end % if

z1 = y + err;
z2 = y - err;

p = plot(x,y,x,z1,x,z2);    YLIM = get(gca,'YLim');    delete(p);

X=[x fliplr(x)]';Y=[z1 fliplr(z2)]';
pa1=patch(X,Y,c,'FaceAlpha',0.4,'LineStyle','none');set(gca,'XLim',[min(x) max(x)]);
set(get(get(pa1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
hold on;
p = plot(x,y,'Color',c,'LineWidth',2);
hold off;

set(gca,'Layer','top');

H = [p, pa1];

if (nargout>0) varargout{1} = H; end;
