function txthilite(strcell, varargin)
% color and embolden text on matlab figures (default:red).
% strcell is a cell array of strings to be hilighted
% vargin is an optional colorspec to selection to color of the hilight
% (default:red).
% example:
%txthilight({'G2' 'G7' 'G8' 'G20'}, 'g') 
% -BQR 02/02/12
if ~isempty(varargin)
    color = varargin{1};
else
    color = 'r';
end

figH = findall(0,'Type','figure');
for ifig = 1:length(figH)
    textH = findall(figH(ifig),'type','text');
       for itxt = 1:length(textH)
           H = get(textH(itxt));
           if isfield(H,'String')
              if ismember(H.String,strcell)
                    set(textH(itxt),'Color',color,'FontWeight','bold')
              end
           end
       end
end
       

