function subhandles = btq_panels(nvert,nhoriz,pos)
% function subhandles = btq_panels(nvert,nhoriz,pos)
% This function creates several subplots in a window that are directly
% attached to each other (no gap inbetween as in the standard subplot
% function). Input parameters are NVERT, the number of axes systems in the
% vertical direction, and NHORIZ, the number of axes systems in the
% horizontal direction. %
% These handles are stored in the vector SUBHANDLES. The
% numbering starts in the top left subplot (No 1) and continues row wise
% until the last subplot at the bottom right (No NVERT*NHORIZ). The
% following script illustrates the numbering scheme:
% 
% figure
% nv = 3;
% nh = 4;
% subhandles = btq_panels(nv,nh,'top');
% % To number the plots:
% for i=1:nv
%     for j=1:nh
%         subplot(subhandles((i-1)*nh+j))
%         text(0.5,0.5,num2str((i-1)*nh+j))
%     end
% end
% 

% Set the bottom position of the plots:
if nargin ==2
    bot = 0.05
else
    if pos=='top'
        bot =0.05;
    else
        bot = 0.05;
    end
end
left = 0.05;

% The width and height of the subplots
width = 0.9/nhoriz;
height = 0.9/nvert;

% Place the subplots
for i=1:nvert
    for j=1:nhoriz
        subhandles((i-1)*nhoriz + j)=subplot('position',...
            [left+(j-1)*width bot+(nvert-i)*height width height]);
    end
end
