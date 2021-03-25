% This script downloads survey results from a survey of some 250,000
% participants who described colors on their computer monitors, using 
% their own words. After you run this script (only takes a second), you'll
% forever be able to make plots like this 
% 
% x = -pi:.1:pi;
% y = sin(x);
% plot(x,y,'linewidth',4,'color',rgb('lavender'))
% 
% 
% A description of the XKCD survey results can be found here:
% http://blog.xkcd.com/2010/05/03/color-survey-results/ 
% 

if ~exist('rgb.txt','file')
    try
    urlwrite('http://xkcd.com/color/rgb.txt','rgb.txt');
    catch 
        disp('Having trouble downloading the data file. You''ll need to download it manually')
        disp('from http://xkcd.com/color/rgb.txt and place it in your current folder.')
        return
    end
end
    
fid = fopen('rgb.txt'); 
RGB = textscan(fid,'%s %s','delimiter','\t');
fclose(fid);

colorlist = RGB{1}; 
hex = RGB{2};

rgblist = hex2rgb(char(hex));

save('xkcd_rgb_data.mat','colorlist','rgblist')

