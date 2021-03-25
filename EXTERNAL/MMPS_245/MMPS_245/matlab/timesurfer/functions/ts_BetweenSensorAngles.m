function theta = BetweenSensorAngles(sens,zshift)
% Purpose: calculate the angle between sensors
% Input: sens - TimeSurfer sensor_info structure
% Output: theta - the angle between each pair of sensors
%
% This script calculates the angle (in 3d) between MEG sensors.
% written nima dehghani MD , june 29th 2009, nima@ucsd.edu
% 
% last modified on 07-Sep-2009 by jason sherfey
if nargin < 2, zshift = 0.08; end
x = arrayfun(@(x)(x.loc(1,4)),sens,'UniformOutput',true)';
y = arrayfun(@(x)(x.loc(2,4)),sens,'UniformOutput',true)';
z = arrayfun(@(x)(x.loc(3,4)),sens,'UniformOutput',true)';

z2 = z - zshift; % shifts the Z center to almost the middle of the helmet
a2 = [x y z2];
% plot3(x,y,z2,'.')
% for i=1:length(x)
%     text(x(i), y(i), z2(i),num2str(i))
% end
% grid on
% vline(0,'g')
% hline(0,'r')

for i=1:length(x)
  a = a2(i,:);
  for j =1:length(x)
      b = a2(j,:);
      if i== j
         theta(i,j) = nan;
      else
         theta(i,j) =  acos(dot(a,b)/(norm(a)*norm(b)))*180/pi;
      end	  
  end 
end
theta = real(theta);
