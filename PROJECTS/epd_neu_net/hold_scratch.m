 
%% Supporting functions
%
function save_image_real_vs_random_dist(path,rand_accuracy,real_accuracy,jj,tt,zz, higher_than, p)

f= figure;

set(f,'color',[1 1 1])

h = histogram(rand_accuracy,10);

h.FaceAlpha = 0.4;

xl = xline(real_accuracy,'k:');

xl.LineWidth = 5;

title(sprintf('Real accuracy = %0.1f%%\nGreater than %0.1f%% of random accuracy (p = %0.3f)',100*real_accuracy,higher_than, p))

print(gcf,...
    
fullfile(path,['Compared_with_random_Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)]),...
    
'-dpng','-r300')

close all

end
