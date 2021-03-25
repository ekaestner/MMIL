
load matfiles/avg_data_post.mat

badchanfile = '060325LP_badchans.txt';
badchan_i = ts_read_txt_badchans(badchanfile,{avg_data.sensor_info.label});

grad1 = 1:3:306;
grad2 = 2:3:306;
grad = union(grad1,grad2);

grad = 1:306;
grad = grad1;

grad = setdiff(grad,badchan_i);
data = avg_data.averages(3).data(grad,:);

eegplot(data,'dispchans',20,'winlength',1,'srate',1000)


