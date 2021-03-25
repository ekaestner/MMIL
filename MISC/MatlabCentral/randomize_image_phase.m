function [ImScrambled RandomPhase] = randomize_image_phase(Im, likeOriginal,randPhase2Use,phase_mask)
%close all
%Im = mat2gray(double(imread('C:\Users\mmil\Desktop\celeb_best4_front\crop_noise_rotated_lessface_edge/celeb_aniston_1.jpg')));
%% saliency map may help randomize phase
% saliency_map_resolution = [size(Im,1) size(Im,2)];                 
% the target saliency map resolution; the most important parameter for spectral saliency approaches
% smap_smoothing_filter_params = {'gaussian',20,5}; 
% filter parameters for the final saliency map
% cmap_smoothing_filter_params = {};                 
% optionally, you can also smooth the conspicuity maps
% cmap_normalization = 1;                             
% specify the normalization of the conspicuity map here
% extended_parameters = {};                           
% @note: here you can specify advanced algorithm parameters for the selected algorithm, e.g. the quaternion axis
% do_figures = false;                                
% enable/disable spectral_saliency_multichannel's integrated visualizations
% saliency_map = spectral_saliency_multichannel(Im,saliency_map_resolution,'fft',smap_smoothing_filter_params,cmap_smoothing_filter_params,cmap_normalization,extended_parameters,do_figures);

Im = mat2gray(double(Im));

%read and rescale (0-1) image
num_rands = 1;
ImSize = size(Im);
%[b1_orig g2_orig] = weibull_constrast_image_fit(rgb2gray(Im));
for rands = 1:num_rands
if ~exist('randPhase2Use', 'var') || isempty(randPhase2Use)
    RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
%     RandomPhase(saliency_map < 1e-6) = 0;
else
    RandomPhase = randPhase2Use;
end
if exist('phase_mask', 'var')
    RandomPhase(~phase_mask) = 0;
end
% RandomPhase(saliency_map < .01) = 0;
%generate random phase structure

for layer = 1:ImSize(3)
    ImFourier(:,:,layer) = fft2(Im(:,:,layer));       
%Fast-Fourier transform
    Amp(:,:,layer) = abs(ImFourier(:,:,layer));       
%amplitude spectrum
    Phase(:,:,layer) = angle(ImFourier(:,:,layer));   
%phase spectrum
    Phase(:,:,layer) = Phase(:,:,layer) + RandomPhase - RandomPhase*likeOriginal; 
    % likeOriginal between 0 and 1, 1 is most like the original image
%add random phase to original phase
    ImScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(Phase(:,:,layer))));   
%combine Amp and Phase then perform inverse Fourier
end
% imagesc(RandomPhase)
ImScrambled = real(ImScrambled); 
%get rid of imaginery part in image (due to rounding error)
% imtool(Phase)
% imtool(ImScrambled)

%imwrite(ImScrambled,'PhaseScrambled.jpg','jpg');
%     if plot_flag
%         figure(1)
%         sfPlot(Im,1);
%         hold on;
%         sfPlot(ImScrambled,1);
% 
%         %imshow(Im)
%         [b1_scramble g2_scramble] = weibull_constrast_image_fit(rgb2gray(ImScrambled));
% 
%         figure(2);
%         scatter(b1_orig, g2_orig);
%         text(b1_orig, g2_orig,'ORIG');
%         hold all;
%         scatter(b1_scramble, g2_scramble);
%         hold on;           
%             colormap('gray')%# Add to the plot
%             figure(3)
%             imshow(ImScrambled);
%         %text(b1_scramble, g2_scramble,'SCRAMBLED');
%     end
end

