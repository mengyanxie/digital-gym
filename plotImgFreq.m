clear;
clc;
%% read the data
num_frames = 300;

%for N = 21:1:23;
    vidMat = load('vid/h41.mat');%(['vid/c' num2str(N) '.mat'],'gvid');
    vidMat = vidMat.new;

    vidMat = vidMat(:,:,1:1000);
    length = size(vidMat,3);

    
    
    n = length;
    Fs = 30;                    % Sampling frequency
    T = 1/Fs;                     % Sampling period
    L = n;                     % Length of ynal
    
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    
    f = Fs/2*linspace(0,1,NFFT/2+1); % frequency axis
    
    vidCrop = vidMat; %crop the video
    
    intensity_vals = reshape(vidCrop,[numel(vidCrop(:,:,1)), size(vidCrop,3)]); %reshape such that intensity variation of each pixel is along each row
    intensity_vals = intensity_vals'; %transpose - intensity variation of each pixel is along each column
    DC = mean(intensity_vals,1); %Mean of intensities
    intensity_vals = bsxfun(@minus,intensity_vals,DC); %Subtracting the mean which removes the DC
    
    Y = fft(intensity_vals,NFFT)/L; %taking the FFT (fourier)
    Y = 2*abs(Y(1:1 + NFFT/2,:)); %Magnitude of fourier
    
    [M,I] = max(Y,[],1);
    f_peaks = f(I);
    f_img = reshape(f_peaks,size(vidCrop(:,:,1)));
    figure, imagesc(f_img); axis image; colormap jet
    
    %%
    %f_img_1 = f_img./max(f_img(:));
    f_img_1 = f_img;
    vidNew = permute(vidCrop,[1 2 4 3]);
    vidNew = cat(3, repmat(f_img_1,[1 1 1 1000]), vidNew, vidNew);
    
    implay(vidNew);