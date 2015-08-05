clear;
clc;
%% read the data
N=39;

%for N = 21:1:23;
    vidMat = load(['vid/h' num2str(N) '.mat'],'new');
    vidMat = vidMat.new;
%    weight_fft = cell(size(vidMat,1), size(vidMat,2));
    length = size(vidMat,3)-1000;
%% get average intensity
    int_average = mean(vidMat(:,:,1:length),3);

%     int_average=zeros(size(vidMat,1),size(vidMat,2));
%     for j = 1:1:size(vidMat,1),
%         for k = 1:1:size(vidMat,2),
%             int_average(j,k) = mean(vidMat(j,k,1:length));
%         end
%     end
%% second method to plot fft
 
    i = 1;
    
    n = length;
    Fs = 30;                    % Sampling frequency
    T = 1/Fs;                     % Sampling period
    L = n;                     % Length of ynal
    
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    
    f = Fs/2*linspace(0,1,NFFT/2+1); % frequency axis
    
    vidCrop = vidMat(40:120,40:85,500:4000); %crop the video
    
    intensity_vals = reshape(vidCrop,[numel(vidCrop(:,:,1)), size(vidCrop,3)]); %reshape such that intensity variation of each pixel is along each row
    intensity_vals = intensity_vals'; %transpose - intensity variation of each pixel is along each column
    DC = mean(intensity_vals,1); %Mean of intensities
    intensity_vals = bsxfun(@minus,intensity_vals,DC); %Subtracting the mean which removes the DC
    
    Y = fft(intensity_vals,NFFT)/L; %taking the FFT (fourier)
    Y = 2*abs(Y(1:1 + NFFT/2,:)); %Magnitude of fourier
    
    y_norms = sqrt(sum(intensity_vals.^2,1)); %norm of signal
    Y_weighted = bsxfun(@rdivide,Y,y_norms); % divide by norm to weight the fft
    
    Y_sum = sum(Y_weighted,2); % sum of weighter fft
    
    figure, plot(f, Y_sum) 
 
%   % select the leg area
%     for iy = 58:105,
%        for ix = 120:181,
%         %n = size(vidMat,3);
%         I = zeros(30,30,n,'uint8');
%         t = zeros(1,n);
%            
%         y = vidMat(ix,iy,1:length);
%         y = y-mean(y);
%         y = squeeze(y);
%         y = single(y);
%                  
%         Y = fft(y,NFFT)/L;  
%         f = Fs/2*linspace(0,1,NFFT/2+1);
%         fft_result = 2*abs(Y(1:NFFT/2+1));
%         
%         fft_matrix = [f; fft_result'];
%       if i == 1 
%          fft_tmp = fft_result / int_average(ix,iy);
%          i = i+1;
%       else
%           fft_tmp = fft_tmp + fft_result / int_average(ix,iy);
%       end
%        % final_fft(1,:) = f;
%         %final_fft(2,:) = final_fft(2,:) + fft_matrix(2,:) / int_average(ix,iy);
%           final_fft = [f;fft_tmp'];
%  %       h = figure;
%  %       figure, plot(f,2*abs(Y(1:NFFT/2+1)));
%  %       saveas(h,['figures/' num2str(6) num2str(iy) num2str(ix) 'pointfft.png']);
%        end
%     
%     end
%     
% % plot FFT result
% x = final_fft(1,:);
% y = final_fft(2,:);
% figure
% plot(x,y)
% ylim([0 700]);
 title('Adithya c, Machine 27, Speed = 4');
 xlabel('Frequency/Hz');
 ylabel('FFT result');


