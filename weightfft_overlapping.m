clear;
clc;
%% read the data
N = 40;

%for N = 21:1:23;
    vidMat = load(['vid/h' num2str(N) '.mat'],'new');
    vidMat = vidMat.new;
    l = floor(size(vidMat,3));
    length = ((l / 100) - 10) * 100;

%% get average intensity
    int_average = zeros(size(vidMat,1),size(vidMat,2));
    for j = 1:1:size(vidMat,1),
        for k = 1:1:size(vidMat,2),
            int_average(j,k) = mean(vidMat(j,k,1:length));
        end
    end


%% overlapping window fft
n = 200;
j = 1;
I = zeros(30,30,n,'uint8');
t = zeros(1,n);
Fs = 30;                    % Sampling frequency
T = 1/Fs;                   % Sample time
L = n;                      % Length of ynal
NFFT = 2^nextpow2(L);       % Next power of 2 from length of y
Len = NFFT/2+1;
f = Fs/2*linspace(0,1,Len);
window_fft = zeros(Len, 1);

% scan 200 frames, every 100 frames
for it = 500:100:length
  % select the leg area
    fft_tmp = zeros(Len, 1);
    for iy = 40:1:85,
       for ix = 40:1:120,
           
            y = vidMat(ix,iy,it:it+200);
            y = y-mean(y);
            y = squeeze(y);
            y = single(y);

            Y = fft(y,NFFT) / L;

            fft_result = 2 * abs(Y(1:Len));
        
    %    fft_matrix = [f; fft_result'];
            if int_average(ix,iy) == 0
                % avoid to be nan
                continue
            end
            fft_tmp = fft_tmp + fft_result / int_average(ix,iy);
        end
    end
    
    window_fft = window_fft + fft_tmp;
    j = j + 1;
end

% get the average
window_fft_mean = window_fft / (j - 1);
final_fft = [f; window_fft_mean'];
                       
% plot FFT result
x = final_fft(1,:);
y = final_fft(2,:);
figure
plot(x,y)
ylim([0 900]);
title('overlapping, Yoshitaka h, Machine 27, Speed = 5');
xlabel('Frequency/Hz');
ylabel('FFT result');


