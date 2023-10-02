clear all;
close all;
clc;

%Test Data 10Hz

%Dataset 1
%filename = 'accel_dataset1run1.mat'; %30 steps
filename = 'dataset1take1.mat'; %50 steps

%Dataset 2
%filename = 'dataset2_take2.mat'; %50 steps
%filename = 'dataset2_take3.mat'; %60 steps

%Workshop 
%filename = 'sample_dataset2.mat';
%filename = "sample_dataset1.mat";
    

load(filename); 
%Acceleration is name of imported timetable
T = timetable2table(Acceleration);

time_data = table2array(T(:, 1));
ax_data = table2array(T(:, 2)); %x dimensional data
ay_data = table2array(T(:, 3)); %y dimensional data
az_data = table2array(T(:, 4)); %z dimensional data

a_norm = sqrt(ax_data.^2+ay_data.^2+az_data.^2); %norm of acceleration vector
a_norm = a_norm-mean(a_norm); %removes gravity 

% Lines below taken from MATLAB Documentation https://au.mathworks.com/help/matlab/ref/fft.html
X = a_norm;
Y = fft(X); %Computed DFT

Fs = 10;            % Sampling frequency (Hz)                  
T = 1/Fs;             % Sampling period       
L = length(time_data);  % Length of signal = no. datapoints
t = (0:L-1)*T;        % Time vector

% Double sided amplitude spectrum

%ffshift 
Yshift2 = abs(fftshift(Y));
fshift2 = linspace(-Fs/2, Fs/2, L);

%Single sided amplitude spectrum

Yshift1 = Yshift2(L/2+1:L);
fshift1 = fshift2(L/2+1:L);

% Get location of where highest peak occurs
highest_val = 0;
location = 0;
for i = 1:length(Yshift1)
    if (Yshift1(i) > highest_val)
        highest_val = Yshift1(i);   
        location = i;
    end
end
highest_val_freq = fshift1(location);


n = 2; % order of butterworth filter used
%Picking lower and upper bound of bandpass filter according to data
difference = 0.25; %making algorithm intellegent by always creating a filter with respect to walking frequency
lowbound = highest_val_freq-difference; %Hz
upperbound = highest_val_freq+difference; %Hz

%Conversion
w1 = lowbound*(2/Fs); 
w2 = upperbound*(2/Fs); 
[b,a] = butter(n, [w1,w2], "bandpass");
y = filter(b, a, a_norm);
        

% Count the number of steps taken

%Lines below taken from https://au.mathworks.com/help/matlabmobile/ug/counting-steps-by-capturing-acceleration-data.html

beta = 4;
minPeakHeight = std(y)-std(y)/beta; %Minimum height to be counted as a "step"

[pks,locs] = findpeaks(y, 'MINPEAKHEIGHT' ,minPeakHeight); %Return the location of peaks "loc" and height of each peak "pks"

%[pks, locs] = findpeaks(y, "MinPeakProminence", 0.4);

numSteps = numel(pks); % gets length of array pks which returns the number of steps
disp("The number of steps taken is: " + numSteps)

%Graphs

figure(1);
subplot(3, 1, 1);
plot(t, a_norm);
xlabel("Relative Time (s)");
ylabel({"Acceleration Magnitude", "Accounting For Gravity - m/s^2"});
title("Acceleration vs Relative Time - Corrupted Data");

subplot(3, 1, 2);
plot(t, y);
title("Filtered Data");
xlabel("Relative Time (s)"); 
ylabel({"Acceleration Magnitude", "Accounting For Gravity - m/s^2"});

subplot(3, 1, 3);
plot(t, y);
hold on;
plot(t(locs), pks, 'r', 'Marker', 'v', 'LineStyle', 'none');
title('Counting Steps');
xlabel('Relative Time (s)');
ylabel({"Acceleration Magnitude", "Accounting For Gravity - m/s^2"});
hold off;


figure(2);
subplot(2, 1, 1);
plot(fshift1, Yshift1);
hold on;
plot(highest_val_freq, highest_val, 'r*')
title("Single-Sided Amplitude Spectrum of Data - Before Filter");
xlabel("f (Hz)");
ylabel("YShift");

Fs = 10;            % Sampling frequency (Hz)                  
T = 1/Fs;             % Sampling period       
L = length(time_data);  % Length of signal = no. datapoints
t = (0:L-1)*T;        % Time vector


%ffshift 

Y = fft(y);
Yshift2 = abs(fftshift(Y));
fshift2 = linspace(-Fs/2, Fs/2, L);

%Single sided amplitude spectrum

Yshift1 = Yshift2(L/2+1:L);
fshift1 = fshift2(L/2+1:L);

subplot(2, 1, 2);
plot(fshift1, Yshift1);
title("Single-Sided Amplitude Spectrum of Data - After Filter");
xlabel("f (Hz)");
ylabel("YShift");

figure(3);
plot(t, ax_data, "b");
hold on;
plot(t, ay_data, "r");
hold on;
plot(t, az_data, "y");
xlabel("Time - Seconds");
ylabel("Acceleration - m/s^2");
title("Acceleration From 3 Axes");
legend("X Axis", "Y Axis", "Z Axis");
hold off;
