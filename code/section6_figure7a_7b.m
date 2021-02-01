%This Matlab script can be used to reproduce Figures 6.7(a) and 6.7(b) in the monograph:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Empty workspace and close figures
close all;
clear;

%Number of antennas
N = 2;

SNR = 1;

%Number of channel realizations
nbrOfRealizations = 10000000;

%Generate CN(0,1) realizations
h1 = (randn(N,nbrOfRealizations) + 1i*randn(N,nbrOfRealizations))/sqrt(2);
h2 = (randn(N,nbrOfRealizations) + 1i*randn(N,nbrOfRealizations))/sqrt(2);


%Prepare to save results
signal_gain_MR = zeros(nbrOfRealizations,1);
signal_gain_LP_MMSE = zeros(nbrOfRealizations,1);
interference_gain_MR = zeros(nbrOfRealizations,1);
interference_gain_LP_MMSE = zeros(nbrOfRealizations,1);
norm2_MR = zeros(nbrOfRealizations,1);
norm2_LP_MMSE = zeros(nbrOfRealizations,1);



for n = 1:nbrOfRealizations
    
    %Display the progress
    if mod(n,10000) == 0
        n        
    end
    
    %Compute precoding vector to UE 1
    w_MR = h1(:,n); %MR without normalization
    w_LP_MMSE = (h1(:,n)*h1(:,n)'+h2(:,n)*h2(:,n)'+(1/SNR)*eye(N))\h1(:,n); %LP-MMSE without normalization
    
    %Compute the squared norm to be used for normalization
    norm2_MR(n) = norm(w_MR)^2;
    norm2_LP_MMSE(n) = norm(w_LP_MMSE)^2;
    
    %Compute signal power (the imaginary part is zero in theory)
    signal_gain_MR(n) = SNR*abs(h1(:,n)'*w_MR)^2;
    signal_gain_LP_MMSE(n) = SNR*abs(h1(:,n)'*w_LP_MMSE)^2;
    
    %Compute interference power (the imaginary part is zero in theory)
    interference_gain_MR(n) = SNR*abs(h2(:,n)'*w_MR)^2;
    interference_gain_LP_MMSE(n) = SNR*abs(h2(:,n)'*w_LP_MMSE)^2;
    
end

signal_gain_MR = signal_gain_MR/mean(norm2_MR);
signal_gain_LP_MMSE = signal_gain_LP_MMSE/mean(norm2_LP_MMSE);
interference_gain_MR = interference_gain_MR/mean(norm2_MR);
interference_gain_LP_MMSE = interference_gain_LP_MMSE/mean(norm2_LP_MMSE);


%% Plot the simulation results
% Plot Figure 6.7(a)
figure; box on;
h1 = histogram(signal_gain_MR);
hold on
h2 = histogram(signal_gain_LP_MMSE);
h1.Normalization = 'pdf';
h1.BinWidth = 0.05;
h2.Normalization = 'pdf';
h2.BinWidth = 0.05;
set(gca,'fontsize',16);
xlabel('Signal-to-noise ratio','Interpreter','Latex');
ylabel('PDF','Interpreter','Latex');
xlim([0,10]);

% Plot Figure 6.7(b)
figure; box on;
h1 = histogram(interference_gain_MR);
hold on
h2 = histogram(interference_gain_LP_MMSE);
h1.Normalization = 'pdf';
h1.BinWidth = 0.025;
h2.Normalization = 'pdf';
h2.BinWidth = 0.025;
set(gca,'fontsize',16);
xlabel('Interference-to-noise ratio','Interpreter','Latex');
ylabel('PDF','Interpreter','Latex');
xlim([0,5]);
