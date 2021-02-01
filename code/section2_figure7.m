%This Matlab script can be used to reproduce Figure 2.7 in the monograph:
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

%Set the side length of the simulation area
squareLength = 400;

%Total number of APs
nbrOfAPs = 64;

%Number of realizations of the random AP locations
nbrOfSetups = 100000;

%Define a function to compute the large-scale fading as function of the
%horizontal distance  measured in meter. The AP is 10 meter above the UE.
largescaleFading = @(hor_dist) db2pow(-30.5-36.7*log10(sqrt(hor_dist.^2+10^2)));

%Set the standard deviation of the shadow fading
shadowFadingStd = 4;


%The UE is in the center of the simulation area
UElocation = (squareLength/2)*(1+1i);


%Prepare to save the simulation results
variance = zeros(nbrOfSetups,length(nbrOfAPs));



%% Go through all random realizations of the AP locations
for n = 1:nbrOfSetups
    
    %Generate random AP locations with uniform distribution
    APcellfree = squareLength*(rand(nbrOfAPs,1)+1i*rand(nbrOfAPs,1));
    
    %Compute distances to the APs
    distances = abs(APcellfree - UElocation);
    
    %Compute large-scale fading based on the model in (2.16)
    beta = largescaleFading(distances) .* 10.^(shadowFadingStd*randn(nbrOfAPs,1)/10);
    
    %Compute the variance value in (2.26) without the 1/N term
    variance(n) = sum(beta)/(sum(sqrt(beta)))^2;
    
end

%Range of number of AP antennas
Nrange = 1:10;

%Compute the median with respect to the AP locations
medianVariance = median(variance);

%Compute an interval that contains 90% of the realizations
varianceSorted = sort(variance,'ascend');
lowLimit = medianVariance - varianceSorted(round(0.05*nbrOfSetups));
highLimit = varianceSorted(round(0.95*nbrOfSetups)) - medianVariance;


%% Plot simulation results
figure;
hold on; box on;
errorbar(Nrange,medianVariance./Nrange,lowLimit./Nrange,highLimit./Nrange,'b--','LineWidth',2);
plot(Nrange,(1/64)*ones(size(Nrange)),'k:','LineWidth',2);
xlabel('Number of antennas ($N$)','Interpreter','latex');
ylabel('Variance of channel hardening','Interpreter','latex');
legend({'Cell-free','Reference case'},'Interpreter','latex','Location','NorthEast');
set(gca,'fontsize',16);
ylim([0 0.3]);
