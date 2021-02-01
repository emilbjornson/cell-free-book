%This Matlab script can be used to reproduce Figure 4.6 in the monograph:
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


%Number of AP antennas
N = 8;

%Range of angular standard deviations in the local scattering model (in radians)
ASD_varphi_deg = [1e-5 1e-4 1e-3 1e-2 1e-1:1e-1:9e-1 1:1:9 10:5:60];
ASD_varphi = deg2rad(ASD_varphi_deg);
ASD_thetaRange = deg2rad([0 10]);

%Define the range of nominal angles of arrival
varphiRange = linspace(-pi,+pi,100);
theta = -pi/12;


%Define the effective SNR in (4.15), including processing gain
SNRdB = 10;
SNR = 10.^(SNRdB/10);

%Preallocate matrices for storing the simulation results
NMSE_localscattering = zeros(length(ASD_varphi),length(varphiRange),length(ASD_thetaRange));


%% Go through the range of nominal angles
for r = 1:length(varphiRange)
    
    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);    
    
    %Go through all ASDs
    for n = 1:length(ASD_varphi)
        
        for m = 1:length(ASD_thetaRange)
        
            %Compute the spatial correlation matrix using local scattering
            %model in (2.18) with Gaussian angular distribution
            R = functionRlocalscattering(N,varphiRange(r),theta,ASD_varphi(n),ASD_thetaRange(m));
            
            %Compute the NMSE according to (4.28)
            NMSE_localscattering(n,r,m) = real(trace(R - SNR*R*((SNR*R+eye(N))\R))/trace(R));
            
        end
               
    end
    
end

%Compute the NMSE for the uncorrelated fading case according to (4.28)
NMSE_uncorrelated = 1/(SNR+1);


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(ASD_varphi_deg,mean(NMSE_localscattering(:,:,2),2),'k-','LineWidth',2);
plot(ASD_varphi_deg,mean(NMSE_localscattering(:,:,1),2),'k:','LineWidth',2);
plot(ASD_varphi_deg,NMSE_uncorrelated*ones(size(ASD_varphi)),'r-.','LineWidth',2);
plot(ASD_varphi_deg,NMSE_localscattering(1)*ones(size(ASD_varphi)),'b--','LineWidth',2);

xlabel('ASD $(\sigma_{\varphi})$ [degree]','Interpreter','latex');
ylabel('NMSE','Interpreter','latex');
set(gca,'YScale','log');

legend({'Local scattering model: $\sigma_{\theta}=10^\circ$','Local scattering model: $\sigma_{\theta}=0^\circ$','Limit: Uncorrelated', 'Limit: Fully correlated'},'Interpreter','latex','Location','NorthEast');
ylim([0.01 1]);
