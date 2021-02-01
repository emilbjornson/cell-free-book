%This Matlab script can be used to reproduce Figure 4.7 in the monograph:
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


%Define the range of AP antennas
Nvalues = 1:16;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(10);
ASD_theta = deg2rad(10);

%Define the range of nominal angles of arrivals
varphiRange = linspace(-pi,+pi,100);
theta = -pi/12;


%Define the range of effective SNRs in (4.15), including processing gain
SNRdB = [0, 10, 20];
SNR = 10.^(SNRdB/10);


%Preallocate matrices for storing the simulation results
NMSE = zeros(length(Nvalues),length(varphiRange),length(SNR));


%% Go through the range of nominal angles
for r = 1:length(varphiRange)

    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);    
    
    %Compute the spatial correlation matrix using local scattering
    %model in (2.18) with Gaussian angular distribution
    R = functionRlocalscattering(max(Nvalues),varphiRange(r),theta,ASD_varphi,ASD_theta);
    
    
    %Go through all SNRs
    for n = 1:length(SNR)
        
        %Go through all number of antennas
        for mind = 1:length(Nvalues)
            
            %Pick out the correlation matrix for current number of antennas
            Rmatrix = R(1:Nvalues(mind),1:Nvalues(mind));
            
            %Compute the NMSE according to (4.28)
            NMSE(mind,r,n) = real(trace(Rmatrix - SNR(n)*Rmatrix*((SNR(n)*Rmatrix+eye(Nvalues(mind)))\Rmatrix)))/trace(Rmatrix);
            
        end
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'fontsize',16);
plot(Nvalues,mean(NMSE(:,:,1),2),'k-','LineWidth',2);
plot(Nvalues,mean(NMSE(:,:,2),2),'r-.','LineWidth',2);
plot(Nvalues,mean(NMSE(:,:,3),2),'b--','LineWidth',2);
xlabel('Number of AP antennas ($N$)','Interpreter','Latex');
ylabel('NMSE','Interpreter','Latex');
set(gca,'YScale','log');
legend({'SNR = 0 dB','SNR = 10 dB','SNR = 20 dB'},'Interpreter','Latex','Location','SouthWest');
xlim([1 16]);
ylim([1e-3 1]);
