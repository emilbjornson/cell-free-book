%This Matlab script can be used to reproduce Figure 1.2 in the monograph:
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

%Length of the coverage area in meter (as a square)
squareLength = 1000;

%Total number of APs in the area
nbrAPs = 4;

%Range of number of antennas and antenna gains
numberOfAntennas = [1 64];
antennaGains = [8 1];

%APs per dimension
nbrAPsPerDim = sqrt(nbrAPs);

%Set the cell edge SNR to 0 dB) in case (a), which becomes -9 dB (1/8) when
%not accounting for the antenn gain
SNRedge = 1/8;

%Number of UE locations in horizontal/vertical dimension
Kdims = 20;

%Number of UEs per cell
K = Kdims*Kdims;

%Pathloss exponent
alpha = 4;

%Maximum spectral efficiency
maxSE = 8;

%Bandwidth in MHz
bandwidth = 10;


%Distance between APs in horizontal/vertical dimension
interSiteDistance = squareLength/nbrAPsPerDim;

%Put out the APs on a square grid
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);


%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrAPs 1]);


%Prepare to put out the UEs in the cells
UEpositions = zeros(K,nbrAPs);
perBS = zeros(nbrAPs,1);

%Prepare to compute SEs
SEs = zeros(K,nbrAPs,length(numberOfAntennas));


%Go through all the cells
for l = 1:nbrAPs
    
    %Put out K users in each cell
    while min(perBS(l))<K
        
        posTmp = repmat(linspace(0,1,Kdims),[Kdims 1]);
        posTmp2 = posTmp';
        
        posX = interSiteDistance*posTmp(:) + real(BSpositions(l)) - interSiteDistance/2;
        posY = interSiteDistance*posTmp2(:) + imag(BSpositions(l)) - interSiteDistance/2;
        
        UEpositions(:,l) = posX  + 1i*posY;
        perBS(l) = K;
        
    end
    
    %Compute the distance from the UEs to AP l
    distancesSquaredBSl = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(l,:),[K 1]) ),[],2);
    
    %Compute the signal power for each UE
    signalPower = ((interSiteDistance/2)./distancesSquaredBSl).^(alpha);
    
    %Compute the interference power for each UE
    interferencePower = zeros(size(distancesSquaredBSl));
    
    for j = 1:nbrAPs
        
        %Compute the distance from the users to UE j
        distancesSquaredBSj = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Add interference from non-serving APs
        if j~=l
            
            interferencePower = interferencePower + ((interSiteDistance/2)./distancesSquaredBSj).^(alpha);
            
        end
        
    end
    
    %Compute SE for different number of antennas, assuming i.i.d. Rayleigh
    %fading and perfect CSI
    for m = 1:length(numberOfAntennas)
        
        SEs(:,l,m) = log2(1+numberOfAntennas(m)*antennaGains(m)*signalPower./(antennaGains(m)*interferencePower+1/SNRedge));
        SEs(SEs(:,l,m)>maxSE,l,m) = maxSE;
        
    end
    
end


%Plot the simulation results
for m = 1:length(numberOfAntennas)
    
    figure;
    hold on; box on;
    
    for l = 1:nbrAPs
        
        surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs(:,l,m),[Kdims Kdims]));
        
    end
    
    xlim([0 1000]);
    ylim([0 1000]);
    zlim([0 bandwidth*maxSE]);
    xlabel('Position [m]','Interpreter','latex');
    ylabel('Position [m]','Interpreter','latex');
    zlabel('Data rate [Mbit/s/user]','Interpreter','latex');
    caxis([0,bandwidth*maxSE]);
    view(-52,28);
    set(gca,'fontsize',16);
    
end