function [gainOverNoisedB_AP,gainOverNoisedB_BS,R_AP,R_BS,pilotIndex,BSassignment,D,D_small,APpositions,UEpositions,BSpositions,distances] = generateSetup_with_cellular(L,K,N,M,nbrOfSetups,seed,ASD_varphi,ASD_theta)
%This function generates realizations of the simulation setup described in
%Sections 5.3 and 5.4.3.
%
%INPUT:
%L                  = Number of APs for cell-free and small-cell systems
%K                  = Number of UEs in the network
%N                  = Number of antennas per AP for cell-free and small-cell
%                     systems
%M                  = Number of antennas per AP for Cellular Massive MIMO
%nbrOfSetups        = Number of setups with random UE and AP locations
%seed               = Seed number of pseudorandom number generator
%ASD_varphi         = Angular standard deviation in the local scattering model 
%                     for the azimuth angle (in radians)
%ASD_theta          = Angular standard deviation in the local scattering model
%                     for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB_AP = Matrix with dimension L x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by the
%                     noise variance) between AP l and UE k in setup n for
%                     cell-free and small-cell setups
%gainOverNoisedB_BS = Matrix with dimension nbrBSs x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by the
%                     noise variance) between AP l and UE k in setup n for
%                     Cellular Massive MIMO setup
%R_AP               = Matrix with dimension N x N x L x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%                     variance for cell-free and small-cell setups
%R_BS               = Matrix with dimension N x N x nbrBSs x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%                     variance for Cellular Massive MIMO setup
%pilotIndex         = Matrix with dimension K x nbrOfSetups containing the
%                     pilots assigned to the UEs
%BSassignment       = Matrix with dimension K x nbrOfSetups containing the
%                     AP indices assigned to the UEs for Cellular Massive
%                     MIMO system
%D                  = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                     is one if AP l serves UE k in setup n and zero otherwise
%                     for cell-free setup
%D_small            = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                     is one if AP l serves UE k in setup n and zero otherwise
%                     for small-cell setup
%APpositions        = Vector of length L with the AP locations, where the real
%                     part is the horizontal position and the imaginary part
%                     is the vertical position for cell-free and small-cell
%                     systems
%UEpositions        = Vector of length K with UE positions, measured in the
%                     same way as APpositions
%distances          = Matrix with same dimension as gainOverNoisedB_AP containing
%                     the distances in meter between APs and UEs for
%                     cell-free nad small-cell setups
%
%This Matlab function was developed to generate simulation results to:
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

%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>5)&&(seed>0)
    rng(seed)
end

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance



%Number of cellular APs (BSs)
nbrBSs = 4;


%Number of pilots is equal to the number of UEs per cell
tau_p = K/nbrBSs;

%Number of cellular BSs per dimension on the grid
nbrBSsPerDim = sqrt(nbrBSs);

%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);



%Number of APs per dimension on the grid for cell-free and small-cell
%setups
nbrAPsPerDim = sqrt(L);

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

%Deploy APs on the grid
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

%Compute alternative BS and AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Prepare to save results
gainOverNoisedB_AP = zeros(L,K,nbrOfSetups);
gainOverNoisedB_BS = zeros(nbrBSs,K,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
R_BS = zeros(M,M,nbrBSs,K,nbrOfSetups);
R_AP = zeros(N,N,L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
BSassignment = zeros(K,nbrOfSetups);
masterAPs = zeros(K,1);
D = zeros(L,K,nbrOfSetups);
D_small = zeros(L,K,nbrOfSetups);



%% Go through all setups
for n = 1:nbrOfSetups
    
    
    
    %Prepare to store number of UEs per cellular BS
    nbrOfUEsPerBS = zeros(nbrBSs,1);
    %Initialize the number of dropped UEs in the network
    nbrOfUEs = 0;
    %Prepare to compute UE locations
    UEpositions = zeros(K,1);
    
    
    
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);
    shadowBSrealizations = zeros(K,nbrBSs);
    
    
    
    %Add UEs
    while nbrOfUEs<K
        
        %Generate a random UE location in the area
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichposAP] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances(:,nbrOfUEs+1,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        
        %If this is not the first UE
        if nbrOfUEs>0
            
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(nbrOfUEs,1);
            
            for i = 1:nbrOfUEs
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end
            
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:nbrOfUEs,1:nbrOfUEs);
            meanvaluesAP = term1*shadowAPrealizations(1:nbrOfUEs,:);
            meanvaluesBS = term1*shadowBSrealizations(1:nbrOfUEs,:);
            
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            
        else %If this is the first UE
            
            %Add the UE and begin to store shadow fading correlation values
            meanvaluesAP = 0;
            meanvaluesBS = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
            
        end
        
        %Generate the shadow fading realizations for cell-free and
        %small-cell setups
        shadowing = meanvaluesAP + stdvalue*randn(1,L);
        
        %Compute the channel gain divided by noise power (in dB)
        gainsAP = constantTerm - alpha*log10(distances(:,nbrOfUEs+1,n)) + shadowing' - noiseVariancedBm;
        
        %Determine the master AP for dropped UE  by looking for AP with best
        %channel condition
        [~,mastertemp] = max(gainsAP);
        
        %Assign orthogonal pilots to the first tau_p UEs
        if nbrOfUEs+1 <= tau_p
            
            pilotIndextemp = nbrOfUEs+1;
            
        else %Assign pilot for remaining UEs
            
            %Compute received power from to the master AP from each pilot
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB_AP(mastertemp,pilotIndex(1:nbrOfUEs,n)==t,n)));
                
            end
            
            %Find the pilot with the least receiver power according to
            %Algorithm 4.1
            [~,bestpilot] = min(pilotinterference);
            pilotIndextemp = bestpilot;
            
        end
        
        %Generate the shadow fading realizations for cellular setup
        
        shadowingBS = meanvaluesBS + stdvalue*randn(1,nbrBSs);
        
        
        %Compute distances from the UE to each of the BSs
        [distanceBSstoUE,whichposBS] = min(abs(BSpositionsWrapped - repmat(UEposition,size(BSpositionsWrapped))),[],2);
        distancesBS = sqrt(distanceVertical^2+distanceBSstoUE.^2);
        
        %Compute the channel gain divided by noise power (in dB)
        gainsBS = constantTerm - alpha*log10(distancesBS) + shadowingBS' - noiseVariancedBm;
        
        %Find which BS the UE would like to connect to in cellular setup
        [~,bestBS] = max(gainsBS);
        
        %If the BS doesn't have tau_p UEs yet
        if nbrOfUEsPerBS(bestBS)<tau_p
            
            if sum(BSassignment(pilotIndex(:,n)==pilotIndextemp,n)==bestBS)<1
                
                
                %Add the UE to the preferred cell
                nbrOfUEsPerBS(bestBS) = nbrOfUEsPerBS(bestBS) + 1;
                
                %Compute the number of dropped UEs until now and and store the UE index
                nbrOfUEs = sum(nbrOfUEsPerBS);
                k = nbrOfUEs;
                
                %Store the channel gains divided by noise power
                gainOverNoisedB_AP(:,k,n) = gainsAP;
                gainOverNoisedB_BS(:,k,n) = gainsBS;

                
                
                %Update shadowing correlation matrix and store realizations
                shadowCorrMatrix(1:k-1,k) = newcolumn;
                shadowCorrMatrix(k,1:k-1) = newcolumn';
                shadowAPrealizations(k,:) = shadowing;
                shadowBSrealizations(k,:) = shadowingBS;
                
                %Assign the UE to the cellular BS
                BSassignment(k,n) = bestBS;
                
                %Store the UE position
                UEpositions(k) = UEposition;
                
                
                %Determine the master AP for UE k by looking for AP with best
                %channel condition
                D(mastertemp,k,n) = 1;
                masterAPs(k) = mastertemp;
                
                %Store pilot assignment
                pilotIndex(k,n) = pilotIndextemp;
                %Go through all APs
                for l = 1:L
                  
                    %Compute nominal angle between UE k and AP l
                    angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichposAP(l)));
                    angletoUE_theta = asin(distanceVertical/distances(l,k,n));
                    %Generate spatial correlation matrix using the local
                    %scattering model in (2.18) and Gaussian angular distribution
                    %by scaling the normalized matrices with the channel gain
                    if nargin>6
                        R_AP(:,:,l,k,n) = db2pow(gainOverNoisedB_AP(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                    else
                        R_AP(:,:,l,k,n) = db2pow(gainOverNoisedB_AP(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
                    end
                end
                
                %Go through all BSs
                for l = 1:nbrBSs
                    
                    %Compute nominal angle between the new UE k and BS l
                    angletoUE_varphi = angle(UEpositions(k)-BSpositionsWrapped(l,whichposBS(l)));
                    angletoUE_theta = asin(distanceVertical/distancesBS(l));
                    %Generate spatial correlation matrices using the
                    %local scattering model
                    if nargin>6
                        R_BS(:,:,l,k,n) = db2pow(gainOverNoisedB_BS(l,k,n))*functionRlocalscattering(M,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                    else
                        R_BS(:,:,l,k,n) = db2pow(gainOverNoisedB_BS(l,k,n))*eye(M);
                    end
                    
                end
                
                
                
            end
        end
        
    end
    
    
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots in the cell-free setup
    for l = 1:L
        
        for t = 1:tau_p
            
            pilotUEs = find(t==pilotIndex(:,n));
            [~,UEindex] = max(gainOverNoisedB_AP(l,pilotUEs,n));
            D(l,pilotUEs(UEindex),n) = 1;
            
            
            
        end
        
    end
    
    %Determine the AP serving each UE in the small-cell setup according to
    %(5.47) by considering only the APs from the set M_k for UE k, i.e.,
    %where D(:,k,n) is one.
    for k=1:K
        
        tempmat = -inf*ones(L,1);
        tempmat(D(:,k,n)==1,1) = gainOverNoisedB_AP(D(:,k,n)==1,k,n);
        [~,servingAP] = max(tempmat);
        D_small(servingAP,k,n) = 1;
        
    end
    
    
end

