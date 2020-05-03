function [RFOutput]= RFMLOC(ContrailsInfo, CloudsInfo, AtmosData)
% RFMLOC Calculate RF from multiple-layer overlapping clouds
%
% Input data: contrails and clouds properties & atmospheric data
% 
% Output data: RF for only clouds, only contrails, all layers and 
% independent contrails as if no overlapping

RFOutput = struct;


% Setup layers characteristics
Contrails_RF_data = zeros(1,5);
for ilayer = 1:size(ContrailsInfo,1)
    Contrails_RF_data(ilayer,:) = [ContrailsInfo(ilayer).level ContrailsInfo(ilayer).area ...
        ContrailsInfo(ilayer).T ContrailsInfo(ilayer).od ContrailsInfo(ilayer).g];
end
Clouds_RF_data = zeros(1,5);
for ilayer = 1:size(CloudsInfo,1)
    Clouds_RF_data(ilayer,:) = [CloudsInfo(ilayer).level CloudsInfo(ilayer).area ...
        CloudsInfo(ilayer).T CloudsInfo(ilayer).od CloudsInfo(ilayer).g];
end

% Setup atmospheric data
OLR = AtmosData.OLR;
SDR = AtmosData.SDR;
t = AtmosData.t;
alpha = AtmosData.alpha;
mu = AtmosData.mu;

% Setup parameters
Sigma=1.607e-4;
delta = 0.75;
k = 2.528;

% Remove empty layers
Clouds_RF_data = Clouds_RF_data(Clouds_RF_data(:,4)>0,:);
Contrails_RF_data = Contrails_RF_data(Contrails_RF_data(:,4)>0,:);

% Remove inconsistent effective optical depth when Sun in the horizon
if abs(mu)<0.01
    mu = 1;
end


% Start of calculations
% 4 cases are evaluated (4 output values):
%   1. Only clouds overlapping
%   2. Contrails overlapping without clouds
%   3. Clouds and contrails overlapping
%   4. Clouds and independent contrails

for icase = 1:4
    data = zeros(3,1); RF_LW = 0; RF_SW = 0;
    if icase == 1
        % overlapping clouds
        RF_data = Clouds_RF_data;
        Ntimes = 1;
    elseif icase == 2
        % overlapping contrails
        RF_data = Contrails_RF_data;
        Ntimes = 1;
    elseif icase == 3
        % all layers
        RF_data = [Clouds_RF_data; Contrails_RF_data];
        Ntimes = 1;
    elseif icase == 4
        % independent contrails in all-sky
        Ntimes = size(Contrails_RF_data,1);
    end
    
    if ~isempty(RF_data)

    for i = 1:Ntimes
        
    if icase == 4
        RF_data = [Clouds_RF_data; Contrails_RF_data(i,:)];
    end
        
    RF_data_ordbA = sortrows(RF_data,2);
    A_smaller = 0;
    for p=1:size(RF_data,1)
        RF_data_perStep = RF_data_ordbA(p:size(RF_data,1),:);
        Aov_perStep = min(RF_data_perStep(:,2)) - A_smaller;
        A_smaller = min(RF_data_perStep(:,2));
        
        RF_data_ordbH = sortrows(RF_data_perStep,1);
        
        TotTau = 0; TotTauEff = 0; F1 = 1; F3 = 0; OLRn = OLR;
        
        % go through every layer
        GammaAll = 1./(1-RF_data_ordbH(:,5));
        for q=1:size(RF_data_ordbH,1)
            % LW
            Tau = RF_data_ordbH(q,4); T = RF_data_ordbH(q,3); g = RF_data_ordbH(q,5);
            Eps = 1-exp(-delta*Tau);
            O = Sigma*T^k;
            
            OLRn = (1-Eps)*OLRn + Eps*O;
            OLRn = OLRn/1;

            % SW
            TotTau = TotTau + Tau;
            TotTauEff = TotTauEff + Tau/mu;
            
            Gamma = 1/(1-g);
            
            F1 = F1 * Gamma;
            F3 = F3 + Tau*prod(GammaAll(1:end ~= q));
        end
        RF_LWAov = (OLR-OLRn);
        RF_LW = RF_LW + RF_LWAov*Aov_perStep;
        
        F2 = TotTau;
        GammaW = F1*F2/F3;

        R = TotTauEff/(GammaW + TotTauEff);
        R_prime = 2*TotTau/(GammaW + 2*TotTau);
        
        RF_SWAov = -SDR*t*(1-alpha)*(R-alpha*R_prime)/(1-alpha*R_prime); 
        RF_SW = RF_SW + RF_SWAov*Aov_perStep;
                
    end
    end
        
    data(1) = RF_LW; data(2) = RF_SW; data(3) = max(RF_data(:,2));
    end

    % setup results
    if icase == 1
        RFOutput(icase).LW = data(1);
        RFOutput(icase).SW = data(2);
        RFOutput(icase).Area = data(3);
    elseif icase == 2
        RFOutput(icase).LW = data(1);
        RFOutput(icase).SW = data(2);
        RFOutput(icase).Area = data(3);
    elseif icase == 3
        RFOutput(icase).LW = data(1);
        RFOutput(icase).SW = data(2);
        RFOutput(icase).Area = data(3);
    elseif icase == 4
        data(1) = data(1) - Ntimes.*RFOutput(1).LW;
        data(2) = data(2) - Ntimes.*RFOutput(1).SW;
        data(3) = RFOutput(3).Area;
        if data(1) == 0; data = zeros(1,3); end
        RFOutput(icase).LW = data(1);
        RFOutput(icase).SW = data(2);
        RFOutput(icase).Area = data(3);
    end

end

            
end
