function Data = vonMisesFit(X,Ang)

% calculating von mises fit for firing Xs/reaction times etc.
% adpated from the sess von mises function, but just takes in linear and
% cirucluar variable nad does the same thing

% x is a linear variable (firing X, reaction time, etc)
% ang is the ciruclar variable in RADIANS (direction, phase of lfp etc)

%% stuff from sess VonMisesTuningCurve

doBaselineComparison = 0; % i dont think i need this, but we can add it later

options = statset('mlecustom');
options.GradObj = 'on';
options.MaxIter = 2e3;
options.MaxFunEvals = 4e3;


% Construct a negative log likelihood function for each xi
nloglf_alt = @(params,data,cens,freq) nloglfVonMisesPoisson(params,data,cens,freq,Ang);
nloglf_null = @(params,data,cens,freq) nloglfConstantPoisson(params,data,cens,freq,Ang);

% Choose parameter starting points
Astart = median(X);
Bstart = 1;
Brange = [ 0 Inf ];

% This (initial) approach to finding inverted tuning curves is basically useless.
% The Poisson distribution has a long tail, so posDiff always ends up greater than negDiff.
% A single outlier biases the comparison.
posDiff = max(X)-Astart;
negDiff = Astart-min(X);
if negDiff>posDiff, Bstart = -1; Brange = [ -Inf 0 ]; end

% This is a more rigorous (albeit computationally intensive) approach to
% finding inverted tuning curves.
if doBaselineComparison
    % Identify inverted tuning curves by applying constraints from (Wang et al, PNAS 2004)
    % - the unit must be inhibited in response to at least 1 target
    % - the unit must not be excited by any of the remaining targets
    % These criteria are tested by comparing target-specific
    % responses to the baseline firing X.
    %     Target = [ Trials{iTaskComp}.Target ]; this is becasue we have target
    %     numbers instead of angles
    Target = Directions;
    uniqueDirs = unique(Directions);
    p = zeros(1,8);
    for k=1:8
        tInd = find(Target == uniqueDirs(k));
        %         tInd = find(Target==k);
        p(k) = calcTimeSeriesPermutationTest(BaseX',X(tInd)',1e3);
    end
    posDiffTotal = numel(find(p<0.05));
    negDiffTotal = numel(find(p>0.95));
    if negDiffTotal>0 && posDiffTotal==0, Bstart = -1; Brange = [ -Inf 0 ]; end
end

Kstart = 1;
S = X*sin(Ang');
C = X*cos(Ang');
Mustart = atan2(S,C); % trigonometric moment


if Bstart < 0
    if Mustart < 0 Mustart = Mustart+pi; elseif Mustart > 0; Mustart = Mustart-pi; end
end
   
% Bstart
% Mustart
% % Brange
% [phat_alt,pci_alt] = mle(X,'nloglf',nloglf_alt,'options',options,...
%     'lowerbound',[0 Brange(1) 0 -pi],'upperbound',[Inf Brange(2) Inf pi],...
%     'start',[Astart Bstart Kstart Mustart]);
% phat_alt
% Astart = mean(X);
% %Pretty unnecessary since Astart is the MLE of lambda for constant poisson model
% [phat_null,pci_null] = mle(X,'nloglf',nloglf_null,'options',options,...
%     'lowerbound',[0],'upperbound',[Inf],...
%     'start',[Astart]);
% 
[phat_alt,pci_alt] = mle(X,'nloglf',nloglf_alt,'options',options,...
    'lowerbound',[0 Brange(1) 0 -pi],'upperbound',[Inf Brange(2) Inf pi],...
    'start',[Astart Bstart Kstart Mustart]);


Astart = mean(X);
%Pretty unnecessary since Astart is the MLE of lambda for constant poisson model
[phat_null,pci_null] = mle(X,'nloglf',nloglf_null,'options',options,...
    'lowerbound',[0],'upperbound',[Inf],...
    'start',[Astart]);

val_alt = nloglfVonMisesPoisson(phat_alt,X,[],[],Ang);
val_null = nloglfConstantPoisson(phat_null,X,[],[],Ang);
D = 2*(val_null-val_alt);

p = 1 - chi2cdf(real(D),3);

Data.phat_alt = phat_alt;
Data.pci_alt = pci_alt;
Data.phat_null = phat_null;
Data.pci_null = pci_null;
Data.p = p;
Data.X = X;
Data.Ang = Ang;
xi = linspace(-pi,pi,100);
Data.Vm.X = fnVonMises(xi, phat_alt);
Data.Gs.X = repmat(phat_null,1,length(xi));
Data.Fit.Angle = xi;
% call the "Fit" which ever best describes the data.
if p < 0.05
Data.Fit.X = Data.Vm.X;
else
Data.Fit.X = Data.Gs.X;
end
% 
% bins{1}(1,:) = [-pi./8,pi./8];
% bins{2}(1,:) = [pi./8,3*pi./8];
% bins{3}(1,:) = [3*pi./8,5*pi./8];
% bins{4}(1,:) = [5*pi./8,7*pi./8];
% bins{5}(1,:) = [7*pi./8,pi];
% bins{5}(2,:) = [2*pi-pi,2*pi-7*pi./8];
% bins{6}(1,:) = [2*pi-7*pi./8,2*pi-5*pi./8];
% bins{7}(1,:) = [2*pi-5*pi./8,2*pi-3*pi./8];
% bins{8}(1,:) = [2*pi-3*pi./8,2*pi-pi./8];
% x = [0, pi./4, pi./2, 3.*pi./4, pi, -3*pi./4, -pi./2, -pi./4];
% Data.Obs.X = calcTuningCurve(Data.Ang, Data.X, bins);
% Data.Obs.Angle = x;


% figure; plot(Data.Obs.Angle,Data.Obs.X,'or'); hold on;
% plot(Data.Fit.Angle,Data.Fit.X,'-b')

end