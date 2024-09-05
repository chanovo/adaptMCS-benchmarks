function [estPf, nSamples, KLvec, sigmaVec, param, sampTot, gTot, LRTot, nStates] = BiCE_MC(LSF, dim, category, inpDist, Ns0, tarCV1, tarCV2, K)
%% 
%{ 
(1)this code is programmed by chano on 11.05.2022
(2)this code employs adaptive ICE method (Geyer et.al.,2019) with a mixture of independent Categorical model 
(3)input:
   LSF: LSF that maps samples(Ns, dim) to the performance(Ns, 1).
   n: dimension of the problem
   category{d}: a vector of size (numberOfStates_d, 1) that storages the possible state of the d-th dimension
   inpDist{d}: a vector of size (numberOfStates_d, 1) that storages the initial probability that is assigned to each category of the d-th dimension
   Ns: number of samples entering into the optimization problem to get the optimal IS distribution
   targetCV: target CV
   maxLev: maximum intermediate level
   K: number of clusters in the mixture parametic family

(4)output:
   estPf    : estimated failure probability
   compCost : computational cost 
   param    : parameters of the mixture model
%}

% model selection parameters
% K   = 10; % there are 10 mixtures in each component
Rho = 0;  % screening out the dimensions whose IM is smaller equal than Rho

% initialization  
l        = 1;     
nSamples = 0;    
sigma    = 10000;
maxLev   = 30;

nMix     = 1;
alpha    = 1;
theta{1} = inpDist;   

% define the memory
sigmaVec       = sigma;
sampCVofW      = [];
visitedStates  = [];
gVisitedStates = [];
KLvec          = [];

param   = cell(maxLev, 1); param{1,1} = {alpha; theta};
sampTot = cell(maxLev, 1);
gTot    = cell(maxLev, 1);
LRTot   = cell(maxLev, 1);

%% iCE iterations
while true
    %% sampling
    Ns = Ns0;
    while true 
        % Sampling from the mixture model and calculating the correspoding LSFs
        samples = zeros(Ns, dim);
        g       = zeros(Ns, 1);  
        probNom = ones(Ns, 1); 
        probRef = ones(Ns, 1); 
        
        z = nMix+1-sum(rand(Ns, 1) <= repmat(cumsum(alpha'), Ns, 1), 2);
                
        for i = 1:Ns
            
            u  = rand(dim, 1);             
            probRefMat = zeros(nMix, dim);
            
            for d = 1:dim
                
                s      = category{d};
                D_id   = theta{z(i)}{d}; 
                inpD_d = inpDist{d};
                
                ind = find(u(d) <= cumsum(D_id), 1, 'first'); % sampling
                
                samples(i, d) = s(ind);
                probNom(i)    = probNom(i)*inpD_d(ind);  
                
                for k = 1:nMix
                    D_kd = theta{k}{d};
                    probRefMat(k, d) = D_kd(ind);
                end
                               
            end 
            
            probRef(i) = alpha'*prod(probRefMat, 2);
            
            % calculate the performance
            if isempty(visitedStates)
                flag = 0;
            else
                [flag, idx] = ismember(samples(i, :), visitedStates, 'rows');
            end
            
            if flag == 0
                g(i) = LSF(samples(i, :));
                visitedStates  = [visitedStates; samples(i, :)];
                gVisitedStates = [gVisitedStates; g(i)];             
            else
                g(i) = gVisitedStates(idx);
            end

            disp(['lev: ', num2str(l), '...sample: ', num2str(i)]);
        end
        
        g(g<=0)  = 0;
        nSamples = nSamples+Ns;

        if length( unique(g) ) == 1; Ns = 2*Ns; else; break; end        
    
    end
    
    %% calculate the likelhood ratio, LR    
    % Likelihood ratio   
    L = probNom./probRef; % of size (Ns, 1)
    
    %% check if current distribution is close enough to the optimal ISD
    q   = (g<=0)./normcdf(-1./sigma.*g);
    cov = std(q)/mean(q);
    if cov < tarCV1 || sum(isnan(q)) > 0 || l >= maxLev; break; end
    
    %% calculate the estimate of the KL(p_target||p_ref)
    w = normcdf(-1./sigma.*g).*L; % of size (Ns, 1)
    
    w( w == 0 ) = 10^-30;
    ita    = mean(w);
    KL_est = mean(w/ita.*log(w/ita));
    
    %% update the sigma (determine the target distribution)
    % define the objective function to minimize    
    % f = @(x) normcdf(-1./x.*g).*L;    
    f = @(x) normcdf(-1./x.*g)./normcdf(-1./sigma.*g);
    
    objFunc = @(x) abs( std(f(x))/mean(f(x))-tarCV2 );
    
    % get the new sigma
    %{
    grid    = 10^-3:10^-3:sigma; grid = grid';
    gridVal = zeros(length(grid), 1);
    for i = 1:length(grid)
        gridVal(i) = objFunc(grid(i));
    end
    [~, ind] = min(gridVal);
    sigma    = grid(ind);
    
    % plot the objFuncVal with respect to sigma    
    opt_pic = 0; % plotting option 
    if opt_pic
        plot(grid, gridVal, 'linewidth', 1);
        xlim([0,10]); hold on;   
    end
    %}
    
    sigma = fminbnd(objFunc, 0, sigma);
    
    %% update the parameters of the mixture model  
    w = normcdf(-1./sigma.*g).*L; % of size (Ns, 1)    
    
    % EM algorithm for weighted MLE estimator considering model selection
    [alpha, theta, nMix, ~] = EM_MS(samples, w, category, inpDist, K, Rho);
    
    % increase the counter and store the parameters of the reference distribution 
    param{l+1, 1} = {alpha; theta};
    sampTot{l}  = samples;
    gTot{l}     = g;
    LRTot{l}    = L;
    
    sigmaVec    = [sigmaVec; sigma];
    KLvec       = [KLvec; KL_est];
    sampCVofW   = [sampCVofW; std(w)/mean(w)];
    
    l           = l+1;

end 

w = normcdf(-1./sigma.*g).*L; % of size (Ns, 1)    
w( w == 0 ) = 10^-30;

ita = mean(w);                   
KL_est = mean(w/ita.*log(w/ita));
KLvec  = [KLvec; KL_est];        
    
w = (g<=0).*L; % of size (Ns, 1)
w( w == 0 ) = 10^-30;

ita    = mean(w);
KL_est = mean(w/ita.*log(w/ita));

% Calculate the failure probability
estPf = ita;   

% Calculate the computationEffort 
nSamples = nSamples + Ns;

% storage 
sampCVofW  = [sampCVofW; std(w)/mean(w)];
KLvec      = [KLvec; KL_est];
sampTot{l} = samples;
gTot{l}    = g;
LRTot{l}   = L;
nStates    = size(gVisitedStates, 1);

% delete empty entrances
param(cellfun('length', param) == 0)     = [];
sampTot(cellfun('length', sampTot) == 0) = [];
gTot(cellfun('length', gTot) == 0)   = [];
LRTot(cellfun('length', LRTot) == 0) = [];

end

%===========================================================================
%===========================NESTED FUNCTIONS================================
%===========================================================================
%%
function [alpha, theta, nMix, maxLL] = EM_MS(samps, w, category, inpDist, K, Rho)
%{
 this code is used for EM algorithm for updating mixture of independent categorical distributions
 input: 
 output:
%}
% hyperparameters
nPilotRun = 50; % for choosing a good starting point of EM algorithm

% initialization 
Ns  = size(samps, 1); 
dim = size(samps, 2);

% compute the score for each dimension (relative Bayesian IM)
BM = cell(dim, 1);
s  = ones(dim, 1);
%{c
for d = 1:dim
    nsd   = length(category{d});
    BM{d} = zeros(nsd, 1);
    for i = 1:nsd
        BM{d}(i) = w'*double(samps(:, d) == category{d}(i))/sum(w);
    end
    s(d) = BM{d}(1)/inpDist{d}(1);
end
%}

% storage
BIC_mat    = zeros(length(K), length(Rho));
Dim_mat    = zeros(length(K), length(Rho));
maxLL_mat  = zeros(length(K), length(Rho));
param_cell = cell(length(K), length(Rho)); 
% color       = rand(length(Rho), 3);

for i = 1:length(K)
    
    nMix = K(i);
     
    for j = 1:length(Rho) 
        
        % choose the "important" comp.s(dimensions), I, based on their average failure probabilities        
        % t = prctile( s, 100*(1-Rho(j)) );          
        t = Rho(j);        
        I = find( s >= t );
                
        % EM algorithm for fitting the mixture model M(k, rho)
        % pilot run determining the starting point of EM algorithm    
               
        if nMix == 1
            probOfZgivenXk = ones(Ns, 1);
        else          
            maxLL_plt = -inf;
            for ii = 1:nPilotRun
                
                t  = sort(rand(Ns, nMix-1), 2); 
                t1 = [zeros(Ns, 1), t]; t2 = [t, ones(Ns, 1)];
                randomStart = t2 - t1;

                [~,~,LL_ii] = EM(samps, w, category, inpDist, nMix, I, randomStart, 'pilot'); 
                if LL_ii > maxLL_plt
                    maxLL_plt = LL_ii; 
                    probOfZgivenXk = randomStart; 
                end               
            end     
        end
        
        % offcial run
        [alpha_ij, theta_ij, maxLL_ij] = EM(samps, w, category, inpDist, nMix, I, probOfZgivenXk, 'official'); 
        
        % model dimension and BIC
        Dim_ij = nMix-1;  % #free parameters in alpha
        for d = 1:dim
            if ismember(d, I)
                nd = length(category{d});
                Dim_ij = Dim_ij+nMix*(nd-1);
            end
        end
        
        BIC_ij = Dim_ij*log(Ns)-2*maxLL_ij;
        
        % storage
        % scatter(D_ij, maxLL_ij, 50, color(j, :), 'filled'); hold on;
        BIC_mat(i,j)     = BIC_ij;
        Dim_mat(i, j)    = Dim_ij; 
        maxLL_mat(i, j)  = maxLL_ij;
        param_cell{i, j} = {alpha_ij, theta_ij};
            
    end    
end

[ind1, ind2] = find(BIC_mat == min(min(BIC_mat)));

%{
% a new criterion for chosing appropriate models
ind = []; cand = [];
for i = 1:length(K)
    for j = 1:length(Rho)
        t = sum(sum( (D_mat <= D_mat(i,j)).*(maxLL_mat >= maxLL_mat(i,j)) ));
        if t-1 == 0 % exlude (i,j) itself
            ind = [ind;[i, j]]; cand = [cand; [D_mat(i,j), maxLL_mat(i, j)]]; 
        end
    end
end
c     = sortrows([cand, ind]);     
slope = (c(2:end, 2)-c(1:end-1, 2))./(c(2:end, 1)-c(1:end-1, 1)); slope = [log(Ns)/2; slope];
ind1  = c(slope == max(slope), 3);
ind2  = c(slope == max(slope), 4); 
%}

nMix  = K(ind1); rho = Rho(ind2);
maxLL = maxLL_mat(ind1, ind2);
alpha = param_cell{ind1, ind2}{1};
theta = param_cell{ind1, ind2}{2};

end
%%
function [alpha, theta, maxLL] = EM(samps, w, category, inpDist, nMix, I, probOfZgivenXk, opt)

% hyper parameters
switch opt
    case 'pilot'
        maxIt = 20; text = 0;
    case 'official'
        maxIt = 500; text = 1;        
end
a   = 10^-8;
b   = 10;

% initalization
convg = 0;
it    = 1; 
Ns    = size(samps, 1); 
n     = size(samps, 2);
wn    = Ns*w./sum(w); % nomalizing the weight
tol   = 0.1/Ns;

% define the storage      
LP      = zeros(maxIt, 1);
LL      = zeros(maxIt, 1);
paramEM = cell(maxIt, 1); 
gammaEM = cell(maxIt+1, 1); gammaEM{1} = probOfZgivenXk;

while it <= maxIt && ~convg
           
    % (M step) update the parameters of the mixture model
    t     = probOfZgivenXk'*wn+a;
    alpha = t/sum(t); 
    
    theta    = cell(nMix, 1);
    logPrior = a*sum(log(alpha));
    probOfXkGivenZ = ones(Ns, nMix);
    
    for i = 1:nMix
        
        theta_i = cell(n, 1);
        
        for d = 1:n
            
            nd = length(category{d});
            
            for r = 1:nd
                
                ind_r = double(samps(:, d) == category{d}(r)); 
                
                if ismember(d, I)
                    t1 = (wn.*probOfZgivenXk(:, i))' * ind_r;
                    t2 = sum(wn.*probOfZgivenXk(:, i));
                    theta_idr = (t1+b)/(t2+nd*b);   
                    logPrior  = logPrior + b*log(theta_idr);
                else
                    theta_idr = inpDist{d}(r);    
                end
                
                theta_i{d} = [theta_i{d}; theta_idr];
                probOfXkGivenZ(:, i) = probOfXkGivenZ(:, i).*(repmat(theta_idr, Ns, 1).^ind_r);
                
            end                 
        end
        
        theta{i} = theta_i;
    
    end
   
    % (E step) compute the probability of the hiden variable Z given the sample x_k
    % under the current parameters, alpha and theta
    probOfXkAndZ   = probOfXkGivenZ.*alpha';             % of size (Ns, nMix)  
    probOfZgivenXk = probOfXkAndZ./sum(probOfXkAndZ, 2); % of size (Ns, nMix)
    
    % calculate the observed data log posterior
    LP(it) = wn'*log(sum(probOfXkAndZ,2)) + logPrior;
    LL(it) = wn'*log(sum(probOfXkAndZ,2));
    
    % check the convergence of the EM algorithm after the 3rd iteration 
    if it >= 2
        t = LP(it) - LP(it-1); t = abs( t/LP(it) );
        if t <= tol; convg = 1; end            
    end
    
    % update the storage
    paramEM{it}   = {alpha, theta};
    gammaEM{it+1} = probOfZgivenXk;
    
    % increase the counter
    it = it + 1;    
    
end

if text == 1
    if it == maxIt+1
        disp('The EM algorithm does not converge'); 
    else
        disp(['The EM algorithm converges in ', num2str(it-1), ' iterations.']); 
    end
end

maxLL = LL(it-1);

end
