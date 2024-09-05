function [Pf, comp_cost, b_tot, condiProb_tot, samples_tot, f_samples_tot, seeds_tot, f_seeds_tot, lambda_tot] = aE_SuS_RA(LSF, dim, N0, tol, p0, scheme)
%{
This code appends new samples in conditional levels randomly from the exsist samples, named aE-SuS-RA.
%}
% check the input 
if mod(N0*p0, 1) ~= 0; error('N0*p0 should be an integral'); end
% if tol >= 1 || tol <= 0; error('tol should be in (0, 1)'); end

% initialization
max_l      = 50;
lambda     = 0.6;
l          = 0;  
b          = [0, inf]; 
comp_cost  = 0; 
% define the memory
samples_tot        = cell(max_l, 1); 
f_samples_tot      = cell(max_l, 1);
seeds_tot          = cell(max_l, 1); 
f_seeds_tot        = cell(max_l, 1);
b_tot              = zeros(max_l, 2); 
lambda_tot         = zeros(max_l, 1); 
condiProb_tot      = zeros(max_l, 1);

while b(1, 2) > 0 && l < max_l
    
    % initialization of the l-th level
    iter = 1; Ns = 0; N = N0; 
    samples = []; f_samples = [];
    
    while Ns < tol*N0*p0        
        % Generate the samples
        if l == 0
            samples_iter = zeros(N, dim); f_samples_iter = zeros(N, 1); 
            rng('shuffle');
            for i = 1:N 
                samples_iter(i, :) = normrnd(zeros(1, dim), ones(1, dim)); 
                f_samples_iter(i) = LSF(samples_iter(i, :)); 
                disp(['No:', num2str(i), '/', num2str(N), ' Level:', num2str(l), ' Iteration:',num2str(iter)]);
            end 
            % append new samples
            samples = [samples; samples_iter]; f_samples = [f_samples; f_samples_iter];
        else 
            % append new samples to the seeds
            [samples, f_samples, lambda] = my_aCS(seeds, f_seeds, LSF, b, N, lambda, adapt, l, iter);                      
        end
                       
        % sort the samples by their LSF-s
        temp = [samples, f_samples]; temp = sortrows(temp, dim+1);
        samples_sorted = temp(:, 1:dim); f_sorted = temp(:, dim+1);
        N_cur = size(samples_sorted, 1); % current sample size
        
        % determine the samples inside the failure domain 
        if iter == 1
            b_temp = f_sorted(p0*N0);
            if b_temp <= 0
                b_temp = 0; 
                Ns = sum(f_sorted <= 0); sign = 1;
                break;
            elseif b_temp ~= f_sorted(end)
                Ns = sum(f_sorted <= b_temp); sign = 1;
                break;
            end
        end 
        
        Ns = sum(f_sorted < b_temp); sign = 0;
        
        % update the parameters if necessary
        if Ns < tol*N0*p0

            switch scheme
                % N_new is the sample size for the next iteration.
                case 'geometric'
                    N_nex = 2*N_cur; 
                case 'arithmetic'
                    N_nex = N_cur+N0;
                case 'adaptive'
                    N_nex = ceil(tol*N0*p0/max(1, Ns)*N_cur);
            end
            
            if l == 0 % MCS level
                N = N_nex-N_cur;  
            else % conditional level
                N = N_nex;  
                adapt = 0; % whether adapt the lambda in the following iteration
                seeds = samples_sorted; f_seeds = f_sorted; % reset the seeds and the corresponding value
            end
            
            iter = iter+1;
            
        end

    end
    
    % initialization for the (l+1) level(next conditional level)  
    b     = [sign, b_temp]; 
    adapt = 1;
    seeds = samples_sorted(1:Ns, :); f_seeds = f_sorted(1:Ns);
    
    % calculate the total comp. cost and the conditional probability
    comp_cost = comp_cost + N_cur;
    p_l       = Ns/N_cur;
    
    % update the memory
    samples_tot(l+1)        = {samples_sorted}; 
    f_samples_tot(l+1)      = {f_sorted};
    seeds_tot(l+1)          = {seeds}; 
    f_seeds_tot(l+1)        = {f_seeds};
    b_tot(l+1, :)           = b;
    condiProb_tot(l+1)      = p_l;
    lambda_tot(l+1)         = lambda;
    
    % increase the counter
    l = l+1;
end

% delete the empty elements
samples_tot(cellfun('length', samples_tot) == 0) = [];
f_samples_tot(cellfun('length', f_samples_tot) == 0) = [];
seeds_tot(cellfun('length', seeds_tot) == 0) = [];
f_seeds_tot(cellfun('length', f_seeds_tot) == 0) = [];
b_tot(all(b_tot == 0, 2), :) = []; 
condiProb_tot(condiProb_tot == 0) = [];
lambda_tot(lambda_tot == 0) = [];

% calculate the failure probability
Pf = prod(condiProb_tot);end

%% define the aCS function
function [samples, f_samples, lambda_updated] = my_aCS(seeds, f_seeds, LSF, b, N, lambda, adapt, lev, it)
   
    % check the input
    Ns = size(seeds, 1); 
    dim = size(seeds, 2);
    if N < Ns; error('the sample size should be larger than the number of seeds'); end
    % if there is no adaptation, then permutating the seeds and adaptively updating the lambda are unnecessary
    if adapt == 0
        permut = 0; update = 0; 
    else
        permut = 1; update = 1; 
    end
    
    % randomly permutate the seeds
    if permut == 0
    else
        temp = [seeds, f_seeds]; temp = temp(randperm(Ns), :);
        seeds = temp(:, 1:dim); f_seeds = temp(:, dim+1);
    end
    
    % determine the number of groups, in [1, Ns];
    Ng = max(1, min(Ns, ceil((N-Ns)/100)));    % update the lambda for approximately every 100 new samples
    
    % determine the length of each group
    % len_group =  floor(Ns/Ng)*ones(Ng, 1); 
    % index = randperm(Ng, mod(Ns, Ng)); len_group(index) =  len_group(index)+1;
    
    % determine the length of each chain
    len_chain = floor(N/Ns)*ones(Ns, 1); 
    index = randperm(Ns, mod(N, Ns)); len_chain(index) = len_chain(index)+1;
    
    % initialization
    optim_accepted_rate = 0.44;
    choice = 'b';
    switch choice
        case 'a'; sigma_init = std(seeds, 1); % (1, dim)
        case 'b'; sigma_init = ones(1, dim);         % (1, dim)
    end
    
    samples = zeros(N, dim);
    f_samples = zeros(N, 1);
    accepted_rate_tot = zeros(Ng, 1);
    lambda_tot = zeros(Ng, 1);
    
    accepted_times = 0;        % counts the accepted times in each group
    sample_times = 0;          % counts the conditional sampling times in each group
    invokeLSF_times = 0;       % counts the times of invoking the LSF
    iter = 1;                  % counts the No. of iteration(group)
    counter = 1;               % counts the No. of sample
    
    
    for i = 1:Ns
        
        % set the i-th seeds as the first sample of the i-th chain
        samples(counter, :) = seeds(i, :); f_samples(counter) = f_seeds(i); counter = counter+1;
        
        % generate the rest (len_chain-1) samples using CS algorithm
        sigma = min(1, lambda*sigma_init);
        rou = sqrt(1-sigma.^2);  % calculate the cross-correlation parameter, rou, which is a row vector
        
        if len_chain(i) == 1            
        else
            for j = 1:(len_chain(i)-1)
  
                candidate = normrnd(samples(counter-1, :)*diag(rou), sqrt(1-rou.^2)); % componentwise sampling scheme
                
                f_candidate = LSF(candidate); invokeLSF_times = invokeLSF_times+1;
                
                if b(1, 1) == 1; condition = f_candidate <= b(1, 2); else; condition = f_candidate < b(1, 2); end
                
                if condition
                    samples(counter, :) = candidate; f_samples(counter) = f_candidate;
                    accepted_times = accepted_times+1;
                else
                    samples(counter, :) = samples(counter-1, :); f_samples(counter) = f_samples(counter-1);
                end
                
                sample_times = sample_times+1;
                counter = counter+1; 
                disp(['No: ', num2str(invokeLSF_times), '/', num2str(N-Ns), ' Level:', num2str(lev), ' Iteration:',num2str(it)]);
            
            end          
        end
        
        % update the lambda when entering into a new group if adapt ~= 0
        
        if update == 0
        else       
            if mod(i, floor(Ns/Ng)) == 0
            % if any(i == cumsum(len_group))
                % calculate the estimate of the accepted rate for the last group
                if sample_times == 0
                else
                    estimate_accepted_rate = accepted_times/sample_times;
                    accepted_rate_tot(iter) = estimate_accepted_rate;
                    lambda_tot(iter) = lambda;
                    accepted_times = 0; sample_times = 0;
                    % update the lambda and counter
                    lambda = lambda*exp(iter^(-1/2)*(estimate_accepted_rate-optim_accepted_rate));
                    iter = iter+1;
                end               
            end
        end
        
    end
      
lambda_updated = lambda;    

end 