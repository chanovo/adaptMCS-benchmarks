clear, clc, close all
dbstop if error
%% aE-SuS for different dimensions

% choose the parameters
dim_series = 300;
counter = 1;
appending_method = 'EC'; % 'EC' or 'RA'
scheme           = 'adaptive'; % 'geometric', 'arithmetic' or 'adaptive'
p0 = 0.1; N0 = 1000; tol = 0.5; 

for dim = dim_series 
    
    % define the parameters
    num_item = dim/2; p_ber = 10^-3; mean_norm = -10; std_norm = 1; target_Pf = 10^-3;
    
    % calculate the threshold for current dimension
    func = @(x) FailureProb(x, num_item, mean_norm, std_norm, p_ber)-target_Pf;
    options = optimoptions(@fsolve,'OptimalityTolerance', 10^-16);
    threshold = fsolve(func, mean_norm, options);
    
    % define the LSF
    LSF = @(x) comb_disRv_conRv(x, threshold, mean_norm, std_norm, p_ber);

    % run aE-SuS
    num_run = 500; 
    est_Pf = zeros(num_run,1); comp_cost = zeros(num_run,1);
    
    for i = 1:num_run
          
        switch appending_method
        case 'RA'
            [est_Pf(i), comp_cost(i)] = aE_SuS_RA(LSF, dim, N0, tol, p0);
        case 'EC'
            [est_Pf(i), comp_cost(i)] = aE_SuS_EC(LSF, dim, N0, tol, p0);
        end     
        
        % storage
        if i/20 == fix(i/20)
            mean_est_Pf(counter) = mean(est_Pf(1:i));
            std_est_Pf(counter) = std(est_Pf(1:i));
            mean_comp_cost(counter) = mean(comp_cost(1:i));
        else
        end 
        
        disp([counter, i]);
    end
    
    counter = counter+1;
    save 'temp.mat';
end 


%% subfunctions
function [y] = comb_disRv_conRv(x, threshold, mean_norm, std_norm, p_ber)
% x is a matrix with size(m, n), m is the number of samples while n represents the dimension. Each dimension obeys independent standard normal distribution.
% y is the output with size(m, 1) which is the sum of a series of product of Bernoulli distribution~B(p_ber) and normal distribution~N(mean,var) 
num_sam = size(x, 1); 
dim     = size(x, 2);
if dim/2 ~= fix(dim/2)
    disp('The dimension must be an even');
end

y = -1*threshold*ones(num_sam, 1);
for j = 1:dim/2
    y = y+(std_norm.*x(:, j)+mean_norm).*( x(:, dim/2+j)>norminv(1-p_ber));
end

end

function Pf = FailureProb(threshold, nitem, mean_norm, std_norm, p_ber)
Pf = 0;

for i = 1:nitem
    Pf = Pf+binopdf(i, nitem, p_ber)*normcdf((threshold-i*mean_norm)/sqrt(i*std_norm^2), 0, 1);
end

Pf = Pf + (threshold>=0)*binopdf(0, nitem, p_ber);

end