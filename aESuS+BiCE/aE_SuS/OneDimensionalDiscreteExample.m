clear, clc
dbstop if error 

% define the parameters
func_choice      = 3; 
appending_method = 'EC'; % 'EC' or 'RA'
scheme           = 'adaptive'; % 'geometric', 'arithmetic' or 'adaptive'
dim = 1; p0  = 0.1; N0 = 1000; tol = 0.5;

% test functions£¨one-dimensional problem
switch func_choice
    case 1
        % continuous function
        LSF = @(x) x+5;
        Pf_exactly = 3.167*10^-5;
    case 2
        % discrete function without big jump
        CDF_vec = [10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 0.5, 1];
        state = [-6, -4, -3, -2, -1, 0, 1];
        LSF = @(x) state(find((x-norminv(CDF_vec))<=0, 1))+5;
        Pf_exactly = 10^-5;
    case 3
        % discrete function with big jump
        CDF_vec = [10^-5, 3*10^-5, 5*10^-5, 3*10^-2, 10^-1, 0.5, 1];
        state = [-6, -4, -3, -2, -1, 0, 1];
        LSF = @(x) state(find((x-norminv(CDF_vec))<=0, 1))+5;
        Pf_exactly = 10^-5;    
end

% one-run of Subset Simulation(different appending methods)
%{
switch appending_method
    case 'RA'
        [Pf, comp_cost, b_tot, condiProb_tot, samples_tot, f_samples_tot, seeds_tot, f_seeds_tot, lamda_tot] = aE_SuS_RA(LSF, dim, N0, tol, p0, scheme);
    case 'EC'
        [Pf, comp_cost, b_tot, condiProb_tot, samples_tot, f_samples_tot, seeds_tot, f_seeds_tot, lamda_tot] = aE_SuS_EC(LSF, dim, N0, tol, p0, scheme);
end
%}

% multi-runs of Subset Simulation(different appending methods)
%{c
num_run = 1000;
est_Pf = zeros(num_run,1); 
num_eval = zeros(num_run,1);
for i = 1:num_run
    switch appending_method
        case 'RA'
            [est_Pf(i), num_eval(i)] = aE_SuS_RA(LSF, dim, N0, tol, p0, scheme);
        case 'EC'
            [est_Pf(i), num_eval(i)] = aE_SuS_EC(LSF, dim, N0, tol, p0, scheme);
    end
    disp(i); 
end 
relative_bias = (mean(est_Pf)-Pf_exactly)/Pf_exactly;
cov_est_Pf = std(est_Pf)/mean(est_Pf);
mean_num_eval = mean(num_eval);
%}