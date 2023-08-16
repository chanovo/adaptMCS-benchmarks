function mpc = add_branchCapacity(mpc, alpha)

if all( mpc.branch(:, 6) == 0 )
    % run powerflow and set the edge capacity twice the power flow 
    mpopt = mpoption('out.all', 0);

    results_dcpf = rundcpf(mpc, mpopt);

    mpc.branch(:, 6)   = alpha*0.5*abs( results_dcpf.branch(:, 14) - results_dcpf.branch(:, 16) ); 
end

end