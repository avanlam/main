function [sync_param_full] = eval_circ_phys(N, factors)
%% Parametrisation

    % System and reduced system
    FinalTime=500;          
    Nsteps=2*FinalTime;
    dt=FinalTime/Nsteps;    % Time step
    Nasympt=Nsteps/2;       % Final time steps (used to compute QoIs)

    Ntrials=1;             

    % Initial point on the stable limit cycle, for each neuron
    x0=repmat([0.0998 0.2468 2.0151 0.0339].',N,1);

    % Fixed parameters (Michaelis-Menten + Light)
    params.nu=[factors(1) factors(1) NaN factors(1) NaN factors(1) NaN factors(1)]';
    params.Ki=[factors(2) factors(2) NaN factors(2) NaN factors(2) NaN factors(2)]';
    params.k=[NaN NaN factors(3) NaN factors(3) NaN factors(3)]';
    params.nuc=factors(4);
    params.Kc=factors(5);
    params.omega=2*pi/factors(6);

    rng('default');
    % Intrinsic periods: Uniform distribution
    a=0.8; b=1.2; tau=a+(b-a)*rand(N,Ntrials);
    %tau_all=sort(tau_all,1);

    % Parameter snapshots
    K=0.6; L0=0.02;

    %% Construct training set

    [t_out,x_out]=RK3 (@(t,x) circadian_rhs(t,x,params,tau,K,L0),[0 FinalTime], x0, dt);
    
    % Estimate asymptotic synchronisation order
    [synchrony_full,sync_param_full,avg_gene_conc_full,spectral_ampl_fact_full,estimated_period_full,order_param_X_full] = compute_QoIs(t_out,x_out,N,Nasympt,params.omega,L0);

end