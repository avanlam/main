function [out] = eval_circ_all(factors)
%% Parametrisation

    % System and reduced system
    N = factors(1);         % Number of oscillators
    factors(1) = [];
        
    chosen_o = factors(end);% Chosen type of indicator for the model output

    FinalTime=500;          
    Nsteps=2*FinalTime;
    dt=FinalTime/Nsteps;    % Time step
    Nasympt=Nsteps/2;       % Final time steps (used to compute QoIs)

    Ntrials=1;             

    % Initial point on the stable limit cycle, for each neuron
    x0=repmat([0.0998 0.2468 2.0151 0.0339].',N,1);
    st = 0.05;
    
    % Fixed parameters (Michaelis-Menten + Light)
    params.nu=[normrnd(factors(1), st, N,1) normrnd(factors(2), st, N,1) NaN(N,1) normrnd(factors(3), st, N,1) NaN(N,1) normrnd(factors(4), st, N,1) NaN(N,1) normrnd(factors(5), st, N,1)]';
    params.Ki=[normrnd(factors(6), st, N,1) normrnd(factors(7), st, N,1) NaN(N,1) normrnd(factors(8), st, N,1) NaN(N,1) normrnd(factors(9), st, N,1) NaN(N,1) normrnd(factors(10), st, N,1)]';
    params.k=[NaN(N,1) NaN(N,1) normrnd(factors(11), st, N,1) NaN(N,1) normrnd(factors(12), st, N,1) NaN(N,1) normrnd(factors(13), st, N,1)]';
    params.nuc=normrnd(factors(14), st, N,1);
    params.Kc=normrnd(factors(15), st, N,1);
    params.omega=repelem(2*pi/factors(16), N,1);

    rng('default');
    % Intrinsic periods: Uniform distribution
    a=0.8; b=1.2; tau=a+(b-a)*rand(N,Ntrials);
    %tau_all=sort(tau_all,1);

    % Parameter snapshots
    K=0.6; L0=0.02;

    %% Construct training set

    [t_out,x_out]=RK3 (@(t,x) circadian_rhs(t,x,params,tau,K,L0),[0 FinalTime], x0, dt);
    
    if N > 1
        % Estimate asymptotic synchronisation order
        [synchrony_full,sync_param_full,avg_gene_conc_full,spectral_ampl_fact_full,estimated_period_full,order_param_X_full, order_param_factor] = compute_QoIs(t_out,x_out,N,Nasympt,factors(16),L0);
        possible_outputs = [sync_param_full, spectral_ampl_fact_full, order_param_factor];
        out = possible_outputs(chosen_o);
    else 
        % Estimate asymptotic period and phase
        [~,locs]=findpeaks(x_out(end-Nasympt:end));
        estimated_period = mean(diff(t_out(locs)));
        if isnan(estimated_period)
            plot(x_out(end-Nasympt:end))
            pause
        end
        out = estimated_period;
    end

end