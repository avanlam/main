function [out] = eval_circ(factors)
% EVAL_CIRC : Simulate the four-variable modified Goodwin model for the
% circadian clock, given as input the eighteen constitutive parameters as
% well as the number of oscillators and the chosen model output

%% Parametrisation

    % System and reduced system
    N = factors(1);         % Number of oscillators
    factors(1) = [];
        
    chosen_o = factors(end);% Chosen model output to measure synchronisation

    FinalTime=500;          
    Nsteps=2*FinalTime;
    dt=FinalTime/Nsteps;    % Time step
    Nasympt=Nsteps/2;       % Final time steps (used to compute QoIs)

    Ntrials=1;             

    % Initial point on the stable limit cycle, for each neuron
    x0=repmat([0.0998 0.2468 2.0151 0.0339].',N,1);

    % Fixed parameters (Michaelis-Menten + Light)
     params.nu=[factors(1) factors(2) NaN factors(3) NaN factors(4) NaN factors(5)]';
     params.Ki=[factors(6) factors(7) NaN factors(8) NaN factors(9) NaN factors(10)]';
     params.k=[NaN NaN factors(11) NaN factors(12) NaN factors(13)]';
     params.nuc=factors(14);
     params.Kc=factors(15);
     params.K=factors(16); params.L0=factors(17);
     params.omega=2*pi/factors(18);
    
    rng('default');
    % Intrinsic periods: Uniform distribution
    a=0.8; b=1.2; tau=a+(b-a)*rand(N,Ntrials);
    %tau_all=sort(tau_all,1);

    %% Construct training set

    [t_out,x_out]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 FinalTime], x0, dt);
    
    if N > 1
        % Estimate asymptotic synchronisation order
        [synchrony_full,sync_param_full,avg_gene_conc_full,spectral_ampl_fact_full,estimated_period_full,order_param_X_full, order_param_factor] = compute_QoIs(t_out,x_out,N,Nasympt,params.omega,params.L0);
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