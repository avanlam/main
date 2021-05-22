function [synchrony,sync_param,avg_gene_conc,spectral_ampl_fact,estimated_period,order_param_X, order_param_factor] = compute_QoIs(t,x,N,Nasympt,omega,L0)

    %F(t)^2/expectation(V_i(t)^2)
    V=x(:,4:4:end);
    synchrony=mean(V,2).^2./mean(V.^2,2);
    
    %sqrt of time-averaged synchrony variable
    sync_param=sqrt(mean(synchrony(end-Nasympt:end)));
    
    %expectation(X_i(t))
    X=x(:,1:4:end-3);
    avg_gene_conc=mean(X,2);
    %time avg of exp(-i*omega*t)*avg_gene_conc, absolute value squared and normalized by 4/L0^2
    spectral_ampl_fact=4/L0^2*abs(mean(exp(-1i*omega*t(end-Nasympt:end)).*avg_gene_conc(end-Nasympt:end)))^2;
    
    %find period of the avg_gene_conc variable
    [peaks,locs]=findpeaks(avg_gene_conc(end-Nasympt:end));
    estimated_period=mean(diff(t(locs)));
    
    %find phase for each oscillator, using the variable X_i
    for i=1:N
        %X_i
        X=x(:,1:4:end-3);
        %estimate period
        [peaks_i,locs_i]=findpeaks(X(end-Nasympt:end,i));
        %phase grows linearly in the period, and is equal to 2*k*pi at peaks
        if numel(locs_i)<2
            estimated_phase(:,i) = zeros(size(t(end-Nasympt:end)));
        else
            estimated_phase(:,i)=interp1(t(end-Nasympt+locs_i),2*pi*(0:length(locs_i)-1)',t(end-Nasympt:end),'linear','extrap');
        end
    end
    %kuramoto order parameter for the X_i variables
    order_param_X = abs(mean(exp(1i*estimated_phase),2));
    order_param_factor = mean(order_param_X);
end