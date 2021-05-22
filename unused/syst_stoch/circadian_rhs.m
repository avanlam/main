function [out] = circadian_rhs(t,x,params,tau,K,L0)

    %extract variables
    X=x(1:4:end-3,:); Y=x(2:4:end-2,:); Z=x(3:4:end-1,:); V=x(4:4:end,:); F=mean(V,1);
    
    out=zeros(size(x));
    % x_dot
    if numel(tau) > 1
        out(1:4:end-3)=params.nu(1).*params.Ki(1).^4./(params.Ki(1).^4+Z.^4)-params.nu(2).*X./(params.Ki(2)+X)+params.nuc.*K.*F./(params.Kc+K.*F)+L0/2.*(1+sin(params.omega.*t));
    else
        out(1:4:end-3)=params.nu(1)*params.Ki(1)^4./(params.Ki(1)^4+Z.^4)-params.nu(2)*X./(params.Ki(2)+X)+L0/2*(1+sin(params.omega*t));
    end
    % y_dot
    out(2:4:end-2)=params.k(3).*X-params.nu(4).*Y./(params.Ki(4)+Y);
    % z_dot
    out(3:4:end-1)=params.k(5).*Y-params.nu(6).*Z./(params.Ki(6)+Z);
    % v_dot
    out(4:4:end)=params.k(7).*X-params.nu(8).*V./(params.Ki(8)+V);
    
    % divide rhs by frequency of each oscillator (same quantity for all variables)
    out=out./repelem(tau,4,1);

end