function out=circadian_rhs_colloc(t,x,params,tau,U,meanU_v)
%CIRCADIAN_RHS_COLLOC 
% Author: Niccolo Discacciati

    % Compute reconstruction for deim indices only
    reconstr=U*x;
    X=reconstr(1:4:end-3,:); Y=reconstr(2:4:end-2,:); Z=reconstr(3:4:end-1,:); V=reconstr(4:4:end,:); 
    
    % Coupling term depends on F and can be precomputed
    F=meanU_v*x(3*length(meanU_v)+1:4*length(meanU_v));
    
    % Define collocated output
    out=zeros(size(U,1),1);
    out(1:4:end-3)=params.nu(1)*params.Ki(1)^4./(params.Ki(1)^4+Z.^4)-params.nu(2)*X./(params.Ki(2)+X)+params.nuc*params.K*F./(params.Kc+params.K*F)+params.L0/2*(1+sin(params.omega*t));
    out(2:4:end-2)=params.k(3)*X-params.nu(4)*Y./(params.Ki(4)+Y);
    out(3:4:end-1)=params.k(5)*Y-params.nu(6)*Z./(params.Ki(6)+Z);
    out(4:4:end)=params.k(7)*X-params.nu(8)*V./(params.Ki(8)+V);
    out=out./repelem(tau,4,1);
    
end
