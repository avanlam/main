function training = train_pod(input, system, base_sample, params, tau_all)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

data=zeros(input.n_trials*input.n_samples*(system.n_steps+1),4*input.n_osc);
data_rhs=zeros(input.n_trials*(input.n_samples)*(system.n_steps+1),4*input.n_osc);

for j = 1 : input.n_samples
    params.omega = 2*pi/base_sample(j, 1);
    params.L0 = base_sample(j, 2);
    params.nu(6) = base_sample(j, 3);

    for i=1:input.n_trials
        [t_train,x_train]=RK3 (@(t,x) circadian_rhs(t,x,params,tau_all(:,i)),[0 system.final_time], system.x0, system.dt);
        data((j-1)*(system.n_steps+1)*input.n_trials+(system.n_steps+1)*(i-1)+1:(j-1)*(system.n_steps+1)*input.n_trials+(system.n_steps+1)*i,:)=x_train;

        x_rhs=zeros(size(x_train));
        for t=1:length(t_train)
            x_rhs(t,:)=circadian_rhs(t_train(t),x_train(t,:)',params,tau_all(:,i))';
        end
        data_rhs((j-1)*(system.n_steps+1)*input.n_trials+(system.n_steps+1)*(i-1)+1:(j-1)*(system.n_steps+1)*input.n_trials+(system.n_steps+1)*i,:)=x_rhs;
    end
end

training = {data,data_rhs};

end

