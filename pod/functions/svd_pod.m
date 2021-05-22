function [U, S, U_pod, U_old] = svd_pod(data,num_n, n_latent_pod)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

for i=1:4
    % SVD for X,Y,Z,V
    [U_old{i},S{i},~]=svd(data(:,i:4:end-4+i).','econ'); 
    
    U{i}=U_old{i}(:,1:n_latent_pod); S{i}=diag(S{i});
end

% The matrix is blockdiagonal, but I have to reorder the indexes 
ids=[1:4:4*num_n-3, 2:4:4*num_n-2, 3:4:4*num_n-1, 4:4:4*num_n];
U_pod(ids,:)=(blkdiag(U{:}));

end

