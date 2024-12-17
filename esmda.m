
function [xall,yall] = esmda(xf,yf,range,obs,sd,N_Iter)
currentdir = pwd;
Nobs = length(obs);     % number of the measurements
Npar = size(xf,1);      % number of the parameters
Ne = size(xf,2);        % ensemble size

Cd = eye(Nobs);         
for i = 1:Nobs
    Cd(i,i) = sd(i)^2;  % covariance of the measurement errors   
end
factor = ones(N_Iter,1)*sqrt(N_Iter);       % inflate Cd in the multiple data assimilation scheme

xall = xf; yall = yf;                       % store results at each iteration
xa = nan(size(xf)); ya = nan(size(yf));     % define the updated ensemble

for n_i = 1:N_Iter
    beta = factor(n_i);
    xa = es_1(xf,yf,range,sd*beta,obs);
    save xa.mat xa;
    for k = 1:Ne 
        ya(:,k) = model_H(xa(:,k)); 
    end
    % Enhance the performance of ILUES in problems with many parameters 
    if Npar > 10    
        likf = Cal_Log_Lik(yf,obs,sd);   
        lika = Cal_Log_Lik(ya,obs,sd);
        cc = (exp(lika - likf)) < rand(Ne,1);
        xa(:,cc) = xf(:,cc);
        ya(:,cc) = yf(:,cc);
    end
    
    xall = [xall xa]; 
    yall = [yall ya];
    xf = xa; yf = ya;
    save(['xa',num2str(n_i),'.mat'],'xa');
    save(['ya',num2str(n_i),'.mat'],'ya');    
    save ESMDA.mat;
    
end

delete ESMDA.mat;

end