function [xall,yall] = ilues(xf,yf,range,obs,sd,alpha,N_Iter)
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

meanxf = repmat(mean(xf,2),1,Ne);           % mean of the prior parameters
Cm = (xf - meanxf)*(xf - meanxf)'/(Ne - 1); % auto-covariance of the prior parameters

for n_i = 1:N_Iter
      
    J1 = nan(Ne,1);
    for i = 1:Ne
        J1(i,1) = (yf(:,i)-obs)'/Cd*(yf(:,i)-obs);
    end
    
    beta = factor(n_i);
    parfor j = 1:Ne
        xa(:,j) = local_update(xf,yf,Cm,sd,range,obs,alpha,beta,J1,j);
    end

    
	parfor k = 1:Ne
        ya(:,k) = model_H(xa(:,k),k);
    end
    
    % Enhance the performance of ILUES in problems with many parameters 
    if Npar > 10    
        likf = Cal_Log_Lik(yf,obs,sd);   
        lika = Cal_Log_Lik(ya,obs,sd);
        cc = (exp(lika - likf)) < rand(Ne,1);
        xa(:,cc) = xf(:,cc);
        ya(:,cc) = yf(:,cc);
    end
    
    xall = [xall xa]; %#ok<*AGROW>
    yall = [yall ya];
    xf = xa; yf = ya;
        
    save ILUES.mat
    
end

delete ILUES.mat

end