%% computing probability of virus spread by high-contact people
function probA=infect_f(taup)
mu0=5.34;                               % mean of latency time
sigma0=2.7249^2;                        % variance of latency time
mu=log(mu0^2/sqrt(sigma0+mu0^2));       % finding correponding mean in log-normal
sigma=sqrt(log(sigma0/mu0^2+1));        % finding correponding std in log-normal
Nsample=1000000;                         % number of samples
s=0;                                    % event of virus spread    
Tg=7+taup;                              % cycle+latency
for i=1:Nsample
        taui=rand*(Tg-taup);            % infection time
        ns=lognrnd(mu,sigma,1,1);
        deltat=Tg-(taui+ns);
       if deltat>0
            s=s+1;                     % event of virus spread                
        end        
end
probA=s/Nsample;
end
