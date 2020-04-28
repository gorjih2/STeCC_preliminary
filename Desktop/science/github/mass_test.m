%% Mass Testing
% created by Hossein Gorji April 2020
% based on "STeCC: Smart Testing with Contact Counting Enhances Covid-19 Mitigation by Bluetooth App Based Contact Tracing"
clear all
close all
%% latency time
mu0=5.34;                               % mean of latency time
sigma0=2.7249^2;                        % variance of latency time
mu=log(mu0^2/sqrt(sigma0+mu0^2));       % corresponding log-normal mean of latency time
sigma=sqrt(log(sigma0/mu0^2+1));        % corresponding log-normal std of latency time  
Nsample=100000;                         % number of samples
ns=lognrnd(mu,sigma,1,Nsample);


%% Testing 
lN=200;                                 % discretization parameter
Nlin=linspace(0.001,0.5,lN);             % testing frequency
N=1./Nlin;
tauvec=[2 1.5 1 0.5];                   % test processing time
Nvec=tauvec;
eta=0.05;                               % ratio of false negatives
%% loop over test processing times
for jjj=1:4
taup=tauvec(jjj);

%% computing detection rate from exposed compartment

TimeL=N;
for j=1:lN
    NN=N(j);
    xs=1+(rand(1,Nsample))*(NN-1);
    LL=0;
    jj=0;
    for i=1:Nsample
        if ((NN-xs(i)+taup)<ns(i))
           LL=LL+1/(NN-xs(i)+taup);
           jj=jj+1;
        end
    end
    TimeL(j)=LL/Nsample;
end
ke=TimeL;

%% computing detection rate from asymptomatic compartment     

TimeL=N;
for j=1:lN
    NN=N(j);
    xs=1+(rand(1,Nsample))*(NN-1);
    LL=0;
    jj=0;
    for i=1:Nsample
         if ((NN-xs(i)+taup)>ns(i))
            if ((NN-xs(i)+taup)<(ns(i)+3/2))
                if rand<1/3
                    LL=LL+1/(NN-xs(i)+taup);
                       jj=jj+1;
                end
            else
         
                if ((NN-xs(i)+taup)<(ns(i)+11.5))
                         LL=LL+1/(NN-xs(i)+taup);
                         jj=jj+1;                 
             
                end  
            end
          end
        
        
    end
    TimeL(j)=LL/Nsample;
end
ka=TimeL;

%% computing detection rate from symptomatic compartment     


TimeL=N;
for j=1:lN
    NN=N(j);
    xs=1+(rand(1,Nsample))*(NN-1);
    LL=0;
    jj=0;
    for i=1:Nsample
          if ((NN-xs(i)+taup)>ns(i))
            if ((NN-xs(i)+taup)<(ns(i)+3/2))
                if rand>1/3
                    LL=LL+1/(NN-xs(i)+taup);
                       jj=jj+1;
                end            
            end
          end
        
        
    end
    TimeL(j)=LL/Nsample;
end


ks=TimeL;

%% Computing R_0

alpha=0.6711;                   % infection rate of symptomatic
epsilon=0.1;                    % infection rate of self-isolated
betaa=0.078;                    % rate of asymptomatic
betas=0.1573;                   % rate of symptomatic
gammaa=0.087;                   % recovery rate of asymptomatic
xims=0.6667;                    % recovery rate of mild symptomatic
gammams=0.08;                   % hospitalization rate of mild symptomatic
xiss=0.02;                      % rate coefficient for successively stronger symptoms
gammass=0.0727;                 % recovery rate of hospitalized




ks=(1-eta)*ks;
ka=(1-eta)*ka;
ke=(1-eta)*ke;

R01=alpha*betas./(betaa+betas).*(betaa./(2*gammaa*betas)+1./xims+epsilon./(gammams+xiss));                              % R0 without testing
R0=alpha*betas./(betaa+betas+ke).*(betaa./(2.*(gammaa+ka).*betas)+1./(xims+ks).*(1+epsilon*xims./(gammams+xiss)));      % R0 with testing


Nvec(jjj)=interp1(R0/R01,N,1/R01);                                                                                       % R0^{wt}=1                          
Nvec(jjj)=1./Nvec(jjj)*100000;                                                                                           % number of tests per 100'000 per day

%% visualization

switch jjj
    case 1
                plot(1./N*100000,R0/R01,'Color',[0.75, 0, 0.75],'LineWidth',2)   
    case 2
                plot(1./N*100000,R0/R01,'-b','LineWidth',2)
    case 3
                plot(1./N*100000,R0/R01,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
    case 4
                plot(1./N*100000,R0/R01,'-k','LineWidth',2)

end
hold on
set(gca,'FontSize',15)

set(gca,'FontSize',15)
set(gca,'FontName','Helvetica Neue')
hold on
xlabel('number of tests per 100''000 per day')
ylabel('R^{wt}_0/R_0')
end

AX=legend('2','1.5','1 (realistic for next gen. seq.)','0.5 (realistic for current tests)');
AX.FontSize = 15;
AX.Title.String = 'Test processing time (days)';
AX.Title.FontSize = 12;
hold on 
h2=plot(linspace(Nvec(4),4500,100),ones(1,100)/2.4,'Color',[0, 0.5, 0],'LineWidth',1);
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


