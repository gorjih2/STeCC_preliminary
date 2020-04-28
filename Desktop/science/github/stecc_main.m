%% STeCC reduction in R0
% created by Hossein Gorji April 2020
% based on "STeCC: Smart Testing with Contact Counting Enhances Covid-19 Mitigation by Bluetooth App Based Contact Tracing"
close all
clear all
%% scenario selection
%scenario=1;              % STeCC-A 
%scenario=2;              % STeCC-B
scenario=3;               % STeCC-B with contact tracing

%% genetal parameters
R0=2.4;       % base R0
m=4;          % minimum contacts
gamma=0.3;    % power law k^{-(2+gamma)}
kc=552;       % maximum contact
%% contact tracing
tps=0.5;      % pre-symptomatic time
ts=1;         % symptomatic time
ta=11.5;      % asymptomatic recovery time
tss=10;       % self-quarantine symptomatic recovery time
r1=1/3;       % fraction of asymptomatic
r2=0.5;       % infectiousness of asymptomatic
theta_trace=((1-r1)*(tps+ts)+(1-r1)*tss/10)./(r1*r2*ta+(1-r1)*(tps+ts)+(1-r1)*tss/10);    % R_0^{sym}/R_0
%% network
lN=100;                                                                         % number of discretized points in test and user spaces
alphab=(1+gamma)*m^(1+gamma)./(1-(kc/m)^(-(1+gamma)));                          % pre-factor of the power law distribution
Fcap=@(a) alphab/(1+gamma)*(a^(-(1+gamma))-kc^(-(1+gamma)));                    % fraction of people with contacts above a
meank=-alphab/(gamma).*(kc.^(-(gamma))-m^(-(gamma))); 
vark=alphab/(1-gamma).*(kc.^((1-gamma))-m^((1-gamma)));
%% simulation parameters 
xetavec=linspace(0.5,1,lN);             % percentage of users
Ntests=linspace(10,100,lN)*5;           % number of tests
R0new=zeros(lN,lN);
cycle=7;                                % testing cycle
eta=0.05;                               % false negatives
taup=1;                                 % test processing time
infect=infect_f(taup);                  % probability of virus spread by high-contact people

%% R0 reduction estimate
    for sss=1:lN                            % loop over percentage of users
        for i=1:lN                          % loop over number of tests
            switch logical(scenario)
                case scenario==1
                 efficacy=1-(1-xetavec(sss)+xetavec(sss)*(eta+(1-eta)*infect));                            % efficacy of STeCC-A
                case scenario==2
                 efficacy=1-(1-xetavec(sss)+xetavec(sss)*(eta+(1-eta)*(1-xetavec(sss))*infect));           % efficacy of STeCC-B                    
                case scenario==3
                 efficacy=1-(1-xetavec(sss)+xetavec(sss)*(eta+(1-eta)*(1-xetavec(sss))*infect));           % efficacy of STeCC-B
            end
                  
                 kcn=exp(-1/(gamma+1)*log((gamma+1)*Ntests(i)*cycle/(100000*alphab*(xetavec(sss)))+kc^(-(1+gamma)))); % computing cut-off
                 mean1=-alphab/(gamma).*(kcn.^(-(gamma))-m^(-(gamma)));                                    % weighted average contact of core network
                 mean2=-alphab/(gamma).*(kc.^(-(gamma))-kcn.^(-(gamma)));                                  % weighted average contact of high-contact network
                 factor1=mean1+(1-efficacy)*mean2;                                                         % average contact in new network


                 var1=alphab/(1-gamma).*(kcn.^((1-gamma))-m^((1-gamma)))-mean1;                             % weighted <k^2-k> of core network
                 var2=alphab/(1-gamma).*(kc.^((1-gamma))-kcn.^((1-gamma)))-mean2;                           % weighted <k^2-k> of high-contact network
                 factor2=var1+(1-efficacy)*var2;                                                            % <k^2-k> in new network
 
                switch logical(scenario)
                     case scenario==1
                             factor1=factor1/meank;         
                             factor2=factor2/(vark-meank);
                             R0new(i,sss)=R0*factor2/factor1;                                         % reduction in R0 due to network modification
                    case scenario==2
                            factor1=factor1/meank;         
                             factor2=factor2/(vark-meank);
                             R0new(i,sss)=R0*factor2/factor1;                                         % reduction in R0 due to network modification
                    case scenario==3             
                            R0ct=(1-xetavec(sss)*theta_trace)*R0;                                   % reduction in R0 due to contact tracing
                            factor1=factor1/meank;
                            factor2=factor2/(vark-meank);
                            R0new(i,sss)=R0ct*factor2/factor1;                                        % reduction in R0 due to network modification                       
                end
        end
    end

%% visualization

h1=figure;
set(h1, 'Position', [160 303 606 459]);
pcolor(Ntests,xetavec*100,R0new')
h=colorbar; 
set(get(h,'label'),'string','achievable R_0 value','FontSize',15.4,'FontName', 'Helvetica Neue','Color','black');

caxis([0.5 2]);
shading interp;

hold on 
C=contour(Ntests,xetavec*100 , R0new', [1 1],'-k','LineWidth',2);

clabel(C,'FontSize',15,'FontName', 'Helvetica Neue','Color','black')


C=contour(Ntests,xetavec*100 , R0new',[1.7 1.7],'--k','LineWidth',1);
clabel(C,'FontSize',15, 'FontName','Helvetica Neue','Color','black')

 C=contour(Ntests,xetavec*100 , R0new',[0.8 0.8],'--k','LineWidth',1);
 clabel(C,'FontSize',15,'FontName', 'Helvetica Neue','Color','black')
 
 C=contour(Ntests,xetavec*100 , R0new',[0.6 0.6],'-.k','LineWidth',1);
 clabel(C,'FontSize',15,'FontName', 'Helvetica Neue','Color','black')
 axis square
 hold on

set(gca,'FontName','Helvetica Neue')
set(gca,'FontSize',15)

HH=sprintf('percentage of app users  \n among smart-phone owners');
 ylabel(HH)
 xlabel('number of tests per 100''000 per day')
 
 hold on 
