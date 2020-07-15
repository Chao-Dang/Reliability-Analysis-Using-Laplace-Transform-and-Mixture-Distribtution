% Author: Chao Dang; E-mail: chaodang@outlook.com    Date: 15/07/2020 
% Dang C., Xu J. Unified reliability assessment for problems with low- to high-dimensional random inputs 
% using the Laplace transform and a mixture distribution. Reliability Engineering & System Safety (2020) https://doi.org/10.1016/j.ress.2020.107124
%% To show how the proposed method works, in the following we just take the Exmaple 3 in the paper as an illustration 

clearvars; close all; clc

%% Define the inputs of your model
% dimension of random variables
Num_rv = 8;
% distribution types of the random variables 
Dis_type={'Lognormal','Lognormal','Lognormal','Lognormal','Lognormal','Lognormal','Lognormal','Lognormal'};
% mean value of the random variables 
Mean = [1.5 0.01 1 0.01 0.05 0.02 15 100];  
% cofficient of variation of the random variables
Cov = [0.1 0.1 0.2 0.2 0.4 0.5 0.1 0.1];
% standard deviation of the random variables
Std = Cov.*Mean;

% the parameters of the random variables
for i = 1:Num_rv
        switch Dis_type{i}
        case 'Normal' 
             Para(1,i) = Mean(i);
             Para(2,i) = Std(i);
        case 'Lognormal'
             Para(1,i) = log((Mean(i)^2)/sqrt(Std(i)^2+Mean(i)^2));
             Para(2,i) = sqrt(log(Std(i)^2/(Mean(i)^2)+1));            
        otherwise
            disp('Please add more distributions!')
    end     
end

%% Define your model (Limit state function, or performance function)
per_fun = @TDOF_dynamic_syetem;

%% Monte Carlo simulation
% sample size
Smaple_size = 1e7;

% simple random sampling
for i = 1:Num_rv
    X(:,i) = lognrnd(Para(1,i),Para(2,i),Smaple_size,1);
end   

% evaluate the limit state function 
G = per_fun(X);

% failure probability
Pf_mcs = sum(G<=0)./Smaple_size;

%% Proposed method
% generate the integration points and weights by the cubature formaula II in standard normal space, 
% if 4< Num_rv <30
[ weights , points ] = Five_degree_cubature_II( Num_rv );

% if 1=< Num_rv <=3
% [ weights , points ] = Five_degree_cubature_I( n );

% if Num_rv >= 30，please first download the LSS.m and LPSS.m at: https://ww2.mathworks.cn/matlabcentral/profile/authors/6129139-michael-shields
% num = 31^2;

% Num_rv is odd
% lpss_design = [2*ones((Num_rv-1)/2,1);1];
% lpss_strata = [sqrt(num)*ones((Num_rv-1)/2,1);num] ;

% Num_rv is even
% lpss_design = [2*ones((Num_rv)/2,1)];
% lpss_strata = [sqrt(num)*ones((Num_rv)/2,1)] ;
% 
% x_lpss = LPSS(lpss_design,lpss_strata);
% 
% weights = ones(1,num)./num;
% points = norminv(x_lpss);


% transform the integration points into the original random-variate space
% by e.g. isoprobabilitic transformation 
for j = 1:Num_rv
    XX(:,j) = logninv(normcdf(points(:,j)),Para(1,j),Para(2,j));
end


% evaluate the performance function
GX = per_fun(XX);

% specify the discretization
g = min(GX)-10:0.01:max(GX)+10;

% recover the PDF and CDF via the mixture distribution with two components
k = [0.1, 0.2, 0.3, -0.1, -0.2, -0.3, 0.4].*1;

[ xx ,Pdf,Cdf ] = mixture_of_two_skew_normal_distributions( k,weights,GX,g );


PDF = @(para,g) 2./(para(2)).*normpdf((g-para(1))./para(2),0,1).*normcdf(para(3).*(g-para(1))./para(2),0,1);

PDF1 = xx(1).*PDF(xx(2:4), g );
PDF2 = (1-xx(1)).*PDF(xx(5:7), g );


% failure probability by the proposed method
pf_est = interp1(g,Cdf,0);


%% comparsion 
figure(1)
[j,i]=hist(G,60);
j=j/length(G)/mean(diff(i));
b=bar(i,j,1);
% xlim([-500 700])
hold on
plot(g,Pdf,'r-','LineWidth',1.5)
 plot(g,PDF1,'g--','LineWidth',1.5)
 plot(g,PDF2,'g--','LineWidth',1.5)
h=legend('MCS','Proposed method','Mixture components');
set(h,'Interpreter','latex','FontSize',14)
xlabel('$z$','interpreter','latex','FontSize',25)
ylabel('$\rm PDF$','interpreter','latex','FontSize',14)
set(gca,'FontSize',12);
set(gca,'FontName','Timesnewroman');

figure(2)
gg = min(G):0.05:max(G);
h_mcs = hist(G,gg);
cdf_mcs = cumsum(h_mcs)/sum(h_mcs);
semilogy(gg,cdf_mcs,'b-','LineWidth',1.5)
hold on
semilogy(g,Cdf,'r--','LineWidth',1.5)
h=legend('MCS','Proposed method');
set(h,'Interpreter','latex','FontSize',14)
xlabel('$z$','interpreter','latex','FontSize',16)
ylabel('$\rm CDF$','interpreter','latex','FontSize',14)
ylim([1e-6 1])
set(gca,'FontSize',12);
set(gca,'FontName','Timesnewroman');
grid on

