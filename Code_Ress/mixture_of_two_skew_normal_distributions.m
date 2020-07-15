function [ xx ,Pdf,Cdf ] = mixture_of_two_skew_normal_distributions( k,weights,GX,g )
%MIXTURE_OF_TWO_SKEW_NORMAL_DISTRIBUTIONS 

central_m(1) = weights*GX;  % mean value
central_m(2) = sqrt(weights*(GX-central_m(1)).^2);  % standard deviation 
central_m(3) = weights*(GX-central_m(1)).^3./central_m(2).^3; % skewness

para(1) = central_m(1) - sign(central_m(3))*central_m(2)*abs(central_m(3))^(1/3)*(2/(4-pi)).^(1/3) ;
para(2) = central_m(2)*sqrt(1+abs(central_m(3))^(2/3)*(2/(4-pi)).^(2/3)) ;
para(3) = sign(central_m(3))*abs(central_m(3))^(1/3)*(2/(4-pi)).^(1/3)*(2/pi+abs(central_m(3))^(2/3)*((2/(4-pi)).^(2/3))*(2/pi-1) )^(-1/2);

% initial values for the mixture
PA1 = [para(1),para(2),para(3)];
PA2 = [para(1),para(2),para(3)];


x0 = [0.40,PA1,PA2];

%% object function -- by system of equations
fun = @(x)[ analytical_raw_moments(x,k(1))-estimated_raw_moments(weights,GX,k(1));
    analytical_raw_moments(x,k(2))-estimated_raw_moments(weights,GX,k(2));
    analytical_raw_moments(x,k(3))-estimated_raw_moments(weights,GX,k(3));
    analytical_raw_moments(x,k(4))-estimated_raw_moments(weights,GX,k(4));
    analytical_raw_moments(x,k(5))-estimated_raw_moments(weights,GX,k(5));
    analytical_raw_moments(x,k(6))-estimated_raw_moments(weights,GX,k(6));
    analytical_raw_moments(x,k(7))-estimated_raw_moments(weights,GX,k(7));
];

AlGO = {'trust-region-dogleg','trust-region','levenberg-marquardt','trust-region-reflective'};
% 
% opts = optimoptions('fsolve','Algorithm',AlGO{1},'MaxFunctionEvaluations',2e4,'FunctionTolerance',1e-100,'StepTolerance',1e-100,'MaxIterations',2e4,'Display','iter','OptimalityTolerance',1e-18);
% xx = fsolve(fun,x0,opts);

options  = optimoptions('lsqnonlin','Display','iter','Algorithm',AlGO{4},'MaxFunctionEvaluations',2e4,'FunctionTolerance',1e-100,'StepTolerance',2e-100,'MaxIterations',1e4,'OptimalityTolerance',1e-18);
xx = lsqnonlin(fun,x0,[0,-inf,0,-inf,-inf,0,-inf],[1,inf,inf,inf,inf,inf,inf],options);

%% object function -- by optiminzation
% fun = @(x) (analytical_raw_moments(x,k(1))./estimated_raw_moments(weights,GX,k(1)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(2))./estimated_raw_moments(weights,GX,k(2)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(3))./estimated_raw_moments(weights,GX,k(3)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(4))./estimated_raw_moments(weights,GX,k(4)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(5))./estimated_raw_moments(weights,GX,k(5)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(6))./estimated_raw_moments(weights,GX,k(6)) - 1)^2  + ...
%            (analytical_raw_moments(x,k(7))./estimated_raw_moments(weights,GX,k(7)) - 1)^2  ;
% opt = optimset('Display','iter','TolFun',1e-16,'TolX',1e-16, 'MaxFunEvals',1e5,'MaxIter',1e5,'FunValCheck','on'); 
% xx = fminsearch(fun,x0,opt);

% options = saoptimset('Display','iter','TolFun',1e-16,'TemperatureFcn',@temperatureexp,'AnnealingFcn',@annealingboltz,'ReannealInterval',100);
% % lb = [0 -inf 0 -inf 0 -inf 0 -inf -inf 0 -inf];
% % ub = [1 -inf 0 -inf 1 -inf 0 -inf -inf 0 -inf];
% lb = [0 -inf 0 -inf  -inf 0 -inf ];
% ub = [1 inf inf inf  inf inf inf ];
% 
% xx = simulannealbnd(fun,x0,lb,ub,options);

Pdf = xx(1).*skew_normal(xx(2:4), g ) + (1-xx(1)).*skew_normal(xx(5:7), g);
Cdf = cumsum(Pdf)./sum(Pdf);



end


function [raw_a] = analytical_raw_moments(x,k)

% raw_a = x(1)*exp(x(3)/x(2))*sqrt(2*x(3)/pi)*x(2)^(k-1/2)*besselk(1/2-k,x(3)/x(2))+(1-x(1))*exp(k*x(4)+k^2*x(5)^2/2);

raw_a = x(1)*2.*exp(-x(2).*k+x(3)^2*k^2./2).*normcdf(-x(3)*x(4)./sqrt(1+x(4)^2).*k,0,1) + (1-x(1))*2.*exp(-x(5).*k+x(6)^2*k^2./2).*normcdf(-x(6)*x(7)./sqrt(1+x(7)^2).*k,0,1);



end

function [raw_e] = estimated_raw_moments(weights,GX,k)

raw_e = weights*exp(-k.*GX);

end

function [pdf ] = skew_normal(para, g)
%SKEW_NORMAL_DISTIRBUTION 
%   central_m - -  the first-three central moments
pdf = 2./(para(2)).*normpdf((g-para(1))./para(2),0,1).*normcdf(para(3).*(g-para(1))./para(2),0,1);
end


