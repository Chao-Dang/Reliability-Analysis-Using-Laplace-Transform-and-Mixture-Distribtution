function [ G ] = TDOF_dynamic_syetem( X )
%TDOF_DYNAMIC_SYETEM: A two-degree-of-freedom primary-secondary system subjected to white noise 

kesia = (X(:,5)+X(:,6))./2;
omigas = sqrt(X(:,4)./X(:,2));
omigap = sqrt(X(:,3)./X(:,1));
omigaa = (omigap+omigas)./2;
v = X(:,2)./X(:,1);
eta = (omigap-omigas)./omigaa;

G = X(:,7)-X(:,4).*3.*sqrt((pi.*X(:,8)./(4.*X(:,6).*omigas.^3)).*((kesia.*X(:,6))./(X(:,5).*X(:,6).*(4.*kesia.^2+eta.^2)+v.*kesia.^2)).*((X(:,5).*omigap.^3+X(:,6).*omigas.^3).*omigap./(4.*kesia.*omigaa.^4)));


end

