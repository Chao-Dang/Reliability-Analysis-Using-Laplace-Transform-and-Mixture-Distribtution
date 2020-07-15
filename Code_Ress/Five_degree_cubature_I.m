function [ points, weights ] = Five_degree_cubature_I( n)
%   n -- dimension of the problem 
%   k -- the free parameter (>0): k = sqrt(2)./2

k = sqrt(2)/2; % this setting corresponds to the formula given in the paper of cubature I


w0 = ((2*k^2-1)^2*n^2-(24*k^4-20*k^2+3)*n+2*(4*k^2-1)^2)/(2*((1-k^2)*n+4*k^2-1).^2);
w1 = (k^4*(4-n))/(2*((1-k^2)*n+4*k^2-1).^2);
w2 = 1./(4*((1-k^2)*n+4*k^2-1).^2);

u0 = 0;
u1 = sqrt((1-k^2)*n+4*k^2-1)/k;
u2 = sqrt((1-k^2)*n+4*k^2-1);

weights = zeros(1,2*n^2+1);
weights(1) = w0;
weights(2:2*n+1) = w1;
weights(2*n+2:end) = w2;


points = zeros(2*n^2+1 ,n );
for i = 2:n+1
    points(i,i-1) = u1;
    points(i+n,i-1) = -u1;
end

A = nchoosek([1:n],2);
for j = 1:nchoosek(n,2)
    points(2*n+1+j,[A(j,:)]) = [u2,u2];
    points(2*n+1+nchoosek(n,2)+j,[A(j,:)]) = [u2,-u2];
    points(2*n+1+2*nchoosek(n,2)+j,[A(j,:)]) = [-u2,u2];
    points(2*n+1+3*nchoosek(n,2)+j,[A(j,:)]) = [-u2,-u2];
end





end

