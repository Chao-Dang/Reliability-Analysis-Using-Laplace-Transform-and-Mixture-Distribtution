function [ weights , points ] = Five_degree_cubature_II( n )
%THREE_COMBINED_WITH_SEVEN_DEGREE_CUBATURE 
%   n^2+3*n+3

w0 = 2/(n+2);
w1 = n^2*(7-n)/(2*(n+1)^2*(n+2)^2);
w2 = 2*(n-1)^2/((n+1)^2*(n+2)^2);


R1 = sqrt(n+2);

%% r = 0
for j = 1:1
    p0(j,:) = zeros(1,n);
    ww0 = w0;
end



%% r = 1
for j = 1:n+1
    p1(j,:) = R1.*fun_A(n,j);
    p1(j+n+1,:) = -R1.*fun_A(n,j);
    ww1(j) = w1;
    ww1(j+n+1) = w1;
end

  temp = nchoosek([1:n+1],2);
for  j = 1:n*(n+1)/2
    b(j,:) = sqrt(n/(2*(n-1))).*(fun_A(n,temp(j,1))+fun_A(n,temp(j,2)));   
    p2(j,:) = R1.*b(j,:);
    p2(j+n*(n+1)/2,:) = -R1.*b(j,:);
    ww2(j) = w2;
    ww2(j+n*(n+1)/2) = w2;  
end


points = [ p0;p1;p2];

weights = [ ww0,ww1,ww2];


end

function [ai] = fun_A(n,k)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
  
    for i = 1:n+1
        for j = 1:n
            if j<i
                A(j,i) = -sqrt((n+1)/(n*(n-j+2)*(n-j+1)));
            elseif i == j
                A(j,i) = sqrt((n+1)*(n-i+1)/(n*(n-i+2)));
            else
                A(j,i) = 0;
            end
        end
    end
    ai = A(:,k)';
    
end
