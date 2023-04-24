function y = gausswinup(L,alpha)
%GAUSSWINUP Summary of this function goes here
%   Detailed explanation goes here
y = gausswin(L,alpha);
y = y(1:(L/2));
end

