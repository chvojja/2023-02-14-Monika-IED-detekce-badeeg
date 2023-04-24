function y = gausswindw(L,alpha)
%GAUSSWINUP Summary of this function goes here
%   Detailed explanation goes here
y = gausswin(L,alpha);
%y = flipud(y);
y = flipud(y(1:(L/2)));
end

