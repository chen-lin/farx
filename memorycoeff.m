function [ mc ] = memorycoeff( forder, mem )
%MEMORYCOEFF calculates the coefficients in evaluating the fractional order
% derivatives
%  arguments
%   forder = fractional order
%   memory = time shift vector of past values considered
%     forder = [.1 0.2];
%     mem    = 0:1:50;
    norder = length(forder);
    mc     = zeros(norder,length(mem));
    for j = 1:norder
        c1     = (-1).^mem;
        c2     = (gamma(forder(j)+1))./times(gamma(mem+1), ...
            gamma(-mem+forder(j)+1));
        mc(j,:)= times(c1,c2);
    end
end

