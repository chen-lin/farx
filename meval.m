function [ fit, mae, mse, maxae ] = meval( y, y_hat )
%MEVAL Summary of this function goes here
%   Detailed explanation goes here
    fit   = 100*(1-norm(y_hat-y)/norm(y-mean(y)));
    mae   = mean(abs(y-y_hat));
    maxae = max(abs(y-y_hat));
    mse   = mean((y-y_hat).^2);
end

