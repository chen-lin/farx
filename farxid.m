% Perform FARX model identification using the simulation data and also
% compare it to other models
%
% Variables:
%   y_tr  - output for training
%   u_tr  - inputs for training
%   dtr - corresponding time
% .........................................................................
%  Author: L. Chen, Research Fellow, Trinity College Dublin, Dublin 2
%  Email: l.chen.tj@gmail.com, or chenl1@tcd.ie
%  Date: Mar 1, 2016
%  Latest Revision: Feb 10, 2019
% -------------------------------------------------------------------------
clear; clc
load data/Qf.mat; load data/Qh.mat; load data/Ta.mat; load data/Tz.mat;
t = 0:length(Ta)-1;
day = t/240;
strt = 240*120;
lens = 240*92;

st   = .1;
dtr  = day(strt+1:strt+lens);
Qff  = Qf(strt+1:strt+lens);
Qhh  = Qh(strt+1:strt+lens);
Taa  = Ta(strt+1:strt+lens);
y_tr = Tz(strt+1:strt+lens);
Ts   = 0:st:(lens-1)*st;
u_tr = [Qhh Taa Qff];
data_tr = iddata(y_tr, u_tr, st);
y    = Tz;
u    = [Qh Ta Qf];
data = iddata(y, u, st);

%% Use ARX
ord = 4;
tic
sys = arx(data_tr,[1 ones(1,3)*ord ones(1,3)*0]);
toc

[yp, ~, ~] = compare(data_tr, sys);
y1      = yp.Outputdata;
[fit_arx, mae_arx, mse_arx, maxae_arx] = meval(y_tr, y1);
err_arx = y_tr - y1;

%% Use FARX
af = 1;
bt = (0:5)*0.28;   %
L = 100;           % the length for fractional order approximation 
na = length(af);
nb = length(bt);

tic
farx1 = farx(y_tr,u_tr,af,bt,L,'LS');
toc
[~, mae_farx, mse_farx, maxae_farx] = meval(y_tr', farx1.yhat);
fit_farx  = farx1.fit;
err_farx  = farx1.err;

%% Plot results
figure(1)

% set (gcf,'Units','centimeters ','Position',[1,15,24,28], 'color','w');
subplot(211)
hold on
plot(dtr, y_tr, 'r')
plot(dtr, y1, 'k')
plot(dtr, farx1.yhat, 'b--')
hold off
box on

subplot(212)
hold on
plot(dtr, err_arx, 'k')
plot(dtr, err_farx, 'b')
hold off
box on
ylim([-2 2])
display([fit_arx mae_arx mse_arx maxae_arx; fit_farx mae_farx mse_farx ...
    maxae_farx])