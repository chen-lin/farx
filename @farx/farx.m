classdef farx
    %FARX Creates an object of fractional order autoregressive model with
    %exogenous input
    %  
    % .....................................................................
    %  Author: L. Chen, Research Associate, Trinity College Dublin, Dublin 2
    %  Email: l.chen.tj@gmail.com, or chenl1@tcd.ie
    %  Date: 01-03-2016
    %  Latest Revision: 12-03-2016
    % ---------------------------------------------------------------------
    properties 
        y          % the output
        u          % the input
        af         % alpha: the predefined orders of the autoregressive model
        na         % number of frational orders on output
        bt         % beta:  the predefined orders of the exogenous input        
        nb         % number of fractional orders on inputs
        L          % number of past values for approximating the GL derivative
        N          % number of samples for system identification
        method     % the method for parameter estimation: LS, RLS, or RIV
        nu         % number of inputs
        afmc       % memory coefficients corresponding to fractional order af
        btmc       % memory coefficients corresponding to fractional order bt
        Y          % memory terms for output
        U          % memory terms for inputs
        a          % identified coefficients corresponding to alpha
        b          % identified coefficients corresponding to beta
        yhat       % predicted output based on the identified model
        err        % prediction error
        fit        % factor to evaluate the modeling accuracy
%         yiv        % instrumental variable in RIV method
%         YIV        % variable to account for long memory effect of yiv
    end

    methods
        % constructor  
        function obj = farx(y,u,af,bt,L,method)
            if nargin < 6; method = 'LS';end;
            obj.method = method;
            [my, Ny]   = size(y);
            [mu, Nu]   = size(u);
            [ma, Na]   = size(af);
            [mb, Nb]   = size(bt);
            if my>Ny; y = y'; [my, Ny]=size(y);end;
            if mu>Nu; u = u'; [mu, Nu]=size(u);end;
            if ma~=Ny; af = af'; [ma, Na] = size(af);end
            if mb~=1;  bt = bt'; [mb, Nb] = size(bt);end
            if ~(Ny==Nu && my==1 && ma==1 && mb==1)
                error('Failed to construct a FARX model, check the data!');
            end
            obj.na = Na;
            obj.nb = Nb;
            obj.y  = y;
            obj.u  = u;
            obj.af = af;
            obj.bt = bt;
            obj.L  = L;
            obj.N  = Ny;
            obj.nu = mu;
            [af_mc, bt_mc] = mc(obj);
            obj.afmc = af_mc;
            obj.btmc = bt_mc;
            [YY, UU] = YU(obj);
            obj.Y    = YY;
            obj.U    = UU;
            switch method
                case 'LS'
                    [af_coef, bt_coef] = lse(obj);
                otherwise
                    [af_coef, bt_coef] = rls(obj);
            end
            obj.a  = af_coef;
            obj.b  = bt_coef;
            [y_hat, errs, yfit] = yest(obj);
            obj.yhat = y_hat;
            obj.err  = errs;
            obj.fit  = yfit;
        end 
        
        % prepare coefficients for evaluate Y and U;
        function [af_mc, bt_mc] = mc(obj) 
            af_mc = memorycoeff(obj.af, 0:1:obj.L);
            bt_mc = memorycoeff(obj.bt, 0:1:obj.L);
        end
       
        % evaluate the parameters depending on past values, Eqs.(9,10)
        function [Y, U] = YU(obj)
            af_mc          = obj.afmc;
            bt_mc          = obj.btmc;
            Y              = zeros(obj.na, obj.N); 
            U              = zeros(obj.nu*obj.nb, obj.N);
            for k = 1:obj.N
                L1   = min(obj.L,k);
                for i = 1:obj.na
                    yL     = obj.y(k+1-1:-1:k-L1+1);
                    Y(i,k) = af_mc(i,2:L1+1)*yL';
                end
                k1 = k-1;
                L1 = min(obj.L,k1);
                for j = 1:obj.nb
                    uL = obj.u(:,k1+1:-1:k1+1-L1);
                    U(j:obj.nb:obj.nu*obj.nb,k1+1) = (bt_mc(j,1:L1+1)*uL')';
                end    
            end
        end
        
        % least squares estimation of the model parameters
        function [af_coef, bt_coef] = lse(obj)
            YY        = obj.Y;
            UU        = obj.U;
            phi       = [-YY(:,1:obj.N-1); UU(:,2:obj.N)]';
            yk        = obj.y(2:obj.N)';
            theta     = (phi'*phi)\(phi'*yk);
            af_coef   = theta(1:obj.na);
            bt_coef   = reshape(theta(obj.na+1:end),length(obj.bt(1,:)),obj.nu);
        end
        
        % recursive least squares estimation of the model parameters
        function [af_coef, bt_coef] = rls(obj)
            YY        = obj.Y;
            UU        = obj.U;
            phi       = [-YY(:,1:obj.N-1); UU(:,2:obj.N)];
            theta     = zeros(obj.na+obj.nb*obj.nu,obj.N);
            % initial parameters for recursive least squares estimation
            F0        = 1000*eye(obj.na+obj.nb*obj.nu);
            lambda    = 1;
            yk        = obj.y(1:obj.N)';
            theta0    = theta(:,1);
            if strcmp(obj.method, 'RLS')
                for k = 1:obj.N-1
                    phik  = phi(:,k);
                    epsk1  = (yk(k+1) - theta0'*phik)/(1+phik'*F0*phik);
                    theta1 = theta0 + F0*phik*epsk1;
                    F1     = 1/lambda*(F0 - (F0*(phik*phik')*F0)/(1+phik'*F0*phik));
                    theta(:,k+1) = theta1;
                    theta0 = theta1;
                    F0     = F1;
                end
            else
                af_mc     = obj.afmc;
                yiv       = zeros(size(obj.y));      % intrumental variable
                Yiv       = zeros(obj.na, obj.N);
                for k = 1:obj.N-1
                    L1   = min(obj.L,k);
                    for i = 1:obj.na
                        yivL     = yiv(k+1-1:-1:k-L1+1);
                        Yiv(i,k) = af_mc(i,2:L1+1)*yivL';
                    end
                    phik  = phi(:,k);
                    psik  = [-Yiv(:,k);UU(:,k+1)];
                    epsk1  = (yk(k+1) - theta0'*phik)/(1+phik'*F0*psik);
                    theta1 = theta0 + F0*psik*epsk1;
                    F1     = 1/lambda*(F0 - (F0*(psik*phik')*F0)/(1+phik'*F0*psik));
                    theta(:,k+1) = theta1;
                    theta0   = theta1;
                    F0       = F1;
                    yiv(k+1) = theta0'*psik;
                end
            end           
            af_coef    = theta(1:obj.na,:);
            bt_coef    = theta(obj.na+1:end, :);
        end
        
        % prediction based on the farx model
        function [y_hat, err, fit] = yest(obj)
            % reconstruction
            af_mc     = obj.afmc;
            y_hat     = zeros(size(obj.y));
            y_hat(1)  = obj.y(1);
            Y_hat     = zeros(obj.na, obj.N); 
            if ~strcmp(obj.method, 'LS')
                theta     = [obj.a(:,obj.N); obj.b(:,obj.N)];
            else
                theta     = [obj.a; reshape(obj.b,[],1)];
            end
            for k = 1:obj.N-1
                L1      = min(k,obj.L);
                for i = 1:obj.na
                    yL_hat      = y_hat(k+1-1:-1:k-L1+1);
                    Y_hat(i,k)  = af_mc(i,2:L1+1)*yL_hat';
                end
                phik_hat   = [-Y_hat(:,k); obj.U(:,k+1)]; 
                y_hat(k+1) = theta'*phik_hat;
            end
            err = obj.y-y_hat;
            fit = 100*(1-norm(obj.y-y_hat)/norm(obj.y-mean(obj.y)));            
        end
        
        % for prediction
        function [yp, err, fit] = pred(obj, y_p, u_p)
            [my, Ny]   = size(y_p);
            [mu, Nu]   = size(u_p);
            if my>Ny; y_p = y_p';[~, Np]   = size(y_p);end;
            if mu>Nu; u_p = u_p';end;
            yp        = zeros(size(y_p));
            yp(:,1)     = y_p(:,1);
            if ~strcmp(obj.method, 'LS')
                theta     = [obj.a(:,obj.N); obj.b(:,obj.N)];
            else
                theta     = [obj.a; reshape(obj.b,[],1)];
            end    
            af_mc          = obj.afmc;
            bt_mc          = obj.btmc;
            YY             = zeros(obj.na, Np); 
            UU             = zeros(obj.nu*obj.nb, Np);
            for k = 1:Np
                k1 = k-1;
                L1 = min(obj.L,k1);
                for j = 1:obj.nb
                    upL = u_p(:,k1+1:-1:k1+1-L1);
                    UU(j:obj.nb:obj.nu*obj.nb,k1+1) = (bt_mc(j,1:L1+1)*upL')';
                end 
            end
            for k = 1:Np-1
                L1   = min(obj.L,k);
                for i = 1:obj.na
                    ypL     = yp(k+1-1:-1:k-L1+1);
                    YY(i,k) = af_mc(i,2:L1+1)*ypL';
                end
                phik_hat  = [-YY(:,k); UU(:,k+1)]; 
                yp(k+1)   = theta'*phik_hat;
            end
            err = y_p-yp;
            fit = 100*(1-norm(yp-y_p)/norm(y_p-mean(y_p)));            
        end
    end
end