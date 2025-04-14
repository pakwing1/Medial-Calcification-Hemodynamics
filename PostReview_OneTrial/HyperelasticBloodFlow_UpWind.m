% Media Sclerosis Simulation June 1st 2022
% Units: mmHg - cm - sec
% Uses Upwind scheme to solve 1D blood flow equations
%
% MS_main.m calls this function

%%
PLOT = 0; % PLOT = 1 means plot as you go, PLOT = 0 means don't plot
CLASSICAL = 0; % CLASSICAL = 1 means use classical tube law

%%
if CLASSICAL == 1
    % analytic forms for classical tube law
    E = 3750; % mmHg conversion of Young's modulus from 500 kPa
    beta = (3/2)*sqrt(pi)*E*h;
    
    ffun = @(A,A0) f_classical(A,A0,Pext,beta);
    f_Afun = @(A,A0) f_classical_A(A,A0,beta);
    f_A0fun = @(A,A0) f_classical_A0(A,A0,beta);
    cfun = @(A,A0) c_classical(A,A0,beta,rho);
    nufun = @(A,A0) nu_classical(A,A0,beta,rho);
    nu_invfun = @(y,A0) nu_inv_classical(y,A0,beta,rho);
else
    % interpolants from hyperelastic tube law
    ffun = f_in;
    f_Afun = f_A_in;
    f_B0fun = f_B0_in;
    cfun = c_in;
    nufun = nu_in;
    nu_invfun = nu_inv_in;
end

%% Use f, f_A, f_B0 etc. directly without interpolation - VERY SLOW!
% ffun = @(A,B0)f(A,B0);
% f_Afun = @(A,B0)f_A(A,B0);
% f_B0fun = @(A,B0)f_B0(A,B0);
% cfun = @(A,B0)c(A,B0);
% nufun = @(A,B0)nu(A,B0);
% nu_invfun = @(A,B0)nu_inv(A,B0);


c0 = cfun(A0,B0); %sqrt( beta/2/rho ) .* A0.^(-1/4);
dx = x(2) - x(1);
LAMBDA = max(c0);
dt = 0.4*dx/LAMBDA; % time step
tvec = [0:dt:T];
M = length(tvec);
P_set1 = zeros(1,M);
P_set2 = zeros(1,M);
Q_set1 = zeros(1,M);
Q_set2 = zeros(1,M);
A_set1 = zeros(1,M);
A_set2 = zeros(1,M);
F_set1 = zeros(1,M);
F_set2 = zeros(1,M);
% initial conditions

t = 0;
A = A0;
U0 = 0; %Q_in(0)/A0(1);
Pvec = Pext; % pressure after first resistor ["P_b" for terminal bcs]

U = U0*(1-x/L); %U(1) = Q_in(0)/A0(1);
A_new = zeros(1,N);
U_new = zeros(1,N);

if CLASSICAL == 1
    Wf = U + nufun(A,A0);
    Wb = U - nufun(A,A0);
else    
    Wf = U + nufun(A,B0);
    Wb = U - nufun(A,B0);
end

Wf_new = zeros(1,N);
Wb_new = zeros(1,N);

Wf_temp = zeros(1,N);
Wb_temp = zeros(1,N);

fsolve_options = optimoptions('fsolve','FunctionTolerance',1e-14,'Display','off','Algorithm','levenberg-marquardt');
% fsolve_options = optimoptions('fsolve','Display','off');
% fsolve_options = optimoptions('fsolve','Display','iter');
i = 1;
% while t <= T
while i <= length(tvec)    
    
        t = tvec(i);
        P = ffun(A,B0);
        
        % record (P,Q,A) at x=0 and x=L
        P_set1(i) = P(1);
        P_set2(i) = P(end);
        Q_set1(i) = A(1)*U(1);
        Q_set2(i) = A(end)*U(end);
        A_set1(i) = A(1);
        A_set2(i) = A(end);
        
        Phi = -8*mu*pi*U;
        
        if CLASSICAL == 1
            A0_p = my_gradient(A0,dx);
            dA0 = 1e-6;
            nu_A0 = ( nufun(A,A0+dA0) - nufun(A,A0-dA0) )/2/dA0;
            lambdaf = U+cfun(A,A0);
            lambdab = U-cfun(A,A0);
            RHS1 = Phi/rho./A + ( -1/rho*f_A0fun(A,A0) +  lambdaf.* nu_A0 ) .* A0_p;
            RHS2 = Phi/rho./A + ( -1/rho*f_A0fun(A,A0) -  lambdab.* nu_A0 ) .* A0_p;
        else
            B0_p = my_gradient(B0,dx);
            dB0 = 1e-6;
            nu_B0 = ( nufun(A,B0+dB0) - nufun(A,B0-dB0) )/2/dB0;
                        
            lambdaf = U+cfun(A,B0);
            lambdab = U-cfun(A,B0);
            RHS1 = Phi/rho./A + ( -1/rho*f_B0fun(A,B0) +  lambdaf.* nu_B0 ) .* B0_p;
            RHS2 = Phi/rho./A + ( -1/rho*f_B0fun(A,B0) -  lambdab.* nu_B0 ) .* B0_p;
        end     
            
        F = Phi - my_gradient(rho*U.^2.*A,dx) - A.*my_gradient(P,dx); % vis-a-tergo
        F_set1(i) = F(1);
        F_set2(i) = F(end);
    
        
    if mod(i-1,20) == 0
        sprintf('time t = %f, i = %d',t,i)
        if PLOT == 1
        
        subplot(4,2,1); plot(x,A,x,B0,'k--'); xlabel('x (cm)'); ylabel('A (cm^2)'); %axis([0 L 0 0.8]); 
        subplot(4,2,2); plot(x,U); xlabel('x (cm)'); ylabel('u (cm/s)'); % axis([0 L -1 3]); 
        subplot(4,2,3); plot(x,Wf); xlabel('x (cm)'); ylabel('Wf'); %axis([0 L -600 600]); 
        subplot(4,2,4); plot(x,Wb); xlabel('x (cm)'); ylabel('Wb'); %axis([0 L -600 600]); 
        subplot(4,2,5); plot(x,cfun(A,B0)); xlabel('x'); ylabel('PWV, c (cm/s)');
        subplot(4,2,6); plot(x,A.*U); xlabel('x (cm)'); ylabel('Flow Rate, Q (cm^3/s)'); %axis([0 L -50 50]);
        subplot(4,2,7); plot(x,P); xlabel('x (cm)'); ylabel('P (mmHg)'); % axis([0 L 70 80]);
        subplot(4,2,8); plot(x,F); xlabel('x (cm)'); ylabel('Force p.u. length (mmHg cm)');
        tit = sprintf('time t = %f sec',t);
        title(tit);
        drawnow;
        
        %pause();
        end
    end
        
    % upwind
        if CLASSICAL == 1
            %lambdaf = U+cfun(A,A0);
            %lambdab = U-cfun(A,A0);
            Wf_new(2:N) = Wf(2:N) + ( -lambdaf(2:N) .* (Wf(2:N) - Wf(1:N-1))/dx + RHS1(2:N) )*dt;
            Wb_new(1:N-1) = Wb(1:N-1) + ( -lambdab(1:N-1) .* (Wb(2:N) - Wb(1:N-1))/dx + RHS2(1:N-1) )*dt;
            [~,~,Wf_new(1)] = LeftBC(Wb_new(1),A0,P_in,ffun,nufun,A0(1),t,dt,fsolve_options);
            [~,~,Wb_new(N),Pvec_new] = RightBC(Wf_new(N),A0,A(N),U(N),C,Res0,Res,P0,h,dt,ffun,f_Afun,nufun,Pvec,fsolve_options);           
            A_new = nu_invfun( (Wf_new - Wb_new)/2, A0);
        else
            %lambdaf = U+cfun(A,B0);
            %lambdab = U-cfun(A,B0);
            Wf_new(2:N) = Wf(2:N) + ( -lambdaf(2:N) .* (Wf(2:N) - Wf(1:N-1))/dx + RHS1(2:N) )*dt;
            Wb_new(1:N-1) = Wb(1:N-1) + ( -lambdab(1:N-1) .* (Wb(2:N) - Wb(1:N-1))/dx + RHS2(1:N-1) )*dt;
            [~,~,Wf_new(1)] = LeftBC(Wb_new(1),B0,P_in,ffun,nufun,A(1),t,dt,fsolve_options);
            [~,~,Wb_new(N),Pvec_new] = RightBC(Wf_new(N),B0,A(N),U(N),C,Res0,Res,P0,h,dt,ffun,f_Afun,nufun,Pvec,fsolve_options);                   
            A_new = nu_invfun( (Wf_new - Wb_new)/2, B0);
        end
        
    U_new = (Wf_new + Wb_new)/2;
    A = A_new; U = U_new; Wf = Wf_new; Wb = Wb_new; Pvec = Pvec_new;

    %t = t + dt;
    i = i+1; 

end



function v_x = my_gradient(v,dx)

N = length(v);
v_x(2:N-1) = (v(3:N) - v(1:N-2))/2/dx;
v_x(1) = (-(3/2)*v(1) +2*v(2) -(1/2)*v(3))/dx;
v_x(N) = ((1/2)*v(N-2) -2*v(N-1) + (3/2)*v(N))/dx;


end
function [A_new_1,U_new_1,Wf_new_1] = LeftBC(Wb_new_1,B0,P_in,ffun,nufun,A_1,t,dt,fsolve_options)

    [temp,~,exitflag] = fsolve(@(A) P_in(t+dt) - ffun(A,B0(1)),A_1,fsolve_options); 
    if exitflag <= 0
        sprintf('there was a problem with fsolve in the inlet bc. Exit flag = %d',exitflag)
        pause;
    end
    A_new_1 = temp;
    Wf_new_1 = Wb_new_1 + 2*nufun(A_new_1,B0(1));
    U_new_1 = nufun(A_new_1,B0(1)) + Wb_new_1;

end
function [A_new_N,U_new_N,Wb_new_N,Pvec_new] = RightBC(Wf_new_end,B0,Aend,Uend,C,Res0,Res,P0,h,dt,ffun,f_Afun,nufun,Pvec,fsolve_options)

% Pvec = [Pb]

    %[temp,~,exitflag] = fsolve(@(x) myfun2(x,Wf_new_end,B0(end),Aend,Aend*Uend,C,Res1,Res2,Pout,h,dt,ffun,f_Afun,nufun),[Aend*Uend;Aend],fsolve_options);
    [temp,~,exitflag] = fsolve(@(x) myfun3(x,Wf_new_end,B0(end),Pvec,C,Res0,Res,P0,dt,ffun,nufun),[Aend*Uend;Aend;Pvec],fsolve_options);   
    if exitflag <= 0       
        sprintf('there was a problem with fsolve in the outlet bc. Exit flag = %d',exitflag)
        pause;
    end
    Q_new = temp(1);
    A_new_N = temp(2);
    Pvec_new = temp(3);
    U_new_N = Q_new/A_new_N;
    Wb_new_N = U_new_N - nufun(A_new_N,B0(end));

    
end
function out = myfun3(x,Wf,B0_end,Pvec,C,Res0,Res,Pout,dt,f,nu)
% function for outlet condition
Q_new = x(1); A_new = x(2); Pvec_new = x(3); % = Pb; Pvec and Pvec_new are just a scalar

% Q_end = A_end*U_end;
% Pb = Pvec(1);
   
% new version
        out = [
            A_new*Wf - Q_new - A_new * nu(A_new,B0_end) ;
            %
            C*(Pvec_new - Pvec)/dt + Pvec_new/Res - Pout/Res - Q_new;
            %
            f(A_new,B0_end) - Pvec_new(1) - Q_new*Res0;
            ];   
        
        
        
end










