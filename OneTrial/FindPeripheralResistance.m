function CR_params = FindPeripheralResistance

% computes the WindKessel parameters [C Res0 Res]
% that produces an outlet pressure from the simulation
% code HEBF_UpWind_fun.m which is closest in norm to data from distal popliteal artery 
% in Raines (1974) 

MM = 200; % # pts in discretization of X/A that underlies interpolants
NN = 200; % # pts in discretization of B0 that underlies interpolants
MMp = 200; % ceil((ymax-ymin)/(Xmax-Xmin) * MM);

lambda = 0;
s = 0;

[P_in_interp,P_out_interp] = process_inflow_outflow;
P_in = @(t) fnval(P_in_interp,t);

[N,~,T,L,x,Pout,~,~,~,rho,mu,h,MechanicalParameters,tapering_factor,ref_rad,~] = LoadParameters(lambda,s);

Pext = P_in(0);
lumen_radius = ref_rad*(1-x/L) + ref_rad*tapering_factor*x/L; 
B0 = pi*lumen_radius.^2;
A0 = PressurizeVessel(B0,Pext,h,MechanicalParameters);

B0min = 0.95*min(B0); B0max = 1.05*max(B0);
% B0min = min(B0); B0max = max(B0);

P = [linspace(0,200,50) linspace(200,1000,10)]; % calculate up to P=1000 mmHg but you don't need to plot the whole range!
A1 = zeros(1,length(P));
A2 = zeros(1,length(P));
for i=1:length(P)
    sprintf('Calculating Pressure-Area Curves, %d/%d',i,length(P));
    if i==1
        Aguess1 = B0min+B0max; Aguess2 = B0min+B0max;
    else
        Aguess1 = A1(i-1); Aguess2 = A2(i-1);
    end
    A1(i) = fzero(@(A) f_Kun(A, B0max,h, MechanicalParameters) - P(i),Aguess1); %inlet
    A2(i) = fzero(@(A) f_Kun(A, B0min,h, MechanicalParameters) - P(i),Aguess2); %outlet
    
end
Amax = max([A1 A2]);

Xmin = 0; Xmax = Amax-B0max; %Xmin must be ZERO
ymin = 0; ymax = 1500; % ymin must be ZERO

Xgrid = linspace(Xmin,Xmax,MM);
B0grid = linspace(B0min,B0max,NN);
ygrid = linspace(ymin,ymax,MMp);

[f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in] = make_interpolants(B0grid,Xgrid,ygrid,MM,NN,MMp,h,MechanicalParameters,rho);

CR_params0 = [6.4 10 0.02]; % [C Res0 Res]
options = optimset('Display','iter','OutputFcn', @outfun);
lb = [0 0 0]; ub = 30*ones(1,3);
CR_params = fmincon(@(CR_params) ObjectiveFunction(CR_params,P_in,f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in,A0,B0,T,L,N,x,Pout,Pext,rho,mu,h,P_out_interp),CR_params0,[],[],[],[],lb,ub,[],options);

end

function stop = outfun(x, optimValues, state)

stop = false;
fprintf('Current state: [%1.5f, %1.5f, %1.5f]\n',x(1),x(2),x(3))

end
function norm_squared = ObjectiveFunction(CR_params,P_in,f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in,A0,B0,T,L,N,x,Pout,Pext,rho,mu,h,P_out_interp)

[tvec,P_set2] = HEBF_UpWind_fun(CR_params,...
    P_in,f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in,...
    A0,B0,T,L,N,x,Pout,Pext,rho,mu,h);

% a0 = 75.2615;
% a = [-0.0131  -16.4786   -3.9808    3.0870    3.7009    1.1620    0.1862   -0.6369   -0.7611   -0.4006    0.0050    0.1389    0.2444];
% b = [17.4895    7.9162  -11.2509   -4.0684   -0.3901    1.4035    1.2086    0.9641    0.1586   -0.3973   -0.4373   -0.4674    -0.1369];
% P_outlet = @(t) trig_poly(t,a0,a,b); % pressure outlet data that we must match to

dt = tvec(2)-tvec(1);

norm_squared = dt*norm( P_set2 - fnval(P_out_interp,tvec) )^2;

end
function [t_short,P_set2_short] = HEBF_UpWind_fun(CR_params,P_in,f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in,A0,B0,T,L,N,x,Pout,Pext,rho,mu,h)

% Media Sclerosis Simulation June 1st 2022
% Units: mmHg - cm - sec
% Uses Upwind scheme

%%
PLOT = 0; % PLOT = 1 means plot as you go, PLOT = 0 means don't plot
CLASSICAL = 0; % CLASSICAL = 1 means use classical tube law

%%  CR_params = [C Res0 Res]
C = CR_params(1);
Res0 = CR_params(2);
Res = CR_params(3);

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

%%
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
Pvec = Pext; % pressures at junctions in the arterial tree [for terminal bcs]

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
        % sprintf('time t = %f',t)
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
            [~,~,Wb_new(N),Pvec_new] = RightBC(Wf_new(N),A0,A(N),U(N),C,Res0,Res,Pout,h,dt,ffun,f_Afun,nufun,Pvec,fsolve_options);           
            A_new = nu_invfun( (Wf_new - Wb_new)/2, A0);
        else
            %lambdaf = U+cfun(A,B0);
            %lambdab = U-cfun(A,B0);
            Wf_new(2:N) = Wf(2:N) + ( -lambdaf(2:N) .* (Wf(2:N) - Wf(1:N-1))/dx + RHS1(2:N) )*dt;
            Wb_new(1:N-1) = Wb(1:N-1) + ( -lambdab(1:N-1) .* (Wb(2:N) - Wb(1:N-1))/dx + RHS2(1:N-1) )*dt;
            [~,~,Wf_new(1)] = LeftBC(Wb_new(1),B0,P_in,ffun,nufun,A(1),t,dt,fsolve_options);
            [~,~,Wb_new(N),Pvec_new] = RightBC(Wf_new(N),B0,A(N),U(N),C,Res0,Res,Pout,h,dt,ffun,f_Afun,nufun,Pvec,fsolve_options);                   
            A_new = nu_invfun( (Wf_new - Wb_new)/2, B0);
        end
        
    U_new = (Wf_new + Wb_new)/2;
    A = A_new; U = U_new; Wf = Wf_new; Wb = Wb_new; Pvec = Pvec_new;

    %t = t + dt;
    i = i+1;
end

t_short = tvec .* (tvec>T-1); ii = find(t_short == 0); t_short(ii) = []; % only keep times t>1 [only take the average when the system exhibits regular periodic behavior]
% t_short = tvec; ii = [];
P_set1_short = P_set1; P_set1_short(ii) = [];
P_set2_short = P_set2; P_set2_short(ii) = [];

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
function [P_in_interp,P_out_interp] = process_inflow_outflow

% inlet time-pressure data 
A = [    0.0045   65.0829
    0.0145   66.4088
    0.0246   68.6188
    0.0302   72.1547
    0.0369   76.7956
    0.0425   80.5525
    0.0515   85.4144
    0.0615   90.2762
    0.0660   95.3591
    0.0739  101.3260
    0.0806  106.4088
    0.0885  111.2707
    0.0974  116.3536
    0.1119  120.3315
    0.1287  124.7514
    0.1443  126.5193
    0.1587  124.9724
    0.1731  122.7624
    0.1897  119.2265
    0.1997  115.2486
    0.2118  109.9448
    0.2228  104.1989
    0.2305   98.0110
    0.2460   91.1602
    0.2637   86.0773
    0.2859   81.4365
    0.3125   77.0166
    0.3369   75.0276
    0.3713   73.2597
    0.4091   73.9227
    0.4469   76.3536
    0.4880   77.6796
    0.5247   77.9006
    0.5613   76.7956
    0.5869   75.4696
    0.6202   73.7017
    0.6535   72.3757
    0.6846   71.2707
    0.7301   69.7238
    0.7679   69.0608
    0.8046   68.1768
    0.8468   67.2928
    0.8890   66.8508
    0.9301   66.6298
    0.9623   66.4088
    0.9901   66.8508];

B = [    0.0211   61.1680
    0.0468   60.9414
    0.0760   61.4133
    0.1041   63.7454
    0.1158   69.0970
    0.1205   75.8422
    0.1287   81.4255
    0.1357   88.4039
    0.1415   93.5216
    0.1474  100.0345
    0.1544  106.7804
    0.1637  113.0616
    0.1731  119.1103
    0.1836  124.4616
    0.1930  128.8824
    0.2082  130.7464
    0.2257  129.1225
    0.2386  125.4046
    0.2503  121.6864
    0.2620  117.0379
    0.2690  113.0861
    0.2772  109.1345
    0.2842  105.4152
    0.2912  101.4633
    0.3006   97.0469
    0.3099   92.8630
    0.3216   87.9820
    0.3310   82.8679
    0.3450   77.9875
    0.3556   73.5713
    0.3708   69.6214
    0.3860   65.6714
    0.4047   62.4200
    0.4246   59.8664
    0.4526   58.4776
    0.4842   59.4152
    0.5029   62.6754
    0.5181   65.0045
    0.5380   67.7998
    0.5602   69.8980
    0.5813   72.2285
    0.6035   72.9313
    0.6257   72.7039
    0.6515   71.7797
    0.6795   70.1583
    0.7088   68.5372
    0.7310   67.3796
    0.7556   65.9899
    0.7883   64.6022
    0.8164   64.3762
    0.8503   62.9887
    0.8795   62.7630
    0.9064   62.7692
    0.9392   62.0792
    0.9661   62.5505
    0.9977   62.5579];

t1 = A(:,1); p1 = A(:,2); t2 = B(:,1); p2 = B(:,2);

% periodically extend each dataset to [-1,4]

t1 = [t1'-1 t1' t1'+1 t1'+2 t1'+3];
t2 = [t2'-1 t2' t2'+1 t2'+2 t2'+3];
p1 = [p1' p1' p1' p1' p1'];
p2 = [p2' p2' p2' p2' p2'];

%plot(t1,p1,'r.',t2,p2,'b.');

P_in_interp = csaps(t1,p1);
P_out_interp = csaps(t2,p2);

end
function A0 = PressurizeVessel(B0,Pext,h,MechanicalParameters)

% Initial configuration under pressure distribution Pext

% Data from Raines: 0.1963 cm^2 - 0.2827 cm^2.  Deformed areas at 100 mmHg
% B0_fun = @(x) 0.028 - 0.008*(x/L); % cm^2
% B0_p_fun = @(x) -0.01/L * x.^0;
% 
% B0 = B0_fun(x); % cm^2 [Ref: Lorbeer gives max and min radii to be 0.4 and 0.3 cm]
% B0_p = B0_p_fun(x); 
% U0 = 0*ones(1,N); % cm/s
% A0 = zeros(1,N); % cm^2

sprintf('Calculating deformed vessel areas...')
[m,n] = size(B0);
for i=1:m
    for j=1:n
        % compute inflated lumen areas
        A0(i,j) = fzero(@(A) f_Kun(A, B0(i,j),h,MechanicalParameters) - Pext,2*B0(i,j));
    end
end
end
function Pmat = f_Kun(Avec, A0vec, h, MechanicalParameters) % use this line for MS_main.m

%Pmat is the P matrix; pressure unit: mmHg;
%Important: Input A and A0 as a row vector

n1=4000; %n1 is the number of nodes for the wall

kappaip = MechanicalParameters.kappaip;
kappaop = MechanicalParameters.kappaop;
mu = MechanicalParameters.mu;
kappai1 = MechanicalParameters.kappai1;
kappai2 = MechanicalParameters.kappai2;
varphii = MechanicalParameters.varphii;
lamb = MechanicalParameters.lamb;
phi1 = MechanicalParameters.phi1;

Ap=2*kappaip*kappaop;
Bp=2*kappaop*(1-2*kappaip);

R1vec=sqrt(A0vec/pi);


A1c = length(Avec);
A0c = length(A0vec);
Pmat=zeros(A0c,A1c);

for rn=1:A0c %row number
    R1=R1vec(rn);
    R2=R1+h;
    h1=(R2-R1)/(n1-1); %subinterval length
    Rv1=linspace(R1, R2, n1);

    for cn=1:A1c %col number
        A=Avec(cn);
        r2=sqrt( A/pi + R2^2-R1^2);
        Pmat(rn, cn)=Find_P();

    end

end


%........(begin) Find P using bisection method..................................
    function final_P=Find_P()
%         FUN=@f_function;
%         final_P=fzero(FUN, 50);
[Trr_MINUS_Ttt]=FindT_no_p();
[rv1]=r_vector(r2);
rv1_1dev=Derivative(rv1); %first derivative
Integrands=1./rv1.*rv1_1dev.*Trr_MINUS_Ttt;
final_P= - Integration(1, n1, Integrands);
    end
%........(end)Find P..........................................


% %...............(begin) T2rr=T3rr at r3 equation, find the left and right sides of the equation.........................................................
%     function [diff]=f_function(Pressure)
%         [Tr1, Tt1, Tz1]=FindT(Pressure);
%         f_inner=Tr1(n1);
%         f_outer=0;
%         diff=f_inner-f_outer;
%     end
%..................(end) T2rr=T3rr at r3 equation, find the left and right sides of the equation.........................................................



%...................(start) Find all Cauchy Stress Tensor.....
    function [Trr_MINUS_Ttt]=FindT_no_p()

        rv1=r_vector(r2);
        rv1_1dev=Derivative(rv1); %first derivative

        I1 = rv1_1dev.^2  +  pi^2*rv1.^2./( (pi-phi1)^2 * Rv1.^2 )  +  lamb^2;
        I4 = pi^2*rv1.^2*(cos(varphii))^2./((pi-phi1)^2*Rv1.^2)  +  lamb^2*(sin(varphii))^2;
        Ei=Ap*I1 + Bp*I4 + (1-3*Ap-Bp)*rv1_1dev.^2 - 1;

        Trr_no_p = mu*rv1_1dev.^2  + 4*kappai1*exp(kappai2*Ei.^2).*Ei*(1-2*Ap-Bp).*rv1_1dev.^2;
        Ttt_no_p = mu*pi^2*rv1.^2./((pi-phi1)^2*Rv1.^2)  + 4*kappai1*exp(kappai2*Ei.^2).*Ei*pi^2.*rv1.^2./((pi-phi1)^2*Rv1.^2)*( Ap+Bp*(cos(varphii))^2 );
        Trr_MINUS_Ttt=Trr_no_p-Ttt_no_p;

%         inner_vec=1./rv1.*rv1_1dev.*(Trr_no_p-Ttt_no_p);
%         Trr_with_p=0*Trr_no_p;
% 
%         for i=1:n1
%             Trr_with_p(i)=  -Integration(1, i, inner_vec) - Pressure;
%         end
% 
%         Pv1=Trr_no_p-Trr_with_p; % the hydrostatic parameter vector p
% 
%         Tr1=Trr_with_p;
%         Tt1 = -Pv1 + Ttt_no_p;
%         Tz1= -Pv1 + mu*lamb^2 + 4*kappai1*exp(kappai2*Ei.^2).*Ei.*lamb^2*( Ap+Bp*(sin(varphii))^2 );

    end
%...................(end) Find all Cauchy Stress Tensor...



%.................(start)Find rv1 as a vector of r in the three layers
    function [rv1]=r_vector(r2)

        rv1=sqrt(   r2^2 - (pi-phi1)/(pi*lamb)*(R2^2 - Rv1.^2)  );
    end
%.................(end)Find rv1 as a vector of r in the three layers



%.................(begin) Integration Function..........................................................
    function int_val=Integration(lower_value, upper_value, vector)%This is for integration.
        %int_val is only a number, not a vector
        %lower_value is the number of lower node starting for integration,
        %upper_value is the number of upper most node we finish integration.
        %vector is the integrand vector;

        hl1=h1;

        if lower_value==upper_value
            int_val=0;
        elseif upper_value==lower_value+1
            int_val=1/2*(vector(lower_value)+vector(upper_value));
        else
            int_val=1/2*(vector(lower_value)+vector(upper_value));
            for gg=1:upper_value-lower_value-1
                int_val=int_val+vector(lower_value+gg);
            end
        end

        int_val=hl1*int_val;
    end
%.................(end) Integration Function..........................................................



%.................(begining) derivative function...............................
    function der_vec=Derivative(vector)%This is for derivative.
        %der_vec is a vector of the derivative of the input vector.
        %vector is the original vector needing differentiation;
        hl1=h1; %subinterval length;
        noden=n1;
        der_vec=0*vector;
        der_vec(1)=(vector(2)-vector(1))/hl1;
        der_vec(noden)=(vector(noden)-vector(noden-1))/hl1;

        for i1=2:noden-1
            der_vec(i1)=(vector(i1+1)-vector(i1-1))/(2*hl1);
        end
    end
%.................(end) derivative function...............................


end
function [f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in] = make_interpolants(B0,X,y,M,N,Mp,h,MechanicalParameters,rho)
    
    %y_min = yvec(1); y_max = yvec(2);
    %yvals = linspace(y_min,y_max,MMp);
    
    PLOT = 0;
    DISPLAY = 0;
    
    max_lumen_pressure = 250; % mmHg
    
    fhat_vals = zeros(M,N);
    fhat_A_vals = zeros(M,N);
    fhat_B0_vals = zeros(M,N);
    chat_vals = zeros(M,N);
    nuhat_vals = zeros(M,N);
    nuhat_inv_vals = zeros(Mp,N);
    
    if DISPLAY == 1
        sprintf('computing fhat, fhat_A, fhat_B0 and c_hat on 2D grid of (X,B0)')
    end
    parfor i=1:M
        for j=1:N
            [i j];
            fhat_vals(i,j) = f(B0(j) + X(i),B0(j),h,MechanicalParameters);
            fhat_A_vals(i,j) = f_A(B0(j) + X(i),B0(j),h,MechanicalParameters);
            fhat_B0_vals(i,j) = f_B0(B0(j) + X(i),B0(j),h,MechanicalParameters);
            chat_vals(i,j) = c(B0(j) + X(i),B0(j),h,MechanicalParameters,rho);
        end
    end
    
    
    method = 'cubic';
%     [Xvals,B0vals] = meshgrid(X,B0);
%     [yvals,B0vals] = meshgrid(y,B0);
%     fhat_in = @(X,B0) interp2(Xvals,B0vals,fhat_vals,X,B0,method); 
    

    % compute nuhat
    if DISPLAY == 1
        sprintf('computing nuhat on 2D grid of (X,B0)')
    end
    options = odeset('AbsTol',1e-14,'RelTol',1e-12);
    parfor j=1:N
        j;
        [~,nuhat]=ode15s(@(Xp,y) c(Xp+B0(j),B0(j),h,MechanicalParameters,rho)/(Xp+B0(j)),X,0,options);
        nuhat_vals(:,j) = nuhat;
    end
    
    % compute inverse of nuhat
    if DISPLAY == 1
        sprintf('computing nuhat inverse on 2D grid of (y,B0)')
    end
    parfor j=1:N
        j;
        [~,Xfun]=ode15s(@(y,X) (X+B0(j))/c(X+B0(j),B0(j),h,MechanicalParameters,rho),y,0,options);
        nuhat_inv_vals(:,j) = B0(j) + Xfun;
    end
    
        
    % interpolate
    [Xvals,B0vals1] = meshgrid(X,B0);
    [yvals,B0vals2] = meshgrid(y,B0);
%     fhat_in = interp2(Xvals,B0vals,fhat_vals,method);
%     fhat_A_in = interp2(Xvals,B0vals,fhat_A_vals,method);
%     fhat_B0_in = interp2(Xvals,B0vals,fhat_B0_vals,method);
%     chat_in = interp2(Xvals,B0vals,chat_vals,method);
%     nuhat_in = interp2(Xvals,B0vals,nuhat_vals,method);
%     nuhat_inv_in = interp2(yvals,B0vals,nuhat_inv_vals,method);
    
    % interpolants
    f_in = @(A,B0) interp2(Xvals,B0vals1,fhat_vals',A-B0,B0,method); 
    f_A_in = @(A,B0) interp2(Xvals,B0vals1,fhat_A_vals',A-B0,B0,method);
    f_B0_in = @(A,B0) interp2(Xvals,B0vals1,fhat_B0_vals',A-B0,B0,method);
    c_in = @(A,B0) interp2(Xvals,B0vals1,chat_vals',A-B0,B0,method);
    nu_in = @(A,B0) interp2(Xvals,B0vals1,nuhat_vals',A-B0,B0,method);
    nu_inv_in = @(y,B0) interp2(yvals,B0vals2,nuhat_inv_vals',y,B0,method);
    
    %f_in = @(A,B0) fhat_in(A-B0,B0);
    %f_in = @(A,B0) interp2(Xvals,B0vals,fhat_vals,A-B0,B0,method);
%     f_A_in = @(A,B0) interp2(Xvals,B0vals,fhat_A_vals,A-B0,B0,method);
%     f_B0_in = @(A,B0) interp2(Xvals,B0vals,fhat_B0_vals,A-B0,B0,method);
%     c_in = @(A,B0) interp2(Xvals,B0vals,chat_vals,A-B0,B0,method);
%     nu_in = @(A,B0) interp2(Xvals,B0vals,nuhat_vals,A-B0,B0,method);
%     nu_inv_in = @(y,B0) interp2(yvals,B0vals,nuhat_inv_vals,y,B0,method);
    


    if PLOT == 1
    % plot over a coarse mesh
    figure(2);
    subplot(2,3,1); surf(X,B0,fhat_vals'); xlabel('X'); ylabel('B_0'); zlabel('f_{hat}'); zlim([0 max_lumen_pressure]);
    subplot(2,3,2); surf(X,B0,fhat_A_vals'); xlabel('X'); ylabel('B_0'); zlabel('f_{hatA}'); 
    subplot(2,3,3); surf(X,B0,fhat_B0_vals'); xlabel('X'); ylabel('B_0'); zlabel('f_{hatB0}'); 
    subplot(2,3,4); surf(X,B0,chat_vals'); xlabel('X'); ylabel('B_0'); zlabel('c_{hat}'); 
    subplot(2,3,5); surf(X,B0,nuhat_vals'); xlabel('X'); ylabel('B_0'); zlabel('nu_{hat}'); 
    subplot(2,3,6); surf(y,B0,nuhat_inv_vals'); xlabel('y'); ylabel('B_0'); zlabel('nu^{-1}_{hat}');
    
    % plot over a dense mesh
    figure(3);
    [Xdense,B0dense] = meshgrid(linspace(min(X),max(X),100),linspace(min(B0),max(B0),100));
    [ydense,~] = meshgrid(linspace(min(y),max(y),100),linspace(min(B0),max(B0),100));
    fhat_dense = f_in(Xdense+B0dense,B0dense);
    fhat_Adense = f_A_in(Xdense+B0dense,B0dense);
    fhat_B0dense = f_B0_in(Xdense+B0dense,B0dense);
    chat_dense = c_in(Xdense+B0dense,B0dense);
    nuhat_dense = nu_in(Xdense+B0dense,B0dense);
    nuhat_inv_dense = nu_inv_in(ydense,B0dense);
    %
    subplot(2,3,1); surf(Xdense,B0dense,fhat_dense); xlabel('X'); ylabel('B_0'); zlabel('f_{hat}'); zlim([0 max_lumen_pressure]);
    subplot(2,3,2); surf(Xdense,B0dense,fhat_Adense,'EdgeColor','None'); xlabel('X'); ylabel('B_0'); zlabel('f_{hatA}'); 
    subplot(2,3,3); surf(Xdense,B0dense,fhat_B0dense,'EdgeColor','None'); xlabel('X'); ylabel('B_0'); zlabel('f_{hatB0}');
    subplot(2,3,4); surf(Xdense,B0dense,chat_dense,'EdgeColor','None'); xlabel('X'); ylabel('B_0'); zlabel('c_{hat}'); 
    subplot(2,3,5); surf(Xdense,B0dense,nuhat_dense,'EdgeColor','None'); xlabel('X'); ylabel('B_0'); zlabel('nu_{hat}');
    subplot(2,3,6); surf(ydense,B0dense,nuhat_inv_dense,'EdgeColor','None'); xlabel('y'); ylabel('B_0'); zlabel('nu^{-1}_{hat}');
    end
end
function out = f(A,B0,h,MechanicalParameters)

if length(A) ~= length(B0)
    sprintf('A and B0 do not have the same length...aborting in f')
    pause;
end

    [m,n] = size(A);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            out(i,j) = f_Kun(A(i,j),B0(i,j),h,MechanicalParameters);
            % out(i) = ( sqrt(A(i)) - sqrt(B0(i)) )/ B0(i);
        end
    end
end
function out = f_A(A,B0,h,MechanicalParameters) 

if length(A) ~= length(B0)
    sprintf('A and B0 do not have the same length...aborting')
    pause;
end

    dA = 1e-6;
    [m,n] = size(A);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            out(i,j) = ( f(A(i,j)+dA,B0(i,j),h,MechanicalParameters) - f(A(i,j)-dA,B0(i,j),h,MechanicalParameters) )/2/dA;
        end
    end

end
function out = f_B0(A,B0,h,MechanicalParameters)

if length(A) ~= length(B0)
    sprintf('A and B0 do not have the same length...aborting')
    pause;
end

    dB0 = 1e-6;
    [m,n] = size(A);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            out(i,j) = ( f(A(i,j),B0(i,j)+dB0,h,MechanicalParameters) - f(A(i,j),B0(i,j)-dB0,h,MechanicalParameters) )/2/dB0;
        end
    end

end
function out = c(A,B0,h,MechanicalParameters,rho)

if length(A) ~= length(B0)
    sprintf('A and B0 do not have the same length...aborting in c')
    pause;
end

%global rho
    % analytic solution for usual tube law: 
    % out = sqrt(beta/2/rho./A0) .* 4* ( A.^(1/4) - A0^(1/4) );
    [m,n] = size(A);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            out(i,j) = sqrt(A(i,j)/rho.* f_A(A(i,j),B0(i,j),h,MechanicalParameters) );
        end
    end

end






