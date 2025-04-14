% MS_loop1.m, MS_loop2.m, MS_loop3.m and MS_loop4.m all run the
% hemodynamics code HyperelasticBloodFlow_UpWind.m but iterate over
% different pairs of parameters.
%
% - MS_loop1.m: iterates over (lambda, s) [degree of calcification (DoC) and
% fraction of arteriole closures]
%
%  - MS_loop2.m: iterates over DoC, lambda. 
% But reference radius, pressure
% amplification factor and closure fraction are also functions of lambda
%
% - MS_loop3.m: iterates over DoC and reference radius
%
% - MS_loop4.m: iterates over DoC and pressure amplification factor, beta
%
% All the MS_loop functions call the following (separate) m files:
% --> LoadParameters.m
% --> HyperelasticBloodFlow_UpWind.m
% --> PlotResults.m




% Compute max, min and mean of (Pressure, Flow Rate, Area) as functions of (lambda,s)
% produces colormaps

% size of lambda_vec and s_vec
MMM = 4; 
NNN = 4;

%% mesh sizes for interpolants
MM = 200; % # pts in discretization of X/A that underlies interpolants
NN = 200; % # pts in discretization of B0 that underlies interpolants
MMp = 200; % ceil((ymax-ymin)/(Xmax-Xmin) * MM);

lambda_vec = linspace(0,1,MMM);
pressure_ampl_factor_vec = linspace(1,1,MMM);
ref_rad_vec = linspace(0.15,0.15,MMM);
% pressure_ampl_factor_vec = linspace(1,1.6,MMM);
% ref_rad_vec = linspace(0.15,0.21,MMM);
s_vec = linspace(0.3,0.85,NNN);





DELTA = zeros(MMM,NNN); MEAN_deltaP = zeros(MMM,NNN);
MAX_Q = zeros(MMM,NNN); MEAN_Q = zeros(MMM,NNN);
MEAN_A1 = zeros(MMM,NNN); MEAN_A2 = zeros(MMM,NNN);



for iii=1:length(lambda_vec)

        %% Load Parameters (except for WK parameters)
        [N,nm,T,L,x,P0,~,~,~,rho,mu,...
            h,MechanicalParameters,tapering_factor,~,~] = LoadParameters(lambda_vec(iii),NaN);
            
        pressure_ampl_factor = pressure_ampl_factor_vec(iii);
        ref_rad = ref_rad_vec(iii);

        
        %% set up inlet and hydrostatic pressures
        P_in = LoadInPressure;
        P_in = @(t) pressure_ampl_factor*P_in(t);
        Pext = P_in(0);

        %% set up arterial geometry
        lumen_radius = ref_rad*(1-x/L) + ref_rad*tapering_factor*x/L; 
        B0 = pi*lumen_radius.^2; B0min = 0.95*min(B0); B0max = 1.05*max(B0);
        A0 = PressurizeVessel(B0,Pext,h,MechanicalParameters);

        %% Make interpolants
        Pmaximum = 1000; % mmHg
        Amax = fzero(@(A) f_Kun(A, B0max,h,MechanicalParameters) - Pmaximum,A0(1));
        Xmin = 0; Xmax = Amax-B0max; % Xmin must be ZERO
        ymin = 0; ymax = 1500; % ymin must be ZERO
        Xgrid = linspace(Xmin,Xmax,MM);
        B0grid = linspace(B0min,B0max,NN);
        ygrid = linspace(ymin,ymax,MMp);
        [f_in,f_A_in,f_B0_in,c_in,nu_in,nu_inv_in] = make_interpolants(B0grid,Xgrid,ygrid,MM,NN,MMp,h,MechanicalParameters,rho);
    
        for jjj=1:length(s_vec)
            sprintf('iii = %d/%d, jjj = %d/%d, (lambda,s) = (%d,%d)',iii,MMM,jjj,NNN,lambda_vec(iii),s_vec(jjj))
            % modify WK parameters
            [~,~,~,~,~,~,Res0,Res,C,~,~,~,~,~,~,~] = LoadParameters(lambda_vec(iii),s_vec(jjj));
            
            HyperelasticBloodFlow_UpWind;
            t_short = tvec .* (tvec>T-1); ii = find(t_short == 0); % only keep times t>1 [only take the average when the system exhibits regular periodic behavior]

            P_set1_short = P_set1; P_set1_short(ii) = [];
            Q_set1_short = Q_set1; Q_set1_short(ii) = [];
            A_set1_short = A_set1; A_set1_short(ii) = [];
            
            P_set2_short = P_set2; P_set2_short(ii) = [];
            Q_set2_short = Q_set2; Q_set2_short(ii) = [];
            A_set2_short = A_set2; A_set2_short(ii) = [];

            % flow metrics to plot
            DELTA(iii,jjj) = max(P_set2_short) - max(P_set1_short); % difference in systolic pressures
            MEAN_deltaP(iii,jjj) = mean(P_set2_short) - mean(P_set1_short); 
            MAX_Q(iii,jjj) = max(Q_set2_short);
            MEAN_Q(iii,jjj) = mean(Q_set2_short);
            MEAN_A2(iii,jjj) = mean(A_set2_short); % mean area at outlet
            MEAN_A1(iii,jjj) = mean(A_set1_short); % mean area at inlet
        end
end

%%
% subplot(2,3,1); surf(lambda_vec,s_vec,DELTA','LineStyle','None'); title('\Delta = Outlet Max Press. - Inlet Max Press. (mmHg)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);
% subplot(2,3,4); surf(lambda_vec,s_vec,MEAN_deltaP','LineStyle','None'); title('Outlet Mean Press. - Inlet Mean Press. (mmHg)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);
% 
% subplot(2,3,2); surf(lambda_vec,s_vec,MAX_Q','LineStyle','None'); title('Max Outlet Flow Rate (cm^3/s)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);
% subplot(2,3,5); surf(lambda_vec,s_vec,MEAN_Q','LineStyle','None'); title('Mean Outlet Flow Rate (cm^3/s)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);
% 
% subplot(2,3,3); surf(lambda_vec,s_vec,MEAN_A1','LineStyle','None'); title('Mean Inlet Lumen Area (cm^2)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);
% subplot(2,3,6); surf(lambda_vec,s_vec,MEAN_A2','LineStyle','None'); title('Mean Outlet Lumen Area (cm^2)'); colorbar; xlabel('\lambda'); ylabel('s'); view(2);

subplot(1,3,1); surf(pressure_ampl_factor_vec,ref_rad_vec,DELTA');  title('\Delta = Outlet Max Press. - Inlet Max Press. (mmHg)'); colorbar; xlabel('\alpha'); ylabel('R_a (cm)'); view(2);
subplot(1,3,2); surf(pressure_ampl_factor_vec,ref_rad_vec,MEAN_deltaP'); title('Outlet Mean Press. - Inlet Mean Press. (mmHg)'); colorbar; xlabel('\alpha'); ylabel('R_a (cm)'); view(2);
subplot(1,3,3); surf(pressure_ampl_factor_vec,ref_rad_vec,MEAN_Q'); title('Mean Outlet Flow Rate (cm^3/s)'); colorbar; xlabel('\alpha'); ylabel('R_a (cm)'); view(2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF FUNCTIONS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
function P_in = LoadInPressure

% Case 1) Fourier interpolation - using Figure 1 from Mohajer - see MohajerData.fig
% a0 = 59.6327;
% a =  [6.7699   -1.8622   -5.8850   -2.7353   -0.7109   -0.1563   -0.2699];
% b = [6.0590   11.3548    1.3139   -1.0616   -0.9991   -0.3302   -0.1777];

% Case 2) Fourier interpolation - using Figure 12 of Raines - see
% InletOutletPressures.fig
%
%a0 = 78.5471;
%a = [3.6185   -9.0953  -10.1179   -1.6387    0.8341    1.3620    0.5244    0.5791    0.4257    0.2691];
%b = [16.2956   13.2335   -2.8412   -4.9046   -2.7385   -0.6774   -0.0717    -0.0434    0.3606    0.4084];
% 
% % normal pressure input
% P_in = @(t) trig_poly(t,a0,a,b); % Pressure conditions (mmHg)

% Case 3) constant pressure at inlet
% P_in = @(t) 71.2179*t.^0;

% Case 4) "impulse" pressure at inlet
% a = 0.1;
%P_in = @(t) 71.2179 + 10*( sin(t/a) .*(t<a) );

% Case 5) Spline interpolation
[P_in_interp,~] = process_inflow_outflow;
P_in = @(t) fnval(P_in_interp,t); % remember: t < 4 !

% sprintf('Inflow function P_in(t) loaded')
%Q_in = @(t) 0;
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

% sprintf('Calculating deformed vessel areas...')
[m,n] = size(B0);
for i=1:m
    for j=1:n
        % compute inflated lumen areas
        A0(i,j) = fzero(@(A) f_Kun(A, B0(i,j),h,MechanicalParameters) - Pext,2*B0(i,j));
    end
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
function out = nu(A,B0,h,MechanicalParameters,rho)

opts = odeset('AbsTol',1e-14,'RelTol',1e-12);

if length(A) ~= length(B0)
    sprintf('A and B0 do not have the same length...aborting in nu')
    pause;
end

%global rho
    % analytic solution for usual tube law:
%     out = sqrt(beta/2/rho./A0) .* 4* ( A.^(1/4) - A0^(1/4) );

    [m,n] = size(A);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            % dnu/dA = c(A)/A   s.t. nu(B0) = 0
            if A(i,j) == B0(i,j)
                out(i,j) = 0;
            else
                [~,temp] = ode45(@(Ap,y) c(Ap,B0(i,j),h,MechanicalParameters,rho)/Ap,[B0(i,j),A(i,j)],0,opts);
                out(i,j) = temp(end);
            end
        end
    end

end
function out = nu_inv(y,B0,h,MechanicalParameters,rho)

opts = odeset('AbsTol',1e-14,'RelTol',1e-12);

if length(y) ~= length(B0)
    sprintf('y and B0 do not have the same length...aborting in nu_inv')
    pause;
end
    [m,n] = size(y);
    out = zeros(m,n);
    for i=1:m
        for j=1:n
            if y(i,j) == 0
                out(i,j) = B0(i,j);
            else
                [~,temp] = ode15s(@(yp,A) A/c(A,B0(i,j),h,MechanicalParameters,rho),[0,y(i,j)],B0(i,j),opts);
                out(i,j) = temp(end);
            end
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

