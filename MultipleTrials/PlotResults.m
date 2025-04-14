

%%
L1 = sprintf('x=%.2f cm',x_c_inlet); L2 = sprintf('outlet, x=%.2f cm',x_c_outlet);


figure;
subplot(5,1,1); plot(tvec,P_set1,tvec,P_set2,'LineWidth',2); 
xlabel('time (s)'); ylabel('Pressure (mmHg)'); legend(L1,L2);
subplot(5,1,2); plot(tvec,Q_set1,tvec,Q_set2,'LineWidth',2); 
xlabel('time (s)'); ylabel('Flow Rate (cm^3/s)'); legend(L1,L2);
subplot(5,1,3); plot(tvec,A_set1,tvec,A_set2,'LineWidth',2); 
xlabel('time (s)'); ylabel('Area (cm^2)'); legend(L1,L2);
subplot(5,1,4); plot(tvec,Q_set1./A_set1,tvec,Q_set2./A_set2,'LineWidth',2); 
xlabel('time (s)'); ylabel('Velocity (cm/s)'); legend(L1,L2);
subplot(5,1,5); plot(tvec,F_set1,tvec,F_set2,'LineWidth',2);
xlabel('time (s)'); ylabel('Force p.u. length (mmHg.cm)'); legend(L1,L2);


%%

t_short = tvec .* (tvec>T-1); ii = find(t_short == 0); t_short(ii) = []; % only keep times t>1 [only take the average when the system exhibits regular periodic behavior]
% t_short = tvec; ii = [];
P_set1_short = P_set1; P_set1_short(ii) = [];
P_set2_short = P_set2; P_set2_short(ii) = [];
Q_set1_short = Q_set1; Q_set1_short(ii) = [];
Q_set2_short = Q_set2; Q_set2_short(ii) = [];
A_set1_short = A_set1; A_set1_short(ii) = [];
A_set2_short = A_set2; A_set2_short(ii) = [];
U_set1_short = Q_set1./A_set1; U_set1_short(ii) = [];
U_set2_short = Q_set2./A_set2; U_set2_short(ii) = [];
F_set1_short = F_set1; F_set1_short(ii) = [];
F_set2_short = F_set2; F_set2_short(ii) = [];

P_set1_mean = mean(P_set1_short);
P_set2_mean = mean(P_set2_short);
Q_set1_mean = mean(Q_set1_short);
Q_set2_mean = mean(Q_set2_short);
A_set1_mean = mean(A_set1_short);
A_set2_mean = mean(A_set2_short);
U_set1_mean = mean(U_set1_short);
U_set2_mean = mean(U_set2_short);
F_set1_mean = mean(F_set1_short);
F_set2_mean = mean(F_set2_short);


subplot(5,1,1); hold on;

param_tit = sprintf('lambda = %d\n T = %.2f, L = %.0f cm, R_0 = %1.2f mmHg.s/cm^3, R = %1.2f mmHg.s/cm^3, C = %1.3f cm^3/mmHg,\nP_{ext} = %1.2f mmHg, P_0 = %1.2f mmHg, U_0 = %1.2f cm/s, Tapering = %.2f \ndx = %1.2e, dt = %.2e, MM=%d, NN = %d, MMp = %d\n',lambda,T,L,Res0,Res,C,Pext,P0,U0,tapering_factor,dx,dt,MM,NN,MMp);
tit = sprintf('Mean proximal: %1.2f mmHg, Mean distal: %1.2f mmHg',P_set1_mean,P_set2_mean);
plot(t_short,P_set1_short,t_short,P_set2_short,'LineWidth',2); title(param_tit,tit); %axis([0 T 0 200]); 
plot(t_short,P_set1_mean*t_short.^0,t_short,P_set2_mean*t_short.^0);
xlabel('time (s)'); ylabel('Pressure (mmHg)'); legend(L1,L2);

subplot(5,1,2); hold on;
tit = sprintf('Mean proximal: %1.2f cm^3/s, Mean distal: %1.2f cm^3/s',Q_set1_mean,Q_set2_mean);
plot(t_short,Q_set1_short,t_short,Q_set2_short,'LineWidth',2); title('',tit); %axis([0 T -2 4]);
plot(t_short,Q_set1_mean*t_short.^0,t_short,Q_set2_mean*t_short.^0);
xlabel('time (s)'); ylabel('Flow Rate (cm^3/s)'); legend(L1,L2);


subplot(5,1,3); hold on;
tit = sprintf('Mean proximal: %1.2f cm^2, Mean distal: %1.2f cm^2',A_set1_mean,A_set2_mean);
plot(t_short,A_set1_short,t_short,A_set2_short,'LineWidth',2); title('',tit); %axis([0 T 0.69 0.71]);
plot(t_short,A_set1_mean*t_short.^0,t_short,A_set2_mean*t_short.^0);
xlabel('time (s)'); ylabel('Area (cm^2)'); legend(L1,L2);

subplot(5,1,4); hold on;
tit = sprintf('Mean proximal: %1.2f cm/s, Mean distal: %1.2f cm/s',U_set1_mean,U_set2_mean);
plot(t_short,U_set1_short,t_short,U_set2_short,'LineWidth',2); title('',tit); %axis([0 T -5 5]);
plot(t_short,U_set1_mean*t_short.^0,t_short,U_set2_mean*t_short.^0);
xlabel('time (s)'); ylabel('Velocity (cm/s)'); legend(L1,L2);

subplot(5,1,5); hold on;
tit = sprintf('Mean proximal: %1.2f mmHg.cm, Mean distal: %1.2f mmHg.cm',F_set1_mean,F_set2_mean);
plot(t_short,F_set1_short,t_short,F_set2_short,'LineWidth',2); title('',tit); %axis([0 T -5 5]);
plot(t_short,F_set1_mean*t_short.^0,t_short,F_set2_mean*t_short.^0);
xlabel('time (s)'); ylabel('Force p.u. length (mmHg.cm)'); legend(L1,L2);

figure;

subplot(4,2,1); plot(x_c,A,x_c,B0,'k--'); xlabel('x (cm)'); ylabel('A (cm^2)'); %axis([0 L 0 0.8]); 
subplot(4,2,2); plot(x_c,U); xlabel('x (cm)'); ylabel('u (cm/s)'); %axis([0 L -100 100]); 
subplot(4,2,3); plot(x_c,Wf); xlabel('x (cm)'); ylabel('Wf'); %axis([0 L -600 600]); 
subplot(4,2,4); plot(x_c,Wb); xlabel('x (cm)'); ylabel('Wb'); %axis([0 L -600 600]); 
%subplot(4,2,5); plot(x,(Wf - Wb)/2); xlabel('x'); ylabel('(Wf-Wb)/2 (cm/s)'); %axis([0 L 800 1000]); 
subplot(4,2,5); plot(x_c,cfun(A,B0)); xlabel('x'); ylabel('PWV, c (cm/s)');
subplot(4,2,6); plot(x_c,A.*U); xlabel('x (cm)'); ylabel('Flow Rate, Q (cm^3/s)'); %axis([0 L -50 50]);
subplot(4,2,7); plot(x_c,P); xlabel('x (cm)'); ylabel('P (mmHg)'); %axis([0 L 0 400]);
% subplot(4,2,8); plot(x,2*sqrt(A/pi)); xlabel('x (cm)'); ylabel('lumen diameter (cm)'); %axis([0 L 0 1]);
subplot(4,2,8); plot(x_c,F); xlabel('x (cm)'); ylabel('Force p.u. length (mmHg cm)');
% subplot(4,2,8); plot(x,my_gradient(A,dx)); xlabel('x'); ylabel('A_x'); xlim([-1 L+1]);
tit = sprintf('time t = %f sec',t);
title(tit);


%%
figure; hold on;

subplot(2,2,1); hold on;
param_tit = sprintf('lambda = %d, s = %d, T = %.2f, L = %.0f cm\nR_0 = %1.2f mmHg.s/cm^3, R = %1.2f mmHg.s/cm^3\n C = %1.2f cm^3/mmHg, h = %1.2f cm\nP_{ext} = %1.2f mmHg, P_{out} = %1.2f mmHg, U_0 = %1.2f cm/s, Tapering = %.2f, ref-rad = %.2f cm\ndx = %1.2e, dt = %.2e, MM=%d, NN = %d, MMp = %d\n',lambda,s,T,L,Res0,Res,C,h,Pext,P0,U0,tapering_factor,ref_rad,dx,dt,MM,NN,MMp);
tit = sprintf('Mean Inlet: %1.2f mmHg, Mean Outlet: %1.2f mmHg',P_set1_mean,P_set2_mean);
plot(t_short,P_set1_short,'b',t_short,P_set2_short,'r','LineWidth',2); xlim([T-1 T]); 
plot(t_short,P_set1_mean*t_short.^0,'b--',t_short,P_set2_mean*t_short.^0,'r--','LineWidth',2);
xlabel('time (s)'); ylabel('Pressure (mmHg)'); legend(L1,L2);
h1=gca; set(h1,'FontSize',16);
box on; grid on;
title(param_tit,'FontSize',10); 

subplot(2,2,2); hold on;
tit = sprintf('Mean Inlet: %1.2f cm^3/s, Mean Outlet: %1.2f cm^3/s',Q_set1_mean,Q_set2_mean);
plot(t_short,Q_set1_short,'b',t_short,Q_set2_short,'r','LineWidth',2); xlim([T-1 T]); 
plot(t_short,Q_set1_mean*t_short.^0,'b--',t_short,Q_set2_mean*t_short.^0,'r--','LineWidth',2);
xlabel('time (s)'); ylabel('Flow Rate (cm^3/s)'); legend(L1,L2);
h1=gca; set(h1,'FontSize',16);
box on; grid on;

subplot(2,2,3); hold on;
tit = sprintf('Mean Inlet: %1.2f cm^2, Mean Outlet: %1.2f cm^2',A_set1_mean,A_set2_mean);
plot(t_short,A_set1_short,'b',t_short,A_set2_short,'r','LineWidth',2); xlim([T-1 T]); 
plot(t_short,A_set1_mean*t_short.^0,'b--',t_short,A_set2_mean*t_short.^0,'r--','LineWidth',2);
xlabel('time (s)'); ylabel('Area (cm^2)'); legend(L1,L2);
h1=gca; set(h1,'FontSize',16);
box on; grid on;

subplot(2,2,4); hold on;
tit = sprintf('Mean Inlet: %1.2f cm^2, Mean Outlet: %1.2f cm^2',U_set1_mean,U_set2_mean);
plot(t_short,U_set1_short,'b',t_short,U_set2_short,'r','LineWidth',2); xlim([T-1 T]); 
plot(t_short,U_set1_mean*t_short.^0,'b--',t_short,U_set2_mean*t_short.^0,'r--','LineWidth',2);
xlabel('time (s)'); ylabel('Velocity (cm/s)'); legend(L1,L2);
h1=gca; set(h1,'FontSize',16);
box on; grid on;