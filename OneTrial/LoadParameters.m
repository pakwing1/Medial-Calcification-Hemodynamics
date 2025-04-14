function [N,Nm,T,L,x,P0,Res0,Res,C,rho,mu,h,MechanicalParameters,tapering_factor,ref_rad,pressure_ampl_factor] = LoadParameters(lambda,s)
%
% INPUTS:
% s: the fraction of blocked vessels
% lambda: degree of calcification 
%
% lambda = 0 -- completely healthy
% lambda = 1 -- completely calcified
%
% OUTPUTS:
%
% N: # discretization points along artery
% Nm: not used (obsolete variable)
% T: final time
% L: length of artery in cm
% x: discretization of artery in axial direction (a vector)
% P0: pressure at the end of microcirculation (mmHg)
% Res0, Res, C: parameters of Windkessel Model
% rho: density of blood, in mmHg.s^2/cm^2.
% h: thickness of arterial wall, in cm
% MechanicalParameters: data structure for mechanical parameters of
% arterial wall
% tapering_factor: degree of tapering (1: no tapering, 0: outlet radius is
% 0)
% ref_rad: radius of artery with no pressure
% pressure_ampl_factor: dimensionless scalar that multiplies driving
% pressure



C_T = 0.03301; % cm^3/mmHg
Res0 = 1.05145; % mmHg s/cm^3
Res_T = 3.02434; % mmHg s/cm^3

Res = Res_T/(1-s);
C = C_T*(1-s);

kappa_ip_healthy = 0.12; % dimensionless
kappa_op_healthy = 0.39; % dimensionless
mu_healthy = 114.6; % mmHg
kappai1_healthy = 86.9; % mmHg
kappai2_healthy = 3.54; % dimensionless
varphii_healthy = 49/180*pi; % radians; 0<varphii<pi/2

kappa_ip_calcified = 0.13; % dimensionless
kappa_op_calcified = 0.47; % dimensionless
mu_calcified = 894.4; % mmHg
kappai1_calcified = 649.9; % mmHg
kappai2_calcified = 171.23; % dimensionless
varphii_calcified = 44/180*pi; % radians; 0<varphii<pi/2
   
ref_rad = 0.15; % cm
h = 0.074; % cm; wall thickness
pressure_ampl_factor = 1.0;

MechanicalParameters.kappaip = (1-lambda)*kappa_ip_healthy + lambda*kappa_ip_calcified;
MechanicalParameters.kappaop = (1-lambda)*kappa_op_healthy + lambda*kappa_op_calcified; 
MechanicalParameters.mu = (1-lambda)*mu_healthy + lambda*mu_calcified; % mu has unit mmHg
MechanicalParameters.kappai1 = (1-lambda)*kappai1_healthy + lambda*kappai1_calcified; % kappai1 has unit mmHg
MechanicalParameters.kappai2 = (1-lambda)*kappai2_healthy + lambda*kappai2_calcified; % no unit
MechanicalParameters.varphii = (1-lambda)*varphii_healthy + lambda*varphii_calcified; % varphii has unit radians; 0<varphii<pi/2

%%%%%%%%%%%%%%%%%%%%%%Other parameters not important
MechanicalParameters.lamb=1; % axial stretch; no unit
MechanicalParameters.phi1=0; %opening angle for releasing residual stress. We don't assume residual stress here
%%%%%%%%%%%%%%%%%%%%%%%%%Other parameters not important

T = 2; % Final time in seconds
% mmHg = 133; % number of Pascals in 1 mmHg. And 1 Pa = 1/133 mmHg

L = 40; % cm
N = 30; % Increasing N to about 500 helps boundary spatial oscillations [but does not remove them completely]
Nm = NaN;
x = linspace(0,L,N);

% L = 21; % cm, length of artery
% Nm = 101;
% Lm = 20; % cm, approximate measurement point
% x1 = linspace(0,Lm,Nm); dx = x1(2) - x1(1);
% x2 = [Lm+dx:dx:L];
% x = [x1 x2];
% N = length(x);

% terminal bc parameters
Pa2mmHg = 0.0075; % 1 Pa = 0.0075 mmHg
P0 = 0; % mmHg, pressure at the end of microcirculation
rho = 7.5e-4; % density of blood in mmHg.s^2/cm^2. Conversion from 1000 kg/m^3. Pa = kg /m/s^2
mu = 0 * 0.0075; % dynamic viscosity of blood, mmHg s
tapering_factor = 0.75; % tapering_factor = 1 means no tapering

end



















