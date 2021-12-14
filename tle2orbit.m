clear n
clc
close all

muE = 398600;

% CP7 (DAVE)              
% 1 43615U 18070C   20311.73749609  .00004108  00000-0  10343-3 0  9994
% 2 43615  93.0291 213.4989 0016794 138.0316 222.2229 15.39816900120365

% ISS (ZARYA)             
% 1 25544U 98067A   20312.62797447  .00001164  00000-0  28917-4 0  9997
% 2 25544  51.6472 358.7040 0002003  92.6677   9.1598 15.49392194254253

% INTELSAT 2-F2           
% 1 02639U 67001A   20312.75954111 -.00000059  00000-0  00000+0 0  9999
% 2 02639   1.9465 287.8917 0009103 316.2065  67.8614  1.00312972 98739
    
CPDAVE = [24946  86.3843 127.9418 0008700 151.1544 209.0134 14.33702974211628];
CPEXO = [33776  86.4036 138.4324 0015334 156.4007 214.2811 14.34129899613840];
Intel = [02639   1.9486 287.8968 0009089 316.0382 116.9384  1.00312936 98719];
GLONAS = [13610  64.0303 137.2256 0008118 199.9412 343.5184  2.14005188297957];

[RDave,VDave] = tle(CPDAVE,muE);
[Rexo, Vexo] = tle(CPEXO, muE);
[Rintel,Vintel] = tle(Intel,muE);
[RGLONAS, VGLONAS] = tle(GLONAS, muE);


radE = 6378;

nT = 5;
step = 50;
format short

[DaveComp] = n2comp(RDave,VDave,nT,muE,step);
[ExoComp] = n2comp(Rexo,Vexo,nT,muE,step);
[IntelComp] = n2comp(Rintel,Vintel,nT,muE,step);
[GLONASComp] = n2comp(RGLONAS,VGLONAS,nT,muE,step);

[row,col] = size(DaveComp);


options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

tspan = [0 24*3600];

state = [RDave;VDave];
[timenew, newstate] = ode45(@TwoBody, tspan, state, options, muE);

state = [Rexo;Vexo];
[timexo, newexo] = ode45(@TwoBody, tspan, state, options, muE);

state = [Rintel;Vintel];
[timeintel, newintel] = ode45(@TwoBody, tspan, state, options, muE);

state = [RGLONAS;VGLONAS];
[timeglo, newglo] = ode45(@TwoBody, tspan, state, options, muE);

% display
fprintf('Radius of CP7 [km]: %f\n',norm(RDave))
fprintf('Velocity of CP7 [km/s]: %f\n',norm(VDave))
fprintf('\nRadius of CP10 [km]: %f\n', norm(Rexo))
fprintf('Velocity of CP10 [km/s]: %f\n',norm(Vexo))
fprintf('\nRadius of Intelsat [km]: %f\n',norm(Rintel))
fprintf('Velocity of Intelsat [km/s]: %f\n',norm(Vintel))


%%
[X,Y,Z] = sphere; % creating sphere (of rad 1) to approximate Earth
r = 6378; % radius of Earth [km]
XE = X*r; % scales sphere to rad 6371 [km]
YE = Y*r;
ZE = Z*r;


%for i = 1:row
i = 69;
figure(1)
plot3(newstate(:,1),newstate(:,2),newstate(:,3),'m','linewidth',2);
hold on
plot3(newexo(:,1),newexo(:,2),newexo(:,3),'b','linewidth',2);
plot3(newintel(:,1),newintel(:,2),newintel(:,3),'r','linewidth',2);
plot3(newglo(:,1),newglo(:,2),newglo(:,3),'g','linewidth',2);
plot3(DaveComp(i,1),DaveComp(i,2),DaveComp(i,3),'-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
plot3(ExoComp(i,1),ExoComp(i,2),ExoComp(i,3),'-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
plot3(IntelComp(i,1),IntelComp(i,2),IntelComp(i,3),'-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
plot3(GLONASComp(i,1),GLONASComp(i,2),GLONASComp(i,3),'-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
surf(XE,YE,ZE);
alpha 0.1
axis equal
hold off
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

legend('CP7 - DAVE Orbit Path (LEO)','Iridium 33 DEB','IntelSat Orbit Path (GEO)','GLONAS Orbit Path (MEO)','CP7 - DAVE','ISS','IntelSat','GLONAS','Earth Model');

%end

function [R,V] = tle(CP7, muE);

inc = CP7(2);
raan = CP7(3);
ecc = CP7(4) /(10^7);
arg = CP7(5);
Me = CP7(6);
n = CP7(7);

Me = deg2rad(Me);

if Me < pi
    E_0 = Me - ecc;
else
    E_0 = Me + ecc;
end

f = @(E) Me - E + ecc*sin(E);
fp = @(E) -1 + ecc*sin(E);

E_1 = E_0 - (f(E_0)/fp(E_0));

err = abs(E_1 - E_0);

while err > 1*10^-8
    E_0 = E_1;
    E_1 = E_0 - (f(E_0)/fp(E_0));
    err = abs(E_1 - E_0);
end

TA = 2*atand((sqrt((1+ecc)/(1-ecc)) * tan(E_1/2)));
TA = mod(TA,360);

T = (n/(24*3600))^-1;

a = (T*sqrt(muE)/(2*pi))^(2/3);

r = a*(1-ecc^2)/(1+ecc*cosd(TA));

h = sqrt(a*muE*(1-ecc^2));

Rmatr = r*[cosd(TA);sind(TA);0];
Vmatr = muE/h * [-sind(TA); ecc + cosd(TA); 0];

Q1 = [cosd(arg) sind(arg) 0; -sind(arg) cosd(arg) 0; 0 0 1];
Q2 = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
Q3 = [cosd(raan) sind(raan) 0; -sind(raan) cosd(raan) 0; 0 0 1];

Q = Q1*Q2*Q3;

R = Q'*Rmatr;
V = Q'*Vmatr;

end

function [Rcomp] = n2comp(R,V,nT,muE,step)
hbar = cross(R,V);
h = norm(hbar);

eccbar = 1/muE*(cross(V,hbar) - muE*R/norm(R));
ecc = norm(eccbar);

rp = h^2/muE * (1/(1+ecc));
ra = h^2/muE * (1/(1-ecc)) ;

a = (ra+rp)/2;

T = 2*pi*a^1.5 / sqrt(muE);

tspan = [0:step:nT*T]';

[row,col] = size(tspan);

Rcomp = zeros(row,4);

for i = 1:row

Me = muE^2 * (1-ecc^2)^1.5 * tspan(i) / h^3;

f = @(E) Me - E + ecc*sin(E);
fp = @(E) -1 + ecc*cos(E);

if Me < pi
    E_0 = Me - ecc;
else
    E_0 = Me + ecc;
end

E_1 = E_0 - (f(E_0)/fp(E_0));
err = abs(E_1 - E_0);

while err > 1*10^-8
    E_0 = E_1;
    E_1 = E_0 - (f(E_0)/fp(E_0));
    err = abs(E_1 - E_0);
end

TA = 2*atand(sqrt((1+ecc)/(1-ecc)) * tan(E_1/2));
TA = mod(TA,360);

rad = h^2/muE * (1/(1+ecc*cosd(TA)));

inc = acosd(hbar(3)/h);

k = [0 0 1];
N = cross(k,hbar);

raan = acosd(N(1)/norm(N));
if N(2) < 0
    raan = 360-raan;
end

arg = acosd(dot(N,eccbar)/(norm(N)*ecc));

Q1 = [cosd(arg) sind(arg) 0; -sind(arg) cosd(arg) 0; 0 0 1];
Q2 = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
Q3 = [cosd(raan) sind(raan) 0; -sind(raan) cosd(raan) 0; 0 0 1];
Q = Q1*Q2*Q3;

Rmatr = rad* [cosd(TA); sind(TA); 0];
Vmatr = muE/h * [-sind(TA); ecc+ cosd(TA); 0];


comp = Q'*Rmatr;
VCOMP = Q'*Vmatr;

Rcomp(i,1) = comp(1);
Rcomp(i,2) = comp(2);
Rcomp(i,3) = comp(3);
Rcomp(i,4) = norm(VCOMP);
end
end

function dstate = TwoBody(time, state, mu)

x = state(1); % defining position elements in state vector
y = state(2);
z = state(3);
vx = state(4); % defining velocity elements in state vector
vy = state(5);
vz = state(6);

rad = norm([x y z]); % def. of "radius"

ax = -mu*x/rad^3; % two body equation
ay = -mu*y/rad^3;
az = -mu*z/rad^3;

dstate = [vx; vy; vz; ax; ay; az]; % new state vector

end