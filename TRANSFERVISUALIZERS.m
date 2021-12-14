%% AERO 351 FINAL PROJECT

% HAYDEN BUSS
% FRANCISCO LEON-GOMEZ
% HECTOR DELGADO MARQUEZ
% GAGANDEEP THAPAR

% FALL QUARTER 2020

%% HOUSEKEEPING
clear
clc
close all
muE = 398600;

%% OBJECTS OF INTEREST
% IRIDIUM 33 (LEO-1) - Iridium debris
% 1 24946U 97051C   20312.77658151  .00000096  00000-0  27434-4 0  9994
% 2 24946  86.3843 127.9418 0008700 151.1544 209.0134 14.33702974211628
% 
% IRIDIUM 33 (LEO - 2) - Iridium debris
% 1 33776U 97051P   20312.80323369  .00000235  00000-0  76976-4 0  9999
% 2 33776  86.4036 138.4324 0015334 156.4007 214.2811 14.34129899613840
% 
% GLONAS (MEO) - Rocket Body
% 1 13610U 82100H   20312.07351556  .00000096  00000-0  00000-0 0  9997
% 2 13610  64.0303 137.2256 0008118 199.9412 343.5184  2.14005188297957
% 
% INTELSAT 2-F2 (GEO) - Defunt GTO sat          
% 1 02639U 67001A   20312.75954111 -.00000059  00000-0  00000+0 0  9999
% 2 02639   1.9465 287.8917 0009103 316.2065  67.8614  1.00312972 98739

iri1tle =   [24946  86.3843 127.9418 0008700 151.1544 209.0134 14.33702974211628];
iri2tle =   [33776  86.4036 138.4324 0015334 156.4007 214.2811 14.34129899613840];
glonastle = [13610  64.0303 137.2256 0008118 199.9412 343.5184  2.14005188297957];
inteltle =  [02639   1.9465 287.8917 0009103 316.2065  67.8614  1.0031297298739];

%% FINDING R, V VECTORS; COES FROM TLE

[iri1R, iri1V] =        tle2RV(iri1tle, muE);
[iri2R, iri2V] =        tle2RV(iri2tle, muE);
[glonasR, glonasV] =    tle2RV(glonastle, muE);
[intelR, intelV] =      tle2RV(inteltle, muE);

[iri1h, iri1vr, iri1inc, iri1raan, iri1ecc,...
    iri1arg, iri1TA, iri1ra, iri1rp, iri1a, iri1T] =    RV2COE(iri1R, iri1V, muE);
[iri2h, iri2vr, iri2inc, iri2raan, iri2ecc,...
    iri2arg, iri2TA, iri2ra, iri2rp, iri2a, iri2T] =    RV2COE(iri2R, iri2V, muE);
[glonash, glonasvr, glonasinc, glonasraan,...
    glonasecc, glonasarg, glonasTA, glonasra,...
    glonasrp, glonasa, glonasT] =                       RV2COE(glonasR, glonasV, muE);
[intelh, intelvr, intelinc, intelraan,...
    intelecc, intelarg, intelTA, intelra,...
    intelrp, intela, intelT] =                          RV2COE(intelR, intelV, muE);

[iri1tsp] = TA2t(iri1ecc, iri1TA, iri1h, muE);
    
%% STATE VECTORS

iri1State =     [iri1R;iri1V];
iri2State =     [iri2R;iri2V];
glonasState =   [glonasR;glonasV];
intelState =    [intelR;intelV];

%% ODE45 CALL FOR ORBIT TRAJECTORY

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

iri1tspan =     [0 iri1T];
iri2tspan =     [0 iri2T];
glonastspan =   [0 glonasT];
inteltspan =    [0 intelT];

[iri1TIME, iri1New] =      ode45(@TwoBody, iri1tspan, iri1State, options, muE);
[iri2TIME, iri2New] =      ode45(@TwoBody, iri2tspan, iri2State, options, muE);
[glonasTIME, glonasNew] =  ode45(@TwoBody, glonastspan, glonasState, options, muE);
[intelTIME, intelNew] =    ode45(@TwoBody, inteltspan, intelState, options, muE);

%% POSITION OF SATELLITES AT START

tStart = 181.78619214501 + 365.25; % days

tsp1 = TA2t(iri1ecc, iri1TA, iri1h, muE);
tsp1 = tsp1/(24*3600); %days
tsp2 = TA2t(iri2ecc, iri2TA, iri2h, muE);
tsp2 = tsp2/(24*3600);
tsp3 = TA2t(glonasecc, glonasTA, glonash, muE);
tsp3 = tsp3/(24*3600);
tsp4 = TA2t(intelecc, intelTA, intelh, muE);
tsp4 = tsp4/(24*3600);

dT1 = tStart - (312.77658151 - tsp1); % days
dT1 = dT1*24*3600; %seconds

dT2 = tStart - (312.80323369 - tsp2);
dT2 = dT2 * 24*3600;

dT3 = tStart - (312.07351556 - tsp3);
dT3 = dT3 * 24*3600;

dT4 = tStart - (312.75954111 - tsp4);
dT4 = dT4 * 24*3600;

revs1 = dT1/iri1T;
frac1 = (revs1 - floor(revs1))*iri1T;
TA1atStart = t2TA(frac1, iri1h, muE, iri1ecc);

revs2 = dT2/iri2T;
frac2 = (revs2 - floor(revs2))*iri2T;
TA2atStart = t2TA(frac2, iri2h, muE, iri2ecc);

revs3 = dT3/glonasT;
frac3 = (revs3 - floor(revs3))*glonasT;
TA3atStart = t2TA(frac3, glonash, muE, glonasecc);

revs4 = dT4/intelT;
frac4 = (revs4 - floor(revs4))*intelT;
TA4atStart = t2TA(frac4, intelh, muE, intelecc);

[iri1rX, iri1rY, iri1rZ] = TA2Pos(TA1atStart, iri1arg, iri1inc, iri1raan, iri1h, muE, iri1ecc);
[iri2rX, iri2rY, iri2rZ] = TA2Pos(TA2atStart, iri2arg, iri2inc, iri2raan, iri2h, muE, iri2ecc);
[glonasrX, glonasrY, glonasrZ] = TA2Pos(TA3atStart, glonasarg, glonasinc, glonasraan, glonash, muE, glonasecc);
[intelrX, intelrY, intelrZ] = TA2Pos(TA4atStart, intelarg, intelinc, intelraan, intelh, muE, intelecc);

%% LEO1 TO LEO2 TRANSFER

% coast in LEO1 until apogee
tsincep = TA2t(iri1ecc, iri1TA, iri1h, muE);
tsp1 = tsincep;
t2Peri = iri1T - tsp1;
t2Apo = t2Peri + iri1T/2;

delT1 = t2Apo;  % delT for first phase in transfer (coast time)
delV1 = 0;      % delV for first phase in transfer (no burn)      

% instant burn to circularize with rad = ra of LEO1
va1 = iri1h/iri1ra;
vcirc = sqrt(muE/iri1ra);

delT2 = 0;          % delT for second phase in transfer (instant burn)
delV2 = vcirc-va1;  % delV for second phase in transfer (circularize) 

% instant burn to perform inc+raan change
calpha = cosd(iri1inc)*cosd(iri2inc) + sind(iri1inc)*sind(iri2inc)*cosd(iri2raan-iri1raan);
alpha = acos(calpha);

delT3 = 0;                      % delT for third phase in transfer (instand burn)
delV3 = 2*vcirc*sin(alpha/2);    % delV for third phase in transfer (inc+raan change)

% hohmann tfr to circular orbit with rad = ra of LEO2\
rptfr = iri1ra;
ratfr = iri2ra;
vp = sqrt(muE/iri1ra);
va = sqrt(muE/iri2ra);

ecctfr = (ratfr - rptfr)/(ratfr+rptfr);
htfr = sqrt(rptfr*muE*(1+ecctfr));
atfr = (ratfr+rptfr)/2;
Ttfr = 2*pi*atfr^1.5/sqrt(muE);

vDtfr = htfr/rptfr;
vAtfr = htfr/ratfr;
delV4a = vDtfr - vp;
delV4b = va - vAtfr;

delT4 = Ttfr/2;             % delT for fourth phase in transfer (hohmann period)
delV4 = delV4a + delV4b;    % delV for fourth phase in transfer (hohmann burns)

% coast from peri-side apse line 1 to apo-side apse line 2

TAi = 0;
TAf = 180+(iri2arg - iri1arg);
delTA = TAf-TAi;
hcirc2 = sqrt(iri2ra*muE);

[tToTAf] = TA2t(0, delTA, hcirc2, muE);

delT5 = tToTAf;     % delT for fifth phase in transfer (coast time)
delV5 = 0;          % delV for fifth phase in transfer (no burn)

% burn to decircularize into LEO2 orbit at apogee
va2 = iri2h/iri2ra;
vcirc2 = sqrt(muE/iri2ra);

delT6 = 0;              % delT for sixth phase in transfer (instant burn)
delV6 = vcirc2 - va2;   % delV for sixth phase in transfer (decircularize)

% coast into perigee of LEO2 orbit

delT7 = iri2T/2; % delT for seventh phase in transfer (coast time)
delV7 = 0;       % delV for seventh phase in transfer (no burn)

% calculating TA of object 2 since start

Tpast = delT1+delT2+delT3+delT4+delT5+delT6+delT7;
revs = Tpast/iri2T;
tInNew = (revs - floor(revs))*iri2T;
TANew = t2TA(tInNew, iri2h, muE, iri2ecc);

% choose to rendezvous with object 2 after it completes orbit + 1 extra
% orbit

tRemain = iri2T - tInNew;
tRendez = tRemain + (1*iri2T);

% find phasing orbit such that 1 period = time til rendezvous (tRendez)

rpPhase = iri2rp;
aPhase = (tRendez*sqrt(muE)/(2*pi))^(2/3); %need ra; rpPhase = rpLeo2
raPhase = 2*aPhase - rpPhase;
eccPhase = (raPhase - rpPhase)/(raPhase + rpPhase);
hPhase = sqrt(rpPhase*muE*(1+eccPhase));

% get onto Phase Orbit

vDphase = hPhase/rpPhase;
vi = iri2h/iri2rp;

delT8 = 0;              % delT for eigth phase in transger (instant burn)
delV8 = vDphase - vi;   % delV for eigth phase in transfer (get onto phase orbit)

% stay on Phase Orbit for 1 period

 TPhase = 2*pi*aPhase^1.5 / sqrt(muE);
 
 delT9 = TPhase;    % delT for ninth phase in transfer (coast time) 
 delV9 = 0;         % delV for ninth phase in transfer (no burn)
 
% get off Phase Orbit back on to LEO2 (Rendezvous'd with Object 2!)

vAphase = hPhase/rpPhase;
vf = iri2h/iri2rp;

delT10 = 0;             % delT for tenth phase in transfer (instant burn)
delV10 = vAphase - vf;  % delV for tenth phase in transfer (get off phase orbit)

% stick with object 2 for 5 periods

delT11 = 5*iri2T;   % stay in orbit with object for 5 periods
delV11 = 0;         % no burn; coasting

delT = [delT1 delT2 delT3 delT4 delT5 delT6 delT7 delT8 delT9 delT10 delT11];
DeltaT = sum(delT);

delV = [delV1 delV2 delV3 delV4 delV5 delV6 delV7 delV8 delV9 delV10 delV11];
DeltaV = sum(delV);

% final JD
% JD = 21 182.4137503600938
[iri2FrX, iri2FrY, iri2FrZ] = TA2Pos(0,iri2arg,iri2inc,iri2raan, iri2h, muE, iri2ecc);
[iri2FvX, iri2FvY, iri2FvZ] = TA2Vel(0, iri2arg, iri2inc, iri2raan, iri2h, muE, iri2ecc);

% R, V vectors after rendezvous1
RpostTFR1 = [iri2FrX, iri2FrY, iri2FrZ];
VpostTFR1 = [iri2FvX, iri2FvY, iri2FvZ];

% calc position of Object 3 at JD = 21 182.4137503600938

tsincep = TA2t(glonasecc, glonasTA, glonash, muE);
tperirendez = tsincep +((181.78619214501 + 365.25 - 312.07351556)*24*3600)  + DeltaT;
revs = tperirendez/glonasT;
tInNew = (revs-floor(revs))*glonasT;
TANEW = t2TA(tInNew, glonash, muE, glonasecc);
[glonasrXpost1, glonasrYpost1, glonasrZpost1] = TA2Pos(TANEW,glonasarg,glonasinc,glonasraan, glonash, muE, glonasecc);
[glonasvXpost1, glonasvYpost1, glonasvZpost1] = TA2Vel(TANEW,glonasarg,glonasinc,glonasraan, glonash, muE, glonasecc);

%% LEO1 TO LEO2 TRANSFER ORBIT PLOTTING

% circularize at orbit 1 ra
timespan = [0 2*iri1T];
[X Y Z] = sphere;
X = X*6378;
Y = Y*6378;
Z = Z*6378;
circ1h = vcirc*iri1ra;
[circ1rX, circ1rY, circ1rZ] = TA2Pos(180,iri1arg,iri1inc,iri1raan, circ1h, muE, 0);
[circ1vX, circ1vY, circ1vZ] = TA2Vel(180,iri1arg,iri1inc,iri1raan, circ1h, muE, 0);
statecirc1= [circ1rX;circ1rY;circ1rZ;circ1vX;circ1vY;circ1vZ];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,circrF] = ode45(@TwoBody, timespan, statecirc1, options, muE);

% circular orbit w/ r = orbit2ra
circ2h = vcirc2*iri2ra;
[circ2rX, circ2rY, circ2rZ] = TA2Pos(180,iri2arg,iri2inc,iri2raan, circ2h, muE, 0);
[circ2vX, circ2vY, circ2vZ] = TA2Vel(180,iri2arg,iri2inc,iri2raan, circ2h, muE, 0);
statecirc2= [circ2rX;circ2rY;circ2rZ;circ2vX;circ2vY;circ2vZ];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,circ2rF] = ode45(@TwoBody, timespan, statecirc2, options, muE);

%inc+ raan change into plane of orbit 2
[circ3rX, circ3rY, circ3rZ] = TA2Pos(180,iri2arg,iri2inc,iri2raan, circ1h, muE, 0);
[circ3vX, circ3vY, circ3vZ] = TA2Vel(180,iri2arg,iri2inc,iri2raan, circ1h, muE, 0);
statecirc3= [circ3rX;circ3rY;circ3rZ;circ3vX;circ3vY;circ3vZ];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,circ3rF] = ode45(@TwoBody, timespan, statecirc3, options, muE);

% hohmann from neworbit1 to orbit2
[hohrX, hohrY, hohrZ] = TA2Pos(0,iri2arg,iri2inc,iri2raan, htfr, muE, ecctfr);
[hohvX, hohvY, hohvZ] = TA2Vel(0,iri2arg,iri2inc,iri2raan, htfr, muE, ecctfr);
hoh= [hohrX;hohrY;hohrZ;hohvX;hohvY;hohvZ];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,hohrF] = ode45(@TwoBody, timespan, hoh, options, muE);

% phasing orbit
[phaserX, phaserY, phaserZ] = TA2Pos(0,iri2arg,iri2inc,iri2raan, hPhase, muE, eccPhase);
[phasevX, phasevY, phasevZ] = TA2Vel(0,iri2arg,iri2inc,iri2raan, hPhase, muE, eccPhase);
phase= [phaserX;phaserY;phaserZ;phasevX;phasevY;phasevZ];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,phaserF] = ode45(@TwoBody, timespan, phase, options, muE);

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(iri1New(:,1),iri1New(:,2),iri1New(:,3), 'g', 'linewidth', 4);
hold on
plot3(circrF(:,1),circrF(:,2),circrF(:,3),'c', 'linewidth', 2)
plot3(circ3rF(:,1),circ3rF(:,2),circ3rF(:,3),'r','linewidth',2)
plot3(hohrF(:,1),hohrF(:,2),hohrF(:,3),'m','linewidth',2)
plot3(circ2rF(:,1),circ2rF(:,2),circ2rF(:,3),'k', 'linewidth', 2)
plot3(phaserF(:,1),phaserF(:,2),phaserF(:,3),'g','linewidth',2)
plot3(iri2New(:,1),iri2New(:,2),iri2New(:,3),'m', 'linewidth', 4);
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('Iridium - 33 (LEO-1) Inital Orbit',...
    'Pt1: Circularize at LEO-1 Apogee',...
    'Pt2: Inc + RAAN change into plane of LEO-2 Orbit',...
    'Pt3: Hohmann Transfer Orbit into LEO-2 Orbit',...
    'Pt4: Circuarize about LEO-2 Apogee',...
    'Pt5: Phasing Orbit',...
    'Pt6: Burn into LEO-2 Orbit',...
    'Blue Marble')
xlabel('X-Direction [km]');
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Orbit Trajectories at Start Time: JD 21 181.7861921')

%% LEO2 TO MEO TRANSFER
Rchase2 = RpostTFR1';
Vchase2 = VpostTFR1';
Rtarget3 = [glonasrXpost1;glonasrYpost1;glonasrZpost1];
Vtarget3 = [glonasvXpost1;glonasvYpost1;glonasvZpost1];

[dVTfr2, dTTfr2, Vbounce2] = twoImpulse(Rchase2, Vchase2, Rtarget3, Vtarget3);

timespan = [0 dTTfr2]; %sec

state = [Rtarget3;Vtarget3]; %state vectorR3
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,Rf] = ode45(@TwoBody, timespan, state, options, muE);

dTpers = 5*glonasT;

dTtotal = dTTfr2 + dTpers;

%% MEO TO GEO TRANSFER

% dTTfr2 measured starting at end of TFR1
% Need pos,vel of GLONAS, INTEL at end of TFR + 5 periods

% Pos of MEO post-(tfr2+periods)

timespan = [0 dTtotal]; %sec
state = [Rtarget3;Vtarget3]; %state vectorR3
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,statenew] = ode45(@TwoBody, timespan, state, options, muE);
RGloPreTFR3 =[statenew(end,1);statenew(end,2);statenew(end,3)];
VGloPreTFR3 = [statenew(end,4);statenew(end,5);statenew(end,6)];

% Pos of GEO post-(tfr2 + periods)

dTtotal = dTtotal + DeltaT; %Time since start of mission

[intelvX, intelvY, intelvZ] = TA2Vel(TA4atStart, intelarg, intelinc, intelraan, intelh, muE, intelecc);

Rintel = [intelrX;intelrY;intelrZ];
Vintel = [intelvX;intelvY;intelvZ];

timespan = [0 dTtotal]; %sec
state = [Rintel;Vintel]; %state vectorR3
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,statenew] = ode45(@TwoBody, timespan, state, options, muE);
RIntelPreTFR3 =[statenew(end,1);statenew(end,2);statenew(end,3)];
VIntelPreTFR3 = [statenew(end,4);statenew(end,5);statenew(end,6)];

% function
Rchase = RGloPreTFR3;
Vchase = VGloPreTFR3;
Rtarget = RIntelPreTFR3;
Vtarget = VIntelPreTFR3;
[dVTfr3, dTTfr3, Vbounce3] = twoImpulse(Rchase, Vchase, Rtarget, Vtarget);

timespan = [0 dTTfr3]; %sec
state = [Rtarget;Vtarget]; %state vectorR3
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tnew,Rf4] = ode45(@TwoBody, timespan, state, options, muE);

%% FINAL NUMBERS

DVFINAL = DeltaV+dVTfr2+dVTfr3;
DTFINAL = DeltaT+dTTfr2+(5*glonasT)+dTTfr3+(5*intelT);

JDFINAL = 181.7861921 + (DTFINAL/(24*3600));

%% MISSION DEBRIEF WITH TWO IMPULSE
fprintf('*******************************************\n\n')
fprintf('END OF MISSION DEBRIEF WITH 2 IMPULSE \n\n')
fprintf('TOTAL DELTA-V REQUIRED [km/s]: %f\n\n',DVFINAL);
fprintf('DELTA-V REQUIRED FOR LEO1 TO LEO2 [km/s]: %f\n',DeltaV);
fprintf('DELTA-V REQUIRED FOR LEO2 TO MEO [km/s]: %f\n',dVTfr2);
fprintf('DELTA-V REQUIRED FOR MEO TO GEO [km/s]: %f\n\n',dVTfr3);
fprintf('TOTAL TIME REQUIRED [days]: %f\n\n',DTFINAL/(24*3600));
fprintf('DELTA-T REQUIRED FOR LEO1 TO LEO2 TFR [days]: %f\n',DeltaT/(24*3600));
fprintf('DELTA-T REQUIRED FOR 5 LEO2 PERIODS [days]: %f\n',(5*iri2T)/(24*3600));
fprintf('DELTA-T REQUIRED FOR LEO2 TO MEO TFR [days]: %f\n',dTTfr2/(24*3600));
fprintf('DELTA-T REQUIRED FOR 5 MEO PERIODS [days]: %f\n',(5*glonasT)/(24*3600));
fprintf('DELTA-T REQUIRED FOR MEO TO GEO TFR [days]: %f\n',dTTfr3/(24*3600));
fprintf('DELTA-T REQUIRED FOR 5 GEO PERIODS [days]: %f\n\n',(5*intelT)/(24*3600));
fprintf('INITAL JD WHILE ON LEO-1 OBJECT: JD 21 181.7861921\n')
fprintf('FINAL JD AFTER 5 GEO PERIODS: JD 21 %f\n',JDFINAL)
fprintf('\n*******************************************\n')

%% FIGURE 1: ALL ORBITS AND SATS AT START TIME
[X,Y,Z] = sphere;
X = X*6378;
Y = Y*6378;
Z = Z*6378;

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(iri1New(:,1),iri1New(:,2),iri1New(:,3), 'g', 'linewidth', 2);
hold on
plot3(iri1rX, iri1rY, iri1rZ, '-o','markeredgecolor','k','markerfacecolor','g', 'markersize', 10)
plot3(iri2New(:,1),iri2New(:,2),iri2New(:,3),'m', 'linewidth', 2);
plot3(iri2rX, iri2rY, iri2rZ, '-o','markeredgecolor','k','markerfacecolor','m', 'markersize', 10)
plot3(glonasNew(:,1),glonasNew(:,2),glonasNew(:,3),'b', 'linewidth', 2);
plot3(glonasrX, glonasrY, glonasrZ, '-o','markeredgecolor','k','markerfacecolor','b', 'markersize', 10)
plot3(intelNew(:,1),intelNew(:,2),intelNew(:,3),'r', 'linewidth', 2);
plot3(intelrX, intelrY, intelrZ, '-o','markeredgecolor','k','markerfacecolor','r', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('Iridium - 33 (LEO-1) Trajectory', 'Iridium - 33 (LEO-1) Satellite',...
    'Iridium - 33 (LEO-2) Trajectory', 'Iridium - 33 (LEO-2) Satellite',...
    'GLONAS (MEO) Trajectory', 'GLONAS (MEO) Satellite',...
    'IntelSAT (GEO) Trajectory', 'IntelSAT (GEO) Satellite',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Orbit Trajectories at Start Time: JD 21 181.7861921');

%% FIGURE 2: LEO1 AND LEO2 

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(iri1New(:,1),iri1New(:,2),iri1New(:,3), 'g', 'linewidth', 2);
hold on
plot3(iri1rX, iri1rY, iri1rZ, '-o','markeredgecolor','k','markerfacecolor','g', 'markersize', 10)
plot3(iri2New(:,1),iri2New(:,2),iri2New(:,3),'m', 'linewidth', 2);
plot3(iri2rX, iri2rY, iri2rZ, '-o','markeredgecolor','k','markerfacecolor','m', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('Iridium - 33 (LEO-1) Trajectory','Iridium - 33 (LEO-1) Satellite','Iridium - 33 (LEO-2) Trajectory', 'Iridium - 33 (LEO-2) Satellite','Blue Marble')
xlabel('X-Direction [km]');
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Orbit Trajectories at Start Time: JD 21 181.7861921')

%% FIGURE 3: RENDEZVOUS AT LEO 2
figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(iri2New(:,1),iri2New(:,2),iri2New(:,3),'m', 'linewidth', 2);
hold on
plot3(iri2FrX, iri2FrY, iri2FrZ, '-o','markeredgecolor','k','markerfacecolor','m', 'markersize', 10)
plot3(iri2rX, iri2rY, iri2rZ,'-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('Iridium - 33 (LEO-2) Trajectory', 'Iridium - 33 (LEO-2) Satellite',...
    'Iridium - 33 (LEO-2) at Initial Time',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Iridium - 33 Satellite Post-Rendezvous/Transfer 1')

%% FIGURE 4: LEO2 AND MEO 

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(iri2New(:,1),iri2New(:,2),iri2New(:,3),'m', 'linewidth', 2);
hold on
plot3(iri2FrX, iri2FrY, iri2FrZ, '-o','markeredgecolor','k','markerfacecolor','m', 'markersize', 10)
plot3(glonasNew(:,1),glonasNew(:,2),glonasNew(:,3),'b', 'linewidth', 2);
plot3(glonasrXpost1, glonasrYpost1, glonasrZpost1, '-o','markeredgecolor','k','markerfacecolor','b', 'markersize', 10)
plot3(Rf(end,1),Rf(end,2),Rf(end,3), '-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
plot3(iri2rX, iri2rY, iri2rZ,'-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
plot3(glonasrX, glonasrY, glonasrZ, '-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('Iridium - 33 (LEO-2) Trajectory', 'Iridium - 33 (LEO-2) Satellite',...
    'GLONAS (MEO) Trajectory', 'GLONAS (MEO) Satellite',...
    'GLONAS (MEO) Position at Rendezvous',...
    'Iridium - 33 (LEO-2) at Initial Time',...
    'GLONAS (MEO) at Initial Time',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Orbit Trajectories at Start Time: JD 21 182.4137504');

%% FIGURE 5: RENDEZVOUS AT MEO

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(glonasNew(:,1),glonasNew(:,2),glonasNew(:,3),'b', 'linewidth', 2);
hold on
plot3(Rf(end,1), Rf(end,2), Rf(end,3), '-o','markeredgecolor','k','markerfacecolor','b', 'markersize', 10)
plot3(glonasrX, glonasrY, glonasrZ, '-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('GLONAS (MEO) Trajectory', 'GLONAS (MEO) Satellite',...
    'GLONAS (MEO) Satellite at Initial Time',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('GLONAS Rocket Body Post-Rendezvous/Transfer 2');

%% FIGURE 6: MEO AND GEO

figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(glonasNew(:,1),glonasNew(:,2),glonasNew(:,3),'b', 'linewidth', 2);
hold on
plot3(RGloPreTFR3(1), RGloPreTFR3(2), RGloPreTFR3(3), '-o','markeredgecolor','k','markerfacecolor','b', 'markersize', 10)
plot3(intelNew(:,1),intelNew(:,2),intelNew(:,3),'r', 'linewidth', 2);
plot3(RIntelPreTFR3(1), RIntelPreTFR3(2), RIntelPreTFR3(3),'-o','markeredgecolor','k','markerfacecolor','r', 'markersize', 10)
plot3(Rf4(end,1),Rf4(end,2),Rf4(end,3),'-o','markeredgecolor','k','markerfacecolor','k', 'markersize', 10)
plot3(glonasrX, glonasrY, glonasrZ, '-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
plot3(intelrX, intelrY, intelrZ, '-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('GLONAS (MEO) Trajectory', 'GLONAS (MEO) Satellite',...
    'IntelSAT (GEO) Trajectory', 'IntelSAT (GEO) Satellite',...
    'IntelSat (GEO) Position at Rendezvous',...
    'GLONAS (MEO) Satellite at Initial Time',...
    'IntelSAT (GEO) Satellite at Initial Time',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('Orbit Trajectories at Start Time: JD 21 185.1668090');

%% FIGURE 7: RENDEZVOUS AT GEO
figure('units','normalized','outerposition',[0.25 0.25 0.75 0.75])
plot3(intelNew(:,1),intelNew(:,2),intelNew(:,3),'r', 'linewidth', 2);
hold on
plot3(Rf4(end,1), Rf4(end,2), Rf4(end,3),'-o','markeredgecolor','k','markerfacecolor','r', 'markersize', 10)
plot3(intelrX, intelrY, intelrZ, '-o','markeredgecolor','k','markerfacecolor','w', 'markersize', 10)
surf(X,Y,Z);
alpha 0.1;
axis equal;
hold off

legend('IntelSAT (GEO) Trajectory', 'IntelSAT (GEO) Satellite',...
    'IntelSAT (GEO) Satellite at Initial Time',...
    'Blue Marble')
xlabel('X-Direction [km]'); 
ylabel('Y-Direction [km]');
zlabel('Z-Direction [km]');
title('IntelSat Post-Rendezvous/Transfer 3');

%% functions
function [R,V] = tle2RV(CP7, muE)

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

E_1 = mod(E_1, 2*pi);

TA = 2*atand((sqrt((1+ecc)/(1-ecc)) * tan(E_1/2)));

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

function [h,vr,inc,raan,ecc,arg,TA,ra,rp,a,T] = RV2COE(R,V, muE)

hbar = cross(R,V);
h = norm(hbar);
vr = dot(R,V)/norm(R);

inc = acosd(hbar(3)/h);

N = cross([0 0 1], hbar);
raan = acosd(N(1)/norm(N));

if N(2) < 0
    raan = 360-raan;
end

eccbar = cross(V,hbar)/muE - R/norm(R);
ecc = norm(eccbar);

arg = acosd((dot(N,eccbar))/(norm(N)*ecc));

TA = acosd((dot(eccbar, R))/(norm(R)*ecc));
if vr < 0
    TA = 360-TA;
end

ra = h^2/muE * (1/(1-ecc));
rp = h^2/muE * (1/(1+ecc));
a = (ra+rp)/2;

T = 2*pi*a^1.5 / sqrt(muE);

end

function [rX, rY, rZ] = TA2Pos(TA,arg,inc,raan, h, muE, ecc)

Q1 = [cosd(arg) sind(arg)   0;...
    -sind(arg)  cosd(arg)   0;...
    0           0           1];

Q2 = [1 0           0;...
    0   cosd(inc)   sind(inc);...
    0   -sind(inc)  cosd(inc)];

Q3 = [cosd(raan)    sind(raan)  0;...
    -sind(raan)     cosd(raan)  0;...
    0               0           1];

Q = Q1*Q2*Q3;

Rmatr = h^2/(muE*(1+ecc*cosd(TA))) * [cosd(TA);sind(TA);0];

R = Q'*Rmatr;
rX = R(1);
rY = R(2);
rZ = R(3);


end

function [vX, vY, vZ] = TA2Vel(TA, arg, inc, raan, h, muE, ecc)

Q1 = [cosd(arg) sind(arg)   0;...
    -sind(arg)  cosd(arg)   0;...
    0           0           1];

Q2 = [1 0           0;...
    0   cosd(inc)   sind(inc);...
    0   -sind(inc)  cosd(inc)];

Q3 = [cosd(raan)    sind(raan)  0;...
    -sind(raan)     cosd(raan)  0;...
    0               0           1];

Q = Q1*Q2*Q3;

Vmatr = muE/h * [-sind(TA); ecc + cosd(TA); 0];
V = Q'*Vmatr;

vX = V(1);
vY = V(2);
vZ = V(3);

end

function [tsp] = TA2t(ecc, TA, h, muE)

E = 2*atan((sqrt((1-ecc)/(1+ecc)) * tand(TA/2)));

E = mod(E,2*pi);

Me = E - ecc*sin(E);

tsp = Me*h^3 / (muE^2 * (1-ecc^2)^1.5);

end

function [TA] = t2TA(t, h, muE, ecc)

Me = muE^2 * (1-ecc^2)^1.5 * t /h^3;

f = @(E) Me - E+ecc*sin(E);
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

E_1 = mod(E_1, 2*pi);

inner = tan(E_1/2);
inner = sqrt((1+ecc)/(1-ecc)) * inner;

TA = 2*atand(inner);

TA = mod(TA, 360);

end

function [V1, V2] = Lambert(R1, R2, delT, tol, muE)

r1 = norm(R1);
r2 = norm(R2);

traj = cross(R1,R2);
delTA = acosd((dot(R1,R2))/(r1*r2));

if traj(3) < 0
    delTa = 360-delTA;
end

A = sind(delTA)*sqrt((r1*r2)/(1-cosd(delTA)));
Z_0 = 1;

C = @(z) (1-cos(sqrt(z)))/z;
S = @(z) (sqrt(z) - sin(sqrt(z)))/((sqrt(z))^3);

y = @(z) r1 + r2 + (A*z*S(z)-A)/(sqrt(C(z)));
T = @(z) A*sqrt(y(z)/muE);

F = @(z) ((y(z)/C(z))^1.5)*S(z) + A*sqrt(y(z)) - sqrt(muE)*delT;
Fp = @(z) (y(z)/C(z))^1.5 * (1/(2*z)*(C(z) - (3*S(z)/(2*C(z)))) + (3*(S(z))^2 / (4* ((C(z)))))) + A/8 * (3*S(z)/C(z) * sqrt(y(z)) + A*(sqrt(C(z)/y(z))));

Z_1 = Z_0 - (F(Z_0)/Fp(Z_0));
err = abs(Z_1 - Z_0);

while err > tol
    Z_0 = Z_1;
    Z_1 = Z_0 - (F(Z_0)/Fp(Z_0));

    err = abs(Z_1 - Z_0);
end

f = @(z) 1-y(z)/r1;
fdot = @(z) sqrt(muE)/(r1*r2) * sqrt(y(z)/C(z))*(z*S(z) - 1);
g = @(z) A*sqrt(y(z)/muE);
gdot = @(z) 1- y(z)/r2;

V1 = 1/(g(Z_1)) * (R2 - (f(Z_1)*R1));
V2 = 1/(g(Z_1)) * (gdot(Z_1)*R2 - R1);

end

function [dVTfr, dTTfr, Vbounce] = twoImpulse(Rchase, Vchase, Rtarget, Vtarget)

ihat = Rtarget/norm(Rtarget);
jhat = Vtarget/norm(Vtarget);
khat = cross(ihat,jhat);

Q = [ihat';jhat';khat'];

delR = Rchase - Rtarget;

nTarget = norm(Vtarget)/norm(Rtarget);
raantarget =  nTarget*khat;

delV = Vchase - Vtarget - cross(raantarget, delR);

delR0 = Q*delR;
delV0min = Q*delV;

dv = zeros(500,5);

for i = 1:500

    dv(i,1) = i;
    
    t = i*3600;

Phirr = [4-3*cos(nTarget*t) 0 0; 6*(sin(nTarget*t) - nTarget*t) 1 0; 0 0 cos(nTarget*t)];

Phirv = [sin(nTarget*t)/nTarget 2*(1-cos(nTarget*t))/nTarget 0; 2*(cos(nTarget*t) - 1)/nTarget (4*sin(nTarget*t) - 3*nTarget*t)/nTarget 0; 0 0 sin(nTarget*t)/nTarget];

Phivr = [3*nTarget *sin(nTarget*t) 0 0; 6*nTarget*(cos(nTarget*t) -1) 0 0; 0 0 -nTarget*sin(nTarget*t)];

Phivv = [cos(nTarget*t) 2*sin(nTarget*t) 0; -2*sin(nTarget*t) -4*cos(nTarget*t)-3 0; 0 0 cos(nTarget*t)];

delvI = -inv(Phirv)*Phirr*delR0;

delv0plus = Phivr*delR0 + Phivv*delvI;

delvFmin = Phivr*delR0 + Phivv*delv0plus;

deltavI = delv0plus - delV0min;
deltavF = [0;0;0] - delvFmin;

DELTAV = norm(deltavI) + norm(deltavF);

dv(i,2) = deltavI(1);
dv(i,3) = deltavI(2);
dv(i,4) = deltavI(3);
dv(i,5) = DELTAV;
end

[val,idx] = min(dv(:,5));% t = 59hour

dVTfr = val; % km/s
dTTfr = idx*3600; % sec
delVF = [dv(idx,2);dv(idx,3);dv(idx,4)];
Vbounce = Vchase + delVF; 

end