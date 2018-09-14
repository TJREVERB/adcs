% Thomas Baron and Justin Zhou
% Updated 11/21/2017 by S. Hur-Diaz, Emergent
% CubeSat Detumble Simulation Code
% Credit to Mirko Leomanni, [Add citation] 
% Assumes spacecraft is in a circular orbit and is at the equatoral crossing
% at time zero

% SIMULATION PARAMETERS

%% SPACECRAFT
% Original Spacecraft properties of the paper
% sc.inertia=diag([0.33,0.37,0.35]); %Spacecraft inertia (kg*m^2)
% sc.rc=7021; %Orbit radius (km) (LEO)
% sc.In=65*pi/180; %Inclination (rad)
% sc.Torb=5855; %Orbital period (s)

% TJ REVERB's estimated properties
sc.inertia = diag([.213 .213 .213]); % Spacecraft inertia matrix (kg*m^2), based on uniform density 
sc.rc = 6771; % Orbit radius (km) (LEO), based on ISS average altitude
sc.In = 51.6*pi/180; % Orbit inclination (rad), based on ISS inclination
sc.Torb = 2*pi/sqrt(3.986e5/sc.rc^3); % Orbital period (s), based on calculated orbital period from ISS altitude, note decay 

sc.wo = 2*pi/sc.Torb; % Angular orbital speed (rad/s)

%% Magtorquer  properties
ct.mmax=.2; % Maximum magnetic moment (A*m^2), based on proposed mag spec sheet 
% ct.mmax=2; % Maximum magnetic moment (A*m^2), paper 
% TODO: where is this located on the spacecraft?

%% Magnetic field/environment properties
mg.Mt = 7.8379e6; % Magnetic moment Earth T*km^3
mg.Gta = 11.44*pi/180; % Geomagnetic tilt angle (rad)
mg.kmg = mg.Mt/sc.rc^3; % Dipole Magnitude (T)
mg.we = 7.2921150e-5; % Earth rotation speed (rad/s)
mg.MeanB2 = 1.0131e-009; 
% mg.Beta1 = 0;

%% Simulation parameters
Ts = 40000; % Simulation time (sec)
tcamp = 10; % Plot sampling time
tint = .1; % Fixed integration step in RK4

%% GAIN CALCULATION of the torquer bar controller 
% Calculate the inclination of the orbital plane with respect to the
% geomagnetic equator (epsm)  
[~, sinepsm] = magField(0,sc,mg);

% Gain formula based on reference #3
% Can try other values to see if performance is improved
% If the gain is too high, then you start losing finer control
ct.kw = 2*sc.wo*(1+sinepsm)*min(diag(sc.inertia));

%% INITIAL CONDITION
q0 = [-0.062;0.925;-0.007;0.375];  % initial quaternion [vec;sca]
w0 = [.604;-1.760;-.384]; % initial angular vector (rad/s)   
w0 = [.604;-.760;-.384]; % initial angular vector (rad/s)   

% Preallocating matrix with initial condition and zeroes 
Yt=[q0' w0';zeros(round(Ts/tcamp),7)];

% Time vector for simulation
tt = [0:tcamp:Ts];

for Js=1:round(Ts/tcamp)

    % SIMULATION WITH FIXED INTEGRATON STEP (tint)

    Ys=ode4(@(t,Y)AttDyn(t,Y,sc,mg,ct),tt(Js):tint:(tt(Js)+tcamp),Yt(Js,:)');
    Yt(Js+1,:)=Ys(end,:);
    Yt(Js+1,:) = Yt(Js+1,:)/norm(Yt(Js+1,:));

end

%% GET INPUT FOR PLOTTING (Only for angular speed feedback)

usim = zeros(3,length(tt));
bv = zeros(3,length(tt));
att = zeros(3,length(tt));

for i=1:length(tt)

    [~,bv(:,i),usim(:,i)] = AttDyn(tt(i),Yt(i,:)',sc,mg,ct);
    % Transform inertial attitude to orbital attitude
    % TBI = TBO*TOI ==> TBI*(TOI)^-1 = TBO
    % If we have inertial position and velocity from GPS, we can use a different formula
    q = qmult(Yt(i,1:4)',qinv(qmult(getq(1,-pi/2),qmult(getq(3,sc.wo*tt(i)+pi/2),getq(1,sc.In)))));
    att(:,i) = 2*q(1:3);
    
end

% Plotting
% Actuator Graphs
r2d = 180/pi;
figure
subplot(3,3,1)
plot(tt/60, usim(1,:))
title("Torque Bar 1(A*m^2)")

subplot(3,3,2)
plot(tt/60, usim(2,:))
title("Torque Bar 2 (A*m^2)")

subplot(3,3,3)
plot(tt/60, usim(3,:))
title("Torque Bar 3 (A*m^2)")

% Angular Velocity Graphs
subplot(3,3,4)
plot(tt/60, Yt(:,5)*r2d)
title("Angular Velocity 1 (deg/s)")

subplot(3,3,5)
plot(tt/60, Yt(:,6)*r2d)
title("Angular Velocity 2 (deg/s)")

subplot(3,3,6)
plot(tt/60, Yt(:,7)*r2d)
title("Angular Velocity 3 (deg/s)")

% Magnetic field graphs
subplot(3,3,7)
plot(tt/60,bv(1,:))
title("Magnetic field in body-x (T)")
xlabel('min')

subplot(3,3,8)
plot(tt/60,bv(2,:))
title("Magnetic field in body-y (T)")
xlabel('min')

subplot(3,3,9)
plot(tt/60,bv(3,:))
title("Magnetic field in body-z (T)")
xlabel('min')

% Plot attitude
figure
subplot(3,1,1)
plot(tt/60,att(1,:)*r2d)
title('Roll (deg)')
xlabel('min')

subplot(3,1,2)
plot(tt/60,att(2,:)*r2d)
title('Pitch (deg)')
xlabel('min')

subplot(3,1,3)
plot(tt/60,att(3,:)*r2d)
title('Yaw (deg)')
xlabel('min')


