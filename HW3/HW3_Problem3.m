%comparing CW and nonlinear
clear all
clc
close all

%Constants%%%
G = 6.6742*10^(-20); %gravitational constant [km^3/(kg*s^2)]
m1 = 5.972*10^24; %mass of Earth [kg]
mu = G*m1; %gravitational parameter [km^3/s^2]

R1 = 6786; %[km] radius of ISS is semimajor axis b/c circular
n1 = sqrt(mu/(R1^3)); %[rad/s] angular velocity of ISS
T_Iss = 2*pi/n1; %[sec] orbital periond of ISS

%Given transfer time: half of ISS period
t12 = 1/2*T_Iss;

global theta_dot_c
%%%%%%%%%%%%

%%%%%%%%%%%%%%CHIEF%%%%%%%%%%%%%%%%%%%
%Chief orbital elements at initial time
a=R1; %semi-major axis [kilometers] (circular so a=r)
enorm=0; %eccentricity (circular so 0)
inclination=51.6; %inclination [degrees] of ISS and other sattelite bc coplanar
RAAN=0; %  right ascension of the ascending node [degrees]
arg_per=0; % argument of perigee [degrees] circular
true_anomaly=0; %True anomaly [degrees] Consider zero for optimal transfer angle later

%converting orbital parameters to ECI cartsian coordinates 
r=zeros(3,1); %fill r vector with zeroes
v=zeros(3,1); %fill v vector with zeroes

p=a*(1-enorm^2); %intermediate variable
q=p/(1+enorm*cosd(true_anomaly));%intermediate variable

% Creating r vector in pqw coordinates
R_pqw(1,1) = q*cosd(true_anomaly);
R_pqw(2,1) = q*sind(true_anomaly);
R_pqw(3,1) = 0;
    
% Creating v vector in pqw coordinates    
V_pqw(1,1) =-(mu/p)^.5*sind(true_anomaly);
V_pqw(2,1) =((mu/p)^.5)*(enorm + cosd(true_anomaly));
V_pqw(3,1) =   0;

%Solving for 313 rotation matrices
R1_i = [1 0 0; 0 cosd(inclination) sind(inclination); 0 -sind(inclination) cosd(inclination)];
R3_Om = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
R3_om = [cosd(arg_per) sind(arg_per) 0; -sind(arg_per) cosd(arg_per) 0; 0 0 1];

support_var = R3_om*R1_i*R3_Om; %Intermediate variable
support_var=support_var'; %Transposed

r_c(:,1)=support_var*R_pqw; %Radius r [km] in ECI Cartesian
v_c(:,1)=support_var*V_pqw; %Velocity v [km/s] in ECI Cartesian
theta_dot_c = sqrt(mu/norm(r_c)^3);
%%%%%%%%%%%%%%CHIEF%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%DEPUTY%%%%%%%%%%%%%%%%%%%
%Deputy must be 
R2 = 12000; %[km]radius of outer sattelite
n2 = sqrt(mu/(R2^3)); %[rad/sec]radial velocity of outer sattelite
T_sat = 2*pi/n2; %[sec]period of outer sattelite

%Can only leave for Hohmann transfer when THIS is satisfied
PHI_0 = pi - n2*t12; %[rad]

%Deputy orbital elements at initial time
a=R2; %semi-major axis [kilometers] 
enorm=0; %eccentricity. 0
inclination=51.6; %inclination [degrees]
RAAN=0; %  right ascension of the ascending node [degrees]
arg_per=0; % argument of perigee [degrees]
true_anomaly=0 + PHI_0; %True anomaly [degrees], add optimal transfer to initial condition, which I chose as 0


%converting orbital parameters to ECI cartsian coordinates 
r=zeros(3,1); %fill r vector with zeroes
v=zeros(3,1); %fill v vector with zeroes

p=a*(1-enorm^2); %intermediate variable
q=p/(1+enorm*cosd(true_anomaly));%intermediate variable

% Creating r vector in pqw coordinates
R_pqw(1,1) = q*cosd(true_anomaly);
R_pqw(2,1) = q*sind(true_anomaly);
R_pqw(3,1) = 0;
    
% Creating v vector in pqw coordinates    
V_pqw(1,1) =-(mu/p)^.5*sind(true_anomaly);
V_pqw(2,1) =((mu/p)^.5)*(enorm + cosd(true_anomaly));
V_pqw(3,1) =   0;

%Solving for 313 rotation matrices
R1_i = [1 0 0; 0 cosd(inclination) sind(inclination); 0 -sind(inclination) cosd(inclination)];
R3_Om = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
R3_om = [cosd(arg_per) sind(arg_per) 0; -sind(arg_per) cosd(arg_per) 0; 0 0 1];

support_var = R3_om*R1_i*R3_Om; %Intermediate variable
support_var=support_var'; %Transposed

r_d(:,1)=support_var*R_pqw; %Radius r [km] in ECI Cartesian
v_d(:,1)=support_var*V_pqw; %Velocity v [km/s] in ECI Cartesian
%%%%%%%%%%%%%%DEPUTY%%%%%%%%%%%%%%%%%%%

%CHIEF numerical integration of Cartesian State Vector
Y0 = [r_c; v_c];
tf = 10*3600; %s 
TSPAN = [0 tf]; %duration of integration in seconds
OPTIONS = odeset('Maxstep', 1);
[TOUT_c,YOUT_c] = ode45(@relacc,TSPAN,Y0,OPTIONS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CHIEF numerical integration of Cartesian State Vector
Y0 = [r_d; v_d];
tf = 10*3600; %s 
TSPAN = [0 tf]; %duration of integration in seconds
OPTIONS = odeset('Maxstep', 1);
[TOUT_d,YOUT_d] = ode45(@relacc,TSPAN,Y0,OPTIONS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b] = size(TOUT_d);

%converting into relative position and velocity in LVLH
for i = 1:a
    i_lvlh = YOUT_c(i, 1:3)/norm(YOUT_c(i, 1:3));
    j_lvlh = YOUT_c(i, 4:6)/norm(YOUT_c(i, 4:6));
    k_lvlh = cross(i_lvlh, j_lvlh);
    rho(1:3, i) = [i_lvlh; j_lvlh; k_lvlh]*(YOUT_d(i, 1:3)'-YOUT_c(i, 1:3)');
    vrel(1:3, i) = [i_lvlh; j_lvlh; k_lvlh]*(YOUT_d(i, 4:6)'-YOUT_c(i, 4:6)') - cross(theta_dot_c*[0; 0; 1], rho(1:3, i));
end

figure(1)
%plot(rho(2,:), rho(1,:))
grid on
hold on
plot(rho(2,1), rho(1,1),'r*')
plot(rho(2,end), rho(1,end),'*')
xlabel('along-track y [km]')
ylabel('radial x [km]')

%integrating CW eqq using initial conditions in LVLH computed above
TSPAN = [0 tf]; %duration of integration in seconds
X0 = [rho(1:3,1);  vrel(1:3, 1)];
OPTIONS = odeset('Maxstep', 1);
[TOUT,YOUT] = ode45(@CWeqq_HW3,TSPAN,X0,OPTIONS);

figure(1)
plot(YOUT(:,2), YOUT(:,1), 'r')
plot(YOUT(1,2), YOUT(1,1), 'k*')

X0

%% Part B: Find Approach Velocity 

t = t12;
n = n1;
PHIrr = [(4-3*cos(n*t)) (0) (0); (6*(sin(n*t) - n*t)) (1) (0); (0) (0) (cos(n*t))];
PHIrv = [((1/n)*sin(n*t)) ((2/n)*(1-cos(n*t))) (0); ((2/n)*(cos(n*t)-1)) ((1/n)*(4*sin(n*t)-3*n*t)) (0); (0) (0) ((1/n)*sin(n*t))];
PHIvr = [(3*n*sin(n*t)) (0) (0); (6*n*(cos(n*t)-1)) (0) (0); (0) (0) (-n*sin(n*t))];
PHIvv = [(cos(n*t)) (2*sin(n*t)) (0); (-2*sin(n*t)) (4*cos(n*t)-3) (0); (0) (0) (cos(n*t))];

%Set initial x y and z to zero, solve using matrix inverse
ini_pos = [X0(1); X0(2); X0(3)];
ini_vel = [X0(4); X0(5); X0(6)];
invPHIrv = inv(PHIrv); %inverse for ezpz use
ini_pos %[km] matrix to show initial position
ini_vel %[km/s] matrix to show initial velocity
approach_vel = -invPHIrv*PHIrr*ini_pos %[km/sec] matrix of approach velocity
