% Carlos Ortiz, EAS 4510 Homework 1
%5 Feb, 2021

clc;clear;close;

%%Define variables and initialize x-y plot for later
mu_e = 398600; %[km^2/s^2], gravitational constant of earth

Re = 6378;      %[km], radius of earth
h_min = 8622;   %[km], minimum height above earth
h_max = 40622;   %[km], maxumum height above earth

rp = Re + h_min; %[km], perigee
ra = Re + h_max; %[km], apogee

figure(1), hold on, axis equal, grid on
xlabel('ECI X - Position [km]')
ylabel('ECI Y - Position [km]')
t_plot = linspace(0,2*pi,100);            %[rad], values of theta to plot earth
plot(Re*cos(t_plot) , Re*sin(t_plot),'b', 'linewidth',2) %plot a circle to represent earth on x y plane


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Question 1%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A) Eccentricity, e
e = (ra - rp)/(ra + rp); %[-]

%B)Semi-Major Axis, a 
a = ra/(1+e); %[km], also works with a = 1/2(ra+rp)

%C)Orbital Period, T
T = (2*pi)/(mu_e)^(1/2)*a^(3/2); %[sec]

%D)Velocity at Perigee, Vp
Vp = sqrt(mu_e*(2/rp-1/a)); %[km/sec]

%E)Velocity at Apogee, Va
Va = sqrt(mu_e*(2/ra-1/a)); %[km/sec]

%F)Specific Angular Momentum at Perigee, hp
hp = rp*Vp; %[km^2/sec]

%G)Specific Angular Momentum at Apogee, ha
ha = ra*Va; %[km^2/sec]
%They are equal becuase in keplerain orbits, h is constant
%h = hp = ha = const

%Display Solution:
fprintf("Question 1:\nA) e = %g\nB) a = %d [km]\nC) T = %e [s]\nD) Vp = %f [km/s]\nE) Va = %f [km/s]\nF) hp = %e [km^2/s]\nG) ha = %e [km^2/s]\n",e,a,T,Vp,Va,hp,ha)
fprintf("\nVp is much faster than Va. It is %g times faster.\n", Vp/Va)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Question 2%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use ODE45 to integrate and plot trajectory
r = [rp; 0; 0]; %[km], initial position vector
v = [0; Vp; 0]; %[km/s], initial velocity vector
s = [r; v];     %[km;km/s], state vector

tspan = 0:T; %[s], time span for one orbital period
OPTIONS = odeset('Maxstep',10); %set a minimum time span of 10 sec
[~,S] = ode45(@Ortiz_HW1_2021_P2_EOM, tspan, s, OPTIONS);%integrate Equations of Motion for Problem 2
%S contains position and velocity for every second

%Split S matrix into its components
Xa = S(:,1)';%[km] - ECI X coordinate
Ya = S(:,2)';%[km] - ECI Y coordinate
Za = S(:,3)';%[km] - ECI Z coordinate

%Plot results
plot(Xa, Ya, 'r', 'linewidth', 2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Question 3%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True anomoly approach

%First define angular momentum vector (establishes plane of orbit)
%get unit vector in direction of h and of e, find the vector in J with
%cross product >>>> Establish basis in true anomoly reference frame

h_vec = cross(r,v); %[km/s], specific angular momentum vector
e_vec = ((Vp^2-mu_e/rp)*r-(r'*v)*v)/mu_e; %eccentricity Vector

%Establish unit vector basis
z_hat = h_vec / norm(h_vec); % unit vector normal to the plane of rotation
x_hat = e_vec / norm(e_vec); % unit vector in the direction of e
y_hat = cross(z_hat, x_hat); % unit vector in y direction

theta = (pi/180 : pi/180 : 2*pi);
R = (norm(h_vec)^2/mu_e) * (1./(1+norm(e_vec)*cos(theta))); % satellite position based on true anomaly
 
r = zeros(3, length(theta)); % create a matrix to fill with values of position
for i = 1:length(theta)
    r(:,i) = R(i)*cos(theta(i))*x_hat + R(i)*sin(theta(i))*y_hat;
end

%Plot Results
plot(r(1,:),r(2,:), 'k', 'linewidth', 2)
legend('Earth','Problem 2 Trajectory','Problem 3 Trajectory')
fprintf("\nQuestion 3:\nNote:\nThe trajectory from problem 3 matches that of problem 2.\n")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Question 4%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part A done using ODE45 to integrate

figure(2), hold on, axis equal, grid on
xlabel('ECI X - Position [km]')
ylabel('ECI Y - Position [km]')
t_plot = linspace(0,2*pi,100);%[rad], values of theta to plot earth
plot(Re*cos(t_plot) , Re*sin(t_plot),'b', 'linewidth',2) %plot a circle to represent earth on x y plane

%Note: Given Semilatus rectum, find value of a from: p = a(1-e^2)
%a) p = 19000 e = 0.6
%Find semimajor axis, a from given information
a4_a = 19000/(1-0.6^2);
%Find radius at periapse for r and v
rp4a = a4_a*(1-0.6);

r = [rp4a;0;0];           %[km], position vector
v = [0;sqrt(2*mu_e*(1/rp4a-1/(2*a4_a)));0]; %[km/s], velocity vector
s = [r; v];            %state vector, [km;km/s]

[~,S] = ode45(@Ortiz_HW1_2021_P2_EOM, tspan, s, OPTIONS);%integrate Equations of Motion for Problem 2
%S contains position and velocity for every second

%Split S matrix into its components
X4a = S(:,1)';%[km] - ECI X coordinate
Y4a = S(:,2)';%[km] - ECI Y coordinate
Z4a = S(:,3)';%[km] - ECI Z coordinate

%Plot results
plot(X4a, Y4a, 'c', 'linewidth', 2)
fprintf("\nQuestion 4a:\nNote:\nThe trajectory shown in problem 4a is elliptical.\n")

%-------------------------------------------------------------------------------------------------%

%Part B done with true anomaly method

%b) p = 35000 e = 4

%For hyperbolic, solve for a, r, for v to find h
a4_b = 35000/(4^2 - 1);
rp4b = a4_b*(4-1);
Vp4b = sqrt(mu_e/a4_b + 2*mu_e/rp4b);
r = [rp4b;0;0];
v = [0;Vp4b;0];

%with known parameters, find and normalize vectors for basis
h_vec = cross(r,v);
e_vec = ((Vp4b^2-mu_e/rp4b)*r-(r'*v)*v)/mu_e;

z_hat2 = h_vec/norm(h_vec);
x_hat2 = e_vec/norm(e_vec);
y_hat2 = cross(z_hat2, x_hat2);

%for hyperbolic case, only consider theta from -pi/2 to pi/2
theta2 = (-pi/2 : pi/180 : pi/2);
R2 = (norm(h_vec)^2/mu_e) * (1./(1+norm(e_vec)*cos(theta2))); % satellite position based on true anomaly
 
r2 = zeros(3, length(theta2)); % create a matrix to fill with values of position
for i = 1:length(theta2)
    r2(:,i) = R2(i)*cos(theta2(i))*x_hat2 + R2(i)*sin(theta2(i))*y_hat2;
end

%Plot Results
plot(r2(1,:),r2(2,:), 'm', 'linewidth', 2)
fprintf("\nQuestion 4b:\nNote:\nThe trajectory shown in problem 4b is hyperbolic.\n")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Question 5%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use initial conditions from 4a

V_esc = sqrt(2*mu_e/rp4a); %[km/s], escape velocity from elliptical periapsis
r = [rp4a;0;0];        %[km], position vector
v = [0;V_esc;0];       %[km/s], velocity vector
s = [r; v];            %state vector, [km;km/s]

[~,S] = ode45(@Ortiz_HW1_2021_P2_EOM, tspan, s, OPTIONS);%integrate Equations of Motion for Problem 2
%S contains position and velocity for every second

%Split S matrix into its components
X5 = S(:,1)';%[km] - ECI X coordinate
Y5 = S(:,2)';%[km] - ECI Y coordinate
Z5 = S(:,3)';%[km] - ECI Z coordinate

%Plot results
plot(X5, Y5, 'g', 'linewidth', 2)
fprintf("\nQuestion 5:\nNote:\nThe trajectory shown in problem 5 is parabolic. The escape velocity required to reach this orbit is %g [km/sec].\n",V_esc)
legend('Earth','Problem 4a Trajectory','Problem 4b Trajectory','Problem 5 Trajectory')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Extra Credit%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
esp_5 = (V_esc^2/2)-(mu_e/rp4a); %[km^2/s^2], The specific energy of a parabolic orbit is zero
esp_4a = -(mu_e/(2*rp4a));%[km^2/s^2], specific energy for elliptical orbit
delta_esp = esp_5 - esp_4a; %[km^2/s^2]
fprintf("\nExtra Credit:\nThe difference in orbital specific energies between case 4a and case 5 is %g [km^2/s^2].\n",delta_esp)