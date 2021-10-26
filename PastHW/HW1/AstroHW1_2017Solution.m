% Carlos Ortiz Astrodynamics, HW 1 2017 Practice
clc;clear all;close all;

%%Define constants and initialize figure
mu_e = 398600; %[km^3/s^2], Gravitational Parameter assuming me>>>>>m
RE = 6378; %[km], earth's radius

[EX,EY,EZ] = sphere;
figure(1), hold on, axis equal, grid on
xlabel('ECI X - Position [km]')
ylabel('ECI Y - Position [km]')
zlabel('ECI Z - Position [km]')
earth =surf(EX*RE, EY*RE, EZ*RE);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Question 1%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%   Given Perigee, rp = 10000km, and apogee, ra = 30000 find:
rp = 10000; %[km]
ra = 30000; %[km]

%A)Eccentricity, e
e = (ra - rp)/(ra + rp); %[]

%B)Semi-Major Axis, a 
a = ra/(1+e); %[km] also works with a = 1/2(ra+rp)

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
fprintf("Question 1:\nA) e = %g\nB) a = %d [km]\nC) T = %e [sec]\nD) Vp = %f [km/sec]\nE) Va = %f [km/s]\nF) hp = %e [km^2/s]\nG) ha = %e [km^2/s]\n",e,a,T,Vp,Va,hp,ha)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Question 2%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use ODE45 to integrate EOM for one orbital period
%using a min time step of 10 seconds

%Initial conditions:
r = [rp;0;0];   %[km], given
v = [0; Vp; 0]; %[km/s], solved for previously
s = [r;v];      %State vector: [km; km/sec]

tspan = 0:T;                    %[sec], Time span in one orbital period
OPTIONS = odeset('Maxstep',10); %set min time span of 10 seconds between takes
[~,S] = ode45(@HW1_2017_P2_EOM, tspan, s, OPTIONS);%integrate Equations of Motion for Problem 2
%S contains position and velocity for every second

%Split S matrix into its components
Xa = S(:,1)';%[km] - ECI X coordinate
Ya = S(:,2)';%[km] - ECI Y coordinate
Za = S(:,3)';%[km] - ECI Z coordinate

%Plot results
plot3(Xa, Ya, Za, 'r', 'linewidth', 2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Question 3%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeat, but plot using true anomoly (theta) approach

%First define angular momentum vector (establishes plane of orbit)
%get unit vector in direction of h and of e, find the vector in J with
%cross product >>>> Establish basis in true anomoly reference frame

h_vec = cross (r,v); %specific angular momentum vector, [km^2/sec]
e_vec = ((Vp^2-mu_e/rp)*r-(r'*v)*v)/mu_e; %eccentricity Vector

z_hat = h_vec/norm(h_vec); %Unit vector normal to orbital plane
x_hat = e_vec/norm(e_vec); %Unit vector along e
y_hat = cross(z_hat, x_hat); %Complete right handed basis

theta = (1*pi/180 : pi/180 : 360*pi/180); %set a range of true anomoly [rads]
R = (norm(h_vec)^2/mu_e) * (1./(1+norm(e_vec)*cos(theta)));%satellite position based on true anomaly
%MAYBE put e instead of e_vec if it doesn't work

r = zeros(3, length(theta));
for i = 1:length(theta)
    r(:,i) = R(i)*cos(theta(i))*x_hat + R(i)*sin(theta(i))*y_hat;
end

%Plot Results
plot3(r(1,:),r(2,:),r(3,:), 'k')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Question 4%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%redo problem 2 using same semimajor axis, but circular.
r = [a;0;0];           %[km], position vector
v = [0;sqrt(mu_e/a);0]; %[km/s], velocity vector
s = [r, v];            %state vector, [km;km/s]

[~,S] = ode45(@HW1_2017_P2_EOM, tspan, s, OPTIONS);
Xb = S(:,1)';%[km] - ECI X coordinate
Yb = S(:,2)';%[km] - ECI Y coordinate
Zb = S(:,3)';%[km] - ECI Z coordinate
plot3(Xb, Yb, Zb, 'b', 'linewidth', 2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Question 5%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%What is minimum velocity of for spacecraft in circular problem 4
%to escape earth orbit?
% > Parabolic orbit

v_esc = sqrt(2*mu_e/a);%[m/s], escape velocity from current location

r = [a;0;0]; %[km], position vector
v = [0;v_esc;0]; %[km/s], velocity vector to escape orbit
s = [r, v];%state vector, [km;km/s]

[~,S] = ode45(@HW1_2017_P2_EOM, tspan, s, OPTIONS);
Xc = S(:,1)';%[km] - ECI X coordinate
Yc = S(:,2)';%[km] - ECI Y coordinate
Zc = S(:,3)';%[km] - ECI Z coordinate
plot3(Xc, Yc, Zc, 'g', 'linewidth', 2)
fprintf("Question 5:\nTo escape into a PARABOLIC trajectory, the escape velocity is %f [km/s].\n",v_esc)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Extra Credit%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the difference in the orbital specific energies in the craft in part
%4 and in part 5

esp4 = (-mu_e)/(2*a);%[km^2/s^2], orbital specific energy from problem 4
esp5 = (-mu_e)/a+v_esc^2/2; %[km^2/s^2], orbital specific energy from problem 5 (should be 0)
delta_esp = esp4-esp5; %[km^2/s^2], difference in energies
fprintf("Extra Credit:\nThe difference in orbital specific energies is %f [km^2/s^2]",delta_esp)
%%
legend('Earth','Problem 2 Trajectory','Problem 3 Trajectory','Problem 4 Trajectory','Problem 5 Trajectory')