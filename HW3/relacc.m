function dy = relacc(t, u)

%constants of the problem START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 6.6742*10^(-20); %km^3/(kg*s^2) gravitational constant
m1 = 5.972*10^24; %kg mass of Earth
mu = G*m1; %km^3/s^2 gravitational parameter
%constants of the problem END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=u(1);
y=u(2);
z=u(3);

dy = [u(4); u(5); u(6); (-mu/((sqrt(x^2+y^2+z^2))^3))*[x; y; z]];

end
