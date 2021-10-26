function dy = relacc_with_j2(t, u)

%constants of the problem START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 6.6742*10^(-20); %m^3/(kg*s^2) gravitational constant
m1 = 5.972*10^24; %kg mass of Earth
mu = G*m1; %km^3/s^2 gravitational parameter
J2 = 0.00108263; % [-] second zonal harmonic
Re = 6378; %km Earth's radius
%constants of the problem END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=u(1);
y=u(2);
z=u(3);

r = sqrt(x^2+y^2+z^2);

dy = [u(4); u(5); u(6); (-mu/r^3)*[x; y; z] + ... % Keplerian dynamics
    (3/2)*(J2*mu*Re^2/r^4)*[(x/r)*(5*(z^2/r^2)-1); (y/r)*(5*(z^2/r^2)-1); (z/r)*(5*(z^2/r^2)-3)]]; % Earth's oblateness (J2)

end

