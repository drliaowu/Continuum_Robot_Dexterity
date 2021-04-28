function cfg = Config (s1, phi2, kappa2, s2, phi3, kappa3, s3, Ntheta, Nh, deltaX, deltaZ)
%Config  Calculate the configuration of the tip of the robot
%
%   cfg = Config (s1, phi2, kappa2, s2, phi3, kappa3, s3, Ntheta, Nh, deltaX, deltaZ)
%   s1:     curve length of the first segment
%   phi2:   angle of the second segment
%   kappa2: curvature of the second segment
%   s2:     curve length of the second segment
%   phi3:   angle of the third segment
%   kappa3: curvature of the third segment
%   s3:     curve length of the third segment
%   Ntheta: number of intervals for theta
%   Nh:     number of intervals for h
%   deltaX: interval of position X
%   deltaZ: interval of position Z
%   cfg:    [xx, z, theta, psi] - the index of the tip point in terms of
%   position and orientation, after discretization of the workspace

g = ForwardKinematics(s1, phi2, kappa2, s2, phi3, kappa3, s3);

x = g(1, 4); y = g(2, 4); z = g(3, 4);
xx = sqrt(x^2 + y^2);

alpha = atan2(y, x);

oz1 = g(1, 3); oz2 = g(2, 3); oz3 = g(3, 3);

theta = atan2(oz2, oz1) - alpha;
theta = atan2(sin(theta), cos(theta));

xx = floor(xx / deltaX);

z = floor(z / deltaZ);

theta = floor(theta / (2 * pi / Ntheta)) + 0.5;

psi = floor(oz3 / (2 / Nh)) + 0.5;

if psi > Nh / 2
    psi = Nh / 2 - 0.5;
end

cfg = [xx, z, theta, psi];