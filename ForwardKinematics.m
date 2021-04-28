function g = ForwardKinematics (s1, phi2, kappa2, s2, phi3, kappa3, s3)
%ForwardKinematics  Forward kinemtaics of three segments
%   
%   The frist segment is always straight in the paper
%
%   g = ForwardKinematics (s1, phi2, kappa2, s2, phi3, kappa3, s3)
%   s1:     curve length of the first segment
%   phi2:   angle of the second segment
%   kappa2: curvature of the second segment
%   s2:     curve length of the second segment
%   phi3:   angle of the third segment
%   kappa3: curvature of the third segment
%   s3:     curve length of the third segment

g = TF(phi2,kappa2,s2) * TF(phi3,kappa3,s3);

g(3,4) = g(3,4) + s1;