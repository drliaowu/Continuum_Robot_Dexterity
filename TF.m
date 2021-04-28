function g = TF (phi, kappa, s)
%TF Homogeneous transformation for one segment
%
%   Eq.(1) in the paper
%   
%   g = TF (phi, kappa, s)
%   phi:    angle phi as defined in Fig.3
%   kappa:  curvature kappa as defined in Fig.3
%   s:      curve length s as defined in Fig.3
%   g:      homogeneous transformation, 4 x 4

g = [ cos(phi)^2*(cos(kappa*s)-1)+1,      sin(phi)*cos(phi)*(cos(kappa*s)-1),         cos(phi)*sin(kappa*s),  cos(phi)*(1-cos(kappa*s))/kappa;
      sin(phi)*cos(phi)*(cos(kappa*s)-1), cos(phi)^2*(1-cos(kappa*s))+cos(kappa*s),   sin(phi)*sin(kappa*s),  sin(phi)*(1-cos(kappa*s))/kappa;
      -cos(phi)*sin(kappa*s),             -sin(phi)*sin(kappa*s),                     cos(kappa*s),           sin(kappa*s)/kappa;
      0,                                  0,                                          0,                      1];