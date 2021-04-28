% Code for calculating the dexterity of three continuum robots combining
% concentric tube mechanisms and cable driven mechanisms.
%
% Please Refer to:
%
% Liao Wu, Ross Crawford, Jonathan Roberts. Dexterity analysis of three 
% 6-DOF continuum robots combining concentric tube mechanisms and cable 
% driven mechanisms. IEEE Robotics and Automation Letters. 2017, 2(2): 
% 514-521.

clear
clc
close all

%% set up parameters for simulation
N       = 10000;        % number of samples, for quick run
% N       = 100000000;    % number of samples, used in the paper
Ntheta  = 60;           % number of intervals for theta
Nh      = 30;           % number of intervals for h
deltaX  = 5;            % interval of position X
deltaZ  = 5;            % interval of position Z

% select case number
casenumber = 1; % Robot 1 in Group I
% casenumber = 2; % Robot 2 in Group I
% casenumber = 3; % Robot 3 in Group I
% casenumber = 4; % Robot 1 in Group II
% casenumber = 5; % Robot 2 in Group II
% casenumber = 6; % Robot 3 in Group II

% the following parameters are shared by all the cases
s1      = random('unif',0,1,N,1) * 100;

phi2    = random('unif',0,1,N,1) * 2 * pi;

phi3    = random('unif',0,1,N,1) * 2 * pi;

% the following parameters are distinguished case by case
switch casenumber
    case 1  % Robot 1 in Group I
        kappa2  =   random('unif',0,1,N,1)  *   pi / 200; 
        s2      =   random('unif',0,1,N,1)  *   100;         
        kappa3  =   ones(N,1)               *   pi / 200;                
        s3      =   random('unif',0,1,N,1)  *   100;          
    case 2  % Robot 2 in Group I
        kappa2  =   ones(N,1)               *   pi / 200;                
        s2      =   random('unif',0,1,N,1)  *   100;         
        kappa3  =   random('unif',0,1,N,1)  *   pi / 200;  
        s3      =   random('unif',0,1,N,1)  *   100;          
    case 3  % Robot 3 in Group I
        kappa2  =   random('unif',0,1,N,1)  *   pi / 200; 
        s2      =   [zeros(N/2,1);  random('unif',0,1,N/2,1) * 100]; 
        kappa3  =   random('unif',0,1,N,1)  *   pi / 200;  
        s3      =   [random('unif',0,1,N/2,1) * 100;  ones(N/2,1) * 100];  
    case 4  % Robot 1 in Group II
        kappa2  =   random('unif',0,1,N,1)  *   pi / 100; 
        s2      =   random('unif',0,1,N,1)  *   100;         
        kappa3  =   ones(N,1)               *   pi / 200;               
        s3      =   random('unif',0,1,N,1)  *   100;          
    case 5  % Robot 2 in Group II
        kappa2  =   ones(N,1)               *   pi / 200;                
        s2      =   random('unif',0,1,N,1)  *   100;         
        kappa3  =   random('unif',0,1,N,1)  *   pi / 100;  
        s3      =   random('unif',0,1,N,1)  *   100;          
    case 6  % Robot 3 in Group II
        kappa2  =   random('unif',0,1,N,1)  *   pi / 100; 
        s2      =   [zeros(N/2,1);  random('unif',0,1,N/2,1) * 100]; 
        kappa3  =   random('unif',0,1,N,1)  *   pi / 100;  
        s3      =   [random('unif',0,1,N/2,1) * 100;  ones(N/2,1) * 100];  
end

%% calculate the statistics of the configurations of all the samples
C=zeros(N,4);

for i=1:N
    C(i,:)=Config(s1(i),phi2(i),kappa2(i),s2(i),phi3(i),kappa3(i),s3(i),Ntheta,Nh,deltaX,deltaZ);
end

D=unique(C,'rows');

[E,F]=unique(D(:,1:2),'rows');

G=zeros(size(F,1),1);
for i=1:size(F,1)-1
    G(i)=F(i+1)-F(i);
end
G(size(F,1))=size(D,1)-F(size(F,1))+1;

H=zeros(size(G,1),4);
H(:,1)=G/Ntheta/Nh;
for i=1:size(G,1)-1
    H(i,2)=sum(abs(cos(D(F(i):F(i+1)-1,3)*2*pi/Ntheta).*sqrt(1-(D(F(i):F(i+1)-1,4)*2/Nh).^2)))/Ntheta/Nh;
    H(i,3)=sum(abs(sin(D(F(i):F(i+1)-1,3)*2*pi/Ntheta).*sqrt(1-(D(F(i):F(i+1)-1,4)*2/Nh).^2)))/Ntheta/Nh;
    H(i,4)=sum(abs(D(F(i):F(i+1)-1,4)*2/Nh))/Ntheta/Nh;
end
H(size(G,1),2)=sum(abs(cos(D(F(size(G,1)):end,3)*2*pi/Ntheta).*sqrt(1-(D(F(size(G,1)):end,4)*2/Nh).^2)))/Ntheta/Nh;
H(size(G,1),3)=sum(abs(sin(D(F(size(G,1)):end,3)*2*pi/Ntheta).*sqrt(1-(D(F(size(G,1)):end,4)*2/Nh).^2)))/Ntheta/Nh;
H(size(G,1),4)=sum(abs(D(F(size(G,1)):end,4)*2/Nh))/Ntheta/Nh;

%% draw workspace and dexterity
%total dexterity
figure();

hold on;
if casenumber<=3
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,25,H(:,1),'s','filled');
    axis equal;
    axis([-40,180,-20,310]);
else
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,15,H(:,1),'s','filled');
    axis equal;
    axis([-40,180,-110,310]);
end

colorbar;

[a,b]=max(H(:,1));

if casenumber>3
    scatter((E(b,1))*deltaX,(E(b,2))*deltaZ,25,'m','s','filled');
else
    scatter((E(b,1))*deltaX,(E(b,2))*deltaZ,15,'m','s','filled');
end

xlabel('X /mm');
ylabel('Z /mm');
hold off;

%radial dexterity
figure();

hold on;
if casenumber<=3
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,25,H(:,2),'s','filled');
    axis equal;
    axis([-40,180,-20,310]);
else
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,15,H(:,2),'s','filled');
    axis equal;
    axis([-40,180,-110,310]);
end

colorbar;

[ax,bx]=max(H(:,2));

if casenumber>3
    scatter((E(bx,1))*deltaX,(E(bx,2))*deltaZ,25,'m','s','filled');
else
    scatter((E(bx,1))*deltaX,(E(bx,2))*deltaZ,15,'m','s','filled');
end

xlabel('X /mm');
ylabel('Z /mm');
hold off;

%circumferential dexterity
figure();

hold on;
if casenumber<=3
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,25,H(:,3),'s','filled');
    axis equal;
    axis([-40,180,-20,310]);
else
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,15,H(:,3),'s','filled');
    axis equal;
    axis([-40,180,-110,310]);
end

colorbar;

[ay,by]=max(H(:,3));

if casenumber>3
    scatter((E(by,1))*deltaX,(E(by,2))*deltaZ,25,'m','s','filled');
else
    scatter((E(by,1))*deltaX,(E(by,2))*deltaZ,15,'m','s','filled');
end

xlabel('X /mm');
ylabel('Z /mm');
hold off;

%axial dexterity
figure();

hold on;
if casenumber<=3
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,25,H(:,4),'s','filled');
    axis equal;
    axis([-40,180,-20,310]);
else
    scatter(E(:,1)*deltaX,E(:,2)*deltaZ,15,H(:,4),'s','filled');
    axis equal;
    axis([-40,180,-110,310]);
end

colorbar;

[az,bz]=max(H(:,4));

if casenumber>3
    scatter((E(bz,1))*deltaX,(E(bz,2))*deltaZ,25,'m','s','filled');
else
    scatter((E(bz,1))*deltaX,(E(bz,2))*deltaZ,15,'m','s','filled');
end

xlabel('X /mm');
ylabel('Z /mm');
hold off;

%% draw orientation
figure();

hold on;
theta = linspace(0,2*pi,Ntheta+1);
phi = asin(linspace(-1,1,Nh+1));
[theta,phi] = meshgrid(theta,phi);
[xs,ys,zs] = sph2cart(theta,phi,1);
surf(xs,ys,zs);
colormap([1,1,1]);
for j=F(b):F(b+1)-1
    patchx1=cos((D(j,3)-0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)-0.5)*2/Nh)^2);
    patchx2=cos((D(j,3)+0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)-0.5)*2/Nh)^2);
    patchx3=cos((D(j,3)-0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)+0.5)*2/Nh)^2);
    patchx4=cos((D(j,3)+0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)+0.5)*2/Nh)^2);
    patchy1=sin((D(j,3)-0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)-0.5)*2/Nh)^2);
    patchy2=sin((D(j,3)+0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)-0.5)*2/Nh)^2);
    patchy3=sin((D(j,3)-0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)+0.5)*2/Nh)^2);
    patchy4=sin((D(j,3)+0.5)*2*pi/Ntheta)*sqrt(1-((D(j,4)+0.5)*2/Nh)^2);
    patchz1=(D(j,4)-0.5)*2/Nh;
    patchz2=(D(j,4)-0.5)*2/Nh;
    patchz3=(D(j,4)+0.5)*2/Nh;
    patchz4=(D(j,4)+0.5)*2/Nh;
    
    color=zeros(2,2,3);
    color(1,1,:)=[0;1;0];
    color(1,2,:)=color(1,1,:);
    color(2,1,:)=color(1,1,:);
    color(2,2,:)=color(1,1,:);
    surf([patchx1,patchx2;patchx3,patchx4;],[patchy1,patchy2;patchy3,patchy4;],[patchz1,patchz2;patchz3,patchz4;],color);
end
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off;

%% calculate the dexterity indices

% display simulation case
switch casenumber
    case 1
        disp ("Group I - Robot 1:")
    case 2
        disp ("Group I - Robot 2:")
    case 3
        disp ("Group I - Robot 3:")
    case 4
        disp ("Group II - Robot 1:")
    case 5
        disp ("Group II - Robot 2:")
    case 6
        disp ("Group II - Robot 3:")
end

% dexterity
D_t = sum(H(:,1))/size(H,1);
D_r = sum(H(:,2))/size(H,1);
D_c = sum(H(:,3))/size(H,1);
D_a = sum(H(:,4))/size(H,1);

disp ("total dexterity: " + D_t)                % total dexterity
disp ("radial dexterity: " + D_r)               % radial dexterity
disp ("circumferential dexterity: " + D_c)      % circumferential dexterity
disp ("axial dexterity: " + D_a)                % axial dexterity

% workspace
W=zeros(1,floor(a/0.05)+1);

for k=0:floor(a/0.05)
    W(k+1)=size(find(H(:,1)>=k*0.05),1)*deltaX*deltaZ;
end
disp ("workspace: " + W(1))                     % workspace

% dexterous workspace
W_D = round(D_t * W(1), 0);

disp ("dexterous workspace: " + W_D)            % dexterous workspace

% max dexterity
disp ("max total dexterity: " + a)              % max total dexterity
disp ("max radial dexterity: " + ax)            % max radial dexterity
disp ("max circumferential dexterity: " + ay)   % max circumferential dexterity
disp ("max axial dexterity: " + az)             % max axial dexterity