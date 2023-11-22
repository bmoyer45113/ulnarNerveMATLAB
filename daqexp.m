clc; clear variables; close all;
d = daq('ni');

% Initialize standard transformational matrices
syms theta

Xrot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Yrot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Zrot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

% Initialize system with 6 different theta values and set constant values
% and pre set size of positional vectors
x0 = 0; y0 = 0; z0 = 0;

counter = 1;

%number of data points
time = 400;

%{
general equation for data: seconds of data = 0.135 * data points
or 7.41 Hz
%}

%pre allocation of angles and end affector points
datapts = time;
theta1 = zeros(datapts,1);
theta2 = zeros(datapts,1);
theta3 = zeros(datapts,1);
theta4 = zeros(datapts,1);
theta5 = zeros(datapts,1);
theta6 = zeros(datapts,1);
x = zeros(datapts,1);
y = zeros(datapts,1);
z = zeros(datapts,1);
botright = zeros(datapts,1);
topright = zeros(datapts,1);
topleft = zeros(datapts,1);
thetax = zeros(datapts,1);
thetay = zeros(datapts,1);
thetaz = zeros(datapts,1);
xo = zeros(datapts,1);
yo = zeros(datapts,1);
zo = zeros(datapts,1);

% DAQ variable setup
adtheta1 = addinput(d,"Dev1","ai7","Voltage");
adtheta2 = addinput(d,"Dev1","ai6","Voltage");
adtheta3 = addinput(d,"Dev1","ai5","Voltage");
adtheta4 = addinput(d,"Dev1","ai4","Voltage");
adtheta5 = addinput(d,"Dev1","ai3","Voltage");
adtheta6 = addinput(d,"Dev1","ai2","Voltage");
adtheta1.TerminalConfig ="SingleEnded";
adtheta2.TerminalConfig ="SingleEnded";
adtheta3.TerminalConfig ="SingleEnded";
adtheta4.TerminalConfig ="SingleEnded";
adtheta5.TerminalConfig ="SingleEnded";
adtheta6.TerminalConfig ="SingleEnded";
adtheta1.Name = "Joint 1";
adtheta2.Name = "Joint 2";
adtheta3.Name = "Joint 3";
adtheta4.Name = "Joint 4";
adtheta5.Name = "Joint 5";
adtheta6.Name = "Joint 6";

% read into a numerical array
T = read(d);
B = table2array(T);

% initial positions of all joints
T1 = B(1,1) - 1.1534; 
T2 = B(1,2); 
T3 = B(1,3); 
T4 = B(1,4);
T5 = B(1,5);
T6 = B(1,6);


%Link lengths
z1 = 5.4; %5.4 const arm  base to first clevis length
x2 = 33.05; %33.05 const arm shaft to shaft of second joint arm length
z3 = 17.5; %17.5 const arm elbow to (end of orange/beginning of black) motor housing length
z4 = 13.1; %13.1 const arm black rotation to end chamber shaft length
z5 = 4.45; %4.45 const arm end chamber shaft to bottom of sensor clamp length
z6 = 14.2; %14.2 chisel pointer length
y7 = 2.2; %2.2 theoretical sensor width

% main loop for one set of data points at one time
while time>0
    %read voltage data and convert to degrees
    D = read(d);
    C = table2array(D);
    theta1(counter) = (-C(1,1) / 3.333 - T1) * 1.16957; %1.16957
    theta2(counter) = (C(1,2) - T2) * 1.33498; %1.33498
    theta3(counter) = (C(1,3) - T3) * 1.06297; %1.06297
    theta4(counter) = (C(1,4) - T4) * 1.62039; %1.62039
    theta5(counter) = (C(1,5) - T5) * 1.2535; %1.2535
    theta6(counter) = (C(1,6) - T6) * 1.09; %1.09
 
    
    % while variable decrement
    time = time - 1;

    %Substitute matrices with actual theta values and input positional
    %vectors creating each joints transfomational matrix
    H1 = [subs(Zrot,theta,theta1(counter)) [0;0;0]; 0 0 0 1];
    H2 = [subs(Yrot,theta,-theta2(counter)) [0;0;z1]; 0 0 0 1];
    H3 = [subs(Yrot,theta,theta3(counter)) [x2;0;0]; 0 0 0 1];
    H4 = [subs(Zrot,theta,theta4(counter)) [0;0;z3]; 0 0 0 1];
    H5 = [subs(Yrot,theta,theta5(counter)) [0;0;z4]; 0 0 0 1];
    H6 = [subs(Zrot,theta,theta6(counter)) [0;0;(z5+z6)]; 0 0 0 1];
    
    % orientation point
    H7 = [subs(Zrot,theta,theta6(counter)) [0;y7;0]; 0 0 0 1];
    
    %Getting joint 1 positional data
    M1 = H1*H2;
    X1 = M1(1,4);
    Y1 = M1(2,4);
    Z1 = M1(3,4);
    
    %Getting joint 2 positional data
    M2 = M1*H3;
    X2 = M2(1,4);
    Y2 = M2(2,4);
    Z2 = M2(3,4);
    
    %Getting joint 3 positional data
    M3 = M2*H4;
    X3 = M3(1,4);
    Y3 = M3(2,4);
    Z3 = M3(3,4);
    
    %Getting joint 4 positional data
    M4 = M3*H5;
    X4 = M4(1,4);
    Y4 = M4(2,4);
    Z4 = M4(3,4);
    
    %Getting joint 5 positional data
    M5 = M4*H6;
    X5 = M5(1,4);
    Y5 = M5(2,4);
    Z5 = M5(3,4);
    
    %Getting orientation data
    orien = M5*H7;
    X6 = orien(1,4);
    Y6 = orien(2,4);
    Z6 = orien(3,4);
    xo(counter) = orien(1,4);
    yo(counter) = orien(2,4);
    zo(counter) = orien(3,4);
    
    %getting the end effector positional data
	FK = H1*H2*H3*H4*H5*H6;
    botright(counter) = orien(3,3);
    topright(counter) = orien(1,3);
    topleft(counter) = orien(1,1);
    x(counter) = FK(1,4);
    y(counter) = FK(2,4);
    z(counter) = FK(3,4);
   
    %plot data
    p = plot3([x0 X1],[y0 Y1],[z0 Z1],[X1 X2],[Y1 Y2],[Z1 Z2],[X2 X3],[Y2 Y3],[Z2 Z3],[X3 X4],[Y3 Y4],[Z3 Z4],[X4 X5],[Y4 Y5],[Z4 Z5],[X5 X6],[Y5 Y6],[Z5 Z6],x,y,z,'linewidth',4);
    p(1).Color = [128 128 128]/255;
    p(2).Color = [0.7 0.7 0.7];
    p(3).Color = [0.91 0.41 0.17];
    p(4).Color = [0.2 0.2 0.2];
    p(5).Color = [0.91 0.41 0.17];
    p(6).Color = [0.576 0.878 0.961];
    p(7).LineWidth = 1;
    p(7).Color = 'black';
    set(gcf,'position',[450 100 1000 1000]);
    
    xlim([-70 70])
    ylim([-70 70])
    zlim([-10 70])
    xlabel('X-Axis');
    ylabel('Y-Axis');
    zlabel('Z-Axis');
    grid on
    drawnow()
    pause(0.005);

    % internal counting variable increment
    counter = counter+1;
    
end

% end affector 3 angles position
for i = 1:datapts
    thetay(i) = asin(topright(i));
    thetaz(i) = acos(topleft(i)/cos(thetay(i)));
    thetax(i) = acos(botright(i)/cos(thetay(i)));
end

% data organization and write to excel
mat = horzcat(x,y,z); %xo,yo,zo should be added for orientation
delete('endaffector.xlsx');
writematrix(mat,'endaffector.xlsx');
