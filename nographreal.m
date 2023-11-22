clc; clear variables; close all;
d = daq('ni');

% Initialize standard transformational matrices
syms theta

Xrot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Yrot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Zrot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

                    % total sample time (variables)
                            sampleTime = 45;
                             daqRate = 40;

%number of data points
sampleArray = sampleTime*daqRate;
d.Rate = daqRate;

%pre allocation of angles and end affector points
theta1 = zeros(sampleArray,1);
theta2 = zeros(sampleArray,1);
theta3 = zeros(sampleArray,1);
theta4 = zeros(sampleArray,1);
theta5 = zeros(sampleArray,1);
theta6 = zeros(sampleArray,1);
x = zeros(sampleArray,1);
y = zeros(sampleArray,1);
z = zeros(sampleArray,1);
xo = zeros(sampleArray,1);
yo = zeros(sampleArray,1);
zo = zeros(sampleArray,1);

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

%Link lengths
z1 = 5.4; %5.4 const arm  base to first clevis length
x2 = 33.05; %33.05 const arm shaft to shaft of second joint arm length
z3 = 17.5; %17.5 const arm elbow to (end of orange/beginning of black) motor housing length
z4 = 13.1; %13.1 const arm black rotation to end chamber shaft length
z5 = 4.45; %4.45 const arm end chamber shaft to bottom of sensor clamp length
z6 = 14.2; %14.2 chisel pointer length
y7 = 2.2; %2.2 theoretical sensor width

P = imread("green.jpg");
imshow(P)

T = read(d);
B = table2array(T);

D = read(d, seconds(sampleTime));
C = table2array(D);

K = imread("Red.jpg");
imshow(K)

for i = 1:sampleArray
    %read voltage data and convert to degrees
    
    theta1(i) = (-C(i,1) / 3.333 - (B(1,1)-1.1534)) * 1.16957; %1.16957
    theta2(i) = (C(i,2) - B(1,2)) * 1.33498; %1.33498
    theta3(i) = (C(i,3) - B(1,3)) * 1.06297; %1.06297
    theta4(i) = (C(i,4) - B(1,4)) * 1.62039; %1.62039
    theta5(i) = (C(i,5) - B(1,5)) * 1.2535; %1.2535
    theta6(i) = (C(i,6) - B(1,6)) * 1.09; %1.09
    
    %Substitute matrices with actual theta values and input positional
    %vectors creating each joints transfomational matrix
    H1 = [subs(Zrot,theta,theta1(i)) [0;0;0]; 0 0 0 1];
    H2 = [subs(Yrot,theta,-theta2(i)) [0;0;z1]; 0 0 0 1];
    H3 = [subs(Yrot,theta,theta3(i)) [x2;0;0]; 0 0 0 1];
    H4 = [subs(Zrot,theta,theta4(i)) [0;0;z3]; 0 0 0 1];
    H5 = [subs(Yrot,theta,theta5(i)) [0;0;z4]; 0 0 0 1];
    H6 = [subs(Zrot,theta,theta6(i)) [0;0;(z5+z6)]; 0 0 0 1];
    
    % orientation point
    H7 = [subs(Zrot,theta,theta6(i)) [0;y7;0]; 0 0 0 1];
    
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
    xo(i) = orien(1,4);
    yo(i) = orien(2,4);
    zo(i) = orien(3,4);
    
    %getting the end effector positional data
	FK = H1*H2*H3*H4*H5*H6;
    x(i) = FK(1,4);
    y(i) = FK(2,4);
    z(i) = FK(3,4);
    
end

mat = horzcat(x,y,z,xo,yo,zo);
delete('endaffector2.xlsx');
writematrix(mat,'endaffector2.xlsx');