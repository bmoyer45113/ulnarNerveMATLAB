clc; clear variables; close all;
d = daq('ni');

                    % total sample time (variables)
                            sampleTime = 10;
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

%daq channel setup
adtheta1 = addinput(d,"Dev1","ai7","Voltage");
adtheta2 = addinput(d,"Dev1","ai6","Voltage");
adtheta3 = addinput(d,"Dev1","ai5","Voltage");
adtheta4 = addinput(d,"Dev1","ai4","Voltage");
adtheta5 = addinput(d,"Dev1","ai3","Voltage");
adtheta6 = addinput(d,"Dev1","ai2","Voltage");
adtheta1.TerminalConfig = "SingleEnded";
adtheta2.TerminalConfig = "SingleEnded";
adtheta3.TerminalConfig = "SingleEnded";
adtheta4.TerminalConfig = "SingleEnded";
adtheta5.TerminalConfig = "SingleEnded";
adtheta6.TerminalConfig = "SingleEnded";
adtheta1.Name = "Joint 1";
adtheta2.Name = "Joint 2";
adtheta3.Name = "Joint 3";
adtheta4.Name = "Joint 4";
adtheta5.Name = "Joint 5";
adtheta6.Name = "Joint 6";

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
end

mat = horzcat(theta1,theta2,theta3,theta4,theta5,theta6);
delete('armData1.xlsx');
writematrix(mat,'armData1.xlsx');