clc; clear variables; close all;

%frontend command window UI
prompt = "Live Plotting at 7.41 Hz with no path coords, no Live Plotting at 40 Hz with path coords [1/2]: ";
plotState = input(prompt);

if plotState == 1

    % just does live plotting and ends, not for data processing just for
    % explanation and demonstration, as well as a good wellness check
    daqExp();
end
if plotState == 2
    
    %this defines the name for the intermediate data and allows for a new
    %file name every iteration loop
    increment = 0; 

    %dataAquisition is the loop variable that allows the user to take data
    %without running path coordinates
    dataAquisition = 0;

    %main loop for the data aquisition when you want to take data but dont
    %want to run path coordinates
    while dataAquisition == 0

        %dataAq UI
        prompt2 = "how many seconds?: ";
        samTime = input(prompt2);
        [theta1, theta2, theta3, theta4, theta5, theta6] = dataAq(samTime);
    
        %forwardKinematics UI
        prompt3 = "Data Aquisition completed. Please type the name of the ultrasound data file:";
        ultrasoundDataFile = input(prompt3,"s");
        [x,y,z,xu,yu,zu,xo,yo,zo] = forwardKinematics(theta1,theta2,theta3,theta4,theta5,theta6,ultrasoundDataFile);
    
        %rangefinder UI
        
        %finds the size of one column of the data so we can define the end
        %point for rangefinder
        sizex = size(x);
        rfsize = sizex(1);

        %all of the default values of the inputs of rangefinder
        startpt = 1;
        endpt = rfsize;
        nFilter = 1;
        aFilter = 1;

        %function call for rangefinder the first time
        [nerveXYZ] = rangeFinderButter(x,y,z,xu,yu,zu,xo,yo,zo,startpt,endpt,aFilter,nFilter);
        promptA = "current run completed. would you like to save this data? [y/n]";
          answer = input(promptA,'s');
                if answer == 'y'
                    promptB = "would you like to continue to path coordinates? [y/n]";
                    answer2 = input(promptB,'s');
                    if answer2 == 'y'
                        writematrix(nerveXYZ, "nerveXYZ" + increment + ".xlsx");
                        writematrix([theta1,theta2,theta3,theta4,theta5,theta6],'xyzDataTest' + increment + '.xlsx');
                        dataAquisition = 1;
                        increment = increment+1;
                        continue;
                    end
                    if answer2 == 'n'
                        fprintf("data saved\n")
                        writematrix(nerveXYZ, "nerveXYZ" + increment + ".xlsx");
                        writematrix([theta1,theta2,theta3,theta4,theta5,theta6],'xyzDataTest' + increment + '.xlsx');
                        increment = increment+1;
    
                    end
                end
                if answer == 'n'
                    continue;
                end
        %rangefinder iteration loop for variable adjustment
        rangeFinder = 0;
        while rangeFinder == 0
            
            %function call and confirmation query
            [nerveXYZ] = rangeFinderButter(x,y,z,xu,yu,zu,xo,yo,zo,startpt,endpt,aFilter,nFilter);
            prompt4 = "Forward Kinematics completed. Please look at the plot and when it is in the right range and filtered, enter 1, otherwise, enter 0: ";
            rangeFinder = input(prompt4);

            if rangeFinder == 1
                continue
            end

            %start point query
            prompt5 = "new start point?: ";
            interme = input(prompt5);
            if isempty(interme)
                startpt = 1;
            else
                startpt = interme;
            end

            %end point query
            prompt6 = "new end point?: ";
            interme2 = input(prompt6);
            if isempty(interme2)
                endpt = rfsize;
            else
                endpt = interme2;
            end

            %arm filter query
            prompt7 = "new arm data filter?: ";
            interme3 = input(prompt7);
            if isempty(interme3)
                aFilter = 1;
            else
                aFilter = interme3;
            end

            %nerve filter query
            prompt8 = "new ultrasound data filter?: ";
            interme4 = input(prompt8);
            if isempty(interme4)
                nFilter = 1;
            else
                nFilter = interme4;
            end
        end
        dataAquisition = 1;
    end

    %pathCoord UI

    %function call for path coordinates
    [eb,en,et,tau,rho] = pathCoords(nerveXYZ);

    %exit or save query
    prompt9 = "pathCoords completed. type 1 to write to excel spreadsheet, type anything else to exit program: ";
    interme5 = input(prompt9);
    if interme5 == 1
        mat = [eb,en,et,tau,rho];
        delete('UlnarNerveData.xlsx');
        writematrix(mat,'UlnarNerveData.xlsx');
    end
end


%all functions

%Data Aquisition function
%inputs: 
%samTime: amount of time you want to take data for at 40 Hz
%outputs: 
%theta1-6: voltage arrays for each joint

function [theta1,theta2,theta3,theta4,theta5,theta6] = dataAq(samTime)
d = daq('ni');

% total sample time wanted from user (samTime)
sampleTime = samTime;
daqRate = 40;

%number of data points and rate change in the structure
sampleArray = sampleTime*daqRate;
d.Rate = daqRate;

%pre allocation of angles and end affector points (theta1 being base joint)
theta1 = zeros(sampleArray,1);
theta2 = zeros(sampleArray,1);
theta3 = zeros(sampleArray,1);
theta4 = zeros(sampleArray,1);
theta5 = zeros(sampleArray,1);
theta6 = zeros(sampleArray,1);

%daq channel setup

%initial channel setup
adtheta1 = addinput(d,"Dev1","ai7","Voltage");
adtheta2 = addinput(d,"Dev1","ai6","Voltage");
adtheta3 = addinput(d,"Dev1","ai5","Voltage");
adtheta4 = addinput(d,"Dev1","ai4","Voltage");
adtheta5 = addinput(d,"Dev1","ai3","Voltage");
adtheta6 = addinput(d,"Dev1","ai2","Voltage");

%terminal configuration of the new inputs in the structure d
adtheta1.TerminalConfig = "SingleEnded";
adtheta2.TerminalConfig = "SingleEnded";
adtheta3.TerminalConfig = "SingleEnded";
adtheta4.TerminalConfig = "SingleEnded";
adtheta5.TerminalConfig = "SingleEnded";
adtheta6.TerminalConfig = "SingleEnded";

%labels for each new channel
adtheta1.Name = "Joint 1";
adtheta2.Name = "Joint 2";
adtheta3.Name = "Joint 3";
adtheta4.Name = "Joint 4";
adtheta5.Name = "Joint 5";
adtheta6.Name = "Joint 6";

%used to show to user when data acquisition has started
P = imread("green.jpg");
imshow(P)

%reads in the initial voltage positions of the arm (1 x 6)
T = read(d);
B = table2array(T);

%reads in all of the data (sampleArray x 6)
D = read(d, seconds(sampleTime));
C = table2array(D);

%used to show to user when data acquisition has ended
K = imread("Red.jpg");
imshow(K)

    %for loop for converting all of the new data to degrees (comments are
    %most recent extrapolated values)
    for i = 1:sampleArray

        %read voltage data and convert to degrees
        theta1(i) = (-C(i,1) / 3.333 - (B(1,1)-1.1534)) * 1.16957; %1.16957
        theta2(i) = (C(i,2) - B(1,2)) * 1.33498; %1.33498
        theta3(i) = (C(i,3) - B(1,3)) * 1.06297; %1.06297
        theta4(i) = (C(i,4) - B(1,4)) * 1.62039; %1.62039
        theta5(i) = (C(i,5) - B(1,5)) * 1.2535; %1.2535
        theta6(i) = (C(i,6) - B(1,6)) * 1.09; %1.09
    end
end

%robot arm Kinematics function
%inputs:
%theta1-6: angle arrays from dataAq
%USData: premade rows x 4 sheet with YZ position of the ultrasound sensor
%outputs:
%x,y,z: position of the center of the ultrasound sensor, sampleArray x 3
%xo,yo,zo: position of one end of the ultrasound sensor, sampleArray x 3
%xu,yu,zu: position of the Ulnar Nerve, rows x 3

function [x,y,z,xu,yu,zu,xo,yo,zo] = forwardKinematics(theta1,theta2,theta3,theta4,theta5,theta6,USData)

% Initialize and define standard transformational matrices
syms theta

%not needed for this particular robot configuration but leaving it in for
%completeness
%Xrot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Yrot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Zrot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

%create size variable for robot data
size1 = size(theta1);
sampleArray = size1(1);

%define intermediate robot variable sizes

%FK is the main DH matrix after all of the transformations have been made,
%and is a list of every DH matrix one after another
FK = zeros(sampleArray*4,4);

%xyz position of the center of the US sensor
x = zeros(sampleArray,1);
y = zeros(sampleArray,1);
z = zeros(sampleArray,1);

%xyz position of one edge of the US sensor
xo = zeros(sampleArray,1);
yo = zeros(sampleArray,1);
zo = zeros(sampleArray,1);

%read in the ultrasound data and make it into an array
dataImport = readtable(USData);
ultrasound = table2array(dataImport);

%define intermediate ultrasound variable sizes and initial calculations
size2 = size(dataImport);
rows = size2(1);
nerve = zeros(rows*4,4);
nsp = zeros(rows*4,4);
robotTime = zeros(sampleArray,1);
sub = zeros(sampleArray,1);
rightRow = zeros(rows,1);
sampleRate = 40;
xu = zeros(rows,1);
yu = zeros(rows,1);
zu = zeros(rows,1);

%for loop that creates robot time increments
for j = 1:sampleArray, robotTime(j,1) = j/sampleRate; end

%Link lengths
z1 = 5.4; %5.4 const arm  base to first clevis length
x2 = 33.05; %33.05 const arm shaft to shaft of second joint arm length
z3 = 17.5; %17.5 const arm elbow to (end of orange/beginning of black) motor housing length
z4 = 13.1; %13.1 const arm black rotation to end chamber shaft length
z5 = 4.45; %4.45 const arm end chamber shaft to bottom of sensor clamp length
z6 = 10.9; %10.9 sensor length
y7 = 2.2; %2.2 theoretical sensor width

%main robot kinematics loop for everything except ultrasound
for i = 1:sampleArray

    %Substitute matrices with actual theta values and input positional
    %vectors creating each joints transfomational matrix
    H1 = [subs(Zrot,theta,theta1(i,:)) [0;0;0]; 0 0 0 1];
    H2 = [subs(Yrot,theta,-theta2(i,:)) [0;0;z1]; 0 0 0 1];
    H3 = [subs(Yrot,theta,theta3(i,:)) [x2;0;0]; 0 0 0 1];
    H4 = [subs(Zrot,theta,theta4(i,:)) [0;0;z3]; 0 0 0 1];
    H5 = [subs(Yrot,theta,theta5(i,:)) [0;0;z4]; 0 0 0 1];
    H6 = [subs(Zrot,theta,theta6(i,:)) [0;0;(z5+z6)]; 0 0 0 1];
    
    %orientation point
    H7 = [subs(Zrot,theta,theta6(i)) [0;y7;0]; 0 0 0 1];
    
    %Getting joint 1 positional data
    M1 = H1*H2;
    
    %Getting joint 2 positional data
    M2 = M1*H3;
    
    %Getting joint 3 positional data
    M3 = M2*H4;
    
    %Getting joint 4 positional data
    M4 = M3*H5;
    
    %Getting joint 5 positional data
    M5 = M4*H6;
    
    %Getting orientation position
    orien = M5*H7;
    xo(i) = orien(1,4);
    yo(i) = orien(2,4);
    zo(i) = orien(3,4);
    
    %getting the end effector positional data
	FK(1+(4*(i-1)):4+(4*(i-1)),1:4) = H1*H2*H3*H4*H5*H6;
    x(i) = FK(1+(4*(i-1)),4);
    y(i) = FK(2+(4*(i-1)),4);
    z(i) = FK(3+(4*(i-1)),4);
    
end

%for loop for each image iteration
for i = 1:rows

    %for loop that goes through every part of the sample array and
    %calculates the difference between the ultrasound sample time and robot
    %time
    for j = 1:sampleArray
        sub(j) = ultrasound(i,4) - robotTime(j,1);
    end

    %takes absolute value of the current sub to eliminate negatives 
    subabs = abs(sub);

    %finds the smallest value of subabs
    minVal = min(subabs);

    %goes through the subabs and puts the index of the smallest number into
    %rightRow so we know where we are close to the right point in time
    for k = 1:sampleArray
        if subabs(k) == minVal
            rightRow(i,1) = k;
        end
    end
    
    %uses rightRow to place the point and then interpolates every value in
    %the DH matrix individually 
    for j = 2:sampleArray-1

        %special case: we are assuming the first point of the arm is at the
        %first picture of the ultrasound, so instead of interpolating, we
        %are just forcing the first position of the arm to be the position
        %at the first image
        if rightRow(i) == 1
            for n = 1:4
                for m = 1:4
                    nsp(n+(4*(i-1)),m) = FK(n+(4*(i-1)),m);
                end
            end
        end

        %first case: this if case will pass true if we are at the right row
        %and the current image is closer to the point at t-1 than the point
        %at t+1
        if rightRow(i) == j && subabs(j-1) - subabs(j) <= subabs(j+1) - subabs(j)
            %interpolate every cell using FK at t and t-1
            for n = 1:4
                for m = 1:4
                    nsp(n+(4*(i-1)),m) = interp1([robotTime(j-1,1) robotTime(j,1)], [FK(n+(4*(j-2)),m) FK(n+(4*(j-1)),m)], ultrasound(i,4));
                end
            end
        end

        %second case: this if case will pass true if we are at the right row
        %and the current image is closer to the point at t+1 than the point
        %at t-1
        if rightRow(i) == j && subabs(j-1) - subabs(j) > subabs(j+1) - subabs(j)
            %interpolate every cell of FK using FK at t and t+1
            for n = 1:4
                for m = 1:4
                    nsp(n+(4*(i-1)),m) = interp1([robotTime(j,1) robotTime(j+1,1)], [FK(n+(4*(j-1)),m) FK(n+(4*j),m)], ultrasound(i,4));
                end
            end
        end
    end
end

%now that we have a number of new interpolated arm points at the time of
%the images, we use the positional data from the ultrasound as a new static 
%joint offset and calculate the nerves position as the matrix
%multiplication of the interpolated arm position times the static offset
    for i = 1:rows

        %ulnar nerve offset
        H8 = [eye(3) [0;ultrasound(i,2);ultrasound(i,3)]; 0 0 0 1];

        %Getting nerve position
        nerve(1+(4*(i-1)):4+(4*(i-1)),1:4) = nsp(1+(4*(i-1)):4+(4*(i-1)),1:4)*H8;
        xu(i) = nerve(1+(4*(i-1)),4);
        yu(i) = nerve(2+(4*(i-1)),4);
        zu(i) = nerve(3+(4*(i-1)),4);
    end
end

%rangeFinder Function: used to cut off unusable filtered data and verify the
%filtering before path coordinates are found
%inputs:
%x,y,z: position of the center of the ultrasound sensor, sampleArray x 3
%xo,yo,zo: position of one end of the ultrasound sensor, sampleArray x 3
%xu,yu,zu: position of the Ulnar Nerve, rows x 3
%startpt: the user defined value of the arm data to start the path
%endpt: the user defined value of the arm data that ends the path
%butterFilter: the user defined value used in the butterworth filter,
%generally lower values increase smoothing but increase wasted initial
%values
%outputs:


function [nerveXYZ] = rangeFinderButter(x,y,z,xu,yu,zu,xo,yo,zo,startpt,endpt,butterFilter,NerveFilter)

%arm concatenated array
A = [x,y,z,xo,yo,zo];

%nerve concatenated array
B = [xu,yu,zu];

%data smoothing and misc variable declaration
ulnSize = size(B);
preFilterCnt = A(startpt:endpt,1:3);
preFilterUln = B(1:ulnSize(1),1:3);
preFilterSide = A(startpt:endpt,4:6);

%nerve filtering

%cutoff frequency for nerve
Fu = NerveFilter;

%sample frequency for nerve
Fsu = 154/19.6;

%filter creation
[yu, xu] = butter(4,Fu/(Fsu/2));

%filtering individual dimensions
inputSignalxU = preFilterUln(:,1);
outSignalxU = filter(yu, xu, inputSignalxU);
inputSignalyU = preFilterUln(:,2);
outSignalyU = filter(yu, xu, inputSignalyU);
inputSignalzU = preFilterUln(:,3);
outSignalzU = filter(yu, xu, inputSignalzU);

%robot arm filtering

%cutoff frequency for robot
F = butterFilter;

%sample frequency for robot
Fs = 40;

%filter creation
[y, x] = butter(4,F/(Fs/2));

%filtering each dimension
inputSignalxC = preFilterCnt(:,1);
outSignalxC = filter(y, x, inputSignalxC);
inputSignalyC = preFilterCnt(:,2);
outSignalyC = filter(y, x, inputSignalyC);
inputSignalzC = preFilterCnt(:,3);
outSignalzC = filter(y, x, inputSignalzC);
inputSignalxS = preFilterSide(:,1);
outSignalxS = filter(y, x, inputSignalxS);
inputSignalyS = preFilterSide(:,2);
outSignalyS = filter(y, x, inputSignalyS);
inputSignalzS = preFilterSide(:,3);
outSignalzS = filter(y, x, inputSignalzS);

%cutting off the beginning of the data due to the filter
filtered = horzcat(outSignalxC,outSignalyC,outSignalzC,outSignalxS,outSignalyS,outSignalzS);
filtered(1:(175/(sqrt(butterFilter))),:) = [];
nerveFiltered = horzcat(outSignalxU,outSignalyU,outSignalzU);
nerveFiltered(1:25,:) = [];
nerveXYZ = nerveFiltered;

%plotting the two points on the arm and the nerve
plot3(filtered(:,1),filtered(:,2),filtered(:,3));
hold on
plot3(filtered(:,4),filtered(:,5),filtered(:,6));
scatter3(nerveFiltered(:,1),nerveFiltered(:,2),nerveFiltered(:,3),'filled');
hold off
title('Filtered Data')
axis square
set(gcf,'position',[100 100 1800 1000]);

end

%path coordinates Function: takes the point data and turns it into a path
%using a moving 4th order spline, and then calculating all path unit
%vectors and radius of curvature/torsion
%inputs: 
%nerveData: the nerve point data from rangefinder rows,3
%outputs:
%eb,en,et: these are the binormal, normal, and tangential vectors for the
%nerve path, one set of vectors for each point
%tau: this is the torsion of the nerve at all points along the path
%rho: this is the radius of curvature at all points along the path

function [eb,en,et,tau,rho] = pathCoords(nerveData)

%user input
D = nerveData(:,1:3);

%intermediate variable declaration
order = 3;
pts = order + 2;
e = 0;

%declarations for the hat vectors
ihat = [double('i'),770];
jhat = [double('j'),770];
khat = [double('k'),770];
ihatstr = char(ihat);
jhatstr = char(jhat);
khatstr = char(khat);

%gets the number of data points for preallocation
sizepath = size(D);
rows = sizepath(1);

%variable preallocation
xpoints = zeros(rows,pts);
ypoints = zeros(rows,pts);
zpoints = zeros(rows,pts);
xto = zeros(rows,order + 1);
yto = zeros(rows,order + 1);
zto = zeros(rows,order + 1);

%time vector
td = zeros(rows,1);

% evaluated S at every point
sf = zeros(rows-(pts-1),1);

%r'(t)
rtp = zeros(rows-(pts-1),3);

%r''(t)
rtdp = zeros(rows-(pts-1),3);

%r'''(t)
rttp = zeros(rows-(pts-1),3);

%tangential vector
et = zeros(rows-(pts-1),3);

%normal vector
en = zeros(rows-(pts-1),3);

%binormal vector
eb = zeros(rows-(pts-1),3);

%rho (radius of curvature)
rho = zeros(rows-(pts-1),1);

%tau (torsion)
tau = zeros(rows-(pts-1),1);

%various intermediate dot products
dot1 = zeros(rows-(pts-1),1);
dot2 = zeros(rows-(pts-1),1);
dot3 = zeros(rows-(pts-1),1);

%intermediate calculations
int1 = zeros(rows-(pts-1),1);
int5 = zeros(rows-(pts-1),3);
int6 = zeros(rows-(pts-1),3);
int7 = zeros(rows-(pts-1),1);

%symbolic variables for calculation
syms t xt yt zt xtp ytp ztp xtdp ytdp ztdp stp xttp yttp zttp;

%defining time vector
for j = 1:rows, td(j) = j/(154/19.6); end

%main data process for loop
for i = 1:(rows-(pts-1))

    %gets the current 5 point snippet of each dimension
    xpoints(i,1:pts) = D(i:i+(pts-1),1);
    ypoints(i,1:pts) = D(i:i+(pts-1),2);
    zpoints(i,1:pts) = D(i:i+(pts-1),3);

    %fits the data to the cubic spline (nth order polynomial)
    xto(i,1:order+1) = polyfit(td(i:i+(pts-1)),xpoints(i,:),order);
    yto(i,1:order+1) = polyfit(td(i:i+(pts-1)),ypoints(i,:),order);
    zto(i,1:order+1) = polyfit(td(i:i+(pts-1)),zpoints(i,:),order);

    %creates the variable form of the equations
    xt(i) = poly2sym(xto(i,1:order+1),t);
    yt(i) = poly2sym(yto(i,1:order+1),t);
    zt(i) = poly2sym(zto(i,1:order+1),t);

    %takes the first, second, and third derivatives of the functions of time 
    xtp(i) = diff(xt(i),t);
    ytp(i) = diff(yt(i),t);
    ztp(i) = diff(zt(i),t);
    xtdp(i) = diff(xtp(i),t);
    ytdp(i) = diff(ytp(i),t);
    ztdp(i) = diff(ztp(i),t);
    xttp(i) = diff(xtdp(i),t);
    yttp(i) = diff(ytdp(i),t);
    zttp(i) = diff(ztdp(i),t);

    %gets the function ds
    stp(i) = ((xtp(i)^2)+(ytp(i)^2)+(ztp(i)^2))^0.5;

    %calculates stp, rtp, rtdp, rttp at the current point
    sf(i,1) = subs(stp(i),td(i+2));
    rtp(i,1) = subs(xtp(i),td(i+2));
    rtp(i,2) = subs(ytp(i),td(i+2));
    rtp(i,3) = subs(ztp(i),td(i+2));
    rtdp(i,1) = subs(xtdp(i),td(i+2));
    rtdp(i,2) = subs(ytdp(i),td(i+2));
    rtdp(i,3) = subs(ztdp(i),td(i+2));
    rttp(i,1) = subs(xttp(i),td(i+2));
    rttp(i,2) = subs(yttp(i),td(i+2));
    rttp(i,3) = subs(zttp(i),td(i+2));

    %calculates et
    et(i,1) = rtp(i,1)/sf(i);
    et(i,2) = rtp(i,2)/sf(i);
    et(i,3) = rtp(i,3)/sf(i);

    %calculates en and rho
    dot1(i) = (rtp(i,1)*rtdp(i,1)) + (rtp(i,2)*rtdp(i,2)) + (rtp(i,3)*rtdp(i,3));
    dot2(i) = (rtdp(i,1)*rtdp(i,1)) + (rtdp(i,2)*rtdp(i,2)) + (rtdp(i,3)*rtdp(i,3));
    int1(i) = ((dot2(i)*(sf(i)^2)) - (dot1(i)^2))^0.5;
    rho(i) = (sf(i)^3)/int1(i);
    int2x = rtdp(i,1)*(sf(i)^2);
    int2y = rtdp(i,2)*(sf(i)^2);
    int2z = rtdp(i,3)*(sf(i)^2);
    int3x = rtp(i,1)*dot1(i);
    int3y = rtp(i,2)*dot1(i);
    int3z = rtp(i,3)*dot1(i);
    int4x = int2x - int3x;
    int4y = int2y - int3y;
    int4z = int2z - int3z;
    en(i,1) = int4x/(sf(i)*int1(i));
    en(i,2) = int4y/(sf(i)*int1(i));
    en(i,3) = int4z/(sf(i)*int1(i));

    % calculates eb and tau
    eb(i,:) = cross(et(i,:),en(i,:));
    int5(i,:) = cross(rtdp(i,:),rttp(i,:));
    dot3(i) = (rtp(i,1)*int5(i,1)) + (rtp(i,2)*int5(i,2)) + (rtp(i,3)*int5(i,3));
    int6(i,:) = cross(rtp(i,:),rtdp(i,:));
    int7(i) = (int6(i,1)^2) + (int6(i,2)^2) + (int6(i,3)^2);
    tau(i) = dot3(i)/int7(i);

end

    % overall prompt while loop
    while e == 0
    
        %asks for data point and stores it in x
        prompt = "data point of interest: ";
        x = input(prompt);
    
        %increments L for row searching
        for l = 1:(rows-(pts-1))
    
            %checks if the given point matches the current row
            if x == l
    
                %prints variables and names and user specified point
                fprintf("\nRadius of Curvature:   %g cm\n\n",rho(l));
                fprintf("Torsion:   %g cm\n\n",tau(l));
                fprintf("Unit Tangential:   %g",et(l,1));
                fprintf(ihatstr);
                fprintf("   %g",et(l,2));
                fprintf(jhatstr);
                fprintf("   %g",et(l,3));
                fprintf(khatstr);
                fprintf("\n\n");
                fprintf("Unit Normal:   %g",en(l,1));
                fprintf(ihatstr);
                fprintf("   %g",en(l,2));
                fprintf(jhatstr);
                fprintf("   %g",en(l,3));
                fprintf(khatstr);
                fprintf("\n\n");
                fprintf("Unit Binormal:   %g",eb(l,1));
                fprintf(ihatstr);
                fprintf("   %g",eb(l,2));
                fprintf(jhatstr);
                fprintf("   %g",eb(l,3));
                fprintf(khatstr);
                fprintf("\n\n");
            end           
        end
    
        %prompts and stores if the user wants another point
        prompt2 = "again?: [y/n]";
        txt = input(prompt2,'s');
    
        %changes the while variable based on answer
        if txt == 'n', e = 1; end
        if txt == 'y', e = 0; end
    end
end

%daqExp function: used to give a live plot of the arm kinematic data for
%use of finding dead zones and demonstration
%Inputs:
%none
%Outputs:
%none

function daqExp()

%declares daq structure d
d = daq('ni');

% Initialize standard transformational matrices
syms theta

%not needed for this particular robot configuration but leaving it in for
%completeness
%Xrot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
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
   
    %plot robot arm data
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

    %graph setup
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

end
