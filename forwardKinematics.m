clc; clear variables; close all;

% Initialize and define standard transformational matrices
syms theta
Xrot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Yrot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Zrot = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

%read in dataAq data and make it into an array
armDataIn = readtable("xyzDataTest10.xlsx");
thetas = table2array(armDataIn);

%sort robot voltage data
theta1 = thetas(:,1);%1
theta2 = thetas(:,2);%2
theta3 = thetas(:,3);%3
theta4 = thetas(:,4);%4
theta5 = thetas(:,5);%5
theta6 = thetas(:,6);%6

%create size variable for robot data
size1 = size(thetas);
sampleArray = size1(1);

%define intermediate robot variable sizes
FK = zeros(sampleArray*4,4);
x = zeros(sampleArray,1);
y = zeros(sampleArray,1);
z = zeros(sampleArray,1);
xo = zeros(sampleArray,1);
yo = zeros(sampleArray,1);
zo = zeros(sampleArray,1);

%read in the ultrasound data and make it into an array
dataImport = readtable("UN Ultrasound data test 10.xlsx");
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
z6 = 11.0; %11.0 chisel pointer length
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

%concatenates the positions of the center sensor, nerve, and edge sensor
%respectively
delete('xyzDataSnr.xlsx');
writematrix([x,y,z,xo,yo,zo],'xyzDataSnr.xlsx');
delete('xyzDataNrv.xlsx');
writematrix([xu,yu,zu],'xyzDataNrv.xlsx');
