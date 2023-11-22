clc; clear; close all;

%importing data
A = readtable("endaffector5.xlsx");

                            %user input
                           B = A(:,1:3);

%4th order butterworth filtering on each dimension
C = table2array(B);
F = 3;
Fs = 40;
[y, x] = butter(4, F/(Fs/2));
inputSignalx = C(:,1);
outSignalx = filter(y, x, inputSignalx);
inputSignaly = C(:,2);
outSignaly = filter(y, x, inputSignaly);
inputSignalz = C(:,3);
outSignalz = filter(y, x, inputSignalz);

%concatenation of filtered data and removal of startup data
D = horzcat(outSignalx,outSignaly,outSignalz);
D(1:50,:) = [];

%plotting of both filtered and unfiltered data for observation
subplot(121)
plot3(D(:,1),D(:,2),D(:,3));
title('Filtered Data')
subplot(122)
plot3(C(:,1),C(:,2),C(:,3));
title('Raw Data')
set(gcf,'position',[450 400 1000 600]);

%intermediate variable declaration
order = 3;
pts = order + 2;
e = 0;
ihat = [double('i'),770];
jhat = [double('j'),770];
khat = [double('k'),770];
ihatstr= char(ihat);
jhatstr = char(jhat);
khatstr = char(khat);

%gets the number of data points for preallocation
size = size(D);
rows = size(1);

%variable preallocation
xpoints = zeros(rows,pts);
ypoints = zeros(rows,pts);
zpoints = zeros(rows,pts);
xto = zeros(rows,order + 1);
yto = zeros(rows,order + 1);
zto = zeros(rows,order + 1);
td = zeros(rows,1);
sf = zeros(rows-(pts-1),1);
rtp = zeros(rows-(pts-1),3);
rtdp = zeros(rows-(pts-1),3);
rttp = zeros(rows-(pts-1),3);
et = zeros(rows-(pts-1),3);
en = zeros(rows-(pts-1),3);
eb = zeros(rows-(pts-1),3);
rho = zeros(rows-(pts-1),1);
tau = zeros(rows-(pts-1),1);
dot1 = zeros(rows-(pts-1),1);
dot2 = zeros(rows-(pts-1),1);
dot3 = zeros(rows-(pts-1),1);
int1 = zeros(rows-(pts-1),1);
int5 = zeros(rows-(pts-1),3);
int6 = zeros(rows-(pts-1),3);
int7 = zeros(rows-(pts-1),1);

syms t xt yt zt xtp ytp ztp xtdp ytdp ztdp stp xttp yttp zttp;

%defining time vector
for j = 1:rows, td(j) = 0.135*j; end

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

    %takes the first and second derivatives of the functions of time 
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

    %calculates et, en, rho at the current point
    et(i,1) = rtp(i,1)/sf(i);
    et(i,2) = rtp(i,2)/sf(i);
    et(i,3) = rtp(i,3)/sf(i);
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

    % takes the cross product to find the binormal vector
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

            %prints variables and names
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

