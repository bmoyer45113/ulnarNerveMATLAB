clc; clear; close all;

%importing data
Ba = readtable("xyzDataNrv.xlsx");
Aa = readtable("xyzDataSnr.xlsx");
B = table2array(Ba);
A = table2array(Aa);

                      %user inputs for center
                        startpt = 100;
                        endpt = 570;

%data smoothing and misc variable declaration
ulnSize = size(B);
preFilterCnt = A(startpt:endpt,1:3);
preFilterUln = B(35:ulnSize(1)-40,1:3);
preFilterSide = A(startpt:endpt,4:6);

%nerve filtering

%cutoff frequency for nerve
Fu = 1;

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
F = 1;

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
filtered(1:(200/(sqrt(1))),:) = [];
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