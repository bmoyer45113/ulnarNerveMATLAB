nographreal: this code is used for taking data from a 6-DOF Robot arm with a potentiometer on each joint. Its process is that it takes the data, calculates the end positions of each joint, and writes that 3D XYZ point data to an excel sheet. This script is used over the other when you need more than 7.41 Hz of positional data, but you have to graph separately. 
IMPORTANT VARIABLES:
      Inputs:
      SampleTIme: how long you would like to sample for
      DaqRate: controls the Hz of the daq
      Outputs:
      x,y,z: end positions of the 5th joint (end of the arm without the 6th joints rotation)
      xo,yo,zo: end positions of the 6th joint (with the final rotation)

daqexp: this script is used for taking data from a 6-DOF Robot arm with a potentiometer on each joint. Its process is that it takes one point of data, calculated the end positions of each joint, and then graphs the result real time for the specified number of points. This script is used when 7.41 Hz is enough or when you are troubleshooting your system.
IMPORTANT VARIABLES:
      Inputs:
      Time: this specifies the number of points to be taken, at a rate with my particular DAQ of 7.41 Hz
      Outputs:
      x,y,z: end positions of the 5th joint (end of the arm without the 6th joints rotation)
      xo,yo,zo: end positions of the 6th joint (with the final rotation)

pathcoords: This code is used to calculate the Unit tangential, the Unit normal, the unit binormal, the torsion and the radius of curvature of each point of a path. Its process is to take in XYZ positional data, create a running 5 point set, fit a third order spline, take all necessary derivates, and calculate the final values for each running set. It also allows for the user to look at the data of one specified point at a time.  This script is used in conjunction with the data collection scripts above to find all of the path coordinate attributes for each point. 
IMPORTANT VARIABLES:
      Inputs:
      excel spreadsheet of XYZ data points
      Outputs: 
      rho: radius of curvature for each point
      tau: torsion for each point
      et: tangential unit vector for each point
      en: normal unit vector for each point
      eb: binormal unit vector for each point

# ulnarNerveMATLAB
