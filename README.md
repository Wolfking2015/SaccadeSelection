# SaccadeSelection
## The repository contains the SaccadeSelection function script (SelectSaccade.m),an example eyetrace dataset(EyeX, EyeY), and an example visualization script to illustrate how the function was applied (SaccadeSelectionExample).

### This SaccadeSelection Function is to get saccades interval and gaze interval from raw eye trace based on velocity (Check the program for more details)

### **Input**: 
  
  - *EyeDataX, EyeDataY:* The x,y components of eye trace (in array); Each row is the x/y for one trial; Each column is the x/y value   for one sampling bin;
 - *EyeBinWidth:* The bin size for each eye sample (in ms);
 - *V_Threshold:* Threshold on velocity to screen out the saccades (in degree/second);
 - *ContinueBin:* The minimum number of countinous bins each saccade should last;
 - For example: Suppose the ContinueBin is 5. If one saccade segment which was screened out by exceeding the V_threshold only last for 4 bins continously, it will not be taken into account.
 - *ISI_Threshold:* The minimum inter-saccade-interval allowed for saccades  (in ms). If multiple saccades have intervals within this        ISI_threshold, they will be merged into one.

### **Output**:
 ### - Saccade; a structure array contains all information about the screened out saccades; Each item in the array indicates all saccades screened for each trial.
 - *Saccade(i).NumberOfSaccade:* Number of saccades for trial i;
 - *Saccade(i).SaccadeIntervalIndex:* A n by 2 matrix to store the start bin
                                  and end bin for each saccade for trial i; Each row is for one saccade;
                                  The first colume is the start interval;The second is the end interval;
 - *Saccade(i).Time:* A array matrix to store all the time stamp for each
                  saccade (each row is for one saccade);
 - *Saccade(i).Ecc:* The eccentricity from the origin(0,0); Each row:each saccade;
 - *Saccade(i).X:* The X component for each screened saccade;
 - *Saccade(i).Y:* The Y component for each screened saccade;
 - *Saccade(i).V:* The velocity for each screened saccade;
 - *Saccade(i).A:* The acceloration for each screened saccade;
 - *Saccade(i).PeakV:* The peak velocity  for each screened saccade;
 - *Saccade(i).PeakVTime:* The time stamp reaching the peak velocity  for each screened saccade;
 - *Saccade(i).StartStartPoint:* The start coordinate(x,y) for each saccade;
 - *Saccade(i).StartEndPoint:* The end coordinate(x,y) for each saccade;
 - *Saccade(i).SaccadeVector:* The saccade vector for each saccade (end-start);
 - *Saccade(i).SaccadeAngle:* The angle for each saccade(0 to 359 degree; 0: upward direction, counterclockwise);
 - *Saccade(i).SaccadeAmplitude:* The amplitude for each saccade;
 - *Saccade(i).SaccadeStartTime:* The start time stamp for each saccade;
 - *Saccade(i).SaccadeEndTime:* The end time stamp for each saccade;
 - *Saccade(i).SaccadeDuration:* The lasting duration for each saccade;
