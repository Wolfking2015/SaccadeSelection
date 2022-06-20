function Saccade=SelectSaccade(EyeDataX,EyeDataY,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%This function is to get saccades interval and gaze interval from raw eye trace
%Xuefei Yu 02202020
%Input: 
% EyeDataX, EyeDataY: The x,y components of eye trace (in array); Each row is the x/y for one trial; Each column is the x/y value for one sampling bin;
% EyeBinWidth: The bin size for each eye sample (in ms);
% V_Threshold: Threshold on velocity to screen out the saccades (in degree/second);
% ContinueBin: The minimum number of countinous bins each saccade should last;
% For example: Suppose the ContinueBin is 5. If one saccade segment which was screened out by exceeding the V_threshold only last for 4 bins continously, it will not be taken into account.
% ISI_Threshold: The minimum inter-saccade-interval allowed for saccades
% (in ms). If multiple saccades have intervals within this ISI_threshold, they will be merged into one.

%Output:
% Saccade; a structure array contains all information about the screened out
% saccades; Each item in the array indicates all saccades screened for each trial.
% Saccade(i).NumberOfSaccade: Number of saccades for trial i;
% Saccade(i).SaccadeIntervalIndex:A n by 2 matrix to store the start bin
%                                  and end bin for each saccade for trial i; Each row is for one saccade;
%                                  The first colume is the start interval;The second is the end interval;
% Saccade(i).Time: A array matrix to store all the time stamp for each
%                  saccade (each row is for one saccade);
% Saccade(i).Ecc: The eccentricity from the origin(0,0); Each row:each saccade;
% Saccade(i).X: The X component for each screened saccade;
% Saccade(i).Y: The Y component for each screened saccade;
% Saccade(i).V: The velocity for each screened saccade;
% Saccade(i).A: The acceloration for each screened saccade;
% Saccade(i).PeakV: The peak velocity  for each screened saccade;
% Saccade(i).PeakVTime: The time stamp reaching the peak velocity  for each screened saccade;
% Saccade(i).StartStartPoint: The start coordinate(x,y) for each saccade;
% Saccade(i).StartEndPoint: The end coordinate(x,y) for each saccade;
% Saccade(i).SaccadeVector: The saccade vector for each saccade (end-start);
% Saccade(i).SaccadeAngle: The angle for each saccade(0 to 359 degree; 0: upward direction, counterclockwise);
% Saccade(i).SaccadeAmplitude: The amplitude for each saccade;
% Saccade(i).SaccadeStartTime: The start time stamp for each saccade;
% Saccade(i).SaccadeEndTime: The end time stamp for each saccade;
% Saccade(i).SaccadeDuration: The lasting duration for each saccade;



%Get the eccentricity of the eye trace

Ecc=sqrt(EyeDataX.^2+EyeDataY.^2);

%Get velocity from the Eyedata
Ecc_V=diff(Ecc,1,2)/EyeBinWidth*1000;
Ecc_V=[Ecc_V(:,1) ,Ecc_V];

%Get the acceleration from the Eyedata
Ecc_Acc=diff(Ecc_V,1,2)/EyeBinWidth*1000;
Ecc_Acc=[Ecc_Acc(:,1),Ecc_Acc];

%{
%Find the base line
BaseLineTime=EyeTime>-50 & EyeTime<50;
BaseLineSelect=abs(Ecc_V)<20;

BaseLine=BaseLineTime.*BaseLineSelect;

BaseLine(BaseLine==0)=NaN;

Ecc_Base=Ecc.*BaseLine;
Ecc_Base_reshape=reshape(Ecc_Base,1,numel(Ecc_Base));


[NumInEachBin,Ecc_Base_Distribute]=hist(Ecc_Base_reshape);

Ecc_Base_Mean=Ecc_Base_Distribute(NumInEachBin==max(NumInEachBin));

Ecc_Normed=Ecc-Ecc_Base_Mean;                           


Ecc_A_Base=Ecc_Acc.*BaseLine;

Ecc_A_Base_reshape=reshape(Ecc_A_Base,1,numel(Ecc_A_Base));


[NumInEachBin,Ecc_Base_A_Distribute]=hist(Ecc_A_Base_reshape);

Ecc_A_Base_Mean=Ecc_Base_A_Distribute(NumInEachBin==max(NumInEachBin));
%}

%% First find out the trace crossing the velocity threshold 
%Selection_Ecc=abs(Ecc_Normed)>3;
Selection_V_Low=abs(Ecc_V)>V_Threshold;
%Selection_A_Low=abs(Ecc_Acc)>10000;

%SelectionSaccade_Basic=Selection_Ecc;%.*Selection_V_Low;%.*Selection_A_Low;
SelectionSaccade_Basic=Selection_V_Low;

%% Loop over each trial(row) for more information about the screened saccades
for i=1:size(EyeDataX,1)
    
    currEcc=Ecc(i,:);
    
    currEcc_X=EyeDataX(i,:);
    currEcc_Y=EyeDataY(i,:);
    
    
   if isempty(EyeDataX)
       %If there is no eye data for this trial/block
        Saccade(i).NumOfSaccade=0;
        Saccade(i).SaccadeIntervalIndex=NaN;
    
        Saccade(i).Time=NaN;
        Saccade(i).Ecc=NaN;
    
         Saccade(i).X=NaN;
        Saccade(i).Y=NaN;
        Saccade(i).V=NaN;
        Saccade(i).A=NaN;
        Saccade(i).PeakV= NaN;
        Saccade(i).PeakVTime=NaN;
        Saccade(i).SaccadeStartPoint=NaN;
        Saccade(i).SaccadeEndPoint=NaN;
        Saccade(i).SaccadeVector=NaN;
        Saccade(i).SaccadeAngle=NaN;
        Saccade(i).SaccadeAmplitude=NaN;
        Saccade(i).SaccadeStartTime=NaN;
        Saccade(i).SaccadeEndTime=NaN;
        Saccade(i).SaccadeDuration=NaN;
   else
    
    
        currEcc_V=abs(Ecc_V(i,:));
        currEcc_A=abs(Ecc_Acc(i,:));
    
        % currEyeTime=EyeTime(i,:);
        currEyeTime=(0:length(EyeDataX)-1)*EyeBinWidth;
    
        currSelection_Basic=SelectionSaccade_Basic(i,:)==1;
    
    
    
        SaccadeStartTime=currEyeTime(currSelection_Basic);
        diffStartTime=diff(SaccadeStartTime);
        SeparationPoint=find(diffStartTime>(EyeBinWidth*1.5));
        SeparationPoint=[0,SeparationPoint];
        SaccadeIntervalIndex=[(SeparationPoint(1:end-1)+1)',[SeparationPoint(2:end)]'];
  
        SaccadeIntervalIndex((SaccadeIntervalIndex(:,2)-SaccadeIntervalIndex(:,1))<ContinueBin,:)=[];
        NumberOfSaccade=size(SaccadeIntervalIndex,1);
    
        if NumberOfSaccade>0
    
        %Select out saccade
        SaccadeSelection_Time=repmat(SaccadeStartTime,NumberOfSaccade,1);
 
        SaccadeSelection_Time=SelectMatrixData(SaccadeSelection_Time,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
    
        %Find out the saccade end time
    
        SaccadeStartTime=SaccadeSelection_Time(:,1);
    
        temp=~isnan(SaccadeSelection_Time);
    
        Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(SaccadeSelection_Time, 1));
        SaccadeEndTime=arrayfun(@(x,y) SaccadeSelection_Time(x,y),1:size(SaccadeSelection_Time, 1), Indices)';
    
    
    
    
     
        %Check saccade interval if more than one saccade was found

            if NumberOfSaccade>1
               InterSaccadeInterval=SaccadeStartTime(2:end)-(SaccadeEndTime(1:end-1));
         
                while (InterSaccadeInterval<ISI_Threshold) & (NumberOfSaccade>1)
                    CombinedSaccade=find(InterSaccadeInterval<ISI_Threshold);
                    CombineIndex=CombinedSaccade(1);
             
                    SaccadeStartTime(CombineIndex+1)=SaccadeStartTime(CombineIndex);
                    SaccadeStartTime(CombineIndex)=[];
                    SaccadeEndTime(CombineIndex)=[];
             
            %{
             mat_combine=NaN*ones(size(SaccadeSelection_Time));
             mat_combine(CombineIndex,:)=SaccadeSelection_Time(CombineIndex+1,:);
             
             SaccadeSelection_Time=[SaccadeSelection_Time,mat_combine];
             SaccadeSelection_Time(CombineIndex+1,:)=[];
             %}
             
                    InterSaccadeInterval=SaccadeStartTime(2:end)-(SaccadeEndTime(1:end-1));
                    NumberOfSaccade=NumberOfSaccade-1;
             
             %{
             SaccadeIntervalIndex(CombineIndex+1,1)= SaccadeIntervalIndex(CombineIndex,1);
            SaccadeIntervalIndex(CombineIndex,:)=[];
             %}
           
             
             
                end %End of while
    
 
            end %End of if saccade interval was over 1
     
     
     %Select relavent parameters
     
            SaccadeSelection_Time=fillInterval(SaccadeStartTime,SaccadeEndTime,EyeBinWidth);
            SaccadeIntervalIndex=[ceil(SaccadeStartTime/EyeBinWidth)+1,ceil(SaccadeEndTime/EyeBinWidth)+1];
     
            SaccadeSelection_Ecc=repmat(currEcc,NumberOfSaccade,1);
            SaccadeSelection_Ecc=SelectMatrixData(SaccadeSelection_Ecc,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
    
     
            SaccadeSelection_V=repmat(currEcc_V,NumberOfSaccade,1);
            SaccadeSelection_V=SelectMatrixData(SaccadeSelection_V,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
            SaccadeSelection_X=repmat(currEcc_X,NumberOfSaccade,1);
            SaccadeSelection_X=SelectMatrixData(SaccadeSelection_X,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
            SaccadeSelection_Y=repmat(currEcc_Y,NumberOfSaccade,1);
            SaccadeSelection_Y=SelectMatrixData(SaccadeSelection_Y,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
            SaccadeSelection_A=repmat(currEcc_A,NumberOfSaccade,1);
            SaccadeSelection_A=SelectMatrixData(SaccadeSelection_A,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
    
     %{
     SaccadeSelection_Ecc=repmat(currEcc(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_Ecc=SelectMatrixData(SaccadeSelection_Ecc,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_V=repmat(currEcc_V(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_V=SelectMatrixData(SaccadeSelection_V,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_X=repmat(currEcc_X(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_X=SelectMatrixData(SaccadeSelection_X,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_Y=repmat(currEcc_Y(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_Y=SelectMatrixData(SaccadeSelection_Y,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
     SaccadeSelection_A=repmat(currEcc_A(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_A=SelectMatrixData(SaccadeSelection_A,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     %}
     
    %Find out the peak velocity and time
 
            [PeakV,PeakV_Loc] = max(SaccadeSelection_V,[],2);
     
             PeakV_Time=arrayfun(@(x,y) SaccadeSelection_Time(x,y),1:size(SaccadeSelection_Time, 1), PeakV_Loc');
     
     
    %Find out the saccade amplitude, and anlge
   
            SaccadeStartPoint_X=SaccadeSelection_X(:,1);
            SaccadeStartPoint_Y=SaccadeSelection_Y(:,1);
 
            temp=~isnan(SaccadeSelection_X);
    
            Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(SaccadeSelection_X, 1));
            SaccadeEndPoint_X=arrayfun(@(x,y) SaccadeSelection_X(x,y),1:size(SaccadeSelection_X, 1), Indices)';
            SaccadeEndPoint_Y=arrayfun(@(x,y) SaccadeSelection_Y(x,y),1:size(SaccadeSelection_Y, 1), Indices)'; 
   
            SaccadeVector=[SaccadeEndPoint_X,SaccadeEndPoint_Y]-[SaccadeStartPoint_X,SaccadeStartPoint_Y];
  
    %Angle and amplitude
    
            [SaccadeAngle, SaccadeAmplitude]=LocToAngle(SaccadeVector(:,1),SaccadeVector(:,2));
    
   
   
   
            Saccade(i).NumOfSaccade=NumberOfSaccade;
   
            Saccade(i).SaccadeIntervalIndex=SaccadeIntervalIndex;
    
            Saccade(i).Time=SaccadeSelection_Time;
            Saccade(i).Ecc=SaccadeSelection_Ecc;
    
            Saccade(i).X=SaccadeSelection_X;
            Saccade(i).Y=SaccadeSelection_Y;
            Saccade(i).V=SaccadeSelection_V;
            Saccade(i).A= SaccadeSelection_A;
            Saccade(i).PeakV= PeakV;
            Saccade(i).PeakVTime=PeakV_Time;
            Saccade(i).SaccadeStartPoint=[SaccadeStartPoint_X,SaccadeStartPoint_Y];
            Saccade(i).SaccadeEndPoint=[SaccadeEndPoint_X,SaccadeEndPoint_Y];
            Saccade(i).SaccadeVector=SaccadeVector;
            Saccade(i).SaccadeAngle=SaccadeAngle;
            Saccade(i).SaccadeAmplitude=SaccadeAmplitude;
            Saccade(i).SaccadeStartTime=SaccadeStartTime;
            Saccade(i).SaccadeEndTime=SaccadeEndTime;
            Saccade(i).SaccadeDuration=SaccadeEndTime-SaccadeStartTime;
        else
        %If no saccade was found
            Saccade(i).NumOfSaccade=NumberOfSaccade;
            Saccade(i).SaccadeIntervalIndex=NaN;
    
            Saccade(i).Time=NaN;
            Saccade(i).Ecc=NaN;
    
            Saccade(i).X=NaN;
            Saccade(i).Y=NaN;
            Saccade(i).V=NaN;
            Saccade(i).A=NaN;
            Saccade(i).PeakV= NaN;
            Saccade(i).PeakVTime=NaN;
            Saccade(i).SaccadeStartPoint=NaN;
            Saccade(i).SaccadeEndPoint=NaN;
            Saccade(i).SaccadeVector=NaN;
            Saccade(i).SaccadeAngle=NaN;
            Saccade(i).SaccadeAmplitude=NaN;
            Saccade(i).SaccadeStartTime=NaN;
            Saccade(i).SaccadeEndTime=NaN;
            Saccade(i).SaccadeDuration=NaN;
        
        end %End of if number of saccade above 0
   end %End of if EyeDataX is empty

     
end %End of the trial loop


end
function mat=fillInterval(Start,End,Bin)
MaxCol=max(ceil((End-Start)/Bin)+1);
mat=NaN*ones(size(Start,1),MaxCol);

 for i=1:size(Start,1)
     Number=ceil((End(i,:)-Start(i,:))/Bin)+1;
     mat(i,1:Number)=Start(i,:):Bin:End(i,:);
  
 end



end

function MatrixSelected = SelectMatrixData(Matrix,StartMat,EndMat)
%This function is for matrix data selection from start index to the end
%index for all trials; Xuefei Yu; 08112019
%   First check the dimension consistency
if (size(Matrix,1)~=size(StartMat,1))||(size(Matrix,1)~=size(EndMat,1))||(size(StartMat,1)~=size(EndMat,1))
    error('The rows of the inputs are not consistent!');
end
% Set up a matrix large enough to store all the selections 
Col_Max=max(EndMat-StartMat+1);
RowNum=size(Matrix,1);

MatrixSelected=ones(RowNum,Col_Max)*NaN;


for i=1:RowNum
    SelectLength=EndMat(i)-StartMat(i)+1;
    CurrentMat=Matrix(i,:);
    if StartMat(i)<0
        MatrixSelected(i,1:-StartMat(i)+1)=NaN;
        MatrixSelected(i,-StartMat(i)+1+1:EndMat(i)-StartMat(i)+1)=CurrentMat(1:EndMat(i));
        
    else
    
    MatrixSelected(i,1:SelectLength)=CurrentMat(StartMat(i):EndMat(i));
    end
    
end

end
%This function is to change location into angle
%i.e. cartesian into polar
%Define the upward direction as zero
function [Angle, Amp]=LocToAngle(X,Y);
Angle=NaN*ones(size(X));
Amp=NaN*ones(size(X));

for i=1:size(X,1)

 AngleCurr=atan2d(Y(i,:),X(i,:));%Transform into angle in degree
 AngleCurr=AngleCurr+90+180;
if sum(AngleCurr>=360)
                 AngleCurr(AngleCurr>=360)=AngleCurr(AngleCurr>=360)-360;
end
            
AmpCurr=sqrt(X(i,:).^2+Y(i,:).^2);
Angle(i,:)=AngleCurr;
Amp(i,:)=AmpCurr;

end


end
