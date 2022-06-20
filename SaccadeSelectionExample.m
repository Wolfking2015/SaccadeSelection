% An example script for Saccade Selection
path = '/Users/yux8/WorkRelated_YXF/Program_Matlab_Local/DataAnalysis/DataFromSystem/PackedData/';%Your path for the eye data
EyeX_file = 'EyeX.mat';
EyeY_file = 'EyeY.mat';
cd(path);
load(EyeX_file);%load EyeX
load(EyeY_file);%load EyeY

% Define eye sample bin and criterium for saccade selection
EyeBinWidth = 1; %in ms

%Criteria for selecting saccades (check the Select Saccade Function for more details)
V_Threshold=40; % in degree/s
ContinueBin=5; % Number of continous bin for each saccde
ISI_Threshold=20; % Inter-spike interval for mulitple saccades

%% Time to select saccades
Saccades=SelectSaccade(EyeX,EyeY,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%Get the information needed from the output
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime_ori=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime_ori=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';


%% The following futher screen and re-arrangements are optional, but I always add them to only include
% saccades over 2 degree, and happend after sometime (not necessary here for the free viewing task).

AmplitudeThreshold=2;
TimeThreshold=0;


SaccadeAngle=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)))',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeAngle,'UniformOutput' ,0));
SaccadeEndPoint=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeEndPoint,'UniformOutput' ,0));
SaccadeStartPoint=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeStartPoint,'UniformOutput' ,0));


%Rearrange the results into matrix (only my personal preference)
SaccadeNumInEachTrial=cellfun(@(x) size(x,2),SaccadeEndPoint);

SaccadeEndPoint_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeEndPoint_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));

SaccadeStartPoint_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeStartPoint_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));


SaccadeEndPoint_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeEndPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeEndPoint_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeEndPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeStartPoint_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeStartPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeStartPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));



SaccadeStartTime=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeStartTime_ori,'UniformOutput' ,0));
SaccadeEndTime=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeEndTime_ori,'UniformOutput' ,0));


SaccadeAmplitude=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeAmplitude,'UniformOutput' ,0));


%% Finally, for some example visualization 

% 1: For Polar plot
% Segment the data into groups
SeparationPoint=[45:45:360];

%Reshape the saccade angle into a line 
SaccadeAngle = reshape(SaccadeAngle,1,numel(SaccadeAngle));
SaccadeAngle = SaccadeAngle(~isnan(SaccadeAngle));
 
SaccadeVectorTuning=PoloarHistGroup(SaccadeAngle,ones(size(SaccadeAngle,2),1),SeparationPoint);
SaccadeVector=SaccadeVectorTuning.Center;
ProOfSaccades=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades

 
PrefVector=SaccadeVectorTuning.PrefVectorNum;
PrefVector_Amp=SaccadeVectorTuning.PrefVector_AmpNum/sum(SaccadeVectorTuning.DataNum)*100;


% Plot
figure(100);
set(gcf,'color','w');
subplot(1,2,1)
polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[ProOfSaccades ProOfSaccades(1)],...
    '-b','LineWidth',3,'Marker','o');
hold on
h_pc = polarplot([PrefVector PrefVector]/180*pi,[0 PrefVector_Amp],'b','LineWidth',3);


set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);


% 2.For grill plot

SaccadeStartPoint_X_Line = reshape(SaccadeStartPoint_X,[],1);
SaccadeStartPoint_Y_Line = reshape(SaccadeStartPoint_Y,[],1);

SaccadeEndPoint_X_Line = reshape(SaccadeEndPoint_X,[],1);
SaccadeEndPoint_Y_Line = reshape(SaccadeEndPoint_Y,[],1);


X_Comp = SaccadeEndPoint_X_Line - SaccadeStartPoint_X_Line;
Y_Comp = SaccadeEndPoint_Y_Line - SaccadeStartPoint_Y_Line;

%Reorder for the plot according to the y component
X_Comp = X_Comp(~isnan(X_Comp));
Y_Comp = Y_Comp(~isnan(Y_Comp));

[val,Index]= sort(Y_Comp,'ascend');
 

X_Comp = X_Comp(Index);
Y_Comp = Y_Comp(Index);


%Grill plot

subplot(1,2,2)

quiver(X_Comp,Y_Comp,'autoscale','off','color','r','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off




%% Some functions used for visualization
function mat=fillinNaNs(cell_ori)
%Function to convert cell to matrix by filling in NaN at the end

maxnumel=max(cellfun(@(x) size(x,2),cell_ori));
maxrow=max(cellfun(@(x) size(x,1),cell_ori));
mat=cell2mat(cellfun(@(x) [x, NaN*ones(max(maxrow-size(x,1),size(x,1)),maxnumel-size(x,2))],cell_ori,'uniform',0));


end


function DataGrouped=PoloarHistGroup(Angle,Data,Edge)
%Function to bin on the polar coordinates
%Xuefei Yu
%Angles: in degree
DataGrouped=[];
DataLine=[];
GroupLine=[];

%Columnwise;
for j=1:length(Edge)
    if j==length(Edge)
          Interval=[Edge(j),Edge(1)];
          SelectAngle=(Angle>=Interval(1) & Angle<359.999)|(Angle>=0 & Angle<Interval(2));
          Center(j)=(Edge(1)+(Edge(j)-360))/2;
            
    else
           Interval=[Edge(j),Edge(j+1)];
           SelectAngle=Angle>=Interval(1) & Angle<Interval(2);
          Center(j)=mean(Interval);
    end
    if sum(SelectAngle)<3
        DataInBin=NaN;
        DataInBin_mean(j,:)=NaN*ones(1,size(Data,2));
        DataInBin_sem(j,:)=NaN*ones(1,size(Data,2));
        DataInBin_Num(j)=sum(SelectAngle);
        DataInBin_Raw{j}=NaN;
       % Center(j)=NaN;
        
    else
    
         DataInBin=Data(SelectAngle,:);
         DataInBin_mean(j,:)=nanmean(DataInBin,1);
         DataInBin_sem(j,:)=nanstd(DataInBin,[],1)./sqrt(sum(~isnan(DataInBin),1));
       
         
         DataInBin_Num(j)=sum(SelectAngle);
          DataInBin_Raw{j}=DataInBin;
         
         %For Stats
         DataLine=[DataLine;DataInBin];
         GroupLine=[GroupLine;repmat(Center(j),DataInBin_Num(j),1)];
        
        
     end  
end
    
    
%Stats via ANOVA

  try  
      ANOVA_p=anova1_group(DataLine,GroupLine,'off');
  catch
      ANOVA_p=NaN;
  end
 
 %Prefer direction via Vector Summation
 for ii=1:size(DataInBin_mean,2)
     try
        [PrefVector(ii), ~, PrefVector_Amp(ii)] = vectorsumAngle(DataInBin_mean(:,ii),  Center, 0*ones(length(Center),1));
     catch
         PrefVector(ii)=NaN;
         PrefVector_Amp(ii)=NaN;
     end
 
 end
 
 %Prefer direction based on condition numbers
 
 
     try
      [PrefVectorNum, ~, PrefVector_AmpNum] = vectorsumAngle(DataInBin_Num,  Center, 0*ones(length(Center),1));
     catch
         PrefVectorNum(ii)=NaN;
         PrefVector_AmpNum(ii)=NaN;
     end
 
 
  
 

DataGrouped.DataMean=DataInBin_mean;
DataGrouped.DataSEM=DataInBin_sem;
DataGrouped.DataNum=DataInBin_Num;
DataGrouped.Center=Center;

DataGrouped.DataLine=DataLine;
DataGrouped.GroupLine=GroupLine;
DataGrouped.DataInBin_Raw=DataInBin_Raw;
DataGrouped.ANOVA_p=ANOVA_p;
DataGrouped.PrefVector=PrefVector;
DataGrouped.PrefVector_Amp=PrefVector_Amp;

DataGrouped.PrefVectorNum=PrefVectorNum;
DataGrouped.PrefVector_AmpNum=PrefVector_AmpNum;




end
