% Load Moco libraries
import org.opensim.modeling.*;
table1=MocoTrajectory('110001_1.sto')
table2=TimeSeriesTable('110001_1_GRF.sto')
time=table1.getTime
aaa=time.getAsMat
time1=aaa*1.1739
table1.setTime(time1)
%opensimTable = TimeSeriesTable(table2);
data = table2.getMatrix.getAsMat; 
timeVec_temp = table2.getIndependentColumn;
for i=1:length(time1)
     timeVec_temp.set(i-1,time1(i))
end
% timeseries=osimTableFromStruct(table1)
% STOFileAdapter.write(timeseries,"145001timescaled.sto") 
table1.write("0800initial.sto")
%table2.wtite("080001timescaled_grf.sto")
STOFileAdapter.write(table2,"0800initial_grf.sto")