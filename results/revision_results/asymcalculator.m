clear
close all
clc

% Compute the actual step time asymmetry.
[stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength,pelvisSpeed] = computeStepAsymmetryValues(...
    '14500normalsecond.sto', ...
    '14500normalsecond_GRF.sto');
fprintf('\n')
fprintf(['Step Time Asymmetry = ' num2str(stepTimeAsymmetry, '%4.2f') '%%'])
fprintf('\n')
fprintf(['Step Length Asymmetry = ' num2str(stepLengthAsymmetry, '%4.2f') '%%'])
fprintf('\n')
fprintf(['Right Step Length = ' num2str(rightStepLength, '%4.2f') '%'])
fprintf('\n')
fprintf(['left Step Length = ' num2str(leftStepLength, '%4.2f') '%'])
fprintf('\n')
%STOFileAdapter.write(pelvisSpeed, 'pelvisSpeed.sto')
opensimMoco.writeTableToFile(pelvisSpeed, 'pelvisSpeed.sto')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    function : computeStepAsymmetryValues  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength, pelvisSpeed] = computeStepAsymmetryValues(...
            solutionFile, grfsFile)
import org.opensim.modeling.*;
model = Model('Rajagopal2015_torque1.osim');
solution = TimeSeriesTable(solutionFile);
grfs = TimeSeriesTable(grfsFile);
% Get time vector
nrow = grfs.getNumRows();
timeVec = grfs.getIndependentColumn();
time = zeros(nrow, 1);
for i = 1:nrow
    time(i) = timeVec.get(i-1);
end
% Find the time of the left and right heelstrikes
contactForceThreshold = 25; % N
rightVerticalGRF = grfs.getDependentColumn('ground_force_r_vy').getAsMat();
leftVerticalGRF = grfs.getDependentColumn('ground_force_l_vy').getAsMat();
rightHeelStrikeIndex = findHeelStrikeIndex(rightVerticalGRF, ... 
    contactForceThreshold);
leftHeelStrikeIndex = findHeelStrikeIndex(leftVerticalGRF, ...
    contactForceThreshold);
rightHeelStrike = time(rightHeelStrikeIndex);
leftHeelStrike = time(leftHeelStrikeIndex);
% Compute step time asymmetry
if rightHeelStrike < leftHeelStrike
    leftStepTime = leftHeelStrike - rightHeelStrike;
    rightStepTime = time(end) - leftHeelStrike + rightHeelStrike;
else 
    rightStepTime = rightHeelStrike - leftHeelStrike;
    leftStepTime = time(end) - rightHeelStrike + leftHeelStrike;
end
stepTimeAsymmetry = ((rightStepTime - leftStepTime) / ... 
                    (rightStepTime + leftStepTime)) * 100.0;
% Create StatesTrajectory from solution             
statesTraj = StatesTrajectory().createFromStatesTable(model, solution, ... 
    false, true, true);

stateRHS = statesTraj.get(rightHeelStrikeIndex-1);
stateLHS = statesTraj.get(leftHeelStrikeIndex-1);
                
rightStepLength = computeStepLength(model, stateRHS);
leftStepLength = computeStepLength(model, stateLHS);

stepLengthAsymmetry = ((rightStepLength - leftStepLength) / ...
                       (rightStepLength + leftStepLength)) * 100.0;
bodySet = model.getBodySet()
pelvis = bodySet.get('pelvis')
pelvisSpeed = pelvis.getOutput('velocity')
                 
end
function [index] = findHeelStrikeIndex(verticalGRF, forceThreshold)

contactIndices = find(verticalGRF > forceThreshold);
nonContactIndices = find(verticalGRF < forceThreshold);

if nonContactIndices(1) > 1
    index = nonContactIndices(end) + 1;
else
    index = contactIndices(1);
end

if index > length(verticalGRF)
    index = 1;
end

end

function [stepLength] = computeStepLength(model, state)

model.initSystem();
model.realizePosition(state);

leftContactGeometry = model.getContactGeometrySet.get('heel_r');
rightContactGeometry = model.getContactGeometrySet.get('heel_l');

rightHeelPosition = leftContactGeometry.getFrame().getPositionInGround(state);
leftHeelPosition = rightContactGeometry.getFrame().getPositionInGround(state);

stepLength = abs(rightHeelPosition.get(0) - leftHeelPosition.get(0));

end




