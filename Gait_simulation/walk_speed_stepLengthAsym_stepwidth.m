clear
close all
clc
% Load Moco libraries
import org.opensim.modeling.*;
%% tracking reference, initial guess, and solution file names
reference_track="08000normalsecond"
initialguess="08015normal5second"
solutionname="newCodetest"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Options and settings
%% asymmetry target and weight
asymmetrygoal=1    % set this to 1 for applying step length asymmetry goal
targetasym=0.15;   % asymmetry level 
asyweight=100;     % weight of the asymmetry goal
strideLen=0.74;    % this should be assigned when applying step length asymmetry the value depends on the walking speed

%% average speed
AvgSpeed = 0.8;    % Average walking speed in speed constraint (m/s)      
stepwidth=0;        % set this to 1 for applying step width constraint 
stpwidth=0.2;      % step width (m)
SolveTrackingProblem = 1; % Set to 1 to run an optimization
w = (50.0/95)-0.03;              % Weight on control effort term in cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% making the report file:
data0(1,1)=solutionname
data0(1,2)=solutionname;
data0(2,1)="track_ref"
data0(2,2)=reference_track;
data0(3,1)="initial_guess";
data0(3,2)=initialguess;
data0(4,1)="avg speed";
data0(4,2)=AvgSpeed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the motion tracking problem
track = MocoTrack();
track.setName('gaitTracking');            
tableStatesProcessor = TableProcessor(strcat(reference_track,'.sto'));   % Target kinematic states for tracking
tableStatesProcessor.append(TabOpLowPassFilter(10));                     % Smooth the target kinematic data (helpful for shoulder and elbow angles, which were manually digitized)
modelProcessor = ModelProcessor('Rajagopal2015_torque1.osim');                    % Load the OpenSim model
modelProcessor.append(ModOpTendonComplianceDynamicsModeDGF('implicit')); % Use muscle contractile dynamics
%modelProcessor.append(ModOpIgnorePassiveFiberForcesDGF());              % Set passive muscle fiber forces to zero
track.setModel(modelProcessor);                                          % Apply the model to the tracking problem
track.setStatesReference(tableStatesProcessor);                          % Apply the target state data to the tracking problem
track.set_states_global_tracking_weight(1.0);                            % Default tracking weight (is changed below)
track.set_allow_unused_references(true);                                 % Target data can include DoF not in this model
track.set_track_reference_position_derivatives(true);                    % Track speed trajectories
track.set_apply_tracked_states_to_guess(true);                           % Use target data if generating a new initial guess
track.set_initial_time(0.0);                                             % Initial time [s]
track.set_final_time(1.00);                                              % Final time [s], comment out to optimize movement duration

% Specify tracking weights as standard deviations averaged over the gait cycle from Miller et al. (2014)
pq = 1.0/29; % Divided by 29 because there are 23 tracking targets 6 GRF components 

stateWeights = MocoWeightSet();
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tx/value',       pq/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_ty/value',       pq/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tz/value',       pq/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tilt/value',     pq/(1*0.0585)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_list/value',     pq/(1*0.0254)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_rotation/value', pq/(1*0.0474)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_extension/value',            0/(1*0.1745)^2)); % Global trunk posture is tracking in another goal
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_bend/value',           0/(1*0.1745)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_rotation/value',           0/(1*0.1745)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_flexion_r/value',          pq/(1*0.0647)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_adduction_r/value',        pq/(1*0.0410)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_rotation_r/value',         pq/(3*0.0728)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/knee_r/knee_angle_r/value',          pq/(1*0.0889)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/ankle_r/ankle_angle_r/value',        pq/(1*0.0574)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/subtalar_r/subtalar_angle_r/value',  pq/(3*0.0595)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/mtp_r/mtp_angle_r/value',            pq/(3*0.0873)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_flexion_l/value',          pq/(1*0.0647)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_adduction_l/value',        pq/(1*0.0410)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_rotation_l/value',         pq/(3*0.0728)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/knee_l/knee_angle_l/value',          pq/(1*0.0889)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/ankle_l/ankle_angle_l/value',        pq/(1*0.0574)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/subtalar_l/subtalar_angle_l/value',  pq/(3*0.0595)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/mtp_l/mtp_angle_l/value',            pq/(3*0.0873)^2));
pw = pq*0.0001; % Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tx/speed',       pw/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_ty/speed',       pw/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tz/speed',       pw/(3*0.1000)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_tilt/speed',     pw/(1*0.0585)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_list/speed',     pw/(1*0.0254)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/groundPelvis/pelvis_rotation/speed', pw/(1*0.0474)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_extension/speed',            0/(1*0.1745)^2)); % Trunk velocities are not tracked;
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_bend/speed',           0/(1*0.1745)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/back/lumbar_rotation/speed',           0/(1*0.1745)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_flexion_r/speed',          pw/(1*0.0647)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_adduction_r/speed',        pw/(1*0.0410)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_r/hip_rotation_r/speed',         pw/(3*0.0728)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/knee_r/knee_angle_r/speed',          pw/(1*0.0889)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/ankle_r/ankle_angle_r/speed',        pw/(1*0.0574)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/subtalar_r/subtalar_angle_r/speed',  pw/(3*0.0595)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/mtp_r/mtp_angle_r/speed',            pw/(3*0.0873)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_flexion_l/speed',          pw/(1*0.0647)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_adduction_l/speed',        pw/(1*0.0410)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/hip_l/hip_rotation_l/speed',         pw/(3*0.0728)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/knee_l/knee_angle_l/speed',          pw/(1*0.0889)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/ankle_l/ankle_angle_l/speed',        pw/(1*0.0574)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/subtalar_l/subtalar_angle_l/speed',  pw/(3*0.0595)^2));
stateWeights.cloneAndAppend(MocoWeight('/jointset/mtp_l/mtp_angle_l/speed',            pw/(3*0.0873)^2));
track.set_states_weight_set(stateWeights);
% Define the Moco study and problem
study = track.initialize();
problem = study.updProblem();

% Define the periodicity goal
periodicityGoal = MocoPeriodicityGoal('symmetryGoal');
problem.addGoal(periodicityGoal);
model = modelProcessor.process();
model.initSystem();

% All states are periodic except pelvis anterior-posterior translation
for i = 1:model.getNumStateVariables()
    currentStateName = string(model.getStateVariableNames().getitem(i-1));
    if (~contains(currentStateName,'pelvis_tx/value'))
       periodicityGoal.addStatePair(MocoPeriodicityGoalPair(currentStateName));
    end
end

% All controls are periodic
for i = 1:model.getNumControls()
    currentControlName = string(problem.createRep().createControlInfoNames().get(i-1));
    periodicityGoal.addControlPair(MocoPeriodicityGoalPair(currentControlName));
end

% Prescribed average walking speed
speedGoal = MocoAverageSpeedGoal('speed');
speedGoal.set_desired_average_speed(AvgSpeed);
problem.addGoal(speedGoal);

%applying step width goal and Preventing feet penetration
distanceConstraint = MocoFrameDistanceConstraint();
if stepwidth==1
distanceConstraint.setName("stepWidth")
distanceConstraint.setProjection("vector")
distanceConstraint.setProjectionVector(Vec3(0,0,1));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/calcn_l','/bodyset/calcn_r',stpwidth,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/toes_r', stpwidth,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/calcn_l','/bodyset/toes_r', stpwidth,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/calcn_r',stpwidth,100));
else
distanceConstraint.setName('distance_constraint');
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/calcn_l','/bodyset/calcn_r',0.10,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/toes_r', 0.10,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/calcn_l','/bodyset/toes_r', 0.10,100));
distanceConstraint.addFramePair(MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/calcn_r',0.10,100));
end
problem.addPathConstraint(distanceConstraint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define upright torso goal
%changed the weight from /37 to /29 because i removed several tracking
%terms and i want it to be consistent with the tracking weight reduction
torsoGoal = MocoOrientationTrackingGoal('torsoGoal',(3/(0.0873^2))/29);
torsoTable = TableProcessor('modified_torso.sto');
torsoGoal.setStatesReference(torsoTable);
paths = StdVectorString();
paths.add('/bodyset/torso')
torsoGoal.setFramePaths(paths);
problem.addGoal(torsoGoal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control effort term (minimize squared muscle excitations)
effort = MocoControlGoal.safeDownCast(problem.updGoal('control_effort'));
effort.setWeight(w);
effort.setExponent(2);
effort.setDivideByDisplacement(false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   STEP LENGTH ASYMMETRY GOAL    
stepLengthAsymmetry = MocoStepLengthAsymmetryGoal();
stepLengthAsymmetry.setWeight(asyweight);
% Provide the body name for the right foot.
stepLengthAsymmetry.setRightFootFrame('/bodyset/calcn_r');
% Provide the body name for the left foot.
stepLengthAsymmetry.setLeftFootFrame('/bodyset/calcn_l');  
% Value for smoothing term use to compute asymmetry (default is 5). Users may 
% need to adjust this based on convergence and matching the target
% asymmetry.
% Larger smoothing values mean that larger step length errors are required
% for a given step to be counted towards the total asymmetry error minimized
% in the cost.
stepLengthAsymmetry.setAsymmetrySmoothing(4); 
% Target step length asymmetry: positive numbers mean greater right step lengths 
% than left.
stepLengthAsymmetry.setTargetAsymmetry(targetasym);
% Provide the stride length. This in combination with the average walking
% speed determines the stride time.
stepLengthAsymmetry.setStrideLength(strideLen);
% Add  the goal to problem.
if asymmetrygoal==1
   problem.addGoal(stepLengthAsymmetry);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRF tracking
contactTracking = MocoContactTrackingGoal('contact',6/(75.337*9.81)/(0.1073^2)/29);
contactTracking.setExternalLoadsFile('refGRF_3D.xml');
forceNamesRightFoot = StdVectorString();
forceNamesRightFoot.add('/forceset/contactHeel_r');
forceNamesRightFoot.add('/forceset/contactMH1_r');
forceNamesRightFoot.add('/forceset/contactMH3_r');
forceNamesRightFoot.add('/forceset/contactMH5_r');
forceNamesRightFoot.add('/forceset/contactOtherToes_r');
forceNamesRightFoot.add('/forceset/contactHallux_r');
trackRightGRF = MocoContactTrackingGoalGroup(forceNamesRightFoot,'RightGRF');
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r');
contactTracking.addContactGroup(trackRightGRF);
forceNamesLeftFoot = StdVectorString();
forceNamesLeftFoot.add('/forceset/contactHeel_l');
forceNamesLeftFoot.add('/forceset/contactMH1_l');
forceNamesLeftFoot.add('/forceset/contactMH3_l');
forceNamesLeftFoot.add('/forceset/contactMH5_l');
forceNamesLeftFoot.add('/forceset/contactOtherToes_l');
forceNamesLeftFoot.add('/forceset/contactHallux_l');
trackLeftGRF = MocoContactTrackingGoalGroup(forceNamesLeftFoot,'LeftGRF');
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l');
contactTracking.addContactGroup(trackLeftGRF);
problem.addGoal(contactTracking);
%% Bound constraints
%problem.setTimeBounds(0, [0.95, 1.05]); % Uncomment if optimizing movement duration
problem.setStateInfo('/jointset/groundPelvis/pelvis_tilt/value', [-10*pi/180, 10*pi/180]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_list/value', [-10*pi/180, 10*pi/180]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_rotation/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_tx/value', [-0.2, 2]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_ty/value', [0.75, 1.25]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_tz/value', [-0.25, 0.25]);
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25*pi/180, 45*pi/180]);
problem.setStateInfo('/jointset/hip_l/hip_adduction_l/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/hip_l/hip_rotation_l/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25*pi/180, 45*pi/180]);
problem.setStateInfo('/jointset/hip_r/hip_adduction_r/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/hip_r/hip_rotation_r/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-75*pi/180, 0]);
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-75*pi/180, 0]);
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-30*pi/180, 20*pi/180]);
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-30*pi/180, 20*pi/180]);
problem.setStateInfo('/jointset/subtalar_l/subtalar_angle_l/value', [-20*pi/180, 20*pi/180]);
problem.setStateInfo('/jointset/subtalar_r/subtalar_angle_r/value', [-20*pi/180, 20*pi/180]);
problem.setStateInfo('/jointset/mtp_l/mtp_angle_l/value', [-20*pi/180, 70*pi/180]);
problem.setStateInfo('/jointset/mtp_r/mtp_angle_r/value', [-20*pi/180, 70*pi/180]);
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-20*pi/180, 20*pi/180]);
problem.setStateInfo('/jointset/back/lumbar_bending/value', [-15*pi/180, 15*pi/180]);
problem.setStateInfo('/jointset/back/lumbar_rotation/value', [-15*pi/180, 15*pi/180]);

%problem.setStateInfoPattern('/forceset/.*/normalized_tendon_force', [0, 1.8], [], []);
problem.setStateInfoPattern('/forceset/.*/activation',   [0.001, 1.0], [], []);
problem.setStateInfoPattern('/forceset/lumbar_ext/activation',   [-1.0, 1.0], [], []);
problem.setStateInfoPattern('/forceset/lumbar_bend/activation',   [-1.0, 1.0], [], []);
problem.setStateInfoPattern('/forceset/lumbar_rot/activation',   [-1.0, 1.0], [], []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the solver and set its options
solver = MocoCasADiSolver.safeDownCast(study.updSolver());
%solver.set_multibody_dynamics_mode('implicit')
%solver.set_minimize_implicit_multibody_accelerations(true)
%solver.set_implicit_multibody_accelerations_weight(0.000001)
solver.set_optim_max_iterations(3000);
solver.set_num_mesh_intervals(50);
solver.set_optim_constraint_tolerance(1e-04);
solver.set_optim_convergence_tolerance(1e+01);
solver.set_minimize_implicit_auxiliary_derivatives(true)
solver.set_implicit_auxiliary_derivatives_weight(0.00001)
solver.set_enforce_constraint_derivatives(true)
%solver.set_lagrange_multiplier_weight(100.0);
solver.resetProblem(problem);

%% Set the normalized tendon forces if not loading initial guess from file
%guess = solver.createGuess();
% numRows = guess.getNumTimes();
% StateNames = model.getStateVariableNames();
% for i = 1:model.getNumStateVariables()
%     currentStateName = string(StateNames.getitem(i-1));
%     if contains(currentStateName,'normalized_tendon_force')
%         MusName = currentStateName;
%         guess.setState(currentStateName, linspace(0.2,0.2,numRows));
%     end
% end
  
% Uncomment this line if not loading an initial guess
%solver.setGuess(guess);

%% Load the initial guess from this file
solver.setGuessFile(strcat(initialguess,'.sto'));
%% other options for setting the initial guess
%% Insert previous solution from a less complex model
%guess = solver.createGuess();
%prevSolution = MocoTrajectory('30LAsymAless51.sto');
%prevStatesTable = prevSolution.exportToStatesTable();
%prevControlsTable = prevSolution.exportToControlsTable();
%guess.insertStatesTrajectory(prevStatesTable, true);
%guess.insertControlsTrajectory(prevControlsTable, true);
%solver.setGuess(guess);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the problem
if (SolveTrackingProblem == 1)
    gaitTrackingSolution = study.solve();
    % Write the solution to a file
    gaitTrackingSolution.write(strcat(solutionname,'.sto'));
    % Write solution's GRF to a file
    externalForcesTableFlat = opensimMoco.createExternalLoadsTableForGait(model,gaitTrackingSolution,forceNamesRightFoot,forceNamesLeftFoot);
    STOFileAdapter.write(externalForcesTableFlat, strcat(solutionname,'_GRF.sto'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Breakdown of terms in objective function
    disp('   ')
    disp('Breakdown of objective (including weights):')
    costTerms = StdVectorString();
    for i = 1:gaitTrackingSolution.getNumObjectiveTerms()
        termName = gaitTrackingSolution.getObjectiveTermNames().get(i-1);
        termNamestr = char(termName);
        termNamestr=convertCharsToStrings(termNamestr);
        row(1,i)=termNamestr;
        row(2,i)=num2str(gaitTrackingSolution.getObjectiveTermByIndex(i-1));
    end
    row
   
else
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reporting the asymmetry levels acheived 
[stepTimeAsymmetry,stepLengthAsymmetry,rightStepLength,leftStepLength] = computeStepAsymmetryValues(...
    strcat(solutionname,'.sto'), ...
    strcat(solutionname,'_GRF.sto'));
fprintf('\n')
fprintf(['Step Time Asymmetry = ' num2str(stepTimeAsymmetry, '%5.3f') '%%'])
fprintf('\n')
fprintf(['Step Length Asymmetry = ' num2str(stepLengthAsymmetry, '%5.3f') '%%'])
fprintf('\n')
fprintf(['Right Step Length = ' num2str(rightStepLength, '%5.3f') '%'])
fprintf('\n')
fprintf(['left Step Length = ' num2str(leftStepLength, '%5.3f') '%'])
fprintf('\n')

data1(1,1)="target asym";
data1(1,2)=targetasym;
data1(2,1)="asym_weight";
data1(2,2)=asyweight;
data1(3,1)="Step Time Asymmetry";
data1(3,2)=stepTimeAsymmetry;
data1(4,1)="stepLengthAsymmetry";
data1(4,2)=stepLengthAsymmetry;
data1(5,1)="Right Step Length";
data1(5,2)=rightStepLength;
data1(6,1)="left Step Length";
data1(6,2)=leftStepLength;
data1(7,1)="num_iter";
data1(7,2)=gaitTrackingSolution.getNumIterations;
data1(8,1)="solver duration";
data1(8,2)=gaitTrackingSolution.getSolverDuration;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Muscle energy expenditure (Umberger et al., 2003; Koelewijn et al., 2018)
outputPaths = StdVectorString();
outputPaths.add('.*fiber_velocity');
outputPaths.add('.*activation');
outputPaths.add('.*active_fiber_force');
outputPaths.add('.*fiber_length');
outputPaths.add('.*active_force_length_multiplier');
outputTable = study.analyze(gaitTrackingSolution, outputPaths);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OFASym=row';
OFASym1=[data0;OFASym];
xlswrite(strcat(solutionname,'.xlsx'),OFASym1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize the solution
study.visualize(gaitTrackingSolution);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  function : computeStepAsymmetryValues  

function [stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength] = computeStepAsymmetryValues(...
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


