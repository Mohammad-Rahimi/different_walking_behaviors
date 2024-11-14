import opensim as osim
import numpy as np
import pandas as pd

# Tracking reference, initial guess, and solution file names
reference_track = "refQ_31DoF"
initialguess = "145151"
solutionname = "145301"

# Options and settings
asymmetrygoal = 1  # Set this to 1 for applying step length asymmetry goal
targetasym = 0.3  # Asymmetry level
asyweight = 175  # Weight of the asymmetry goal
strideLen = 1.45  # This should be assigned when applying step length asymmetry the value depends on the walking speed

# Average speed
AvgSpeed = 1.45  # Average walking speed in speed constraint (m/s)
study_final_time = 1.0
track_final_time = 1.0
stepwidth = 0  # Set this to 1 for applying step width constraint
stpwidth = 0.1  # Step width (m)
SolveTrackingProblem = 1  # Set to 1 to run an optimization
w = (50.0 / 95) - 0.03  # Weight on control effort term in cost function


def computeStepAsymmetryValues(solutionFile, grfsFile):
    import opensim as osim
    model = osim.Model('Rajagopal2015_torque1.osim')
    solution = osim.TimeSeriesTable(solutionFile)
    grfs = osim.TimeSeriesTable(grfsFile)
    
    # Get time vector
    nrow = grfs.getNumRows()
    timeVec = grfs.getIndependentColumn()
    time = np.array([timeVec[i] for i in range(nrow)])
    
    # Find the time of the left and right heelstrikes
    contactForceThreshold = 25  # N
    rightVerticalGRF = grfs.getDependentColumn('ground_force_r_vy').to_numpy()
    leftVerticalGRF = grfs.getDependentColumn('ground_force_l_vy').to_numpy()
    rightHeelStrikeIndex = findHeelStrikeIndex(rightVerticalGRF, contactForceThreshold)
    leftHeelStrikeIndex = findHeelStrikeIndex(leftVerticalGRF, contactForceThreshold)
    rightHeelStrike = time[rightHeelStrikeIndex]
    leftHeelStrike = time[leftHeelStrikeIndex]
    
    # Compute step time asymmetry
    if rightHeelStrike < leftHeelStrike:
        leftStepTime = leftHeelStrike - rightHeelStrike
        rightStepTime = time[-1] - leftHeelStrike + rightHeelStrike
    else:
        rightStepTime = rightHeelStrike - leftHeelStrike
        leftStepTime = time[-1] - rightHeelStrike + leftHeelStrike
    stepTimeAsymmetry = ((rightStepTime - leftStepTime) / (rightStepTime + leftStepTime)) * 100.0
    
    # Create StatesTrajectory from solution
    statesTraj = osim.StatesTrajectory().createFromStatesTable(model, solution, False, True, True)
    stateRHS = statesTraj[int(rightHeelStrikeIndex)]
    stateLHS = statesTraj[int(leftHeelStrikeIndex)]
    
    rightStepLength = computeStepLength(model, stateRHS)
    leftStepLength = computeStepLength(model, stateLHS)
    
    stepLengthAsymmetry = ((rightStepLength - leftStepLength) / (rightStepLength + leftStepLength)) * 100.0
    
    return stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength

def findHeelStrikeIndex(verticalGRF, forceThreshold):
    contactIndices = np.where(verticalGRF > forceThreshold)[0]
    nonContactIndices = np.where(verticalGRF < forceThreshold)[0]
    
    if nonContactIndices[0] > 1:
        index = nonContactIndices[-1] + 1
    else:
        index = contactIndices[0]
    
    if index >= len(verticalGRF):
        index = 0
    
    return index

def computeStepLength(model, state):
    model.initSystem()
    model.realizePosition(state)
    
    leftContactGeometry = model.getContactGeometrySet().get('heel_r')
    rightContactGeometry = model.getContactGeometrySet().get('heel_l')
    
    rightHeelPosition = leftContactGeometry.getFrame().getPositionInGround(state)
    leftHeelPosition = rightContactGeometry.getFrame().getPositionInGround(state)
    
    stepLength = abs(rightHeelPosition.get(0) - leftHeelPosition.get(0))
    
    return stepLength
# Define the motion tracking problem
track = osim.MocoTrack()
track.setName('gaitTracking')
tableStatesProcessor = osim.TableProcessor(reference_track + '.sto')  # Target kinematic states for tracking
tableStatesProcessor.append(osim.TabOpLowPassFilter(10))  # Smooth the target kinematic data
modelProcessor = osim.ModelProcessor('Rajagopal2015_torque1.osim')  # Load the OpenSim model
modelProcessor.append(osim.ModOpTendonComplianceDynamicsModeDGF('implicit'))  # Use muscle contractile dynamics
# modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())  # Set passive muscle fiber forces to zero
track.setModel(modelProcessor)  # Apply the model to the tracking problem
track.setStatesReference(tableStatesProcessor)  # Apply the target state data to the tracking problem
track.set_states_global_tracking_weight(1.0)  # Default tracking weight
track.set_allow_unused_references(True)  # Target data can include DoF not in this model
track.set_track_reference_position_derivatives(True)  # Track speed trajectories
track.set_apply_tracked_states_to_guess(True)  # Use target data if generating a new initial guess
track.set_initial_time(0.0)  # Initial time [s]
track.set_final_time(track_final_time)  # Final time [s]

# Specify tracking weights as standard deviations averaged over the gait cycle from Miller et al. (2014)
pq = 1.0 / 29  # Divided by 29 because there are 23 tracking targets 6 GRF components

stateWeights = osim.MocoWeightSet()
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tx/value', pq / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_ty/value', pq / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tz/value', pq / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tilt/value', pq / (1 * 0.0585) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_list/value', pq / (1 * 0.0254) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_rotation/value', pq / (1 * 0.0474) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/value', 0 / (1 * 0.1745) ** 2))  # Global trunk posture is tracking in another goal
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_bend/value', 0 / (1 * 0.1745) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_rotation/value', 0 / (1 * 0.1745) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/value', pq / (1 * 0.0647) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_adduction_r/value', pq / (1 * 0.0410) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_rotation_r/value', pq / (3 * 0.0728) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/value', pq / (1 * 0.0889) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/value', pq / (1 * 0.0574) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_r/subtalar_angle_r/value', pq / (3 * 0.0595) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/value', pq / (3 * 0.0873) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/value', pq / (1 * 0.0647) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_adduction_l/value', pq / (1 * 0.0410) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_rotation_l/value', pq / (3 * 0.0728) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/value', pq / (1 * 0.0889) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/value', pq / (1 * 0.0574) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_l/subtalar_angle_l/value', pq / (3 * 0.0595) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/value', pq / (3 * 0.0873) ** 2))
pw = pq * 0.0001  # Scale the generalized speed tracking errors by this constant
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tx/speed', pw / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_ty/speed', pw / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tz/speed', pw / (3 * 0.1000) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_tilt/speed', pw / (1 * 0.0585) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_list/speed', pw / (1 * 0.0254) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/groundPelvis/pelvis_rotation/speed', pw / (1 * 0.0474) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_extension/speed', 0 / (1 * 0.1745) ** 2))  # Trunk velocities are not tracked
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_bend/speed', 0 / (1 * 0.1745) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/back/lumbar_rotation/speed', 0 / (1 * 0.1745) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_flexion_r/speed', pw / (1 * 0.0647) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_adduction_r/speed', pw / (1 * 0.0410) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_r/hip_rotation_r/speed', pw / (3 * 0.0728) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_r/knee_angle_r/speed', pw / (1 * 0.0889) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_r/ankle_angle_r/speed', pw / (1 * 0.0574) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_r/subtalar_angle_r/speed', pw / (3 * 0.0595) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_r/mtp_angle_r/speed', pw / (3 * 0.0873) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_flexion_l/speed', pw / (1 * 0.0647) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_adduction_l/speed', pw / (1 * 0.0410) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/hip_l/hip_rotation_l/speed', pw / (3 * 0.0728) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/knee_l/knee_angle_l/speed', pw / (1 * 0.0889) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/ankle_l/ankle_angle_l/speed', pw / (1 * 0.0574) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/subtalar_l/subtalar_angle_l/speed', pw / (3 * 0.0595) ** 2))
stateWeights.cloneAndAppend(osim.MocoWeight('/jointset/mtp_l/mtp_angle_l/speed', pw / (3 * 0.0873) ** 2))
track.set_states_weight_set(stateWeights)

# Define the Moco study and problem
study = track.initialize()
problem = study.updProblem()

# Define the periodicity goal
periodicityGoal = osim.MocoPeriodicityGoal('symmetryGoal')
problem.addGoal(periodicityGoal)
model = modelProcessor.process()
model.initSystem()

# All states are periodic except pelvis anterior-posterior translation
for i in range(model.getNumStateVariables()):
    currentStateName = str(model.getStateVariableNames().getitem(i))
    if 'pelvis_tx/value' not in currentStateName:
        periodicityGoal.addStatePair(osim.MocoPeriodicityGoalPair(currentStateName))

# All controls are periodic
for i in range(model.getNumControls()):
    currentControlName = str(problem.createRep().createControlInfoNames()[i])
    periodicityGoal.addControlPair(osim.MocoPeriodicityGoalPair(currentControlName))

# Prescribed average walking speed
speedGoal = osim.MocoAverageSpeedGoal('speed')
speedGoal.set_desired_average_speed(AvgSpeed)
problem.addGoal(speedGoal)

# Applying step width goal and preventing feet penetration
distanceConstraint = osim.MocoFrameDistanceConstraint()
if stepwidth == 1:
    distanceConstraint.setName("stepWidth")
    distanceConstraint.setProjection("vector")
    distanceConstraint.setProjectionVector(osim.Vec3(0, 0, 1))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/calcn_l', '/bodyset/calcn_r', stpwidth, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/toes_r', stpwidth, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/calcn_l', '/bodyset/toes_r', stpwidth, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/calcn_r', stpwidth, 100))
else:
    distanceConstraint.setName('distance_constraint')
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/calcn_l', '/bodyset/calcn_r', 0.10, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/toes_r', 0.10, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/calcn_l', '/bodyset/toes_r', 0.10, 100))
    distanceConstraint.addFramePair(osim.MocoFrameDistanceConstraintPair('/bodyset/toes_l', '/bodyset/calcn_r', 0.10, 100))

problem.addPathConstraint(distanceConstraint)

# Define upright torso goal
torsoGoal = osim.MocoOrientationTrackingGoal('torsoGoal', (3 / (0.0873 ** 2)) / 29)
torsoTable = osim.TableProcessor('modified_torso.sto')
torsoGoal.setStatesReference(torsoTable)
paths = osim.StdVectorString()
paths.append('/bodyset/torso')
torsoGoal.setFramePaths(paths)
problem.addGoal(torsoGoal)

# Control effort term (minimize squared muscle excitations)
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
effort.setWeight(w)
effort.setExponent(2)
effort.setDivideByDisplacement(False)

# Step length asymmetry goal
stepLengthAsymmetry = osim.MocoStepLengthAsymmetryGoal()
stepLengthAsymmetry.setWeight(asyweight)
stepLengthAsymmetry.setRightFootFrame('/bodyset/calcn_r')
stepLengthAsymmetry.setLeftFootFrame('/bodyset/calcn_l')
stepLengthAsymmetry.setAsymmetrySmoothing(4)
stepLengthAsymmetry.setTargetAsymmetry(targetasym)
stepLengthAsymmetry.setStrideLength(strideLen)
if asymmetrygoal == 1:
    problem.addGoal(stepLengthAsymmetry)

# GRF tracking
contactTracking = osim.MocoContactTrackingGoal('contact', 6 / (75.337 * 9.81) / (0.1073 ** 2) / 29)
contactTracking.setExternalLoadsFile('refGRF_3D.xml')

forceNamesRightFoot = osim.StdVectorString()
forceNamesRightFoot.append('/forceset/contactHeel_r')
forceNamesRightFoot.append('/forceset/contactMH1_r')
forceNamesRightFoot.append('/forceset/contactMH3_r')
forceNamesRightFoot.append('/forceset/contactMH5_r')
forceNamesRightFoot.append('/forceset/contactOtherToes_r')
forceNamesRightFoot.append('/forceset/contactHallux_r')
trackRightGRF = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'RightGRF')
trackRightGRF.append_alternative_frame_paths('/bodyset/toes_r')
contactTracking.addContactGroup(trackRightGRF)

forceNamesLeftFoot = osim.StdVectorString()
forceNamesLeftFoot.append('/forceset/contactHeel_l')
forceNamesLeftFoot.append('/forceset/contactMH1_l')
forceNamesLeftFoot.append('/forceset/contactMH3_l')
forceNamesLeftFoot.append('/forceset/contactMH5_l')
forceNamesLeftFoot.append('/forceset/contactOtherToes_l')
forceNamesLeftFoot.append('/forceset/contactHallux_l')
trackLeftGRF = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'LeftGRF')
trackLeftGRF.append_alternative_frame_paths('/bodyset/toes_l')
contactTracking.addContactGroup(trackLeftGRF)

problem.addGoal(contactTracking)

# Bound constraints
problem.setStateInfo('/jointset/groundPelvis/pelvis_tilt/value', [-10 * osim.SimTK_PI / 180, 10 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/groundPelvis/pelvis_list/value', [-10 * osim.SimTK_PI / 180, 10 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/groundPelvis/pelvis_rotation/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/groundPelvis/pelvis_tx/value', [-0.2, 2])
problem.setStateInfo('/jointset/groundPelvis/pelvis_ty/value', [0.75, 1.25])
problem.setStateInfo('/jointset/groundPelvis/pelvis_tz/value', [-0.25, 0.25])
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-25 * osim.SimTK_PI / 180, 45 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/hip_l/hip_adduction_l/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/hip_l/hip_rotation_l/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-25 * osim.SimTK_PI / 180, 45 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/hip_r/hip_adduction_r/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/hip_r/hip_rotation_r/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-75 * osim.SimTK_PI / 180, 0])
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-75 * osim.SimTK_PI / 180, 0])
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-30 * osim.SimTK_PI / 180, 20 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-30 * osim.SimTK_PI / 180, 20 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/subtalar_l/subtalar_angle_l/value', [-20 * osim.SimTK_PI / 180, 20 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/subtalar_r/subtalar_angle_r/value', [-20 * osim.SimTK_PI / 180, 20 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/mtp_l/mtp_angle_l/value', [-20 * osim.SimTK_PI / 180, 70 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/mtp_r/mtp_angle_r/value', [-20 * osim.SimTK_PI / 180, 70 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/back/lumbar_extension/value', [-20 * osim.SimTK_PI / 180, 20 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/back/lumbar_bending/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])
problem.setStateInfo('/jointset/back/lumbar_rotation/value', [-15 * osim.SimTK_PI / 180, 15 * osim.SimTK_PI / 180])

problem.setStateInfoPattern('/forceset/.*/activation', [0.001, 1.0], [], [])
problem.setStateInfoPattern('/forceset/lumbar_ext/activation', [-1.0, 1.0], [], [])
problem.setStateInfoPattern('/forceset/lumbar_bend/activation', [-1.0, 1.0], [], [])
problem.setStateInfoPattern('/forceset/lumbar_rot/activation', [-1.0, 1.0], [], [])

# Problem time bounds
#problem.setTimeBounds(0, study_final_time)

# Define the solver and set its options
solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
solver.set_optim_max_iterations(3000)
solver.set_num_mesh_intervals(50)
solver.set_optim_constraint_tolerance(1e-04)
solver.set_optim_convergence_tolerance(1e+01)
solver.set_minimize_implicit_auxiliary_derivatives(True)
solver.set_implicit_auxiliary_derivatives_weight(0.00001)
solver.set_enforce_constraint_derivatives(True)
solver.resetProblem(problem)

# Set the initial guess
solver.setGuessFile(initialguess + '.sto')

# Solve the problem


gaitTrackingSolution = study.solve()
# Write the solution to a file
gaitTrackingSolution.write(solutionname + '.sto')
# Write solution's GRF to a file
externalForcesTableFlat = osim.createExternalLoadsTableForGait(model, gaitTrackingSolution, forceNamesRightFoot, forceNamesLeftFoot)
osim.STOFileAdapter.write(externalForcesTableFlat, solutionname + '_GRF.sto')

# Breakdown of terms in objective function
print('   ')
print('Breakdown of objective (including weights):')
costTerms = []
for i in range(gaitTrackingSolution.getNumObjectiveTerms()):
   termName = gaitTrackingSolution.getObjectiveTermNames()[i]
   termNamestr = str(termName)
   costTerms.append((termNamestr, gaitTrackingSolution.getObjectiveTermByIndex(i)))
   print(costTerms)



# Reporting the asymmetry levels achieved
stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength = computeStepAsymmetryValues(
    solutionname + '.sto',
    solutionname + '_GRF.sto'
)
print(f'\nStep Time Asymmetry = {stepTimeAsymmetry:.3f}%')
print(f'Step Length Asymmetry = {stepLengthAsymmetry:.3f}%')
print(f'Right Step Length = {rightStepLength:.3f}%')
print(f'Left Step Length = {leftStepLength:.3f}%')

# Making the report file

data0 = pd.DataFrame({
    0: ["solutionname","track_ref","initial_guess","avg speed","target asym", "asym_weight", "Step Time Asymmetry", "stepLengthAsymmetry", "Right Step Length", "left Step Length", "num_iter", "solver duration", "track final time", "study final time"],
    1: [solutionname, reference_track, initialguess, AvgSpeed, targetasym, asyweight, stepTimeAsymmetry, stepLengthAsymmetry, rightStepLength, leftStepLength, gaitTrackingSolution.getNumIterations(), gaitTrackingSolution.getSolverDuration(), track_final_time, study_final_time]
})

# Muscle energy expenditure (Umberger et al., 2003; Koelewijn et al., 2018)
outputPaths = osim.StdVectorString()
outputPaths.append('.*fiber_velocity')
outputPaths.append('.*activation')
outputPaths.append('.*active_fiber_force')
outputPaths.append('.*fiber_length')
outputPaths.append('.*active_force_length_multiplier')
outputTable = study.analyze(gaitTrackingSolution, outputPaths)

# Save results to Excel
OFASym = pd.DataFrame(costTerms)
OFASym1 = pd.concat([data0, OFASym], axis=0)
OFASym1.to_excel(solutionname + ".xlsx", index=False)

# Visualize the solution
study.visualize(gaitTrackingSolution)


