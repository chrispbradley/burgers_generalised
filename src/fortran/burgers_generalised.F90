PROGRAM BurgersGeneralised

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifndef NOMPIMOD
  USE MPI
#endif
  IMPLICIT NONE
#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=14

  !Program variables
  INTEGER(CMISSIntg) :: dynamicSolverOutputType
  INTEGER(CMISSIntg) :: nonlinearSolverOutputType
  INTEGER(CMISSIntg) :: linearSolverOutputType
  INTEGER(CMISSIntg) :: equationsOutputType
  INTEGER(CMISSIntg) :: numberGlobalXElements
  REAL(CMISSRP) :: length
  REAL(CMISSRP) :: startTime
  REAL(CMISSRP) :: stopTime
  REAL(CMISSRP) :: timeIncrement

  !Program types
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) ::  analyticField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: dynamicSolver,nonlinearSolver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex,err
  LOGICAL :: directLinearSolverFlag,exportField

  !Intialise OpenCMISS
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Initialise(context,err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  !Set number of elements for FEM discretization
  numberGlobalXElements=6
  length=3.0_CMISSRP

  !Set output parameters
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  !dynamicSolverOutputType=CMFE_SOLVER_NO_OUTPUT
  dynamicSolverOutputType=CMFE_SOLVER_MATRIX_OUTPUT
  !nonlinearSolverOutputType=CMFE_SOLVER_NO_OUTPUT
  nonlinearSolverOutputType=CMFE_SOLVER_MATRIX_OUTPUT
  !linearSolverOutputType=CMFE_SOLVER_NO_OUTPUT
  linearSolverOutputType=CMFE_SOLVER_MATRIX_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  !equationsOutputType=CMFE_EQUATIONS_NO_OUTPUT
  equationsOutputType=CMFE_EQUATIONS_MATRIX_OUTPUT

  !Set time parameter
  startTime=0.0_CMISSRP
  stopTime=0.100001_CMISSRP
  timeIncrement=0.01_CMISSRP

  !Solver parameters
  directLinearSolverFlag=.FALSE.

  !Export parameters
  exportField=.FALSE.

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"BurgersRegion",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL cmfe_Basis_NumberOfXiSet(basis,1,err)
  !Set the basis xi interpolation and number of Gauss points
  CALL cmfe_Basis_InterpolationXiSet(basis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[3],err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[length],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberGlobalXElements],err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region  
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_NO_SCALING,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_set for a Generalised Burgers's equation
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMFE_EQUATIONS_SET_BURGERS_EQUATION_TYPE,CMFE_EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set material field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  !Initialise materials field
  !Set A
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,1.0_CMISSRP,err)
  !Set B
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,-1.0_CMISSRP,err)
  !Set C
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,1.0_CMISSRP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(analyticField,err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSet,CMFE_EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2, &
    & ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type (Sparse/Full)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output (NoOutput/TimingOutput/MatrixOutput/SolverMatrix/ElementMatrixOutput)
  CALL cmfe_Equations_OutputTypeSet(equations,equationsOutputType,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS, &
    & CMFE_PROBLEM_BURGERS_EQUATION_TYPE,CMFE_PROBLEM_DYNAMIC_BURGERS_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(controlLoop,startTime,stopTime,timeIncrement,err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(controlLoop,1,err)
  CALL cmfe_ControlLoop_OutputTypeSet(controlLoop,CMFE_CONTROL_LOOP_NO_OUTPUT,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(dynamicSolver,err)
  CALL cmfe_Solver_Initialise(nonlinearSolver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)

  !Get the dymamic solver
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,dynamicSolver,err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(dynamicSolver,dynamicSolverOutputType,err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(dynamicSolver,0.5_CMISSRP,err)

  !Get the dynamic nonlinear solver
  CALL cmfe_Solver_DynamicNonlinearSolverGet(dynamicSolver,nonlinearSolver,err)
  !Set the nonlinear Jacobian type
  !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  !Set the line search
  CALL cmfe_Solver_NewtonLineSearchTypeSet(nonlinearSolver,CMFE_SOLVER_NEWTON_LINESEARCH_LINEAR,err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(nonlinearSolver,nonlinearSolverOutputType,err)
  !Get the dynamic nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(nonlinearSolver,linearSolver,err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(linearSolver,linearSolverOutputType,err)
  !Set the solver settings

  IF(directLinearSolverFlag) THEN
    CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL cmfe_Solver_LibraryTypeSet(linearSolver,CMFE_SOLVER_MUMPS_LIBRARY,err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(linearSolver,50,err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)


  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the dynamic solver equations
  CALL cmfe_Solver_Initialise(dynamicSolver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,dynamicSolver,err)
  CALL cmfe_Solver_SolverEquationsGet(dynamicSolver,solverEquations,err)
  !Set the solver equations sparsity (Sparse/Full)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,EquationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(problem,err)
  WRITE(*,'(A)') "Problem solved!"

  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(dependentField,"GeneralisedBurgersAnalytics",err)

  !export fields
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"BurgersGeneralised","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"BurgersGeneralised","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF

  CALL cmfe_Finalise(context,err)
  WRITE(*,'(A)') "Program successfully completed."

END PROGRAM BurgersGeneralised
