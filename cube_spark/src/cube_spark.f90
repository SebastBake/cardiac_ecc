!> \file
!> $Id: cube_spark.f90 1528 2010-09-21 01:32:29Z  $
!> \author Chris Bradley
!> \brief This is test example of calcium spark kinetics and diffusion in cardiac cells.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example cube_spark.f90
!! Example program which sets up a single calcium release site inside a cubic cardiac cell geometry .
!! \par Latest Builds:
!<

!> Main program
PROGRAM CUBE_SPARK

  USE OPENCMISS
  USE MPI
  USE FIELDML_API

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  ! program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=6


  INTEGER(CMISSIntg), PARAMETER :: CaEquationsSetUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: CaMaterialsFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: CaFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: CaEquationsSetFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: iCaFieldUserNumber=11


  INTEGER(CMISSIntg), PARAMETER :: CaTnCFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: NumRyRFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=23

  INTEGER(CMISSIntg), PARAMETER :: FCaEquationsSetUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: FCaMaterialsFieldUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: FCaFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: FCaEquationsSetFieldUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: iFCaFieldUserNumber=28

  INTEGER(CMISSIntg), PARAMETER :: FEquationsSetUserNumber=29
  INTEGER(CMISSIntg), PARAMETER :: FMaterialsFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: FFieldUserNumber=31
  INTEGER(CMISSIntg), PARAMETER :: FEquationsSetFieldUserNumber=32
  INTEGER(CMISSIntg), PARAMETER :: iFFieldUserNumber=33


  !CMISS variables
  !setting up fields for Ca and CaM equation sets.
  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSFieldType) :: GeometricField,CaMaterialsField,FCaMaterialsField,FMaterialsField,CaField,FCaField,FField
  TYPE(CMISSFieldType) :: CaEquationsSetField,FCaEquationsSetField,FEquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSMeshElementsType) :: MeshElements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSEquationsType) :: CaEquations,FCaEquations,FEquations
  TYPE(CMISSEquationsSetType) :: CaEquationsSet,FCaEquationsSet,FEquationsSet
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldType) :: iCaField,CaTnCField,NumRyRField,iFCaField,iFField
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh    

  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES
  INTEGER :: st,i,NUMBER_OF_COORDS,NODES_PER_ELE
  INTEGER :: NUMBER_OF_NODES,node,NUMBER_OF_ELEMENTS,element,NUMBER_OF_RYRS,NODE_BD_LABEL
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:) :: RyRDensity
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  REAL(CMISSDP) :: nodex,nodey,nodez,sphere_constant,RELEASE_RADIUS
  LOGICAL :: EXPORT_FIELD=.FALSE.
  INTEGER(CMISSIntg) :: CELL_TYPE,NumRyRsPerCluster,RELEASE_NODE_RANK
  INTEGER(CMISSIntg) :: ryrModelIndex,GeometricMeshComponent
  INTEGER(CMISSIntg) :: Err,EquationsSetIndex,NODE_NUMBER,RYR_NODE_NUMBER,CONDITION,CellMLIndex,ELEM_NUMBER,NUMBER_RELEASE_NODES
  INTEGER(CMISSIntg),DIMENSION(166) :: CELLBOUNDARYNODES
  INTEGER(CMISSIntg) :: SL_BD_MARKER,MITO_BD_MARKER,MITO_REGION_MARKER,WITH_MITO_ELEMENTS,ELEM_LABEL
  REAL(CMISSDP) :: startT,endT,Tstep,ODE_TIME_STEP,VALUE,init_Ca, init_FCa, init_F, ryr_nodex,ryr_nodey,ryr_nodez, &
    & caDiffx, caDiffy,caDiffz,fcaDiffx, fcaDiffy,fcaDiffz,fDiffx,fDiffy,fDiffz,store_coeff,iCa,init_CaTnC,NodeRyRDensity
  INTEGER(CMISSIntg) :: NonZeroNodes
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,CELLPATH,RyRModel,RYRDENSITYFILE
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NodeDomain,ElementDomain
  INTEGER :: MPI_IERROR
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
    
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!_________________________________________________________________________________________________
  !Problem INPUTS. PARAMETERS FROM SOELLER et.al. 2009
  !MESH FILES
  open(unit=9,file='inputs.txt',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening inputs file',st
    STOP
  ELSE
    PRINT *,'inputs file opened correctly'
    READ(9,*) !file info
    READ(9,*) !comment on next set of variables - reading in file info
    READ(9,*) RYR_NODE_NUMBER
    READ(9,*)
    READ(9,*) NumRyRsPerCluster,NodeRyRDensity
    READ(9,*)
    READ(9,*) RyRModel
    READ(9,*)
    READ(9,*) init_Ca,caDiffx,caDiffy,caDiffz,iCa
    READ(9,*)
    READ(9,*) init_F,fDiffx,fDiffy,fDiffz
    READ(9,*)
    READ(9,*) init_FCa,fcaDiffx,fcaDiffy,fcaDiffz
    READ(9,*)
    READ(9,*) init_CaTnC
    READ(9,*) 
    READ(9,*) startT,endT,Tstep,ODE_TIME_STEP
  ENDIF
  CLOSE(9)
  !Write the params out to screen for double checking
  WRITE(*,*) 'RyR Cluster Node:', RYR_NODE_NUMBER
  WRITE(*,*) 'Number of RyRs per Cluster:', NumRyRsPerCluster
  WRITE(*,*) 'RyR Clusters per Unit Volume At Node:', NodeRyRDensity
  WRITE(*,*) 'CellML Model File:', RyRModel
  !cell initial conditions
  WRITE(*,*) 'Initial [Ca]i = ',init_Ca !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'Ca Diff Coeff in x = ',caDiffx
  WRITE(*,*) 'Ca Diff Coeff in y = ',caDiffy
  WRITE(*,*) 'Ca Diff Coeff in z = ',caDiffz

  WRITE(*,*) 'Initial [F]i = ',init_F !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'F Diff Coeff in x = ',fDiffx
  WRITE(*,*) 'F Diff Coeff in y = ',fDiffy
  WRITE(*,*) 'F Diff Coeff in z = ',fDiffz

  WRITE(*,*) 'Initial [FCa]i = ',init_FCa !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'FCa Diff Coeff in x = ',fcaDiffx
  WRITE(*,*) 'FCa Diff Coeff in y = ',fcaDiffy
  WRITE(*,*) 'FCa Diff Coeff in z = ',fcaDiffz

  WRITE(*,*) 'Initial Equil. [CaTnC] - uniform over myofibril region = ',init_CaTnC
  store_coeff = 1.0_CMISSDP

  WRITE(*,*) 'Tstart=',startT
  WRITE(*,*) 'Tend=',endT
  WRITE(*,*) 'Tstep=',Tstep
  WRITE(*,*) 'ODE_Tstep=',ODE_TIME_STEP
  
  EXPORT_FIELD=.FALSE.
!_________________________________________________________________________________________________
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  !get computational nodes for parallel processing
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

   !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemUserNumber,3,Err)
  !The coordinate system is 3D by default;set it to be 3D in above command.
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_LabelSet(Region,"Cube",Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

!_________________________________________________________________________________________________
  !Start the creation of a trilinear-simplex basis
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a trilinear simplex  basis
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  !set interpolation to be linear
  CALL CMISSBasis_InterpolationXiSet(Basis,(/CMISS_Basis_Linear_Simplex_Interpolation, &
   &   CMISS_Basis_Linear_Simplex_Interpolation, CMISS_Basis_Linear_Simplex_Interpolation/),Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)

!_________________________________________________________________________________________________  
  !Time to create a mesh - wohoo!
  !Read in nodes (set up RyRDensity array with column 
  !of zeros for later updating).
  NODEFILE="cuboid.1.node"
  ELEMFILE="cuboid.1.ele"
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  !Read in elements
  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,5))
    DO i = 1,NUMBER_OF_ELEMENTS
      READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
    ENDDO
  ENDIF 
  CLOSE(11)

  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL CMISSMeshElements_Initialise(MeshElements,Err)
  CALL CMISSMeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    CALL CMISSMeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
     &   ElemMap(i,4),ElemMap(i,5)/),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElements,Err)
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 3 field components in 1 mesh component
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      nodez = NodeCoords(i,3)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,3,nodez,Err)
     ENDIF
    ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Cube_Geom","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Cube_Geom","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF 

!______________________________________________________________________________________________________________
  WRITE(*,*) 'Create the cellml reaction with split reaction diffusion equations_set - 1 for each species' 
  !Ca equations
  CALL CMISSEquationsSet_Initialise(CaEquationsSet,Err)
  CALL CMISSField_Initialise(CaEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(CaEquationsSetUserNumber,Region, & 
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & CaEquationsSetFieldUserNumber,CaEquationsSetField,CaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(CaEquationsSet,Err)


  !Create the equations set dependent field variables for Ca
  CALL CMISSField_Initialise(CaField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(CaEquationsSet,CaFieldUserNumber,CaField,Err)
  CALL CMISSField_VariableLabelSet(CaField,CMISS_FIELD_U_VARIABLE_TYPE,"Ca Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(CaEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(CaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_Ca,Err)


  !Create the equations set material field variables - Ca
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMISSField_Initialise(CaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(CaEquationsSet,CaMaterialsFieldUserNumber,CaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"Ca Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(CaEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 1,caDiffx,Err) !ca diff coeff in x
  CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 2,caDiffy,Err) !ca diff coeff in y
  CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 3,caDiffz,Err) !ca diff coeff in z
  CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMISSField_ParameterSetUpdateStart(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaField
  !Might use the field for CellML input of elementary RyR calcium release
  CALL CMISSField_Initialise(iCaField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(CaEquationsSet,iCaFieldUserNumber,iCaField,Err)
  CALL CMISSField_VariableLabelSet(iCaField,CMISS_FIELD_U_VARIABLE_TYPE,"iCa Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(CaEquationsSet,Err)
  !Initialising the iCaField to iCa everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,iCa,Err)


  !NumRyRField
  !Set up number of RyRs in each cluster as a spatial distribution field
  !set up intensity field
  CALL CMISSField_Initialise(NumRyRField,Err)
  CALL CMISSField_CreateStart(NumRyRFieldUserNumber,Region,NumRyRField,Err)
  CALL CMISSField_TypeSet(NumRyRField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(NumRyRField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(NumRyRField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(NumRyRField,1,Err)
  CALL CMISSField_VariableTypesSet(NumRyRField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMISSField_DataTypeSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
  CALL CMISSField_DimensionSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,"Num RyR Field",Err)
  CALL CMISSField_ComponentMeshComponentGet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMISSField_ComponentMeshComponentSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMISSField_ComponentInterpolationSet(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMISSField_CreateFinish(NumRyRField,Err)
  !Initialise Num RyR Field
  CALL CMISSField_ComponentValuesInitialise(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !Define a single RyR release node, multiply the intensity by number of ryrs per cluster
  CALL CMISSDecomposition_NodeDomainGet(Decomposition,RYR_NODE_NUMBER,1,NodeDomain,Err)
  RELEASE_NODE_RANK = NodeDomain

  IF(NodeDomain.EQ.ComputationalNodeNumber) THEN
    RELEASE_NODE_RANK = ComputationalNodeNumber
    CALL CMISSField_ParameterSetUpdateNode(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,1,1,RYR_NODE_NUMBER,1,(NumRyRsPerCluster*NodeRyRDensity),Err)

    !noting the coordinates of the ryr release node.
    CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,1,ryr_nodex,Err)
    CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,2,ryr_nodey,Err)
    CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,3,ryr_nodez,Err)

  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERROR)

  CALL MPI_BCAST(ryr_nodex,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(ryr_nodey,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(ryr_nodez,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)

  CALL CMISSField_ParameterSetUpdateStart(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Now assign ryr density to nodes that fall within a 100 nm (0.1 micron) radius of the center of the RyR_Node_Number
  RELEASE_RADIUS = 0.1_CMISSDP
  NUMBER_RELEASE_NODES=1
  DO node=1,NUMBER_OF_NODES
    NODE_NUMBER=NodeNums(node,1)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain.EQ.ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,1,nodex,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,2,nodey,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,3,nodez,Err)
      sphere_constant = SQRT((ryr_nodex-nodex)**2+(ryr_nodey-nodey)**2+(ryr_nodez-nodez)**2)
      IF(sphere_constant.LE.RELEASE_RADIUS) THEN
        WRITE(*,*) 'Node Number',NODE_NUMBER
        WRITE(*,*) 'COORDS:', nodex,nodey,nodez
        WRITE(*,*) 'REL CENTRE:',ryr_nodex,ryr_nodey,ryr_nodez
        WRITE(*,*) 'DISTANCE:',sphere_constant 
        CALL CMISSField_ParameterSetUpdateNode(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(NumRyRsPerCluster*NodeRyRDensity),Err)
        NUMBER_RELEASE_NODES = NUMBER_RELEASE_NODES+1
      ENDIF
    ENDIF
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  WRITE(*,*) "The final number of nodes which will be releasing Ca2+ =", NUMBER_RELEASE_NODES
  !Set up the fields for the other buffers which will store concentrations of the Ca-Buffer complex
  !F equations
  CALL CMISSEquationsSet_Initialise(FEquationsSet,Err)
  CALL CMISSField_Initialise(FEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(FEquationsSetUserNumber,Region, & 
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & FEquationsSetFieldUserNumber,FEquationsSetField,FEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(FEquationsSet,Err)


  !Create the equations set dependent field variables for F
  CALL CMISSField_Initialise(FField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(FEquationsSet,FFieldUserNumber,FField,Err)
  CALL CMISSField_VariableLabelSet(FField,CMISS_FIELD_U_VARIABLE_TYPE,"F Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(FEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(FField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_F,Err)


  !Create the equations set material field variables - F
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMISSField_Initialise(FMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(FEquationsSet,FMaterialsFieldUserNumber,FMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"F Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(FEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 1,fDiffx,Err) !f diff coeff in x
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 2,fDiffy,Err) !f diff coeff in y
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 3,fDiffz,Err) !f diff coeff in z
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMISSField_ParameterSetUpdateStart(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFField
  CALL CMISSField_Initialise(iFField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(FEquationsSet,iFFieldUserNumber,iFField,Err)
  CALL CMISSField_VariableLabelSet(iFField,CMISS_FIELD_U_VARIABLE_TYPE,"iF Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(FEquationsSet,Err)
  !Initialising the iField to zero everywhere. Might modify for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iFField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !FCa equations
  CALL CMISSEquationsSet_Initialise(FCaEquationsSet,Err)
  CALL CMISSField_Initialise(FCaEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(FCaEquationsSetUserNumber,Region, & 
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & FCaEquationsSetFieldUserNumber,FCaEquationsSetField,FCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(FCaEquationsSet,Err)


  !Create the equations set dependent field variables for FCa
  CALL CMISSField_Initialise(FCaField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(FCaEquationsSet,FCaFieldUserNumber,FCaField,Err)
  CALL CMISSField_VariableLabelSet(FCaField,CMISS_FIELD_U_VARIABLE_TYPE,"FCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(FCaEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(FCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_FCa,Err)


  !Create the equations set material field variables - FCa
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMISSField_Initialise(FCaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(FCaEquationsSet,FCaMaterialsFieldUserNumber,FCaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"FCa Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(FCaEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 1,fcaDiffx,Err) !fca diff coeff in x
  CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 2,fcaDiffy,Err) !fca diff coeff in y
  CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 3,fcaDiffz,Err) !fca diff coeff in z
  CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMISSField_ParameterSetUpdateStart(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFCaField
  CALL CMISSField_Initialise(iFCaField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(FCaEquationsSet,iFCaFieldUserNumber,iFCaField,Err)
  CALL CMISSField_VariableLabelSet(iFCaField,CMISS_FIELD_U_VARIABLE_TYPE,"iFCa Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(FCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iFCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !CaTnC
  CALL CMISSField_Initialise(CaTnCField,Err)
  CALL CMISSField_CreateStart(CaTnCFieldUserNumber,Region,CaTnCField,Err)
  CALL CMISSField_TypeSet(CaTnCField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(CaTnCField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(CaTnCField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(CaTnCField,1,Err)
  CALL CMISSField_VariableTypesSet(CaTnCField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMISSField_DataTypeSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
  CALL CMISSField_DimensionSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,"CaTnC Field",Err)
  CALL CMISSField_ComponentMeshComponentGet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMISSField_ComponentMeshComponentSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMISSField_ComponentInterpolationSet(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMISSField_CreateFinish(CaTnCField,Err)
  !Initialise CaTnC concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMISSField_ComponentValuesInitialise(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,init_CaTnC,Err)
  CALL CMISSField_ParameterSetUpdateStart(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  
!_____________________________________________________________________________________________________________________

  WRITE(*,*) 'Start to set up CellML Fields'

  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import ryr release and buffer source model from a file
  CALL CMISSCellML_ModelImport(CellML,RyRModel,ryrModelIndex,Err)
  ! set iCa as known so that it can be set as spatially varying in opencmiss.
  !CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  ! set RyRDensity as known so that it can be set as spatially varying in opencmiss.
  CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"Dyad/iCa",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"Dyad/NumRyR",Err)

  !to get from the CellML side. variables in cellml model that are not state variables, but are dependent on independent and state variables. 
  !- components of intermediate field
  !fluxes of the different buffers and CaRUs that I want to get out as intermediate variables
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"Dyad/J_ryr",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"Cytoplasm/J_f4",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"Cytoplasm/J_tnc",Err)
  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  !Mapping free calcium in opencmiss to that in cellml.
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/Ca_i",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/Ca_i",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping iCaField to iCa in the cellml model
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,iCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/iCa",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/iCa",CMISS_FIELD_VALUES_SET_TYPE, &
    & iCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping NumRyRField to RyRDensity in the cellml model
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/NumRyR",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/NumRyR",CMISS_FIELD_VALUES_SET_TYPE, &
    & NumRyRField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping Buffer-Complex resting values of cellml model to appropriate fields set up above

   !Mapping F
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,FField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/F4",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/F4",CMISS_FIELD_VALUES_SET_TYPE, &
    & FField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping FCa
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,FCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/F4Ca",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/F4Ca",CMISS_FIELD_VALUES_SET_TYPE, &
    & FCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping CaTnC
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/CaiTnC",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/CaiTnC",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellML_FieldMapsCreateFinish(CellML,Err)


  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL CMISSField_Initialise(CellMLModelsField,Err)
  CALL CMISSCellML_ModelsFieldCreateStart(CellML, &
    & CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellML_ModelsFieldCreateFinish(CellML,Err)

  !By default all field parameters have default model value of 1, i.e. the first model. 
  ! assigning the bufferNryr cellml model (model 1) for all nodes.
  CALL CMISSField_ComponentValuesInitialise(CellMLModelsField, & 
    & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)
  !Start the creation of the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML, &
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)


  !Create the equations set equations for Ca
  WRITE(*,*) 'Creating the equations for the equation sets'
  CALL CMISSEquations_Initialise(CaEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(CaEquationsSet,CaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(CaEquations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(CaEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(CaEquationsSet,Err)

  !Create the equations set equations for F
  CALL CMISSEquations_Initialise(FEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(FEquationsSet,FEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(FEquations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(FEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(FEquationsSet,Err)

  !Create the equations set equations for FCa
  CALL CMISSEquations_Initialise(FCaEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(FCaEquationsSet,FCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(FCaEquations,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(FCaEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(FCaEquationsSet,Err)
!____________________________________________________________________________________________________________

  WRITE(*,*) 'Create the problem'
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a strang split reaction diffusion problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMISS_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,startT,endT,Tstep,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,1,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,5.00_CMISSDP,0.01_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

!______________________________________________________________________________________________________________
  WRITE(*,*) 'Set up the problem solvers for Strang splitting'
  !note, this example is contrived to have strang splitting, when it could be solved as a simple evaluation (as opposed to integration) of source and diffusion
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL CMISSSolver_DynamicThetaSet(Solver,1.0_CMISSDP,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL CMISSSolver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  !set linear solver to be direct solver. Note, I found this stuff in fluidmechanics/darcy/dynamic/src example
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverCMISSLibrary,Err)
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverMUMPSLibrary,Err)
  CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)


  !Third solver is another DAE solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL CMISSSolver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err) !set the third solver's integration time step
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_NO_OUTPUT,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)


  !Start the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSSolver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL CMISSCellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMISSSolver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set for Ca, F and FCa
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,CaEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,FEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,FCaEquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

!_________________________________________________________________________________________________________
  WRITE(*,*) 'Set up boundary conditions'
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  DO node=1,NUMBER_OF_NODES
    NODE_BD_LABEL = NodeNums(node,2)
    NODE_NUMBER = NodeNums(node,1)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NODE_BD_LABEL.EQ.10) THEN !set no flux bc if node is a boundary node.
        CONDITION = CMISS_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSDP
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaField, &
         & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FField, &
         & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FCaField, &
         & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

      ENDIF
    ENDIF
  ENDDO
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
!__________________________________________________________________________________________________________
  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)

  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Ca_Cube_Solution","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Ca_Cube_Solution","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF 
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP


END PROGRAM CUBE_SPARK

