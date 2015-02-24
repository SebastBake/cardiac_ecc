!> \file
!> $Id: cardiac_ecc.f90 2014-12-16 vraj004 $
!> \author Vijay Rajagopal
!> \brief Main program file to simulate cardiac ECC using opencmiss library routines
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

!> \example cardiac_ecc.f90
!! Example program which sets up a singpe species diffusion problem.
!! \par Latest Builds:
!<

!> Main program
PROGRAM CARDIAC_ECC

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
  INTEGER(CMISSIntg), PARAMETER :: RyRDenseFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: RyRReleaseLagFieldUserNumber=34


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

  INTEGER(CMISSIntg), PARAMETER :: CaMEquationsSetUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: CaMMaterialsFieldUserNumber=35
  INTEGER(CMISSIntg), PARAMETER :: CaMFieldUserNumber=36
  INTEGER(CMISSIntg), PARAMETER :: CaMEquationsSetFieldUserNumber=37
  INTEGER(CMISSIntg), PARAMETER :: iCaMFieldUserNumber=38

  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetUserNumber=39
  INTEGER(CMISSIntg), PARAMETER :: ATPMaterialsFieldUserNumber=40
  INTEGER(CMISSIntg), PARAMETER :: ATPFieldUserNumber=41
  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetFieldUserNumber=42
  INTEGER(CMISSIntg), PARAMETER :: iATPFieldUserNumber=43

  INTEGER(CMISSIntg), PARAMETER :: CaMCaEquationsSetUserNumber=44
  INTEGER(CMISSIntg), PARAMETER :: CaMCaMaterialsFieldUserNumber=45
  INTEGER(CMISSIntg), PARAMETER :: CaMCaFieldUserNumber=46
  INTEGER(CMISSIntg), PARAMETER :: CaMCaEquationsSetFieldUserNumber=47
  INTEGER(CMISSIntg), PARAMETER :: iCaMCaFieldUserNumber=48


  INTEGER(CMISSIntg), PARAMETER :: ATPCaEquationsSetUserNumber=49
  INTEGER(CMISSIntg), PARAMETER :: ATPCaMaterialsFieldUserNumber=50
  INTEGER(CMISSIntg), PARAMETER :: ATPCaFieldUserNumber=51
  INTEGER(CMISSIntg), PARAMETER :: ATPCaEquationsSetFieldUserNumber=52
  INTEGER(CMISSIntg), PARAMETER :: iATPCaFieldUserNumber=53

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22


  !Defining CMISS-type variables
  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSFieldType) :: GeometricField,CaMaterialsField,FCaMaterialsField,FMaterialsField,CaField,FCaField,FField
  TYPE(CMISSFieldType) :: CaEquationsSetField,FCaEquationsSetField,FEquationsSetField
  TYPE(CMISSFieldType) :: CaMMaterialsField, CaMEquationsSetField,CaMField
  TYPE(CMISSFieldType) :: CaMCaMaterialsField, CaMCaEquationsSetField,CaMCaField
  TYPE(CMISSFieldType) :: ATPField,ATPMaterialsField,ATPEquationsSetField
  TYPE(CMISSFieldType) :: ATPCaField,ATPCaMaterialsField,ATPCaEquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSMeshElementsType) :: MeshElements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSEquationsType) :: CaEquations,FCaEquations,FEquations,CaMEquations,CaMCaEquations,ATPEquations
  TYPE(CMISSEquationsType) :: ATPCaEquations
  TYPE(CMISSEquationsSetType) :: CaEquationsSet,FCaEquationsSet,FEquationsSet,ATPEquationsSet,CaMEquationsSet
  TYPE(CMISSEquationsSetType) :: ATPCaEquationsSet,CaMCaEquationsSet
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldType) :: iCaField,CaTnCField,RyRDenseField,iFCaField,iFField,RyRReleaseLagField
  TYPE(CMISSFieldType) :: iCaMField,iCaMCaField,iATPField,iATPCaField

  !Defining program-specific fortran variables

  INTEGER(CMISSIntg) :: MPI_IERROR,NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES,outputfreq
  INTEGER :: st,i,NUMBER_OF_COORDS,NODES_PER_ELE,NUMBER_OF_MITOBDFACENODES,FACE_ATTRIBUTES,NUMBER_OF_CELLBDNODES
  INTEGER :: NUMBER_OF_NODES,node,NUMBER_OF_ELEMENTS,element,NUMBER_OF_RYRS,NODE_BD_LABEL
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:,:) :: RyRDensity
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap,MITOBDFaceMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:) :: MITOBDFaceNodes,CELLBDNodes
  REAL(CMISSDP) :: nodex,nodey,nodez
  LOGICAL :: EXPORT_FIELD=.FALSE.
  INTEGER(CMISSIntg) :: CELL_TYPE
  INTEGER(CMISSIntg) :: ryrModelIndex,GeometricMeshComponent
  INTEGER(CMISSIntg) :: Err,EquationsSetIndex,NODE_NUMBER,CONDITION,CellMLIndex,ELEM_NUMBER
  INTEGER(CMISSIntg) :: SL_BD_MARKER,MITO_BD_MARKER,MITO_REGION_MARKER, &
    & WITH_MITO_ELEMENTS,ELEM_LABEL,CYTO_REGION_MARKER,MODEL_ON
  REAL(CMISSDP) :: startT,endT,Tstep,ODE_TIME_STEP,VALUE,init_Ca,init_FCa, init_F, ryr_nodex,ryr_nodey,ryr_nodez, &
    & caDiffx, caDiffy,caDiffz,fcaDiffx, fcaDiffy,fcaDiffz,fDiffx,fDiffy,fDiffz,store_coeff,iCa,init_CaTnC, &
    & NodeRyRDensity,mitoCaDiffx,mitoCaDiffy,mitoCaDiffz,mito_initCa,mito_initF, mito_initFCa,mito_initCaTnC,&
    & mitoFDiffx,mitoFDiffy,mitoFDiffz,mitoFCaDiffx,mitoFCaDiffy,mitoFCaDiffz,init_CaM,init_ATP,init_CaMCa,init_ATPCa, &
    & mitoCaMDiffx,mitoCaMDiffy,mitoCaMDiffz,mitoCaMCaDiffx,mitoCaMCaDiffy,mitoCaMCaDiffz,mitoATPDiffx,mitoATPDiffy,mitoATPDiffz, &
    & mitoATPCaDiffx,mitoATPCaDiffy,mitoATPCaDiffz,CaMDiffx, CaMDiffy,CaMDiffz,CaMCaDiffx, CaMCaDiffy,CaMCaDiffz, &
    & ATPDiffx, ATPDiffy,ATPDiffz,ATPCaDiffx, ATPCaDiffy,ATPCaDiffz,mito_initCaM, mito_initCaMCa,mito_initATP, mito_initATPCa
  INTEGER(CMISSIntg) :: NumRyRsPerCluster,NonZeroNodes,MITOBDFaceNode_idx,NUMBER_OF_MITOBDFACES
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,CELLPATH,RyRModel,RYRDENSITYFILE,MITOBDFACEFILE,CELLBDNODESFILE
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NodeDomain,ElementDomain




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
    READ(9,*) CELLID
    READ(9,*)
    READ(9,*) WITH_MITO_ELEMENTS
    READ(9,*)
    READ(9,*) NODEFILE,ELEMFILE
    READ(9,*)
    READ(9,*) MITOBDFACEFILE,CELLBDNODESFILE
    READ(9,*) 
    READ(9,*) RYRDENSITYFILE,NumRyRsPerCluster
    READ(9,*)
    READ(9,*) RyRModel,MODEL_ON
    READ(9,*)
    READ(9,*) init_Ca,caDiffx,caDiffy,caDiffz,iCa
    READ(9,*)
    READ(9,*) init_F,fDiffx,fDiffy,fDiffz
    READ(9,*)
    READ(9,*) init_FCa,fcaDiffx,fcaDiffy,fcaDiffz
    READ(9,*)
    READ(9,*) init_CaM,CaMDiffx,CaMDiffy,CaMDiffz
    READ(9,*)
    READ(9,*) init_CaMCa,CaMCaDiffx,CaMCaDiffy,CaMCaDiffz
    READ(9,*)
    READ(9,*) init_ATP,ATPDiffx,ATPDiffy,ATPDiffz
    READ(9,*)
    READ(9,*) init_ATPCa,ATPCaDiffx,ATPCaDiffy,ATPCaDiffz
    READ(9,*)
    READ(9,*) init_CaTnC
    READ(9,*)
    READ(9,*) mito_initCa,mitoCaDiffx,mitoCaDiffy,mitoCaDiffz
    READ(9,*)
    READ(9,*) mito_initF,mitoFDiffx,mitoFDiffy,mitoFDiffz
    READ(9,*)
    READ(9,*) mito_initFCa,mitoFCaDiffx,mitoFCaDiffy,mitoFCaDiffz
    READ(9,*)
    READ(9,*) mito_initCaM,mitoCaMDiffx,mitoCaMDiffy,mitoCaMDiffz
    READ(9,*)
    READ(9,*) mito_initCaMCa,mitoCaMCaDiffx,mitoCaMCaDiffy,mitoCaMCaDiffz
    READ(9,*)
    READ(9,*) mito_initATP,mitoATPDiffx,mitoATPDiffy,mitoATPDiffz
    READ(9,*)
    READ(9,*) mito_initATPCa,mitoATPCaDiffx,mitoATPCaDiffy,mitoATPCaDiffz
    READ(9,*)
    READ(9,*) mito_initCaTnC
    READ(9,*)
    READ(9,*) startT,endT,Tstep,ODE_TIME_STEP
    READ(9,*) 
    READ(9,*) outputfreq
  ENDIF
  CLOSE(9)
  !Write the params out to screen for double checking

  WRITE(*,*) 'Cell ID is:',CELLID
  WRITE(*,*) 'MITO REPRESENTATION, WITH_MITO_ELEMENTS=', WITH_MITO_ELEMENTS
  WRITE(*,*) 'Node file:',NODEFILE
  WRITE(*,*) 'Element file:',ELEMFILE
  WRITE(*,*) 'Boundary Face file:',MITOBDFACEFILE
  WRITE(*,*) 'RyR Intensity File:',RYRDENSITYFILE
  WRITE(*,*) 'Number of RyRs per Cluster:', NumRyRsPerCluster
  WRITE(*,*) 'CellML Model File:', RyRModel
  WRITE(*,*) 'CellML Model is switched on:', MODEL_ON
  !cell initial conditions
  WRITE(*,*) 'Initial [Ca]i = ',init_Ca !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'Ca Diff Coeff in x = ',caDiffx
  WRITE(*,*) 'Ca Diff Coeff in y = ',caDiffy
  WRITE(*,*) 'Ca Diff Coeff in z = ',caDiffz
  WRITE(*,*) 'RyR Ca release current = ',iCa
  WRITE(*,*) 'Initial [F]i = ',init_F 
  WRITE(*,*) 'F Diff Coeff in x = ',fDiffx
  WRITE(*,*) 'F Diff Coeff in y = ',fDiffy
  WRITE(*,*) 'F Diff Coeff in z = ',fDiffz

  WRITE(*,*) 'Initial [FCa]i = ',init_FCa 
  WRITE(*,*) 'FCa Diff Coeff in x = ',fcaDiffx
  WRITE(*,*) 'FCa Diff Coeff in y = ',fcaDiffy
  WRITE(*,*) 'FCa Diff Coeff in z = ',fcaDiffz

  WRITE(*,*) 'Initial Equil. [CaTnC] - uniform over myofibril region = ',init_CaTnC
  store_coeff = 1.0_CMISSDP
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    WRITE(*,*) 'Mito Diff Coeff in x = ',mitoCaDiffx
    WRITE(*,*) 'Mito Diff Coeff in y = ',mitoCaDiffy
    WRITE(*,*) 'Mito Diff Coeff in z = ',mitoCaDiffz
    WRITE(*,*) 'Mito Diff Coeff in x = ',mitoFDiffx
    WRITE(*,*) 'Mito Diff Coeff in y = ',mitoFDiffy
    WRITE(*,*) 'Mito Diff Coeff in z = ',mitoFDiffz
    WRITE(*,*) 'Mito Diff Coeff in x = ',mitoFCaDiffx
    WRITE(*,*) 'Mito Diff Coeff in y = ',mitoFCaDiffy
    WRITE(*,*) 'Mito Diff Coeff in z = ',mitoFCaDiffz


    WRITE(*,*) 'Initial [Ca]mito = ',mito_initCa
    WRITE(*,*) 'Initial [F]mito = ',mito_initF
    WRITE(*,*) 'Initial [FCa]mito = ',mito_initFCa
    WRITE(*,*) 'Initial [CaTnC]mito = ',mito_initCaTnC
  ENDIF
  store_coeff = 1.0_CMISSDP

  WRITE(*,*) 'Tstart=',startT
  WRITE(*,*) 'Tend=',endT
  WRITE(*,*) 'Tstep=',Tstep
  WRITE(*,*) 'ODE_Tstep=',ODE_TIME_STEP
  WRITE(*,*) 'Output Frequency=',outputfreq


  !Boundary Markers
  SL_BD_MARKER=0
  MITO_BD_MARKER=20
  MITO_REGION_MARKER = 1
  CYTO_REGION_MARKER = 2
  
  EXPORT_FIELD=.FALSE.
!_________________________________________________________________________________________________
  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  !get computational nodes for parallel processing
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)


  !set diagnostics
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"SOLVE_DIAGNOSTICS", &
  !  & (/"REACTION_DIFFUSION_PRE_SOLVE"/),Err)
  CALL MPI_BCAST(NumberOfComputationalNodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)


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
  CALL CMISSRegion_LabelSet(Region,"Cell",Err)
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

  !Time to create a mesh - wohoo!
  !Read in nodes (set up RyRDensity array with column 
  !of zeros for later updating).
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(RyRDensity(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  PRINT *, 'Total Nodes',NUMBER_OF_NODES
  !Read in elements
  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    IF(WITH_MITO_ELEMENTS.EQ.1 .AND. ELE_ATTRIBUTES.EQ.1) THEN
      ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,6))
      DO i = 1,NUMBER_OF_ELEMENTS
        READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5),ElemMap(i,6)
        IF(ElemMap(i,6).EQ.MITO_REGION_MARKER) THEN
          DO node = 2,5
            NODE_NUMBER = ElemMap(i,node)
            NodeNums(NODE_NUMBER,2) = MITO_REGION_MARKER
          ENDDO
        ENDIF
      ENDDO

    ELSE
      ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,5))
      DO i = 1,NUMBER_OF_ELEMENTS
        READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
      ENDDO

    ENDIF
  ENDIF 
  CLOSE(11)
  PRINT *,'Total Elements',NUMBER_OF_ELEMENTS
  CALL MPI_BCAST(NUMBER_OF_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)


  !Read in MITO BDFACES
  OPEN(unit=12,file=MITOBDFACEFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening face file',st
    STOP
  ELSE
    PRINT *,'Face file opened successfully'
    READ(12,*) NUMBER_OF_MITOBDFACES,FACE_ATTRIBUTES
    ALLOCATE(MITOBDFaceMap(NUMBER_OF_MITOBDFACES,5))
    NUMBER_OF_MITOBDFACENODES = NUMBER_OF_MITOBDFACES*3
    ALLOCATE(MITOBDFaceNodes(NUMBER_OF_MITOBDFACENODES))
    MITOBDFaceNode_idx = 0
    DO i = 1,NUMBER_OF_MITOBDFACES
      READ(12,*) MITOBDFaceMap(i,1),MITOBDFaceMap(i,2),MITOBDFaceMap(i,3),MITOBDFaceMap(i,4),MITOBDFaceMap(i,5)
      DO node = 2,4
        MITOBDFaceNodes(MITOBDFaceNode_idx+1) = MITOBDFaceMap(i,node)
        MITOBDFaceNode_idx = MITOBDFaceNode_idx+1
      ENDDO
    ENDDO
  ENDIF 
  CLOSE(12)
  PRINT *,'Total Mitochondrial Boundary Faces',NUMBER_OF_MITOBDFACES


  !Read in Cell BDNodes
  OPEN(unit=13,file=CELLBDNODESFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening face file',st
    STOP
  ELSE
    PRINT *,'Cell BD Nodes file opened successfully'
    READ(13,*) NUMBER_OF_CELLBDNODES
    ALLOCATE(CELLBDNodes(NUMBER_OF_CELLBDNODES))
    DO i = 1,NUMBER_OF_CELLBDNODES
      READ(13,*) CELLBDNodes(i)
    ENDDO
  ENDIF
  CLOSE(13)
  PRINT *,'Total Cell Boundary Nodes',NUMBER_OF_CELLBDNODES

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
    CALL CMISSFields_NodesExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF 
!______________________________________________________________________________________________________________
  !Create the cellml reaction with split reaction diffusion equations_set - 1 for each species
!###################
  !Ca equations
!###################
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
  !Initialise Ca dependent field
  CALL CMISSField_ComponentValuesInitialise(CaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_Ca,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL CMISSField_ParameterSetUpdateNode(CaField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF

  !Create the equations set material field variables - Ca
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMISSField_Initialise(CaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(CaEquationsSet,CaMaterialsFieldUserNumber,CaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"Ca Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(CaEquationsSet,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaDiffz,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,caDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,caDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,caDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 1,caDiffx,Err) !ca diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 2,caDiffy,Err) !ca diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 3,caDiffz,Err) !ca diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(CaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaField
  !Might use the field for CellML input of elementary RyR calcium release
  CALL CMISSField_Initialise(iCaField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(CaEquationsSet,iCaFieldUserNumber,iCaField,Err)
  CALL CMISSField_VariableLabelSet(iCaField,CMISS_FIELD_U_VARIABLE_TYPE,"iCa Field",Err)
  !CALL CMISSField_ComponentInterpolationSet(iCaField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
  !  & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(CaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,iCa,Err)

!###################
  !RyRDenseField and RyRReleaseLagField
!###################
  !Set up RyR spatial density field and the RyRReleaseLagFields
  !Read in file containing estimates of spatial density and time lag of release at each node
  !update RyRDensity array with density values
  OPEN(unit=12,file=RYRDENSITYFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening RyR density file',st
    STOP
  ELSE
    PRINT *,'RyR density file opened successfully'
    DO i = 1,NUMBER_OF_NODES
      READ(12,*) RyRDensity(i,1),RyRDensity(i,2)
    ENDDO
  ENDIF 
  CLOSE(12)
  !set up intensity field
  CALL CMISSField_Initialise(RyRDenseField,Err)
  CALL CMISSField_CreateStart(RyRDenseFieldUserNumber,Region,RyRDenseField,Err)
  CALL CMISSField_TypeSet(RyRDenseField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(RyRDenseField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(RyRDenseField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(RyRDenseField,1,Err)
  CALL CMISSField_VariableTypesSet(RyRDenseField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMISSField_DataTypeSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
  CALL CMISSField_DimensionSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,"RyR Density Field",Err)
  CALL CMISSField_ComponentMeshComponentGet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMISSField_ComponentMeshComponentSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMISSField_ComponentInterpolationSet(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMISSField_CreateFinish(RyRDenseField,Err)

  !set up timelag field
  CALL CMISSField_Initialise(RyRReleaseLagField,Err)
  CALL CMISSField_CreateStart(RyRReleaseLagFieldUserNumber,Region,RyRReleaseLagField,Err)
  CALL CMISSField_TypeSet(RyRReleaseLagField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(RyRReleaseLagField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(RyRReleaseLagField,GeometricField,Err)
  CALL CMISSField_NumberOfVariablesSet(RyRReleaseLagField,1,Err)
  CALL CMISSField_VariableTypesSet(RyRReleaseLagField,[CMISS_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMISSField_DataTypeSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DP_TYPE,Err)
  CALL CMISSField_DimensionSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMISSField_VariableLabelSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,"RyR Release Lag Field",Err)
  CALL CMISSField_ComponentMeshComponentGet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMISSField_ComponentMeshComponentSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)
  !Specify the interpolation to be same as geometric interpolation
  CALL CMISSField_ComponentInterpolationSet(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMISSField_CreateFinish(RyRReleaseLagField,Err)

  !Initialise RyR intensity Field, multiply the intensity by number of ryrs per cluster
  ! also Set RyR time lag Field

  NonZeroNodes = 0
  DO node = 1,NUMBER_OF_NODES
    NODE_NUMBER = NodeNums(node,1)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then don't release from there. 
        CALL CMISSField_ParameterSetUpdateNode(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0.0_CMISSDP,Err)
        CALL CMISSField_ParameterSetUpdateNode(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0.0_CMISSDP,Err)

      ELSE
        CALL CMISSField_ParameterSetUpdateNode(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(NumRyRsPerCluster*RyRDensity(node,1)),Err)
        CALL CMISSField_ParameterSetUpdateNode(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(RyRDensity(node,2)),Err)

      ENDIF
      IF(RyRDensity(node,1).GE.0.1_CMISSDP) THEN
        NonZeroNodes=NonZeroNodes+1
      ENDIF
    ENDIF
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  CALL CMISSField_ParameterSetUpdateStart(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  WRITE(*,*) 'Number of Non-Zero RyR Nodes',NonZeroNodes

  !Set up the fields for the other buffers which will store concentrations of the Ca-Buffer complex
!###################
  !F equations
!###################
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
  !Initialise F dependent field
  CALL CMISSField_ComponentValuesInitialise(FField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_F,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL CMISSField_ParameterSetUpdateNode(FField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initF,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(FField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF

  !Create the equations set material field variables - F
  CALL CMISSField_Initialise(FMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(FEquationsSet,FMaterialsFieldUserNumber,FMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"F Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(FEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoFDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoFDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoFDiffz,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,fDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,fDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,fDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 1,fDiffx,Err) !f diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 2,fDiffy,Err) !F diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 3,fDiffz,Err) !f diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient 

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

!###################
  !FCa equations
!###################
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
  !Initialise FCa dependent field
  CALL CMISSField_ComponentValuesInitialise(FCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_FCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL CMISSField_ParameterSetUpdateNode(FCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initFCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(FCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - FCa
  CALL CMISSField_Initialise(FCaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(FCaEquationsSet,FCaMaterialsFieldUserNumber,FCaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"FCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, & 
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(FCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoFCaDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoFCaDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoFCaDiffz,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,fcaDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,fcaDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,fcaDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 1,fcaDiffx,Err) !fca diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 2,fcaDiffy,Err) !fca diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(FCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 3,fcaDiffz,Err) !fca diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient 

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


  !###################
  !CaM equations
  !###################
  CALL CMISSEquationsSet_Initialise(CaMEquationsSet,Err)
  CALL CMISSField_Initialise(CaMEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(CaMEquationsSetUserNumber,Region, &
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & CaMEquationsSetFieldUserNumber,CaMEquationsSetField,CaMEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(CaMEquationsSet,Err)


  !Create the equations set dependent field variables for CaM
  CALL CMISSField_Initialise(CaMField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(CaMEquationsSet,CaMFieldUserNumber,CaMField,Err)
  CALL CMISSField_VariableLabelSet(CaMField,CMISS_FIELD_U_VARIABLE_TYPE,"CaM Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(CaMEquationsSet,Err)
  !Initialise CaM dependent field
  CALL CMISSField_ComponentValuesInitialise(CaMField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_CaM,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL CMISSField_ParameterSetUpdateNode(CaMField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaM,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaMField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - CaM
  CALL CMISSField_Initialise(CaMMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(CaMEquationsSet,CaMMaterialsFieldUserNumber,CaMMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"CaM Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(CaMEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaMDiffx,Err)
        CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaMDiffy,Err)
        CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaMDiffz,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,CaMDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,CaMDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
           & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,CaMDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,CaMDiffx,Err) !CaM diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,CaMDiffy,Err) !CaM diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(CaMMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,CaMDiffz,Err) !CaM diff coeff in z
    ENDIF
    CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 4,store_coeff,Err) ! storage coefficient

   !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
   !iCaMField
   CALL CMISSField_Initialise(iCaMField,Err)
   CALL CMISSEquationsSet_SourceCreateStart(CaMEquationsSet,iCaMFieldUserNumber,iCaMField,Err)
   CALL CMISSField_VariableLabelSet(iCaMField,CMISS_FIELD_U_VARIABLE_TYPE,"iCaM Field",Err)
   !Finish the equations set source field variables
   CALL CMISSEquationsSet_SourceCreateFinish(CaMEquationsSet,Err)
   !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
   CALL CMISSField_ComponentValuesInitialise(iCaMField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
     & 1,0.0_CMISSDP,Err)


  !###################
  !CaMCa equations
  !###################
  CALL CMISSEquationsSet_Initialise(CaMCaEquationsSet,Err)
  CALL CMISSField_Initialise(CaMCaEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(CaMCaEquationsSetUserNumber,Region, &
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & CaMCaEquationsSetFieldUserNumber,CaMCaEquationsSetField,CaMCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(CaMCaEquationsSet,Err)


  !Create the equations set dependent field variables for CaMCa
  CALL CMISSField_Initialise(CaMCaField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(CaMCaEquationsSet,CaMCaFieldUserNumber,CaMCaField,Err)
  CALL CMISSField_VariableLabelSet(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,"CaMCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(CaMCaEquationsSet,Err)
  !Initialise CaMCa dependent field
  CALL CMISSField_ComponentValuesInitialise(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_CaMCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL CMISSField_ParameterSetUpdateNode(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaMCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - CaMCa
  CALL CMISSField_Initialise(CaMCaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(CaMCaEquationsSet,CaMCaMaterialsFieldUserNumber,CaMCaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"CaMCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(CaMCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
        !element based assignment of diffusion properties
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaMCaDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaMCaDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaMCaDiffz,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,CaMCaDiffx,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,CaMCaDiffy,Err)
          CALL CMISSField_ParameterSetUpdateElement(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,CaMCaDiffz,Err)
        ENDIF
      ENDIF
     ENDDO
     CALL CMISSField_ParameterSetUpdateStart(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
     CALL CMISSField_ParameterSetUpdateFinish(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,CaMCaDiffx,Err) !CaMCa diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,CaMCaDiffy,Err) !CaMCa diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(CaMCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,CaMCaDiffz,Err) !CaMCa diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaMCaField
  CALL CMISSField_Initialise(iCaMCaField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(CaMCaEquationsSet,iCaMCaFieldUserNumber,iCaMCaField,Err)
  CALL CMISSField_VariableLabelSet(iCaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,"iCaMCa Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(CaMCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iCaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)


  !###################
  !ATP equations
  !###################
  CALL CMISSEquationsSet_Initialise(ATPEquationsSet,Err)
  CALL CMISSField_Initialise(ATPEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(ATPEquationsSetUserNumber,Region, &
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & ATPEquationsSetFieldUserNumber,ATPEquationsSetField,ATPEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(ATPEquationsSet,Err)


  !Create the equations set dependent field variables for ATP
  CALL CMISSField_Initialise(ATPField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(ATPEquationsSet,ATPFieldUserNumber,ATPField,Err)
  CALL CMISSField_VariableLabelSet(ATPField,CMISS_FIELD_U_VARIABLE_TYPE,"ATP Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(ATPEquationsSet,Err)
  !Initialise ATP dependent field
  CALL CMISSField_ComponentValuesInitialise(ATPField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_ATP,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL CMISSField_ParameterSetUpdateNode(ATPField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initATP,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(ATPField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(ATPField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - ATP
  CALL CMISSField_Initialise(ATPMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(ATPEquationsSet,ATPMaterialsFieldUserNumber,ATPMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"ATP Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(ATPEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoATPDiffx,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoATPDiffy,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoATPDiffz,Err)
      ELSE
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,ATPDiffx,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,ATPDiffy,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,ATPDiffz,Err)
      ENDIF
    ENDIF
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
  CALL CMISSField_ComponentValuesInitialise(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,ATPDiffx,Err) !ATP diff coeff in x
  CALL CMISSField_ComponentValuesInitialise(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,ATPDiffy,Err) !ATP diff coeff in y
  CALL CMISSField_ComponentValuesInitialise(ATPMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,ATPDiffz,Err) !ATP diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPField
  CALL CMISSField_Initialise(iATPField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(ATPEquationsSet,iATPFieldUserNumber,iATPField,Err)
  CALL CMISSField_VariableLabelSet(iATPField,CMISS_FIELD_U_VARIABLE_TYPE,"iATP Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(ATPEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iATPField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)


  !###################
  !ATPCa equations
  !###################
  CALL CMISSEquationsSet_Initialise(ATPCaEquationsSet,Err)
  CALL CMISSField_Initialise(ATPCaEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(ATPCaEquationsSetUserNumber,Region, &
    & GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE, &
    & ATPCaEquationsSetFieldUserNumber,ATPCaEquationsSetField,ATPCaEquationsSet,Err)
   !Set the equations set to be a standard Diffusion no source problem
   !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(ATPCaEquationsSet,Err)


  !Create the equations set dependent field variables for ATPCa
  CALL CMISSField_Initialise(ATPCaField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(ATPCaEquationsSet,ATPCaFieldUserNumber,ATPCaField,Err)
  CALL CMISSField_VariableLabelSet(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,"ATPCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(ATPCaEquationsSet,Err)
  !Initialise ATPCa dependent field
  CALL CMISSField_ComponentValuesInitialise(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,init_ATPCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL CMISSField_ParameterSetUpdateNode(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initATPCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - ATPCa
  CALL CMISSField_Initialise(ATPCaMaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(ATPCaEquationsSet,ATPCaMaterialsFieldUserNumber,ATPCaMaterialsField,Err)
  CALL CMISSField_VariableLabelSet(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,"ATPCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL CMISSField_ComponentInterpolationSet(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,2, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,3, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(ATPCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL CMISSDecomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoATPCaDiffx,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoATPCaDiffy,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoATPCaDiffz,Err)
      ELSE
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,ATPCaDiffx,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,ATPCaDiffy,Err)
        CALL CMISSField_ParameterSetUpdateElement(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,ATPCaDiffz,Err)
      ENDIF
    ENDIF
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL CMISSField_ComponentValuesInitialise(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,ATPCaDiffx,Err) !ATPCa diff coeff in x
    CALL CMISSField_ComponentValuesInitialise(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 2,ATPCaDiffy,Err) !ATPCa diff coeff in y
    CALL CMISSField_ComponentValuesInitialise(ATPCaMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 3,ATPCaDiffz,Err) !ATPCa diff coeff in z
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(FMaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPCaField
  CALL CMISSField_Initialise(iATPCaField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(ATPCaEquationsSet,iATPCaFieldUserNumber,iATPCaField,Err)
  CALL CMISSField_VariableLabelSet(iATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,"iATPCa Field",Err)
  !Finish the equations set source field variables
  CALL CMISSEquationsSet_SourceCreateFinish(ATPCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMISSField_ComponentValuesInitialise(iATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)


!###################
  !CaTnC
!###################
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
  IF(WITH_MITO_ELEMENTS.EQ.0) THEN
    CALL CMISSField_ComponentValuesInitialise(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & 1,init_CaTnC,Err)
  ELSE
    DO node=1,NUMBER_OF_NODES
      NODE_NUMBER=NodeNums(node,1)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(node,2).EQ.MITO_BD_MARKER .OR. NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
          CALL CMISSField_ParameterSetUpdateNode(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaTnC,Err)
        ELSE
          CALL CMISSField_ParameterSetUpdateNode(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE, &
            & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaTnC,Err)
        ENDIF
      ENDIF
    ENDDO  
    CALL CMISSField_ParameterSetUpdateStart(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  ENDIF  
!_____________________________________________________________________________________________________________________
  !Start to set up CellML Fields

  !Create the CellML environment
  CALL CMISSCellML_Initialise(CellML,Err)
  CALL CMISSCellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import ryr release and buffer source model from a file
  CALL CMISSCellML_ModelImport(CellML,RyRModel,ryrModelIndex,Err)
  ! set iCa as known so that it can be set as spatially varying in opencmiss.
  !CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  ! set RyRDensity as known so that it can be set as spatially varying in opencmiss.
  CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/ryrDensity",Err)
  CALL CMISSCellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/timelag",Err)

  !to get from the CellML side. variables in cellml model that are not state variables, but are dependent on independent and state variables. 
  !- components of intermediate field
  !fluxes of the different buffers and CaRUs that I want to get out as intermediate variables
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"CRU/Jryr",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"FluoBuffer/Jfluo",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"TnCBuffer/Jtnc",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"TnCBuffer/JATP",Err)
  CALL CMISSCellML_VariableSetAsWanted(CellML,ryrModelIndex,"TnCBuffer/JCaM",Err)

  !Finish the CellML environment
  CALL CMISSCellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  !Mapping free calcium in opencmiss to that in cellml.
  CALL CMISSCellML_FieldMapsCreateStart(CellML,Err)
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/Ca_free",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/Ca_free",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping iCaField to iCa in the cellml model
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,iCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/iCa",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/iCa",CMISS_FIELD_VALUES_SET_TYPE, &
    & iCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping RyRDenseField to RyRDensity in the cellml model
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/ryrDensity",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/ryrDensity",CMISS_FIELD_VALUES_SET_TYPE, &
    & RyRDenseField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping RyRReleaseLagField to timelag in the cellml model
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/timelag",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/timelag",CMISS_FIELD_VALUES_SET_TYPE, &
    & RyRReleaseLagField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping Buffer-Complex resting values of cellml model to appropriate fields set up above

   !Mapping F
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,FField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"FluoBuffer/Fluo_free",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"FluoBuffer/Fluo_free",CMISS_FIELD_VALUES_SET_TYPE, &
    & FField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping FCa
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,FCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"FluoBuffer/FluoCa",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"FluoBuffer/FluoCa",CMISS_FIELD_VALUES_SET_TYPE, &
    & FCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

   !Mapping CaTnC
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"TnCBuffer/CaTnC",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"TnCBuffer/CaTnC",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaTnCField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaM
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaMField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CaMBuffer/CaM",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CaMBuffer/CaM",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaMField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaMCa
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CaMBuffer/CaMCa",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CaMBuffer/CaMCa",CMISS_FIELD_VALUES_SET_TYPE, &
    & CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATP
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,ATPField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"ATPBuffer/ATP",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"ATPBuffer/ATP",CMISS_FIELD_VALUES_SET_TYPE, &
    & ATPField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATPCa
  CALL CMISSCellML_CreateFieldToCellMLMap(CellML,ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"ATPBuffer/ATPCa",CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSCellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"ATPBuffer/ATPCa",CMISS_FIELD_VALUES_SET_TYPE, &
    & ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_VALUES_SET_TYPE,Err)

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
  IF(MODEL_ON.EQ.1) THEN 
    CALL CMISSField_ComponentValuesInitialise(CellMLModelsField, & 
      & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)
    IF(WITH_MITO_ELEMENTS.EQ.1) THEN
      DO node=1,NUMBER_OF_NODES
        NODE_NUMBER=NodeNums(node,1)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
            CALL CMISSField_ParameterSetUpdateNode(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE, &
              & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0_CMISSIntg,Err)
          ENDIF
        ENDIF
      ENDDO  
      CALL CMISSField_ParameterSetUpdateStart(CellMLModelsField, &
        & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
      CALL CMISSField_ParameterSetUpdateFinish(CellMLModelsField, &
        & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    ENDIF
  ELSE
    CALL CMISSField_ComponentValuesInitialise(CellMLModelsField, & 
      & CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0_CMISSIntg,Err)
  ENDIF
  CALL CMISSField_ParameterSetUpdateStart(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CellMLModelsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  !Start the creation of the CellML state field
  CALL CMISSField_Initialise(CellMLStateField,Err)
  CALL CMISSCellML_StateFieldCreateStart(CellML, &
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellML_StateFieldCreateFinish(CellML,Err)
  CALL CMISSField_ParameterSetUpdateStart(CellMLStateField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CellMLStateField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)


  !Start the creation of the CellML intermediate field
  CALL CMISSField_Initialise(CellMLIntermediateField,Err)
  CALL CMISSCellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMISSCellML_IntermediateFieldCreateFinish(CellML,Err)
  CALL CMISSField_ParameterSetUpdateStart(CellMLIntermediateField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CellMLIntermediateField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)


  !Start the creation of CellML parameters field
  CALL CMISSField_Initialise(CellMLParametersField,Err)
  CALL CMISSCellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellML_ParametersFieldCreateFinish(CellML,Err)

  CALL CMISSField_ParameterSetUpdateStart(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(CellMLParametersField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

!_____________________________________________________________________________________________________________
 
  !Create the equations set equations for Ca
  CALL CMISSEquations_Initialise(CaEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(CaEquationsSet,CaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(CaEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
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
  CALL CMISSEquations_SparsityTypeSet(FEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
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
  CALL CMISSEquations_SparsityTypeSet(FCaEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(FCaEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(FCaEquationsSet,Err)

  !Create the equations set equations for CaM
  CALL CMISSEquations_Initialise(CaMEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(CaMEquationsSet,CaMEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(CaMEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(CaMEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(CaMEquationsSet,Err)

  !Create the equations set equations for CaMCa
  CALL CMISSEquations_Initialise(CaMCaEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(CaMCaEquationsSet,CaMCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(CaMCaEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(CaMCaEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(CaMCaEquationsSet,Err)

  !Create the equations set equations for ATP
  CALL CMISSEquations_Initialise(ATPEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(ATPEquationsSet,ATPEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(ATPEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(ATPEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(ATPEquationsSet,Err)

  !Create the equations set equations for ATPCa
  CALL CMISSEquations_Initialise(ATPCaEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(ATPCaEquationsSet,ATPCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(ATPCaEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(ATPCaEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(ATPCaEquationsSet,Err)

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
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,outputfreq,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,5.00_CMISSDP,0.01_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

!______________________________________________________________________________________________________________
  !Set up the problem solvers for Strang splitting - 
  !note, this example is contrived to have strang splitting, when it could be solved as a simple evaluation (as opposed to integration) of source and diffusion
  WRITE(*,*) 'Set up the problem solvers for Strang splitting'
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)

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
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

!_______________________________________________________________________________________________________________
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

!_______________________________________________________________________________________________________________
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
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set for Ca, F and FCa
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,CaEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,FEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,FCaEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,CaMEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,CaMCaEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,ATPEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,ATPCaEquationsSet,EquationsSetIndex,Err)

  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)
!_________________________________________________________________________________________________________
  WRITE(*,*) 'Set up boundary conditions'  
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Set 0 conc. bc on nodes within mitos
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO node=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(node,1)
      IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CONDITION = CMISS_BOUNDARY_CONDITION_FIXED
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaField, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FField, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FCaField, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)


          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)


          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

        ENDIF
      ENDIF
    ENDDO
  ENDIF
  !set no flux bc if node is a mito boundary node.

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO node=1,NUMBER_OF_MITOBDFACENODES
      NODE_NUMBER = MITOBDFaceNodes(node)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CONDITION = CMISS_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSDP
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)!

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FCaField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        !make sure to unset mito bd node bcs to be free from dirchlet bcs set above.

        CONDITION = CMISS_BOUNDARY_CONDITION_FREE
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_Ca,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_F,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,FCaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_FCa,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_CaM,Err) !(neumann boundary condition - no flux)

         CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_CaMCa,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_ATP,Err) !(neumann boundary condition - no flux)

        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_ATPCa,Err) !(neumann boundary condition - no flux)



        !Set the initial conc. at these mito boundary nodes to be the cytosolic versions
        CALL CMISSField_ParameterSetUpdateNode(CaField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_Ca,Err)

        CALL CMISSField_ParameterSetUpdateNode(FField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_F,Err)

        CALL CMISSField_ParameterSetUpdateNode(FCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_FCa,Err)

        CALL CMISSField_ParameterSetUpdateNode(CaMField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaM,Err)


        CALL CMISSField_ParameterSetUpdateNode(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaMCa,Err)

        CALL CMISSField_ParameterSetUpdateNode(ATPField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_ATP,Err)

        CALL CMISSField_ParameterSetUpdateNode(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_ATPCa,Err)


      ENDIF
    ENDDO
    CALL CMISSField_ParameterSetUpdateStart(CaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(FField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(FCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(FCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(CaMField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(CaMCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(ATPField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(ATPField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    CALL CMISSField_ParameterSetUpdateStart(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
    CALL CMISSField_ParameterSetUpdateFinish(ATPCaField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

    ENDIF

  !Set no flux on cell boundary
  DO node=1,NUMBER_OF_CELLBDNODES
    NODE_NUMBER = CELLBDNodes(node)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
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

      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMField, &
      & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
      & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPField, &
      & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
      & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

    ENDIF
  ENDDO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!__________________________________________________________________________________________________________
  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)
!__________________________________________________________________________________________________________
  IF(EXPORT_FIELD) THEN
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"Cell_Solution","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"Cell_Solution","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
  ENDIF 
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP


END PROGRAM CARDIAC_ECC
