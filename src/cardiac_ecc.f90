!> \file
!> $Id: cardiac_ecc.f90 2014-12-16 vraj004 $
!> \author Vijay Rajagopal
!> \brief Main program file to simulate cardiac ECC using opencmfe_ library routines
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
!> The Original Code is Opencmfe_
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
  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif



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


  !Defining cmfe_-type variables
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField,CaMaterialsField,FCaMaterialsField,FMaterialsField,CaField,FCaField,FField
  TYPE(cmfe_FieldType) :: CaEquationsSetField,FCaEquationsSetField,FEquationsSetField
  TYPE(cmfe_FieldType) :: CaMMaterialsField, CaMEquationsSetField,CaMField
  TYPE(cmfe_FieldType) :: CaMCaMaterialsField, CaMCaEquationsSetField,CaMCaField
  TYPE(cmfe_FieldType) :: ATPField,ATPMaterialsField,ATPEquationsSetField
  TYPE(cmfe_FieldType) :: ATPCaField,ATPCaMaterialsField,ATPCaEquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: MeshElements
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_EquationsType) :: CaEquations,FCaEquations,FEquations,CaMEquations,CaMCaEquations,ATPEquations
  TYPE(cmfe_EquationsType) :: ATPCaEquations
  TYPE(cmfe_EquationsSetType) :: CaEquationsSet,FCaEquationsSet,FEquationsSet,ATPEquationsSet,CaMEquationsSet
  TYPE(cmfe_EquationsSetType) :: ATPCaEquationsSet,CaMCaEquationsSet
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldType) :: iCaField,CaTnCField,RyRDenseField,iFCaField,iFField,RyRReleaseLagField
  TYPE(cmfe_FieldType) :: iCaMField,iCaMCaField,iATPField,iATPCaField

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
  !Intialise Opencmfe_
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(cmfe_ERRORS_TRAP_ERROR,Err)
  !get computational nodes for parallel processing
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)


  !set diagnostics
  !CALL cmfe_DiagnosticsSetOn(cmfe_FROM_DIAG_TYPE,(/1,2,3,4,5/),"SOLVE_DIAGNOSTICS", &
  !  & (/"REACTION_DIFFUSION_PRE_SOLVE"/),Err)
  CALL MPI_BCAST(NumberOfComputationalNodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)


   !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystemUserNumber,3,Err)
  !The coordinate system is 3D by default;set it to be 3D in above command.
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_LabelSet(Region,"Cell",Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

!_________________________________________________________________________________________________
  !Start the creation of a trilinear-simplex basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a trilinear simplex  basis
  CALL cmfe_Basis_TypeSet(Basis,cmfe_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
  !set interpolation to be linear
  CALL cmfe_Basis_InterpolationXiSet(Basis,(/cmfe_Basis_Linear_Simplex_Interpolation, &
   &   cmfe_Basis_Linear_Simplex_Interpolation, cmfe_Basis_Linear_Simplex_Interpolation/),Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

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

  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL cmfe_MeshElements_Initialise(MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    CALL cmfe_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
     &   ElemMap(i,4),ElemMap(i,5)/),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,cmfe_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 3 field components in 1 mesh component
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      nodez = NodeCoords(i,3)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,3,nodez,Err)
     ENDIF
    ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF 
!______________________________________________________________________________________________________________
  !Create the cellml reaction with split reaction diffusion equations_set - 1 for each species
!###################
  !Ca equations
!###################
  CALL cmfe_EquationsSet_Initialise(CaEquationsSet,Err)
  CALL cmfe_Field_Initialise(CaEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(CaEquationsSetUserNumber,Region, & 
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaEquationsSetFieldUserNumber,CaEquationsSetField,CaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(CaEquationsSet,Err)


  !Create the equations set dependent field variables for Ca
  CALL cmfe_Field_Initialise(CaField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(CaEquationsSet,CaFieldUserNumber,CaField,Err)
  CALL cmfe_Field_VariableLabelSet(CaField,cmfe_FIELD_U_VARIABLE_TYPE,"Ca Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(CaEquationsSet,Err)
  !Initialise Ca dependent field
  CALL cmfe_Field_ComponentValuesInitialise(CaField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_Ca,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL cmfe_Field_ParameterSetUpdateNode(CaField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF

  !Create the equations set material field variables - Ca
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL cmfe_Field_Initialise(CaMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(CaEquationsSet,CaMaterialsFieldUserNumber,CaMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"Ca Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(CaEquationsSet,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaDiffz,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,caDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,caDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,caDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 1,caDiffx,Err) !ca diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 2,caDiffy,Err) !ca diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 3,caDiffz,Err) !ca diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(CaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaField
  !Might use the field for CellML input of elementary RyR calcium release
  CALL cmfe_Field_Initialise(iCaField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(CaEquationsSet,iCaFieldUserNumber,iCaField,Err)
  CALL cmfe_Field_VariableLabelSet(iCaField,cmfe_FIELD_U_VARIABLE_TYPE,"iCa Field",Err)
  !CALL cmfe_Field_ComponentInterpolationSet(iCaField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
  !  & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(CaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
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
PRINT *,'______________________________________________________________________________'
  !set up intensity field
  CALL cmfe_Field_Initialise(RyRDenseField,Err)
  CALL cmfe_Field_CreateStart(RyRDenseFieldUserNumber,Region,RyRDenseField,Err)
  CALL cmfe_Field_TypeSet(RyRDenseField,cmfe_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(RyRDenseField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(RyRDenseField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(RyRDenseField,1,Err)
  CALL cmfe_Field_VariableTypesSet(RyRDenseField,[cmfe_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,"RyR Density Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & cmfe_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL cmfe_Field_CreateFinish(RyRDenseField,Err)

  !set up timelag field
  CALL cmfe_Field_Initialise(RyRReleaseLagField,Err)
  CALL cmfe_Field_CreateStart(RyRReleaseLagFieldUserNumber,Region,RyRReleaseLagField,Err)
  CALL cmfe_Field_TypeSet(RyRReleaseLagField,cmfe_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(RyRReleaseLagField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(RyRReleaseLagField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(RyRReleaseLagField,1,Err)
  CALL cmfe_Field_VariableTypesSet(RyRReleaseLagField,[cmfe_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,"RyR Release Lag Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & cmfe_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL cmfe_Field_CreateFinish(RyRReleaseLagField,Err)

  !Initialise RyR intensity Field, multiply the intensity by number of ryrs per cluster
  ! also Set RyR time lag Field

  NonZeroNodes = 0
  DO node = 1,NUMBER_OF_NODES
    NODE_NUMBER = NodeNums(node,1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then don't release from there. 
        CALL cmfe_Field_ParameterSetUpdateNode(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0.0_CMISSDP,Err)
        CALL cmfe_Field_ParameterSetUpdateNode(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0.0_CMISSDP,Err)

      ELSE
        CALL cmfe_Field_ParameterSetUpdateNode(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(NumRyRsPerCluster*RyRDensity(node,1)),Err)
        CALL cmfe_Field_ParameterSetUpdateNode(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(RyRDensity(node,2)),Err)

      ENDIF
      IF(RyRDensity(node,1).GE.0.1_CMISSDP) THEN
        NonZeroNodes=NonZeroNodes+1
      ENDIF
    ENDIF
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_Field_ParameterSetUpdateStart(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  WRITE(*,*) 'Number of Non-Zero RyR Nodes',NonZeroNodes

  !Set up the fields for the other buffers which will store concentrations of the Ca-Buffer complex
!###################
  !F equations
!###################
  CALL cmfe_EquationsSet_Initialise(FEquationsSet,Err)
PRINT *,'______________________________________________________________________________'
  CALL cmfe_Field_Initialise(FEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(FEquationsSetUserNumber,Region, & 
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & FEquationsSetFieldUserNumber,FEquationsSetField,FEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(FEquationsSet,Err)

PRINT *,'______________________________________________________________________________'
  !Create the equations set dependent field variables for F
  CALL cmfe_Field_Initialise(FField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(FEquationsSet,FFieldUserNumber,FField,Err)
  CALL cmfe_Field_VariableLabelSet(FField,cmfe_FIELD_U_VARIABLE_TYPE,"F Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(FEquationsSet,Err)
  !Initialise F dependent field
  CALL cmfe_Field_ComponentValuesInitialise(FField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_F,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL cmfe_Field_ParameterSetUpdateNode(FField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initF,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(FField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF
PRINT *,'______________________________________________________________________________'
  !Create the equations set material field variables - F
  CALL cmfe_Field_Initialise(FMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(FEquationsSet,FMaterialsFieldUserNumber,FMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"F Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(FEquationsSet,Err)
PRINT *,'______________________________________________________________________________'
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoFDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoFDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoFDiffz,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,fDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,fDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,fDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 1,fDiffx,Err) !f diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 2,fDiffy,Err) !F diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 3,fDiffz,Err) !f diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient 
PRINT *,'______________________________________________________________________________'
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFField
  CALL cmfe_Field_Initialise(iFField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(FEquationsSet,iFFieldUserNumber,iFField,Err)
  CALL cmfe_Field_VariableLabelSet(iFField,cmfe_FIELD_U_VARIABLE_TYPE,"iF Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(FEquationsSet,Err)
  !Initialising the iField to zero everywhere. Might modify for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iFField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)
PRINT *,'F equations'
!###################
  !FCa equations
!###################
  CALL cmfe_EquationsSet_Initialise(FCaEquationsSet,Err)
  CALL cmfe_Field_Initialise(FCaEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(FCaEquationsSetUserNumber,Region, & 
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & FCaEquationsSetFieldUserNumber,FCaEquationsSetField,FCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(FCaEquationsSet,Err)


  !Create the equations set dependent field variables for FCa
  CALL cmfe_Field_Initialise(FCaField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(FCaEquationsSet,FCaFieldUserNumber,FCaField,Err)
  CALL cmfe_Field_VariableLabelSet(FCaField,cmfe_FIELD_U_VARIABLE_TYPE,"FCa Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(FCaEquationsSet,Err)
  !Initialise FCa dependent field
  CALL cmfe_Field_ComponentValuesInitialise(FCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_FCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init. 
          CALL cmfe_Field_ParameterSetUpdateNode(FCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initFCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(FCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - FCa
  CALL cmfe_Field_Initialise(FCaMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(FCaEquationsSet,FCaMaterialsFieldUserNumber,FCaMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"FCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, & 
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(FCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
          !element based assignment of diffusion properties
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoFCaDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoFCaDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoFCaDiffz,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,fcaDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,fcaDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,fcaDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 1,fcaDiffx,Err) !fca diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 2,fcaDiffy,Err) !fca diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(FCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 3,fcaDiffz,Err) !fca diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient 

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFCaField
  CALL cmfe_Field_Initialise(iFCaField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(FCaEquationsSet,iFCaFieldUserNumber,iFCaField,Err)
  CALL cmfe_Field_VariableLabelSet(iFCaField,cmfe_FIELD_U_VARIABLE_TYPE,"iFCa Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(FCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iFCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)
PRINT *,'FCa equations'

  !###################
  !CaM equations
  !###################
  CALL cmfe_EquationsSet_Initialise(CaMEquationsSet,Err)
  CALL cmfe_Field_Initialise(CaMEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(CaMEquationsSetUserNumber,Region, &
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaMEquationsSetFieldUserNumber,CaMEquationsSetField,CaMEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(CaMEquationsSet,Err)


  !Create the equations set dependent field variables for CaM
  CALL cmfe_Field_Initialise(CaMField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(CaMEquationsSet,CaMFieldUserNumber,CaMField,Err)
  CALL cmfe_Field_VariableLabelSet(CaMField,cmfe_FIELD_U_VARIABLE_TYPE,"CaM Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(CaMEquationsSet,Err)
  !Initialise CaM dependent field
  CALL cmfe_Field_ComponentValuesInitialise(CaMField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_CaM,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL cmfe_Field_ParameterSetUpdateNode(CaMField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaM,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaMField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - CaM
  CALL cmfe_Field_Initialise(CaMMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(CaMEquationsSet,CaMMaterialsFieldUserNumber,CaMMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"CaM Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(CaMEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaMDiffx,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaMDiffy,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaMDiffz,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,CaMDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,CaMDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
           & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,CaMDiffz,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 1,CaMDiffx,Err) !CaM diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 2,CaMDiffy,Err) !CaM diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(CaMMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 3,CaMDiffz,Err) !CaM diff coeff in z
    ENDIF
    CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 4,store_coeff,Err) ! storage coefficient

   !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
   !iCaMField
   CALL cmfe_Field_Initialise(iCaMField,Err)
   CALL cmfe_EquationsSet_SourceCreateStart(CaMEquationsSet,iCaMFieldUserNumber,iCaMField,Err)
   CALL cmfe_Field_VariableLabelSet(iCaMField,cmfe_FIELD_U_VARIABLE_TYPE,"iCaM Field",Err)
   !Finish the equations set source field variables
   CALL cmfe_EquationsSet_SourceCreateFinish(CaMEquationsSet,Err)
   !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
   CALL cmfe_Field_ComponentValuesInitialise(iCaMField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
     & 1,0.0_CMISSDP,Err)

PRINT *,'CaM equations'
  !###################
  !CaMCa equations
  !###################
  CALL cmfe_EquationsSet_Initialise(CaMCaEquationsSet,Err)
  CALL cmfe_Field_Initialise(CaMCaEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(CaMCaEquationsSetUserNumber,Region, &
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaMCaEquationsSetFieldUserNumber,CaMCaEquationsSetField,CaMCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(CaMCaEquationsSet,Err)


  !Create the equations set dependent field variables for CaMCa
  CALL cmfe_Field_Initialise(CaMCaField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(CaMCaEquationsSet,CaMCaFieldUserNumber,CaMCaField,Err)
  CALL cmfe_Field_VariableLabelSet(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,"CaMCa Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(CaMCaEquationsSet,Err)
  !Initialise CaMCa dependent field
  CALL cmfe_Field_ComponentValuesInitialise(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_CaMCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL cmfe_Field_ParameterSetUpdateNode(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaMCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - CaMCa
  CALL cmfe_Field_Initialise(CaMCaMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(CaMCaEquationsSet,CaMCaMaterialsFieldUserNumber,CaMCaMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"CaMCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(CaMCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

    DO i=1,NUMBER_OF_ELEMENTS
      ELEM_NUMBER = ElemMap(i,1)
      CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
      IF(ElementDomain==ComputationalNodeNumber) THEN
        ELEM_LABEL = ElemMap(i,6)
        IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
        !element based assignment of diffusion properties
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoCaMCaDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoCaMCaDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoCaMCaDiffz,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,CaMCaDiffx,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,CaMCaDiffy,Err)
          CALL cmfe_Field_ParameterSetUpdateElement(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,CaMCaDiffz,Err)
        ENDIF
      ENDIF
     ENDDO
     CALL cmfe_Field_ParameterSetUpdateStart(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
     CALL cmfe_Field_ParameterSetUpdateFinish(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 1,CaMCaDiffx,Err) !CaMCa diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 2,CaMCaDiffy,Err) !CaMCa diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(CaMCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 3,CaMCaDiffz,Err) !CaMCa diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaMCaField
  CALL cmfe_Field_Initialise(iCaMCaField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(CaMCaEquationsSet,iCaMCaFieldUserNumber,iCaMCaField,Err)
  CALL cmfe_Field_VariableLabelSet(iCaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,"iCaMCa Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(CaMCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iCaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

PRINT *,'CaMCa equations'
  !###################
  !ATP equations
  !###################
  CALL cmfe_EquationsSet_Initialise(ATPEquationsSet,Err)
  CALL cmfe_Field_Initialise(ATPEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(ATPEquationsSetUserNumber,Region, &
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ATPEquationsSetFieldUserNumber,ATPEquationsSetField,ATPEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(ATPEquationsSet,Err)


  !Create the equations set dependent field variables for ATP
  CALL cmfe_Field_Initialise(ATPField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(ATPEquationsSet,ATPFieldUserNumber,ATPField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,"ATP Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(ATPEquationsSet,Err)
  !Initialise ATP dependent field
  CALL cmfe_Field_ComponentValuesInitialise(ATPField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_ATP,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL cmfe_Field_ParameterSetUpdateNode(ATPField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initATP,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - ATP
  CALL cmfe_Field_Initialise(ATPMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(ATPEquationsSet,ATPMaterialsFieldUserNumber,ATPMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"ATP Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(ATPEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoATPDiffx,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoATPDiffy,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoATPDiffz,Err)
      ELSE
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,ATPDiffx,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,ATPDiffy,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,ATPDiffz,Err)
      ENDIF
    ENDIF
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
  CALL cmfe_Field_ComponentValuesInitialise(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,ATPDiffx,Err) !ATP diff coeff in x
  CALL cmfe_Field_ComponentValuesInitialise(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 2,ATPDiffy,Err) !ATP diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(ATPMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 3,ATPDiffz,Err) !ATP diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPField
  CALL cmfe_Field_Initialise(iATPField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(ATPEquationsSet,iATPFieldUserNumber,iATPField,Err)
  CALL cmfe_Field_VariableLabelSet(iATPField,cmfe_FIELD_U_VARIABLE_TYPE,"iATP Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(ATPEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iATPField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

PRINT *,'ATP equations'
  !###################
  !ATPCa equations
  !###################
  CALL cmfe_EquationsSet_Initialise(ATPCaEquationsSet,Err)
  CALL cmfe_Field_Initialise(ATPCaEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(ATPCaEquationsSetUserNumber,Region, &
    & GeometricField,[cmfe_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & cmfe_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & cmfe_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ATPCaEquationsSetFieldUserNumber,ATPCaEquationsSetField,ATPCaEquationsSet,Err)
   !Set the equations set to be a standard Diffusion no source problem
   !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(ATPCaEquationsSet,Err)


  !Create the equations set dependent field variables for ATPCa
  CALL cmfe_Field_Initialise(ATPCaField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(ATPCaEquationsSet,ATPCaFieldUserNumber,ATPCaField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,"ATPCa Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(ATPCaEquationsSet,Err)
  !Initialise ATPCa dependent field
  CALL cmfe_Field_ComponentValuesInitialise(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
    & cmfe_FIELD_VALUES_SET_TYPE,1,init_ATPCa,Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO i=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(i,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(i,2).EQ.MITO_REGION_MARKER) THEN !if node is mito-associated node, then set initial conc. to mito_init.
          CALL cmfe_Field_ParameterSetUpdateNode(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initATPCa,Err)
        ENDIF
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  ENDIF


  !Create the equations set material field variables - ATPCa
  CALL cmfe_Field_Initialise(ATPCaMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(ATPCaEquationsSet,ATPCaMaterialsFieldUserNumber,ATPCaMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,"ATPCa Materials Field",Err)
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN ! element based assignment of diffusion properties to distinguish mito diffusion from the rest
    CALL cmfe_Field_ComponentInterpolationSet(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,2, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,3, &
      & cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDIF
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(ATPCaEquationsSet,Err)

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN

  DO i=1,NUMBER_OF_ELEMENTS
    ELEM_NUMBER = ElemMap(i,1)
    CALL cmfe_Decomposition_ElementDomainGet(Decomposition,ELEM_NUMBER,ElementDomain,Err)
    IF(ElementDomain==ComputationalNodeNumber) THEN
      ELEM_LABEL = ElemMap(i,6)
      IF(ELEM_LABEL.EQ.MITO_REGION_MARKER) THEN
      !element based assignment of diffusion properties
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,mitoATPCaDiffx,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,mitoATPCaDiffy,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,mitoATPCaDiffz,Err)
      ELSE
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,1,ATPCaDiffx,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,2,ATPCaDiffy,Err)
        CALL cmfe_Field_ParameterSetUpdateElement(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,ELEM_NUMBER,3,ATPCaDiffz,Err)
      ENDIF
    ENDIF
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 1,ATPCaDiffx,Err) !ATPCa diff coeff in x
    CALL cmfe_Field_ComponentValuesInitialise(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 2,ATPCaDiffy,Err) !ATPCa diff coeff in y
    CALL cmfe_Field_ComponentValuesInitialise(ATPCaMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 3,ATPCaDiffz,Err) !ATPCa diff coeff in z
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(FMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 4,store_coeff,Err) ! storage coefficient

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPCaField
  CALL cmfe_Field_Initialise(iATPCaField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(ATPCaEquationsSet,iATPCaFieldUserNumber,iATPCaField,Err)
  CALL cmfe_Field_VariableLabelSet(iATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,"iATPCa Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(ATPCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL cmfe_Field_ComponentValuesInitialise(iATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

PRINT *,'ATPCa equations'
!###################
  !CaTnC
!###################
  CALL cmfe_Field_Initialise(CaTnCField,Err)
  CALL cmfe_Field_CreateStart(CaTnCFieldUserNumber,Region,CaTnCField,Err)
  CALL cmfe_Field_TypeSet(CaTnCField,cmfe_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(CaTnCField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(CaTnCField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(CaTnCField,1,Err)
  CALL cmfe_Field_VariableTypesSet(CaTnCField,[cmfe_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,"CaTnC Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
    & cmfe_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL cmfe_Field_CreateFinish(CaTnCField,Err)
  !Initialise CaTnC concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  IF(WITH_MITO_ELEMENTS.EQ.0) THEN
    CALL cmfe_Field_ComponentValuesInitialise(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
      & 1,init_CaTnC,Err)
  ELSE
    DO node=1,NUMBER_OF_NODES
      NODE_NUMBER=NodeNums(node,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        IF(NodeNums(node,2).EQ.MITO_BD_MARKER .OR. NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
          CALL cmfe_Field_ParameterSetUpdateNode(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,mito_initCaTnC,Err)
        ELSE
          CALL cmfe_Field_ParameterSetUpdateNode(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaTnC,Err)
        ENDIF
      ENDIF
    ENDDO  
    CALL cmfe_Field_ParameterSetUpdateStart(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  ENDIF  
!_____________________________________________________________________________________________________________________
  !Start to set up CellML Fields

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import ryr release and buffer source model from a file
  CALL cmfe_CellML_ModelImport(CellML,RyRModel,ryrModelIndex,Err)
  ! set iCa as known so that it can be set as spatially varying in opencmfe_.
  !CALL cmfe_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  ! set RyRDensity as known so that it can be set as spatially varying in opencmfe_.
  CALL cmfe_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/ryrDensity",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/timelag",Err)

  !to get from the CellML side. variables in cellml model that are not state variables, but are dependent on independent and state variables. 
  !- components of intermediate field
  !fluxes of the different buffers and CaRUs that I want to get out as intermediate variables
  CALL cmfe_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"CRU/Jryr",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"FluoBuffer/Jfluo",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"TnCBuffer/Jtnc",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"ATPBuffer/JATP",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"CaMBuffer/JCaM",Err)

  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> Opencmfe_ field maps
  !Mapping free calcium in opencmfe_ to that in cellml.
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/Ca_free",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/Ca_free",cmfe_FIELD_VALUES_SET_TYPE, &
    & CaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

   !Mapping iCaField to iCa in the cellml model
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,iCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/iCa",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/iCa",cmfe_FIELD_VALUES_SET_TYPE, &
    & iCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

   !Mapping RyRDenseField to RyRDensity in the cellml model
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/ryrDensity",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/ryrDensity",cmfe_FIELD_VALUES_SET_TYPE, &
    & RyRDenseField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping RyRReleaseLagField to timelag in the cellml model
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CRU/timelag",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CRU/timelag",cmfe_FIELD_VALUES_SET_TYPE, &
    & RyRReleaseLagField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping Buffer-Complex resting values of cellml model to appropriate fields set up above

   !Mapping F
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,FField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"FluoBuffer/Fluo_free",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"FluoBuffer/Fluo_free",cmfe_FIELD_VALUES_SET_TYPE, &
    & FField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

   !Mapping FCa
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,FCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"FluoBuffer/FluoCa",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"FluoBuffer/FluoCa",cmfe_FIELD_VALUES_SET_TYPE, &
    & FCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

   !Mapping CaTnC
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"TnCBuffer/CaTnC",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"TnCBuffer/CaTnC",cmfe_FIELD_VALUES_SET_TYPE, &
    & CaTnCField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaM
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CaMField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CaMBuffer/CaM_free",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CaMBuffer/CaM_free",cmfe_FIELD_VALUES_SET_TYPE, &
    & CaMField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaMCa
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"CaMBuffer/CaMCa",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"CaMBuffer/CaMCa",cmfe_FIELD_VALUES_SET_TYPE, &
    & CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATP
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"ATPBuffer/ATP_free",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"ATPBuffer/ATP_free",cmfe_FIELD_VALUES_SET_TYPE, &
    & ATPField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATPCa
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"ATPBuffer/ATPCa",cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"ATPBuffer/ATPCa",cmfe_FIELD_VALUES_SET_TYPE, &
    & ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Finish the creation of CellML <--> Opencmfe_ field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)


  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML, &
    & CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)

  !By default all field parameters have default model value of 1, i.e. the first model. 
  ! assigning the bufferNryr cellml model (model 1) for all nodes.
  IF(MODEL_ON.EQ.1) THEN 
    CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField, & 
      & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)
    IF(WITH_MITO_ELEMENTS.EQ.1) THEN
      DO node=1,NUMBER_OF_NODES
        NODE_NUMBER=NodeNums(node,1)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,cmfe_FIELD_U_VARIABLE_TYPE, &
              & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0_CMISSIntg,Err)
          ENDIF
        ENDIF
      ENDDO  
      CALL cmfe_Field_ParameterSetUpdateStart(CellMLModelsField, &
        & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
      CALL cmfe_Field_ParameterSetUpdateFinish(CellMLModelsField, &
        & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    ENDIF
  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField, & 
      & cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1,0_CMISSIntg,Err)
  ENDIF
  CALL cmfe_Field_ParameterSetUpdateStart(CellMLModelsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CellMLModelsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML, &
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(CellMLStateField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CellMLStateField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)


  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(CellMLIntermediateField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CellMLIntermediateField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)


  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

  CALL cmfe_Field_ParameterSetUpdateStart(CellMLParametersField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(CellMLParametersField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

!_____________________________________________________________________________________________________________
 
  !Create the equations set equations for Ca
  CALL cmfe_Equations_Initialise(CaEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(CaEquationsSet,CaEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(CaEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(CaEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(CaEquationsSet,Err)

  !Create the equations set equations for F
  CALL cmfe_Equations_Initialise(FEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(FEquationsSet,FEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(FEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(FEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(FEquationsSet,Err)

  !Create the equations set equations for FCa
  CALL cmfe_Equations_Initialise(FCaEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(FCaEquationsSet,FCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(FCaEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(FCaEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(FCaEquationsSet,Err)

  !Create the equations set equations for CaM
  CALL cmfe_Equations_Initialise(CaMEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(CaMEquationsSet,CaMEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(CaMEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(CaMEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(CaMEquationsSet,Err)

  !Create the equations set equations for CaMCa
  CALL cmfe_Equations_Initialise(CaMCaEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(CaMCaEquationsSet,CaMCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(CaMCaEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(CaMCaEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(CaMCaEquationsSet,Err)

  !Create the equations set equations for ATP
  CALL cmfe_Equations_Initialise(ATPEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(ATPEquationsSet,ATPEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(ATPEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(ATPEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(ATPEquationsSet,Err)

  !Create the equations set equations for ATPCa
  CALL cmfe_Equations_Initialise(ATPCaEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(ATPCaEquationsSet,ATPCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(ATPCaEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(ATPCaEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(ATPCaEquationsSet,Err)

!____________________________________________________________________________________________________________

  WRITE(*,*) 'Create the problem'
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[cmfe_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & cmfe_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,cmfe_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],Problem,Err)
  !Set the problem to be a strang split reaction diffusion problem
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,cmfe_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,startT,endT,Tstep,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,outputfreq,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,cmfe_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !CALL cmfe_ControlLoopTimesSet(ControlLoop,0.0_CMISSDP,5.00_CMISSDP,0.01_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

!______________________________________________________________________________________________________________
  !Set up the problem solvers for Strang splitting - 
  !note, this example is contrived to have strang splitting, when it could be solved as a simple evaluation (as opposed to integration) of source and diffusion
  WRITE(*,*) 'Set up the problem solvers for Strang splitting'
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL cmfe_Solver_DynamicThetaSet(Solver,1.0_CMISSDP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  !set linear solver to be direct solver. Note, I found this stuff in fluidmechanics/darcy/dynamic/src example
  !CALL cmfe_SolverLinearTypeSet(LinearSolver,cmfe_SolverLinearDirectSolveType,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_Solvercmfe_Library,Err)
  !CALL cmfe_SolverLinearTypeSet(LinearSolver,cmfe_SolverLinearDirectSolveType,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_SolverMUMPSLibrary,Err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)


  !Third solver is another DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err) !set the third solver's integration time step
  CALL cmfe_Solver_OutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverTimingOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverSolverOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SOLVER_PROGRESS_OUTPUT,Err)

  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

!_______________________________________________________________________________________________________________
  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

!_______________________________________________________________________________________________________________
  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,cmfe_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquationsSparsityTypeSet(SolverEquations,cmfe_SolverEquationsSparseMatrices,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,cmfe_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set for Ca, F and FCa
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CaEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,FEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,FCaEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CaMEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CaMCaEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,ATPEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,ATPCaEquationsSet,EquationsSetIndex,Err)

  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)
!_________________________________________________________________________________________________________
  WRITE(*,*) 'Set up boundary conditions'  
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Set 0 conc. bc on nodes within mitos
  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO node=1,NUMBER_OF_NODES
      NODE_NUMBER = NodeNums(node,1)
      IF(NodeNums(node,2).EQ.MITO_REGION_MARKER) THEN
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CONDITION = cmfe_BOUNDARY_CONDITION_FIXED
          VALUE=0.0_CMISSDP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaField, &
            & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FField, &
            & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FCaField, &
            & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)


          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)


          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(dirchlet boundary condition)

        ENDIF
      ENDIF
    ENDDO
  ENDIF
  !set no flux bc if node is a mito boundary node.

  IF(WITH_MITO_ELEMENTS.EQ.1) THEN
    DO node=1,NUMBER_OF_MITOBDFACENODES
      NODE_NUMBER = MITOBDFaceNodes(node)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CONDITION = cmfe_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSDP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)!

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FCaField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        !make sure to unset mito bd node bcs to be free from dirchlet bcs set above.

        CONDITION = cmfe_BOUNDARY_CONDITION_FREE
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_Ca,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_F,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FCaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_FCa,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_CaM,Err) !(neumann boundary condition - no flux)

         CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_CaMCa,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_ATP,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
          & cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,init_ATPCa,Err) !(neumann boundary condition - no flux)



        !Set the initial conc. at these mito boundary nodes to be the cytosolic versions
        CALL cmfe_Field_ParameterSetUpdateNode(CaField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_Ca,Err)

        CALL cmfe_Field_ParameterSetUpdateNode(FField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_F,Err)

        CALL cmfe_Field_ParameterSetUpdateNode(FCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_FCa,Err)

        CALL cmfe_Field_ParameterSetUpdateNode(CaMField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaM,Err)


        CALL cmfe_Field_ParameterSetUpdateNode(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_CaMCa,Err)

        CALL cmfe_Field_ParameterSetUpdateNode(ATPField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_ATP,Err)

        CALL cmfe_Field_ParameterSetUpdateNode(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE, &
          & cmfe_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,init_ATPCa,Err)


      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateStart(CaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(FField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(FCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(FCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(CaMField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(CaMCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(ATPField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    CALL cmfe_Field_ParameterSetUpdateStart(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(ATPCaField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

    ENDIF

  !Set no flux on cell boundary
  DO node=1,NUMBER_OF_CELLBDNODES
    NODE_NUMBER = CELLBDNodes(node)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CONDITION = cmfe_BOUNDARY_CONDITION_FIXED
      VALUE=0.0_CMISSDP
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaField, &
        & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
        & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FField, &
        & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
        & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FCaField, &
        & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
        & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMField, &
      & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
      & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
      & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)


      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
      & cmfe_FIELD_DELUDELN_VARIABLE_TYPE,1,cmfe_NO_GLOBAL_DERIV, &
      & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

!__________________________________________________________________________________________________________
  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)
!__________________________________________________________________________________________________________
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"Cell_Solution","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"Cell_Solution","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF 
  !Finialise cmfe_
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP


END PROGRAM CARDIAC_ECC
