### Anura3D_2022 ###	
*############################################
*##### dimension                        #####
*############################################
$$DIMENSION
*set var dimension=GenData(dimension_ID,int)
*if(dimension==100)
2D-plane_strain
*end if
*if(dimension==101)
2D-axisymmetric
*end if
*if(dimension==102)
3D-cartesian
*end if
*if(dimension==103)
3D-cylindrical
*end if
*if((ndime==3) && (dimension==100))
*MessageBox **** INPUT ERROR: 3D element type is incompatible with chosen dimension (2D - plane strain). ****
*end if
*if((ndime==3) && (dimension==101))
*MessageBox **** INPUT ERROR: 3D element type is incompatible with chosen dimension (2D - axisymmetric). ****
*end if
*if((ndime==2) && (dimension==102))
*MessageBox **** INPUT ERROR: 2D element type is incompatible with chosen dimension (3D - Cartesian). ****
*end if
*if((ndime==2) && (dimension==103))
*MessageBox **** INPUT ERROR: 2D element type is incompatible with chosen dimension (3D - cylindrical). ****
*end if
*############################################
*##### elementtype                      #####
*############################################
$$ELEMENTTYPE
*loop elems
*set var I=elemstype
*set var J=nnode
*end elems
*if((I==2) && (J==3))
*set elems(triangle)
triangular_3-noded
*#*elseif((I==2) && (J==6))
*#*set elems(triangle)
*#triangular_6-noded
*#*elseif((I==3) && (J==4))
*#*set elems(quadrilateral)
*#quadrilateral_4-noded
*#*elseif((I==3) && (J==8))
*#*set elems(quadrilateral)
*#quadrilateral_8-noded
*#*elseif((I==4) && (J==4))
*#*set elems(tetrahedra)
*#tetrahedral_4-noded
*elseif((I==4) && (J==10))
*set elems(tetrahedra)
tetrahedral_old
*#*elseif((I==5) && (J==8))
*#*set elems(hexahedra)
*#hexahedral_8-noded
*#*elseif((I==5) && (J==20))
*#*set elems(hexahedra)
*#hexahedral_20-noded
*else
*MessageBox **** INPUT ERROR: Element type not properly defined. Only the following element types are supported: triangular 3-noded, tetrahedral 10-noded. ****
*end if
*################################################
*##### formulation (single or double point) #####
*################################################
$$FORMULATION
*if(ndime==2)
*set cond 2D_-_Single-point_formulation
*set var IsSinglePointFormulation2D=0
*loop elems *OnlyIncond
*set var IsSinglePointFormulation2D=1
*end elems
*if((IsSinglePointFormulation2D==1))
single-point
*end if
*set cond 2D_-_Double-point_formulation
*set var IsDoublePointFormulation2D=0
*loop elems *OnlyIncond
*set var IsDoublePointFormulation2D=1
*end elems
*if((IsDoublePointFormulation2D==1))
double-point
*end if
*if((IsDoublePointFormulation2D==1) && (IsSinglePointFormulation2D==1))
*MessageBox **** INPUT ERROR: The Material Point Specification can be for either 'single-point formulation' or 'double-point formulation' and not combined. ****
*end if
*if((IsDoublePointFormulation2D==0) && (IsSinglePointFormulation2D==0))
*MessageBox **** INPUT ERROR: No material points are assigned. Please assign material points for solid and/or liquid! Select <Anura3D><Material Point Specification> to assign the material points. Regenerate the mesh. ****
*end if
*end if
*if(ndime==3)
*set cond 3D_-_Single-point_formulation
*set var IsSinglePointFormulation3D=0
*loop elems *OnlyIncond
*set var IsSinglePointFormulation3D=1
*end elems
*if((IsSinglePointFormulation3D==1))
single-point
*end if
*set cond 3D_-_Double-point_formulation
*set var IsDoublePointFormulation3D=0
*loop elems *OnlyIncond
*set var IsDoublePointFormulation3D=1
*end elems
*if((IsDoublePointFormulation3D==1))
double-point
*end if
*if((IsDoublePointFormulation3D==1) && (IsSinglePointFormulation3D==1))
*MessageBox **** INPUT ERROR: The Material Point Specification can be for either 'single-point formulation' or 'double-point formulation' and not combined. ****
*end if
*if((IsDoublePointFormulation3D==0) && (IsSinglePointFormulation3D==0))
*MessageBox **** INPUT ERROR: No material points are assigned. Please assign material points for solid and/or liquid! Select <Anura3D><Material Point Specification> to assign the material points. Regenerate the mesh. ****
*end if
*end if
*############################################
*##### counters                         #####
*############################################
$$STARTCOUNTERS
*nelem *npoin
*############################################
*##### nodal coordinates                #####
*############################################
$$STARTNODES
*loop nodes
*if(ndime==3)
*NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real)
*elseif(ndime==2)
*NodesCoord(1,real) *NodesCoord(2,real)
*endif
*end nodes
*############################################
*##### element connectivities           #####
*############################################
$$STARTELEMCON
*loop elems
*elemsconec
*end elems
*############################################
*##### fixities                         #####
*############################################
$$START_FIXITY_SURFACE_SOLID
*set Cond Solid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ISS=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_LINE_SOLID
*set Cond Solid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ILS=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_POINT_SOLID
*set Cond Solid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var IPS=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_SURFACE_LIQUID
*set Cond Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ISL=I
*loop nodes*OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_LINE_LIQUID
*set Cond Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ILL=I
*loop nodes *OnlyInCond
*if(ndime ==3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_POINT_LIQUID
*set Cond Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_and_Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var IPL=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_SURFACE_GAS
*set Cond Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ISG=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_LINE_GAS
*set Cond Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var ILG=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
$$START_FIXITY_POINT_GAS
*set Cond Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*set var IPG=I
*loop nodes *OnlyInCond
*if(ndime == 3)
*NodesNum *cond(1) *cond(2) *cond(3)
*else
*NodesNum *cond(1) *cond(2)
*endif
*end nodes
*##### check if any fixities are assigned ###
*set var ITOT=operation(ISS+ILS+IPS+ISL+ILL+IPL+ISG+ILG+IPG)
*if((ITOT==0))
*MessageBox **** INPUT ERROR: No fixities are assigned. Please specify fixities for solid and/or liquid and/or gas! Select <Anura3D><Fixities> to assign the fixities. Regenerate the mesh. ****
*end if
*############################################
*##### remove fixities                  #####
*############################################
$$START_REMOVE_FIXITY_SURFACE_SOLID
*set Cond Remove_Solid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_LINE_SOLID
*set Cond Remove_Solid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_POINT_SOLID
*set Cond Remove_Solid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_SURFACE_LIQUID
*set Cond Remove_Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_LINE_LIQUID
*set Cond Remove_Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_POINT_LIQUID
*set Cond Remove_Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_and_Liquid_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_SURFACE_GAS
*set Cond Remove_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(surface) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_LINE_GAS
*set Cond Remove_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(line) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
$$START_REMOVE_FIXITY_POINT_GAS
*set Cond Remove_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*add Cond Remove_Solid_Liquid_and_Gas_Fixity_(point) *nodes *or(1,int) *or(2,int) *or(3,int)
*set var I=0
*loop nodes *OnlyInCond
*set var I=Operation(I+1)
*end nodes
*I
*loop nodes *OnlyInCond
*if(ndime==3)
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int)) *operation((-1)*cond(3,int))
*else
*NodesNum *operation((-1)*cond(1,int)) *operation((-1)*cond(2,int))
*endif
*end nodes
*############################################
*##### prescribed velocities 2D/3D      #####
*#####   (only written when specified)  ##### 
*############################################
*##### 2D case
*if(ndime==2)
*set Cond 2D_-_Material_point_velocity_(surface) *elems
*set var nPrescribedVelocity2D=0
*loop elems *OnlyInCond
*set var nPrescribedVelocity2D=Operation(nPrescribedVelocity2D+1)
*end elems
*if(nPrescribedVelocity2D>0)
$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE
*nPrescribedVelocity2D
*loop elems *OnlyInCond
*elemsnum *cond(1) *cond(3) *cond(2) *cond(4)
*end
*end if
*set Cond 2D_-_Nodal_velocity_(surface) *nodes
*set var nPrescribedVelocity2D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity2D=Operation(nPrescribedVelocity2D+1)
*end nodes
*if(nPrescribedVelocity2D>0)
$$PRESCRIBED_NODAL_VELOCITY_SURFACE
*nPrescribedVelocity2D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(2) *cond(4)
*end
*end if
*set Cond 2D_-_Nodal_velocity_(line) *nodes
*set var nPrescribedVelocity2D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity2D=Operation(nPrescribedVelocity2D+1)
*end nodes
*if(nPrescribedVelocity2D>0)
$$PRESCRIBED_NODAL_VELOCITY_LINE
*nPrescribedVelocity2D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(2) *cond(4)
*end
*end if
*set Cond 2D_-_Nodal_velocity_(point) *nodes
*set var nPrescribedVelocity2D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity2D=Operation(nPrescribedVelocity2D+1)
*end nodes
*if(nPrescribedVelocity2D>0)
$$PRESCRIBED_NODAL_VELOCITY_POINT
*nPrescribedVelocity2D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(2) *cond(4)
*end
*end if
*end if
*##### 3D case
*if(ndime==3)
*set Cond 3D_-_Material_point_velocity_(volume) *elems
*set var nPrescribedVelocity3D=0
*loop elems *OnlyInCond
*set var nPrescribedVelocity3D=Operation(nPrescribedVelocity3D+1)
*end elems
*if(nPrescribedVelocity3D>0)
$$PRESCRIBED_MATERIAL_POINT_VELOCITY_VOLUME
*nPrescribedVelocity3D
*loop elems *OnlyInCond
*elemsnum *cond(1) *cond(3) *cond(5) *cond(2) *cond(4) *cond(6)
*end
*end if
*set Cond 3D_-_Nodal_velocity_(volume) *nodes
*set var nPrescribedVelocity3D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity3D=Operation(nPrescribedVelocity3D+1)
*end nodes
*if(nPrescribedVelocity3D>0)
$$PRESCRIBED_NODAL_VELOCITY_VOLUME
*nPrescribedVelocity3D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(5) *cond(2) *cond(4) *cond(6)
*end
*end if
*set Cond 3D_-_Nodal_velocity_(surface) *nodes
*set var nPrescribedVelocity3D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity3D=Operation(nPrescribedVelocity3D+1)
*end nodes
*if(nPrescribedVelocity3D>0)
$$PRESCRIBED_NODAL_VELOCITY_SURFACE
*nPrescribedVelocity3D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(5) *cond(2) *cond(4) *cond(6)
*end
*end if
*set Cond 3D_-_Nodal_velocity_(line) *nodes
*set var nPrescribedVelocity3D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity3D=Operation(nPrescribedVelocity3D+1)
*end nodes
*if(nPrescribedVelocity3D>0)
$$PRESCRIBED_NODAL_VELOCITY_LINE
*nPrescribedVelocity3D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(5) *cond(2) *cond(4) *cond(6)
*end
*end if
*set Cond 3D_-_Nodal_velocity_(point) *nodes
*set var nPrescribedVelocity3D=0
*loop nodes *OnlyInCond
*set var nPrescribedVelocity3D=Operation(nPrescribedVelocity3D+1)
*end nodes
*if(nPrescribedVelocity3D>0)
$$PRESCRIBED_NODAL_VELOCITY_POINT
*nPrescribedVelocity3D
*loop nodes *OnlyInCond
*nodesnum *cond(1) *cond(3) *cond(5) *cond(2) *cond(4) *cond(6)
*end
*end if
*end if
*############################################
*##### initial conditions 2D/3D         #####
*#####   (only written when specified)  ##### 
*############################################
*##### Initial Velocity
*##### 2D case
*if(ndime==2)
*set var nInitialVelocity2D=0
*set Cond 2D_-_Initial_velocity_(on_material_points) *elems
*loop elems *OnlyInCond
*set var nInitialVelocity2D=Operation(nInitialVelocity2D+1)
*end elems
*if(nInitialVelocity2D>0)
$$INITIAL_VELOCITY_MATERIAL_POINT
*nInitialVelocity2D
*loop elems *OnlyInCond
*elemsnum *cond(1,real) *cond(2,real)
*end
*end if
*end if
*##### 3D case
*if(ndime==3)
*set var nInitialVelocity3D=0
*set Cond 3D_-_Initial_velocity_(on_material_points) *elems
*loop elems *OnlyInCond
*set var nInitialVelocity3D=Operation(nInitialVelocity3D+1)
*end elems
*if(nInitialVelocity3D>0)
$$INITIAL_VELOCITY_MATERIAL_POINT
*nInitialVelocity3D
*loop elems *OnlyInCond
*elemsnum *cond(1,real) *cond(2,real) *cond(3,real)
*end
*end if
*end if
*##### Initial Water Surface From File
*##### 2D case
*if(ndime==2)
*set cond 2D_-_Water_surface_with_0_pressure_from_file
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var nWaterSurfaceFromFileMaterials=0
*loop elems *OnlyInCond
*set var MaterialIndex=elemsmat
*if(MaterialIndex==tmp1)
*set var tmp1=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp2)
*set var tmp2=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp3)
*set var tmp3=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp4)
*set var tmp4=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*end elems
*if(nWaterSurfaceFromFileMaterials>0)
$$INITIAL_WATER_SURFACE_FROM_FILE
*nWaterSurfaceFromFileMaterials
*endif
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var nWaterSurfaceFromFileMaterials=0
*loop elems onlyincond
*set var MaterialIndex=elemsmat
*if(MaterialIndex==tmp1)
*elemsmat *cond(1) 
*set var tmp1=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp2)
*elemsmat *cond(1) 
*set var tmp2=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp3)
*elemsmat *cond(1) 
*set var tmp3=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp4)
*elemsmat *cond(1) 
*set var tmp4=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*end elems
*endif
*##### 3D case
*if(ndime==3)
*set cond 3D_-_Water_surface_with_0_pressure_from_file
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var nWaterSurfaceFromFileMaterials=0
*loop elems *OnlyInCond
*set var MaterialIndex=elemsmat
*if(MaterialIndex==tmp1)
*set var tmp1=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp2)
*set var tmp2=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp3)
*set var tmp3=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp4)
*set var tmp4=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*end elems
*if(nWaterSurfaceFromFileMaterials>0)
$$INITIAL_WATER_SURFACE_FROM_FILE
*nWaterSurfaceFromFileMaterials
*endif
*set var tmp1=1
*set var tmp2=2
*set var tmp3=3
*set var tmp4=4
*set var nWaterSurfaceFromFileMaterials=0
*loop elems onlyincond
*set var MaterialIndex=elemsmat
*if(MaterialIndex==tmp1)
*elemsmat *cond(1) 
*set var tmp1=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp2)
*elemsmat *cond(1) 
*set var tmp2=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp3)
*elemsmat *cond(1) 
*set var tmp3=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*if(MaterialIndex==tmp4)
*elemsmat *cond(1) 
*set var tmp4=-1
*set var nWaterSurfaceFromFileMaterials=operation(nWaterSurfaceFromFileMaterials+1)
*endif
*end elems
*endif
*############################################
*##### loading conditions 2D/3D         #####
*#####   (only written when specified)  #####
*############################################
*##### 3D traction load nodes
*if(ndime==3)
*set Cond 3D_-_Solid_Traction *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(21),"nodes")==0)
*if(strcasecmp(cond(22),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(22),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_NODES_SOLID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(21),"nodes")==0) && (strcasecmp(cond(22),"A")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real)+LocalAxesDef(3)*cond(5,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real)+LocalAxesDef(6)*cond(5,real))
*set var Z=Operation(LocalAxesDef(7)*cond(3,real)+LocalAxesDef(8)*cond(4,real)+LocalAxesDef(9)*cond(5,real))
*GlobalNodes *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(6,real)
*Set var Y0=cond(7,real)
*Set var Z0=cond(8,real)
*Set var tracX0=cond(9,real)
*Set var tracY0=cond(10,real)
*Set var tracZ0=cond(11,real)
*Set var ixX=cond(12,real)
*Set var iyX=cond(13,real)
*Set var izX=cond(14,real)
*Set var ixY=cond(15,real)
*Set var iyY=cond(16,real)
*Set var izY=cond(17,real)
*Set var ixZ=cond(18,real)
*Set var iyZ=cond(19,real)
*Set var izZ=cond(20,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0)+izX*(Z1-Z0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0)+izY*(Z1-Z0))
*Set var tracZ1=(tracZ0+ixZ*(X1-X0)+iyZ*(Y1-Y0)+izZ*(Z1-Z0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1+LocalAxesDef(3)*tracZ1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1+LocalAxesDef(6)*tracZ1)
*set var tracZ1local=(LocalAxesDef(7)*tracX1+LocalAxesDef(8)*tracY1+LocalAxesDef(9)*tracZ1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0)+izX*(Z2-Z0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0)+izY*(Z2-Z0))
*Set var tracZ2=(tracZ0+ixZ*(X2-X0)+iyZ*(Y2-Y0)+izZ*(Z2-Z0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2+LocalAxesDef(3)*tracZ2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2+LocalAxesDef(6)*tracZ2)
*set var tracZ2local=(LocalAxesDef(7)*tracX2+LocalAxesDef(8)*tracY2+LocalAxesDef(9)*tracZ2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0)+izX*(Z3-Z0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0)+izY*(Z3-Z0))
*Set var tracZ3=(tracZ0+ixZ*(X3-X0)+iyZ*(Y3-Y0)+izZ*(Z3-Z0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3+LocalAxesDef(3)*tracZ3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3+LocalAxesDef(6)*tracZ3)
*set var tracZ3local=(LocalAxesDef(7)*tracX3+LocalAxesDef(8)*tracY3+LocalAxesDef(9)*tracZ3)
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var tracX4=(tracX0+ixX*(X4-X0)+iyX*(Y4-Y0)+izX*(Z4-Z0))
*Set var tracY4=(tracY0+ixY*(X4-X0)+iyY*(Y4-Y0)+izY*(Z4-Z0))
*Set var tracZ4=(tracZ0+ixZ*(X4-X0)+iyZ*(Y4-Y0)+izZ*(Z4-Z0))
*set var tracX4local=(LocalAxesDef(1)*tracX4+LocalAxesDef(2)*tracY4+LocalAxesDef(3)*tracZ4)
*set var tracY4local=(LocalAxesDef(4)*tracX4+LocalAxesDef(5)*tracY4+LocalAxesDef(6)*tracZ4)
*set var tracZ4local=(LocalAxesDef(7)*tracX4+LocalAxesDef(8)*tracY4+LocalAxesDef(9)*tracZ4)
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var tracX5=(tracX0+ixX*(X5-X0)+iyX*(Y5-Y0)+izX*(Z5-Z0))
*Set var tracY5=(tracY0+ixY*(X5-X0)+iyY*(Y5-Y0)+izY*(Z5-Z0))
*Set var tracZ5=(tracZ0+ixZ*(X5-X0)+iyZ*(Y5-Y0)+izZ*(Z5-Z0))
*set var tracX5local=(LocalAxesDef(1)*tracX5+LocalAxesDef(2)*tracY5+LocalAxesDef(3)*tracZ5)
*set var tracY5local=(LocalAxesDef(4)*tracX5+LocalAxesDef(5)*tracY5+LocalAxesDef(6)*tracZ5)
*set var tracZ5local=(LocalAxesDef(7)*tracX5+LocalAxesDef(8)*tracY5+LocalAxesDef(9)*tracZ5)
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var tracX6=(tracX0+ixX*(X6-X0)+iyX*(Y6-Y0)+izX*(Z6-Z0))
*Set var tracY6=(tracY0+ixY*(X6-X0)+iyY*(Y6-Y0)+izY*(Z6-Z0))
*Set var tracZ6=(tracZ0+ixZ*(X6-X0)+iyZ*(Y6-Y0)+izZ*(Z6-Z0))
*set var tracX6local=(LocalAxesDef(1)*tracX6+LocalAxesDef(2)*tracY6+LocalAxesDef(3)*tracZ6)
*set var tracY6local=(LocalAxesDef(4)*tracX6+LocalAxesDef(5)*tracY6+LocalAxesDef(6)*tracZ6)
*set var tracZ6local=(LocalAxesDef(7)*tracX6+LocalAxesDef(8)*tracY6+LocalAxesDef(9)*tracZ6)
*end if
*end nodes
*GlobalNodes *tracX1local *tracY1local *tracZ1local *tracX2local *tracY2local *tracZ2local *tracX3local *tracY3local *tracZ3local *tracX4local *tracY4local *tracZ4local *tracX5local *tracY5local *tracZ5local *tracX6local *tracY6local *tracZ6local
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_NODES_SOLID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(21),"nodes")==0) && (strcasecmp(cond(22),"B")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real)+LocalAxesDef(3)*cond(5,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real)+LocalAxesDef(6)*cond(5,real))
*set var Z=Operation(LocalAxesDef(7)*cond(3,real)+LocalAxesDef(8)*cond(4,real)+LocalAxesDef(9)*cond(5,real))
*GlobalNodes *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(6,real)
*Set var Y0=cond(7,real)
*Set var Z0=cond(8,real)
*Set var tracX0=cond(9,real)
*Set var tracY0=cond(10,real)
*Set var tracZ0=cond(11,real)
*Set var ixX=cond(12,real)
*Set var iyX=cond(13,real)
*Set var izX=cond(14,real)
*Set var ixY=cond(15,real)
*Set var iyY=cond(16,real)
*Set var izY=cond(17,real)
*Set var ixZ=cond(18,real)
*Set var iyZ=cond(19,real)
*Set var izZ=cond(20,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0)+izX*(Z1-Z0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0)+izY*(Z1-Z0))
*Set var tracZ1=(tracZ0+ixZ*(X1-X0)+iyZ*(Y1-Y0)+izZ*(Z1-Z0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1+LocalAxesDef(3)*tracZ1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1+LocalAxesDef(6)*tracZ1)
*set var tracZ1local=(LocalAxesDef(7)*tracX1+LocalAxesDef(8)*tracY1+LocalAxesDef(9)*tracZ1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0)+izX*(Z2-Z0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0)+izY*(Z2-Z0))
*Set var tracZ2=(tracZ0+ixZ*(X2-X0)+iyZ*(Y2-Y0)+izZ*(Z2-Z0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2+LocalAxesDef(3)*tracZ2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2+LocalAxesDef(6)*tracZ2)
*set var tracZ2local=(LocalAxesDef(7)*tracX2+LocalAxesDef(8)*tracY2+LocalAxesDef(9)*tracZ2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0)+izX*(Z3-Z0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0)+izY*(Z3-Z0))
*Set var tracZ3=(tracZ0+ixZ*(X3-X0)+iyZ*(Y3-Y0)+izZ*(Z3-Z0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3+LocalAxesDef(3)*tracZ3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3+LocalAxesDef(6)*tracZ3)
*set var tracZ3local=(LocalAxesDef(7)*tracX3+LocalAxesDef(8)*tracY3+LocalAxesDef(9)*tracZ3)
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var tracX4=(tracX0+ixX*(X4-X0)+iyX*(Y4-Y0)+izX*(Z4-Z0))
*Set var tracY4=(tracY0+ixY*(X4-X0)+iyY*(Y4-Y0)+izY*(Z4-Z0))
*Set var tracZ4=(tracZ0+ixZ*(X4-X0)+iyZ*(Y4-Y0)+izZ*(Z4-Z0))
*set var tracX4local=(LocalAxesDef(1)*tracX4+LocalAxesDef(2)*tracY4+LocalAxesDef(3)*tracZ4)
*set var tracY4local=(LocalAxesDef(4)*tracX4+LocalAxesDef(5)*tracY4+LocalAxesDef(6)*tracZ4)
*set var tracZ4local=(LocalAxesDef(7)*tracX4+LocalAxesDef(8)*tracY4+LocalAxesDef(9)*tracZ4)
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var tracX5=(tracX0+ixX*(X5-X0)+iyX*(Y5-Y0)+izX*(Z5-Z0))
*Set var tracY5=(tracY0+ixY*(X5-X0)+iyY*(Y5-Y0)+izY*(Z5-Z0))
*Set var tracZ5=(tracZ0+ixZ*(X5-X0)+iyZ*(Y5-Y0)+izZ*(Z5-Z0))
*set var tracX5local=(LocalAxesDef(1)*tracX5+LocalAxesDef(2)*tracY5+LocalAxesDef(3)*tracZ5)
*set var tracY5local=(LocalAxesDef(4)*tracX5+LocalAxesDef(5)*tracY5+LocalAxesDef(6)*tracZ5)
*set var tracZ5local=(LocalAxesDef(7)*tracX5+LocalAxesDef(8)*tracY5+LocalAxesDef(9)*tracZ5)
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var tracX6=(tracX0+ixX*(X6-X0)+iyX*(Y6-Y0)+izX*(Z6-Z0))
*Set var tracY6=(tracY0+ixY*(X6-X0)+iyY*(Y6-Y0)+izY*(Z6-Z0))
*Set var tracZ6=(tracZ0+ixZ*(X6-X0)+iyZ*(Y6-Y0)+izZ*(Z6-Z0))
*set var tracX6local=(LocalAxesDef(1)*tracX6+LocalAxesDef(2)*tracY6+LocalAxesDef(3)*tracZ6)
*set var tracY6local=(LocalAxesDef(4)*tracX6+LocalAxesDef(5)*tracY6+LocalAxesDef(6)*tracZ6)
*set var tracZ6local=(LocalAxesDef(7)*tracX6+LocalAxesDef(8)*tracY6+LocalAxesDef(9)*tracZ6)
*end if
*end nodes
*GlobalNodes *tracX1local *tracY1local *tracZ1local *tracX2local *tracY2local *tracZ2local *tracX3local *tracY3local *tracZ3local *tracX4local *tracY4local *tracZ4local *tracX5local *tracY5local *tracZ5local *tracX6local *tracY6local *tracZ6local
*end if
*end if
*end elems
*end if
*##### 2D traction load nodes
*if(ndime==2)
*set Cond 2D_-_Solid_Traction *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(13),"nodes")==0)
*if(strcasecmp(cond(14),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(14),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_NODES_SOLID
*J
*endif
*loop elems *OnlyInCond
*if((strcasecmp(cond(13),"nodes")==0) && (strcasecmp(cond(14),"A")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real))
*GlobalNodes *X *Y *X *Y
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(5,real)
*Set var Y0=cond(6,real)
*Set var tracX0=cond(7,real)
*Set var tracY0=cond(8,real)
*Set var ixX=cond(9,real)
*Set var iyX=cond(10,real)
*Set var ixY=cond(11,real)
*Set var iyY=cond(12,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3)
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local *tracX3local *tracY3local
*elseif(nnode==3)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local
*end if
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_NODES_SOLID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(13),"nodes")==0) && (strcasecmp(cond(14),"B")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real))
*GlobalNodes *X *Y *X *Y
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(5,real)
*Set var Y0=cond(6,real)
*Set var tracX0=cond(7,real)
*Set var tracY0=cond(8,real)
*Set var ixX=cond(9,real)
*Set var iyX=cond(10,real)
*Set var ixY=cond(11,real)
*Set var iyY=cond(12,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3)
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local *tracX3local *tracY3local
*elseif(nnode==3)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local
*end if
*end if
*end if
*end elems
*end if
*##### 3D traction load material points 
*if(ndime==3)
*set Cond 3D_-_Solid_Traction *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(21),"material_points")==0)
*if(strcasecmp(cond(22),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(22),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_MATERIAL_POINTS_SOLID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(21),"material_points")==0) && (strcasecmp(cond(22),"A")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real)+LocalAxesDef(3)*cond(5,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real)+LocalAxesDef(6)*cond(5,real))
*set var Z=Operation(LocalAxesDef(7)*cond(3,real)+LocalAxesDef(8)*cond(4,real)+LocalAxesDef(9)*cond(5,real))
*GlobalNodes *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(6,real)
*Set var Y0=cond(7,real)
*Set var Z0=cond(8,real)
*Set var tracX0=cond(9,real)
*Set var tracY0=cond(10,real)
*Set var tracZ0=cond(11,real)
*Set var ixX=cond(12,real)
*Set var iyX=cond(13,real)
*Set var izX=cond(14,real)
*Set var ixY=cond(15,real)
*Set var iyY=cond(16,real)
*Set var izY=cond(17,real)
*Set var ixZ=cond(18,real)
*Set var iyZ=cond(19,real)
*Set var izZ=cond(20,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0)+izX*(Z1-Z0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0)+izY*(Z1-Z0))
*Set var tracZ1=(tracZ0+ixZ*(X1-X0)+iyZ*(Y1-Y0)+izZ*(Z1-Z0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1+LocalAxesDef(3)*tracZ1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1+LocalAxesDef(6)*tracZ1)
*set var tracZ1local=(LocalAxesDef(7)*tracX1+LocalAxesDef(8)*tracY1+LocalAxesDef(9)*tracZ1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0)+izX*(Z2-Z0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0)+izY*(Z2-Z0))
*Set var tracZ2=(tracZ0+ixZ*(X2-X0)+iyZ*(Y2-Y0)+izZ*(Z2-Z0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2+LocalAxesDef(3)*tracZ2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2+LocalAxesDef(6)*tracZ2)
*set var tracZ2local=(LocalAxesDef(7)*tracX2+LocalAxesDef(8)*tracY2+LocalAxesDef(9)*tracZ2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0)+izX*(Z3-Z0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0)+izY*(Z3-Z0))
*Set var tracZ3=(tracZ0+ixZ*(X3-X0)+iyZ*(Y3-Y0)+izZ*(Z3-Z0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3+LocalAxesDef(3)*tracZ3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3+LocalAxesDef(6)*tracZ3)
*set var tracZ3local=(LocalAxesDef(7)*tracX3+LocalAxesDef(8)*tracY3+LocalAxesDef(9)*tracZ3)
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var tracX4=(tracX0+ixX*(X4-X0)+iyX*(Y4-Y0)+izX*(Z4-Z0))
*Set var tracY4=(tracY0+ixY*(X4-X0)+iyY*(Y4-Y0)+izY*(Z4-Z0))
*Set var tracZ4=(tracZ0+ixZ*(X4-X0)+iyZ*(Y4-Y0)+izZ*(Z4-Z0))
*set var tracX4local=(LocalAxesDef(1)*tracX4+LocalAxesDef(2)*tracY4+LocalAxesDef(3)*tracZ4)
*set var tracY4local=(LocalAxesDef(4)*tracX4+LocalAxesDef(5)*tracY4+LocalAxesDef(6)*tracZ4)
*set var tracZ4local=(LocalAxesDef(7)*tracX4+LocalAxesDef(8)*tracY4+LocalAxesDef(9)*tracZ4)
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var tracX5=(tracX0+ixX*(X5-X0)+iyX*(Y5-Y0)+izX*(Z5-Z0))
*Set var tracY5=(tracY0+ixY*(X5-X0)+iyY*(Y5-Y0)+izY*(Z5-Z0))
*Set var tracZ5=(tracZ0+ixZ*(X5-X0)+iyZ*(Y5-Y0)+izZ*(Z5-Z0))
*set var tracX5local=(LocalAxesDef(1)*tracX5+LocalAxesDef(2)*tracY5+LocalAxesDef(3)*tracZ5)
*set var tracY5local=(LocalAxesDef(4)*tracX5+LocalAxesDef(5)*tracY5+LocalAxesDef(6)*tracZ5)
*set var tracZ5local=(LocalAxesDef(7)*tracX5+LocalAxesDef(8)*tracY5+LocalAxesDef(9)*tracZ5)
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var tracX6=(tracX0+ixX*(X6-X0)+iyX*(Y6-Y0)+izX*(Z6-Z0))
*Set var tracY6=(tracY0+ixY*(X6-X0)+iyY*(Y6-Y0)+izY*(Z6-Z0))
*Set var tracZ6=(tracZ0+ixZ*(X6-X0)+iyZ*(Y6-Y0)+izZ*(Z6-Z0))
*set var tracX6local=(LocalAxesDef(1)*tracX6+LocalAxesDef(2)*tracY6+LocalAxesDef(3)*tracZ6)
*set var tracY6local=(LocalAxesDef(4)*tracX6+LocalAxesDef(5)*tracY6+LocalAxesDef(6)*tracZ6)
*set var tracZ6local=(LocalAxesDef(7)*tracX6+LocalAxesDef(8)*tracY6+LocalAxesDef(9)*tracZ6)
*end if
*end nodes
*GlobalNodes *tracX1local *tracY1local *tracZ1local *tracX2local *tracY2local *tracZ2local *tracX3local *tracY3local *tracZ3local *tracX4local *tracY4local *tracZ4local *tracX5local *tracY5local *tracZ5local *tracX6local *tracY6local *tracZ6local
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_MATERIAL_POINTS_SOLID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(21),"material_points")==0) && (strcasecmp(cond(22),"B")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real)+LocalAxesDef(3)*cond(5,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real)+LocalAxesDef(6)*cond(5,real))
*set var Z=Operation(LocalAxesDef(7)*cond(3,real)+LocalAxesDef(8)*cond(4,real)+LocalAxesDef(9)*cond(5,real))
*GlobalNodes *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z *X *Y *Z
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(6,real)
*Set var Y0=cond(7,real)
*Set var Z0=cond(8,real)
*Set var tracX0=cond(9,real)
*Set var tracY0=cond(10,real)
*Set var tracZ0=cond(11,real)
*Set var ixX=cond(12,real)
*Set var iyX=cond(13,real)
*Set var izX=cond(14,real)
*Set var ixY=cond(15,real)
*Set var iyY=cond(16,real)
*Set var izY=cond(17,real)
*Set var ixZ=cond(18,real)
*Set var iyZ=cond(19,real)
*Set var izZ=cond(20,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0)+izX*(Z1-Z0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0)+izY*(Z1-Z0))
*Set var tracZ1=(tracZ0+ixZ*(X1-X0)+iyZ*(Y1-Y0)+izZ*(Z1-Z0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1+LocalAxesDef(3)*tracZ1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1+LocalAxesDef(6)*tracZ1)
*set var tracZ1local=(LocalAxesDef(7)*tracX1+LocalAxesDef(8)*tracY1+LocalAxesDef(9)*tracZ1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0)+izX*(Z2-Z0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0)+izY*(Z2-Z0))
*Set var tracZ2=(tracZ0+ixZ*(X2-X0)+iyZ*(Y2-Y0)+izZ*(Z2-Z0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2+LocalAxesDef(3)*tracZ2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2+LocalAxesDef(6)*tracZ2)
*set var tracZ2local=(LocalAxesDef(7)*tracX2+LocalAxesDef(8)*tracY2+LocalAxesDef(9)*tracZ2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0)+izX*(Z3-Z0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0)+izY*(Z3-Z0))
*Set var tracZ3=(tracZ0+ixZ*(X3-X0)+iyZ*(Y3-Y0)+izZ*(Z3-Z0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3+LocalAxesDef(3)*tracZ3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3+LocalAxesDef(6)*tracZ3)
*set var tracZ3local=(LocalAxesDef(7)*tracX3+LocalAxesDef(8)*tracY3+LocalAxesDef(9)*tracZ3)
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var tracX4=(tracX0+ixX*(X4-X0)+iyX*(Y4-Y0)+izX*(Z4-Z0))
*Set var tracY4=(tracY0+ixY*(X4-X0)+iyY*(Y4-Y0)+izY*(Z4-Z0))
*Set var tracZ4=(tracZ0+ixZ*(X4-X0)+iyZ*(Y4-Y0)+izZ*(Z4-Z0))
*set var tracX4local=(LocalAxesDef(1)*tracX4+LocalAxesDef(2)*tracY4+LocalAxesDef(3)*tracZ4)
*set var tracY4local=(LocalAxesDef(4)*tracX4+LocalAxesDef(5)*tracY4+LocalAxesDef(6)*tracZ4)
*set var tracZ4local=(LocalAxesDef(7)*tracX4+LocalAxesDef(8)*tracY4+LocalAxesDef(9)*tracZ4)
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var tracX5=(tracX0+ixX*(X5-X0)+iyX*(Y5-Y0)+izX*(Z5-Z0))
*Set var tracY5=(tracY0+ixY*(X5-X0)+iyY*(Y5-Y0)+izY*(Z5-Z0))
*Set var tracZ5=(tracZ0+ixZ*(X5-X0)+iyZ*(Y5-Y0)+izZ*(Z5-Z0))
*set var tracX5local=(LocalAxesDef(1)*tracX5+LocalAxesDef(2)*tracY5+LocalAxesDef(3)*tracZ5)
*set var tracY5local=(LocalAxesDef(4)*tracX5+LocalAxesDef(5)*tracY5+LocalAxesDef(6)*tracZ5)
*set var tracZ5local=(LocalAxesDef(7)*tracX5+LocalAxesDef(8)*tracY5+LocalAxesDef(9)*tracZ5)
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var tracX6=(tracX0+ixX*(X6-X0)+iyX*(Y6-Y0)+izX*(Z6-Z0))
*Set var tracY6=(tracY0+ixY*(X6-X0)+iyY*(Y6-Y0)+izY*(Z6-Z0))
*Set var tracZ6=(tracZ0+ixZ*(X6-X0)+iyZ*(Y6-Y0)+izZ*(Z6-Z0))
*set var tracX6local=(LocalAxesDef(1)*tracX6+LocalAxesDef(2)*tracY6+LocalAxesDef(3)*tracZ6)
*set var tracY6local=(LocalAxesDef(4)*tracX6+LocalAxesDef(5)*tracY6+LocalAxesDef(6)*tracZ6)
*set var tracZ6local=(LocalAxesDef(7)*tracX6+LocalAxesDef(8)*tracY6+LocalAxesDef(9)*tracZ6)
*end if
*end nodes
*GlobalNodes *tracX1local *tracY1local *tracZ1local *tracX2local *tracY2local *tracZ2local *tracX3local *tracY3local *tracZ3local *tracX4local *tracY4local *tracZ4local *tracX5local *tracY5local *tracZ5local *tracX6local *tracY6local *tracZ6local
*end if
*end if
*end elems
*end if
*##### 2D traction load material points
*if(ndime==2)
*set Cond 2D_-_Solid_Traction *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(13),"material_points")==0)
*if(strcasecmp(cond(14),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(14),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_MATERIAL_POINTS_SOLID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(13),"material_points")==0) && (strcasecmp(cond(14),"A")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real))
*GlobalNodes *X *Y *X *Y
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(5,real)
*Set var Y0=cond(6,real)
*Set var tracX0=cond(7,real)
*Set var tracY0=cond(8,real)
*Set var ixX=cond(9,real)
*Set var iyX=cond(10,real)
*Set var ixY=cond(11,real)
*Set var iyY=cond(12,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3)
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local *tracX3local *tracY3local
*elseif(nnode== 3)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local
*end if
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_MATERIAL_POINTS_SOLID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(13),"material_points")==0) && (strcasecmp(cond(14),"B")==0))
*if(strcasecmp(cond(2),"uniform_distribution")==0)
*set var X=Operation(LocalAxesDef(1)*cond(3,real)+LocalAxesDef(2)*cond(4,real))
*set var Y=Operation(LocalAxesDef(4)*cond(3,real)+LocalAxesDef(5)*cond(4,real))
*GlobalNodes *X *Y *X *Y
*end if
*if(strcasecmp(cond(2),"linear_distribution")==0)
*Set var X0=cond(5,real)
*Set var Y0=cond(6,real)
*Set var tracX0=cond(7,real)
*Set var tracY0=cond(8,real)
*Set var ixX=cond(9,real)
*Set var iyX=cond(10,real)
*Set var ixY=cond(11,real)
*Set var iyY=cond(12,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var tracX1=(tracX0+ixX*(X1-X0)+iyX*(Y1-Y0))
*Set var tracY1=(tracY0+ixY*(X1-X0)+iyY*(Y1-Y0))
*set var tracX1local=(LocalAxesDef(1)*tracX1+LocalAxesDef(2)*tracY1)
*set var tracY1local=(LocalAxesDef(4)*tracX1+LocalAxesDef(5)*tracY1)
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var tracX2=(tracX0+ixX*(X2-X0)+iyX*(Y2-Y0))
*Set var tracY2=(tracY0+ixY*(X2-X0)+iyY*(Y2-Y0))
*set var tracX2local=(LocalAxesDef(1)*tracX2+LocalAxesDef(2)*tracY2)
*set var tracY2local=(LocalAxesDef(4)*tracX2+LocalAxesDef(5)*tracY2)
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var tracX3=(tracX0+ixX*(X3-X0)+iyX*(Y3-Y0))
*Set var tracY3=(tracY0+ixY*(X3-X0)+iyY*(Y3-Y0))
*set var tracX3local=(LocalAxesDef(1)*tracX3+LocalAxesDef(2)*tracY3)
*set var tracY3local=(LocalAxesDef(4)*tracX3+LocalAxesDef(5)*tracY3)
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local *tracX3local *tracY3local
*elseif(nnode== 3)
*GlobalNodes *tracX1local *tracY1local *tracX2local *tracY2local
*end if
*end if
*end if
*end elems
*end if
*##### 3D pressure nodes
*if(ndime==3)
*set Cond 3D_-_Liquid_Pressure *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(10),"nodes")==0) 
*if(strcasecmp(cond(11),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(11),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_NODES_LIQUID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(10),"nodes")==0) && (strcasecmp(cond(11),"A")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var Z0=cond(6,real)
*Set var ix=cond(7,real)
*Set var iy=cond(8,real)
*Set var iz=cond(9,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0)+iz*(Z1-Z0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0)+iz*(Z2-Z0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0)+iz*(Z3-Z0))
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var pp4=(pp0+ix*(X4-X0)+iy*(Y4-Y0)+iz*(Z4-Z0))
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var pp5=(pp0+ix*(X5-X0)+iy*(Y5-Y0)+iz*(Z5-Z0))
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var pp6=(pp0+ix*(X6-X0)+iy*(Y6-Y0)+iz*(Z6-Z0))
*end if
*end nodes
*GlobalNodes *pp1 *pp1 *pp1 *pp2 *pp2 *pp2 *pp3 *pp3 *pp3 *pp4 *pp4 *pp4 *pp5 *pp5 *pp5 *pp6 *pp6 *pp6
*end if
*end if
*end elems
*if(K>0)
$START_LOAD_ON_NODES_LIQUID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(10),"nodes")==0) && (strcasecmp(cond(11),"B")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var Z0=cond(6,real)
*Set var ix=cond(7,real)
*Set var iy=cond(8,real)
*Set var iz=cond(9,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0)+iz*(Z1-Z0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0)+iz*(Z2-Z0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0)+iz*(Z3-Z0))
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var pp4=(pp0+ix*(X4-X0)+iy*(Y4-Y0)+iz*(Z4-Z0))
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var pp5=(pp0+ix*(X5-X0)+iy*(Y5-Y0)+iz*(Z5-Z0))
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var pp6=(pp0+ix*(X6-X0)+iy*(Y6-Y0)+iz*(Z6-Z0))
*end if
*end nodes
*GlobalNodes *pp1 *pp1 *pp1 *pp2 *pp2 *pp2 *pp3 *pp3 *pp3 *pp4 *pp4 *pp4 *pp5 *pp5 *pp5 *pp6 *pp6 *pp6
*end if
*end if
*end elems
*end if
*##### 2D pressure nodes
*if(ndime==2)
*set Cond 2D_-_Liquid_Pressure *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(8),"nodes")==0)
*if(strcasecmp(cond(9),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(9),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_NODES_LIQUID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(8),"nodes")==0) && (strcasecmp(cond(9),"A")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var ix=cond(6,real)
*Set var iy=cond(7,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0))
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 *pp3 *pp3
*elseif(nnode==3)
*GlobalNodes *pp1 *pp1 *pp2 *pp2
*end if
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_NODES_LIQUID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(8),"nodes")==0) && (strcasecmp(cond(9),"B")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var ix=cond(6,real)
*Set var iy=cond(7,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0))
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 *pp3 *pp3
*elseif(nnode==3)
*GlobalNodes *pp1 *pp1 *pp2 *pp2
*end if
*end if
*end if
*end elems
*end if
*##### 3D pressure material points
*if(ndime==3)
*set Cond 3D_-_Liquid_Pressure *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(10),"material_points")==0)
*if((strcasecmp(cond(11),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(11),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_MATERIAL_POINTS_LIQUID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(10),"material_points")==0) && (strcasecmp(cond(11),"A")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var Z0=cond(6,real)
*Set var ix=cond(7,real)
*Set var iy=cond(8,real)
*Set var iz=cond(9,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0)+iz*(Z1-Z0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0)+iz*(Z2-Z0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0)+iz*(Z3-Z0))
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var pp4=(pp0+ix*(X4-X0)+iy*(Y4-Y0)+iz*(Z4-Z0))
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var pp5=(pp0+ix*(X5-X0)+iy*(Y5-Y0)+iz*(Z5-Z0))
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var pp6=(pp0+ix*(X6-X0)+iy*(Y6-Y0)+iz*(Z6-Z0))
*end if
*end nodes
*GlobalNodes *pp1 *pp1 *pp1 *pp2 *pp2 *pp2 *pp3 *pp3 *pp3 *pp4 *pp4 *pp4 *pp5 *pp5 *pp5 *pp6 *pp6 *pp6
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_MATERIAL_POINTS_LIQUID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(10),"material_points")==0) && (strcasecmp(cond(11),"B")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var Z0=cond(6,real)
*Set var ix=cond(7,real)
*Set var iy=cond(8,real)
*Set var iz=cond(9,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*Set var Node4=GlobalNodes(4)
*Set var Node5=GlobalNodes(5)
*Set var Node6=GlobalNodes(6)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var Z1=NodesCoord(3,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0)+iz*(Z1-Z0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var Z2=NodesCoord(3,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0)+iz*(Z2-Z0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var Z3=NodesCoord(3,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0)+iz*(Z3-Z0))
*end if
*if((NodesNum==Node4))
*Set var X4=NodesCoord(1,real)
*Set var Y4=NodesCoord(2,real)
*Set var Z4=NodesCoord(3,real)
*Set var pp4=(pp0+ix*(X4-X0)+iy*(Y4-Y0)+iz*(Z4-Z0))
*end if
*if((NodesNum==Node5))
*Set var X5=NodesCoord(1,real)
*Set var Y5=NodesCoord(2,real)
*Set var Z5=NodesCoord(3,real)
*Set var pp5=(pp0+ix*(X5-X0)+iy*(Y5-Y0)+iz*(Z5-Z0))
*end if
*if((NodesNum==Node6))
*Set var X6=NodesCoord(1,real)
*Set var Y6=NodesCoord(2,real)
*Set var Z6=NodesCoord(3,real)
*Set var pp6=(pp0+ix*(X6-X0)+iy*(Y6-Y0)+iz*(Z6-Z0))
*end if
*end nodes
*GlobalNodes *pp1 *pp1 *pp1 *pp2 *pp2 *pp2 *pp3 *pp3 *pp3 *pp4 *pp4 *pp4 *pp5 *pp5 *pp5 *pp6 *pp6 *pp6
*end if
*end if
*end elems
*end if
*##### 2D pressure material points
*if(ndime==2)
*set Cond 2D_-_Liquid_Pressure *elems *CanRepeat
*set var J=0
*set var K=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(8),"material_points")==0)
*if(strcasecmp(cond(9),"A")==0)
*set var J=Operation(J+1)
*elseif(strcasecmp(cond(9),"B")==0)
*set var K=Operation(K+1)
*end if
*end if
*end elems
*if(J>0)
$$START_LOAD_ON_MATERIAL_POINTS_LIQUID
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(8),"material_points")==0) && (strcasecmp(cond(9),"A")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var ix=cond(6,real)
*Set var iy=cond(7,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0))
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 *pp3 *pp3
*elseif(nnode==3)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 
*end if
*end if
*end if
*end elems
*if(K>0)
$$START_LOAD_ON_MATERIAL_POINTS_LIQUID_B
*K
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(8),"material_points")==0) && (strcasecmp(cond(9),"B")==0))
*if(strcasecmp(cond(1),"uniform_distribution")==0)
*GlobalNodes *cond(2,real) *cond(2,real) *cond(2,real) *cond(2,real)
*end if
*if(strcasecmp(cond(1),"linear_distribution")==0)
*Set var pp0=cond(3,real)
*Set var X0=cond(4,real)
*Set var Y0=cond(5,real)
*Set var ix=cond(6,real)
*Set var iy=cond(7,real)
*Set var Node1=GlobalNodes(1)
*Set var Node2=GlobalNodes(2)
*Set var Node3=GlobalNodes(3)
*loop nodes
*if((NodesNum==Node1))
*Set var X1=NodesCoord(1,real)
*Set var Y1=NodesCoord(2,real)
*Set var pp1=(pp0+ix*(X1-X0)+iy*(Y1-Y0))
*end if
*if((NodesNum==Node2))
*Set var X2=NodesCoord(1,real)
*Set var Y2=NodesCoord(2,real)
*Set var pp2=(pp0+ix*(X2-X0)+iy*(Y2-Y0))
*end if
*if((NodesNum==Node3))
*Set var X3=NodesCoord(1,real)
*Set var Y3=NodesCoord(2,real)
*Set var pp3=(pp0+ix*(X3-X0)+iy*(Y3-Y0))
*end if
*end nodes
*if(nnode==6)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 *pp3 *pp3
*elseif(nnode==3)
*GlobalNodes *pp1 *pp1 *pp2 *pp2 
*end if
*end if
*end if
*end elems
*end if
*############################################
*##### Soil surface definition 2D   #####
*#####   (only written when specified)  #####
*############################################
*##### 2D phreatic surface
*if(ndime==2)
*set Cond 2D_-_Soil_Surface_Line *elems *CanRepeat
*set var J=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(1),"nodes")==0)
*set var J=Operation(J+1)
*end if
*end elems
*if(J>0)
$$START_SOIL_SURFACE_NODES
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(1),"nodes")==0))
*GlobalNodes 
*end if
*end elems
*end if
*############################################
*##### phreatic surface definition 2D   #####
*#####   (only written when specified)  #####
*############################################
*##### 2D phreatic surface
*if(ndime==2)
*set Cond 2D_-_Phreatic_Line *elems *CanRepeat
*set var J=0
*loop elems *OnlyInCond
*if(strcasecmp(cond(1),"nodes")==0)
*set var J=Operation(J+1)
*end if
*end elems
*if(J>0)
$$START_PHREATIC_SURFACE_NODES
*J
*end if
*loop elems *OnlyInCond
*if((strcasecmp(cond(1),"nodes")==0))
*GlobalNodes 
*end if
*end elems
*end if
*############################################
*##### hydraulic boundary conditions    #####
*############################################
*if(ndime==3)
*##### 3D hydraulic head
*set Cond Hydraulic_Head
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*set var zmin=0
*set var zmax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_y-_and_z-min_[3D_only]")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*set var zmin = nodescoord(3)
*end if
*if(strcasecmp(cond(2),"x-_y-_and_z-max_[3D_only]")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*set var zmax = nodescoord(3)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_HYDRAULIC_HEAD_AREA
*xmin *xmax
*ymin *ymax
*zmin *zmax
*end if
*end if
*if(ndime==2)
*##### 2D hydraulic head
*set Cond Hydraulic_Head
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_and_y-min")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*end if
*if(strcasecmp(cond(2),"x-_and_y-max")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_HYDRAULIC_HEAD_AREA
*xmin *xmax
*ymin *ymax
*end if
*end if
*if(ndime==3)
*##### 3D seepage face
*set Cond Seepage_Face
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*set var zmin=0
*set var zmax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_y-_and_z-min_[3D_only]")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*set var zmin = nodescoord(3)
*end if
*if(strcasecmp(cond(2),"x-_y-_and_z-max_[3D_only]")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*set var zmax = nodescoord(3)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_SEEPAGE_AREA
*xmin *xmax
*ymin *ymax
*zmin *zmax
*end if
*end if
*if(ndime==2)
*##### 2D seepage face
*set Cond Seepage_Face
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_and_y-min")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*end if
*if(strcasecmp(cond(2),"x-_and_y-max")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_SEEPAGE_AREA
*xmin *xmax
*ymin *ymax
*end if
*end if
*if(ndime==3)
*##### 3D infiltration
*set Cond Infiltration
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*set var zmin=0
*set var zmax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_y-_and_z-min_[3D_only]")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*set var zmin = nodescoord(3)
*end if
*if(strcasecmp(cond(2),"x-_y-_and_z-max_[3D_only]")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*set var zmax = nodescoord(3)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_INFILTRATION_AREA
*xmin *xmax
*ymin *ymax
*zmin *zmax
*end if
*set var J=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(J==1)
$$INFILTRATION_RATE
*cond(4) *cond(5) *cond(6) 
*endif
*end nodes
*end if
*if(ndime==2)
*##### 2D infiltration
*set Cond Infiltration
*set var J=0
*set var xmin=0
*set var xmax=0
*set var ymin=0
*set var ymax=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(strcasecmp(cond(2),"x-_and_y-min")==0)
*set var xmin = nodescoord(1) 
*set var ymin = nodescoord(2)
*end if
*if(strcasecmp(cond(2),"x-_and_y-max")==0)
*set var xmax = nodescoord(1) 
*set var ymax = nodescoord(2)
*end if
*end nodes
*if(J>0)
$$BOUNDARY_INFILTRATION_AREA
*xmin *xmax
*ymin *ymax
*end if
*set var J=0
*loop nodes *OnlyInCond
*set var J=Operation(J+1)
*if(J==1)
$$INFILTRATION_RATE
*cond(4) *cond(5)
*endif
*end nodes
*end if
*############################################
*##### contact properties               #####
*############################################
*if(ndime==3)
*# *##### 3D boundary contact 
*# $$START_CONTACT_SURFACE
*# *set Cond 3D_-_Boundary_contact *elems *CanRepeat
*# *set var J=0
*# *loop elems *OnlyInCond
*# *set var J=Operation(J+1)
*# *end elems
*# *J
*# *loop elems *OnlyInCond
*# *if((cond(1,int)==1))
*# *GlobalNodes *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0 "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==2))
*# *GlobalNodes *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==3))
*# *GlobalNodes *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==4))
*# *GlobalNodes *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "*cond(12)" *cond(13,real) *cond(14,real) 
*# *end if
*# *end elems
*##### 3D body contact 
$$START_CONTACT_VOLUME
*set Cond 3D_-_Body_contact *elems *CanRepeat
*set var J=0
*loop elems *OnlyInCond
*set var J=Operation(J+1)
*end elems
*J
*loop elems *OnlyInCond
*if((cond(1,int)==1))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0 "NAN" 0.0 0.0
*end if
*if((cond(1,int)==2))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0
*end if
*if((cond(1,int)==3))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "NAN" 0.0 0.0
*end if
*if((cond(1,int)==4))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "*cond(12)" *cond(13,real) *cond(14,real) 
*end if
*end elems
*end if
*if(ndime==2)
*##### 2D body contact 
*set Cond 2D_-_Body_contact *elems *CanRepeat
*set var nBodyContact2D=0
*loop elems *OnlyInCond
*set var nBodyContact2D=Operation(nBodyContact2D+1)
*end elems
*if(nBodyContact2D>0)
$$START_BODY_CONTACT_2D
*nBodyContact2D
*end if
*loop elems *OnlyInCond
*if((cond(1,int)==1))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0 "NAN" 0.0 0.0
*end if
*if((cond(1,int)==2))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0
*end if
*if((cond(1,int)==3))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "NAN" 0.0 0.0
*end if
*if((cond(1,int)==4))
*elemsnum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "*cond(12)" *cond(13,real) *cond(14,real) 
*end if
*end elems
*# *##### 2D boundary contact
*# *set Cond 2D_-_Boundary_contact *nodes
*# *set var nBoundaryContact2D=0
*# *loop nodes *OnlyInCond
*# *set var nBoundaryContact2D=Operation(nBoundaryContact2D+1)
*# *end nodes
*# *if(nBoundaryContact2D>0)
*# $$START_BOUNDARY_CONTACT_2D
*# *nBoundaryContact2D
*# *end if
*# *loop nodes *OnlyInCond
*# *if((cond(1,int)==1))
*# *NodesNum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0 "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==2))
*# *NodesNum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "NAN" 0.0 0.0 "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==3))
*# *NodesNum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "NAN" 0.0 0.0
*# *end if
*# *if((cond(1,int)==4))
*# *NodesNum *cond(1,int) "*cond(3)" *cond(4,real) *cond(5,real) "*cond(6)" *cond(7,real) *cond(8,real) "*cond(9)" *cond(10,real) *cond(11,real) "*cond(12)" *cond(13,real) *cond(14,real) 
*# *end if
*# *end nodes
*end if
*############################################
*##### excavation properties            #####
*############################################
*set cond 3D_-_Solid_Excavation *elems 
*set var nExcavation=0
*loop elems *OnlyInCond
*set var nExcavation=operation(nExcavation+1)
*end elems
*if(nExcavation>0)
$$START_EXCAVATION_SOLID
*end if
*loop elems *OnlyInCond
*cond(2) *elemsnum
*end elems
*set cond 2D_-_Solid_Excavation *elems 
*set var nExcavation2D=0
*loop elems *OnlyInCond
*set var nExcavation2D=operation(nExcavation2D+1)
*end elems
*if(nExcavation2D>0)
$$START_EXCAVATION_SOLID
*end if
*loop elems *OnlyInCond
*cond(2) *elemsnum
*end elems
*############################################
*##### output reaction forces 2D/3D     #####
*#####   (only written when specified)  ##### 
*############################################
*##### 2D case
*if(ndime==2)
*set cond 2D_-_Reaction_Forces *elems 
*set var nReactionForces=0
*loop elems *OnlyInCond
*set var nReactionForces=operation(nReactionForces+1)
*end elems
*if(nReactionForces>0)
$$START_OUTPUT_REACTION_FORCES
*nReactionForces
*loop elems *OnlyInCond
"*cond(1)" *elemsnum *GlobalNodes 
*end elems
*end if
*end if
*##### 3D case
*if(ndime==3)
*set cond 3D_-_Reaction_Forces *elems
*set var nReactionForces=0
*loop elems *OnlyInCond
*set var nReactionForces=operation(nReactionForces+1)
*end elems
*if(nReactionForces>0)
$$START_OUTPUT_REACTION_FORCES
*nReactionForces
*loop elems *OnlyInCond
"*cond(1)" *elemsnum *GlobalNodes
*end elems
*end if
*end if
*############################################
*##### absorbing boundaries 2D/3D       #####
*#####   (only written when specified)  ##### 
*############################################
*set var IsAbsorbingBoundary=0
*##### surface conditions
*set Cond Solid_Absorbing_Boundary_(surface)
*set var nABSolidSurface=0
*loop nodes *OnlyInCond
*set var nABSolidSurface=operation(nABSolidSurface+1)
*end nodes
*if(nABSolidSurface>0)
*set var IsAbsorbingBoundary=1
$$ABSORBING_BOUNDARY_SURFACE_SOLID
*nABSolidSurface
*end if
*loop nodes *OnlyInCond
*if(ndime==2)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6)
*end if
*if(ndime==3)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) *cond(9)
*end if
*end nodes
*##### line conditions
*set Cond Solid_Absorbing_Boundary_(line)
*set var nABSolidLine=0
*loop nodes *OnlyInCond
*set var nABSolidLine=operation(nABSolidLine+1)
*end nodes
*if(nABSolidLine>0)
*set var IsAbsorbingBoundary=1
$$ABSORBING_BOUNDARY_LINE_SOLID
*nABSolidLine
*end if
*loop nodes *OnlyInCond
*if(ndime==2)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6)
*end if
*if(ndime==3)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) *cond(9)
*end if
*end nodes
*##### point conditions
*set Cond Solid_Absorbing_Boundary_(point)
*set var nABSolidPoint=0
*loop nodes *OnlyInCond
*set var nABSolidPoint=operation(nABSolidPoint+1)
*end nodes
*if(nABSolidPoint>0)
*set var IsAbsorbingBoundary=1
$$ABSORBING_BOUNDARY_POINT_SOLID
*nABSolidPoint
*loop nodes *OnlyInCond
*if(ndime==2)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6)
*end if
*if(ndime==3)
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) *cond(9)
*end if
*end nodes
*end if
*##### reference material
*if(IsAbsorbingBoundary==1)
*if(ndime==2)
*set Cond 2D_-_Absorbing_Boundary_Reference_Material
*end if
*if(ndime==3)
*set Cond 3D_-_Absorbing_Boundary_Reference_Material
*end if
*set var I=0
$$ABSORBING_BOUNDARY_REFERENCE_MATERIAL_INDEX
*loop elems *OnlyInCond
*set var I=operation(I+1)
*if(I==1)
*ElemsMat 
*end if
*end elems
*if((I==0))
*MessageBox **** INPUT ERROR: No reference material for absorbing boundaries is assigned. Select <Anura3D><Absorbing Boundaries> to assign a reference material. Regenerate the mesh. ****
*end if
*end if
*############################################
*##### moving mesh                      #####
*############################################
*set var MM=0
*set Cond Extending_Mesh
*set var nMMextending=0
*loop nodes *OnlyInCond
*set var nMMextending=operation(nMMextending+1)
*end nodes
*if(nMMextending>0)
*set var MM=1
$$EXTENDING_MESH_CORNER_NODES
*end if
*loop nodes *OnlyInCond
*NodesNum
*end nodes
*set Cond Compressing_Mesh
*set var nMMcompressing=0
*loop nodes *OnlyInCond
*set var nMMcompressing=operation(nMMcompressing+1)
*end nodes
*if(nMMcompressing>0)
*set var MM=1
$$COMPRESSING_MESH_CORNER_NODES
*end if
*loop nodes *OnlyInCond
*NodesNum
*end nodes
*set Cond Moving_Mesh
*set var nMMmoving=0
*loop nodes *OnlyInCond
*set var nMMmoving=operation(nMMmoving+1)
*end nodes
*if(nMMmoving>0)
*set var MM=1
$$MOVING_MESH_CORNER_NODES
*end if
*loop nodes *OnlyInCond
*NodesNum
*end nodes
*if(nMMmoving>0)
$$MOVING_MESH_DIRECTION
*set var nNodes=0
*loop nodes *OnlyInCond
*set var nNodes=operation(nNodes+1)
*if(nNodes==1)
*cond(1)
*end if
*end nodes
*end if
*##### reference material 
*if(MM==1)
*set var I=0
$$MOVING_MESH_REFERENCE_MATERIAL_INDEX
*##### 3D case
*set Cond 3D_-_Moving_Mesh_Reference_Material
*loop elems *OnlyInCond
*set var I=operation(I+1)
*if(I==1)
*ElemsMat 
*end if
*end elems
*##### 3D case
*set Cond 2D_-_Moving_Mesh_Reference_Material
*loop elems *OnlyInCond
*set var I=operation(I+1)
*if(I==1)
*ElemsMat 
*end if
*end elems
*if(I==0)
*MessageBox **** INPUT ERROR: No reference material for moving mesh is assigned. Select <Anura3D><Moving Mesh> to assign a reference material. Regenerate the mesh. ****
*end if
*end if
*##############################################
*##### material properties 2D=3D          #####
*#####   (always specified and written)   #####
*##############################################
$$NUMBER_OF_MATERIALS
*nmats
*loop materials
$$MATERIAL_INDEX
*matnum()
$$MATERIAL_NAME
*MatProp(0)
*########## material type: undefined ##########
*if(MatProp(materialtype_id,int)==200)
$$MATERIAL_TYPE
undefined
*endif
*########## material type: dry_material ##########
*if((MatProp(materialtype_id,int)==201))
$$MATERIAL_TYPE
dry_material
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: saturated_material_drained ##########
*if((MatProp(materialtype_id,int)==202))
$$MATERIAL_TYPE
saturated_material_drained
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: saturated_material_undrained_effective ##########
*if((MatProp(materialtype_id,int)==203))
$$MATERIAL_TYPE
saturated_material_undrained_effective
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: saturated_material_undrained_total ##########
*if((MatProp(materialtype_id,int)==204))
$$MATERIAL_TYPE
saturated_material_undrained_total
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: saturated_material_coupled ##########
*if((MatProp(materialtype_id,int)==205))
$$MATERIAL_TYPE
saturated_material_coupled
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$INTRINSIC_PERMEABILITY_LIQUID
*MatProp(intrinsic_permeability_liquid_)
$$BULK_MODULUS_LIQUID
*MatProp(bulk_modulus_liquid_)
$$DYNAMIC_VISCOSITY_LIQUID
*MatProp(dynamic_viscosity_liquid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: liquid ##########
*if(MatProp(materialtype_id,int)==206)
$$MATERIAL_TYPE
liquid
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$BULK_MODULUS_LIQUID
*MatProp(bulk_modulus_liquid_)
$$DYNAMIC_VISCOSITY_LIQUID
*MatProp(dynamic_viscosity_liquid_)
$$LIQUID_CAVITATION
*MatProp(liquid_cavitation_)
$$APPLY_DETECT_LIQUID_SURFACE
*if(strcasecmp(MatProp(detect_free_surface_),"yes_")==0)
1 *MatProp(free_surface_factor_)
*else
0 *MatProp(free_surface_factor_)
*endif
*endif
*########## material type: unsaturated_material_2phase_suction ##########
*if((MatProp(materialtype_id,int)==207))
$$MATERIAL_TYPE
unsaturated_material_2phase_suction
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$INTRINSIC_PERMEABILITY_LIQUID
*MatProp(intrinsic_permeability_liquid_)
$$BULK_MODULUS_LIQUID
*MatProp(bulk_modulus_liquid_)
$$DYNAMIC_VISCOSITY_LIQUID
*MatProp(dynamic_viscosity_liquid_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material type: unsaturated_material_3phase_coupled ##########
*if((MatProp(materialtype_id,int)==208))
$$MATERIAL_TYPE
unsaturated_material_3phase_coupled
$$POROSITY_SOLID
*MatProp(initial_porosity_)
$$DENSITY_SOLID 
*MatProp(density_solid_)
$$DENSITY_LIQUID
*MatProp(density_liquid_)
$$DENSITY_GAS
*MatProp(density_gas_)
$$INTRINSIC_PERMEABILITY_LIQUID
*MatProp(intrinsic_permeability_liquid_)
$$INTRINSIC_PERMEABILITY_GAS
*MatProp(intrinsic_permeability_gas_)
$$BULK_MODULUS_LIQUID
*MatProp(bulk_modulus_liquid_)
$$BULK_MODULUS_GAS
*MatProp(bulk_modulus_gas_)
$$DYNAMIC_VISCOSITY_LIQUID
*MatProp(dynamic_viscosity_liquid_)
$$DYNAMIC_VISCOSITY_GAS
*MatProp(dynamic_viscosity_gas_)
$$SWELLING_INDEX
*MatProp(elastic_swelling_index_)
$$K0_VALUE_SOLID
*MatProp(K0-value_)
*endif
*########## material model: dry_material / saturated_material_drained / saturated_material_coupled / unsaturated_material_2phase_suction / unsaturated_material_3phase_coupled ##########
*if((MatProp(materialtype_id,int)==201)||(MatProp(materialtype_id,int)==202)||(MatProp(materialtype_id,int)==205)||(MatProp(materialtype_id,int)==207)||(MatProp(materialtype_id,int)==208))
*if(MatProp(model_id,int)==100)
$$MATERIAL_MODEL_SOLID
undefined
*endif
*if(MatProp(model_id,int)==101)
$$MATERIAL_MODEL_SOLID
linear_elasticity
$$YOUNG_MODULUS
*MatProp(effective_Young_modulus_)
$$POISSON_RATIO
*MatProp(effective_Poisson_ratio_)
*endif
*if(MatProp(model_id,int)==102)
$$MATERIAL_MODEL_SOLID
mohr_coulomb
$$YOUNG_MODULUS
*MatProp(effective_Young_modulus_)
$$POISSON_RATIO
*MatProp(effective_Poisson_ratio_)
$$FRICTION_ANGLE
*MatProp(effective_friction_angle_)
$$COHESION
*MatProp(effective_cohesion_)
$$DILATANCY_ANGLE
*MatProp(effective_dilatancy_angle_)
$$TENSILE_STRENGTH
*MatProp(tensile_strength_)
*endif
*if(MatProp(model_id,int)==105)
$$MATERIAL_MODEL_SOLID
rigid_body
$$constraint_XDISPLACEMENT
*MatProp(constraint_x-displacement_[0_or_1]:)
$$constraint_YDISPLACEMENT
*MatProp(constraint_y-displacement_[0_or_1]:)
$$constraint_ZDISPLACEMENT
*MatProp(constraint_z-displacement_[0_or_1]:)
*endif
*if(MatProp(model_id,int)==199)
$$MATERIAL_MODEL_SOLID
external_soil_model
$$MATERIAL_MODEL_DLL
*MatProp(External_Material_Model_DLL_)
$$UMAT_DIMENSION
*MatProp(size_of_stress_tensor_inside_UMAT_)
$$MATERIAL_PARAMETER_SOLID_01
*MatProp(material_parameter_solid_01_)
$$MATERIAL_PARAMETER_SOLID_02
*MatProp(material_parameter_solid_02_)
$$MATERIAL_PARAMETER_SOLID_03
*MatProp(material_parameter_solid_03_)
$$MATERIAL_PARAMETER_SOLID_04
*MatProp(material_parameter_solid_04_)
$$MATERIAL_PARAMETER_SOLID_05
*MatProp(material_parameter_solid_05_)
$$MATERIAL_PARAMETER_SOLID_06
*MatProp(material_parameter_solid_06_)
$$MATERIAL_PARAMETER_SOLID_07
*MatProp(material_parameter_solid_07_)
$$MATERIAL_PARAMETER_SOLID_08
*MatProp(material_parameter_solid_08_)
$$MATERIAL_PARAMETER_SOLID_09
*MatProp(material_parameter_solid_09_)
$$MATERIAL_PARAMETER_SOLID_10
*MatProp(material_parameter_solid_10_)
$$MATERIAL_PARAMETER_SOLID_11
*MatProp(material_parameter_solid_11_)
$$MATERIAL_PARAMETER_SOLID_12
*MatProp(material_parameter_solid_12_)
$$MATERIAL_PARAMETER_SOLID_13
*MatProp(material_parameter_solid_13_)
$$MATERIAL_PARAMETER_SOLID_14
*MatProp(material_parameter_solid_14_)
$$MATERIAL_PARAMETER_SOLID_15
*MatProp(material_parameter_solid_15_)
$$MATERIAL_PARAMETER_SOLID_16
*MatProp(material_parameter_solid_16_)
$$MATERIAL_PARAMETER_SOLID_17
*MatProp(material_parameter_solid_17_)
$$MATERIAL_PARAMETER_SOLID_18
*MatProp(material_parameter_solid_18_)
$$MATERIAL_PARAMETER_SOLID_19
*MatProp(material_parameter_solid_19_)
$$MATERIAL_PARAMETER_SOLID_20
*MatProp(material_parameter_solid_20_)
$$MATERIAL_PARAMETER_SOLID_21
*MatProp(material_parameter_solid_21_)
$$MATERIAL_PARAMETER_SOLID_22
*MatProp(material_parameter_solid_22_)
$$MATERIAL_PARAMETER_SOLID_23
*MatProp(material_parameter_solid_23_)
$$MATERIAL_PARAMETER_SOLID_24
*MatProp(material_parameter_solid_24_)
$$MATERIAL_PARAMETER_SOLID_25
*MatProp(material_parameter_solid_25_)
$$MATERIAL_PARAMETER_SOLID_26
*MatProp(material_parameter_solid_26_)
$$MATERIAL_PARAMETER_SOLID_27
*MatProp(material_parameter_solid_27_)
$$MATERIAL_PARAMETER_SOLID_28
*MatProp(material_parameter_solid_28_)
$$MATERIAL_PARAMETER_SOLID_29
*MatProp(material_parameter_solid_29_)
$$MATERIAL_PARAMETER_SOLID_30
*MatProp(material_parameter_solid_30_)
$$MATERIAL_PARAMETER_SOLID_31
*MatProp(material_parameter_solid_31_)
$$MATERIAL_PARAMETER_SOLID_32
*MatProp(material_parameter_solid_32_)
$$MATERIAL_PARAMETER_SOLID_33
*MatProp(material_parameter_solid_33_)
$$MATERIAL_PARAMETER_SOLID_34
*MatProp(material_parameter_solid_34_)
$$MATERIAL_PARAMETER_SOLID_35
*MatProp(material_parameter_solid_35_)
$$MATERIAL_PARAMETER_SOLID_36
*MatProp(material_parameter_solid_36_)
$$MATERIAL_PARAMETER_SOLID_37
*MatProp(material_parameter_solid_37_)
$$MATERIAL_PARAMETER_SOLID_38
*MatProp(material_parameter_solid_38_)
$$MATERIAL_PARAMETER_SOLID_39
*MatProp(material_parameter_solid_39_)
$$MATERIAL_PARAMETER_SOLID_40
*MatProp(material_parameter_solid_40_)
$$MATERIAL_PARAMETER_SOLID_41
*MatProp(material_parameter_solid_41_)
$$MATERIAL_PARAMETER_SOLID_42
*MatProp(material_parameter_solid_42_)
$$MATERIAL_PARAMETER_SOLID_43
*MatProp(material_parameter_solid_43_)
$$MATERIAL_PARAMETER_SOLID_44
*MatProp(material_parameter_solid_44_)
$$MATERIAL_PARAMETER_SOLID_45
*MatProp(material_parameter_solid_45_)
$$MATERIAL_PARAMETER_SOLID_46
*MatProp(material_parameter_solid_46_)
$$MATERIAL_PARAMETER_SOLID_47
*MatProp(material_parameter_solid_47_)
$$MATERIAL_PARAMETER_SOLID_48
*MatProp(material_parameter_solid_48_)
$$MATERIAL_PARAMETER_SOLID_49
*MatProp(material_parameter_solid_49_)
$$MATERIAL_PARAMETER_SOLID_50
*MatProp(material_parameter_solid_50_)
$$INITIAL_STATE_VARIABLE_SOLID_01
*MatProp(initial_state_variable_solid_01_)
$$INITIAL_STATE_VARIABLE_SOLID_02
*MatProp(initial_state_variable_solid_02_)
$$INITIAL_STATE_VARIABLE_SOLID_03
*MatProp(initial_state_variable_solid_03_)
$$INITIAL_STATE_VARIABLE_SOLID_04
*MatProp(initial_state_variable_solid_04_)
$$INITIAL_STATE_VARIABLE_SOLID_05
*MatProp(initial_state_variable_solid_05_)
$$INITIAL_STATE_VARIABLE_SOLID_06
*MatProp(initial_state_variable_solid_06_)
$$INITIAL_STATE_VARIABLE_SOLID_07
*MatProp(initial_state_variable_solid_07_)
$$INITIAL_STATE_VARIABLE_SOLID_08
*MatProp(initial_state_variable_solid_08_)
$$INITIAL_STATE_VARIABLE_SOLID_09
*MatProp(initial_state_variable_solid_09_)
$$INITIAL_STATE_VARIABLE_SOLID_10
*MatProp(initial_state_variable_solid_10_)
$$INITIAL_STATE_VARIABLE_SOLID_11
*MatProp(initial_state_variable_solid_11_)
$$INITIAL_STATE_VARIABLE_SOLID_12
*MatProp(initial_state_variable_solid_12_)
$$INITIAL_STATE_VARIABLE_SOLID_13
*MatProp(initial_state_variable_solid_13_)
$$INITIAL_STATE_VARIABLE_SOLID_14
*MatProp(initial_state_variable_solid_14_)
$$INITIAL_STATE_VARIABLE_SOLID_15
*MatProp(initial_state_variable_solid_15_)
$$INITIAL_STATE_VARIABLE_SOLID_16
*MatProp(initial_state_variable_solid_16_)
$$INITIAL_STATE_VARIABLE_SOLID_17
*MatProp(initial_state_variable_solid_17_)
$$INITIAL_STATE_VARIABLE_SOLID_18
*MatProp(initial_state_variable_solid_18_)
$$INITIAL_STATE_VARIABLE_SOLID_19
*MatProp(initial_state_variable_solid_19_)
$$INITIAL_STATE_VARIABLE_SOLID_20
*MatProp(initial_state_variable_solid_20_)
$$INITIAL_STATE_VARIABLE_SOLID_21
*MatProp(initial_state_variable_solid_21_)
$$INITIAL_STATE_VARIABLE_SOLID_22
*MatProp(initial_state_variable_solid_22_)
$$INITIAL_STATE_VARIABLE_SOLID_23
*MatProp(initial_state_variable_solid_23_)
$$INITIAL_STATE_VARIABLE_SOLID_24
*MatProp(initial_state_variable_solid_24_)
$$INITIAL_STATE_VARIABLE_SOLID_25
*MatProp(initial_state_variable_solid_25_)
$$INITIAL_STATE_VARIABLE_SOLID_26
*MatProp(initial_state_variable_solid_26_)
$$INITIAL_STATE_VARIABLE_SOLID_27
*MatProp(initial_state_variable_solid_27_)
$$INITIAL_STATE_VARIABLE_SOLID_28
*MatProp(initial_state_variable_solid_28_)
$$INITIAL_STATE_VARIABLE_SOLID_29
*MatProp(initial_state_variable_solid_29_)
$$INITIAL_STATE_VARIABLE_SOLID_30
*MatProp(initial_state_variable_solid_30_)
$$INITIAL_STATE_VARIABLE_SOLID_31
*MatProp(initial_state_variable_solid_31_)
$$INITIAL_STATE_VARIABLE_SOLID_32
*MatProp(initial_state_variable_solid_32_)
$$INITIAL_STATE_VARIABLE_SOLID_33
*MatProp(initial_state_variable_solid_33_)
$$INITIAL_STATE_VARIABLE_SOLID_34
*MatProp(initial_state_variable_solid_34_)
$$INITIAL_STATE_VARIABLE_SOLID_35
*MatProp(initial_state_variable_solid_35_)
$$INITIAL_STATE_VARIABLE_SOLID_36
*MatProp(initial_state_variable_solid_36_)
$$INITIAL_STATE_VARIABLE_SOLID_37
*MatProp(initial_state_variable_solid_37_)
$$INITIAL_STATE_VARIABLE_SOLID_38
*MatProp(initial_state_variable_solid_38_)
$$INITIAL_STATE_VARIABLE_SOLID_39
*MatProp(initial_state_variable_solid_39_)
$$INITIAL_STATE_VARIABLE_SOLID_40
*MatProp(initial_state_variable_solid_40_)
$$INITIAL_STATE_VARIABLE_SOLID_41
*MatProp(initial_state_variable_solid_41_)
$$INITIAL_STATE_VARIABLE_SOLID_42
*MatProp(initial_state_variable_solid_42_)
$$INITIAL_STATE_VARIABLE_SOLID_43
*MatProp(initial_state_variable_solid_43_)
$$INITIAL_STATE_VARIABLE_SOLID_44
*MatProp(initial_state_variable_solid_44_)
$$INITIAL_STATE_VARIABLE_SOLID_45
*MatProp(initial_state_variable_solid_45_)
$$INITIAL_STATE_VARIABLE_SOLID_46
*MatProp(initial_state_variable_solid_46_)
$$INITIAL_STATE_VARIABLE_SOLID_47
*MatProp(initial_state_variable_solid_47_)
$$INITIAL_STATE_VARIABLE_SOLID_48
*MatProp(initial_state_variable_solid_48_)
$$INITIAL_STATE_VARIABLE_SOLID_49
*MatProp(initial_state_variable_solid_49_)
$$INITIAL_STATE_VARIABLE_SOLID_50
*MatProp(initial_state_variable_solid_50_)
*endif
*endif 
*########## material model: saturated_material_undrained_total ##########
*if(MatProp(materialtype_id,int)==204)
*if(MatProp(total_id,int)==400)
$$MATERIAL_MODEL_SOLID
undefined
*endif
*if(MatProp(total_id,int)==401)
$$MATERIAL_MODEL_SOLID
linear_elasticity
$$YOUNG_MODULUS
*MatProp(undrained_Young_modulus_)
$$POISSON_RATIO
*MatProp(undrained_Poisson_ratio_)
*endif
*if(MatProp(total_id,int)==402)
$$MATERIAL_MODEL_SOLID
mohr_coulomb
$$YOUNG_MODULUS
*MatProp(undrained_Young_modulus_)
$$POISSON_RATIO
*MatProp(undrained_Poisson_ratio_)
$$FRICTION_ANGLE
*MatProp(friction_angle_)
$$COHESION
*MatProp(cohesion_)
$$DILATANCY_ANGLE
*MatProp(dilatancy_angle_)
$$TENSILE_STRENGTH
*MatProp(tensile_strength_)
*endif
*if(MatProp(total_id,int)==499)
$$MATERIAL_MODEL_SOLID
external_soil_model
$$MATERIAL_MODEL_DLL
*MatProp(External_Material_Model_DLL_)
$$UMAT_DIMENSION
*MatProp(size_of_stress_tensor_inside_UMAT_)
$$MATERIAL_PARAMETER_SOLID_01
*MatProp(material_parameter_solid_01_)
$$MATERIAL_PARAMETER_SOLID_02
*MatProp(material_parameter_solid_02_)
$$MATERIAL_PARAMETER_SOLID_03
*MatProp(material_parameter_solid_03_)
$$MATERIAL_PARAMETER_SOLID_04
*MatProp(material_parameter_solid_04_)
$$MATERIAL_PARAMETER_SOLID_05
*MatProp(material_parameter_solid_05_)
$$MATERIAL_PARAMETER_SOLID_06
*MatProp(material_parameter_solid_06_)
$$MATERIAL_PARAMETER_SOLID_07
*MatProp(material_parameter_solid_07_)
$$MATERIAL_PARAMETER_SOLID_08
*MatProp(material_parameter_solid_08_)
$$MATERIAL_PARAMETER_SOLID_09
*MatProp(material_parameter_solid_09_)
$$MATERIAL_PARAMETER_SOLID_10
*MatProp(material_parameter_solid_10_)
$$MATERIAL_PARAMETER_SOLID_11
*MatProp(material_parameter_solid_11_)
$$MATERIAL_PARAMETER_SOLID_12
*MatProp(material_parameter_solid_12_)
$$MATERIAL_PARAMETER_SOLID_13
*MatProp(material_parameter_solid_13_)
$$MATERIAL_PARAMETER_SOLID_14
*MatProp(material_parameter_solid_14_)
$$MATERIAL_PARAMETER_SOLID_15
*MatProp(material_parameter_solid_15_)
$$MATERIAL_PARAMETER_SOLID_16
*MatProp(material_parameter_solid_16_)
$$MATERIAL_PARAMETER_SOLID_17
*MatProp(material_parameter_solid_17_)
$$MATERIAL_PARAMETER_SOLID_18
*MatProp(material_parameter_solid_18_)
$$MATERIAL_PARAMETER_SOLID_19
*MatProp(material_parameter_solid_19_)
$$MATERIAL_PARAMETER_SOLID_20
*MatProp(material_parameter_solid_20_)
$$MATERIAL_PARAMETER_SOLID_21
*MatProp(material_parameter_solid_21_)
$$MATERIAL_PARAMETER_SOLID_22
*MatProp(material_parameter_solid_22_)
$$MATERIAL_PARAMETER_SOLID_23
*MatProp(material_parameter_solid_23_)
$$MATERIAL_PARAMETER_SOLID_24
*MatProp(material_parameter_solid_24_)
$$MATERIAL_PARAMETER_SOLID_25
*MatProp(material_parameter_solid_25_)
$$MATERIAL_PARAMETER_SOLID_26
*MatProp(material_parameter_solid_26_)
$$MATERIAL_PARAMETER_SOLID_27
*MatProp(material_parameter_solid_27_)
$$MATERIAL_PARAMETER_SOLID_28
*MatProp(material_parameter_solid_28_)
$$MATERIAL_PARAMETER_SOLID_29
*MatProp(material_parameter_solid_29_)
$$MATERIAL_PARAMETER_SOLID_30
*MatProp(material_parameter_solid_30_)
$$MATERIAL_PARAMETER_SOLID_31
*MatProp(material_parameter_solid_31_)
$$MATERIAL_PARAMETER_SOLID_32
*MatProp(material_parameter_solid_32_)
$$MATERIAL_PARAMETER_SOLID_33
*MatProp(material_parameter_solid_33_)
$$MATERIAL_PARAMETER_SOLID_34
*MatProp(material_parameter_solid_34_)
$$MATERIAL_PARAMETER_SOLID_35
*MatProp(material_parameter_solid_35_)
$$MATERIAL_PARAMETER_SOLID_36
*MatProp(material_parameter_solid_36_)
$$MATERIAL_PARAMETER_SOLID_37
*MatProp(material_parameter_solid_37_)
$$MATERIAL_PARAMETER_SOLID_38
*MatProp(material_parameter_solid_38_)
$$MATERIAL_PARAMETER_SOLID_39
*MatProp(material_parameter_solid_39_)
$$MATERIAL_PARAMETER_SOLID_40
*MatProp(material_parameter_solid_40_)
$$MATERIAL_PARAMETER_SOLID_41
*MatProp(material_parameter_solid_41_)
$$MATERIAL_PARAMETER_SOLID_42
*MatProp(material_parameter_solid_42_)
$$MATERIAL_PARAMETER_SOLID_43
*MatProp(material_parameter_solid_43_)
$$MATERIAL_PARAMETER_SOLID_44
*MatProp(material_parameter_solid_44_)
$$MATERIAL_PARAMETER_SOLID_45
*MatProp(material_parameter_solid_45_)
$$MATERIAL_PARAMETER_SOLID_46
*MatProp(material_parameter_solid_46_)
$$MATERIAL_PARAMETER_SOLID_47
*MatProp(material_parameter_solid_47_)
$$MATERIAL_PARAMETER_SOLID_48
*MatProp(material_parameter_solid_48_)
$$MATERIAL_PARAMETER_SOLID_49
*MatProp(material_parameter_solid_49_)
$$MATERIAL_PARAMETER_SOLID_50
*MatProp(material_parameter_solid_50_)
$$INITIAL_STATE_VARIABLE_SOLID_01
*MatProp(initial_state_variable_solid_01_)
$$INITIAL_STATE_VARIABLE_SOLID_02
*MatProp(initial_state_variable_solid_02_)
$$INITIAL_STATE_VARIABLE_SOLID_03
*MatProp(initial_state_variable_solid_03_)
$$INITIAL_STATE_VARIABLE_SOLID_04
*MatProp(initial_state_variable_solid_04_)
$$INITIAL_STATE_VARIABLE_SOLID_05
*MatProp(initial_state_variable_solid_05_)
$$INITIAL_STATE_VARIABLE_SOLID_06
*MatProp(initial_state_variable_solid_06_)
$$INITIAL_STATE_VARIABLE_SOLID_07
*MatProp(initial_state_variable_solid_07_)
$$INITIAL_STATE_VARIABLE_SOLID_08
*MatProp(initial_state_variable_solid_08_)
$$INITIAL_STATE_VARIABLE_SOLID_09
*MatProp(initial_state_variable_solid_09_)
$$INITIAL_STATE_VARIABLE_SOLID_10
*MatProp(initial_state_variable_solid_10_)
$$INITIAL_STATE_VARIABLE_SOLID_11
*MatProp(initial_state_variable_solid_11_)
$$INITIAL_STATE_VARIABLE_SOLID_12
*MatProp(initial_state_variable_solid_12_)
$$INITIAL_STATE_VARIABLE_SOLID_13
*MatProp(initial_state_variable_solid_13_)
$$INITIAL_STATE_VARIABLE_SOLID_14
*MatProp(initial_state_variable_solid_14_)
$$INITIAL_STATE_VARIABLE_SOLID_15
*MatProp(initial_state_variable_solid_15_)
$$INITIAL_STATE_VARIABLE_SOLID_16
*MatProp(initial_state_variable_solid_16_)
$$INITIAL_STATE_VARIABLE_SOLID_17
*MatProp(initial_state_variable_solid_17_)
$$INITIAL_STATE_VARIABLE_SOLID_18
*MatProp(initial_state_variable_solid_18_)
$$INITIAL_STATE_VARIABLE_SOLID_19
*MatProp(initial_state_variable_solid_19_)
$$INITIAL_STATE_VARIABLE_SOLID_20
*MatProp(initial_state_variable_solid_20_)
$$INITIAL_STATE_VARIABLE_SOLID_21
*MatProp(initial_state_variable_solid_21_)
$$INITIAL_STATE_VARIABLE_SOLID_22
*MatProp(initial_state_variable_solid_22_)
$$INITIAL_STATE_VARIABLE_SOLID_23
*MatProp(initial_state_variable_solid_23_)
$$INITIAL_STATE_VARIABLE_SOLID_24
*MatProp(initial_state_variable_solid_24_)
$$INITIAL_STATE_VARIABLE_SOLID_25
*MatProp(initial_state_variable_solid_25_)
$$INITIAL_STATE_VARIABLE_SOLID_26
*MatProp(initial_state_variable_solid_26_)
$$INITIAL_STATE_VARIABLE_SOLID_27
*MatProp(initial_state_variable_solid_27_)
$$INITIAL_STATE_VARIABLE_SOLID_28
*MatProp(initial_state_variable_solid_28_)
$$INITIAL_STATE_VARIABLE_SOLID_29
*MatProp(initial_state_variable_solid_29_)
$$INITIAL_STATE_VARIABLE_SOLID_30
*MatProp(initial_state_variable_solid_30_)
$$INITIAL_STATE_VARIABLE_SOLID_31
*MatProp(initial_state_variable_solid_31_)
$$INITIAL_STATE_VARIABLE_SOLID_32
*MatProp(initial_state_variable_solid_32_)
$$INITIAL_STATE_VARIABLE_SOLID_33
*MatProp(initial_state_variable_solid_33_)
$$INITIAL_STATE_VARIABLE_SOLID_34
*MatProp(initial_state_variable_solid_34_)
$$INITIAL_STATE_VARIABLE_SOLID_35
*MatProp(initial_state_variable_solid_35_)
$$INITIAL_STATE_VARIABLE_SOLID_36
*MatProp(initial_state_variable_solid_36_)
$$INITIAL_STATE_VARIABLE_SOLID_37
*MatProp(initial_state_variable_solid_37_)
$$INITIAL_STATE_VARIABLE_SOLID_38
*MatProp(initial_state_variable_solid_38_)
$$INITIAL_STATE_VARIABLE_SOLID_39
*MatProp(initial_state_variable_solid_39_)
$$INITIAL_STATE_VARIABLE_SOLID_40
*MatProp(initial_state_variable_solid_40_)
$$INITIAL_STATE_VARIABLE_SOLID_41
*MatProp(initial_state_variable_solid_41_)
$$INITIAL_STATE_VARIABLE_SOLID_42
*MatProp(initial_state_variable_solid_42_)
$$INITIAL_STATE_VARIABLE_SOLID_43
*MatProp(initial_state_variable_solid_43_)
$$INITIAL_STATE_VARIABLE_SOLID_44
*MatProp(initial_state_variable_solid_44_)
$$INITIAL_STATE_VARIABLE_SOLID_45
*MatProp(initial_state_variable_solid_45_)
$$INITIAL_STATE_VARIABLE_SOLID_46
*MatProp(initial_state_variable_solid_46_)
$$INITIAL_STATE_VARIABLE_SOLID_47
*MatProp(initial_state_variable_solid_47_)
$$INITIAL_STATE_VARIABLE_SOLID_48
*MatProp(initial_state_variable_solid_48_)
$$INITIAL_STATE_VARIABLE_SOLID_49
*MatProp(initial_state_variable_solid_49_)
$$INITIAL_STATE_VARIABLE_SOLID_50
*MatProp(initial_state_variable_solid_50_)
*endif
*endif
*########## material model: saturated_material_undrained_effective ##########
*if(MatProp(materialtype_id,int)==203)
*if(MatProp(effective_id,int)==500)
$$MATERIAL_MODEL_SOLID
undefined
*endif
*if(MatProp(effective_id,int)==501)
$$MATERIAL_MODEL_SOLID
linear_elasticity
$$YOUNG_MODULUS
*MatProp(effective_Young_modulus_)
$$POISSON_RATIO
*MatProp(effective_Poisson_ratio_)
$$POISSON_RATIO_UNDRAINED
*MatProp(_undrained_Poisson_ratio_)
*endif
*if(MatProp(effective_id,int)==502)
$$MATERIAL_MODEL_SOLID
mohr_coulomb
$$YOUNG_MODULUS
*MatProp(effective_Young_modulus_)
$$POISSON_RATIO
*MatProp(effective_Poisson_ratio_)
$$POISSON_RATIO_UNDRAINED
*MatProp(_undrained_Poisson_ratio_)
$$FRICTION_ANGLE
*MatProp(effective_friction_angle_)
$$COHESION
*MatProp(effective_cohesion_)
$$DILATANCY_ANGLE
*MatProp(effective_dilatancy_angle_)
$$TENSILE_STRENGTH
*MatProp(tensile_strength_)
*endif
*if(MatProp(effective_id,int)==599)
$$MATERIAL_MODEL_SOLID
external_soil_model
$$MATERIAL_MODEL_DLL
*MatProp(External_Material_Model_DLL_)
$$UMAT_DIMENSION
*MatProp(size_of_stress_tensor_inside_UMAT_)
$$MATERIAL_PARAMETER_SOLID_01
*MatProp(material_parameter_solid_01_)
$$MATERIAL_PARAMETER_SOLID_02
*MatProp(material_parameter_solid_02_)
$$MATERIAL_PARAMETER_SOLID_03
*MatProp(material_parameter_solid_03_)
$$MATERIAL_PARAMETER_SOLID_04
*MatProp(material_parameter_solid_04_)
$$MATERIAL_PARAMETER_SOLID_05
*MatProp(material_parameter_solid_05_)
$$MATERIAL_PARAMETER_SOLID_06
*MatProp(material_parameter_solid_06_)
$$MATERIAL_PARAMETER_SOLID_07
*MatProp(material_parameter_solid_07_)
$$MATERIAL_PARAMETER_SOLID_08
*MatProp(material_parameter_solid_08_)
$$MATERIAL_PARAMETER_SOLID_09
*MatProp(material_parameter_solid_09_)
$$MATERIAL_PARAMETER_SOLID_10
*MatProp(material_parameter_solid_10_)
$$MATERIAL_PARAMETER_SOLID_11
*MatProp(material_parameter_solid_11_)
$$MATERIAL_PARAMETER_SOLID_12
*MatProp(material_parameter_solid_12_)
$$MATERIAL_PARAMETER_SOLID_13
*MatProp(material_parameter_solid_13_)
$$MATERIAL_PARAMETER_SOLID_14
*MatProp(material_parameter_solid_14_)
$$MATERIAL_PARAMETER_SOLID_15
*MatProp(material_parameter_solid_15_)
$$MATERIAL_PARAMETER_SOLID_16
*MatProp(material_parameter_solid_16_)
$$MATERIAL_PARAMETER_SOLID_17
*MatProp(material_parameter_solid_17_)
$$MATERIAL_PARAMETER_SOLID_18
*MatProp(material_parameter_solid_18_)
$$MATERIAL_PARAMETER_SOLID_19
*MatProp(material_parameter_solid_19_)
$$MATERIAL_PARAMETER_SOLID_20
*MatProp(material_parameter_solid_20_)
$$MATERIAL_PARAMETER_SOLID_21
*MatProp(material_parameter_solid_21_)
$$MATERIAL_PARAMETER_SOLID_22
*MatProp(material_parameter_solid_22_)
$$MATERIAL_PARAMETER_SOLID_23
*MatProp(material_parameter_solid_23_)
$$MATERIAL_PARAMETER_SOLID_24
*MatProp(material_parameter_solid_24_)
$$MATERIAL_PARAMETER_SOLID_25
*MatProp(material_parameter_solid_25_)
$$MATERIAL_PARAMETER_SOLID_26
*MatProp(material_parameter_solid_26_)
$$MATERIAL_PARAMETER_SOLID_27
*MatProp(material_parameter_solid_27_)
$$MATERIAL_PARAMETER_SOLID_28
*MatProp(material_parameter_solid_28_)
$$MATERIAL_PARAMETER_SOLID_29
*MatProp(material_parameter_solid_29_)
$$MATERIAL_PARAMETER_SOLID_30
*MatProp(material_parameter_solid_30_)
$$MATERIAL_PARAMETER_SOLID_31
*MatProp(material_parameter_solid_31_)
$$MATERIAL_PARAMETER_SOLID_32
*MatProp(material_parameter_solid_32_)
$$MATERIAL_PARAMETER_SOLID_33
*MatProp(material_parameter_solid_33_)
$$MATERIAL_PARAMETER_SOLID_34
*MatProp(material_parameter_solid_34_)
$$MATERIAL_PARAMETER_SOLID_35
*MatProp(material_parameter_solid_35_)
$$MATERIAL_PARAMETER_SOLID_36
*MatProp(material_parameter_solid_36_)
$$MATERIAL_PARAMETER_SOLID_37
*MatProp(material_parameter_solid_37_)
$$MATERIAL_PARAMETER_SOLID_38
*MatProp(material_parameter_solid_38_)
$$MATERIAL_PARAMETER_SOLID_39
*MatProp(material_parameter_solid_39_)
$$MATERIAL_PARAMETER_SOLID_40
*MatProp(material_parameter_solid_40_)
$$MATERIAL_PARAMETER_SOLID_41
*MatProp(material_parameter_solid_41_)
$$MATERIAL_PARAMETER_SOLID_42
*MatProp(material_parameter_solid_42_)
$$MATERIAL_PARAMETER_SOLID_43
*MatProp(material_parameter_solid_43_)
$$MATERIAL_PARAMETER_SOLID_44
*MatProp(material_parameter_solid_44_)
$$MATERIAL_PARAMETER_SOLID_45
*MatProp(material_parameter_solid_45_)
$$MATERIAL_PARAMETER_SOLID_46
*MatProp(material_parameter_solid_46_)
$$MATERIAL_PARAMETER_SOLID_47
*MatProp(material_parameter_solid_47_)
$$MATERIAL_PARAMETER_SOLID_48
*MatProp(material_parameter_solid_48_)
$$MATERIAL_PARAMETER_SOLID_49
*MatProp(material_parameter_solid_49_)
$$MATERIAL_PARAMETER_SOLID_50
*MatProp(material_parameter_solid_50_)
$$INITIAL_STATE_VARIABLE_SOLID_01
*MatProp(initial_state_variable_solid_01_)
$$INITIAL_STATE_VARIABLE_SOLID_02
*MatProp(initial_state_variable_solid_02_)
$$INITIAL_STATE_VARIABLE_SOLID_03
*MatProp(initial_state_variable_solid_03_)
$$INITIAL_STATE_VARIABLE_SOLID_04
*MatProp(initial_state_variable_solid_04_)
$$INITIAL_STATE_VARIABLE_SOLID_05
*MatProp(initial_state_variable_solid_05_)
$$INITIAL_STATE_VARIABLE_SOLID_06
*MatProp(initial_state_variable_solid_06_)
$$INITIAL_STATE_VARIABLE_SOLID_07
*MatProp(initial_state_variable_solid_07_)
$$INITIAL_STATE_VARIABLE_SOLID_08
*MatProp(initial_state_variable_solid_08_)
$$INITIAL_STATE_VARIABLE_SOLID_09
*MatProp(initial_state_variable_solid_09_)
$$INITIAL_STATE_VARIABLE_SOLID_10
*MatProp(initial_state_variable_solid_10_)
$$INITIAL_STATE_VARIABLE_SOLID_11
*MatProp(initial_state_variable_solid_11_)
$$INITIAL_STATE_VARIABLE_SOLID_12
*MatProp(initial_state_variable_solid_12_)
$$INITIAL_STATE_VARIABLE_SOLID_13
*MatProp(initial_state_variable_solid_13_)
$$INITIAL_STATE_VARIABLE_SOLID_14
*MatProp(initial_state_variable_solid_14_)
$$INITIAL_STATE_VARIABLE_SOLID_15
*MatProp(initial_state_variable_solid_15_)
$$INITIAL_STATE_VARIABLE_SOLID_16
*MatProp(initial_state_variable_solid_16_)
$$INITIAL_STATE_VARIABLE_SOLID_17
*MatProp(initial_state_variable_solid_17_)
$$INITIAL_STATE_VARIABLE_SOLID_18
*MatProp(initial_state_variable_solid_18_)
$$INITIAL_STATE_VARIABLE_SOLID_19
*MatProp(initial_state_variable_solid_19_)
$$INITIAL_STATE_VARIABLE_SOLID_20
*MatProp(initial_state_variable_solid_20_)
$$INITIAL_STATE_VARIABLE_SOLID_21
*MatProp(initial_state_variable_solid_21_)
$$INITIAL_STATE_VARIABLE_SOLID_22
*MatProp(initial_state_variable_solid_22_)
$$INITIAL_STATE_VARIABLE_SOLID_23
*MatProp(initial_state_variable_solid_23_)
$$INITIAL_STATE_VARIABLE_SOLID_24
*MatProp(initial_state_variable_solid_24_)
$$INITIAL_STATE_VARIABLE_SOLID_25
*MatProp(initial_state_variable_solid_25_)
$$INITIAL_STATE_VARIABLE_SOLID_26
*MatProp(initial_state_variable_solid_26_)
$$INITIAL_STATE_VARIABLE_SOLID_27
*MatProp(initial_state_variable_solid_27_)
$$INITIAL_STATE_VARIABLE_SOLID_28
*MatProp(initial_state_variable_solid_28_)
$$INITIAL_STATE_VARIABLE_SOLID_29
*MatProp(initial_state_variable_solid_29_)
$$INITIAL_STATE_VARIABLE_SOLID_30
*MatProp(initial_state_variable_solid_30_)
$$INITIAL_STATE_VARIABLE_SOLID_31
*MatProp(initial_state_variable_solid_31_)
$$INITIAL_STATE_VARIABLE_SOLID_32
*MatProp(initial_state_variable_solid_32_)
$$INITIAL_STATE_VARIABLE_SOLID_33
*MatProp(initial_state_variable_solid_33_)
$$INITIAL_STATE_VARIABLE_SOLID_34
*MatProp(initial_state_variable_solid_34_)
$$INITIAL_STATE_VARIABLE_SOLID_35
*MatProp(initial_state_variable_solid_35_)
$$INITIAL_STATE_VARIABLE_SOLID_36
*MatProp(initial_state_variable_solid_36_)
$$INITIAL_STATE_VARIABLE_SOLID_37
*MatProp(initial_state_variable_solid_37_)
$$INITIAL_STATE_VARIABLE_SOLID_38
*MatProp(initial_state_variable_solid_38_)
$$INITIAL_STATE_VARIABLE_SOLID_39
*MatProp(initial_state_variable_solid_39_)
$$INITIAL_STATE_VARIABLE_SOLID_40
*MatProp(initial_state_variable_solid_40_)
$$INITIAL_STATE_VARIABLE_SOLID_41
*MatProp(initial_state_variable_solid_41_)
$$INITIAL_STATE_VARIABLE_SOLID_42
*MatProp(initial_state_variable_solid_42_)
$$INITIAL_STATE_VARIABLE_SOLID_43
*MatProp(initial_state_variable_solid_43_)
$$INITIAL_STATE_VARIABLE_SOLID_44
*MatProp(initial_state_variable_solid_44_)
$$INITIAL_STATE_VARIABLE_SOLID_45
*MatProp(initial_state_variable_solid_45_)
$$INITIAL_STATE_VARIABLE_SOLID_46
*MatProp(initial_state_variable_solid_46_)
$$INITIAL_STATE_VARIABLE_SOLID_47
*MatProp(initial_state_variable_solid_47_)
$$INITIAL_STATE_VARIABLE_SOLID_48
*MatProp(initial_state_variable_solid_48_)
$$INITIAL_STATE_VARIABLE_SOLID_49
*MatProp(initial_state_variable_solid_49_)
$$INITIAL_STATE_VARIABLE_SOLID_50
*MatProp(initial_state_variable_solid_50_)
*endif
*endif 
*########## material model: liquid ##########
*if(MatProp(materialtype_id,int)==206)
*if(MatProp(liquid_id,int)==300)
$$MATERIAL_MODEL_LIQUID
undefined
*endif
*if(MatProp(liquid_id,int)==301)
$$MATERIAL_MODEL_LIQUID
newtonian_liquid
*endif
*if(MatProp(liquid_id,int)==302)
$$MATERIAL_MODEL_LIQUID
bingham_liquid
$$BINGHAM_YIELD_STRESS
*MatProp(Bingham_yield_stress_)
$$BINGHAM_YOUNG_MODULUS
*MatProp(elastic_Young_modulus_)
$$BINGHAM_POISSON_RATIO
*MatProp(elastic_Poisson_ratio_)
*endif
*if(MatProp(liquid_id,int)==303)
$$MATERIAL_MODEL_LIQUID
frictional_liquid
$$LIQUID_FRICTION_ANGLE
*MatProp(friction_angle_liquid_)
$$LIQUID_YOUNG_MODULUS
*MatProp(elastic_Young_modulus_)
$$LIQUID_POISSON_RATIO
*MatProp(elastic_Poisson_ratio_)
*endif
*endif
*########## material model: unsaturated_material_2phase_suction / unsaturated_material_3phase_coupled ##########
*if((MatProp(materialtype_id,int)==207)||(MatProp(materialtype_id,int)==208))
*if(MatProp(retention_id,int)==600)
$$WATER_RETENTION_CURVE
undefined
*endif
*if(MatProp(retention_id,int)==601)
$$WATER_RETENTION_CURVE
linear
$$av
*MatProp(linear_coefficient)
*endif
*if(MatProp(retention_id,int)==602)
$$WATER_RETENTION_CURVE
van_genuchten
$$Smin
*MatProp(minimum_degree_of_saturation)
$$Smax
*MatProp(maximum_degree_of_saturation)
$$P0
*MatProp(reference_pressure)
$$L
*MatProp(lambda)
*endif
*if(MatProp(conductivity_id,int)==700)
$$HYDR_CONDUCTIVITY_CURVE
undefined
*endif
*if(MatProp(conductivity_id,int)==701)
$$HYDR_CONDUCTIVITY_CURVE
constant
*endif
*if(MatProp(conductivity_id,int)==702)
$$HYDR_CONDUCTIVITY_CURVE
hillel
$$r
*MatProp(r_exponent)
*endif
*if(MatProp(conductivity_id,int)==703)
$$HYDR_CONDUCTIVITY_CURVE
mualem
$$Smin
*MatProp(minimum_degree_of_saturation_)
$$Smax
*MatProp(maximum_degree_of_saturation_)
$$P0
*MatProp(reference_pressure_)
$$L
*MatProp(lambda_)
*endif
*endif
*end materials
*############################################
*##### element materials 2D=3D          #####
*#####   (always specified and written) #####
*############################################
$$STARTELMMAT
*loop elems
*ElemsMat 
*end elems
*############################################
*##### element local damping 2D/3D      #####
*#####   (only written when specified)  ##### 
*############################################
*##### 2D case
*set var IsLocalDamping2D=0
*if(ndime==2)
*set Cond 2D_-_Local_damping_per_surface_entity *elems 
*loop elems *OnlyInCond
*set var IsLocalDamping2D=1
*end elems
*if(IsLocalDamping2D==1)
$$STARTDAMPING
*tcl(Anura3D::WriteDamping 2D_-_Local_damping_per_surface_entity)
*end if
*end if
*##### 3D case
*set var IsLocalDamping3D=0
*if(ndime==3)
*set Cond 3D_-_Local_damping_per_volume_entity *elems 
*loop elems *OnlyInCond
*set var IsLocalDamping3D=1
*end elems
*if(IsLocalDamping3D==1)
$$STARTDAMPING
*tcl(Anura3D::WriteDamping 3D_-_Local_damping_per_volume_entity)
*end if
*end if
*##################################################
*##### material point specification 2D/3D     #####
*#####   (always specified and written)       #####
*##################################################
$$START_NUMBER_OF_MATERIAL_POINTS
*##### 2D case 
*if(ndime==2)
*##### single-point
*set Cond 2D_-_Single-point_formulation *elems
*set var IsSinglePoint2D=0
*loop elems *OnlyIncond
*set var IsSinglePoint2D=1
*end elems
*if(IsSinglePoint2D==1)
*tcl(Anura3D::WriteNumberOfMaterialPointsSP 2D_-_Single-point_formulation)
*end if
*##### double-point
*set Cond 2D_-_Double-point_formulation *elems
*set var IsDoublePoint2D=0
*loop elems *OnlyIncond
*set var IsDoublePoint2D=1
*end elems
*if(IsDoublePoint2D==1)
*tcl(Anura3D::WriteNumberOfMaterialPointsDP 2D_-_Double-point_formulation)
*end if
*end if
*##### 3D case 
*if(ndime==3)
*##### single-point
*set Cond 3D_-_Single-point_formulation *elems
*set var IsSinglePoint3D=0
*loop elems *OnlyIncond
*set var IsSinglePoint3D=1
*end elems
*if(IsSinglePoint3D==1)
*tcl(Anura3D::WriteNumberOfMaterialPointsSP 3D_-_Single-point_formulation)
*end if
*##### double-point
*set Cond 3D_-_Double-point_formulation *elems
*set var IsDoublePoint3D=0
*loop elems *OnlyIncond
*set var IsDoublePoint3D=1
*end elems
*if(IsDoublePoint3D==1)
*tcl(Anura3D::WriteNumberOfMaterialPointsDP 3D_-_Double-point_formulation)
*end if
*end if
$$FINISH