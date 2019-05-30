file(GLOB NICE_EXAMPLES_TO_INSTALL 
  thermique1d transport1d transport2d_ns transport2d_s #Old examples
  HeatEquation1DExplicit HeatEquation1DImplicit BurgersEquation1D # 1D Heat and Burgers' equations
  TransportEquation1DUpwindExplicit TransportEquation1DUpwindImplicit TransportEquation1DCenteredImplicit TransportEquation1DCenteredExplicit # 1D Transport equation 
  Poisson2DEF Poisson2DEF_disk Poisson2DVF Poisson2DVF_disk Poisson3DEF Poisson3DEF_BALL Poisson3DVF Poisson3DVF_BALL Poisson3DSphereEF Poisson3DTorusEF  
  WaveSystem1DUpwind WaveSystem1DFV_RiemannProblem 
WaveSystem_Stationary/WaveSystemUpwind WaveSystem_Stationary/WaveSystemCentered WaveSystem_Stationary/WaveSystemPStag WaveSystem_Stationary/WaveSystemStaggered 
WaveSystem_Shock/WaveSystemUpwind WaveSystem_Shock/WaveSystemCentered WaveSystem_Shock/WaveSystemPStag WaveSystem_Shock/WaveSystemStaggered 
)

install(DIRECTORY ${NICE_EXAMPLES_TO_INSTALL} DESTINATION share/examples)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
    ADD_SUBDIRECTORY(TransportEquation1DUpwindExplicit)
    ADD_SUBDIRECTORY(TransportEquation1DUpwindImplicit)
    ADD_SUBDIRECTORY(TransportEquation1DCenteredImplicit)
    ADD_SUBDIRECTORY(TransportEquation1DCenteredExplicit)
    ADD_SUBDIRECTORY(HeatEquation1DExplicit)
    ADD_SUBDIRECTORY(HeatEquation1DImplicit)
    ADD_SUBDIRECTORY(BurgersEquation1D)
    ADD_SUBDIRECTORY(Poisson2DVF)
    ADD_SUBDIRECTORY(Poisson2DVF_DISK)
    ADD_SUBDIRECTORY(Poisson2DEF)
    ADD_SUBDIRECTORY(Poisson2DEF_DISK)
    ADD_SUBDIRECTORY(Poisson3DSphereEF)
    ADD_SUBDIRECTORY(Poisson3DTorusEF)
    ADD_SUBDIRECTORY(Poisson3DEF)
    ADD_SUBDIRECTORY(Poisson3DEF_BALL)
    ADD_SUBDIRECTORY(Poisson3DVF)
    ADD_SUBDIRECTORY(Poisson3DVF_BALL)
    ADD_SUBDIRECTORY(WaveSystem1DUpwind)
    ADD_SUBDIRECTORY(WaveSystem1DFV_RiemannProblem)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemUpwind)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemPStag)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemStaggered)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemCentered)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemUpwind)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemPStag)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemStaggered)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemCentered)
ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


