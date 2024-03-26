    module PM4Sand
    !----------------------------------------------
    ! PM4Sand stress-strain model to capture the cyclic
    ! plastic deformation of soils. 
    !
    !----------------------------------------------
   
    ! Modules to be used in the PM4Sand material --> Temporary
    use ModElementEvaluationTETRA
    use ModElementEvaluationTRI
    use ModElementEvaluationQUAD
    !use ModElementEvaluation
    use ModString
    use ModFeedback
    use ModGlobalConstants
    use ModCounters
    use ModGeometryMath
    use ModNURBS

    implicit none

    contains
    
    subroutine PM4SandMaterial()
    !----------------------------------------------
    ! To have the variables that are needed as constants 
    ! in the PM4Sand implementation. 
    !----------------------------------------------
    implicit none 
    
    real(REAL_TYPE), dimension(24) :: oData
    
    ! initialization of list of parameters
    oData(1) = 101.3    !// P_atm       ! Numbering according to the manual
	oData(2) = -1       !// h0          ! parameter 1
	oData(3) = 0.8      !// emax        ! parameter 2
	oData(4) = 0.5      !// emin        ! parameter 2
	oData(5) = 0.5      !// nb          ! parameter 3
	oData(6) = 0.1      !// nd          ! parameter 4
	oData(7) = -1       !// Ado         ! parameter 5
	oData(8) = -1       !// z_max       ! parameter 6
	oData(9) = 250.0    !// cz          ! parameter 7 
	oData(10) = -1       !// ce         ! parameter 8
	oData(11) = 33.0    !// phi_cv      ! parameter 9
	oData(12) = 0.3     !// nu          ! parameter 10
	oData(13) = 2.0     !// Cgd         ! parameter 11
	oData(14) = -1      !// C_DR        ! parameter 
	oData(15) = -1      !// Ckaf        ! parameter 12
	oData(16) = 10.0    !// Q           ! parameter 13
	oData(17) = 1.5     !// R           ! parameter 14
	oData(18) = 0.01    !// m           ! parameter 15
	oData(19) = -1      !// Fsed_min    ! 
	oData(20) = -1      !//p_sdeo       !
	oData(21) = 5		 !// IntScheme  !
	oData(22) = 0		 !// TanType    !
	oData(23) = 1.0e-8	 !// TolF       !
	oData(24) = 1.0e-8	 !// TolR       !
    
    
    
    end subroutine PM4SandMaterial
    
    
    
    
    
    subroutine FullConstructor(Dr, G0, hp0, mDen, P_atm, &
                                       h0, emax, emin, nb, nd, &
                                       Ado, z_max, cz, ce, &
                                       phi_cv, nu, &
                                       Cgd, Cdr, Ckaf, &
                                       QQ, RR, &
                                       mm, Fsed_min, p_sdeo, &
                                       integrationScheme, tangentType, TolF, TolR)
    !----------------------------------------------
    ! To populate the stress-strain variables needed 
    ! 
    !---------------------------------------------- 
    implicit none 
    
    ! PM4SAND input variables 
    real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1
    real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2
    real(REAL_TYPE), intent(in) :: hp0      ! primiary parameter 3
    real(REAL_TYPE), intent(in) :: mDen     
    real(REAL_TYPE), intent(in) :: P_atm    ! primiary parameter 4
    real(REAL_TYPE), intent(in) :: h0                               ! secondary parameter 1
    real(REAL_TYPE), intent(in) :: emax                             ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: emin                             ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: nb                               ! secondary parameter 3
    real(REAL_TYPE), intent(in) :: nd                               ! secondary parameter 4
    real(REAL_TYPE), intent(in) :: Ado                              ! secondary parameter 5
    real(REAL_TYPE), intent(in) :: z_max                            ! secondary parameter 6
    real(REAL_TYPE), intent(in) :: cz                               ! secondary parameter 7
    real(REAL_TYPE), intent(in) :: ce                               ! secondary parameter 8
    real(REAL_TYPE), intent(in) :: phi_cv                           ! secondary parameter 9
    real(REAL_TYPE), intent(in) :: nu                               ! secondary parameter 10
    real(REAL_TYPE), intent(in) :: Cgd                              ! secondary parameter 11
    real(REAL_TYPE), intent(in) :: Cdr                              ! secondary parameter 
    real(REAL_TYPE), intent(in) :: Ckaf                             ! secondary parameter 12
    real(REAL_TYPE), intent(in) :: QQ                               ! secondary parameter 13
    real(REAL_TYPE), intent(in) :: RR                               ! secondary parameter 14
    real(REAL_TYPE), intent(in) :: mm                               ! secondary parameter 15
    real(REAL_TYPE), intent(in) :: Fsed_min                         ! secondary parameter 
    real(REAL_TYPE), intent(in) :: p_sdeo                           ! secondary parameter 1
    integer(INTEGER_TYPE), intent(in) :: integrationScheme          ! secondary parameter 1
    integer(INTEGER_TYPE), intent(in) :: tangentType
    real(REAL_TYPE), intent(in) :: TolF
    real(REAL_TYPE), intent(in) :: TolR
    
    
    
    
    !PM4SAND material variables 
    real(REAL_TYPE) :: mEpsilon(3)                  ! new strain ! --> this needs to be an input to this subroutine
	real(REAL_TYPE) :: mEpsilon_n(3)                ! current strain 
	real(REAL_TYPE) :: mEpsilonE(3)                 ! new elastic strain 
	real(REAL_TYPE) :: mEpsilonE_n(3)               ! current elastic strain 
	real(REAL_TYPE) :: mSigma(3)                    ! new stress
	real(REAL_TYPE) :: mSigma_b(3)                  ! new 
	real(REAL_TYPE) :: mSigma_n(3)                  ! current strain      
	real(REAL_TYPE) :: mAlpha(3)                    ! new alpha vector 
	real(REAL_TYPE) :: mAlpha_n(3)                  ! current alpha vector
	real(REAL_TYPE) :: mAlpha_in(3)                 ! new initial alpha vector
	real(REAL_TYPE) :: mAlpha_in_n(3)               ! current initial alpha vector
	real(REAL_TYPE) :: mAlpha_in_p(3)               ! 
	real(REAL_TYPE) :: mAlpha_in_p_n(3)             ! current
	real(REAL_TYPE) :: mAlpha_in_true(3)            !
	real(REAL_TYPE) :: mAlpha_in_true_n(3)          ! current
	real(REAL_TYPE) :: mAlpha_in_max(3)             !
	real(REAL_TYPE) :: mAlpha_in_max_n(3)           ! current
	real(REAL_TYPE) :: mAlpha_in_min(3)             !
	real(REAL_TYPE) :: mAlpha_in_min_n(3)           ! current
	real(REAL_TYPE) :: mFabric(3)                   ! new fabric vector
	real(REAL_TYPE) :: mFabric_n(3)                 ! current fabric vector 
	real(REAL_TYPE) :: mFabric_in(3)                ! new initial fabric vector 
	real(REAL_TYPE) :: mFabric_in_n(3)              ! current initial fabric vector
	real(REAL_TYPE) :: mCe(3, 3)                    ! stiffness matrix
	real(REAL_TYPE) :: mCep(3, 3)                   ! stiffness matrix
	real(REAL_TYPE) :: mCep_Consistent(3, 3)        ! stiffness matrix
	real(REAL_TYPE) :: mTracker(3)
    
    !real(REAL_TYPE) :: mEpsilonE(3), mEpsilonE_n(3)
    
    
    !PM4SAND material variables 
    real(REAL_TYPE) :: m_Dr
    real(REAL_TYPE) :: m_G0
    real(REAL_TYPE) :: m_hpo
    real(REAL_TYPE) :: massDen
    real(REAL_TYPE) :: m_P_atm
    real(REAL_TYPE) :: m_h0_1
    real(REAL_TYPE) :: m_h0_2
    real(REAL_TYPE) :: m_h0
    
    real(REAL_TYPE) :: m_emax
    real(REAL_TYPE) :: m_emin
    real(REAL_TYPE) :: m_nb
    real(REAL_TYPE) :: m_nd
    
    real(REAL_TYPE) :: m_Ado
    real(REAL_TYPE) :: m_z_max
    
    real(REAL_TYPE) :: m_cz
    real(REAL_TYPE) :: m_ce
    real(REAL_TYPE) :: m_Mc
    real(REAL_TYPE) :: m_nu
    real(REAL_TYPE) :: m_Cgd
    real(REAL_TYPE) :: m_Cdr
    real(REAL_TYPE) :: m_Ckaf
    real(REAL_TYPE) :: m_QQ
    real(REAL_TYPE) :: m_RR
    
    real(REAL_TYPE) :: m_m
    real(REAL_TYPE) :: m_Fsed_min
    real(REAL_TYPE) :: m_p_sedo
    real(REAL_TYPE) :: m_FirstCall
    real(REAL_TYPE) :: m_PostShake
    
    real(REAL_TYPE) :: mScheme
    real(REAL_TYPE) :: mTangtype
    real(REAL_TYPE) :: mTolF
    real(REAL_TYPE) :: mTolR
    
    real(REAL_TYPE) :: m_e_init
    real(REAL_TYPE) :: miter
    
    ! new local variables 
    real(REAL_TYPE) :: mMb
    real(REAL_TYPE) :: mMd
    real(REAL_TYPE) :: mKp
    real(REAL_TYPE), dimension(3) :: rrr
    real(REAL_TYPE), dimension(3) :: initStress
    real(REAL_TYPE) :: mzxp
    real(REAL_TYPE) :: mzcum
    real(REAL_TYPE) :: mpzp
    real(REAL_TYPE) :: mzpeak
    real(REAL_TYPE) :: mK
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mMcur
    
    
    logical :: m_pzpFlag
    
    m_Dr = Dr
    m_G0 = G0
    m_hpo = hp0
    massDen = mDen
    
    if (P_atm<0) then
    m_P_atm = 101.3
    end if
    
    if (h0<0) then 
    m_h0_1 = 0.3
    m_h0_2 = (0.25+m_Dr)/2
    m_h0 = max(m_h0_1, m_h0_2)
    end if 
    
    if (emax<0) then 
        m_emax = 0.8    
    else 
        m_emax = emax
    end if 
    
    if (emin<0) then 
        m_emin = 0.5 
    else 
        m_emin = emin
    end if 
    
    if (nb<0) then 
        m_nb = 0.5
    else
        m_nb = nb
    end if 
    
    if (nd<0) then 
        m_nd = 0.5
    else
        m_nd = nd
    end if
    
    m_Ado = Ado
    m_z_max = z_max
    
    if (cz<0) then 
        m_cz = 250
    else 
        m_cz = cz
    end if 
    
    
    if (ce>0) then 
        m_ce = ce
    else 
        if (m_Dr > 0.75) then 
			m_ce = 0.2
		else if (m_Dr < 0.55) then 
			m_ce = 0.5
		else
			m_ce = 0.5 - ((m_Dr - 0.55) * 1.5)
        end if
                                       
    end if 
                                       
    if (phi_cv<0) then 
       m_Mc = 2.0 * sin(33.0/180.0 * 3.14159265359)      
    else 
       m_Mc = 2.0 * sin(phi_cv/180.0 * 3.14159265359)
    end if 
    
    if (nu<0) then 
        m_nu = 0.3
    else 
        m_nu = nu 
    end if 
    
    if (Cgd < 0) then 
        m_Cgd = 2.0
    else 
        m_Cgd = Cgd
    end if 
    
    if (Cdr < 0.0) then 
        m_Cdr = (5 + 25 * (m_Dr - 0.35))
    else 
        m_Cdr = Cdr
    end if 
    m_Cdr = min(m_Cdr, 10.0)
    
    if (Ckaf < 0) then 
        m_Ckaf = (5.0 + 220.0 * (m_Dr - 0.26)**3)
    else 
        m_Ckaf = Ckaf
    end if 
    
    if (m_Ckaf > 35) then 
        m_Ckaf = 35.0
    else 
        m_Ckaf = m_Ckaf
    end if 
    
    if (m_Ckaf < 4) then 
        m_Ckaf = 4.0
    else 
        m_Ckaf = m_Ckaf
    end if 
    
    if (QQ < 0) then 
        m_QQ = 10.0
    else 
        m_QQ = QQ
    end if 
    
    if (RR < 0) then 
        m_RR = 1.5
    else 
        m_RR = RR
    end if  
    
    if (mm < 0) then 
        m_m = 0.01
    else 
        m_m = mm
    end if
    
    if (Fsed_min<0.0) then 
        m_Fsed_min = (0.03 * exp(2.6 * m_Dr))
    else 
        m_Fsed_min = Fsed_min
    end if 
    m_Fsed_min = min(m_Fsed_min, 0.99)
    
    if (p_sdeo < 0.0) then 
        m_p_sedo = (m_P_atm / 5.0)
    else 
        m_p_sedo = p_sdeo
    end if 
    
    m_FirstCall = 0
    m_PostShake = 0
    
    mScheme = integrationScheme
    mTangType = tangentType
	mTolF = TolF
	mTolR = TolR
    
    m_e_init = m_emax - (m_emax - m_emin) * m_Dr
    mIter = 0
    
    call initialize(initStress, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m)
    
    ! we obtain the strain, mEpsilon, from the element and then input it into the integrate matrix
    ! we will find the strain increment inside the the integrate subroutine
    call setTrialStrain(mEpsilon, mEpsilon_n, &
        mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
        mEpsilonE, mEpsilonE_n) !initStress, 
    
    end subroutine FullConstructor
    
    
 !   subroutine CommitState() 
 !   !-------------------------------------------------
 !   ! To commit stress state and state parameters
 !   !
 !   !-------------------------------------------------
 !   implicit none 
 !   
 !   call GetElasticModuli_(mSigma, mK, mG, mMcur, mzcum) !--> temporary subroutine
 !   
 !   call GetTrace(mSigma, GetTrace_mSigma) !--> temporary subroutine
 !   
 !   if ((mMcur>mMb) .and. me2p) then
 !       pp = 0.5 * mSigma_Trace ! real type 
 !       r = (mSigma - (p*mI1)) * (1/pp) * (mMb/mMcur)
 !       mSigma = (p*mI1) + (r*p)
 !       mAlpha = r * (mMb - m_m)/mMb
 !   end if 
 !   mAlpha_in_n = mAlpha_in
	!mAlpha_n = mAlpha
	!mAlpha_in_p_n = mAlpha_in_p
	!mAlpha_in_true_n = mAlpha_in_true
	!mAlpha_in_max_n = mAlpha_in_max
	!mAlpha_in_min_n = mAlpha_in_min
	!mSigma_n = mSigma
	!mEpsilon_n = mEpsilon
	!mEpsilonE_n = mEpsilonE
	!dFabric = mFabric - mFabric_n
 !   
 !   ! update cumulated fabric 
 !   call DoubleDot2_2_Contr(dFabric, dFabric, DoubleDot2_2_Contr_dFabric_dFabric) !--> temporary subroutine
 !   mzcum = mzcum + sqrt(DoubleDot2_2_Contr_dFabric_dFabric/2.0)
 !   call DoubleDot2_2_Contr(mFabric, mFabric, DoubleDot2_2_Contr_mFabric_mFabric) !--> temporary subroutine
 !   mzpeak = max( sqrt(DoubleDot2_2_Contr_mFabric_mFabric/2.0), mzpeak)
 !   
 !   mFabric_n = mFabric
 !   mFabric_in_n = mFabric_in
 !   mDGamma_n = mDGamma 
 !   call GetTrace(mEpsilon, GetTrace_mEpsilon) !--> temporary subroutine
 !   mVoidRatio = m_e_init - ( (1 + m_e_init) * GetTrace_mEpsilon )
 !   
 !   call GetStiffness(mK, mG, mCe) !--> temporary subroutine
 !   call GetElastoPlasticTangent(mSigma_n, mCe, R, n, mKp, mCep) !--> temporary subroutine
 !   mCep_Consistent = mCe
 !   
 !   
 !   end subroutine CommitState
                                       
                                       
              
                                       
    subroutine commitState()
    implicit none
    ! Local variables
    real(REAL_TYPE) :: pp
    real(REAL_TYPE), dimension(3) :: nn, RR, dFabric
    real(REAL_TYPE), dimension(3) :: rrr
    
    ! Local variables 
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mK 
    real(REAL_TYPE) :: mMcur
    real(REAL_TYPE), dimension(3) :: mSigma
    real(REAL_TYPE), dimension(3) :: mSigma_n
    integer(INTEGER_TYPE) :: me2p
    real(REAL_TYPE) :: mMb
    real(REAL_TYPE) :: mzcum
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE), dimension(3) :: mI1
    real(REAL_TYPE), dimension(3) :: mAlpha
    real(REAL_TYPE) :: m_m
    real(REAL_TYPE) :: mAlpha_in
    real(REAL_TYPE) :: mAlpha_in_n
    real(REAL_TYPE), dimension(3) :: mAlpha_n
    real(REAL_TYPE) :: mAlpha_in_p_n
    real(REAL_TYPE) :: mAlpha_in_p
    real(REAL_TYPE) :: mAlpha_in_true_n
    real(REAL_TYPE) :: mAlpha_in_true
    real(REAL_TYPE) :: mAlpha_in_max_n
    real(REAL_TYPE) :: mAlpha_in_max
    real(REAL_TYPE) :: mAlpha_in_min
    real(REAL_TYPE) :: mAlpha_in_min_n
    
    !real(REAL_TYPE), dimension(3) :: mSigma_n
    real(REAL_TYPE), dimension(3) :: mEpsilonE_n
    real(REAL_TYPE), dimension(3) :: mEpsilon
    real(REAL_TYPE), dimension(3) :: mEpsilonE
    real(REAL_TYPE), dimension(3) :: mEpsilon_n
    
    real(REAL_TYPE), dimension(3) :: mFabric_n
    real(REAL_TYPE), dimension(3) :: mFabric
    real(REAL_TYPE), dimension(3) :: mFabric_in_n
    real(REAL_TYPE), dimension(3) :: mFabric_in
    real(REAL_TYPE) :: mvoidratio
    real(REAL_TYPE) :: m_e_init
    real(REAL_TYPE), dimension(3,3) :: CC
    real(REAL_TYPE), dimension(3,3) :: mCe
    real(REAL_TYPE) :: mKp
    real(REAL_TYPE), dimension(3,3) :: mCep
    real(REAL_TYPE), dimension(3,3) :: mCep_consistent
    !real(REAL_TYPE), dimension(3) :: mFabric_in_n
    real(REAL_TYPE) :: m_pmin
    real(REAL_TYPE) :: m_p_atm
    real(REAL_TYPE) :: m_pmin2
    !real(REAL_TYPE) :: mzcum 
    real(REAL_TYPE) :: mzpeak 
    real(REAL_TYPE) :: m_z_max
    !real(REAL_TYPE) :: mK
    !real(REAL_TYPE) :: mG
    
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
  
    !real(REAL_TYPE) :: mzpeak
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    ! Get elastic moduli
    call GetElasticModuli_(mSigma, mK, mG, mMcur, mzcum)

    ! Check condition for plasticity
    if ((mMcur > mMb) .and. (me2p==1)) then
        call GetTrace(mSigma, GetTrace_result)
        pp = 0.5 * GetTrace_result!(mSigma)
        rrr = (mSigma - pp * mI1) * (mMb / mMcur / pp)
        mSigma = pp * mI1 + rrr * pp
        mAlpha = rrr * (mMb - m_m) / mMb
    endif

    ! Update state variables
    mAlpha_in_n = mAlpha_in
    mAlpha_n = mAlpha
    mAlpha_in_p_n = mAlpha_in_p
    mAlpha_in_true_n = mAlpha_in_true
    mAlpha_in_max_n = mAlpha_in_max
    mAlpha_in_min_n = mAlpha_in_min
    mSigma_n = mSigma
    mEpsilon_n = mEpsilon
    mEpsilonE_n = mEpsilonE

    ! Calculate change in fabric
    dFabric = mFabric - mFabric_n

    ! Update cumulated fabric
    call DoubleDot2_2_Contr(dFabric, dFabric, DoubleDot2_2_Contr_result)
    mzcum = mzcum + sqrt(DoubleDot2_2_Contr_result / 2.0)
    call DoubleDot2_2_Contr(mFabric, mFabric, DoubleDot2_2_Contr_result)
    mzpeak = max(sqrt(DoubleDot2_2_Contr_result / 2.0), mzpeak)

    ! Update fabric variables
    mFabric_n = mFabric
    mFabric_in_n = mFabric_in

    ! Update void ratio
    call GetTrace(mEpsilon, GetTrace_result)
    mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace_result!GetTrace(mEpsilon)

    ! Calculate stiffness matrices
    call GetStiffness(mK, mG, CC)
    mCe = CC!GetStiffness(mK, mG)
    call GetElastoPlasticTangent(mSigma_n, mCe, RR, nn, mKp, mCep)
    mCep_Consistent = mCe
    
    end subroutine commitState


    !---------------------------------------------------------------------
    


    subroutine revertToLastCommit()
    ! Need to be added
    end subroutine revertToLastCommit

!---------------------------------------------------------------------


    !subroutine revertToStart()
    !! Added: C.McGann, U.Washington for InitialStateAnalysis
    !if (ops_InitialStateAnalysis) then
    !    ! Do nothing, keep state variables from last step
    !else
    !    ! Normal call for revertToStart (not initialStateAnalysis)
    !    call initialize(mSigma)
    !endif
    !
    !end subroutine revertToStart

!---------------------------------------------------------------------


    !subroutine getCopy(clone)
    !! Outputs
    !type(PM4Sand), pointer :: clone
    !! Create a clone
    !clone => PM4Sand()
    !clone = this
    !
    !end subroutine getCopy

!---------------------------------------------------------------------

    subroutine getType(materialType)
    ! Outputs
    character(LEN=*) :: materialType
    ! Set material type
    materialType = "PlaneStrain"
    end subroutine getType

!---------------------------------------------------------------------

    subroutine getOrder(materialOrder)
    ! Outputs
    integer, intent(out) :: materialOrder
    ! Set material order
    materialOrder = 3
    end subroutine getOrder

!---------------------------------------------------------------------

    !subroutine setResponse(argv, argc, output, responseObj)
    !! Inputs
    !character(LEN=*), dimension(*), intent(in) :: argv
    !integer, intent(in) :: argc
    !type(OPS_Stream), intent(inout) :: output
    !! Outputs
    !type(Response), pointer :: responseObj
    !! Set response based on input arguments
    !if (strcmp(argv(1), "stress") == 0 .or. strcmp(argv(1), "stresses") == 0) then
    !    responseObj => MaterialResponse(this, 1, this%getStress())
    !elseif (strcmp(argv(1), "strain") == 0 .or. strcmp(argv(1), "strains") == 0) then
    !    responseObj => MaterialResponse(this, 2, this%getStrain())
    !elseif (strcmp(argv(1), "state") == 0) then
    !    responseObj => MaterialResponse(this, 3, this%getState())
    !elseif (strcmp(argv(1), "alpha") == 0 .or. strcmp(argv(1), "backstressratio") == 0) then
    !    responseObj => MaterialResponse(this, 4, this%getAlpha())
    !elseif (strcmp(argv(1), "fabric") == 0) then
    !    responseObj => MaterialResponse(this, 5, this%getFabric())
    !elseif (strcmp(argv(1), "alpha_in") == 0 .or. strcmp(argv(1), "alphain") == 0) then
    !    responseObj => MaterialResponse(this, 6, this%getAlpha_in())
    !elseif (strcmp(argv(1), "trackers") == 0 .or. strcmp(argv(1), "tracker") == 0) then
    !    responseObj => MaterialResponse(this, 7, this%getTracker())
    !else
    !    responseObj => NULL()
    !endif
    !end subroutine setResponse
                                       
                                       
                                       
    
    !------------------------------------------------------------------------
    
!    subroutine getResponse(responseID, matInfo)
!    use Information
!    implicit none
!    ! Inputs
!    integer, intent(in) :: responseID
!    type(Information), intent(inout) :: matInfo
!    ! Local variables
!    integer :: res
!
!    select case (responseID)
!    case (-1)
!        res = -1
!    case (1)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getStress()
!        endif
!        res = 0
!    case (2)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getStrain()
!        endif
!        res = 0
!    case (3)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getState()
!        endif
!        res = 0
!    case (4)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getAlpha()
!        endif
!        res = 0
!    case (5)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getFabric()
!        endif
!        res = 0
!    case (6)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getAlpha_in()
!        endif
!        res = 0
!    case (7)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getTracker()
!        endif
!        res = 0
!    case default
!        res = -1
!    end select
!
!end subroutine getResponse
    
        
    !------------------------------------------------------------------------
    
    !subroutine sendSelf(commitTag, theChannel)
    !use Channel
    !implicit none
    !! Inputs
    !integer, intent(in) :: commitTag
    !type(Channel), intent(in) :: theChannel
    !! Local variables
    !integer :: res
    !integer, parameter :: dataSize = 101
    !real(REAL_TYPE) :: data(dataSize)
    !
    !! Assign data to vector
    !data(1) = real(this%getTag())
    !
    !data(2) = this%m_Dr
    !data(3) = this%m_G0
    !data(4) = this%m_hpo
    !data(5) = this%massDen
    !data(6) = this%m_P_atm
    !data(7) = this%m_h0
    !data(8) = this%m_emax
    !data(9) = this%m_emin
    !data(10) = this%m_e_init
    !data(11) = this%m_nb
    !data(12) = this%m_nd
    !data(13) = this%m_Ado
    !data(14) = this%m_cz
    !data(15) = this%m_ce
    !data(16) = this%m_Mc
    !data(17) = this%m_nu
    !data(18) = this%m_Cgd
    !data(19) = this%m_Cdr
    !data(20) = this%m_Ckaf
    !data(21) = this%m_Q
    !data(22) = this%m_R
    !data(23) = this%m_m
    !data(24) = this%m_z_max
    !data(25) = this%m_Fsed_min
    !data(26) = this%m_p_sedo
    !data(27) = this%m_FirstCall
    !data(28) = this%m_PostShake
    !
    !data(29) = this%mTolF
    !data(30) = this%mTolR
    !data(31) = this%mScheme
    !data(32) = this%mTangType
    !data(33) = this%m_Pmin
    !data(34) = this%m_Pmin2
    !data(35) = this%m_pzpFlag
    !data(36) = this%me2p
    !
    !data(37) = this%mDGamma
    !data(38) = this%mDGamma_n
    !data(39) = this%mK
    !data(40) = this%mG
    !data(41) = this%mVoidRatio
    !data(42) = this%mKp
    !data(43) = this%mzcum
    !data(44) = this%mzpeak
    !data(45) = this%mpzp
    !data(46) = this%mzxp
    !data(47) = this%mMb
    !data(48) = this%mMd
    !data(49) = this%mMcur
    !
    !data(50) = this%mEpsilon(1);     data(53) = this%mEpsilon_n(1);      data(56) = this%mSigma(1);     data(59) = this%mSigma_n(1);    data(62) = this%mSigma_b(1)
    !data(51) = this%mEpsilon(2);     data(54) = this%mEpsilon_n(2);      data(57) = this%mSigma(2);     data(60) = this%mSigma_n(2);    data(63) = this%mSigma_b(2)
    !data(52) = this%mEpsilon(3);     data(55) = this%mEpsilon_n(3);      data(58) = this%mSigma(3);     data(61) = this%mSigma_n(3);    data(64) = this%mSigma_b(3)
    !
    !data(65) = this%mEpsilonE(1);    data(68) = this%mEpsilonE_n(1);     data(71) = this%mAlpha(1);     data(74) = this%mAlpha_n(1);    data(77) = this%mAlpha_in_n(1)
    !data(66) = this%mEpsilonE(2);    data(69) = this%mEpsilonE_n(2);     data(72) = this%mAlpha(2);     data(75) = this%mAlpha_n(2);    data(78) = this%mAlpha_in_n(2)
    !data(67) = this%mEpsilonE(3);    data(70) = this%mEpsilonE_n(3);     data(73) = this%mAlpha(3);     data(76) = this%mAlpha_n(3);    data(79) = this%mAlpha_in_n(3)
    !
    !data(80) = this%mAlpha_in_p_n(1); data(83) = this%mAlpha_in_true_n(1); data(86) = this%mAlpha_in_max_n(1);   data(89) = this%mAlpha_in_min_n(1)
    !data(81) = this%mAlpha_in_p_n(2); data(84) = this%mAlpha_in_true_n(2); data(87) = this%mAlpha_in_max_n(2);   data(90) = this%mAlpha_in_min_n(2)
    !data(82) = this%mAlpha_in_p_n(3); data(85) = this%mAlpha_in_true_n(3); data(88) = this%mAlpha_in_max_n(3);   data(91) = this%mAlpha_in_min_n(3)
    !
    !data(92) = this%mFabric(1);        data(95) = this%mFabric_n(1);      data(98) = this%mFabric_in_n(1)
    !data(93) = this%mFabric(2);        data(96) = this%mFabric_n(2);      data(99) = this%mFabric_in_n(2)
    !data(94) = this%mFabric(3);        data(97) = this%mFabric_n(3);      data(100) = this%mFabric_in_n(3)
    !
    !! Send data vector via the channel
    !!res
    !
    !
    !end subroutine 
    !------------------------------------------------------------------------
    
!    subroutine recvSelf(commitTag, theChannel, theBroker)
!    use Channel
!    use FEM_ObjectBroker
!    implicit none
!    ! Inputs
!    integer, intent(in) :: commitTag
!    type(Channel), intent(in) :: theChannel
!    type(FEM_ObjectBroker), intent(in) :: theBroker
!    ! Local variables
!    integer :: res
!    integer, parameter :: dataSize = 101
!    real(REAL_TYPE) :: data(dataSize)
!
!    ! Receive data from the channel
!    res = theChannel%recvVector(this%getDbTag(), commitTag, data)
!    if (res < 0) then
!        print *, "WARNING: PM4Sand::recvSelf - failed to receive vector from channel"
!        return
!    endif
!
!    ! Set received data to class attributes
!    this%setTag(int(data(1)))
!
!    ! Set other attributes using received data
!    this%m_Dr = data(2)
!    this%m_G0 = data(3)
!    this%m_hpo = data(4)
!    this%massDen = data(5)
!    this%m_P_atm = data(6)
!    this%m_h0 = data(7)
!    this%m_emax = data(8)
!    this%m_emin = data(9)
!    this%m_e_init = data(10)
!    this%m_nb = data(11)
!    this%m_nd = data(12)
!    this%m_Ado = data(13)
!    this%m_cz = data(14)
!    this%m_ce = data(15)
!    this%m_Mc = data(16)
!    this%m_nu = data(17)
!    this%m_Cgd = data(18)
!    this%m_Cdr = data(19)
!    this%m_Ckaf = data(20)
!    this%m_Q = data(21)
!    this%m_R = data(22)
!    this%m_m = data(23)
!    this%m_z_max = data(24)
!    this%m_Fsed_min = data(25)
!    this%m_p_sedo = data(26)
!    this%m_FirstCall = data(27)
!    this%m_PostShake = data(28)
!
!    this%mTolF = data(29)
!    this%mTolR = data(30)
!    this%mScheme = data(31)
!    this%mTangType = data(32)
!    this%m_Pmin = data(33)
!    this%m_Pmin2 = data(34)
!    this%m_pzpFlag = data(36)
!    this%me2p = data(37)
!
!    this%mDGamma = data(38)
!    this%mDGamma_n = data(39)
!    this%mK = data(40)
!    this%mG = data(41)
!    this%mVoidRatio = data(42)
!    this%mKp = data(43)
!    this%mzcum = data(44)
!    this%mzpeak = data(45)
!    this%mpzp = data(46)
!    this%mzxp = data(47)
!    this%mMb = data(48)
!    this%mMd = data(49)
!    this%mMcur = data(50)
!
!    this%mEpsilon(1) = data(51);     this%mEpsilon_n(1) = data(52);      this%mSigma(1) = data(53);     this%mSigma_n(1) = data(54);    this%mSigma_b(1) = data(55)
!    this%mEpsilon(2) = data(56);     this%mEpsilon_n(2) = data(57);      this%mSigma(2) = data(58);     this%mSigma_n(2) = data(59);    this%mSigma_b(2) = data(60)
!    this%mEpsilon(3) = data(61);     this%mEpsilon_n(3) = data(62);      this%mSigma(3) = data(63);     this%mSigma_n(3) = data(64);    this%mSigma_b(3) = data(65)
!
!    this%mEpsilonE(1) = data(66);    this%mEpsilonE_n(1) = data(67);     this%mAlpha(1) = data(68);     this%mAlpha_n(1) = data(69);    this%mAlpha_in_n(1) = data(70)
!    this%mEpsilonE(2) = data(71);    this%mEpsilonE_n(2) = data(72);     this%mAlpha(2) = data(73);     this%mAlpha_n(2) = data(74);    this%mAlpha_in_n(2) = data(75)
!    this%mEpsilonE(3) = data(76);    this%mEpsilonE_n(3) = data(77);     this%mAlpha(3) = data(78);     this%mAlpha_n(3) = data(79);    this%mAlpha_in_n(3) = data(80)
!
!    this%mAlpha_in_p_n(1) = data(81); this%mAlpha_in_true_n(1) = data(82); this%mAlpha_in_max_n(1) = data(83); this%mAlpha_in_min_n(1) = data(84)
!    this%mAlpha_in_p_n(2) = data(85); this%mAlpha_in_true_n(2) = data(86); this%mAlpha_in_max_n(2) = data(87); this%mAlpha_in_min_n(2) = data(88)
!    this%mAlpha_in_p_n(3) = data(89); this%mAlpha_in_true_n(3) = data(90); this%mAlpha_in_max_n(3) = data(91); this%mAlpha_in_min_n(3) = data(92)
!
!    this%mFabric(1) = data(93);       this%mFabric_n(1) = data(94);       this%mFabric_in_n(1) = data(95)
!    this%mFabric(2) = data(96);       this%mFabric_n(2) = data(97);       this%mFabric_in_n(2) = data(98)
!    this%mFabric(3) = data(99);       this%mFabric_n(3) = data(100);      this%mFabric_in_n(3) = data(101)
!end subroutine recvSelf
    
    
    !------------------------------------------------------------------------
    
!    subroutine Print(s, flag)
!    use OPS_Stream
!    implicit none
!    ! Inputs
!    type(OPS_Stream), intent(inout) :: s
!    integer, intent(in) :: flag
!
!    ! Print information to the stream
!    call s%write("PM4Sand Material, tag: ")
!    call s%write(this%getTag())
!    call s%write("Type: ")
!    call s%write(this%getType())
!    call s%write(endln)
!
!end subroutine Print
    
    !------------------------------------------------------------------------
    
!    subroutine SetParameter(argv, argc, param, result)
!    implicit none
!    ! Inputs
!    character(len=*), dimension(*), intent(in) :: argv
!    integer, intent(in) :: argc
!    type(Parameter), intent(inout) :: param
!    ! Output
!    integer, intent(out) :: result
!
!    ! Local variables
!    integer :: theMaterialTag
!
!    ! Check if argc is less than 2
!    if (argc < 2) then
!        result = -1
!        return
!    end if
!
!    ! Convert argv[1] to integer
!    read(argv(1), *) theMaterialTag
!
!    ! Check if theMaterialTag matches the current object's tag
!    if (theMaterialTag == this%getTag()) then
!        ! Check argv[0] for different cases and add objects accordingly
!        if (trim(argv(0)) == "updateMaterialStage") then
!            print *, this%getTag(), "update Material Stage"
!            result = param%addObject(1, this)
!        else if (trim(argv(0)) == "materialState") then
!            result = param%addObject(5, this)
!        else if (trim(argv(0)) == "IntegrationScheme") then
!            result = param%addObject(2, this)
!        else if (trim(argv(0)) == "refShearModulus" .or. trim(argv(0)) == "ShearModulus") then
!            result = param%addObject(6, this)
!        else if (trim(argv(0)) == "poissonRatio") then
!            result = param%addObject(7, this)
!        else if (trim(argv(0)) == "FirstCall") then
!            result = param%addObject(8, this)
!        else if (trim(argv(0)) == "voidRatio") then
!            result = param%addObject(9, this)
!        else if (trim(argv(0)) == "PostShake") then
!            result = param%addObject(13, this)
!        end if
!    else
!        result = -1
!    end if
!
!end subroutine SetParameter
    
    !------------------------------------------------------------------------
!    subroutine UpdateParameter(responseID, info)
!    implicit none
!    ! Inputs
!    integer, intent(in) :: responseID
!    type(Information), intent(in) :: info
!
!    ! Local variables
!    real(REAL_TYPE) :: eps_v
!
!    ! Update material parameters based on responseID
!    if (responseID == 1) then
!        me2p = info%theInt
!    else if (responseID == 5) then
!        me2p = int(info%theDouble)
!    else if (responseID == 2) then
!        mScheme = int(info%theDouble)
!    else if (responseID == 6) then
!        m_G0 = info%theDouble
!    else if (responseID == 7) then
!        m_nu = info%theDouble
!    else if (responseID == 8) then
!        m_FirstCall = info%theInt
!        call initialize(mSigma_n)
!        print *, this%getTag(), "initialize"
!    else if (responseID == 9) then
!        eps_v = GetTrace(mEpsilon)
!        m_e_init = (info%theDouble + eps_v) / (1.0 - eps_v)
!    else if (responseID == 13) then
!        m_PostShake = 1
!        ! mElastFlag = 1
!        call GetElasticModuli(mSigma, mK, mG, mMcur, mzcum)
!        print *, this%getTag(), "activate post shaking reconsolidation"
!    else
!        ! Invalid responseID
!        return -1
!    end if
!
!    ! Return success
!    return 0
!
!end subroutine UpdateParameter
    
    
    
    !------------------------------------------------------------------------
    subroutine initialize(initStress, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, MmCur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m)
    ! Initialize PM4Sand Material
    
    ! input variables 
    real(REAL_TYPE), dimension(3), intent(in) :: initStress  ! parameter
    real(REAL_TYPE), dimension(3), intent(inout) :: mSigma_n

    real(REAL_TYPE), intent(inout), dimension(3) :: rrr
    
    real(REAL_TYPE), intent(inout) :: mMb
    real(REAL_TYPE), intent(inout) :: mMd
    real(REAL_TYPE), intent(inout) :: mKp
    
    real(REAL_TYPE) :: ksi !intent(inout)
    real(REAL_TYPE) :: Mcut !intent(inout)
    real(REAL_TYPE) :: Mfin !intent(inout)
    
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_n

    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min_n
    
    real(REAL_TYPE), intent(inout) :: mzpeak
    real(REAL_TYPE), intent(inout) :: mzcum
    real(REAL_TYPE), intent(inout) :: mpzp
    real(REAL_TYPE), intent(inout) :: mzxp
    real(REAL_TYPE), intent(inout) :: m_z_max     ! parameter 
    real(REAL_TYPE), intent(inout) :: m_Ado
    
    
    real(REAL_TYPE), intent(inout) :: mK
    real(REAL_TYPE), intent(inout) :: mG
    real(REAL_TYPE), intent(inout) :: mMcur

    logical, intent(inout) :: m_pzpFlag
    
    real(REAL_TYPE), intent(in) :: m_P_atm
    real(REAL_TYPE), intent(in) :: m_Dr
    real(REAL_TYPE), intent(in) :: m_Mc
    real(REAL_TYPE), intent(in) :: m_nd
    real(REAL_TYPE), intent(in) :: m_nb
    real(REAL_TYPE), intent(in) :: m_m

    
    real(REAL_TYPE) :: m_Pmin
    real(REAL_TYPE) :: m_Pmin2
    
    real(REAL_TYPE), dimension(3, 3) :: CC
    real(REAL_TYPE), dimension(3,3) :: mCe
    real(REAL_TYPE), dimension(3,3) :: mCep
    real(REAL_TYPE), dimension(3,3) :: mCep_Consistent
    real(REAL_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE), dimension(3) :: mSigma_b
    

    real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3,3) :: GetStiffness_result
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: GetKsi_result
    real(REAL_TYPE) :: p0
    
    ! initialize identity matrix 
    mI1(1) = 1 
    mI1(2) = 1
    mI1(3) = 0
    
    
    ! find p0 
    call GetTrace(initStress, GetTrace_result)
    p0 = 0.5 * p0 ! mean effective stress 
    
    ! minimum p' 
    m_Pmin = max(p0/200, m_P_atm/200)
    
    ! p_min for stress
    m_Pmin2 = m_Pmin * 10
    
    
    ! check if p0 is less than m_Pmin 
    if (p0 < m_Pmin) then 
        
        ! initial p is small, set p to p_min and store the difference(mSigmab), the difference
		! will be added to the stress returned to element
        
        mSigma_n = m_Pmin * mI1
        mSigma_b = initStress - mSigma_n
        p0 = m_Pmin
        mAlpha = 0
        mAlpha_n = 0
        
    else 

        mSigma_n = initStress
        mSigma_b = 0.0
        call GetDevPart(initStress, GetDevPart_result)
        mAlpha_n = GetDevPart_result/p0
    
    end if 
    
    
    call GetKsi(m_Dr, p0, GetKsi_result)
    ksi = GetKsi_result
    
    if (m_z_max < 0) then 
        m_z_max = min(0.7*exp(-6.1*ksi), 20.0)
    end if
    
    
    ! bounding and dilatancy surface variations depends on whether it is dense of loose
    if (ksi < 0) then 
        
        ! dense of critical
        mMb = m_Mc * exp(-1.0 * m_nb * ksi)
        mMd = m_Mc * exp(m_nd * ksi)
        
        if (m_Ado < 0) then 
            
            if (mMb > 2.0) then 
                ! Warning, Mb is larger than 2, using Ado = 1.5.
                m_Ado = 1.5
            else 
                m_Ado = 2.5 * (asin(mMb/2.0) - asin(m_Mc/2.0)) / (mMb - mMd)
            end if 
            
        end if 
        
    else 
        
        
        ! loose of critical 
        mMb = m_Mc * exp(-1.0 * (m_nb/4.0) * ksi)
        mMd = m_Mc * exp(       (m_nd*4.0) * ksi)
        
        if (m_Ado < 0) then 
            m_Ado = 1.24
        end if 
        
    end if 
    
    
    ! check if initial stresses are inside bounding and dilatancy surface
    Mcut = max(mMb, mMd)
    
    call GetDevPart(mSigma_n, GetDevPart_result)
    call GetNorm_Contr(GetDevPart_result, GetNorm_Contr_result)
    Mfin = sqrt(2.0) * GetNorm_Contr_result
    Mfin = Mfin / p0
    
    if (Mfin > Mcut) then
        
        rrr = (mSigma_n - (p0*mI1)) * (1/p0) * (Mcut/Mfin)
        ! initial stress outside bounding/dilatancy surface, scale shear stress and store the difference(mSigma_b),
		! the difference will be added to the stress returned to element to maintain global equilibrium
        
        mSigma_n = (p0*mI1) + (rrr*p0)
        mSigma_b = initStress - mSigma_n
        mAlpha_n = rrr * (Mcut - m_m) * (1/Mcut)
        
    end if
    
    mzcum = 0.0
    
    call GetElasticModuli_(mSigma_n, mK, mG, mMcur, mzcum)
    
    call GetStiffness(mK, mG, GetStiffness_result)
    mCe = GetStiffness_result
    mCep = GetStiffness_result
    mCep_consistent = GetStiffness_result
    
    mKp = 100 * mG
    mAlpha = mAlpha_n
    
    mAlpha_in = 0.0
    mAlpha_in_n = 0.0
    mAlpha_in_p = 0.0
    mAlpha_in_p_n = 0.0
    
    mAlpha_in_true = mAlpha_n
    mAlpha_in_true_n = mAlpha_n
	mAlpha_in_max = mAlpha_n
	mAlpha_in_max_n = mAlpha_n
	mAlpha_in_min = mAlpha_n
	mAlpha_in_min_n = mAlpha_n
    
    mFabric = 0.0
    mFabric_in = 0.0
    mFabric_in_n = 0.0
    mFabric_n = 0.0
    
    !mTracker = 0.0 ! --> not sure what this is doing 
    
    mzpeak = m_z_max / 100000.0
    mpzp = max(p0, m_Pmin) / 100.0
    mzxp = 0.0
    m_pzpFlag = .true.
    
    
    
    end subroutine initialize
    
    
    
    !------------------------------------------------------------------------
    
    subroutine setTrialStrain(strain_from_element, mEpsilon_n, &
        mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
        mEpsilonE, mEpsilonE_n)
    ! Input variables
    real(REAL_TYPE), dimension(3), intent(in) :: strain_from_element 
    ! --> Epsilon_xx, Epsilon_yy, Epsilon_xy
    
    ! variables that need to be an input?
    real(REAL_TYPE), dimension(3) :: mEpsilon
    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilon_n
    
    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilonE, mEpsilonE_n
    
    
    ! variables 
    real(REAL_TYPE), dimension(3), intent(inout) :: mSigma
    real(REAL_TYPE), dimension(3), intent(inout) :: mSigma_n
    
    real(REAL_TYPE), dimension(3), intent(inout) :: rrr
    real(REAL_TYPE), intent(inout) :: mMb
    real(REAL_TYPE), intent(inout) :: mMd
    real(REAL_TYPE), intent(inout) :: mKp
    
    
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_n

    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min_n
    
    real(REAL_TYPE), intent(inout) :: mzpeak
    real(REAL_TYPE), intent(inout) :: mzcum
    real(REAL_TYPE), intent(inout) :: mpzp
    real(REAL_TYPE), intent(inout) :: mzxp
    real(REAL_TYPE), intent(inout) :: m_z_max     ! parameter 
    real(REAL_TYPE), intent(inout) :: m_Ado
    
    
    real(REAL_TYPE), intent(inout) :: mK
    real(REAL_TYPE), intent(inout) :: mG
    real(REAL_TYPE), intent(inout) :: mMcur

    logical, intent(inout) :: m_pzpFlag
    
    real(REAL_TYPE), intent(in) :: m_P_atm
    real(REAL_TYPE), intent(in) :: m_Dr
    real(REAL_TYPE), intent(in) :: m_Mc
    real(REAL_TYPE), intent(in) :: m_nd
    real(REAL_TYPE), intent(in) :: m_nb
    real(REAL_TYPE), intent(in) :: m_m
    
    ! Set trial strain and call integrate subroutine
    mEpsilon = strain_from_element
    ! note that now compression is positive, tension is negative
    ! but I am not sure if this is appropriate
    mEpsilon = mEpsilon * (-1.0) ! -1.0 is for geotechnical sign convention
    
    call integrate(mEpsilon, mEpsilon_n, mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
        mEpsilonE, mEpsilonE_n)
    
    !integrate(mEpsilon, mEpsilon_n, mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
    !    mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
    !    mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
    !    mEpsilonE, mEpsilonE_n)
    
    end subroutine setTrialStrain
    
    
    
    
    
    !------------------------------------------------------------------------

    subroutine GetState(State)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(16), intent(out) :: State

    ! variables
    real(REAL_TYPE), dimension(3) :: mEpsilonE
    real(REAL_TYPE), dimension(3) :: mAlpha_n
    real(REAL_TYPE), dimension(3) :: mFabric_n
    real(REAL_TYPE), dimension(3) :: mAlpha_in_n
    real(REAL_TYPE) :: mVoidRatio
    real(REAL_TYPE) :: mDGamma_n
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mKp
    
    ! Assemble state parameters into the result vector
    State(1:3) = mEpsilonE
    State(4:6) = mAlpha_n
    State(7:9) = mFabric_n
    State(10:12) = mAlpha_in_n
    State(13) = mVoidRatio
    State(14) = mDGamma_n
    State(15) = mG
    State(16) = mKp


    end subroutine GetState

    !------------------------------------------------------------------------

subroutine GetAlpha(Alpha)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha
    
    ! variables 
    real(REAL_TYPE) :: mAlpha_n

    ! Return alpha tensor
    Alpha = mAlpha_n

end subroutine GetAlpha

    !------------------------------------------------------------------------

    subroutine GetFabric(Fabric)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Fabric
    
    real(REAL_TYPE) :: mFabric_n 

    ! Return fabric tensor
    Fabric = mFabric_n

    end subroutine GetFabric

    !------------------------------------------------------------------------

    subroutine GetAlpha_in(Alpha_in)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in

    real(REAL_TYPE) :: mAlpha_in_n
    
    ! Return alpha_in tensor
    Alpha_in = mAlpha_in_n

    end subroutine GetAlpha_in

    !------------------------------------------------------------------------

    subroutine GetTracker(Tracker)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Tracker

    ! variables 
    real(REAL_TYPE) :: mTracker
    
    ! Return tracker vector
    Tracker = mTracker

    end subroutine GetTracker

    !------------------------------------------------------------------------


    subroutine GetKp(Kp)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: Kp

    ! variables 
    real(REAL_TYPE) :: mKp
    
    ! Return Kp
    Kp = mKp

    end subroutine GetKp

    !------------------------------------------------------------------------


    subroutine GetG(G)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: G

    ! variables 
    real(REAL_TYPE):: mG
    
    ! Return shear modulus
    G = mG

    end subroutine GetG

    !------------------------------------------------------------------------


    subroutine GetAlpha_in_p(Alpha_in_p)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_p

    real(REAL_TYPE), dimension(3) :: mAlpha_in_p_n
    
    ! Return previous alpha_in tensor
    Alpha_in_p = mAlpha_in_p_n


    end subroutine GetAlpha_in_p

    !------------------------------------------------------------------------


    subroutine GetDGamma(DGamma)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: DGamma
    
    real(REAL_TYPE) :: mDGamma_n

    ! Return previous L
    DGamma = mDGamma_n

    end subroutine GetDGamma
    
    !------------------------------------------------------------------------

    subroutine GetTangent(TangType, Ce, Cep, Cep_Consistent, Tangent)
    implicit none
    ! Inputs
    integer, intent(in) :: TangType
    real(REAL_TYPE), dimension(3,3), intent(in) :: Ce, Cep, Cep_Consistent
    ! Output
    real(REAL_TYPE), dimension(3,3), intent(out) :: Tangent

    ! Determine which tangent to return based on TangType
    if (TangType == 0) then
        Tangent = Ce
    else if (TangType == 1) then
        Tangent = Cep
    else
        Tangent = Cep_Consistent
    endif


    end subroutine GetTangent


    !------------------------------------------------------------------------

    subroutine getInitialTangent(initialTangent)
    ! Input/output variables
    real(REAL_TYPE), dimension(6,6), intent(out) :: initialTangent
    
    ! variables 
    real(REAL_TYPE), dimension(6,6) :: mCe
    
    ! Return initial tangent matrix
    initialTangent = mCe
    
    end subroutine getInitialTangent

    

    !------------------------------------------------------------------------

    subroutine getStress(stress)
    ! Input/output variables
    real(REAL_TYPE), dimension(3), intent(out) :: stress
    
    ! variables 
    real(REAL_TYPE), dimension(3) :: mSigma
    real(REAL_TYPE), dimension(3) :: mSigma_b
    
    ! Calculate and return stress
    stress = mSigma + mSigma_b
    stress = stress * (-1.0)  ! -1.0 is for geotechnical sign convention
    
    end subroutine getStress

    
        
    !------------------------------------------------------------------------

    subroutine getStrain(strain)
    ! Input/output variables
    real(REAL_TYPE), dimension(3), intent(out) :: strain
    
    
    ! variables 
    real(REAL_TYPE), dimension(3) :: mEpsilon
    
    ! Calculate and return strain
    strain = mEpsilon
    strain = strain * (-1.0)  ! -1.0 is for geotechnical sign convention
    
    end subroutine getStrain
    
    !------------------------------------------------------------------------

    
    !subroutine GetStress(Stress, Sigma_r)
    !implicit none
    !! Inputs
    !real(REAL_TYPE), dimension(6), intent(in) :: Stress
    !! Output
    !real(REAL_TYPE), dimension(6), intent(out) :: Sigma_r
    !
    !! Compute stress
    !Sigma_r = Stress + mSigma_b
    !Sigma_r = Sigma_r * (-1.0)
    !
    !
    !end subroutine GetStress
    
    
    
    !-----------------------------------------------------------------------
    subroutine GetElasticStrain(Epsilon, EpsilonE_r)
    implicit none
    ! Inputs
    real(REAL_TYPE), dimension(6), intent(in) :: Epsilon
    ! Output
    real(REAL_TYPE), dimension(6), intent(out) :: EpsilonE_r

    ! Compute elastic strain
    EpsilonE_r = Epsilon * (-1.0)

    end subroutine GetElasticStrain
    !subroutine getElasticStrain(elasticStrain)
    !! Input/output variables
    !real*8, dimension(3), intent(out) :: elasticStrain
    !
    !! Calculate elastic strain
    !elasticStrain = mEpsilon
    !elasticStrain = elasticStrain * (-1.0)
    !
    !end subroutine getElasticStrain
    
    
    
    
    !------------------------------------------------------------------------------------
    !// ---------------------------------------------------------------------------------
    !/*************************************************************/
    !// Plastic Integrator
    !/*************************************************************/
    
    subroutine integrate(mEpsilon, mEpsilon_n, mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
        mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
        mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
        mEpsilonE, mEpsilonE_n)
    implicit none
    ! Local variables
    real(REAL_TYPE) :: pp, zxpTemp
    real(REAL_TYPE), dimension(3) :: rrr
    integer(INTEGER_TYPE) :: ii
    real(REAL_TYPE), dimension(3) :: n_tr, tmp0, tmp1, mAlpha_mAlpha_in_true
    
    ! variables --> interconnection
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_n

    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min
    real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min_n
    
    real(REAL_TYPE), dimension(3), intent(inout) :: mSigma_n    
    real(REAL_TYPE), dimension(3), intent(inout) :: mSigma

    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilon
    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilon_n
    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilonE
    real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilonE_n
    real(REAL_TYPE), dimension(3,3) :: mCe
    
    
    real(REAL_TYPE), intent(inout) :: mMb
    real(REAL_TYPE), intent(inout) :: mMd
    real(REAL_TYPE), intent(inout) :: mKp
    
    real(REAL_TYPE), intent(inout) :: mzpeak
    real(REAL_TYPE), intent(inout) :: mzcum
    real(REAL_TYPE), intent(inout) :: mpzp
    real(REAL_TYPE), intent(inout) :: mzxp
    real(REAL_TYPE), intent(inout) :: m_z_max     ! parameter 
    real(REAL_TYPE), intent(inout) :: m_Ado
    
    
    real(REAL_TYPE), intent(inout) :: mK
    real(REAL_TYPE), intent(inout) :: mG
    real(REAL_TYPE), intent(inout) :: mMcur
    
    real(REAL_TYPE), intent(in) :: m_P_atm
    real(REAL_TYPE), intent(in) :: m_Dr
    real(REAL_TYPE), intent(in) :: m_Mc
    real(REAL_TYPE), intent(in) :: m_nd
    real(REAL_TYPE), intent(in) :: m_nb
    real(REAL_TYPE), intent(in) :: m_m
    
    
    ! Local variables 
    real(REAL_TYPE), dimension(3) :: GetNormalToYield_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    integer(INTEGER_TYPE) :: me2p
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: m_pmin
    real(REAL_TYPE) :: GetNorm_Contr_result
    !real(REAL_TYPE), intent(inout) :: mzxp
    !real(REAL_TYPE), intent(inout) :: mpzp
    logical, intent(inout) :: m_pzpflag !--> this might need to be logical !integer(INTEGER_TYPE)
    !real(REAL_TYPE), dimension(3) :: mEpsilon_n !--> dimension(3)?
    !real(REAL_TYPE), dimension(3) :: mEpsilon
    real(REAL_TYPE) :: mVoidratio
    !real(REAL_TYPE), intent(inout) :: mG
    !real(REAL_TYPE), intent(inout) :: mK
    real(REAL_TYPE), dimension(3,3) :: mCep
    real(REAL_TYPE), dimension(3,3) :: mCep_consistent
    real(REAL_TYPE) :: mdGamma
    real(REAL_TYPE) :: m_e_init
    
    
    
    ! Assignments, new = old
    mAlpha = mAlpha_n
    mAlpha_in = mAlpha_in_n
    mAlpha_in_true = mAlpha_in_true_n
    mAlpha_in_p = mAlpha_in_p_n
    mAlpha_in_max = mAlpha_in_max_n
    mAlpha_in_min = mAlpha_in_min_n
    mFabric = mFabric_n
    mFabric_in = mFabric_in_n

    ! Calculate trial stress
    tmp0 = mSigma_n
    tmp1 = mEpsilon
    tmp1 = tmp1 - mEpsilon_n
    tmp0 = tmp0 + matmul(mCe, tmp1)!(mCe * tmp1)
    call GetNormalToYield(tmp0, mAlpha, GetNormalToYield_result)
    n_tr = GetNormalToYield_result!GetNormalToYield(tmp0, mAlpha)

    ! Check loading reversal condition
    ! if ((DoubleDot2_2_Contr(mAlpha - mAlpha_in_true, n_tr) < 0.0) && me2p) {
    mAlpha_mAlpha_in_true = mAlpha - mAlpha_in_true
    call DoubleDot2_2_Contr(mAlpha_mAlpha_in_true, n_tr, DoubleDot2_2_Contr_result)
    if ( (DoubleDot2_2_Contr_result < 0.0) .and. (me2p /= 0) ) then !--> me2p will actually be 1
        mAlpha_in_p = mAlpha_in
        mAlpha_in_true = mAlpha
        mFabric_in = mFabric
        
        ! Update pzp
        call GetTrace(mSigma_n, GetTrace_result)
        pp = 0.5 * GetTrace_result
        if (pp <= m_Pmin) then
            pp = m_Pmin
        endif
        call GetNorm_Contr(mFabric_n, GetNorm_Contr_result)
        zxpTemp = GetNorm_Contr_result * pp
        if ( ((zxpTemp > mzxp) .and. (pp > mpzp)) .or. m_pzpFlag ) then
            mzxp = zxpTemp
            mpzp = pp
            m_pzpFlag = .false.
        endif
        
        ! Track initial back-stress ratio history
        do ii = 1, 3
            if (mAlpha_in(ii) > 0.0) then
                ! Minimum positive value
                mAlpha_in_min(ii) = min(mAlpha_in_min(ii), mAlpha(ii))
            else
                ! Maximum negative value
                mAlpha_in_max(ii) = max(mAlpha_in_max(ii), mAlpha(ii))
            endif
        enddo
        
        ! Update mAlpha_in based on loading direction
        if (mAlpha(3) * mAlpha_in_p(3) > 0) then
            do ii = 1, 3
                if (n_tr(ii) > 0.0) then
                    ! Positive loading direction
                    mAlpha_in(ii) = max(0.0, mAlpha_in_min(ii))
                else
                    ! Negative loading direction
                    mAlpha_in(ii) = min(0.0, mAlpha_in_max(ii))
                endif
            enddo
        else
            mAlpha_in = mAlpha
        endif
    endif

    ! Force elastic response
    if (me2p == 0) then !--> this will never be the case
        call elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha, &
            mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent)
    else
        ! ElastoPlastic response
        ! Explicit schemes
        call explicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in, &
            mAlpha_in_p, mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, &
            mK, mCe, mCep, mCep_Consistent)
    endif

end subroutine integrate
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !------------------------------------------------------------------------
    !/*************************************************************/
    !// Elastic Integrator
    !/*************************************************************/
    
    subroutine elastic_integrator(CurStress, CurStrain, CurElasticStrain, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, NextVoidRatio, &
        GG, KK, aC, aCep, aCep_Consistent)
    implicit none
    ! Input variables
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), dimension(3), intent(in) :: CurStrain
    real(REAL_TYPE), dimension(3), intent(in) :: CurElasticStrain
    real(REAL_TYPE), dimension(3), intent(in) :: NextStrain

    ! Output variables
    real(REAL_TYPE), dimension(3), intent(out) :: NextElasticStrain
    real(REAL_TYPE), dimension(3), intent(out) :: NextStress
    real(REAL_TYPE), dimension(3), intent(out) :: NextAlpha
    real(REAL_TYPE), intent(out) :: NextVoidRatio
    real(REAL_TYPE), intent(inout) :: GG ! I changed this to inout. It was in.
    real(REAL_TYPE), intent(inout) :: KK ! I changed this to inout. It was in.
    real(REAL_TYPE), dimension(3, 3), intent(in) :: aC
    real(REAL_TYPE), dimension(3, 3), intent(inout) :: aCep ! I changed this to inout. It was in.
    real(REAL_TYPE), dimension(3, 3), intent(inout) :: aCep_Consistent ! I changed this to inout. It was in.

    ! Local variables
    real(REAL_TYPE) :: pp
    real(REAL_TYPE), dimension(3) :: dStrain
    
    real(REAL_TYPE) :: GetTrace_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
	
	real(REAL_TYPE) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: m_e_init
    
    real(REAL_TYPE) :: m_Pmin
    
    real(REAL_TYPE), dimension(3) :: DetDevPart_result
    
    !real(REAL_TYPE) :: GetTrace_result

    
    ! Calculate strain increment
    !// calculate elastic response
	!// dStrain = NextStrain - CurStrain;
    dStrain = NextStrain - CurStrain
    ! Calculate the void ratio
    call GetTrace(NextStrain, GetTrace_result)
    NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace_result!sum(NextStrain)
    ! Calculate the next elastic strain
    NextElasticStrain = CurElasticStrain + dStrain
    ! Calculate elastic moduli
    call GetElasticModuli(CurStress, KK, GG)
    ! Set consistent aCep
    aCep_Consistent = aCep
    aCep = aC
    ! Calculate the next stress
    call DoubleDot4_2(aC, dStrain, DoubleDot4_2_result)
    NextStress = CurStress + DoubleDot4_2_result!DoubleDot4_2(aC, dStrain)
    ! Calculate mean effective stress
    call GetTrace(NextStress, GetTrace_result)
    pp = 0.5 * GetTrace_result !sum(NextStress)
    ! Calculate plastic multiplier alpha if mean effective stress is above the minimum threshold
    
    if (pp > m_Pmin) then
        call GetDevPart(NextStress, DetDevPart_result)
        NextAlpha = DetDevPart_result / pp
    end if


    end subroutine elastic_integrator
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Explicit Integrator
    !/*************************************************************/
    
    
    
    subroutine explicit_integrator(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
    NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, mG, mK, aC, aCep, aCep_Consistent)
    
    implicit none
    
    ! Input variables
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), dimension(3), intent(in) :: CurStrain
    real(REAL_TYPE), dimension(3), intent(in) :: CurElasticStrain
    real(REAL_TYPE), dimension(3), intent(in) :: CurAlpha
    real(REAL_TYPE), dimension(3), intent(in) :: CurFabric
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_p
    real(REAL_TYPE), dimension(3), intent(in) :: NextStrain

    ! Output variables
    real(REAL_TYPE), dimension(3), intent(out) :: NextElasticStrain
    real(REAL_TYPE), dimension(3), intent(out) :: NextStress
    real(REAL_TYPE), dimension(3), intent(out) :: NextAlpha
    real(REAL_TYPE), dimension(3), intent(out) :: NextFabric
    real(REAL_TYPE), intent(out) :: NextL
    real(REAL_TYPE), intent(out) :: NextVoidRatio
    real(REAL_TYPE), intent(in) :: mG
    real(REAL_TYPE), intent(in) :: mK
    real(REAL_TYPE), dimension(3, 3), intent(inout) :: aC !in
    real(REAL_TYPE), dimension(3, 3), intent(inout) :: aCep !in
    real(REAL_TYPE), dimension(3, 3), intent(inout) :: aCep_Consistent !in

    ! Local variables
    real(REAL_TYPE) :: elasticRatio, ff, fn, dVolStrain
    real(REAL_TYPE), dimension(3) :: dStrain, dSigma, dDevStrain, nn, tmp, dElasStrain

    real(REAL_TYPE) :: m_e_init
    !real(REAL_TYPE) :: mG
    !real(REAL_TYPE) :: mK
    real(REAL_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE) :: GetTrace_result
    
    real(REAL_TYPE) :: mTolF
    
    real(REAL_TYPE) :: aa
    
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    
    real(REAL_TYPE) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: a0
    real(REAL_TYPE) :: a1
    
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 0
    
    ! Calculate the next void ratio using the strain
    call GetTrace(NextStrain, GetTrace_result) !--> trace to find volumetric strain
    NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace_result!sum(NextStrain)
    ! Calculate strain increment
    dStrain = NextStrain - CurStrain
    ! Calculate the elastic strain increment
    NextElasticStrain = CurElasticStrain + dStrain
    ! Calculate volumetric strain and deviatoric strain
    call GetTrace(dStrain, GetTrace_result)
    dVolStrain = GetTrace_result!sum(dStrain)
    dDevStrain = dStrain - dVolStrain / 3.0

    ! Calculate the stress increment
    tmp = 2 * mG * dDevStrain
    dSigma = mK * dVolStrain * mI1 + tmp
    ! Update the stress
    NextStress = CurStress + dSigma

    ! Calculate some parameters for plasticity
    !ff = GetF(NextStress, CurAlpha)
    call GetFYieldFunction(NextStress, CurAlpha, ff)
    !fn = GetF(CurStress, CurAlpha)
    call GetFYieldFunction(CurStress, CurAlpha, fn)
    !nn = GetNormalToYield(NextStress, CurAlpha)
    call GetNormalToYield(NextStress, CurAlpha, nn)
    ! Perform plasticity calculations
    if (ff <= mTolF) then
        ! Pure elastic loading/unloading
        NextAlpha = CurAlpha
        NextFabric = CurFabric
        NextL = 0
        aCep_Consistent = aCep
        aCep = aC
    elseif (fn < -mTolF) then
        ! Transition from elastic to plastic
        a0 = 0.0
        a1 = 1.0
        call IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1, aa)
        elasticRatio = aa
        dElasStrain = elasticRatio * (NextStrain - CurStrain)
        call DoubleDot4_2(aC, dElasStrain, DoubleDot4_2_result)
        !dSigma = DoubleDot4_2(aC, dElasStrain)
        dSigma = DoubleDot4_2_result!DoubleDot4_2(aC, dElasStrain)
        ! Temporary commented exp_int
        !call exp_int(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
        !    CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
        !    NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    elseif (abs(fn) < mTolF) then
        ! Pure plastic step or elastic unloading followed by plastic loading
        call GetNormalToYield(CurStress, CurAlpha, nn)
        call DoubleDot2_2_Contr(nn, dSigma, DoubleDot2_2_Contr_result)
        call GetNorm_Contr(dSigma, GetNorm_Contr_result)
        !if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha), dSigma) / (max(1.0, GetNorm_Contr(dSigma))) > (-sqrt(mTolF))) then
        if (DoubleDot2_2_Contr_result / (max(1.0, GetNorm_Contr_result)) > (-sqrt(mTolF))) then
            ! Pure plastic step
            ! Temporary commented exp_int
            !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
            !    NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
        else
            ! Elastic unloading followed by plastic loading
            call IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha, aa)
            elasticRatio = aa
            dElasStrain = elasticRatio * (NextStrain - CurStrain)
            call DoubleDot4_2(aC, dElasStrain, DoubleDot4_2_result)
            dSigma = DoubleDot4_2_result!DoubleDot4_2(aC, dElasStrain)
            ! Temporary commented exp_int
            !call exp_int(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
            !    CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
            !    NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
        end if
    else
        call GetFYieldFunction(CurStress, CurAlpha,ff)
        ! Illegal stress state
        print *, "PM4Sand: Encountered an illegal stress state! Tag: "!, !getTag()
        print *, "            f = ", ff
        ! Temporary commented exp_int
        !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
        !    NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    end if

end subroutine explicit_integrator
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------------
    !// ----------------------------------------------------------------------
    !/*************************************************************/
    !// Forward-Euler Integrator
    !/*************************************************************/
    subroutine ForwardEuler(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
                                 NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    
    implicit none
    
    ! Inputs --> I changed it from intent(in) to intent(inout)
    real(REAL_TYPE), intent(inout) :: CurStress(3), CurStrain(3), CurElasticStrain(3), CurAlpha(3), CurFabric(3), alpha_in(3), alpha_in_p(3), NextStrain(3)
    
    ! Outputs
    real(REAL_TYPE), intent(out) :: NextElasticStrain(3), NextStress(3), NextAlpha(3), NextFabric(3), NextL, NextVoidRatio, GG, KK
    real(REAL_TYPE), intent(out) :: aC(3,3), aCep(3,3), aCep_Consistent(3,3)
    
    ! Local variables
    real(REAL_TYPE) :: CurVoidRatio, CurDr, Cka, hh, pp, dVolStrain, DD, AlphaAlphaBDotN
    real(REAL_TYPE) :: nn(3), RR(3), alphaD(3), dPStrain(3), bb(3), dDevStrain(3), rrr(3), dStrain(3)
    real(REAL_TYPE) :: dSigma(3), dAlpha(3), dFabric(3)
    
    real(REAL_TYPE) :: temp4
    real(REAL_TYPE) :: GetTrace_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
	
	real(REAL_TYPE) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: mMcur
    real(REAL_TYPE) :: mZcum
    real(REAL_TYPE) :: m_e_init
    real(REAL_TYPE) :: m_emax
    real(REAL_TYPE) :: m_emin
    real(REAL_TYPE) :: m_pmin
    !real(REAL_TYPE) :: pp
    real(REAL_TYPE), dimension(3) :: mFabric_in
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mZpeak
    real(REAL_TYPE) :: mpzp
    real(REAL_TYPE) :: mKp
    real(REAL_TYPE) :: mDGamma
    integer(INTEGER_TYPE) :: DebugFlag
    
    integer(INTEGER_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE) :: m_z_max 
    real(REAL_TYPE) :: m_cz
    real(REAL_TYPE) :: two3
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    two3 = 0.666666666666666667
    
    ! Get elastic moduli
    call GetElasticModuli_(NextStress, KK, GG, mMcur, mzcum)
    ! Calculate current void ratio
    call GetTrace(CurStrain, GetTrace_result)
    CurVoidRatio = m_e_init - (1.0 + m_e_init) * GetTrace_result !sum(CurStrain) 
    ! m_e_init needs to be global to the entire PM4Sand module 
    ! Calculate current Dr
    CurDr = (m_emax - CurVoidRatio) / (m_emax - m_emin) 
    ! m_emax and m_emin needs to be global the entire PM4Sand module 
    ! Calculate mean effective stress
    call GetTrace(CurStress, GetTrace_result)
    pp = 0.5 * GetTrace_result!sum(CurStress)
    if (pp < m_Pmin) pp = m_Pmin ! Apply tension cutoff
    ! Calculate void ratio from NextStrain
    call GetTrace(NextStrain, GetTrace_result)
    NextVoidRatio = m_e_init - (1.0 + m_e_init) * GetTrace_result!sum(NextStrain)
    ! Calculate strain increment
    !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
    dStrain = NextStrain - CurStrain
    ! Calculate NextElasticStrain
    NextElasticStrain = CurElasticStrain + dStrain
    !// using NextStress instead of CurStress to get correct n
	!// 1: calculate state parameters and initialize them
    ! Calculate state parameters
    call GetStateDependent(NextStress, CurAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum, &
                           mzpeak, mpzp, mMcur, CurDr, nn, DD, RR, mKp, alphaD, Cka, hh, bb, AlphaAlphaBDotN)
    ! Calculate volumetric strain increment
    call GetTrace(dStrain, GetTrace_result)
    !// dVolStrain = GetTrace(NextStrain - CurStrain);
    dVolStrain = GetTrace_result!sum(dStrain)
    ! Calculate deviatoric strain increment
	!// dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
    dDevStrain = dStrain - dVolStrain / 3.0
    ! Calculate stress ratio
    call GetDevPart(NextStress, GetDevPart_result)
    rrr = GetDevPart_result / pp
    ! Calculate denominator term
    !double precision :: temp4
    call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
    temp4 = mKp + 2.0 * GG - KK * DD * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(nn, rrr)
    !// if (temp4 < 0.0) {
	!// 	mKp = -0.5 * (2 * G - K* D *DoubleDot2_2_Contr(n, r));
	!// 	temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
	!// 	h = 1.5 * mKp / (p * AlphaAlphaBDotN);
	!// }
    
    ! Check denominator
    if (abs(temp4) < small) then
    
        ! Neutral loading
        dSigma = 0.0
        dAlpha = 0.0
        dFabric = 0.0
        ! dPStrain = dDevStrain + dVolStrain * mI1
        dPStrain = dStrain

    else
    
        ! Calculate L
        call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
        call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
        !NextL = (2.0 * GG * DoubleDot2_2_Mixed(nn, dDevStrain) - DoubleDot2_2_Contr(nn, rrr) * KK * dVolStrain) / temp4 ! Equation 31
        NextL = (2.0 * GG * DoubleDot2_2_Mixed_result - DoubleDot2_2_Contr_result * KK * dVolStrain) / temp4 ! Equation 31

        ! Set mDGamma to NextL
        mDGamma = NextL
        if (NextL < 0.0) then
            ! If NextL is negative
            if (debugFlag) then
                write(*,*) "NextL is smaller than 0"
                write(*,*) "NextL = ", NextL
            end if
        
            call ToContraviant(dDevStrain, ToContraviant_result)
            dSigma = 2.0 * GG * ToContraviant_result + KK * dVolStrain * mI1
            dAlpha = 0.0
            dFabric = 0.0
            dPStrain = 0.0
    
        else
            ! If NextL is non-negative
            ! Calculate dSigma
            call ToContraviant(dDevStrain, ToContraviant_result)
            call Macauley(NextL, Macauley_result)
            !dSigma = 2.0 * GG * ToContraviant(dDevStrain) + KK * dVolStrain * mI1 - Macauley(NextL) * &
            !     (2.0 * GG * nn + KK * DD * mI1)
            dSigma = 2.0 * GG * ToContraviant_result + KK * dVolStrain * mI1 - Macauley_result * &
                 (2.0 * GG * nn + KK * DD * mI1)
            ! Update fabric
            call DoubleDot2_2_Contr(alphaD - CurAlpha, nn, DoubleDot2_2_Contr_result)
            if (DoubleDot2_2_Contr_result < 0.0) then
                ! Update fabric according to Equation 57
                dFabric = nn
                dFabric = m_z_max * nn + CurFabric
                !dFabric = (-1.0 * (m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0))) * Macauley(NextL) * &
                !       MacauleyIndex(-DD) * dFabric)
                call Macauley(mzcum / 2.0 / m_z_max - 1.0, Macauley_result)
                call Macauley(NextL, Macauley_result2)
                call MacauleyIndex(-DD, MacauleyIndex_result)
                !//  dz in Equation 57 
                dFabric = (-1.0 * (m_cz / (1 + Macauley_result)) * Macauley_result2 * &
                       MacauleyIndex_result * dFabric)
            end if
            ! Update alpha
            !// dAlpha = two3 * NextL * h * b;
			!// dPStrain = NextL * mIIco * R;
            call ToCovariant(RR, ToCovariant_result)
            dPStrain = ToCovariant_result * NextL
            dAlpha = bb * (two3 * NextL * hh)
    
        end if
    end if


    ! Update NextFabric, NextElasticStrain, NextStress, and NextAlpha
    !// NextFabric = CurFabric + dFabric;
    NextFabric = CurFabric + dFabric
    !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain
    !// NextStress = CurStress + dSigma;
    NextStress = CurStress + dSigma
    !// NextAlpha = CurAlpha + dAlpha;
    NextAlpha = CurAlpha + dAlpha

    ! Perform stress correction
    call Stress_Correction_(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio)



    end subroutine ForwardEuler

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !--------------------------------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Integrator Constraining Maximum Strain Increment
    !/*************************************************************/
    
    subroutine MaxStrainInc(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
        GG, KK, aC, aCep, aCep_Consistent)
    
    implicit none
    
    ! Inputs
    real(REAL_TYPE), intent(in) :: CurStress(3), CurStrain(3), CurElasticStrain(3), CurAlpha(3), CurFabric(3), alpha_in(3), alpha_in_p(3), NextStrain(3)
    ! Outputs
    real(REAL_TYPE), intent(out) :: NextElasticStrain(3), NextStress(3), NextAlpha(3), NextFabric(3), NextL, NextVoidRatio, GG, KK
    real(REAL_TYPE), intent(out) :: aC(3,3), aCep(3,3), aCep_Consistent(3,3)
    
    ! Function pointer to the integration scheme
    ! Temporary commented exp_int
    !procedure(PM4Sand_IntegrationScheme), pointer :: exp_int !--> more work
    
    ! Local variables
    real(REAL_TYPE) :: StrainInc(3), maxInc
    integer(INTEGER_TYPE) :: ii, numSteps
    real(REAL_TYPE) :: cStress(3), cStrain(3), cAlpha(3), cFabric(3), cAlpha_in(3), cAlpha_in_p(3), cEStrain(3)
    real(REAL_TYPE) :: nStrain(3)
    real(REAL_TYPE), dimension(3,3) :: nCe, nCep, nCepC
    
    real(REAL_TYPE) :: maxStrainIncValue
    
    maxStrainIncValue = 1e-5 !maxStrainInc

    ! Assign function pointer based on integration scheme
    ! Temporary commented exp_int
    !select case (mScheme)
    !    case (INT_MAXSTR_FE)
    !        exp_int => ForwardEuler !PM4Sand_
    !    case (INT_MAXSTR_ME)
    !        exp_int => ModifiedEuler !PM4Sand_
    !    case default
    !        exp_int => ModifiedEuler !PM4Sand_
    !end select
    
    ! Compute maximum strain increment
    StrainInc = NextStrain - CurStrain
    maxInc = StrainInc(1)
    do ii = 2, 3
        if (abs(StrainInc(ii)) > abs(maxInc)) maxInc = StrainInc(ii)
    end do
    
    ! Apply maximum strain increment constraint
    if (abs(maxInc) > maxStrainIncValue) then
        numSteps = floor(abs(maxInc) / maxStrainIncValue) + 1 
        ! floor(A) returns the greatest integer less than or equal to A
        StrainInc = (NextStrain - CurStrain) / numSteps !real(, kind=8)

        ! Initialize temporary variables
        cStress = CurStress
        cStrain = CurStrain
        cAlpha = CurAlpha
        cFabric = CurFabric
        cAlpha_in = alpha_in
        cAlpha_in_p = alpha_in_p
        cEStrain = CurElasticStrain

        do ii = 1, numSteps
            nStrain = cStrain + StrainInc

            ! Call integration scheme subroutine --> pointer
            ! Temporary commented exp_int
            !call exp_int(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, cAlpha_in_p, &
            !             nStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
            !             NextL, NextVoidRatio, GG, KK, nCe, nCep, nCepC)

            ! Update temporary variables
            cStress = NextStress
            cStrain = nStrain
            cEStrain = NextElasticStrain
            cAlpha = NextAlpha
            cFabric = NextFabric
        end do
    else
        ! Call integration scheme subroutine without applying maximum strain increment constraint
        ! Temporary commented exp_int
        !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, &
        !             NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
        !             NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    end if


    end subroutine MaxStrainInc !PM4Sand_
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !---------------------------------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Modified-Euler Integrator
    !/*************************************************************/
    
    subroutine ModifiedEuler(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, &
        alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
        GG, KK, aC, aCep, aCep_Consistent)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain
    real(REAL_TYPE), intent(out) :: NextElasticStrain(3), NextStress(3), NextAlpha(3), NextFabric(3)
    real(REAL_TYPE), intent(out) :: NextL, NextVoidRatio, GG, KK
    real(REAL_TYPE), dimension(3, 3), intent(in) :: aC, aCep, aCep_Consistent
    ! Local variables
    real(REAL_TYPE) :: NextDr, dVolStrain, pp, Cka, temp4, curStepError, qq, stressNorm, hh, DD, AlphaAlphaBDotN
    real(REAL_TYPE) :: TT, dT, dT_min, TolE
    
    real(REAL_TYPE) :: small != 1.0e-10
    real(REAL_TYPE) :: tmp0(3), tmp1(3), tmp2(3), nn(3), R1(3), R2(3), alphaD(3), dDevStrain(3), rr(3), bb(3)
    real(REAL_TYPE) :: dSigma1(3), dSigma2(3), dAlpha1(3), dAlpha2(3), dFabric1(3), dFabric2(3), dPStrain1(3), dPStrain2(3)
    real(REAL_TYPE) :: mDGamma, m_Pmin, m_e_init, m_emax, m_emin, m_z_max, m_cz, mKp, mMcur, mzcum, mzpeak, mpzp
    real(REAL_TYPE), dimension(3) :: mFabric_in
    real(REAL_TYPE) :: two3, mIIco ! These should be initialized with actual values

    
    ! Local variables for use with subroutines below 
    real(REAL_TYPE) :: GetTrace_result
	real(REAL_TYPE) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
    
    real(REAL_TYPE), dimension(3) :: nStress
    real(REAL_TYPE), dimension(3) :: nFabric
    real(REAL_TYPE), dimension(3) :: nAlpha
    
    real(REAL_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE), dimension(3) :: alphaD_NextAlpha
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    
    TT = 0.0 
    dT = 1.0 
    dT_min = 1e-4 
    TolE = 1e-5
    
    small = 1.0e-10
    
    two3 = 2.0/3.0
    mIIco = 1.0
    
    ! Initialize NextElasticStrain, NextStress, NextAlpha, NextFabric
    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
    NextStress = CurStress
    NextAlpha = CurAlpha
    NextFabric = CurFabric
    
    !this->GetElasticModuli(NextStress, K, G, mMcur, mzcum); 
    call GetElasticModuli_(NextStress, KK, GG, mMcur, mzcum)

    ! Calculate p
    GetTrace_result = 0.0
    call GetTrace(CurStress, GetTrace_result)
    pp = 0.5 * GetTrace_result

    ! Check if p < m_Pmin / 5.0 and adjust NextStress if needed
    if (pp < m_Pmin / 5.0) then
        !if (debugFlag) then
        !    write(*,*) "Tag = ", this%getTag(), " : p < pmin / 5, should not happen"
        !end if
        GetDevPart_result = 0.0
        call GetDevPart(NextStress, GetDevPart_result)
        NextStress = GetDevPart_result + m_Pmin / 5.0 * mI1
    end if

    ! Main loop
    do while (TT < 1.0)
        ! Calculate NextVoidRatio and NextDr
        tmp0 = NextStrain - CurStrain
        tmp0 = tmp0 * TT
        tmp0 = CurStrain + tmp0
        GetTrace_result = 0.0
        call GetTrace(tmp0, GetTrace_result)
        NextVoidRatio = m_e_init - (1.0 + m_e_init) * GetTrace_result != m_e_init - (1 + m_e_init) * GetTrace(CurStrain + T*(NextStrain - CurStrain));
        NextDr = (m_emax - NextVoidRatio) / (m_emax - m_emin)

        ! Calculate dVolStrain and dDevStrain
        tmp0 = NextStrain - CurStrain
        GetTrace_result = 0.0
        call GetTrace(tmp0, GetTrace_result)
        dVolStrain = dT * GetTrace_result
        
        dDevStrain = (mI1 * (-dVolStrain / 3.0)) + ( (NextStrain - CurStrain) * dT)

        GetTrace_result = 0.0
        call GetTrace(NextStress, GetTrace_result)
        pp = 0.5 * GetTrace_result

        
        ! Calculate Delta 1
        ! Implement the calculation of dSigma1, dAlpha1, dFabric1, and dPStrain1
        call GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, NextFabric, mFabric_in, GG, mzcum, &
                            mzpeak, mpzp, mMcur, NextDr, nn, DD, R1, mKp, alphaD, Cka, hh, bb, AlphaAlphaBDotN)
        ! --> above subroutine is temporary 
        
        GetDevPart_result = 0.0
        call GetDevPart(NextStress, GetDevPart_result)
        rr = GetDevPart_result / pp

        DoubleDot2_2_Contr_result = 0.0
        call DoubleDot2_2_Contr(nn, rr, DoubleDot2_2_Contr_result)
        temp4 = mKp + 2 * GG - KK * DD * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(n, r)

        if (abs(temp4) < small) then
    
            ! neutral loading
            dSigma1 = 0.0
            dAlpha1 = 0.0
            dFabric1 = 0.0
            dPStrain1 = tmp0

        else
            call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
            call DoubleDot2_2_Contr(nn, rr, DoubleDot2_2_Contr_result)


            !NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4
            NextL = (2 * GG * DoubleDot2_2_Mixed_result - DoubleDot2_2_Contr_result * KK * dVolStrain) / temp4
            
            if (NextL < 0) then
                !if (debugFlag) then
                !    write(*,*) "1 NextL is smaller than 0"
                !    write(*,*) "NextL = ", NextL
                !end if
                call ToContraviant(dDevStrain, ToContraviant_result)
        
                dSigma1 = 2 * GG * ToContraviant_result + KK * dVolStrain * mI1
                dAlpha1 = 0.0
                dFabric1 = 0.0
                dPStrain1 = 0.0
    
            else
                tmp0 = nn * (2.0 * GG)
                tmp1 = mI1 * (KK * DD)
                tmp1 = tmp1 + tmp0
                call Macauley(NextL, Macauley_result)
                tmp1 = tmp1 * (-Macauley_result)
                tmp2 = mI1 * (KK * dVolStrain)
                call ToContraviant(dDevStrain, ToContraviant_result)
                dSigma1 = ToContraviant_result * (2.0 * GG)
                dSigma1 = dSigma1 + tmp2 + tmp1

                ! update fabric				
                !// if (DoubleDot2_2_Contr(alphaD - NextAlpha, n) < 0.0) {
                alphaD_NextAlpha = alphaD - NextAlpha
                call DoubleDot2_2_Contr(alphaD_NextAlpha, nn, DoubleDot2_2_Contr_result)
                if (DoubleDot2_2_Contr_result < 0.0) then
                ! Equation 57
                dFabric1 = nn * m_z_max
                dFabric1 = dFabric1 + NextFabric
                                
                !dFabric1 = dFabric1 * (-m_cz / (1 + Macauley(mzcum / (2.0 * m_z_max) - 1.0)) * Macauley(NextL) * MacauleyIndex(-D))
                call Macauley(mzcum / (2.0 * m_z_max) - 1.0, Macauley_result)
                call Macauley(NextL, Macauley_result2)
                call MacauleyIndex(-DD, MacauleyIndex_result)
                dFabric1 = dFabric1 * (-m_cz / (1 + Macauley_result) * Macauley_result2 * MacauleyIndex_result)
                end if
                call ToCovariant(R1, ToCovariant_result)
                dPStrain1 = ToCovariant_result * NextL
                dAlpha1 = bb * (two3 * NextL * hh)
            end if
            
        end if

        ! Calculate Delta 2
        ! Implement the calculation of dSigma2, dAlpha2, dFabric2, and dPStrain2
        ! tmp0 = NextStress; tmp0 += dSigma1;   ! tmp0 is NextStress + dSigma1 until line 1511
        tmp0 = NextStress + dSigma1
        GetTrace_result = 0.0
        call GetTrace(tmp0, GetTrace_result)
        pp = 0.5 * GetTrace_result
        if (pp < 0) then
            if (dT == dT_min) then
        
                !if (debugFlag) then
                !
                !    write(*,*) "Delta 1: p < 0"
                !
                !endif
                NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
                NextStress = CurStress
                NextAlpha = CurAlpha
                NextFabric = CurFabric
        
                return
    
            endif
    
            dT = max(0.1 * dT, dT_min)
    
            cycle

        endif

        ! tmp1.Zero();  tmp1 += NextAlpha; tmp1 += dAlpha1;  ! tmp1 is NextAlpha + dAlpha1
        tmp1 = NextAlpha + dAlpha1
        ! tmp2.Zero();  tmp2 += NextFabric; tmp2 += dFabric1;  ! tmp2 is NextFabric + dFabric1
        tmp2 = NextFabric + dFabric1
        call GetStateDependent(tmp0, tmp1, alpha_in, alpha_in_p, tmp2, mFabric_in, GG, mzcum &
            , mzpeak, mpzp, mMcur, NextDr, nn, DD, R2, mKp, alphaD, Cka, hh, bb, AlphaAlphaBDotN)
        ! r = GetDevPart(NextStress + dSigma1) / p;
        call GetDevPart(tmp0, GetDevPart_result)
        rr = GetDevPart_result / pp
        call DoubleDot2_2_Contr(nn, rr, DoubleDot2_2_Contr_result)
        temp4 = mKp + 2 * GG - KK * DD * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(nn, rr)
        if (abs(temp4) < small) then
        
            ! neutral loading
            dSigma2 = 0
            dAlpha2 = 0
            dFabric2 = 0
            ! dPStrain2 = dDevStrain + dVolStrain * mI1;
            dPStrain2 = dPStrain1

        else
            call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
            call DoubleDot2_2_Contr(nn, rr, DoubleDot2_2_Contr_result)

            !NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4
            NextL = (2 * GG * DoubleDot2_2_Mixed_result - DoubleDot2_2_Contr_result * KK * dVolStrain) / temp4
            mDGamma = NextL
            if (NextL < 0) then
                !if (debugFlag) then
                !    write(*,*) "2 NextL is smaller than 0"
                !    write(*,*) "NextL = ", NextL
                !endif
        
                call ToContraviant(dDevStrain, ToContraviant_result)
                
                dSigma2 = 2 * GG * ToContraviant_result + KK * dVolStrain * mI1
                dAlpha2 = 0
                dFabric2 = 0
                dPStrain2 = 0
                ! dSigma2.Zero();
                ! dPStrain2 = dPStrain1;
    
            
            else
        
                ! dSigma2 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextL) *
                !  (2.0 * G * n + K * D * mI1);
                tmp0 = nn * (2.0 * GG)
                tmp1 = mI1 * (KK * DD)
                tmp1 = tmp1 + tmp0
                
                call Macauley(NextL, Macauley_result)
                
                tmp1 = tmp1 * (-Macauley_result)
                tmp2 = mI1 * (KK * dVolStrain)
                
                call ToContraviant(dDevStrain, ToContraviant_result)
                
                dSigma2 = ToContraviant_result * (2.0 * GG)
                dSigma2 = dSigma2 + tmp2 + tmp1
                ! update fabric
                ! if (DoubleDot2_2_Contr(alphaD - (NextAlpha + dAlpha1), n) < 0.0) then
                alphaD_NextAlpha = alphaD - NextAlpha - dAlpha1
                call DoubleDot2_2_Contr(alphaD_NextAlpha, nn, DoubleDot2_2_Contr_result)
                if (DoubleDot2_2_Contr_result < 0.0) then
                    ! Equation 57
                    ! dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + NextFabric + dFabric1);
                    dFabric2 = nn * m_z_max
                    dFabric2 = dFabric2 + NextFabric + dFabric1
                    Macauley_result = 0.0
                    Macauley_result2 = 0.0
                    MacauleyIndex_result = 0.0
                    call Macauley(mzcum / 2.0 / m_z_max - 1.0, Macauley_result)
                    call Macauley(NextL, Macauley_result2)
                    call MacauleyIndex(-DD, MacauleyIndex_result)
                    !dFabric2 = dFabric2 * (-1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL) * MacauleyIndex(-D))
                    dFabric2 = dFabric2 * (-1.0 * m_cz / (1 + Macauley_result) * Macauley_result2 * MacauleyIndex_result)

                end if
                ! dPStrain2 = NextL * mIIco * R2;
                ! dAlpha2 = two3 * NextL * h * b;
                ToCovariant_result = 0.0
                call ToCovariant(R2,ToCovariant_result)
                dPStrain2 = ToCovariant_result * NextL
                
                dAlpha2 = bb * (two3 * NextL * hh)
            
            
            endif

        endif

        
        ! Update nStress, nFabric, and nAlpha
        nStress = 0.5 * (dSigma1 + dSigma2) + NextStress
        nFabric = 0.5 * (dFabric1 + dFabric2) + NextFabric
        nAlpha = 0.5 * (dAlpha1 + dAlpha2) + NextAlpha
        
        

        ! Check p
        call GetTrace(nStress, GetTrace_result)
        pp = 0.5 * GetTrace_result !GetTrace(nStress)

        ! Perform stress correction if needed
        if (pp < 0) then
            if (dT == dT_min) then
                ! Implement the stress correction
            end if
                
            NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
            NextStress = CurStress
            NextAlpha = CurAlpha
            NextFabric = CurFabric
            return
        end if

        ! Calculate current step error
        call GetNorm_Contr(NextStress, GetNorm_Contr_result)
        stressNorm = GetNorm_Contr_result!GetNorm_Contr(NextStress)
        tmp0 = dSigma2 - dSigma1
        if (stressNorm < 0.5) then
            GetNorm_Contr_result = 0.0
            call GetNorm_Contr(tmp0, GetNorm_Contr_result)
            curStepError = GetNorm_Contr_result!GetNorm_Contr(tmp0)
        else
            GetNorm_Contr_result = 0.0
            call GetNorm_Contr(tmp0, GetNorm_Contr_result)
            !curStepError = GetNorm_Contr(tmp0) / (2 * stressNorm)
            curStepError = GetNorm_Contr_result / (2 * stressNorm)

        end if

        ! Check if the current step error is within tolerance
        if (curStepError > TolE) then
            ! Adjust step size
            qq = max(0.8 * sqrt(TolE / curStepError), 0.1)
            if (dT == dT_min) then
                ! Implement the step size adjustment and stress correction
            end if
            dT = max(qq * dT, dT_min)
        else
            ! Implement the step size adjustment and update NextElasticStrain, NextStress, NextAlpha, and NextFabric
            NextElasticStrain = NextElasticStrain - 0.5 * (dPStrain1 + dPStrain2)
            NextStress = nStress
            NextAlpha = nAlpha
            NextFabric = nFabric
            TT = TT + dT
            qq = max(0.8 * sqrt(TolE / curStepError), 0.5)
            dT = max(qq * dT, dT_min)
            dT = min(dT, 1.0 - TT)
        
        end if
    
    end do

    end subroutine ModifiedEuler
    
    
    
    
    
    
    
    !--------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Runge-Kutta Integrator
    !/*************************************************************/
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !--------------------------------------------------------------
    !/*************************************************************/
    !//            Pegasus Iterations                             //
    !/*************************************************************/
    subroutine IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1, aa)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress, CurStrain, NextStrain, CurAlpha
    real(REAL_TYPE), intent(out) :: aa
    ! Local variables
    real(REAL_TYPE), intent(inout) :: a0, a1
    real(REAL_TYPE) :: ff, f0, f1
    real(REAL_TYPE), dimension(3) :: dSigma, dSigma0, dSigma1, strainInc, tmp
    ! Constants
    integer(INTEGER_TYPE) :: maxIter! = 10
    real(REAL_TYPE) :: small! = 1.0e-10
    real(REAL_TYPE) :: mTolF! = 1.0e-8
    
    integer(INTEGER_TYPE) :: ii
    real(REAL_TYPE) :: GetFYieldFunction_result
    
    real(REAL_TYPE), dimension(3,3) :: mCe
    
    
    maxIter = 10
    small = 1.0e-10
    mTolF = 1.0e-8

    ! Initialize strain increment
    strainInc = NextStrain - CurStrain

    ! Calculate dSigma0 and f0
    dSigma0 = 0.0
    !dSigma0 = dSigma0 + mCe * (a0 * strainInc)
    dSigma0 = dSigma0 + matmul(mCe, (a0 * strainInc))
    tmp = CurStress + dSigma0
    GetFYieldFunction_result = 0 
    call GetFYieldFunction(tmp, CurAlpha, GetFYieldFunction_result)
    f0 = GetFYieldFunction_result

    ! Calculate dSigma1 and f1
    dSigma1 = 0.0
    !dSigma1 = dSigma1 + mCe * (a1 * strainInc)
    dSigma1 = dSigma1 + matmul(mCe, (a1 * strainInc))

    tmp = CurStress + dSigma1
    GetFYieldFunction_result = 0
    call GetFYieldFunction(tmp, CurAlpha, GetFYieldFunction_result)
    f1 = GetFYieldFunction_result!GetFYieldFunction(tmp, CurAlpha)

    ! Iterate using the Pegasus method
    do ii = 1, maxIter
        ! Calculate new a
        aa = a1 - f1 * (a1 - a0) / (f1 - f0)
        !dSigma = mCe * (aa * strainInc)
        dSigma = matmul(mCe , (aa * strainInc))
        tmp = CurStress + dSigma
        call GetFYieldFunction(tmp, CurAlpha, GetFYieldFunction_result)
        ff = GetFYieldFunction_result

        ! Check convergence
        if (abs(ff) < mTolF) then
            exit
        endif

        ! Update bounds for the next iteration
        if (ff * f0 < 0) then
            a1 = aa
            f1 = ff
        else
            f1 = f1 * f0 / (f0 + ff)
            a0 = aa
            f0 = ff
        endif

        ! Check for non-convergence
        if (ii == maxIter) then
            !if (debugFlag) then
            !    write(*, *) "Didn't find alpha!"
            !endif
            aa = 0.0
            exit
        endif
    end do

    ! Ensure a is within bounds
    if (aa > 1.0 - small) aa = 1.0
    if (aa < small) aa = 0.0
    if (isnan(aa)) then
        !if (debugFlag) then
        !    write(*, *) "aa is NaN"
        !endif
        aa = 0.0
    endif

    end subroutine IntersectionFactor
    
    !----------------------------------------------------------------
    !/*************************************************************/
    !//      Pegasus Iterations  (ElastoPlastic Unloading)        //
    !/*************************************************************/
    
    subroutine IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha, factor)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress, CurStrain, NextStrain, CurAlpha
    real(REAL_TYPE), intent(out) :: factor
    ! Local variables
    real(REAL_TYPE) :: aa, a0, a1, da
    real(REAL_TYPE) :: ff, f0, f1, fs
    integer(INTEGER_TYPE) :: nSub! = 20
    real(REAL_TYPE), dimension(3) :: dSigma, dSigma0, dSigma1, strainInc, tmp
    logical :: flag = .false.
    ! Constants
    integer(INTEGER_TYPE) :: maxIter! = 10
    real(REAL_TYPE) :: mTolF! = 1.0e-8
    
    ! local variables 
    real(REAL_TYPE) :: GetFYieldFunction_result
    real(REAL_TYPE), dimension(3,3) :: mCe
    
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    
    integer(INTEGER_TYPE) :: ii, kk
    
    aa = 0.0
    a0 = 0.0 
    a1 = 1.0
    
    nSub = 20
    maxIter = 10
    mTolF = 1.0e-8
    
    ! Initialize strain increment
    strainInc = NextStrain - CurStrain
    ! Compute initial value of f0
    call GetFYieldFunction(CurStress, CurAlpha, GetFYieldFunction_result)
    f0 = GetFYieldFunction_result
    fs = f0
    
    ! Compute initial value of dSigma
    call DoubleDot4_2(mCe, strainInc, DoubleDot4_2_result)
    dSigma = DoubleDot4_2_result!DoubleDot4_2(mCe, strainInc)

    ! Pegasus iterations loop
    do ii = 1, maxIter
        ! Compute increment factor
        da = (a1 - a0) / real(nSub)
        ! Iterate over subintervals
        do kk = 1, nSub
            aa = a0 + da * real(kk)
            ! Compute stress increment
            tmp = aa * dSigma + CurStress
            ! Compute value of f
            call GetFYieldFunction(tmp, CurAlpha, GetFYieldFunction_result)
            ff = GetFYieldFunction_result
            ! Check if f exceeds tolerance
            if (ff > mTolF) then
                a1 = aa
                if (f0 < -mTolF) then
                    f1 = ff
                    flag = .true.
                    exit
                else
                    a0 = 0.0
                    f0 = fs
                    exit
                endif
            else
                a0 = aa
                f0 = ff
            endif
        enddo

        ! Check termination conditions
        if (flag) exit
        if (ii == maxIter) then
            !if (debugFlag) then
            !    write(*, *) "Didn't find alpha! - Unloading, a0 = ", a0, ", a1 = ", a1
            !endif
            factor = 0.0
            return
        endif
    enddo

    ! Output the result
    !if (debugFlag) then
    !    write(*, *) "Found alpha - Unloading, a0 = ", a0, ", a1 = ", a1
    !endif
    call IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1, aa) 
    factor = aa
    !--> temporary


    end subroutine IntersectionFactor_Unloading
    
    
    
    
    
    
    
    
    !--------------------------------------------------------
    !/*************************************************************/
    !//            Stress Correction                              //
    !/*************************************************************/
    subroutine Stress_Correction_(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(inout) :: NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric
    real(REAL_TYPE), intent(inout) :: NextVoidRatio
    
    ! Local variables
    real(REAL_TYPE), dimension(3) :: dSigmaP, dfrOverdSigma, dfrOverdAlpha, nn, RR, alphaD, bb, aBar, rrr
    real(REAL_TYPE), dimension(3) :: nAlpha, nStress, dSigma, tmp0, tmp1
    real(REAL_TYPE) :: lambda, DD, K_p, Cka, hh, pp, fr, AlphaAlphaBDotN
    real(REAL_TYPE) :: maxIter, ii, fr_old, alpha_up, alpha_mid, alpha_down
    integer(INTEGER_TYPE) :: jj
    
    ! Local variables
    real(REAL_TYPE) :: GetTrace_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result2
    
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
	
	real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE) :: GetFYieldFunction_result
    
    real(REAL_TYPE) :: m_Pmin
    
    real(REAL_TYPE) :: mTolF
    
    real(REAL_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE) :: m_Mc, CurDr, m_emax, m_emin
    
    real(REAL_TYPE), dimension(3) :: mFabric_in
    
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mzcum 
    real(REAL_TYPE) :: mzpeak
    real(REAL_TYPE) :: mpzp
    real(REAL_TYPE) :: mMcur
    real(REAL_TYPE) :: mK
    
    real(REAL_TYPE), dimension(3,3) :: GetStiffness_result
    real(REAL_TYPE), dimension(3) :: mDGamma
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    real(REAL_TYPE) :: two3
    
    !real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE), dimension(3) :: mSigma
    

    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    ! Constants
    maxIter = 25
    mTolF = 1.0e-5
    two3 = 0.6666666666667

    
    ! Compute p
    GetTrace_result = 0
    call GetTrace(NextStress, GetTrace_result)
    pp = 0.5 * GetTrace_result!GetTrace(NextStress, GetTrace_result)

    ! Check tension cutoff
    if (pp < m_Pmin / 5.0) then
        GetFYieldFunction_result = 0.0
        call GetFYieldFunction(NextStress, NextAlpha, GetFYieldFunction_result)
        fr = GetFYieldFunction_result
        if (fr < mTolF) then
            ! Stress state inside yield surface
            NextStress = NextStress + (m_Pmin / 5.0 - pp) * mI1
        else
            ! Stress state outside yield surface
            NextStress = m_Pmin / 5.0 * mI1
            NextStress(3) = 0.8 * m_Mc * m_Pmin / 5.0
            NextAlpha = 0.0
            NextAlpha(3) = 0.8 * m_Mc
            return
        endif
    else
        ! Correct the stress
        GetFYieldFunction_result = 0
        call GetFYieldFunction(NextStress, NextAlpha, GetFYieldFunction_result)
        fr = GetFYieldFunction_result
        if (fr < mTolF) return
        ! Compute CurDr
        CurDr = (m_emax - NextVoidRatio) / (m_emax - m_emin)
        nStress = NextStress
        nAlpha = NextAlpha

        ! Stress correction loop
        do ii = 1, maxIter
            ! Compute r
            GetDevPart_result = 0.0
            call GetDevPart(nStress, GetDevPart_result)
            rrr = GetDevPart_result / pp
            ! Call GetStateDependent --> need to be checked
            call GetStateDependent(nStress, nAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum, &
                                    mzpeak, mpzp, mMcur, CurDr, nn, DD, RR, K_p, alphaD, Cka, hh, bb, AlphaAlphaBDotN)
            ! Compute dSigmaP
            call ToCovariant(RR, ToCovariant_result)
            call GetStiffness(mK, mG, GetStiffness_result)
            !dSigmaP = DoubleDot4_2(GetStiffness(mK, mG), mDGamma * ToCovariant(R))
            ToCovariant_result = 0.0
            call DoubleDot4_2(GetStiffness_result, mDGamma * ToCovariant_result, DoubleDot4_2_result)
            dSigmaP = DoubleDot4_2_result !DoubleDot4_2(GetStiffness_result, mDGamma * ToCovariant_result)
            ! Compute aBar
            aBar = two3 * hh * bb
            DoubleDot2_2_Contr_result = 0.0
            call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
            ! Compute dfrOverdSigma
            dfrOverdSigma = -0.5 * DoubleDot2_2_Contr_result * mI1 + nn
            ! Compute dfrOverdAlpha
            dfrOverdAlpha = -pp * nn
            
            ! Compute lambda
            DoubleDot2_2_Contr_result = 0.0
            DoubleDot2_2_Contr_result2 = 0.0
            call DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP, DoubleDot2_2_Contr_result)
            call DoubleDot2_2_Contr(dfrOverdAlpha, aBar, DoubleDot2_2_Contr_result2)
            lambda = fr / (DoubleDot2_2_Contr_result - DoubleDot2_2_Contr_result2)
            ! Update NextStress and NextAlpha
            tmp0 = -lambda * dSigmaP + nStress
            tmp1 = lambda * aBar + nAlpha
            GetFYieldFunction_result = 0.0
            call GetFYieldFunction(tmp0, tmp1, GetFYieldFunction_result)
            if (abs(GetFYieldFunction_result) < abs(fr)) then
                nStress = tmp0
                nAlpha = tmp1
            else
                DoubleDot2_2_Contr_result = 0.0
                call DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma, DoubleDot2_2_Contr_result)
                lambda = fr / DoubleDot2_2_Contr_result !DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma)
                nStress = nStress - lambda * dfrOverdSigma
            end if
            ! Update fr
             call GetFYieldFunction(nStress, nAlpha, GetFYieldFunction_result)
             fr = GetFYieldFunction_result
            if (abs(fr) < mTolF) then
                NextStress = nStress
                NextAlpha = nAlpha
                return
            endif
            ! Update p
            GetTrace_result = 0.0
            call GetTrace(nStress, GetTrace_result)
            pp = max(0.5 * GetTrace_result, m_Pmin)
        end do
        ! Search for the optimum alpha_mid
        dSigma = NextStress - mSigma
        alpha_up = 1.0
        alpha_mid = 0.5
        alpha_down = 0.0
        tmp0 = alpha_mid * dSigma + mSigma
        call GetFYieldFunction(tmp0, NextAlpha, GetFYieldFunction_result)
        fr_old = GetFYieldFunction_result!GetF(tmp0, NextAlpha)
        do jj = 1, maxIter
            if (fr_old < 0.0) then
                alpha_down = alpha_mid
                alpha_mid = 0.5 * (alpha_up + alpha_mid)
            else
                alpha_up = alpha_mid
                alpha_mid = 0.5 * (alpha_down + alpha_mid)
            endif
            tmp0 = alpha_mid * dSigma + mSigma
            GetFYieldFunction_result = 0.0
            call GetFYieldFunction(tmp0, NextAlpha, GetFYieldFunction_result)
            fr_old = GetFYieldFunction_result!GetF(tmp0, NextAlpha)
            if (abs(fr_old) < mTolF) then
                NextStress = tmp0
                exit
            endif
        end do
    endif

end subroutine Stress_Correction_
    
    
    
    
    
    
    
    
    
    !---------------------------------------------------------------
    
    subroutine Stress_Correction(NextStress, NextAlpha, dAlpha, mm, RR, nn, rrr)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(inout) :: NextStress, NextAlpha, dAlpha, RR, nn, rrr
    real(REAL_TYPE) :: mm
    ! Local variables
    real(REAL_TYPE), dimension(3) :: dfrOverdSigma
    real(REAL_TYPE) :: lambda
    integer(INTEGER_TYPE) :: maxIter, ii
    real(REAL_TYPE) :: ff
    
    ! Local variables 
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: GetFYieldFunction_result
    real(REAL_TYPE) :: mTolF
    real(REAL_TYPE), dimension(3) :: mI1
    real(REAL_TYPE) :: root12

    ! Constants
    maxIter = 50
    mTolF = 1.0e-5
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    root12 = 0.70710678118
    
    call GetFYieldFunction(NextStress, NextAlpha, ff)
    !ff = 

    ! Check if correction is needed
    if (ff < mTolF) return

    ! Stress correction loop
    do ii = 1, maxIter
        ! Compute dfrOverdSigma
        call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
        dfrOverdSigma = nn - 0.5 * DoubleDot2_2_Contr_result * mI1

        ! Compute lambda
        call DoubleDot2_2_Contr(dfrOverdSigma, RR, DoubleDot2_2_Contr_result)
        lambda = ff / DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(dfrOverdSigma, R)

        ! Update NextStress and NextAlpha
        NextStress = NextStress - (RR * lambda)
        NextAlpha = NextAlpha - (dAlpha * lambda)

        ! Recompute f
        call GetFYieldFunction(NextStress, NextAlpha, GetFYieldFunction_result)
        ff = GetFYieldFunction_result!GetF(NextStress, NextAlpha)

        ! Check convergence
        if (abs(ff) < mTolF) exit

        ! Check maximum iterations
        !if (ii == maxIter) then
        !    if (debugFlag) print *, "Still outside with f = ", ff
        !endif
    end do

    end subroutine Stress_Correction
    
     
    !---------------------------------------------------------------
    !/*************************************************************/
    !/*************************************************************/
    !//            MATERIAL SPECIFIC METHODS                      //
    !/*************************************************************/
    !/*************************************************************/
    !// Macauley() -------------------------------------------------
    
    subroutine Macauley(xx, xx_Macauley) 
    !--------------------------------------------
    ! To do Macauley bracket. If it is negative, set to zero.
    !--------------------------------------------
    implicit none 
    real(REAL_TYPE), intent(in) :: xx
    real(REAL_TYPE), intent(out) :: xx_Macauley
    
    !// Macauley bracket
    if (xx>0) then 
        xx_Macauley = xx
    else 
        xx_Macauley = 0
    end if 
    
    end subroutine Macauley
    
    
    
    !---------------------------------------------------
    
    subroutine MacauleyIndex(xx, MacauleyIndex_result) 
    !--------------------------------------------
    ! To do Macauley bracket. If it is negative, set to zero.
    !--------------------------------------------
    implicit none 
    real(REAL_TYPE), intent(in) :: xx
    real(REAL_TYPE), intent(out) :: MacauleyIndex_result
    
    !// Macauley bracket
    if (xx>0) then 
        MacauleyIndex_result = 1
    else 
        MacauleyIndex_result = 0
    end if 
    
    end subroutine MacauleyIndex
    
    
    
    
   
    
    !-----------------------------------------------------
    
    subroutine GetFYieldFunction(nStress, nAlpha, ff) 
    !-----------------------------------------------------
    ! To evaluate the yield function. 
    !-----------------------------------------------------
    implicit none 
    
    ! input parameters 
    real(REAL_TYPE), intent(in), dimension(3) :: nAlpha
    real(REAL_TYPE), intent(in), dimension(3) :: nStress
    
    ! local parameters
    real(REAL_TYPE), dimension(3) :: ss ! deviator stress
    real(REAL_TYPE) :: pp
    
    ! output parameters
    real(REAL_TYPE) :: ff ! scalar
    
    ! Local parameters 
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    
    real(REAL_TYPE) :: root12
    real(REAL_TYPE) :: m_m
    
    real(REAL_TYPE) :: m_pmin
    
    
    root12 = 0.70710678118
    
    !// PM4Sand's yield function
	!// s = nStress - p*I
	!//   = [s_xx	s_xy ;	=	[Sigma_xx-p		Sigma_xy   ;
	!//		s_xy	s_yy];		 Sigma_xy		Sigma_yy-p];
	!// f = [(s - p*alpha):(s - p*alpha)]^0.5 - sqrt(1/2)*p*m
    
    call GetTrace(nStress, GetTrace_result)
    
    call GetDevPart(nStress, GetDevPart_result)
    ss = GetDevPart_result !GetDevPart(nStress)
    pp = 0.5 * GetTrace_result
    
    ss = ss - (pp*nAlpha)
    
    call GetNorm_Contr(ss, GetNorm_Contr_result)
    
    ff = GetNorm_Contr_result - (root12 * m_m * pp)
    
    end subroutine GetFYieldFunction
    
    
    !----------------------------------------------------
    
    subroutine GetKsi(Dr, pp, ksi)!, &
                     !m_R, m_Q)
    !-------------------------------------------
    ! To calculate the relative state parameter 
    ! based on Bolton (1986).
    !-------------------------------------------
    implicit none
    
    ! input parameters
    real(REAL_TYPE), intent(in) :: Dr
    real(REAL_TYPE), intent(in) :: pp
    !real(REAL_TYPE), intent(in) :: m_R
    !real(REAL_TYPE), intent(in) :: m_Q
    ! I am not sure how can we make m_R and m_Q being global
    
    ! output parameters 
    real(REAL_TYPE), intent(out) :: ksi
    
    ! local parameters
    real(REAL_TYPE) :: pn
    real(REAL_TYPE) :: m_pmin
    real(REAL_TYPE) :: m_R
    real(REAL_TYPE) :: m_Q
    real(REAL_TYPE) :: m_p_atm
    
    pn = pp
    
    if (pn<=m_Pmin) then 
        pn = m_Pmin
    end if 
    
    !// find ksi_R
    !// calculating the relative density based on Bolton (1986)
    ksi = (m_R / (m_Q - log(100.0 * pn / m_P_atm)) ) - Dr
    
    
                     end subroutine GetKsi
    !-------------------------------------------------
    subroutine GetElasticModuli_(sigma, KK, GG, Mcur, zcum)
    !--------------------------------------------
    ! To calculate the elastic modulus of the soil
    !--------------------------------------------
    
    implicit none 
    
    real(REAL_TYPE) :: two3, Csr0, temp, p, p_sed, Macauley_result, F_sed
    real(REAL_TYPE), intent(in) :: sigma(3) !out
    real(REAL_TYPE), intent(out) :: KK, GG, Mcur, zcum
    real(REAL_TYPE) :: GetTrace_sigma, pn, qn, Csr
    
    real(REAL_TYPE) :: msr, mMb, m_z_max, m_G0, m_Cgd, pp 
    
    integer(INTEGER_TYPE) :: me2p, m_PostShake
    
    real(REAL_TYPE) :: mMd, m_p_sedo, mMcur, m_Fsed_min, m_nu
    real(REAL_TYPE) :: m_pmin
    real(REAL_TYPE) :: mzcum
    !integer(INTEGER_TYPE) :: me2p
    real(REAL_TYPE) :: m_P_atm
    !real(REAL_TYPE) :: m_Pmin
    !
    !real(REAL_TYPE) :: me2p
    
    !real(REAL_TYPE) :: two3
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_Result
    
    
    two3 = 0.666666666667
    msr = 4.0
    Csr0 = 0.5
    
    call GetTrace(sigma, GetTrace_sigma)
    pn = 0.5 * GetTrace_Sigma
    
    if (pn <= m_Pmin) then 
        pn = m_Pmin 
    end if 
    
    qn = 2 * sqrt( (0.5*(sigma(1) - sigma(2)))**2.0 + sigma(3)**2.0)
    
    Mcur = qn/pn
    Csr = 1 - ( Csr0 * min(1.0, (Mcur/mMb)**msr) ) 
    
    temp = zcum / m_z_max
    
    if (me2p == 0) then 
        GG = m_G0 * m_P_atm
    else 
        GG = m_G0 * m_P_atm * sqrt(pn / m_P_atm) * Csr * (1 + temp) / (1 + temp * m_Cgd) 
        !// this equation is in the table 
		!// "Elastic deviatoric strain increment"
        
        if (m_PostShake) then 
			!// reduce elastic shear modulus for post shaking consolidation
            call GetTrace(sigma, GetTrace_sigma)
			pp = 0.5 * GetTrace_sigma
            call Macauley(1 - mMcur / mMd, Macauley_result)
			p_sed = m_p_sedo * (mzcum / (mzcum + m_z_max)) * Macauley_result**0.25 !// Equation 88
			F_sed = min(m_Fsed_min + (1 - m_Fsed_min) * (pp / 20.0 / (p_sed + small)), 1.0) !// Equation 87
			GG = GG * F_sed !// Equation 85
		end if 
        
    end if 
    
        
    if (m_nu==0.5)then 
        m_nu = 0.4999
    end if 
    
    KK = two3 * ( (1+m_nu)/(1 - (2*m_nu)) ) * GG
    !// Equation 86 (kind of...) 
	!// In the table "Elastic volumetric strain increment"
    
    !//
	!// K = (2/3) * ( 2*(1+nu)/3(1-2*nu) ) * G "Elastic volumetric strain increment"
	!//
    
    end subroutine GetElasticModuli_
    
    
    !---------------------------------------------------
    
    
    subroutine GetElasticModuli(sigma, KK, GG)
    !---------------------------------------------------
    ! To get the elastic parameters
    !---------------------------------------------------
    implicit none 
    
    real(REAL_TYPE), intent(in) :: sigma(3)
    real(REAL_TYPE), intent(out) :: KK, GG
    real(REAL_TYPE) :: GetTrace_sigma, pn
    
    ! Local variables 
    real(REAL_TYPE) :: m_pmin
    real(REAL_TYPE) :: me2p
    real(REAL_TYPE) :: m_G0
    real(REAL_TYPE) :: m_p_atm
    real(REAL_TYPE) :: m_nu
    real(REAL_TYPE) :: two3
    
    two3 = 2.0/3.0
    
    call GetTrace(sigma, GetTrace_sigma)
    
    if (pn<=m_Pmin) then 
        pn = m_Pmin
    end if 
    
    if (me2p == 0) then 
        GG = m_G0 * m_P_atm 
    else 
        GG = m_G0 * m_P_atm * sqrt(pn / m_P_atm)
    end if 
    
    if (m_nu==0.5) then 
        m_nu =0.499999
    end if 
    
    
    KK = two3 * ( (1+m_nu)/(1-(2*m_nu)) ) * GG !--> two3 needs to be a global number
    
    
    end subroutine GetElasticModuli                 
                     
                     
                     
                     
    !-------------------------------------------------                
    subroutine GetStiffness(KK, GG, CC)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), intent(in) :: KK, GG
    real(REAL_TYPE), intent(out) :: CC(3, 3)
    ! Local variables
    real(REAL_TYPE) :: aa, bb

    ! Compute coefficients
    aa = KK + 4.0 / 3.0 * GG
    bb = KK - 2.0 / 3.0 * GG

    ! Fill stiffness matrix
    CC(1, 1) = aa
    CC(2, 2) = aa
    CC(3, 3) = GG
    CC(1, 2) = bb
    CC(2, 1) = bb


    end subroutine GetStiffness              
    
    
    
    
    !------------------------------------------------                 
    !subroutine GetCompliance()                 
    !               
    !implicit none
    !
    !
    !
    !
    !
    !end subroutine 
    
    
    !// returns the compliance matrix in its covariant-covariant form
    subroutine GetCompliance(KK, GG, DD)
    implicit none
    ! Inputs
    real(REAL_TYPE), intent(in) :: KK, GG
    ! Output
    real(REAL_TYPE), dimension(3, 3), intent(out) :: DD
    ! Local variables
    real(REAL_TYPE) :: aa, bb, cc

    ! Compute components of compliance matrix
    aa = (KK + 4.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
    bb = (KK - 2.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
    cc = 1.0 / GG

    ! Assign values to the compliance matrix
    DD(1, 1) = aa
    DD(2, 2) = aa
    DD(3, 3) = cc
    DD(1, 2) = bb
    DD(2, 1) = bb

end subroutine GetCompliance
    
    !------------------------------------------------
                     
    !subroutine GetElastoPlasticTangent()                 
    !          
    !implicit none
    !                 
    !              
    !
    !
    !end subroutine GetElastoPlasticTangent
    
    
    
    
    
    subroutine GetElastoPlasticTangent(NextStress, aCe, RR, nn, K_p, aCep)
    implicit none
    ! Inputs
    real(REAL_TYPE), dimension(3), intent(in) :: NextStress
    real(REAL_TYPE), dimension(3, 3), intent(in) :: aCe
    real(REAL_TYPE), dimension(3), intent(in) :: RR
    real(REAL_TYPE), dimension(3), intent(in) :: nn
    real(REAL_TYPE), intent(in) :: K_p
    ! Output
    real(REAL_TYPE), dimension(3, 3), intent(out) :: aCep
    ! Local variables
    real(REAL_TYPE) :: pp, temp3
    real(REAL_TYPE), dimension(3) :: rrr, temp1, temp2

    ! Local variables 
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: m_pmin
    
    real(REAL_TYPE), dimension(3) :: DoubleDot2_4_result
    
    integer(INTEGER_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE), dimension(3,3) :: Dyadic2_2_result
    
    real(REAL_TYPE), dimension(3,3) :: mIIco
    
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    mIIco(1,1) = 1
    mIIco(2,1) = 1
    mIIco(3,1) = 1
    
    mIIco(1,2) = 1
    mIIco(2,2) = 1
    mIIco(3,2) = 1
    
    mIIco(1,3) = 1
    mIIco(2,3) = 1
    mIIco(3,4) = 1
    
    
    
    
    
    
    
    
    ! Calculate mean effective stress
    call GetTrace(NextStress, GetTrace_result)
    pp = 0.5 * GetTrace_result
    if (pp < m_Pmin) pp = m_Pmin

    ! Calculate deviatoric part of stress and normalize
    call GetDevPart(NextStress, GetDevPart_result)
    rrr = GetDevPart_result / pp

    ! Compute temporary vectors and matrices
    call DoubleDot4_2(aCe, RR, DoubleDot4_2_result)
    temp1 = DoubleDot4_2_result!DoubleDot4_2(aCe, RR, DoubleDot4_2_result)
    
    call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
    call DoubleDot2_4((nn - 1.0 / 2.0 * DoubleDot2_2_Contr_result) * mI1, aCe * mIIco, DoubleDot2_4_result)
    temp2 = DoubleDot2_4_result !DoubleDot2_4(n - 1.0 / 2.0 * DoubleDot2_2_Contr(n, r) * mI1, aCe * mIIco)
    
    call DoubleDot2_2_Contr(temp2, RR, DoubleDot2_2_Contr_result)
    temp3 = DoubleDot2_2_Contr_result + K_p

    ! Update aCep based on plasticity condition
    if (temp3 < small) then
        aCep = aCe
    else
        call Dyadic2_2(temp1, temp2, Dyadic2_2_result)
        aCep = aCe - 1.0 / temp3 * Dyadic2_2_result!Dyadic2_2(temp1, temp2)
    end if
end subroutine GetElastoPlasticTangent
    
    
    
    
    
    
    
    
    !------------------------------------------------
    
    
    
    subroutine GetNormalToYield(stress, alpha, nn)
    !-----------------------------------------------
    ! To calculate the normal to the yield surface
    !-----------------------------------------------
    
    implicit none 
    
    ! output variables
    real(REAL_TYPE), intent(out), dimension(3) :: nn
    
    ! input variables 
    real(REAL_TYPE), intent(in), dimension(3) :: stress, alpha
    
    ! local variables 
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: small
    real(REAL_TYPE) :: root12
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_stress
    real(REAL_TYPE) :: GetTrace_stress
    
    real(REAL_TYPE) :: normN !, dimension(3)
    
    
    root12 = 0.7071067811865475
    
    call GetTrace(stress, GetTrace_stress)
    call GetDevPart(stress, GetDevPart_stress)
    pp = 0.5 * GetTrace_stress
    
    if (abs(pp)<small) then 
        nn(2) = root12
        
    else 
        nn = (alpha*-pp) + GetDevPart_stress
        call GetNorm_Contr(nn, normN)
        if (normN<small) then 
            normN = 1
        end if 
        nn = nn * (1/normN)
        
    end if 
    
    
    end subroutine GetNormalToYield
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------
    
    !subroutine GetStateDependent()    
    !!----------------------------------------------------
    !! To calculate the state parameter for PM4Sand. 
    !!----------------------------------------------------   
    !                 
    !implicit none 
    !                 
    !if (p<=m_mPmin) then 
    !    m_Pmin 
    !end if
    !
    !!// relative state parameter index
    !call GetKsi(Dr, pp, ksi, &
    !                 m_R, m_Q)
    !
    !!// get the yield to the normal
    !
    !
    !
    !end subroutine GetStateDependent
    
    
    
    !*************************************************************
    ! GetStateDependent() ----------------------------------------
    ! Alsardi: we are calculating the state parameters for PM4Sand
    subroutine GetStateDependent(stress, alpha, alpha_in, alpha_in_p, fabric, fabric_in, &
    GG, zcum, zpeak, pzp, Mcur, CurDr, nn, DD, RR, K_p, alphaD, Cka, hh, bb, AlphaAlphaBDotN)
    implicit none
    ! Input parameters
    real(REAL_TYPE), dimension(3), intent(in) :: stress
    real(REAL_TYPE), dimension(3), intent(in) :: alpha
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_p
    real(REAL_TYPE), dimension(3), intent(in) :: fabric
    real(REAL_TYPE), dimension(3), intent(in) :: fabric_in
    real(REAL_TYPE), intent(in) :: GG
    real(REAL_TYPE), intent(in) :: zcum
    real(REAL_TYPE), intent(in) :: zpeak
    real(REAL_TYPE), intent(in) :: pzp
    real(REAL_TYPE), intent(in) :: Mcur
    real(REAL_TYPE), intent(in) :: CurDr
    ! Output parameters
    real(REAL_TYPE), dimension(3), intent(out) :: nn
    real(REAL_TYPE), intent(out) :: DD
    real(REAL_TYPE), dimension(3), intent(out) :: RR
    real(REAL_TYPE), intent(out) :: K_p
    real(REAL_TYPE), dimension(3), intent(out) :: alphaD
    real(REAL_TYPE), intent(out) :: Cka
    real(REAL_TYPE), intent(out) :: hh
    real(REAL_TYPE), dimension(3), intent(out) :: bb
    real(REAL_TYPE), intent(out) :: AlphaAlphaBDotN
    ! Local variables
    real(REAL_TYPE) :: alphaD_alpha(3), alphaDr_alpha(3), alpha_mAlpha_in(3), alpha_mAlpha_in_true(3), alpha_mAlpha_p(3), minusFabric(3)
    real(REAL_TYPE) :: Czpk1, Czpk2, Cpzp2, Cg1, Ckp, AlphaAlphaInDotN, AlphaAlphaInTrueDotN, Czin1, Crot1, Mdr
    real(REAL_TYPE) :: pp, ksi
    
    
    ! Local variables --> alphaB is probably not a local variable...?
    real(REAL_TYPE) :: alphaB(3), Macauley_result, GetTrace_result, DoubleDot2_2_Contr_result, DoubleDot2_2_Contr_result2

    real(REAL_TYPE) :: m_Pmin, mMb, m_Mc, m_Nb, mMd, m_Nd, m_m, m_z_max
    real(REAL_TYPE) :: root12
    
    real(REAL_TYPE) :: m_h0
    real(REAL_TYPE) :: malpha_in
    real(REAL_TYPE) :: malpha_in_true
    real(REAL_TYPE) :: m_Ckaf
    
    real(REAL_TYPE) :: Cpzp
    real(REAL_TYPE) :: Cpmin
    
    real(REAL_TYPE) :: m_pmin2
    
    real(REAL_TYPE) :: Czin2
    
    real(REAL_TYPE) :: two3
    
    real(REAL_TYPE) :: temp
    real(REAL_TYPE) :: m_Ado
    real(REAL_TYPE) :: m_Ce
    real(REAL_TYPE) :: Ad
    
    real(REAL_TYPE), dimension(3) :: mI1
    
    real(REAL_TYPE) :: C_pmin2, Cin, Cdz
    
    real(REAL_TYPE) :: one3
    
    real(REAL_TYPE) :: hp, m_hpo, Adc
    
    real(REAL_TYPE) :: Drot, m_Cdr
    
    real(REAL_TYPE) :: Crot2
    
   
    mI1(1) = 1
    mI1(2) = 1
    mI1(3) = 1
    
    two3 = 2.0/3.0
    
    one3 = 1.0/3.0
    
    
    root12 = 0.70710678118
    
    ! Mean effective stress
    call GetTrace(stress, GetTrace_result)
    pp = 0.5d0 * GetTrace_result!GetTrace(stress)
    if (pp <= m_Pmin) pp = m_Pmin ! kind of a tension cutoff
    call GetKsi(CurDr, pp, ksi)
    !ksi = GetKsi(CurDr, p) ! relative state parameter index
    call GetNormalToYield(stress, alpha, nn) ! get the yield to the normal
    if (ksi <= 0.0d0) then ! negative relative state parameter index
        ! dense of critical
        mMb = m_Mc * exp(-1.0d0 * m_nb * ksi) ! bounding ratio based on ksi using m_nb
        mMd = m_Mc * exp(m_nd * ksi) ! dilatancy ratio update based on ksi using m_nd
    else
        ! loose of critical
        mMb = m_Mc * exp(-1.0d0 * m_nb / 4.0d0 * ksi) ! bounding ratio based on ksi using m_nb/4
        mMd = m_Mc * exp(m_nd * 4.0d0 * ksi) ! dilatancy ratio update based on ksi using m_nb*4
    end if
    ! Vector alphaB = root12 * (mMb - m_m) * n;
    alphaB = nn * root12 * (mMb - m_m) ! image back-stress ratio for bounding surface
    ! Vector alphaD = root12 * (mMd - m_m) * n;
    alphaD = nn * root12 * (mMd - m_m) ! image back-stress ratio for dilatancy surface
    ! C_zpk1 = z_peak / (z_cum + z_max/5)
    Czpk1 = zpeak / (zcum + m_z_max / 5.0d0) ! Equation 62 --> Effect of fabric on plastic modulus
    ! C_zpk2 = z_peak / (z_cum + z_max/100)
    Czpk2 = zpeak / (zcum + m_z_max / 100.0d0) ! Equation 63 --> Effect of fabric on plastic modulus
    if (Czpk2 > 1.0d0 - small) Czpk2 = 1.0d0 - small ! there is a limit on this
    ! C_pzp2 = (-<-(p_zp - p)>)/(-<-(p_zp - p)> + p_min)
    call Macauley((pzp - pp), Macauley_result)
    !Cpzp2 = Macauley((pzp - pp)) / (Macauley((pzp - pp)) + m_Pmin) ! Equation 64 --> Effect of fabric on plastic modulus
    Cpzp2 = Macauley_result / (Macauley_result + m_Pmin) ! Equation 64 --> Effect of fabric on plastic modulus

    Cg1 = m_h0 / 200.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Gamma1 = h_o / 200
    Ckp = 2.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Kp = 2
    ! b = alphaB - alpha;
    bb = alphaB - alpha
    call DoubleDot2_2_Contr(bb, nn, DoubleDot2_2_Contr_result)
    AlphaAlphaBDotN = DoubleDot2_2_Contr_result
    alpha_mAlpha_in = alpha - mAlpha_in ! alpha - alpha_in
    
    call DoubleDot2_2_Contr(alpha_mAlpha_in, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
    !AlphaAlphaInDotN = Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in, n)) ! Macauly bracket in the Kp calculation
    AlphaAlphaInDotN = Macauley_result ! Macauly bracket in the Kp calculation
    
    
    alpha_mAlpha_in_true = alpha - mAlpha_in_true ! Alpha - Alpha^true_in
    
    call DoubleDot2_2_Contr(alpha_mAlpha_in_true, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
    AlphaAlphaInTrueDotN = Macauley_result!Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in_true, n)) ! dot product in the denominator of Cka
    
    Cka = 1.0d0 + (m_Ckaf / (1.0d0 + (2.5d0 * AlphaAlphaInTrueDotN)**2)) * Cpzp2 * Czpk1 ! look at table "Plastic deviatoric strain increment"
    ! updataed K_p formulation following PM4Sand V3.1. mAlpha_in is the apparent back-stress ratio.
    alpha_mAlpha_p = alpha - alpha_in_p ! if alpha - alpha_in_p is small
    
    call DoubleDot2_2_Contr(alpha_mAlpha_p, nn, DoubleDot2_2_Contr_result)
    call Macauley(AlphaAlphaBDotN, Macauley_result)
      
    if (abs(AlphaAlphaBDotN) < small) then
        ! adding this condition to avoid division by zero error
        hh = 1.0d10 ! part of Kp calculation
    elseif (DoubleDot2_2_Contr_result <= 0) then
        ! loose of critical states 
        hh = 1.5d0 * GG * m_h0 / pp / (exp(AlphaAlphaInDotN) - 1 + Cg1) / sqrt(abs(AlphaAlphaBDotN)) * &
            Cka / (1 + Ckp * zpeak / m_z_max * Macauley_result * sqrt(1 - Czpk2)) ! part of Kp calculation 
        ! Equation 35 --> Hardening and the update of the back-stress ratio
        hh = hh * (AlphaAlphaInDotN + Cg1) / (AlphaAlphaInTrueDotN + Cg1) ! part of Kp calculation
    else 
        ! 
        hh = 1.5d0 * GG * m_h0 / pp / (exp(AlphaAlphaInDotN) - 1 + Cg1) / sqrt(abs(AlphaAlphaBDotN)) * &
            Cka / (1 + Ckp * zpeak / m_z_max * Macauley_result * sqrt(1 - Czpk2)) ! part of Kp calculation
    end if
    
        
    ! after we evaluate h we plug it in K_p to find the plastic bulk modulus
    call DoubleDot2_2_Contr(bb, nn, DoubleDot2_2_Contr_result)
    K_p = two3 * hh * pp * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(b, n) ! Equation 34 --> Hardening and the update of the back-stress ratio
    
    call DoubleDot2_2_Contr(fabric_in, nn, DoubleDot2_2_Contr_result)
    call DoubleDot2_2_Contr(fabric, nn, DoubleDot2_2_Contr_result2)

    call Macauley(1.0 - exp(-2.0*abs((DoubleDot2_2_Contr_result - DoubleDot2_2_Contr_result2) / m_z_max)), Macauley_result)
    !Czin1 = Macauley(1.0 - exp(-2.0*abs((DoubleDot2_2_Contr(fabric_in, nn) - DoubleDot2_2_Contr(fabric, n)) / m_z_max))) ! Plastic volumetric strain increment... see in table
    Czin1 = Macauley_result! Plastic volumetric strain increment... see in table
    ! rotated dilatancy surface
    minusFabric = fabric
    minusFabric = minusFabric * (-1.0)
    call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
    Crot1 = max(1.0 + 2.0 * Macauley_result / (sqrt(2.0)*m_z_max)*(1 - Czin1), 1.0) ! Plastic volumetric strain increment... see in table
    ! Equation 66 --> Effect of fabric on plastic volumetric dilation 
    Mdr = mMd / Crot1 ! Equation 65 --> Effect of fabric on plastic volumetric dilation 
    ! Vector alphaDr = root12 * (Mdr - m_m) * n;
    alphaDr_alpha = nn
    alphaDr_alpha = alphaDr_alpha * (root12 * (Mdr - m_m))
    alphaDr_alpha = alphaDr_alpha - alpha ! Equation 68
    ! Alpha_DR = sqrt(1/2) * (M^dr - m)n
    alphaD_alpha = alphaD
    alphaD_alpha = alphaD_alpha - alpha
    call DoubleDot2_2_Contr(alphaDr_alpha, nn, DoubleDot2_2_Contr_result)
    if (DoubleDot2_2_Contr_result <= 0) then
        ! dilation
        !double precision Cpzp, Cpmin, Czin2, temp, Ad, Drot
        if (pzp == 0.0) then
            Cpzp = 1.0
        else
            Cpzp = 1.0 / (1.0 + (2.5 * pp / pzp) ** 5.0) ! plastic vol strain increment
        end if
        Cpmin = 1.0 / (1.0 + (m_Pmin2 / pp) ** 2) ! plastic vol strain increment
        Czin2 = (1.0 + Czin1 * (zcum - zpeak) / 3.0 / m_z_max) / (1.0 + 3.0 * Czin1 * (zcum - zpeak) / 3.0 / m_z_max) ! plastic vol strain increment
        ! double temp = (1.0 - Macauley(DoubleDot2_2_Contr(-1.0 * fabric, n)) * root12 / zpeak) ** 3
        
        call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        
        temp = (1.0 - Macauley_result * root12 / zpeak) ** 3.0
        Ad = m_Ado * Czin2 / ((zcum ** 2 / m_z_max) * temp * (m_ce ** 2) * Cpzp * Cpmin * Czin1 + 1.0) ! plastic vol strain increment ! Equation 72
        call DoubleDot2_2_Contr(alphaD_alpha, nn, DoubleDot2_2_Contr_result)
        DD = Ad * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(alphaD_alpha, n) ! D_non-rot ! plastic vol strain increment
        ! double Drot = Ad * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max) * DoubleDot2_2_Contr(alphaDr - alpha, n) / m_Cdr
        
        call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        
        call DoubleDot2_2_Contr(alphaDr_alpha, nn, DoubleDot2_2_Contr_result2)
        !call Macauley(DoubleDot2_2_Contr_result, Macauley_result2)
        !Drot = Ad * Macauley(DoubleDot2_2_Contr(minusFabric, n)) / (sqrt(2.0) * m_z_max) * DoubleDot2_2_Contr(alphaDr_alpha, n) / m_Cdr ! D_rot ! plastic vol strain increment
        Drot = Ad * Macauley_result / (sqrt(2.0) * m_z_max) * DoubleDot2_2_Contr_result2 / m_Cdr ! D_rot ! plastic vol strain increment
        
        ! Equation 69 --> Drot
        if (DD > Drot) then ! plastic vol strain increment --> D_non-rot < D_rot
            call Macauley(mMb - Mcur, Macauley_result)
            DD = DD + (Drot - DD) * Macauley_result / (Macauley_result + 0.01) ! Equation 78
        end if
        if (m_Pmin <= pp .and. pp <= 2 * m_Pmin) then
            call Macauley(mMb - mMd, Macauley_result)
            DD = min(DD, -3.5 * m_Ado * Macauley_result * (2 * m_Pmin - pp) / m_Pmin) ! Equation 78
        end if
    else
        ! contraction
        !double precision hp, Crot2, Cdz, Adc, Cin, C_pmin2
        ! bound K_p to non - negative, following flac practice
        K_p = max(0.0, K_p)
        ! he is using power 2 instead of 2.5 in BZ14
        call Macauley(0.5 - ksi, Macauley_result)
        hp = m_hpo * exp(-0.7 + 7.0 * (Macauley_result) ** 2.0) ! Equation 53 --> Plastic volumetric strains - Contraction
        Crot2 = 1.0 - Czpk2 ! Equation 83
        Cdz = max((1.0 - Crot2 * sqrt(2.0) * zpeak / m_z_max) * (m_z_max / (m_z_max + Crot2 * zcum)), 1.0 / (1.0 + m_z_max / 2.0))
        
        call DoubleDot2_2_Contr(fabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        
        !Adc = m_Ado * (1.0 + Macauley(DoubleDot2_2_Contr(fabric, n))) / hp / Cdz
        Adc = m_Ado * (1.0 + Macauley_result) / hp / Cdz

        !Cin = 2.0 * Macauley(DoubleDot2_2_Contr(fabric, n)) / sqrt(2.0) / m_z_max
        Cin = 2.0 * Macauley_result / sqrt(2.0) / m_z_max

        ! D = min(Adc * (DoubleDot2_2_Contr(alpha - mAlpha_in, n) + Cin) ** 2, 1.5 * m_Ado) * &
        !     DoubleDot2_2_Contr(alphaD - alpha, n) / (DoubleDot2_2_Contr(alphaD - alpha, n) + 0.16)
        call DoubleDot2_2_Contr(alpha_mAlpha_in, nn, DoubleDot2_2_Contr_result)
        call DoubleDot2_2_Contr(alphaD_alpha, nn, DoubleDot2_2_Contr_result2)
        DD = min(Adc * (DoubleDot2_2_Contr_result + Cin) ** 2, 1.5 * m_Ado) * &
            DoubleDot2_2_Contr_result2 / (DoubleDot2_2_Contr_result2 + 0.16)
        ! Apply a factor to D so it doesn't go very big when p is small
        if (pp < m_Pmin * 2.0) then
            C_pmin2 = 0.0
        elseif (pp >= m_Pmin * 18.0) then
            C_pmin2 = 1.0
        else
            C_pmin2 = (pp - 2.0 * m_Pmin) / (16.0 * m_Pmin)
        end if
        
        DD = DD * C_pmin2
        
    end if 
    
    ! R = n + one3 * D * mI1
    RR = mI1
    RR = RR + (one3 * DD)
    RR = RR + nn  ! --> used in plastic deviatoric strain increment
    
    
    
    end subroutine GetStateDependent
    
    
    
    
    
    
    !/*************************************************************/
    !/*************************************************************/
    !//            SYMMETRIC TENSOR OPERATIONS                    //
    !/*************************************************************/
    !/*************************************************************/
    !// In all the functions below, by contravariant tensors, we mean stress-like tensors
    !// and by covariant tensors we mean strain-like tensors
    subroutine GetTrace(vv, GetTrace_result)
    !// computes the trace of the input argument
    implicit none
    ! Arguments
    real, dimension(3), intent(in) :: vv
    real, intent(out) :: GetTrace_result
    
    ! initialize
    GetTrace_result = 0.0

    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetTrace requires vector of size(3)!"
        GetTrace_result = 0.0
        return
    endif
    
    ! Compute the trace
    GetTrace_result = vv(1) + vv(2)

    end subroutine GetTrace
    
    
    
    subroutine GetDevPart(Vector, DevVector)
    !-------------------------------------------
    ! To calculate the deviatoric part of a matrix.
    !-------------------------------------------
    implicit none
    
    ! input variable 
    real(REAL_TYPE), dimension(3), intent(in) :: Vector
    real(REAL_TYPE), dimension(3), intent(out) :: DevVector
    
    ! local variables 
    real(REAL_TYPE) :: pp
    
    ! initialize 
    DevVector = 0.0
    pp = 0.0
    
    ! Check vector size
    if (size(Vector) /= 3) then
        print *, "ERROR! PM4Sand::GetDevPart requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    call GetTrace(Vector, pp)
    
    DevVector = Vector 
    DevVector(1) = DevVector(1) - 0.5*pp
    DevVector(2) = DevVector(2) - 0.5*pp
    
    
    
    end subroutine GetDevPart
    
    
    
    
    
    subroutine DoubleDot2_2_Contr(v1, v2, DoubleDot2_2_Contr_result)
    !--------------------------------------------------
    ! // computes doubledot product for vector-vector arguments,
    !    both "contravariant".
    implicit none 
    
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Contr_result
    integer(INTEGER_TYPE) :: ii
    
    ! initialize 
    DoubleDot2_2_Contr_result = 0.0

    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Contr requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Contr_result = 0 
    
    do ii = 1, 3 ! 3 is for the size of the vector 
        !result = result + v1(i) * v2(i) + (i > 1) * v1(i) * v2(i)
        DoubleDot2_2_Contr_result = DoubleDot2_2_Contr_result + ( v1(ii) * v2(ii) ) !+ (ii>2)*( v1(ii) * v2(ii) )
        if (ii == 3) then 
            ! double count shear component according to cauchy stress tensor
            DoubleDot2_2_Contr_result = DoubleDot2_2_Contr_result + ( v1(ii) * v2(ii) )
        end if 
    end do 
    
    end subroutine DoubleDot2_2_Contr
    
    
    
    
    
    subroutine DoubleDot2_2_Cov(v1, v2, DoubleDot2_2_Cov_result)
    !--------------------------------------------------------
    ! // computes doubledot product for vector-vector arguments, both "covariant"
    !
    implicit none 
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Cov_result
    integer(INTEGER_TYPE) :: ii
    
      
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Cov requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Cov_result = 0
    
    
    ! Compute the doubledot product
    do ii = 1, 3
        !result = result + v1(i) * v2(i) - (i > 2) * 0.5 * v1(i) * v2(i)
        DoubleDot2_2_Cov_result = DoubleDot2_2_Cov_result + (v1(ii) * v2(ii))
        if (ii == 3) then 
            DoubleDot2_2_Cov_result = DoubleDot2_2_Cov_result - (0.5 * v1(ii) * v2(ii))
        end if 
        
    end do 
    
    end subroutine DoubleDot2_2_Cov
    
    
    
    
    subroutine DoubleDot2_2_Mixed(v1, v2, DoubleDot2_2_Mixed_result)
    ! // computes doubledot product for vector-vector arguments, 
    ! one "covariant" and the other "contravariant" 
    implicit none 
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Mixed_result
    integer(INTEGER_TYPE) :: ii
    
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Mixed requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Mixed_result = 0
    
    do ii = 1, 3
        DoubleDot2_2_Mixed_result = DoubleDot2_2_Mixed_result + ( v1(ii) * v2(ii) )
    end do
    
    end subroutine DoubleDot2_2_Mixed
    
    
    
    
    subroutine GetNorm_Contr(vv, GetNorm_Contr_result)
    ! // computes contravariant (stress-like) norm of input 6x1 tensor
    implicit none
    
    real(REAL_TYPE), dimension(3), intent(in) :: vv
    real(REAL_TYPE), intent(out) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    
    
    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetNorm_Contr requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    GetNorm_Contr_result = 0
    DoubleDot2_2_Contr_result = 0 
    
    call DoubleDot2_2_Contr(vv, vv, DoubleDot2_2_Contr_result)
    
    GetNorm_Contr_result = sqrt(DoubleDot2_2_Contr_result)
    
    end subroutine GetNorm_Contr
    
    
    
    
    subroutine GetNorm_Cov(vv, GetNorm_Cov_result)
    ! // computes covariant (strain-like) norm of input 6x1 tensor
    implicit none 
    
    real(REAL_TYPE), dimension(3), intent(in) :: vv
    real(REAL_TYPE), intent(out) :: GetNorm_Cov_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Cov_result
    
    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetNorm_Cov requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    GetNorm_Cov_result = 0
    DoubleDot2_2_Cov_result = 0 
    
    ! evaluate GetNorm_Cov_result using DoubleDot2_2_Cov
    call DoubleDot2_2_Cov(vv, vv, DoubleDot2_2_Cov_result)
    
    ! square root to find GetNorm_Cov
    GetNorm_Cov_result = sqrt(DoubleDot2_2_Cov_result)
    
    
    end subroutine GetNorm_Cov
    
    
    
    
    
    
    subroutine Dyadic2_2(v1, v2, Dyadic2_2_result)
    !// computes dyadic product for two vector-storage arguments
    !// the coordinate form of the result depends on the coordinate form of inputs
    ! https://en.wikipedia.org/wiki/Dyadics
    implicit none 
    
    ! input vectors v1 and v2
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    ! output matrix result Dyadic2_2_result
    real(REAL_TYPE), dimension(3,3), intent(out) :: Dyadic2_2_result
    
    ! local variables 
    integer(INTEGER_TYPE) :: jj, ii
    
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::Dyadic2_2 requires vector of size(3)!"
        return
    endif
    
    ! initialize matrix
    Dyadic2_2_result = 0.0
    
    do ii = 1, 3
        do jj = 1, 3
            
            Dyadic2_2_result(ii,jj) = v1(ii) * v2(jj)
            
        end do 
    end do 
    
    
    
    
    end subroutine Dyadic2_2
    
    
    
    
    subroutine DoubleDot4_2(m1, v1, DoubleDot4_2_result)
    !// computes doubledot product for matrix-vector arguments
    !// caution: second coordinate of the matrix should be in opposite variant form of vector
    implicit none 
    ! input variables 
    real(REAL_TYPE), intent(in), dimension(3,3) :: m1
    real(REAL_TYPE), intent(in), dimension(3) :: v1
    ! output variables 
    real(REAL_TYPE), intent(out), dimension(3) :: DoubleDot4_2_result!m1Timesv1
    ! local variables 
    integer(INTEGER_TYPE) :: ii
    integer(INTEGER_TYPE) :: jj
    
    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_2 requires vector of size(3)!"
        return
    endif

    ! Check matrix size
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_2 requires 3-by-3 matrix "
        return
    endif
    
    
    ! Compute the result
    do ii = 1, 3
        DoubleDot4_2_result(ii) = 0.0
        do jj = 1, 3
            DoubleDot4_2_result(ii) = DoubleDot4_2_result(ii) + ( m1(ii, jj) * v1(jj) )
        enddo
    enddo
    
    end subroutine DoubleDot4_2
    
    
    
    
    
    
    
    
    
    
    subroutine DoubleDot2_4(v1, m1, DoubleDot2_4_result)
    !// computes doubledot product for matrix-vector arguments
    !// caution: first coordinate of the matrix should be in opposite 
    !// variant form of vector
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3, 3), intent(in) :: m1
    real(REAL_TYPE), dimension(3), intent(out) :: DoubleDot2_4_result
    ! Local variables
    integer(INTEGER_TYPE) :: ii, jj

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_4 requires vector of size(3)!"
        return
    endif

    ! Check matrix size
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_4 requires 3-by-3 matrix "
        return
    endif

    ! Compute the result
    do ii = 1, 3
        DoubleDot2_4_result(ii) = 0.0
        do jj = 1, 3
            DoubleDot2_4_result(ii) = DoubleDot2_4_result(ii) + ( m1(jj, ii) * v1(jj) )
        enddo
    enddo

    end subroutine DoubleDot2_4
    
    
    
    
    
    
    
    subroutine DoubleDot4_4(m1, m2, DoubleDot4_4_result)
    !// computes doubledot product for matrix-matrix arguments
    !// caution: second coordinate of the first matrix should be in opposite 
    !// variant form of the first coordinate of second matrix
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3, 3), intent(in) :: m1, m2
    real(REAL_TYPE), dimension(3, 3), intent(out) :: DoubleDot4_4_result
    ! Local variables
    integer(INTEGER_TYPE) :: ii, jj, kk

    ! Check matrix sizes
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3 .or. &
        size(m2, 1) /= 3 .or. size(m2, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_4 requires 3-by-3 matrices"
        return
    endif

    ! Compute the result
    do ii = 1, 3
        do jj = 1, 3
            DoubleDot4_4_result(ii, jj) = 0.0
            do kk = 1, 3
                DoubleDot4_4_result(ii, jj) = DoubleDot4_4_result(ii, jj) + ( m1(ii, kk) * m2(kk, jj) )
            enddo
        enddo
    enddo
    
    
    
    end subroutine DoubleDot4_4
    
    
    
    
    
    subroutine ToContraviant(v1, ToContraviant_result)
    !/*************************************************************/
    !// ToContraviant() ---------------------------------------------
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(out) :: ToContraviant_result
    ! Local variables
    integer :: ii
    
    ! initialize
    ToContraviant_result = 0.0

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::ToContraviant requires vector of size(3)!"
        return
    endif

    ! Copy input vector to result
    ToContraviant_result = v1

    ! Adjust second component of the result
    ToContraviant_result(2) = ToContraviant_result(2) * 0.5

    end subroutine ToContraviant
    
    
    
    
    
    subroutine ToCovariant(v1, ToCovariant_result)
    !/*************************************************************/
    !// ToCovariant() ---------------------------------------------
    implicit none
    ! Arguments
    real, dimension(3), intent(in) :: v1
    real, dimension(3), intent(out) :: ToCovariant_result
    ! Local variables
    integer :: ii

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::ToCovariant requires vector of size(3)!"
        return
    endif

    ! Copy input vector to result
    ToCovariant_result = v1

    ! Adjust second component of the result
    ToCovariant_result(2) = ToCovariant_result(2) * 2.0


    end subroutine ToCovariant
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    end module 