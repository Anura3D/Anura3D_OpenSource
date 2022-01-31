!===============================================================================
! Copyright 2004-2019 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!   Content:
!           Intel(R) Math Kernel Library (Intel(R) MKL) DSS Fortran-90 header file
!
!           Contains main datatypes, routines and constants definition.
!           For CDECL use only.
!
!*******************************************************************************
!DEC$ IF .NOT. DEFINED( __MKL_DSS_F90 )

!DEC$ DEFINE __MKL_DSS_F90


      MODULE MKL_DSS_PRIVATE
        TYPE MKL_DSS_HANDLE; INTEGER(KIND=8) DUMMY; END TYPE
      END MODULE MKL_DSS_PRIVATE

      MODULE MKL_DSS

        USE MKL_DSS_PRIVATE

      INTEGER, PARAMETER :: MKL_DSS_DEFAULTS = 0

!
! Out-of-core level option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_OOC_VARIABLE = 1024
      INTEGER, PARAMETER :: MKL_DSS_OOC_STRONG   = 2048

!
! Refinement steps on / off
!
      INTEGER, PARAMETER :: MKL_DSS_REFINEMENT_OFF = 4096
      INTEGER, PARAMETER :: MKL_DSS_REFINEMENT_ON  = 8192

!
! Solver step's substitution
!

      INTEGER, PARAMETER :: MKL_DSS_FORWARD_SOLVE	=	16384
      INTEGER, PARAMETER :: MKL_DSS_DIAGONAL_SOLVE	=	32768
      INTEGER, PARAMETER :: MKL_DSS_BACKWARD_SOLVE	=	49152
      INTEGER, PARAMETER :: MKL_DSS_TRANSPOSE_SOLVE	=	262144
      INTEGER, PARAMETER :: MKL_DSS_CONJUGATE_SOLVE	=	524288

!
! Single precision
!

      INTEGER, PARAMETER :: MKL_DSS_SINGLE_PRECISION	=	65536

!
! Zero-based indexing
!
      INTEGER, PARAMETER :: MKL_DSS_ZERO_BASED_INDEXING =	131072

!
! Message level option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_SUCCESS = -2147483647
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_DEBUG   = -2147483646
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_INFO    = -2147483645
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_WARNING = -2147483644
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_ERROR   = -2147483643
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_FATAL   = -2147483642

!
! Termination level option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_SUCCESS = 1073741832
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_DEBUG   = 1073741840
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_INFO    = 1073741848
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_WARNING = 1073741856
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_ERROR   = 1073741864
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_FATAL   = 1073741872

!
! Structure option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_SYMMETRIC                   = 536870976
      INTEGER, PARAMETER :: MKL_DSS_SYMMETRIC_STRUCTURE         = 536871040
      INTEGER, PARAMETER :: MKL_DSS_NON_SYMMETRIC               = 536871104
      INTEGER, PARAMETER :: MKL_DSS_SYMMETRIC_COMPLEX           = 536871168
      INTEGER, PARAMETER :: MKL_DSS_SYMMETRIC_STRUCTURE_COMPLEX = 536871232
      INTEGER, PARAMETER :: MKL_DSS_NON_SYMMETRIC_COMPLEX       = 536871296
!
! Reordering option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_AUTO_ORDER         = 268435520
      INTEGER, PARAMETER :: MKL_DSS_MY_ORDER           = 268435584
      INTEGER, PARAMETER :: MKL_DSS_OPTION1_ORDER      = 268435648
      INTEGER, PARAMETER :: MKL_DSS_GET_ORDER          = 268435712
      INTEGER, PARAMETER :: MKL_DSS_METIS_ORDER        = 268435776
      INTEGER, PARAMETER :: MKL_DSS_METIS_OPENMP_ORDER = 268435840

!
! Factorization option definitions
!

      INTEGER, PARAMETER :: MKL_DSS_POSITIVE_DEFINITE           = 134217792
      INTEGER, PARAMETER :: MKL_DSS_INDEFINITE                  = 134217856
      INTEGER, PARAMETER :: MKL_DSS_HERMITIAN_POSITIVE_DEFINITE = 134217920
      INTEGER, PARAMETER :: MKL_DSS_HERMITIAN_INDEFINITE        = 134217984

!
! Return status values
!

      INTEGER, PARAMETER :: MKL_DSS_SUCCESS         = 0
      INTEGER, PARAMETER :: MKL_DSS_ZERO_PIVOT      = -1
      INTEGER, PARAMETER :: MKL_DSS_OUT_OF_MEMORY   = -2
      INTEGER, PARAMETER :: MKL_DSS_FAILURE         = -3
      INTEGER, PARAMETER :: MKL_DSS_ROW_ERR         = -4
      INTEGER, PARAMETER :: MKL_DSS_COL_ERR         = -5
      INTEGER, PARAMETER :: MKL_DSS_TOO_FEW_VALUES  = -6
      INTEGER, PARAMETER :: MKL_DSS_TOO_MANY_VALUES = -7
      INTEGER, PARAMETER :: MKL_DSS_NOT_SQUARE      = -8
      INTEGER, PARAMETER :: MKL_DSS_STATE_ERR       = -9
      INTEGER, PARAMETER :: MKL_DSS_INVALID_OPTION  = -10
      INTEGER, PARAMETER :: MKL_DSS_OPTION_CONFLICT = -11
      INTEGER, PARAMETER :: MKL_DSS_MSG_LVL_ERR     = -12
      INTEGER, PARAMETER :: MKL_DSS_TERM_LVL_ERR    = -13
      INTEGER, PARAMETER :: MKL_DSS_STRUCTURE_ERR   = -14
      INTEGER, PARAMETER :: MKL_DSS_REORDER_ERR     = -15
      INTEGER, PARAMETER :: MKL_DSS_VALUES_ERR      = -16
      INTEGER, PARAMETER :: MKL_DSS_STATISTICS_INVALID_MATRIX = -17
      INTEGER, PARAMETER :: MKL_DSS_STATISTICS_INVALID_STATE  = -18
      INTEGER, PARAMETER :: MKL_DSS_STATISTICS_INVALID_STRING = -19

!
! Function prototypes for DSS routines
!

      INTERFACE
        FUNCTION DSS_CREATE( HANDLE, OPT )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(OUT)   :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER                      :: DSS_CREATE
        END FUNCTION DSS_CREATE
      END INTERFACE

      INTERFACE
        FUNCTION DSS_DEFINE_STRUCTURE( HANDLE, OPT, ROWINDEX, NROWS, RCOLS, COLUMNS, NNONZEROS )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NROWS
          INTEGER,       INTENT(IN)    :: RCOLS
          INTEGER,       INTENT(IN)    :: NNONZEROS
          INTEGER,       INTENT(IN)    :: ROWINDEX( * ) ! * = MIN(NROWS, NCOLS)+1
          INTEGER,       INTENT(IN)    :: COLUMNS( * ) ! * = NNONZEROS
          INTEGER                      :: DSS_DEFINE_STRUCTURE
        END FUNCTION DSS_DEFINE_STRUCTURE
      END INTERFACE

      INTERFACE
        FUNCTION DSS_REORDER( HANDLE, OPT, PERM )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: PERM( * )
          INTEGER                      :: DSS_REORDER
        END FUNCTION DSS_REORDER
      END INTERFACE

      INTERFACE DSS_FACTOR
        FUNCTION DSS_FACTOR_REAL_D( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          REAL(KIND=8),          INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_REAL_D
        END FUNCTION DSS_FACTOR_REAL_D

        FUNCTION DSS_FACTOR_REAL_S( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          REAL(KIND=4),          INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_REAL_S
        END FUNCTION DSS_FACTOR_REAL_S

        FUNCTION DSS_FACTOR_COMPLEX_D( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          COMPLEX(KIND=8),       INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_COMPLEX_D
        END FUNCTION DSS_FACTOR_COMPLEX_D

        FUNCTION DSS_FACTOR_COMPLEX_S( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          COMPLEX(KIND=4),       INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_COMPLEX_S
        END FUNCTION DSS_FACTOR_COMPLEX_S
      END INTERFACE

      INTERFACE DSS_FACTOR_REAL
        FUNCTION DSS_FACTOR_REAL_D_( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          REAL(KIND=8),          INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_REAL_D_
        END FUNCTION DSS_FACTOR_REAL_D_

        FUNCTION DSS_FACTOR_REAL_S_( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          REAL(KIND=4),          INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_REAL_S_
        END FUNCTION DSS_FACTOR_REAL_S_
      END INTERFACE

      INTERFACE DSS_FACTOR_COMPLEX
        FUNCTION DSS_FACTOR_COMPLEX_D_( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          COMPLEX(KIND=8),       INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_COMPLEX_D_
        END FUNCTION DSS_FACTOR_COMPLEX_D_

        FUNCTION DSS_FACTOR_COMPLEX_S_( HANDLE, OPT, RVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          COMPLEX(KIND=4),       INTENT(IN)    :: RVALUES( * )
          INTEGER                      :: DSS_FACTOR_COMPLEX_S_
        END FUNCTION DSS_FACTOR_COMPLEX_S_
      END INTERFACE

      INTERFACE DSS_SOLVE
        FUNCTION DSS_SOLVE_REAL_D( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          REAL(KIND=8),          INTENT(IN)    :: RRHSVALUES( * )
          REAL(KIND=8),          INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_REAL_D
        END FUNCTION DSS_SOLVE_REAL_D

        FUNCTION DSS_SOLVE_REAL_S( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          REAL(KIND=4),          INTENT(IN)    :: RRHSVALUES( * )
          REAL(KIND=4),          INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_REAL_S
        END FUNCTION DSS_SOLVE_REAL_S

        FUNCTION DSS_SOLVE_COMPLEX_D( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          COMPLEX(KIND=8),       INTENT(IN)    :: RRHSVALUES( * )
          COMPLEX(KIND=8),       INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_COMPLEX_D
        END FUNCTION DSS_SOLVE_COMPLEX_D

        FUNCTION DSS_SOLVE_COMPLEX_S( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          COMPLEX(KIND=4),       INTENT(IN)    :: RRHSVALUES( * )
          COMPLEX(KIND=4),       INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_COMPLEX_S
        END FUNCTION DSS_SOLVE_COMPLEX_S
      END INTERFACE

      INTERFACE DSS_SOLVE_REAL
        FUNCTION DSS_SOLVE_REAL_D_( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          REAL(KIND=8),          INTENT(IN)    :: RRHSVALUES( * )
          REAL(KIND=8),          INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_REAL_D_
        END FUNCTION DSS_SOLVE_REAL_D_

        FUNCTION DSS_SOLVE_REAL_S_( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          REAL(KIND=4),          INTENT(IN)    :: RRHSVALUES( * )
          REAL(KIND=4),          INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_REAL_S_
        END FUNCTION DSS_SOLVE_REAL_S_
      END INTERFACE

      INTERFACE DSS_SOLVE_COMPLEX
        FUNCTION DSS_SOLVE_COMPLEX_D_( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          COMPLEX(KIND=8),       INTENT(IN)    :: RRHSVALUES( * )
          COMPLEX(KIND=8),       INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_COMPLEX_D_
        END FUNCTION DSS_SOLVE_COMPLEX_D_

        FUNCTION DSS_SOLVE_COMPLEX_S_( HANDLE, OPT, RRHSVALUES, NRHS, RSOLVALUES )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(INOUT) :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: NRHS
          COMPLEX(KIND=4),       INTENT(IN)    :: RRHSVALUES( * )
          COMPLEX(KIND=4),       INTENT(OUT)   :: RSOLVALUES( * )
          INTEGER                      :: DSS_SOLVE_COMPLEX_S_
        END FUNCTION DSS_SOLVE_COMPLEX_S_
      END INTERFACE

      INTERFACE
        FUNCTION DSS_DELETE( HANDLE, OPT )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(IN)    :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER                      :: DSS_DELETE
        END FUNCTION DSS_DELETE
      END INTERFACE

      INTERFACE
        FUNCTION DSS_STATISTICS( HANDLE, OPT, STAT, RET )
          USE MKL_DSS_PRIVATE
          TYPE(MKL_DSS_HANDLE), INTENT(IN)    :: HANDLE
          INTEGER,       INTENT(IN)    :: OPT
          INTEGER,       INTENT(IN)    :: STAT( * )
          REAL(KIND=8),        INTENT(OUT)   :: RET( * )
          INTEGER                      :: DSS_STATISTICS
        END FUNCTION DSS_STATISTICS
      END INTERFACE

      INTERFACE
        SUBROUTINE MKL_CVT_TO_NULL_TERMINATED_STR(DESTSTR,DESTLEN,SRCSTR)
          INTEGER,       INTENT(OUT)   :: DESTSTR( * )
          INTEGER,       INTENT(IN)    :: DESTLEN
          CHARACTER,          INTENT(IN)    :: SRCSTR(*)
        END SUBROUTINE MKL_CVT_TO_NULL_TERMINATED_STR
      END INTERFACE

      END MODULE MKL_DSS


!DEC$ ENDIF
