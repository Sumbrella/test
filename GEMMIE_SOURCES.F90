MODULE  GLOBE_SOURCES
    USE PRECISION_CONTROL_MODULE
    USE GEMMIE_CONSTANTS
    USE SPHERICAL_HARMONICS
    USE RADIAL_DATA_TYPES_MODULE
    USE LAYERED_MEDIA_MODULE
    USE MPI_UTILITIES_MODULE
    USE SPHERICAL_BESSEL_SINGLE_INTEGRATION_MODULE
    USE RADIAL_INTEGRATION_CTLS
    USE CTL_TYPES_MODULE
    USE S2A_OPERATOR_MODULE
    USE S2O_OPERATOR_MODULE
    USE INTEGRAL_OPERATORS_MODULE
    USE CHECK_MEMORY
    USE TABLE_FUNCTIONS
    USE Lateral_Integration_Utilities
    USE MODIFIED_FUNCTIONS_MODULE


    IMPLICIT NONE
    PRIVATE
    PUBLIC :: EM_FIELDS_FOR_SPHERICAL_HARMONICS_1D
    PUBLIC :: EM_FIELDS_FOR_SPHERICAL_HARMONICS_1D_SITE_WISE
    PUBLIC :: CALCULATE_C1_RESPONSE, CALCULATE_Q_RESPONSES_MANY_FREQ
    PUBLIC :: CALCULATE_RHS_FOR_IE_BY_SH_COEFFS
    PUBLIC :: SCALARS_FOR_TIPPERS_S2A
    PUBLIC :: SCALARS_FOR_TIPPERS_S2O
    PUBLIC :: CALCULATE_RHS_FOR_IE_TIPPERS_SRC
    PUBLIC :: CALCULATE_EM_FIELDS_FOR_TIPPERS_SRC
CONTAINS

    FUNCTION CALCULATE_C1_RESPONSE(bkg, freqs, theta) RESULT (C)
        TYPE(BACKGROUND_TYPE),   INTENT(INOUT) :: bkg
        REAL(DP), INTENT(IN) :: freqs
        REAL(DP), INTENT(IN) :: theta
        COMPLEX(DP) :: C(3)
        COMPLEX(DP) :: E(E_THETA:E_R)
        COMPLEX(DP) :: H(H_THETA:H_R)
        COMPLEX(DP) :: e_coeffs(-1:1, 1)
        INTEGER :: Nr
        Nr=UBOUND(BKG%r,1)
        CALL SET_FREQUENCY(freqs, bkg)
        e_coeffs=C_ZERO
        e_coeffs(0,1)=C_ONE
        CALL EM_FIELDS_FOR_SPHERICAL_HARMONICS_1D(bkg, Nr, e_coeffs,&
            [theta], [R_ZERO], E, H) 
        C=C_ZERO
        C(1)=-BKG%r(Nr)*TAN(theta)/R_TWO*H(H_R)/H(H_THETA)
        C(2)=E(E_PHI)/H(H_THETA)
        !C(3)=E(E_THETA)/H(H_PHI)
    ENDFUNCTION

    SUBROUTINE CALCULATE_Q_RESPONSES_MANY_FREQ(bkg, freqs, Nmax, Q)
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        REAL(DP), INTENT(IN) :: freqs(:)
        INTEGER , INTENT(IN) :: Nmax
        COMPLEX(DP), INTENT(OUT) :: Q(:,:)
        INTEGER :: I
        TYPE(BACKGROUND_TYPE) :: w_bkg
        w_bkg=bkg
        DO I=1, SIZE(freqs)
            CALL CALCULATE_Q_RESPONSES_1D(w_bkg, freqs(I), Nmax, Q(:,I))
        ENDDO
    ENDSUBROUTINE


    SUBROUTINE CALCULATE_Q_RESPONSES_1D(bkg, freq, Nmax, Q) 
        TYPE(BACKGROUND_TYPE),   INTENT(INOUT) :: bkg
        REAL(DP), INTENT(IN) :: freq
        INTEGER , INTENT(IN) :: Nmax
        COMPLEX(DP), INTENT(OUT) :: Q(:)
        COMPLEX(DP) :: e_coeffs(0:0, 1:Nmax)
        COMPLEX(DP) :: t_coeffs(0:0, 1:Nmax, 1:3 )
        INTEGER :: Nr
        Nr=UBOUND(BKG%r,1)
        CALL SET_FREQUENCY(freq, bkg)
        e_coeffs=R_ONE
        t_coeffs=GET_TOTAL_EM_COEFFS(bkg, Nmax, 0,  Nr, e_coeffs )
        Q=t_coeffs(0,:,2)-C_ONE
    ENDSUBROUTINE

    SUBROUTINE EM_FIELDS_FOR_SPHERICAL_HARMONICS_1D(bkg,  Nr, external_coeffs, theta_s, phi_s, E, H) 
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        INTEGER, INTENT(IN) ::   Nr
        COMPLEX(DP), INTENT(IN) :: external_coeffs(:, :)
        REAL(DP) :: theta_s(:), phi_s(:)
        COMPLEX(DP), INTENT(OUT) :: E( SIZE(theta_s), E_THETA:E_R, SIZE(phi_s) )
        COMPLEX(DP), INTENT(OUT) :: H( SIZE(theta_s), H_THETA:H_R, SIZE(phi_s) )

        COMPLEX(DP), ALLOCATABLE :: tot_coeffs(:,:,: )
        INTEGER:: N
        N=SIZE(external_coeffs,2)
        ALLOCATE(tot_coeffs(-N: N, 1:N, 1:3 ))
        tot_coeffs=GET_TOTAL_EM_COEFFS(bkg, N, N,  Nr, external_coeffs )
        IF (SIZE(E)>SIZE(theta_s)*N**2) THEN
            CALL GET_FIELDS_BY_COEFFICIENTS(N,theta_s, phi_s, tot_coeffs, E, H)
        ELSE
            CALL GET_FIELDS_BY_COEFFICIENTS_MEM_OPT(N, theta_s, phi_s, tot_coeffs, E, H) 
        ENDIF

        H=H*bkg%MainSphereRadius 
    ENDSUBROUTINE

    SUBROUTINE EM_FIELDS_FOR_SPHERICAL_HARMONICS_1D_SITE_WISE(bkg,  Nr, external_coeffs, theta_s, phi_s, E, H) 
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        INTEGER, INTENT(IN) ::   Nr
        COMPLEX(DP), INTENT(IN) :: external_coeffs(:, :)
        REAL(DP) :: theta_s(:), phi_s(:)
        COMPLEX(DP), INTENT(OUT) :: E( SIZE(theta_s), E_THETA:E_R )
        COMPLEX(DP), INTENT(OUT) :: H( SIZE(theta_s), H_THETA:H_R )

        COMPLEX(DP), ALLOCATABLE :: tot_coeffs(:,:,: )
        INTEGER:: N
        N=SIZE(external_coeffs,2)
        ALLOCATE(tot_coeffs(-N: N, 1:N, 1:3 ))
        tot_coeffs=GET_TOTAL_EM_COEFFS(bkg, N, N,  Nr, external_coeffs )
        CALL GET_FIELDS_BY_COEFFICIENTS_SITE_WISE(N,theta_s, phi_s, tot_coeffs, E, H)
        H=H*bkg%MainSphereRadius 
    ENDSUBROUTINE

    SUBROUTINE CALCULATE_RHS_FOR_IE_BY_SH_COEFFS(bkg, Anomaly,   external_coeffs,  E_bkg_sph)
        USE GAUSS_INTEGRATION_MODULE
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        TYPE(SPH_DOMAIN_SHAPE_TYPE), INTENT(IN) :: Anomaly
        COMPLEX(DP), INTENT(IN) :: external_coeffs(:, :)
        COMPLEX(DP), INTENT(OUT) :: E_bkg_sph(Anomaly%Nr, Anomaly%Ntheta,&
            &E_THETA:E_R,Anomaly%Nphi_loc(1):Anomaly%Nphi_loc(2))

        COMPLEX(DP), ALLOCATABLE :: tot_coeffs(:,:,:,: )
        INTEGER ::   N
        N=SIZE(external_coeffs,2)
        ALLOCATE(tot_coeffs(Anomaly%Nr, -N: N, 1:N, 1:3 ))
        tot_coeffs=GET_RADIAL_AVERAGED_E_COEFFS(bkg, Anomaly, N, external_coeffs)
        CALL GET_E_FIELD_LAT_INTEGRALS(Anomaly%Nr,  N, Anomaly%dphi, Anomaly%phi0, &
            Anomaly%Nphi_loc, Anomaly%thetas,tot_coeffs, E_bkg_sph)

    ENDSUBROUTINE


    FUNCTION GET_TOTAL_EM_COEFFS(bkg, N, Nm, Nr, external_coeffs) RESULT (RES)
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        INTEGER, INTENT(IN) :: N, Nm,  Nr
        COMPLEX(DP), INTENT(IN) :: external_coeffs(-Nm:Nm, 1:N)
        COMPLEX(DP) :: RES(-Nm:Nm, 1:N, 1:3 )

        TYPE(SPH_DOMAIN_SHAPE_TYPE) :: Domain_shape
        TYPE(RADIAL_OPERATOR_TYPE)::Domain
        TYPE(MPI_STATE_TYPE) :: mpi
        COMPLEX(DP) :: EM_recv_part(3)
        COMPLEX(DP) :: TMP(RC_E_LP:RC_H_RL)
        INTEGER :: I_n, I_m,  Ns, Nr2, J_n
        REAL(DP) :: S_model

        Ns=UBOUND(bkg%r,1)
        S_model=bkg%r(Ns)**2
        IF (bkg%normalized) S_model=S_model*bkg%MainSphereRadius**2

        Domain_shape%Mr1=Nr
        Domain_shape%Mr2=UBOUND(bkg%r,1)+1
        Domain_shape%Nr=Domain_shape%Mr2-Domain_shape%Mr1

        Domain=CREATE_RADIAL_DOMAIN(Domain_shape, bkg, mpi,  N)
        Nr2=1
        EM_recv_part=C_ZERO
        DO I_n=1, N
            TMP= SITE_TO_SITE_RADIAL_D2O_SINGLE_TERM(Domain, Domain_shape%Nr, Nr2,&
                I_n, SOURCE_ABOVE_RECIEVER)

            EM_recv_part(1)=TMP(RC_E_LT)*I_n*PI4*S_model
            EM_recv_part(2)=TMP(RC_H_LT)*I_n*PI4*S_model
            EM_recv_part(3)=TMP(RC_H_RL)*I_n*PI4*S_model
            J_n=MIN(Nm, I_n)
            DO I_m=-J_n, J_n
                RES(I_m, I_n,:)=external_coeffs(I_m, I_n)*EM_recv_part
            ENDDO

        ENDDO
    ENDFUNCTION

    FUNCTION GET_RADIAL_AVERAGED_E_COEFFS(bkg, Anomaly, N,  external_coeffs) RESULT (RES)
        TYPE(BACKGROUND_TYPE),   INTENT(IN) :: bkg
        TYPE(SPH_DOMAIN_SHAPE_TYPE), INTENT(IN) :: Anomaly
        INTEGER, INTENT(IN) :: N
        COMPLEX(DP), INTENT(IN) :: external_coeffs(-N:N, 1:N)
        COMPLEX(DP) :: RES(Anomaly%Nr, -N:N, 1:N, 3)

        TYPE(SPH_DOMAIN_SHAPE_TYPE) :: Domain_shape
        TYPE(BACKGROUND_TYPE)   :: bkg2
        TYPE(RADIAL_OPERATOR_TYPE)::Domain
        TYPE(MPI_STATE_TYPE)     :: mpi
        COMPLEX(DP) :: TMP(RC_E_LP: RC_H_RL)
        COMPLEX(DP) :: A(4,3,1), Et
        REAL(DP) :: S_model
        INTEGER :: I_n,  J,  Ns, Ns2
        Ns=UBOUND(bkg%r,1)

        S_model=bkg%r(Ns)**2
        IF (bkg%normalized) S_model=S_model*bkg%MainSphereRadius**2
        CALL ADD_LAYER_TO_DOMAIN(bkg%r(Ns), Anomaly, bkg, Domain_shape, bkg2)

        Domain=CREATE_RADIAL_DOMAIN(Domain_shape, bkg2, mpi,  N)
        Ns2=GET_OBS_LAYER_INDEX(bkg%r(Ns), Domain_shape%r)

        RES=C_ZERO
        DO I_n=1, N

            DO J=1, Anomaly%Nr

                ASSOCIATE(z=>Domain%Layers(J)%Spherical_Layer%z2,&
                        &s=>Domain%Layers(J)%Spherical_Layer%s,&
                        &r=>Domain%Layers(J)%Spherical_Layer%r2)
                    A =  GET_INTEGRATION_COEFFS_MAIN(z, s, [I_n, I_n])
                    TMP= INTEGRATE_RADIAL_D2O_SINGLE_TERM(domain, J,  Ns2, I_n, A(:,:,1),&
                        RECIEVER_ABOVE_SOURCE ) !We use reciprocity ))
                    Et=TMP(RC_E_LT)*I_n*PI4*S_model
                END ASSOCIATE
                RES(J, -I_n:I_n, I_n, 1 )=external_coeffs(-I_n:I_n, I_n)*Et/&
                    GET_LAYER_VR(Domain%Layers(J)%Spherical_Layer)
            ENDDO
        ENDDO
    ENDFUNCTION

    SUBROUTINE GET_FIELDS_BY_COEFFICIENTS(N,theta_s, phi_s, coeffs, E, H) 
        INTEGER, INTENT(IN) :: N
        REAL(DP), INTENT(IN) :: theta_s(:)
        REAL(DP), INTENT(IN) :: phi_s(:)
        COMPLEX(DP), INTENT(IN) :: coeffs(-N:N, 1:N, 1:3 )
        COMPLEX(DP), INTENT(OUT) :: E(SIZE(theta_s), E_THETA:E_R, SIZE(phi_s)  )
        COMPLEX(DP), OPTIONAL, INTENT(OUT) :: H(SIZE(theta_s), H_THETA:H_R, SIZE(phi_s)  )
        REAL(DP), ALLOCATABLE :: PM(:, :, :)
        REAL(DP), ALLOCATABLE :: PD(:, :, :)
        REAL(DP) :: COS_THETAS(SIZE(theta_s)), SIN_THETAS(SIZE(theta_s))
        COMPLEX(DP) :: EXP_PHI(1:N)
        INTEGER :: I_n, I_m
        INTEGER :: I_theta, I_phi

        ALLOCATE(PM(SIZE(theta_s),0:N, 0:N))
        ALLOCATE(PD(SIZE(theta_s), 0:N, 0:N))
        COS_THETAS=COS(theta_s)
        SIN_THETAS=SIN(theta_s)
        CALL GET_ASSOCIATED_LEGENDRE_POLYNOMAILS_FOR_ANGLES(SIZE(theta_s),&
            N, COS_THETAS, SIN_THETAS,  PM, PD)
        E=C_ZERO
        IF (PRESENT(H)) H=C_ZERO
        DO I_phi=1, SIZE(phi_s)
            DO I_m=1, N
                EXP_PHI(I_m)=EXP(C_IONE*I_m*phi_s(I_phi))
            ENDDO
            DO I_theta=1, SIZE(theta_s)
                DO I_n=1,N

                    E(I_theta, E_PHI, I_phi)=&
                        E(I_theta, E_PHI, I_phi)+coeffs(0,I_n,1)*PD(I_theta, 0,I_n)

                    DO I_m=1, I_n

                        E(I_theta, E_THETA, I_phi)=E(I_theta, E_THETA, I_phi)&
                            -C_IONE*I_m*coeffs(I_m,I_n,1)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                            +C_IONE*I_m*coeffs(-I_m,I_n,1)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                        E(I_theta, E_PHI, I_phi)=E(I_theta, E_PHI, I_phi)&
                            +coeffs(I_m,I_n,1)*PD(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                            +coeffs(-I_m,I_n,1)*PD(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                    ENDDO
                ENDDO
            ENDDO
            IF (PRESENT(H)) THEN
                DO I_theta=1, SIZE(theta_s)
                    DO I_n=1,N


                        H(I_theta, H_THETA, I_phi)=&
                            H(I_theta, H_THETA, I_phi)-coeffs(0,I_n,2)*PD(I_theta, 0,I_n)

                        H(I_theta, H_R, I_phi)=&
                            H(I_theta, H_R, I_phi)-coeffs(0,I_n,3)*PM(I_theta, 0,I_n)
                        DO I_m=1, I_n


                            H(I_theta, H_THETA, I_phi)=H(I_theta, H_THETA, I_phi)&
                                -coeffs(I_m,I_n,2)*PD(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                                -coeffs(-I_m,I_n,2)*PD(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                            H(I_theta, H_PHI, I_phi)=H(I_theta, H_PHI, I_phi)&
                                -C_IONE*I_m*coeffs(I_m,I_n,2)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                                +C_IONE*I_m*coeffs(-I_m,I_n,2)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                            H(I_theta, H_R, I_phi)= H(I_theta, H_R, I_phi)&
                                -coeffs(I_m,I_n,3)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                                -coeffs(-I_m,I_n,3)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                        ENDDO
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
    ENDSUBROUTINE

    SUBROUTINE GET_FIELDS_BY_COEFFICIENTS_MEM_OPT(N,theta_s, phi_s, coeffs, E, H) 
        INTEGER, INTENT(IN) :: N
        REAL(DP), INTENT(IN) :: theta_s(:)
        REAL(DP), INTENT(IN) :: phi_s(:)
        COMPLEX(DP), INTENT(IN) :: coeffs(-N:N, 1:N, 1:3 )
        COMPLEX(DP), INTENT(OUT) :: E(SIZE(theta_s), E_THETA:E_R, SIZE(phi_s)  )
        COMPLEX(DP), OPTIONAL, INTENT(OUT) :: H(SIZE(theta_s), H_THETA:H_R, SIZE(phi_s)  )
        REAL(DP), ALLOCATABLE :: PM(:, :)
        REAL(DP), ALLOCATABLE :: PD(:, :)
        REAL(DP) :: COS_THETAS(SIZE(theta_s)), SIN_THETAS(SIZE(theta_s))
        COMPLEX(DP) :: EXP_PHI(1:N)
        INTEGER :: I_n, I_m
        INTEGER :: I_theta, I_phi

        COS_THETAS=COS(theta_s)
        SIN_THETAS=SIN(theta_s)
        ALLOCATE(PM(0:N, 0:N), PD(0:N, 0:N))
        E=C_ZERO
        IF (PRESENT(H)) H=C_ZERO
        DO I_phi=1, SIZE(phi_s)
            DO I_m=1, N
                EXP_PHI(I_m)=EXP(C_IONE*I_m*phi_s(I_phi))
            ENDDO
            DO I_theta=1, SIZE(theta_s)
                CALL GET_ASSOCIATED_LEGENDRE_POLYNOMAILS_FOR_ANGLES(&
                    N, COS_THETAS(I_theta), SIN_THETAS(I_theta),  PM, PD)
                DO I_n=1,N

                    E(I_theta, E_PHI, I_phi)=&
                        E(I_theta, E_PHI, I_phi)+coeffs(0,I_n,1)*PD(0,I_n)
                    IF (PRESENT(H)) THEN

                        H(I_theta, H_THETA, I_phi)=&
                            H(I_theta, H_THETA, I_phi)-coeffs(0,I_n,2)*PD(0, I_n)

                        H(I_theta, H_R, I_phi)=&
                            H(I_theta, H_R, I_phi)-coeffs(0,I_n,3)*PM(0, I_n)
                    ENDIF

                    DO I_m=1, I_n

                        E(I_theta, E_THETA, I_phi)=E(I_theta, E_THETA, I_phi)&
                            -C_IONE*I_m*coeffs(I_m, I_n, 1)*PM(I_m, I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                            +C_IONE*I_m*coeffs(-I_m ,I_n, 1)*PM(I_m, I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                        E(I_theta, E_PHI, I_phi)=E(I_theta, E_PHI, I_phi)&
                            +coeffs(I_m, I_n, 1)*PD(I_m, I_n)*EXP_PHI(I_m)&
                            +coeffs(-I_m ,I_n, 1)*PD(I_m, I_n)*CONJG(EXP_PHI(I_m))

                        IF (PRESENT(H)) THEN

                            H(I_theta, H_THETA, I_phi)=H(I_theta, H_THETA, I_phi)&
                                -coeffs(I_m,I_n,2)*PD(I_m, I_n)*EXP_PHI(I_m)&
                                -coeffs(-I_m,I_n,2)*PD(I_m, I_n)*CONJG(EXP_PHI(I_m))

                            H(I_theta, H_PHI, I_phi)=H(I_theta, H_PHI, I_phi)&
                                -C_IONE*I_m*coeffs(I_m,I_n,2)*PM(I_m, I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                                +C_IONE*I_m*coeffs(-I_m,I_n,2)*PM(I_m, I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                            H(I_theta, H_R, I_phi)= H(I_theta, H_R, I_phi)&
                                -coeffs(I_m,I_n,3)*PM(I_m, I_n)*EXP_PHI(I_m)&
                                -coeffs(-I_m,I_n,3)*PM(I_m, I_n)*CONJG(EXP_PHI(I_m))

                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDSUBROUTINE

    SUBROUTINE GET_FIELDS_BY_COEFFICIENTS_SITE_WISE(N,theta_s, phi_s, coeffs, E, H) 
        INTEGER, INTENT(IN) :: N
        REAL(DP), INTENT(IN) :: theta_s(:)
        REAL(DP), INTENT(IN) :: phi_s(:)
        COMPLEX(DP), INTENT(IN) :: coeffs(-N:N, 1:N, 1:3 )
        COMPLEX(DP), INTENT(OUT) :: E(SIZE(theta_s), E_THETA:E_R  )
        COMPLEX(DP),  INTENT(OUT) :: H(SIZE(theta_s), H_THETA:H_R )
        REAL(DP) :: PM(SIZE(theta_s),0:N, 0:N)
        REAL(DP) :: PD(SIZE(theta_s), 0:N, 0:N)
        REAL(DP) :: COS_THETAS(SIZE(theta_s)), SIN_THETAS(SIZE(theta_s))
        COMPLEX(DP) :: EXP_PHI(1:N)
        INTEGER :: I_n, I_m
        INTEGER :: I_theta, I_phi

        COS_THETAS=COS(theta_s)
        SIN_THETAS=SIN(theta_s)
        CALL GET_ASSOCIATED_LEGENDRE_POLYNOMAILS_FOR_ANGLES(SIZE(theta_s),&
            N, COS_THETAS, SIN_THETAS,  PM, PD)
        E=C_ZERO
        H=C_ZERO
        !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(I_phi, I_theta, I_m, I_n, &
        !$OMP & EXP_PHI)
        !$OMP DO SCHEDULE (DYNAMIC)
        DO I_phi=1, SIZE(phi_s)
            DO I_m=1, N
                EXP_PHI(I_m)=EXP(C_IONE*I_m*phi_s(I_phi))
            ENDDO
            I_theta=I_phi
            DO I_n=1,N

                E(I_theta, E_PHI)=&
                    E(I_theta, E_PHI)+coeffs(0,I_n,1)*PD(I_theta, 0,I_n)

                H(I_theta, H_THETA)=&
                    H(I_theta, H_THETA)-coeffs(0,I_n,2)*PD(I_theta, 0,I_n)

                H(I_theta, H_R)=&
                    H(I_theta, H_R)-coeffs(0,I_n,3)*PM(I_theta, 0,I_n)

                DO I_m=1, I_n

                    E(I_theta, E_THETA)=E(I_theta, E_THETA)&
                        -C_IONE*I_m*coeffs(I_m,I_n,1)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                        +C_IONE*I_m*coeffs(-I_m,I_n,1)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                    E(I_theta, E_PHI)=E(I_theta, E_PHI)&
                        +coeffs(I_m,I_n,1)*PD(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                        +coeffs(-I_m,I_n,1)*PD(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                    H(I_theta, H_THETA)=H(I_theta, H_THETA)&
                        -coeffs(I_m,I_n,2)*PD(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                        -coeffs(-I_m,I_n,2)*PD(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                    H(I_theta, H_PHI)=H(I_theta, H_PHI)&
                        -C_IONE*I_m*coeffs(I_m,I_n,2)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)/SIN_THETAS(I_theta)&
                        +C_IONE*I_m*coeffs(-I_m,I_n,2)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))/SIN_THETAS(I_theta)

                    H(I_theta, H_R)= H(I_theta, H_R)&
                        -coeffs(I_m,I_n,3)*PM(I_theta, I_m,I_n)*EXP_PHI(I_m)&
                        -coeffs(-I_m,I_n,3)*PM(I_theta, I_m,I_n)*CONJG(EXP_PHI(I_m))

                ENDDO
            ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP ENDPARALLEL
    ENDSUBROUTINE


    SUBROUTINE GET_E_FIELD_LAT_INTEGRALS(Nr,  N, dphi, start_phi,  Nphi, theta_s, coeffs, E)
        USE GAUSS_INTEGRATION_MODULE 
        INTEGER, INTENT(IN) :: Nr, N
        REAL(DP), INTENT(IN) :: theta_s(0:)
        REAL(DP), INTENT(IN) :: dphi, start_phi
        INTEGER, INTENT(IN) :: Nphi(2)
        COMPLEX(DP), INTENT(IN) :: coeffs(Nr, -N:N, 1:N, 1:3 )
        COMPLEX(DP), INTENT(OUT) :: E(Nr,SIZE(theta_s)-1 , E_THETA:E_R, Nphi(1):Nphi(2)  )

        INTEGER, PARAMETER :: Ng=16
        REAL(DP), ALLOCATABLE :: PM(:, :, :)
        REAL(DP), ALLOCATABLE :: PD(:, :, :)
        REAL(DP) :: COS_THETAS(Ng), SIN_THETAS(Ng)
        REAL(DP) :: W(Ng), x(Ng)
        REAL(DP), ALLOCATABLE :: PM_INT(:, :)
        REAL(DP), ALLOCATABLE :: PD_INT(:, :)
        COMPLEX(DP) :: EXP_PHI(1:N, Nphi(1):Nphi(2))
        REAL(DP) :: dphi2, phi0
        REAL(DP) :: Stheta
        INTEGER :: I_n, I_m
        INTEGER :: I_theta, I_phi, J_n

        ALLOCATE(PM(Ng, 0:N, 0:N))
        ALLOCATE(PD(Ng, 0:N, 0:N))
        ALLOCATE(PM_INT(0:N, 0:N))
        ALLOCATE(PD_INT(0:N, 0:N))
        E=C_ZERO
        dphi2=dphi/R_TWO
        DO I_phi=Nphi(1), Nphi(2)
            phi0=dphi*(I_phi-0.5_DP)+start_phi
            DO I_m=1, N
                EXP_PHI(I_m, I_phi)=R_TWO*SIN(dphi2*I_m)*&
                    EXP(C_IONE*I_m*phi0)/I_m/dphi
            ENDDO

        ENDDO

        DO I_theta=1, SIZE(theta_s)-1
            CALL GET_GAUSS_WEIGHTS_NODES(theta_s(I_theta-1),&
                theta_s(I_theta), x, W)
            COS_THETAS=COS(x)
            SIN_THETAS=SIN(x)

            Stheta=R_TWO*SIN((theta_s(I_theta)-theta_s(I_theta-1))/R_TWO)*&
                SIN((theta_s(I_theta)+theta_s(I_theta-1))/R_TWO)

            CALL GET_ASSOCIATED_LEGENDRE_POLYNOMAILS_FOR_ANGLES(Ng,&
                N, COS_THETAS, SIN_THETAS,  PM, PD)
            DO I_n=0, N
                DO J_n=0, N
                    PM_INT(J_n, I_n)=SUM(PM(:,J_n, I_n)*W)/Stheta
                    PD_INT(J_n, I_n)=SUM(PD(:,J_n, I_n)*W*SIN_THETAS)/Stheta
                ENDDO
            ENDDO


            DO I_phi= Nphi(1), Nphi(2)
                DO I_n=1,N
                    E(:, I_theta, E_PHI, I_phi)=&
                        E(:, I_theta, E_PHI, I_phi)+coeffs(:,0,I_n,1)*PD_INT(0, I_n)

                    DO I_m=1, I_n

                        E(:, I_theta, E_THETA, I_phi)=E(:, I_theta, E_THETA, I_phi)&
                            -C_IONE*I_m*coeffs(:,I_m,I_n,1)*PM_INT(I_m, I_n)*EXP_PHI(I_m, I_phi)&
                            +C_IONE*I_m*coeffs(:,-I_m,I_n,1)*PM_INT(I_m, I_n)*CONJG(EXP_PHI(I_m, I_phi))

                        E(:, I_theta,E_PHI,  I_phi)=E(:, I_theta, E_PHI, I_phi)&
                            +coeffs(:, I_m,I_n,1)*PD_INT(I_m, I_n)*EXP_PHI(I_m, I_phi)&
                            +coeffs(:, -I_m,I_n,1)*PD_INT(I_m, I_n)*CONJG(EXP_PHI(I_m, I_phi))

                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDSUBROUTINE

    SUBROUTINE CALCULATE_RHS_FOR_IE_TIPPERS_SRC(run_data, I_src, SCALARS, E_bkg_sph)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_src
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(IN) :: SCALARS
        COMPLEX(DP), INTENT(OUT) :: E_bkg_sph(:,:,:, &
            run_data%SPH_Anomaly%Nphi_loc(1):)

        REAL(DP) :: theta_p, phi_p
        COMPLEX(DP)::v(1, SIZE(E_bkg_sph,1))
        COMPLEX(DP) :: E_t(SIZE(E_bkg_sph,1))
        COMPLEX(DP) :: E_p(SIZE(E_bkg_sph,1))

        INTEGER, PARAMETER :: Ng=8
        REAL(DP) :: thetas(Ng), phis(Ng), Wt(Ng), Wp(Ng), W0, W1
        INTEGER :: Itheta, Iphi
        INTEGER :: Jphi, Jtheta
        INTEGER :: II(1)
        REAL(DP) :: sin_alpha2, sin_alpha(1), cos_alpha2
        REAL(DP) :: sin_beta, cos_beta
        REAL(DP) :: dtheta, stheta, func_corr

        theta_p=DEG2RAD(run_data%Sources(I_src)%Currents_Data%pole_theta)
        phi_p=DEG2RAD(run_data%Sources(I_src)%Currents_Data%pole_phi)

        ASSOCIATE(Nphi_loc=>run_data%SPH_Anomaly%Nphi_loc, Ntheta=>run_data%SPH_Anomaly%Ntheta, &
                theta_r=>run_data%SPH_Anomaly%thetas, phi_r=>run_data%SPH_Anomaly%phis, &
                dphi=>run_data%SPH_Anomaly%dphi)
            E_bkg_sph=C_ZERO
            DO Iphi=Nphi_loc(1), Nphi_loc(2)
                CALL GET_GAUSS_WEIGHTS_NODES(phi_r(Iphi-1),&
                    phi_r(Iphi), phis, Wp)
                DO Itheta=1, run_data%SPH_Anomaly%Ntheta
                    CALL GET_GAUSS_WEIGHTS_NODES(theta_r(Itheta-1),&
                        theta_r(Itheta), thetas, Wt)
                    Wt=Wt*SIN(thetas)
                    dtheta=(theta_r(Itheta)-theta_r(Itheta-1))/R_TWO
                    stheta=(theta_r(Itheta)+theta_r(Itheta-1))/R_TWO
                    W0=dphi*R_TWO*SIN(dtheta)*SIN(stheta)
                    DO Jphi=1, Ng
                        DO Jtheta=1,Ng
                            W1=Wt(Jtheta)*Wp(Jphi)/W0

                            CALL GET_SIN2_COS2(theta_p, thetas(Jtheta),&
                                phi_p, phis(Jphi), sin_alpha2, cos_alpha2)

                            func_corr=SIGN(R_ONE, cos_alpha2-sin_alpha2)
                            sin_alpha=SQRT(MIN(sin_alpha2, cos_alpha2))

                            II=GET_INTERVAL(sin_alpha, SCALARS%TABLE)
                            v=CALCULATE_DERIV_1_BY_TABLE(II, sin_alpha, SCALARS%TF(1))


                            CALL GET_ROTATION_ANGLE_SCALED(theta_p, phi_p, thetas(Jtheta),&
                                phis(Jphi), sin_beta, cos_beta)

                            E_t=v(1,:)*cos_beta*W1
                            E_p=v(1,:)*sin_beta*W1

                            E_bkg_sph(:, Itheta, E_THETA, Iphi)= &
                                E_bkg_sph(:, Itheta, E_THETA, Iphi)+E_t

                            E_bkg_sph(:, Itheta, E_PHI, Iphi)= &
                                E_bkg_sph(:, Itheta, E_PHI, Iphi)+E_p

                            v=CALCULATE_FUNC_BY_TABLE(II, sin_alpha, SCALARS%TF(2))
                            E_bkg_sph(:, Itheta, E_R, Iphi)= &
                                E_bkg_sph(:, Itheta, E_R, Iphi)+v(1,:)*func_corr*W1
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        END ASSOCIATE
    ENDSUBROUTINE

    SUBROUTINE CALCULATE_EM_FIELDS_FOR_TIPPERS_SRC(run_data, I_src, SCALARS, theta_recv, phi_recv, E, H)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_src
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(IN) :: SCALARS
        REAL(DP), INTENT(IN) :: theta_recv(:)
        REAL(DP), INTENT(IN) :: phi_recv(:)
        COMPLEX(DP), INTENT(OUT) :: E(SIZE(theta_recv), E_THETA:E_R, SIZE(phi_recv)  )
        COMPLEX(DP), INTENT(OUT) :: H(SIZE(theta_recv), H_THETA:H_R, SIZE(phi_recv)  )

        REAL(DP) :: theta_p, phi_p
        COMPLEX(DP)::v(1, 1)

        INTEGER :: Itheta, Iphi
        INTEGER :: II(1)
        REAL(DP) :: sin_alpha2, sin_alpha(1), cos_alpha2
        REAL(DP) :: sin_beta, cos_beta
        REAL(DP) ::  func_corr

        theta_p=DEG2RAD(run_data%Sources(I_src)%Currents_Data%pole_theta)
        phi_p=DEG2RAD(run_data%Sources(I_src)%Currents_Data%pole_phi)
        E=C_ZERO
        H=C_ZERO
        DO Iphi=1, SIZE(phi_recv)
            DO Itheta=1, SIZE(theta_recv)
                CALL GET_SIN2_COS2(theta_p, theta_recv(Itheta),&
                    phi_p, phi_recv(Iphi), sin_alpha2, cos_alpha2)

                func_corr=SIGN(R_ONE, cos_alpha2-sin_alpha2)
                sin_alpha=SQRT(MIN(sin_alpha2, cos_alpha2))

                CALL GET_ROTATION_ANGLE_SCALED(theta_p, phi_p, theta_recv(Itheta),&
                    phi_recv(Iphi), sin_beta, cos_beta)
   
                II=GET_INTERVAL(sin_alpha, SCALARS%TABLE)
                II(1)=MAX(II(1),0)

                v=CALCULATE_DERIV_1_BY_TABLE(II, sin_alpha, SCALARS%TF(1))

                E(Itheta, E_THETA, Iphi)=v(1,1)*cos_beta
                E(Itheta, E_PHI, Iphi)=v(1,1)*sin_beta

                v=CALCULATE_FUNC_BY_TABLE(II, sin_alpha, SCALARS%TF(2))
                E(Itheta, E_R, Iphi)=v(1,1)*func_corr

                v=CALCULATE_DERIV_1_BY_TABLE(II, sin_alpha, SCALARS%TF(3))

                H(Itheta, H_THETA, Iphi)=v(1,1)*sin_beta
                H(Itheta, H_PHI, Iphi)=-v(1,1)*cos_beta
            ENDDO
        ENDDO
    ENDSUBROUTINE

    SUBROUTINE GET_SIN2_COS2(theta1, theta2, phi1, phi2, sin_a2, cos_a2)
        REAL(DP), INTENT(IN) :: theta1, theta2, phi1, phi2
        REAL(DP), INTENT(OUT) :: sin_a2, cos_a2
        REAL(DP) :: sin_z2, ss
        sin_z2=SIN((phi1-phi2)/R_TWO)**2
        ss=SIN(theta1)*SIN(theta2)*sin_z2

        sin_a2=SIN((theta1-theta2)/R_TWO)**2+ss
        cos_a2=COS((theta1-theta2)/R_TWO)**2-ss
    ENDSUBROUTINE

    SUBROUTINE GET_ROTATION_ANGLE_SCALED(theta_p, phi_p, theta, phi, sin_alpha, cos_alpha )
        REAL(DP), INTENT(IN) :: theta_p, phi_p, theta, phi
        REAL(DP), INTENT(OUT) :: cos_alpha, sin_alpha
        REAL(DP) :: z
        z=phi-phi_p
        cos_alpha=SIN(theta-theta_p)+R_TWO*COS(theta)*SIN(theta_p)*SIN(z/R_TWO)**2
        sin_alpha=SIN(theta_p)*SIN(z)
    ENDSUBROUTINE

    SUBROUTINE SCALARS_FOR_TIPPERS_S2A(run_data, I_src,  mpi, SCALARS)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_src
        TYPE(MPI_STATE_TYPE), INTENT(IN):: mpi
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(OUT) :: SCALARS

        TYPE(SPH_DOMAIN_SHAPE_TYPE) :: TMP_domain
        TYPE (OPERATOR_D2D) :: S2A_OP
        REAL(DP) :: s
        INTEGER, PARAMETER :: Ng=8
        !CALL MERGE_TWO_DOMAINS(run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain,&
        !   run_data%SPH_Anomaly, &
        !  run_data%work_bkg, TMP_domain)
        CALL MERGE_TWO_DOMAINS(run_data%SPH_Anomaly, &
            run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain,&
            run_data%work_bkg, TMP_domain)

        ASSOCIATE(IE_OP=>S2A_OP%IE_OP)

            S2A_OP%Src=run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain
            S2A_OP%Recv=run_data%SPH_Anomaly

            IE_OP%mpi=mpi
            IE_OP%Domain=TMP_domain

            IE_OP%Bkg=run_data%work_bkg

            CALL ADD_TO_BKG(IE_OP%Domain, IE_OP%Bkg)
            CALL SET_BKG_INDEXES(IE_OP%BKG, IE_OP%Domain)

            IE_OP%Nterms=run_data%cpt%Nterms_RC

            S2A_OP%IE_OP%TENSOR_TYPE=COLLOC_KERNEL

            CALL SET_BKG_INDEXES(S2A_OP%IE_OP%bkg, S2A_OP%Src)
            CALL SET_BKG_INDEXES(S2A_OP%IE_OP%bkg, S2A_OP%Recv)

            S2A_OP%IE_OP%TN_rr=9*S2A_OP%Src%Nr*S2A_OP%Recv%Nr

            s=GET_MIN_S(S2A_OP%IE_OP)    

            CALL CREATE_SYMM_TABLE_STRUCTURE(s, run_data%cpt%STEP_POW, mpi, SCALARS%TABLE)
        END ASSOCIATE
        CALL CALCULATE_SCALARS_FOR_TIPPERS_S2A(run_data, I_src, S2A_OP, SCALARS)


    ENDSUBROUTINE

    SUBROUTINE CALCULATE_SCALARS_FOR_TIPPERS_S2A(run_data, I_source, S2A_OP, SCALARS)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_source
        TYPE (OPERATOR_D2D), INTENT(INOUT) :: S2A_OP
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(INOUT) :: SCALARS

        INTEGER:: I 
        REAL(DP), ALLOCATABLE :: C_MULTS(:)
        TYPE(TENSOR_COMPONENT_SERIES) :: TIPPERS_SERIES(2)


        CALL CALCULATE_SERIES_S2A(run_data, I_source, S2A_OP)  

        ASSOCIATE(IE_OP=>S2A_OP%IE_OP, SERIES=>S2A_OP%IE_OP%SERIES, &
                mpi=>S2A_OP%IE_OP%mpi, GR_TABLE=>S2A_OP%IE_OP%GR_TABLE)


            CALL ALLOCATE_SERIES2(SERIES(1)%Nterms, SERIES(1)%Nrr,&
                SERIES(1)%N_loc, TIPPERS_SERIES(1))

            CALL ALLOCATE_SERIES2(SERIES(1)%Nterms, SERIES(1)%Nrr,&
                SERIES(1)%N_loc, TIPPERS_SERIES(2))

            ALLOCATE(C_MULTS(0:IE_OP%Nterms))

            CALL TIPPERS_SERIES_CORRECTORS(C_MULTS)
            DO I=SERIES(1)%N_loc(1), SERIES(1)%N_loc(2)
                TIPPERS_SERIES(1)%S(I,:)=SERIES(RC_E_LP)%S(I,:)*C_MULTS(I)
                TIPPERS_SERIES(2)%S(I,:)=SERIES(RC_E_RL)%S(I,:)*C_MULTS(I)
            ENDDO

            DEALLOCATE(C_MULTS)

            CALL SERIES_TO_TABLE_FUNC_BY_PARTS(TIPPERS_SERIES, SCALARS, mpi)

            CALL CHECK_MEM("SCALARS FOR SOURCE FOR TIPPERS HAS BEEN  OBTAINED", S2A_OP%IE_OP%MPI)

            DEALLOCATE(TIPPERS_SERIES(2)%A_CORR, TIPPERS_SERIES(2)%S)
            DEALLOCATE(TIPPERS_SERIES(1)%A_CORR, TIPPERS_SERIES(1)%S)

            CALL DEALLOCATE_SERIES(S2A_OP)
            IF (run_data%Sources(I_source)%dump_scalars) CALL DUMP_S2A_SCALARS(SCALARS, mpi)

#ifdef dump_scalar
            CALL FINALIZE_GEMMIE()
            STOP "DUMP SCALARS"
#endif

        END ASSOCIATE
    ENDSUBROUTINE

    SUBROUTINE CALCULATE_SCALARS_FOR_TIPPERS_S2O(run_data, I_source, S2O_OP, SCALARS)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_source
        TYPE (OPERATOR_D2O), INTENT(INOUT) :: S2O_OP
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(INOUT) :: SCALARS

        INTEGER:: I 
        REAL(DP), ALLOCATABLE :: C_MULTS(:)
        TYPE(TENSOR_COMPONENT_SERIES) :: TIPPERS_SERIES(3)


        CALL CALCULATE_SERIES_S2O(run_data, I_source, S2O_OP)  

        ASSOCIATE(IE_OP=>S2O_OP%IE_OP, SERIES=>S2O_OP%IE_OP%SERIES, &
                mpi=>S2O_OP%IE_OP%mpi, GR_TABLE=>S2O_OP%IE_OP%GR_TABLE)

            DO I=1,3
                CALL ALLOCATE_SERIES2(SERIES(1)%Nterms, SERIES(1)%Nrr,&
                    SERIES(1)%N_loc, TIPPERS_SERIES(I))
            ENDDO
            ALLOCATE(C_MULTS(0:IE_OP%Nterms))
        
            CALL TIPPERS_SERIES_CORRECTORS(C_MULTS)

            DO I=SERIES(1)%N_loc(1), SERIES(1)%N_loc(2)
                TIPPERS_SERIES(1)%S(I,:)=SERIES(RC_E_LP)%S(I,:)*C_MULTS(I)
                TIPPERS_SERIES(2)%S(I,:)=SERIES(RC_E_RL)%S(I,:)*C_MULTS(I)
                TIPPERS_SERIES(3)%S(I,:)=SERIES(RC_H_LP)%S(I,:)*C_MULTS(I)&
                    *IE_OP%bkg%MainSphereRadius
            ENDDO

            DEALLOCATE(C_MULTS)

            CALL SERIES_TO_TABLE_FUNC_BY_PARTS(TIPPERS_SERIES, SCALARS, mpi)

            DO I=3,1,-1
                DEALLOCATE(TIPPERS_SERIES(I)%A_CORR, TIPPERS_SERIES(I)%S )
            ENDDO

            CALL DEALLOCATE_SERIES(S2O_OP)

            CALL CHECK_MEM("SCALARS FOR SOURCE FOR TIPPERS HAS BEEN  OBTAINED", S2O_OP%IE_OP%MPI)

            IF (run_data%Sources(I_source)%dump_scalars) CALL DUMP_S2O_SCALARS(GR_TABLE, mpi)

#ifdef dump_scalar
            CALL FINALIZE_GEMMIE()
            STOP "DUMP SCALARS"
#endif

        END ASSOCIATE
    ENDSUBROUTINE


    SUBROUTINE SCALARS_FOR_TIPPERS_S2O(run_data, I_src, I_obs_r,  mpi, SCALARS)
        TYPE(RUN_INPUT), INTENT(IN) :: run_data
        INTEGER, INTENT(IN) :: I_src, I_obs_r
        TYPE(MPI_STATE_TYPE), INTENT(IN):: mpi
        TYPE (GREEN_SCALARS_TABLE_TYPE), INTENT(OUT) :: SCALARS

        TYPE(SPH_DOMAIN_SHAPE_TYPE) :: TMP_domain
        TYPE (OPERATOR_D2O) :: S2O_OP
        REAL(DP) :: s
        REAL(DP) :: anom_bnd(2)
        INTEGER :: Mr
        S2O_OP%r_obs=run_data%sph_recvs%r(I_obs_r)

        CALL ADD_LAYER_TO_DOMAIN(S2O_OP%r_obs,&
            run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain,&
            run_data%work_bkg, tmp_domain)

        anom_bnd(1)=run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain%r(0)
        Mr=run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain%Nr
        anom_bnd(2)=run_data%Sources(I_src)%Currents_Data%SPH_Source_Domain%r(Mr)

        ASSOCIATE(IE_OP=>S2O_OP%IE_OP)

            IE_OP%mpi=mpi
            IE_OP%Domain=tmp_domain
            IE_OP%Bkg=run_data%work_bkg

            CALL ADD_TO_BKG(IE_OP%Domain, IE_OP%Bkg)
            CALL SET_BKG_INDEXES(IE_OP%BKG, IE_OP%Domain)
            IE_OP%Nterms=run_data%cpt%Nterms_RC


            S2O_OP%IE_OP%TENSOR_TYPE=COLLOC_KERNEL

            s=1e-2_DP    

            CALL CREATE_SYMM_TABLE_STRUCTURE( s, run_data%cpt%STEP_POW, mpi, SCALARS%TABLE)

        END ASSOCIATE

        S2O_OP%L_obs=GET_OBS_LAYER_INDEX(S2O_OP%r_obs, S2O_OP%IE_OP%bkg%r)
        S2O_OP%L_src=GET_BACKGROUND_LAYER_INDEX(anom_bnd, S2O_OP%IE_OP%bkg%r)
        IF (S2O_OP%r_obs<=S2O_OP%IE_OP%Domain%r(1)) THEN
            S2O_OP%L_obs_r_op=1
        ELSE
            S2O_OP%L_obs_r_op=GET_OBS_LAYER_INDEX(S2O_OP%r_obs, S2O_OP%IE_OP%Domain%r)
        ENDIF

        S2O_OP%Nr_src=S2O_OP%L_src(2)-S2O_OP%L_src(1)+1

        S2O_OP%IE_OP%TENSOR_TYPE=RC_EH_KERNEL
        S2O_OP%IE_OP%TN_rr=17*S2O_OP%Nr_src

        CALL CALCULATE_SCALARS_FOR_TIPPERS_S2O(run_data, I_src, S2O_OP, SCALARS)

    ENDSUBROUTINE


    SUBROUTINE TIPPERS_SERIES_CORRECTORS(A)
        REAL(DP), INTENT(INOUT) :: A(0:)
        REAL(DP), ALLOCATABLE :: B(:)
        INTEGER :: I, N
        N=UBOUND(A,1)
        ALLOCATE(B(0:N))
        A(0)=R_ZERO
        A(1)=PI**2
        B(0)=-R_TWO*PI**2
        B(1)=R_ZERO
        DO I=2, N
            A(I)=-( (R_TWO*I-R_ONE)*B(I-1)+ (I-R_TWO)*A(I-2))/(I+R_ONE)
            B(I)=-( (R_TWO*I-R_ONE)*A(I-1)+(I-R_ONE)*B(I-2))/I
        ENDDO
    ENDSUBROUTINE

END MODULE  GLOBE_SOURCES
