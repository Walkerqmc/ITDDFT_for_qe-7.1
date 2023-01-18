!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included gate related energy
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !! General self-consistency loop, also for hybrid functionals.  
  !! For non-hybrid functionals it just calls \(\texttt{electron_scf}\).
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE gvect,                ONLY : ecutrho
  USE lsda_mod,             ONLY : nspin, absmag
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, vsol, esol, ef_up, ef_dw
  USE tsvdw_module,         ONLY : EtsvdW
  USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : tr2, niter, conv_elec, restart, lmd, &
                                   do_makov_payne
  USE io_files,             ONLY : iunres, seqopn
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : tefield, etotefield
  USE wvfct,                ONLY : nbnd, wg, et
  USE klist,                ONLY : nks
  USE uspp,                 ONLY : okvan
  USE exx,                  ONLY : aceinit,exxinit, exxenergy2, exxbuff, &
                                   fock0, fock1, fock2, fock3, dexx, use_ace, local_thr 
  USE xc_lib,               ONLY : xclib_dft_is, exx_is_active
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi, gamma_only
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE ions_base,            ONLY : nat
  USE loc_scdm,             ONLY : use_scdm, localize_orbitals
  USE loc_scdm_k,           ONLY : localize_orbitals_k
  !
  USE wvfct_gpum,           ONLY : using_et, using_wg, using_wg_d
  USE scf_gpum,             ONLY : using_vrs
  !
  USE rism_module,          ONLY : lrism, rism_calc3d
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: exxen
  !! used to compute exchange energy
  REAL(DP), EXTERNAL :: exxenergyace
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER :: iter
  !! counter on iterations
  INTEGER :: printout, ik, ios
  !
  REAL(DP) :: tr2_min
  !! estimated error on energy coming from diagonalization
  REAL(DP) :: tr2_final
  !! final threshold for exx minimization 
  !! when using adaptive thresholds.
  LOGICAL :: first, exst
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  LOGICAL :: DoLoc
  !
  !
  DoLoc = local_thr.gt.0.0d0
  exxen = 0.0d0
  iter = 0
  first = .true.
  tr2_final = tr2
  IF ( xclib_dft_is('hybrid') ) THEN
     !printout = 0  ! do not print etot and energy components at each scf step
     printout = 1  ! print etot, not energy components at each scf step
  ELSE IF ( lmd ) THEN
     printout = 1  ! print etot, not energy components at each scf step
  ELSE
     printout = 2  ! print etot and energy components at each scf step
  ENDIF
  IF (xclib_dft_is('hybrid') .AND. adapt_thr ) tr2= tr2_init
  fock0 = 0.D0
  fock1 = 0.D0
  fock3 = 0.D0
  IF (.NOT. exx_is_active () ) fock2 = 0.D0
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%  Iterate hybrid functional  %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  IF ( restart ) THEN
     CALL seqopn (iunres, 'restart_e', 'formatted', exst)
     IF ( exst ) THEN
        ios = 0
        READ (iunres, *, iostat=ios) iter, tr2, dexx
        IF ( ios /= 0 ) THEN
           iter = 0
        ELSE IF ( iter < 0 .OR. iter > niter ) THEN
           iter = 0
        ELSE 
           READ (iunres, *) exxen, fock0, fock1, fock2
           ! FIXME: et and wg should be read from xml file
           READ (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           READ (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           CLOSE ( unit=iunres, status='delete')
           CALL using_et(2); CALL using_wg(2)
#if defined(__CUDA)
           CALL using_wg_d(0)
#endif
           ! ... if restarting here, exx was already active
           ! ... initialize stuff for exx
           first = .false.
           CALL exxinit(DoLoc)
           IF( DoLoc.and.gamma_only) THEN
             CALL localize_orbitals( )
           ELSE IF (DoLoc) THEN
             CALL localize_orbitals_k( )
           ENDIF 
           ! FIXME: ugly hack, overwrites exxbuffer from exxinit
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           IF (exst) READ (iunres, iostat=ios) exxbuff
           IF (ios /= 0) WRITE(stdout,'(5x,"Error in EXX restart!")')
           IF (use_ace) CALL aceinit ( DoLoc )
           !
           CALL v_of_rho( rho, rho_core, rhog_core, &
               ehart, etxc, vtxc, eth, etotefield, charge, v)
           IF (lrism) CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
           IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw,etot_cmp_paw)
           CALL using_vrs(1)
           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
                         nspin, doublegrid )
           !
           WRITE(stdout,'(5x,"Calculation (EXX) restarted from iteration #", &
                        & i6)') iter 
        ENDIF
     ENDIF
     CLOSE ( unit=iunres, status='delete')
  ENDIF
  !
  DO idum=1,niter
     !
     iter = iter + 1
     !
     ! ... Self-consistency loop. For hybrid functionals the exchange potential
     ! ... is calculated with the orbitals at previous step (none at first step)
     !
     CALL electrons_scf ( printout, exxen )
     !
     IF ( .NOT. xclib_dft_is('hybrid') ) RETURN
     !
     ! ... From now on: hybrid DFT only
     !
     IF ( stopped_by_user .OR. .NOT. conv_elec ) THEN
        conv_elec=.FALSE.
        IF ( .NOT. first) THEN
           CALL using_et(0)
           WRITE(stdout,'(5x,"Calculation (EXX) stopped during iteration #", &
                        & i6)') iter
           CALL seqopn (iunres, 'restart_e', 'formatted', exst)
           WRITE (iunres, *) iter-1, tr2, dexx
           WRITE (iunres, *) exxen, fock0, fock1, fock2
           WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
           WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
           CLOSE (unit=iunres, status='keep')
           CALL seqopn (iunres, 'restart_exx', 'unformatted', exst)
           WRITE (iunres) exxbuff
           CLOSE (unit=iunres, status='keep')
        ENDIF
        RETURN
     ENDIF
     !
     first =  first .AND. .NOT. exx_is_active ( )
     !
     ! "first" is true if the scf step was performed without exact exchange
     !
     IF ( first ) THEN
        !
        first = .false.
        !
        ! Activate exact exchange, set orbitals used in its calculation,
        ! then calculate exchange energy (will be useful at next step)
        !
        CALL exxinit(DoLoc)
        IF( DoLoc.and.gamma_only) THEN
          CALL localize_orbitals( )
        ELSE IF (DoLoc) THEN
          CALL localize_orbitals_k( )
        ENDIF 
        IF (use_ace) THEN
           CALL aceinit ( DoLoc ) 
           fock2 = exxenergyace()
        ELSE
           fock2 = exxenergy2()
        ENDIF
        exxen = 0.50d0*fock2 
        etot = etot - etxc 
        !
        ! Recalculate potential because XC functional has changed,
        ! start self-consistency loop on exchange
        !
        CALL v_of_rho( rho, rho_core, rhog_core, &
             ehart, etxc, vtxc, eth, etotefield, charge, v)
        etot = etot + etxc + exxen
        !
        IF (lrism) CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
        !
        IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw,etot_cmp_paw)
        CALL using_vrs(1)
        CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
             nspin, doublegrid )
        !
     ELSE
        !
        ! fock1 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n-1 in the expression of exchange
        !
        IF (use_ace) THEN
           fock1 = exxenergyace()
        ELSE
           fock1 = exxenergy2()
        ENDIF
        !
        ! Set new orbitals for the calculation of the exchange term
        !
        CALL exxinit(DoLoc)
        IF( DoLoc.and.gamma_only) THEN
          CALL localize_orbitals( )
        ELSE IF (DoLoc) THEN
          CALL localize_orbitals_k( )
        ENDIF 
        IF (use_ace) CALL aceinit ( DoLoc, fock3 )
        !
        ! fock2 is the exchange energy calculated for orbitals at step n,
        !       using orbitals at step n in the expression of exchange 
        ! fock0 is fock2 at previous step
        !
        fock0 = fock2
        IF (use_ace) THEN
           fock2 = exxenergyace()
        ELSE
           fock2 = exxenergy2()
        ENDIF
        !
        ! check for convergence. dexx is positive definite: if it isn't,
        ! there is some numerical problem. One such cause could be that
        ! the treatment of the divergence in exact exchange has failed. 
        ! FIXME: to be properly implemented for all cases
        !
!civn 
!       IF (use_ace .AND. (nspin == 1) .AND. gamma_only) THEN
        IF ( DoLoc ) THEN
          dexx =  0.5D0 * ((fock1-fock0)+(fock3-fock2)) 
        ELSE
          dexx = fock1 - 0.5D0*(fock0+fock2)
        ENDIF 
        !
        IF ( dexx < 0.0_dp ) THEN
           !IF( Doloc ) THEN
              WRITE(stdout,'(5x,a,1e12.3)') "BEWARE: negative dexx:", dexx
              dexx = ABS ( dexx )
           !ELSE
           !   CALL errore( 'electrons', 'dexx is negative! &
           !        & Check that exxdiv_treatment is appropriate for the system,&
           !        & or ecutfock may be too low', 1 )
           !ENDIF
        ENDIF
        !
        !   remove the estimate exchange energy exxen used in the inner SCF
        !
        etot = etot + exxen + 0.5D0*fock2 - fock1
        hwf_energy = hwf_energy + exxen + 0.5D0*fock2 - fock1 ! [LP]
        exxen = 0.5D0*fock2 
        !
        IF ( dexx < tr2_final ) THEN
           WRITE( stdout, 9066 ) '!!', etot, hwf_energy
        ELSE
           WRITE( stdout, 9066 ) '  ', etot, hwf_energy
        ENDIF
        IF ( dexx>1.d-8 ) THEN
          WRITE( stdout, 9067 ) dexx
        ELSE
          WRITE( stdout, 9068 ) dexx
        ENDIF
        
        WRITE( stdout, 9062 ) - fock1
        IF (use_ace) THEN
           WRITE( stdout, 9063 ) 0.5D0*fock2
        ELSE
           WRITE( stdout, 9064 ) 0.5D0*fock2
        ENDIF
        !
        IF ( dexx < tr2_final ) THEN
           IF ( do_makov_payne ) CALL makov_payne( etot )
           WRITE( stdout, 9101 )
           RETURN
        ENDIF
        !
        IF ( adapt_thr ) THEN
           tr2 = MAX(tr2_multi * dexx, tr2_final)
           WRITE( stdout, 9121 ) tr2
        ENDIF
     ENDIF
     !
     WRITE( stdout,'(/5x,"EXX: now go back to refine exchange calculation")')
     !
     IF ( check_stop_now() ) THEN
        CALL using_et(0)
        WRITE(stdout,'(5x,"Calculation (EXX) stopped after iteration #", &
                        & i6)') iter
        conv_elec=.FALSE.
        CALL seqopn (iunres, 'restart_e', 'formatted', exst)
        WRITE (iunres, *) iter, tr2, dexx
        WRITE (iunres, *) exxen, fock0, fock1, fock2
        ! FIXME: et and wg are written to xml file
        WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
        WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
        CLOSE (unit=iunres, status='keep')
        RETURN
     ENDIF
     !
  ENDDO
  !
  WRITE( stdout, 9120 ) iter
  FLUSH( stdout )
  !
  RETURN
  !
  ! ... formats
  !
9062 FORMAT( '     - averaged Fock potential =',0PF17.8,' Ry' )
9063 FORMAT( '     + Fock energy (ACE)       =',0PF17.8,' Ry' )
9064 FORMAT( '     + Fock energy (full)      =',0PF17.8,' Ry' )
9066 FORMAT(/,A2,'   total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9067 FORMAT('     est. exchange err (dexx)  =',0PF17.8,' Ry' )
9068 FORMAT('     est. exchange err (dexx)  =',1PE17.1,' Ry' )
9101 FORMAT(/'     EXX self-consistency reached' )
9120 FORMAT(/'     EXX convergence NOT achieved after ',i3,' iterations: stopping' )
9121 FORMAT(/'     scf convergence threshold =',1PE17.1,' Ry' )
  !
END SUBROUTINE electrons
!
!----------------------------------------------------------------------------
SUBROUTINE electrons_scf ( printout, exxen )
  !----------------------------------------------------------------------------
  !! This routine is a driver of the self-consistent cycle.  
  !! It uses the routine \(\texttt{c_bands}\) for computing the bands at fixed
  !! Hamiltonian, the routine \(\texttt{sum_band}\) to compute the charge density,
  !! the routine \(\texttt{v_of_rho}\) to compute the new potential and the routine
  !! \(\text{mix_rho}\) to mix input and output charge densities.
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor, atm, &
                                   ntyp => nsp
  USE basis,                ONLY : starting_pot
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                   two_fermi_energies, tot_charge
  USE fixed_occ,            ONLY : one_atom_occupations
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, ef_up, ef_dw, exdm, ef, &
                                   egrand, vsol, esol
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm, ts_vdw, &
                                   mbd_vdw, use_gpu
  USE control_flags,        ONLY : n_scf_steps, scf_error

  USE io_files,             ONLY : iunmix, output_drho
  USE ldaU,                 ONLY : eth, lda_plus_u, lda_plus_u_kind, &
                                   niter_with_fixed_ns, hub_pot_fix, &
                                   nsg, nsgnew, v_nsg, at_sc, neighood, &
                                   ldim_u, is_hubbard_back
  USE extfield,             ONLY : tefield, etotefield, gate, etotgatefield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report, domag
  USE io_rho_xml,           ONLY : write_scf
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : root_pool, me_pool, my_pool_id, &
                                   inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, get_atomic_number
  USE dftd3_qe,             ONLY : dftd3
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE dfunct,               ONLY : newd
  USE dfunct_gpum,          ONLY : newd_gpu
  USE esm,                  ONLY : do_comp_esm, esm_printpot, esm_ewald
  USE gcscf_module,         ONLY : lgcscf, gcscf_mu, gcscf_ignore_mun, gcscf_set_nelec
  USE clib_wrappers,        ONLY : memstat
  USE fcp_module,           ONLY : lfcp, fcp_mu
  USE rism_module,          ONLY : lrism, rism_calc3d, rism_printpot
  USE iso_c_binding,        ONLY : c_int
  !
  USE plugin_variables,     ONLY : plugin_etot
  USE libmbd_interface,     ONLY : EmbdvdW
  USE add_dmft_occ,         ONLY : dmft, dmft_update, v_dmft, dmft_updated
  !
  USE wvfct_gpum,           ONLY : using_et
  USE scf_gpum,             ONLY : using_vrs
  USE device_fbuff_m,       ONLY : dev_buf, pin_buf
  USE pwcom,                ONLY : report_mag 
  !
#if defined (__ENVIRON)
  USE plugin_flags,         ONLY : use_environ
  USE environ_base_module,  ONLY : calc_environ_energy, print_environ_energies
  USE environ_pw_module,    ONLY : calc_environ_potential
#endif
  !


!! begin modification

  USE ldaU,                 ONLY : eth, lda_plus_u, lda_plus_u_kind, &
                                   niter_with_fixed_ns, hub_pot_fix, &
                                   nsg, nsgnew, v_nsg, at_sc, neighood, &
                                   ldim_u, is_hubbard_back, &
                                   run_ITDDFT, g_dependant_dtau, dtau, ITDDFT_backup



  USE wvfct,                ONLY : wg, et, npwx
  USE control_flags,        ONLY : io_level, gamma_only    
  USE noncollin_module,     ONLY : npol
  USE uspp,                 ONLY : okvan, vkb, nkb, qq_nt
  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE mp_world,             ONLY : root, world_comm


!! end modification



  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: printout
  !! * If printout>0, prints on output the total energy;
  !! * if printout>1, also prints decomposition into energy contributions.
  REAL(DP),INTENT (IN) :: exxen
  !! current estimate of the exchange energy
  !
  ! ... local variables
  !
  REAL(DP) :: dr2
  !! the norm of the diffence between potential
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: deband_hwf
  !! deband for the Harris-Weinert-Foulkes functional
  REAL(DP) :: mag
  !! local magnetization
  INTEGER :: i
  !! counter on polarization
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER :: iter
  !! counter on iterations
  INTEGER :: nt
  !! counter on atomic types
  INTEGER :: ios, kilobytes
  !
  REAL(DP) :: tr2_min
  !! estimated error on energy coming from diagonalization
  REAL(DP) :: descf
  !! correction for variational energy
  REAL(DP) :: en_el=0.0_DP
  !! electric field contribution to the total energy
  REAL(DP) :: eext=0.0_DP
  !! external forces contribution to the total energy
  LOGICAL :: first, exst
  !! auxiliary variables for calculating and storing temporary copies of
  !! the charge density and of the HXC-potential
  !
  TYPE(scf_type) :: rhoin
  !! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  ! 
  REAL(DP) :: latvecs(3,3)
  !! auxiliary variables for grimme-d3
  INTEGER:: atnum(1:nat), na
  !! auxiliary variables for grimme-d3
  LOGICAL :: lhb
  !! if .TRUE. then background states are present (DFT+U)
  !

!! begin modification



  REAL(DP), ALLOCATABLE :: energy_record(:)
  REAL(DP) :: original_tr2, original_mixing_beta
  INTEGER ::  energy_increases, energy_recordings, original_nmix
  LOGICAL :: switched_to_ITDDFT, original_run_ITDDFT, recalc


  CALL mp_bcast(run_ITDDFT, root, world_comm)
  CALL mp_bcast(g_dependant_dtau, root, world_comm)
  CALL mp_bcast(dtau, root, world_comm)
  CALL mp_bcast(ITDDFT_backup, root, world_comm)

   recalc = .TRUE.
   switched_to_ITDDFT = .FALSE.

   ALLOCATE(energy_record(24))
   energy_record = 1e10
   energy_recordings = 0


!! end modification


  lhb = .FALSE.
  IF ( lda_plus_u )  THEN
     DO nt = 1, ntyp
        IF (is_hubbard_back(nt)) lhb = .TRUE.
     ENDDO
  ENDIF
  !
  iter = 0
  dr2  = 0.0_dp
  IF ( restart ) CALL restart_in_electrons( iter, dr2, ethr, et )
  IF ( restart ) CALL using_et(2)
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  CALL memstat( kilobytes )
  IF ( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  CALL start_clock( 'electrons' )
  !
  FLUSH( stdout )
  !
  ! ... calculates the ewald contribution to total energy
  !
  IF ( do_comp_esm ) THEN
     ewld = esm_ewald()
  ELSE
     ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  ENDIF
  !
  IF ( llondon ) THEN
     elondon = energy_london( alat , nat , ityp , at ,bg , tau )
  ELSE
     elondon = 0.d0
  ENDIF
  !
  ! Grimme-D3 correction to the energy
  !
  IF (ldftd3) THEN
     CALL start_clock('energy_dftd3')
     latvecs(:,:)=at(:,:)*alat
     tau(:,:)=tau(:,:)*alat
     DO na = 1, nat
        atnum(na) = get_atomic_number(TRIM(atm(ityp(na))))
     ENDDO
     call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,edftd3)
     edftd3=edftd3*2.d0
     tau(:,:)=tau(:,:)/alat
     CALL stop_clock('energy_dftd3')
  ELSE
     edftd3= 0.0
  ENDIF
  !
  !
  CALL create_scf_type( rhoin )
  !
  WRITE( stdout, 9002 )
  FLUSH( stdout )
  !
  CALL open_mix_file( iunmix, 'mix', exst )
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter
     !
     IF ( check_stop_now() ) THEN
        conv_elec=.FALSE.
        CALL using_et(0)
        CALL save_in_electrons (iter, dr2, ethr, et )
        GO TO 10
     ENDIF
     iter = iter + 1
     !
     IF (use_gpu .and. (iverbosity >= 1) ) CALL dev_buf%print_report(stdout)
     IF (use_gpu .and. (iverbosity >= 1) ) CALL pin_buf%print_report(stdout)
     !
     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
     !
     FLUSH( stdout )
     !
     ! ... Convergence threshold for iterative diagonalization is
     ! ... automatically updated during self consistency
     !
     IF ( iter > 1 ) THEN
        !
        IF ( iter == 2 ) THEN
           !
           IF ( lgcscf ) THEN
              ethr = 1.D-5
           ELSE
              ethr = 1.D-2
           ENDIF
           !
        ENDIF
        !
        ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
        ! ... do not allow convergence threshold to become too small:
        ! ... iterative diagonalization may become unstable
        ethr = MAX( ethr, 1.D-13 )
        !
     ENDIF
     !
     first = ( iter == 1 )
     !
     ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
     ! ... first time here using the input density and potential ( to be
     ! ... used to calculate the Harris-Weinert-Foulkes energy )
     !
     deband_hwf = delta_e()
     !
     ! ... save input density to rhoin
     !
     call scf_type_COPY( rho, rhoin )
     !
     scf_step: DO
        !
        ! ... tr2_min is set to an estimate of the error on the energy
        ! ... due to diagonalization - used only for the first scf iteration
        !
        tr2_min = 0.D0
        !
        IF ( first ) tr2_min = ethr*MAX( 1.D0, nelec ) 
        !
        ! ... diagonalization of the KS hamiltonian

        ! begin mod

        IF( ITDDFT_backup .AND. idum > 2 ) THEN

          energy_recordings = energy_recordings + 1

          DO i = MIN( SIZE(energy_record), energy_recordings ), 2, -1

            energy_record(i) = energy_record(i-1)

          END DO

          energy_record(1) = etot 


          IF( energy_recordings >= SIZE(energy_record) ) THEN

            energy_increases = 0
            DO i = 1, SIZE(energy_record)-1
              IF( energy_record(i) > energy_record(i+1) ) energy_increases = energy_increases + 1
            END DO        

            IF( energy_increases > 10 .AND. .NOT. run_ITDDFT ) THEN

              WRITE(stdout,*) " Switching to it-TDDFT"

              switched_to_ITDDFT  = .TRUE.

              original_run_ITDDFT = run_ITDDFT
              original_tr2 = tr2
              original_mixing_beta = mixing_beta
              original_nmix = nmix         

              run_ITDDFT = .TRUE.                   
              
              IF( g_dependant_dtau ) THEN          
                 tr2 = tr2/5.0d4
              ELSE
                 tr2 = 1.0d-16                
              END IF
              
              mixing_beta = 1
              nmix = 1

            END IF
         
          END IF

        END IF



        IF( run_ITDDFT ) THEN

          WRITE(stdout,*) "  running ITDDFT"

          CALL itddft()

        ELSE

          WRITE(stdout,*) "  running SCF"

          IF ( lelfield ) THEN
             CALL c_bands_efield ( iter )
          ELSE
             CALL c_bands( iter )
          END IF

        END IF


        !
        IF ( stopped_by_user ) THEN
           conv_elec=.FALSE.
           CALL using_et( 0 )
           CALL save_in_electrons( iter-1, dr2, ethr, et )
           GO TO 10
        ENDIF
        !
        IF (one_atom_occupations) CALL new_evc()
        !
        ! ... xk, wk, isk, et, wg are distributed across pools;
        ! ... the first node has a complete copy of xk, wk, isk,
        ! ... while eigenvalues et and weights wg must be
        ! ... explicitly collected to the first node
        ! ... this is done here for et, in sum_band for wg
        !
        CALL using_et(1)
        CALL poolrecover( et, nbnd, nkstot, nks )
        !
        ! ... the new density is computed here. For PAW:
        ! ... sum_band computes new becsum (stored in uspp modules)
        ! ... and a subtly different copy in rho%bec (scf module)
        !
        ! ... for charge self-consistent calculations call update
        ! ... of Fermi weights from DMFT here
        !
        IF ( dmft .AND. idum == 1) THEN
            CALL weights()
            CALL dmft_update()
        ELSEIF ( dmft .AND. idum > 1) THEN
            WRITE( stdout, '(5X,"WARNING: electron_maxstep > 1 not recommended for dmft = .true.")')
        END IF
        !
        IF (.not. use_gpu) CALL sum_band()
        IF (      use_gpu) CALL sum_band_gpu()
        !
        ! ... if DMFT update was made, make sure it is only done in the first iteration
        ! ... (generally in this mode it should only run a single iteration, but just to make sure!)
        !
        IF ( dmft .AND. idum == 1) THEN
            dmft_updated = .TRUE.
            DEALLOCATE(v_dmft)
        END IF
        !
        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        !
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
        If ( okpaw ) hwf_energy = hwf_energy + epaw
        IF ( lda_plus_u ) hwf_energy = hwf_energy + eth
        IF ( lrism ) hwf_energy = hwf_energy + esol + vsol
        !
        IF ( lda_plus_u )  THEN
           !
           ! Write the occupation matrices
           !
           IF ( iverbosity > 0 .OR. first ) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                 CALL write_ns()
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 IF (noncolin) THEN
                    CALL write_ns_nc()
                 ELSE
                    CALL write_ns()
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 CALL write_nsg()
              ENDIF
           ENDIF
           !
           ! Keep the Hubbard potential fixed, i.e. keep the 
           ! occupation matrix equal to the ground-state one.
           ! This is needed for the calculation of Hubbard parameters
           ! in a self-consistent way.
           !
           IF (hub_pot_fix) THEN
             IF (lda_plus_u_kind.EQ.0) THEN
                rho%ns = rhoin%ns ! back to input values
                IF (lhb) rho%nsb = rhoin%nsb
             ELSEIF (lda_plus_u_kind.EQ.1) THEN
                CALL errore('electrons_scf', &
                  & 'hub_pot_fix is not implemented for lda_plus_u_kind=1',1)
             ELSEIF (lda_plus_u_kind.EQ.2) THEN
                nsgnew = nsg
             ENDIF
           ENDIF
           !
           IF ( first .AND. starting_pot == 'atomic' ) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                 CALL ns_adj()
                 rhoin%ns = rho%ns
                 IF (lhb) rhoin%nsb = rho%nsb
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 CALL ns_adj()
                 IF (noncolin) THEN
                    rhoin%ns_nc = rho%ns_nc
                 ELSE
                    rhoin%ns = rho%ns
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 CALL nsg_adj()
              ENDIF
           ENDIF
           IF ( iter <= niter_with_fixed_ns ) THEN
              WRITE( stdout, '(/,5X,"RESET ns to initial values (iter <= mixing_fixed_ns)",/)')
              IF (lda_plus_u_kind.EQ.0) THEN
                 rho%ns = rhoin%ns
                 IF (lhb) rhoin%nsb = rho%nsb
              ELSEIF (lda_plus_u_kind.EQ.1) THEN
                 IF (noncolin) THEN
                    rho%ns_nc = rhoin%ns_nc
                 ELSE
                    rho%ns = rhoin%ns
                 ENDIF
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 nsgnew = nsg
              ENDIF
           ENDIF
           !
        ENDIF
        !
        ! ... calculate total and absolute magnetization
        !
        IF ( lsda .OR. noncolin ) CALL compute_magnetization()
        !
        ! ... eband  = \sum_v \epsilon_v    is calculated by sum_band
        ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
        ! ... eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
        !
        deband = delta_e()
        !
        ! ... mix_rho mixes several quantities: rho in g-space, tauk (for
        ! ... meta-gga), ns and ns_nc (for lda+u) and becsum (for paw)
        ! ... The mixing is done on pool 0 only (image parallelization
        ! ... inside mix_rho => rho_ddot => PAW_ddot is no longer there)
        !
        IF ( my_pool_id == root_pool ) CALL mix_rho( rho, rhoin, &
                mixing_beta, dr2, tr2_min, iter, nmix, iunmix, conv_elec )
        !
        ! ... Results are broadcast from pool 0 to others to prevent trouble
        ! ... on machines unable to yield the same results for the same 
        ! ... calculations on the same data, performed on different procs
        !
        IF ( lda_plus_u )  THEN
           ! ... For DFT+U, ns and ns_nc are also broadcast inside each pool
           ! ... to ensure consistency on all processors of all pools
           IF (noncolin) THEN
              IF (ALLOCATED(rhoin%ns_nc)) CALL mp_bcast( rhoin%ns_nc, root_pool, intra_pool_comm )
           ELSE
              IF (ALLOCATED(rhoin%ns)) CALL mp_bcast( rhoin%ns, root_pool, intra_pool_comm )
           ENDIF
           ! DFT+U+V: this variable is not in "mix-type" variable rhoin
           IF (lda_plus_u_kind.EQ.2) THEN
              IF (ALLOCATED(nsg) ) CALL mp_bcast ( nsg, root_pool, inter_pool_comm)
           ENDIF
        ENDIF
        !
        CALL bcast_scf_type( rhoin, root_pool, inter_pool_comm )
        CALL mp_bcast( dr2, root_pool, inter_pool_comm )
        CALL mp_bcast( conv_elec, root_pool, inter_pool_comm )
        !
        IF (.NOT. scf_must_converge .AND. idum == niter) conv_elec = .TRUE.
        !
        ! ... if convergence is achieved or if the self-consistency error
        ! ... (dr2) is smaller than the estimated error due to diagonalization
        ! ... (tr2_min), rhoin and rho are unchanged: rhoin contains the input
        ! ... density and rho contains the output density.
        ! ... In all other cases, rhoin contains the mixed charge density 
        ! ... (the new input density) while rho is unchanged
        !
        IF ( first .and. nat > 0) THEN
           !
           ! ... first scf iteration: check if the threshold on diagonalization
           ! ... (ethr) was small enough wrt the error in self-consistency (dr2)
           ! ... if not, perform a new diagonalization with reduced threshold
           !
           first = .FALSE.
           !
           IF ( dr2 < tr2_min ) THEN
              !
              WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                               &    "too large:",/,5X,                      &
                               & "Diagonalizing with lowered threshold",/)' )
              !
              ethr = 0.1D0*dr2 / MAX( 1.D0, nelec )
              !
              CYCLE scf_step
              !
           ENDIF
           !
        ENDIF
        !
        IF ( .NOT. conv_elec ) THEN
           !
           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate)
           !
           CALL v_of_rho( rhoin, rho_core, rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v )
           !
           IF (lrism) THEN
              CALL rism_calc3d(rhoin%of_g(:, 1), esol, vsol, v%of_r, dr2)
           ENDIF
           !
           IF (okpaw) THEN
              CALL PAW_potential( rhoin%bec, ddd_paw, epaw,etot_cmp_paw )
              CALL PAW_symmetrize_ddd( ddd_paw )
           ENDIF
           !
           ! ... estimate correction needed to have variational energy:
           ! ... T + E_ion (eband + deband) are calculated in sum_band
           ! ... and delta_e using the output charge density rho;
           ! ... E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
           ! ... above, using the mixed charge density rhoin%of_r.
           ! ... delta_escf corrects for this difference at first order
           !
           descf = delta_escf()
           !
           ! ... now copy the mixed charge density in R- and G-space in rho
           !
           CALL scf_type_COPY( rhoin, rho )
           !
           IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) nsgnew = nsg
           !
           IF ( lgcscf ) THEN
              CALL gcscf_set_nelec( charge )
           END IF
           !
        ELSE 
           !
           ! ... convergence reached:
           ! ... 1) the output HXC-potential is saved in v
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           !
           vnew%of_r(:,:) = v%of_r(:,:)
           !
           IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) nsg = nsgnew
           !
           CALL v_of_rho( rho,rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v )
           !
           IF (lrism) THEN
              CALL rism_calc3d(rho%of_g(:, 1), esol, vsol, v%of_r, tr2)
           END IF
           !
           vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
           !
           IF (okpaw) THEN
              CALL PAW_potential( rho%bec, ddd_paw, epaw, etot_cmp_paw )
              CALL PAW_symmetrize_ddd( ddd_paw )
           ENDIF
           !
           ! ... note that rho is here the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           !
           descf = 0._dp
           !
        ENDIF 
        !
        ! ... if we didn't cycle before we can exit the do-loop
        !
        EXIT scf_step
        !
     ENDDO scf_step
     !
     plugin_etot = 0.0_dp
     !
#if defined (__ENVIRON)
     IF (use_environ) THEN
        CALL calc_environ_energy(plugin_etot, .TRUE.)
        CALL calc_environ_potential(rhoin, conv_elec, dr2, vltot)
     END IF
#endif
     !
     ! ... define the total local potential (external + scf)
     !
     CALL using_vrs(1)
     CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
     !
     ! ... interpolate the total local potential
     !
     CALL using_vrs(1) ! redundant
     CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     ! ... PAW: newd contains PAW updates of NL coefficients
     !
     IF (.not. use_gpu) CALL newd()
     IF (      use_gpu) CALL newd_gpu()
     !
     IF ( lelfield ) en_el =  calc_pol ( )
     !
     IF ( report > 0 ) THEN
        IF ( conv_elec .OR.  MOD(iter,report) == 0 ) CALL report_mag()
     ELSE IF ( report < 0 ) THEN
        IF ( conv_elec ) CALL report_mag(SAVE_LOCALS=.TRUE.)
     END IF
     !
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     !
     IF ( conv_elec ) WRITE( stdout, 9101 )
 
     IF ( conv_elec ) THEN 
           scf_error = dr2
           n_scf_steps = iter
     ENDIF  

     !
     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 .OR. dmft_updated ) THEN
        !
        ! iverbosity == 0 for the PW code
        ! iverbosity >  2 for the HP code
        !
        IF ( lda_plus_u .AND. (iverbosity == 0 .OR. iverbosity > 2) ) THEN
           !
           ! Recompute the occupation matrix:
           ! needed when computing U (and V) in a SCF way
           IF (hub_pot_fix) THEN
              IF (lda_plus_u_kind.EQ.0) THEN
                 CALL new_ns(rho%ns)
                 IF (lhb) CALL new_nsb(rho%nsb)
              ELSEIF (lda_plus_u_kind.EQ.2) THEN
                 CALL new_nsg()
              ENDIF
           ENDIF
           !
           ! Write the occupation matrices
           IF (lda_plus_u_kind == 0) THEN
              CALL write_ns()
           ELSEIF (lda_plus_u_kind == 1) THEN
              IF (noncolin) THEN
                 CALL write_ns_nc()
              ELSE
                 CALL write_ns()
              ENDIF
           ELSEIF (lda_plus_u_kind == 2) THEN
              CALL write_nsg()
           ENDIF
           !
        ENDIF
        CALL print_ks_energies()
        !
     ENDIF
     !
     IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
        WRITE( stdout, 9050 ) charge, nelec
        IF ( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
           IF (.not.lgauss) THEN
              CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
           ELSE
              CALL errore( 'electrons', 'charge is wrong', 1 )
           ENDIF
        ENDIF
     ENDIF
     !
     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
     ! for hybrid calculations, add the current estimate of exchange energy
     ! (it will subtracted later if exx_is_active to be replaced with a better estimate)
     etot = etot - exxen
     hwf_energy = hwf_energy - exxen ! [LP]
     !
     IF (okpaw) etot = etot + epaw
     IF ( lda_plus_u ) etot = etot + eth
     !
     IF ( lelfield ) etot = etot + en_el
     ! not sure about the HWF functional in the above case
     IF( textfor ) THEN
        eext = alat*compute_eextfor()
        etot = etot + eext
        hwf_energy = hwf_energy + eext
     ENDIF
     IF (llondon) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
     ENDIF
     !
     ! grimme-d3 dispersion energy
     IF (ldftd3) THEN
        etot = etot + edftd3
        hwf_energy = hwf_energy + edftd3
     ENDIF
     !
     ! calculate the xdm energy contribution with converged density
     IF (lxdm .and. conv_elec) THEN
        exdm = energy_xdm()  
        etot = etot + exdm
        hwf_energy = hwf_energy + exdm
     ENDIF
     IF (mbd_vdw) THEN
        ! factor 2 converts from Hartree to Ry units
        etot = etot + 2.0d0*EmbdvdW
        hwf_energy = hwf_energy + 2.0d0*EmbdvdW
     ELSE IF (ts_vdw) THEN
        ! factor 2 converts from Hartree to Ry units
        etot = etot + 2.0d0*EtsvdW
        hwf_energy = hwf_energy + 2.0d0*EtsvdW
     ENDIF
     !
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     ENDIF
     ! TB gate energy
     IF ( gate ) THEN
        etot = etot + etotgatefield
        hwf_energy = hwf_energy + etotgatefield
     ENDIF
     !
     IF ( lgcscf .AND. (.NOT. gcscf_ignore_mun) ) THEN
        etot = etot + egrand
        hwf_energy = hwf_energy + egrand
     END IF
     !
     IF ( lfcp ) THEN
        etot = etot + fcp_mu * tot_charge
        hwf_energy = hwf_energy + fcp_mu * tot_charge
     END IF
     !
     IF ( lrism ) THEN
        etot = etot + esol + vsol
        hwf_energy = hwf_energy + esol + vsol
     END IF
     !
     ! ... adds possible external contribution from plugins to the energy
     !
     etot = etot + plugin_etot 
     !
     CALL print_energies ( printout )
     !
     IF ( conv_elec ) THEN
        !
        ! ... if system is charged add a Makov-Payne correction to the energy
        ! ... (not in case of hybrid functionals: it is added at the end)
        !
        IF ( do_makov_payne .AND. printout/= 0 ) CALL makov_payne( etot )
        !
        ! ... print out ESM potentials if desired
        !
        IF ( do_comp_esm ) CALL esm_printpot( rho%of_g )
        !
        ! ... print out 3D-RISM potentials if desired
        !
        IF ( lrism ) CALL rism_printpot()
        !
        WRITE( stdout, 9110 ) iter
        !
        ! ... jump to the end
        !
        GO TO 10
        !
     ENDIF
     !
     ! ... uncomment the following line if you wish to monitor the evolution
     ! ... of the force calculation during self-consistency
     !
     !CALL forces()
     !
     ! ... it can be very useful to track internal clocks during
     ! ... self-consistency for benchmarking purposes
#if defined(__PW_TRACK_ELECTRON_STEPS)
     CALL print_clock_pw()
#endif
     !
  ENDDO
  n_scf_steps = iter
  scf_error = dr2
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
10  FLUSH( stdout )
  !
  ! ... exiting: write (unless disabled) the charge density to file
  ! ... (also write ldaU ns coefficients and PAW becsum)
  !
  IF ( io_level > -1 ) CALL write_scf( rho, nspin )
  !
  ! ... delete mixing info if converged, keep it if not
  !
  IF ( conv_elec ) THEN
     CALL close_mix_file( iunmix, 'delete' )
  ELSE
     CALL close_mix_file( iunmix, 'keep' )
  ENDIF
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  call destroy_scf_type ( rhoin )
  CALL stop_clock( 'electrons' )
  !

!! begin modification

   IF( switched_to_ITDDFT ) THEN
      run_ITDDFT = original_run_ITDDFT
      tr2 = original_tr2
      mixing_beta = original_mixing_beta 
      nmix = original_nmix 
   END IF

  DEALLOCATE(energy_record)


!! end modification


  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2 )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_magnetization()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER :: ir
       !
       !
       IF ( lsda ) THEN
          !
          magtot = 0.D0
          absmag = 0.D0
          !
          DO ir = 1, dfftp%nnr
             !
             mag = rho%of_r(ir,2)
             !
             magtot = magtot + mag
             absmag = absmag + ABS( mag )
             !
          ENDDO
          !
          magtot = magtot * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
          CALL mp_sum( magtot, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
          !
          IF (two_fermi_energies.and.lgauss) bfield(3)=0.5D0*(ef_up-ef_dw)
          !
       ELSEIF ( noncolin ) THEN
          !
          magtot_nc = 0.D0
          absmag    = 0.D0
          !
          DO ir = 1,dfftp%nnr
             !
             mag = SQRT( rho%of_r(ir,2)**2 + &
                         rho%of_r(ir,3)**2 + &
                         rho%of_r(ir,4)**2 )
             !
             DO i = 1, 3
                !
                magtot_nc(i) = magtot_nc(i) + rho%of_r(ir,i+1)
                !
             ENDDO
             !
             absmag = absmag + ABS( mag )
             !
          ENDDO
          !
          CALL mp_sum( magtot_nc, intra_bgrp_comm )
          CALL mp_sum( absmag, intra_bgrp_comm )
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
             !
          ENDDO
          !
          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          !
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE compute_magnetization
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       ! This function computes delta_e, where:
       !
       ! ... delta_e =  - \int rho%of_r(r)  v%of_r(r)
       !                - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                - \sum rho%ns       v%ns       [for DFT+Hubbard]
       !                - \sum becsum       D1_Hxc     [for PAW]
       !
       USE xc_lib,  ONLY : xclib_dft_is
       !
       IMPLICIT NONE
       !
       REAL(DP) :: delta_e
       REAL(DP) :: delta_e_hub
       INTEGER  :: ir, na1, nt1, na2, nt2, m1, m2, equiv_na2, viz, is
       !
       delta_e = 0._DP
       IF ( nspin==2 ) THEN
          !
          DO ir = 1,dfftp%nnr
            delta_e = delta_e - ( rho%of_r(ir,1) + rho%of_r(ir,2) ) * v%of_r(ir,1) &  ! up
                              - ( rho%of_r(ir,1) - rho%of_r(ir,2) ) * v%of_r(ir,2)    ! dw
          ENDDO 
          delta_e = 0.5_DP*delta_e
          !
       ELSE
          delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
       ENDIF
       !
       IF ( xclib_dft_is('meta') ) &
          delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       IF (lda_plus_u .AND. (.NOT.hub_pot_fix)) THEN
         IF (lda_plus_u_kind.EQ.0) THEN
            delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
            IF (lhb) delta_e_hub = delta_e_hub - SUM(rho%nsb(:,:,:,:)*v%nsb(:,:,:,:))
            IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
            delta_e = delta_e + delta_e_hub
         ELSEIF (lda_plus_u_kind.EQ.1) THEN
            IF (noncolin) THEN
              delta_e_hub = - SUM( rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:) )
              delta_e = delta_e + delta_e_hub
            ELSE
              delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
              IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
              delta_e = delta_e + delta_e_hub
            ENDIF
         ELSEIF (lda_plus_u_kind.EQ.2) THEN
            delta_e_hub = 0.0_dp
            DO is = 1, nspin
               DO na1 = 1, nat
                  nt1 = ityp(na1)
                  DO viz = 1, neighood(na1)%num_neigh
                     na2 = neighood(na1)%neigh(viz)
                     equiv_na2 = at_sc(na2)%at
                     nt2 = ityp(equiv_na2)
                     IF ( ANY(v_nsg(:,:,viz,na1,is).NE.0.0d0) ) THEN
                        DO m1 = 1, ldim_u(nt1)
                           DO m2 = 1, ldim_u(nt2)
                              delta_e_hub = delta_e_hub - &
                                  nsgnew(m2,m1,viz,na1,is)*v_nsg(m2,m1,viz,na1,is)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
            delta_e = delta_e + delta_e_hub
         ENDIF
       ENDIF
       !
       IF (okpaw) delta_e = delta_e - SUM( ddd_paw(:,:,:)*rho%bec(:,:,:) )
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       ! This function calculates the difference between the Hartree and XC energy
       ! at first order in the charge density difference delta_rho(r).
       !
       ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
       !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                  - \sum \delta rho%ns       v%ns       [for DFT+Hubbard]
       !                  - \sum \delta becsum       D1         [for PAW] 
       !
       USE xc_lib, ONLY : xclib_dft_is
       !
       IMPLICIT NONE
       REAL(DP) :: delta_escf, delta_escf_hub, rho_dif(2)
       INTEGER  :: ir, na1, nt1, na2, nt2, m1, m2, equiv_na2, viz, is
       !
       delta_escf=0._dp
       IF ( nspin==2 ) THEN
          !
          DO ir=1, dfftp%nnr
             !
             rho_dif = rhoin%of_r(ir,:) - rho%of_r(ir,:)
             !
             delta_escf = delta_escf - ( rho_dif(1) + rho_dif(2) ) * v%of_r(ir,1) &  !up
                                     - ( rho_dif(1) - rho_dif(2) ) * v%of_r(ir,2)    !dw
          ENDDO
          delta_escf = 0.5_dp*delta_escf
          !
       ELSE
         delta_escf = -SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
       ENDIF
       !
       IF ( xclib_dft_is('meta') ) &
          delta_escf = delta_escf - &
                       SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
       !
       delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_escf, intra_bgrp_comm )
       !
       IF (lda_plus_u .AND. (.NOT.hub_pot_fix)) THEN
          IF (lda_plus_u_kind.EQ.0) THEN
             delta_escf_hub = -SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
             IF (lhb) delta_escf_hub = delta_escf_hub - &
                         SUM((rhoin%nsb(:,:,:,:)-rho%nsb(:,:,:,:))*v%nsb(:,:,:,:))
             IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
             delta_escf = delta_escf + delta_escf_hub
          ELSEIF (lda_plus_u_kind.EQ.1) THEN
             IF (noncolin) THEN
                delta_escf_hub = -SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
                delta_escf = delta_escf + delta_escf_hub
             ELSE
                delta_escf_hub = -SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
                IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
                delta_escf = delta_escf + delta_escf_hub
             ENDIF
          ELSEIF (lda_plus_u_kind.EQ.2) THEN
             delta_escf_hub = 0.0_dp
             DO is = 1, nspin
                DO na1 = 1, nat
                   nt1 = ityp(na1)
                   DO viz = 1, neighood(na1)%num_neigh
                      na2 = neighood(na1)%neigh(viz)
                      equiv_na2 = at_sc(na2)%at
                      nt2 = ityp(equiv_na2)
                      IF ( ANY(v_nsg(:,:,viz,na1,is).NE.0.0d0) ) THEN
                         DO m1 = 1, ldim_u(nt1)
                            DO m2 = 1, ldim_u(nt2)
                               delta_escf_hub = delta_escf_hub - &
                                    (nsg(m2,m1,viz,na1,is)-nsgnew(m2,m1,viz,na1,is)) * &
                                     v_nsg(m2,m1,viz,na1,is)
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             IF ( nspin==1 ) delta_escf_hub = 2.d0 * delta_escf_hub
             delta_escf = delta_escf + delta_escf_hub
          ENDIF
       ENDIF

       IF ( okpaw ) delta_escf = delta_escf - &
                                 SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

       RETURN
       !
     END FUNCTION delta_escf
     !
     !-----------------------------------------------------------------------
     FUNCTION calc_pol( ) RESULT ( en_el )
       !-----------------------------------------------------------------------
       !
       USE kinds,     ONLY : DP
       USE constants, ONLY : pi
       USE bp,        ONLY : lelfield, ion_pol, el_pol, fc_pol, l_el_pol_old, &
                             el_pol_acc, el_pol_old, efield, l3dstring, gdir, &
                             transform_el, efield_cart
       !
       IMPLICIT NONE
       REAL (DP) :: en_el
       !
       INTEGER :: i, j 
       REAL(DP):: sca, el_pol_cart(3),  el_pol_acc_cart(3)
       !
       IF (.NOT.l3dstring) THEN
          !
          CALL c_phase_field( el_pol(gdir), ion_pol(gdir), fc_pol(gdir), gdir )
          !
          IF (.NOT.l_el_pol_old) THEN
             l_el_pol_old = .TRUE.
             el_pol_old(gdir) = el_pol(gdir)
             en_el = -efield*(el_pol(gdir)+ion_pol(gdir))
             el_pol_acc(gdir) = 0.d0
          ELSE
             sca = (el_pol(gdir)-el_pol_old(gdir))/fc_pol(gdir)
             IF (sca < - pi) THEN
                el_pol_acc(gdir) = el_pol_acc(gdir)+2.d0*pi*fc_pol(gdir)
             ELSEIF (sca > pi) THEN
                el_pol_acc(gdir) = el_pol_acc(gdir)-2.d0*pi*fc_pol(gdir)
             ENDIF
             en_el = -efield*(el_pol(gdir)+ion_pol(gdir)+el_pol_acc(gdir))
             el_pol_old = el_pol
          ENDIF
          !
       ELSE
          !
          DO i = 1, 3
            CALL c_phase_field( el_pol(i), ion_pol(i), fc_pol(i), i )
          ENDDO
          el_pol_cart(:) = 0.d0
          DO i = 1, 3
             DO j = 1, 3
                !el_pol_cart(i)=el_pol_cart(i)+transform_el(j,i)*el_pol(j)
                el_pol_cart(i) = el_pol_cart(i)+at(i,j)*el_pol(j) / &
                                 (SQRT(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
             ENDDO
          ENDDO
          !
          WRITE( stdout,'( "Electronic Dipole on Cartesian axes" )' )
          DO i = 1, 3
             WRITE(stdout,*) i, el_pol_cart(i)
          ENDDO
          !
          WRITE( stdout,'( "Ionic Dipole on Cartesian axes" )' )
          DO i = 1, 3
             WRITE(stdout,*) i, ion_pol(i)
          ENDDO
          !
          IF (.NOT.l_el_pol_old) THEN
             l_el_pol_old = .TRUE.
             el_pol_old(:) = el_pol(:)
             en_el = 0.d0
             DO i = 1, 3
                en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i))
             ENDDO
             el_pol_acc(:) = 0.d0
          ELSE
             DO i = 1, 3
                sca = (el_pol(i)-el_pol_old(i))/fc_pol(i)
                IF (sca < - pi) THEN
                   el_pol_acc(i) = el_pol_acc(i)+2.d0*pi*fc_pol(i)
                ELSEIF (sca > pi) THEN
                   el_pol_acc(i) = el_pol_acc(i)-2.d0*pi*fc_pol(i)
                ENDIF
             ENDDO
             el_pol_acc_cart(:) = 0.d0
             DO i = 1, 3
                DO j = 1, 3
                   el_pol_acc_cart(i) = el_pol_acc_cart(i)+transform_el(j,i)*el_pol_acc(j)
                ENDDO
             ENDDO
             en_el = 0.d0
             DO i = 1, 3
                en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i)+el_pol_acc_cart(i))
             ENDDO
             el_pol_old(:) = el_pol(:)
          ENDIF
          !
       ENDIF
       !
     END FUNCTION calc_pol
     !
     !-----------------------------------------------------------------------
     SUBROUTINE print_energies ( printout )
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps8, RYTOEV
       INTEGER, INTENT (IN) :: printout
       !
   
       IF ( printout == 0 ) RETURN
       IF ( ( conv_elec .OR. MOD(iter,iprint) == 0 ) .AND. printout > 1 ) THEN
          !
          WRITE( stdout, 9081 ) etot
          IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          !
          IF ( lgcscf )  WRITE( stdout, 9181 ) tot_charge
          !
          IF ( lgauss ) then
             WRITE( stdout, 9070 ) demet
             WRITE( stdout, 9170 ) etot-demet
          END IF
          !
          IF ( lgauss ) then
             WRITE( stdout, 9061 )
          ELSE
             WRITE( stdout, 9060 )
          END IF
          !
          WRITE( stdout, 9062 ) (eband + deband), ehart, ( etxc - etxcc ), ewld
          !
          IF ( llondon ) WRITE ( stdout , 9074 ) elondon
          IF ( ldftd3 )  WRITE ( stdout , 9078 ) edftd3
          IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
          IF ( mbd_vdw ) THEN
             WRITE ( stdout , 9076 ) 2.0d0*Embdvdw
          ELSEIF ( ts_vdw ) THEN
             WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
          ENDIF
          IF ( textfor)  WRITE ( stdout , 9077 ) eext
          IF ( tefield )            WRITE( stdout, 9064 ) etotefield
          IF ( gate )               WRITE( stdout, 9065 ) etotgatefield
          IF ( lda_plus_u )         WRITE( stdout, 9066 ) eth
          IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
          IF ( okpaw ) THEN
            WRITE( stdout, 9067 ) epaw
            ! Detailed printout of PAW energy components, if verbosity is high
            IF(iverbosity>0)THEN
            WRITE( stdout, 9068) SUM(etot_cmp_paw(:,1,1)), &
                                 SUM(etot_cmp_paw(:,1,2)), &
                                 SUM(etot_cmp_paw(:,2,1)), &
                                 SUM(etot_cmp_paw(:,2,2)), &
            SUM(etot_cmp_paw(:,1,1))+SUM(etot_cmp_paw(:,1,2))+ehart, &
            SUM(etot_cmp_paw(:,2,1))+SUM(etot_cmp_paw(:,2,2))+etxc-etxcc
            ENDIF
          ENDIF
          !
          IF ( lrism ) THEN
             WRITE( stdout, 9902 ) esol
             IF ( ABS( vsol ) > eps8 ) WRITE( stdout, 9903 ) vsol
          END IF
          !
          ! ... With Grandcanonical-SCF (GC-SCF) or Fictitious charge particle (FCP),
          ! ... etot is the grand potential energy Omega = E - muN, where -muN is
          ! ... the potentiostat contribution.
          !
          IF ( lgcscf .AND. (.NOT. gcscf_ignore_mun) ) THEN
             WRITE( stdout, 9072 ) egrand
          END IF
          !
          IF ( lfcp ) WRITE( stdout, 9072 ) fcp_mu*tot_charge
          !
       ELSE IF ( conv_elec ) THEN
          !
          WRITE( stdout, 9081 ) etot
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          IF ( lgauss ) then
             WRITE( stdout, 9070 ) demet
             WRITE( stdout, 9170 ) etot-demet
          ENDIF
          !
          IF( lgcscf ) WRITE( stdout, 9181 ) tot_charge
          !
       ELSE
          !
          WRITE( stdout, 9080 ) etot
          IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9083 ) dr2
          ELSE
             WRITE( stdout, 9084 ) dr2
          END IF
          !
          IF( lgcscf ) WRITE( stdout, 9182 ) tot_charge, ef*RYTOEV, &
                                             ABS( ef - gcscf_mu ) * RYTOEV
          !
       ENDIF
       !
#if defined (__ENVIRON)
       IF (use_environ) CALL print_environ_energies('PW')
#endif
       !
       IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
       !
       IF ( noncolin .AND. domag ) &
            WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
       !
       IF ( i_cons == 3 .OR. i_cons == 4 )  &
            WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
       IF ( i_cons /= 0 .AND. i_cons < 4 ) &
            WRITE( stdout, 9073 ) lambda
       !
       FLUSH( stdout )
       !
       RETURN
       !
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9060 FORMAT(/'     The total energy is the sum of the following terms:' )
9061 FORMAT(/'     The total energy is F=E-TS. E is the sum of the following terms:' )
9062 FORMAT( '     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9064 FORMAT( '     electric field correction =',F17.8,' Ry' )
9065 FORMAT( '     gate field correction     =',F17.8,' Ry' ) ! TB
9066 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9068 FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
            /'      -> PAW hartree energy PS =',F17.8,' Ry' &
            /'      -> PAW xc energy AE      =',F17.8,' Ry' &
            /'      -> PAW xc energy PS      =',F17.8,' Ry' &
            /'      -> total E_H with PAW    =',F17.8,' Ry' &
            /'      -> total E_XC with PAW   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' )
9082 FORMAT( '     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9083 FORMAT( '     estimated scf accuracy    <',0PF17.8,' Ry' )
9084 FORMAT( '     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT( '     total all-electron energy =',0PF17.6,' Ry' )
9170 FORMAT( '     internal energy E=F+TS    =',0PF17.8,' Ry' )
9181 FORMAT(                                                  &
            /'!    total charge of GC-SCF    =',0PF17.8,' e' )
9182 FORMAT(                                                  &
            /'     total charge of GC-SCF    =',0PF17.8,' e'  &
            /'     the Fermi energy          =',0PF17.8,' eV' &
            /'                        (error :',0PF17.8,' eV)')

9902 FORMAT( '     solvation energy (RISM)   =',F17.8,' Ry' )
9903 FORMAT( '     level-shifting contrib.   =',F17.8,' Ry' )
  END SUBROUTINE print_energies


!! begin modification



SUBROUTINE itddft()

#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE constants,            ONLY : au_sec
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE wvfct,                ONLY : nbnd, npwx, wg, g2kin, current_k, et
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE uspp,                 ONLY : nkb, vkb, deeq
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE ldaU 
  USE buffers,              ONLY : get_buffer, save_buffer
  USE mp_pools,             ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE mp,                   ONLY : mp_max  
  USE wavefunctions_gpum,   ONLY : evc_d, using_evc_d, using_evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto  
  USE wvfct_gpum,           ONLY : et_d, using_et, using_et_d
  USE uspp_init,            ONLY : init_us_2
  IMPLICIT NONE
  REAL(dp) :: dt
  REAL(dp), SAVE :: max_Ke, last_etot, original_dtau
  INTEGER :: npw, ik, ibnd, jbnd, istep, m, kdmx, i, minindx
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), hpsi_d(:,:)
  REAL(DP), ALLOCATABLE :: dt_d(:)
  LOGICAL :: firstIter = .TRUE.
  external h_psi, s_psi, s_1psi, lr_sm1_psi, h_psi_gpu
#if defined(__CUDA)
  attributes(DEVICE) :: hpsi_d, dt_d
#endif

  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )
  CALL using_becp_auto(2)

    IF(firstIter) THEN
      DO ik = 1, nks
        CALL g2_kin(ik)
        npw = ngk(ik)
        DO m = 1, npw  
          max_Ke = MAX( max_Ke, g2kin(m) )
        END DO
      END DO 
      CALL mp_max( max_Ke, inter_pool_comm )
      original_dtau = dtau      
      firstIter = .FALSE.
    END IF

    IF (use_gpu) THEN
      ALLOCATE(hpsi_d(npol*npwx,nbnd)) 
    ELSE
      ALLOCATE(hpsi(npol*npwx,nbnd))
    END IF  

    IF( idum > 2 ) THEN
      IF( etot > last_etot ) THEN
        dtau = MAX( dtau * 0.5, 0.05 * original_dtau )
      ELSE  
        dtau = MIN( dtau * 1.25, 1 * original_dtau )
      END IF
    END IF
    last_etot = etot


    IF( g_dependant_dtau ) THEN
      ALLOCATE(dt_d(npwx))   
      minindx = max(npwx/20,10)
      dt_d(minindx+1:npwx) = dtau / g2kin(minindx+1:npwx)
      dt_d(1:minindx) = dtau / g2kin(minindx)
    END IF    

    dt = dtau/max_Ke   

    write(stdout,*) "idum, dtau, etot", idum, real(dtau), etot

    DO ik = 1, nks

     npw = ngk(ik)

     current_k = ik
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     IF ( use_gpu ) THEN
        CALL g2_kin_gpu( ik )
     ELSE
        CALL g2_kin( ik )
     END IF

     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb, .true. )

     IF ( nks > 1 .OR. lelfield ) &
          CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     IF ( nks > 1 .OR. lelfield ) CALL using_evc(2)

     IF ( nks > 1 .AND. lda_plus_u .AND. (Hubbard_projectors .NE. 'pseudo') ) &
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )


        ! Propogate

        IF (use_gpu) THEN

           CALL using_evc_d(1)

           CALL h_psi_gpu( npwx, npw, nbnd, evc_d, hpsi_d )

           CALL compute_band_energy_itddft(hpsi,hpsi_d)

           IF( g_dependant_dtau ) THEN

              CALL using_et_d(0)

              !$cuf kernel do(2) <<<*,*>>>
              DO ibnd = 1, nbnd
                 DO i = 1, npw
                    evc_d(i,ibnd) = ( evc_d(i,ibnd) - dt_d(i) * hpsi_d(i,ibnd) ) &
                                    / (1 - dt_d(i) * et_d(ibnd,ik) )
                 END DO
              END DO
                    
           ELSE

              !$cuf kernel do(2) <<<*,*>>>
              DO ibnd = 1, nbnd
                 DO i = 1, npw
                    evc_d(i,ibnd) = evc_d(i,ibnd) - dt * hpsi_d(i,ibnd) 
                 END DO
              END DO

           END IF

           CALL orthonormalize_gpu( evc_d, npw, nbnd )   

        ELSE

           CALL using_evc(1)

           CALL h_psi(npwx, npw, nbnd, evc, hpsi)
           CALL compute_band_energy_itddft(hpsi,hpsi_d)

           IF( g_dependant_dtau ) THEN

              CALL using_et(0)

              DO ibnd = 1, nbnd
                 DO i = 1, npw
                    evc(i,ibnd) = ( evc(i,ibnd) - dt_d(i) * hpsi(i,ibnd) ) &
                                  / ( 1 - dt_d(i) * et(ibnd,ik) )
                 END DO 
              END DO

           ELSE

              evc = evc - dt * hpsi

           END IF

!if(gamma_only) then
           CALL orthonormalize( evc, npw, nbnd )
!else           
!           CALL orthonormalize2(npw)
!end if           

        END IF

      
     CALL using_evc(0)
     IF ( nks > 1 .OR. lelfield ) &
          CALL save_buffer ( evc, nwordwfc, iunwfc, ik )


    END DO


    IF (use_gpu) THEN
      DEALLOCATE(hpsi_d)   
    ELSE
      DEALLOCATE(hpsi)
    END IF  

    IF( g_dependant_dtau ) DEALLOCATE(dt_d) 
  CALL deallocate_bec_type( becp )
  CALL using_becp_auto(2)

END SUBROUTINE itddft




#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )


SUBROUTINE orthonormalize( psi, npw, n_orbitals )

  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wvfct,                ONLY : npwx, current_k
  USE gvect,                ONLY : gstart  
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm
  USE gvect,                ONLY : gstart
  USE mp,                   ONLY : mp_sum  

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npw, n_orbitals
  COMPLEX(DP), INTENT(INOUT) :: psi(:,:)
  INTEGER :: kdim, kdmx, j, m, l, ik, m_start, m_end, npw2, npwx2 
  COMPLEX(DP) :: prod, spsi(npol*npwx), products(n_orbitals)
  REAL(DP) :: rand, psi_norm,  lagrange(n_orbitals)
  external h_psi, s_psi, s_1psi, lr_sm1_psi

  IF(gamma_only) THEN

      npw2 = 2 * npw
      npwx2 = 2 * npwx

      DO m = 1, n_orbitals

      ! ... calculate S|psi>
      CALL s_1psi( npwx, npw, psi(1,m), spsi )

      call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
      lagrange = 0.d0
      if(m_start.le.m_end) &
      CALL DGEMV( 'T', npw2, m_end-m_start+1, 2.D0, psi(1,m_start), npwx2, spsi, 1, 0.D0, lagrange(m_start), 1 )
      CALL mp_sum( lagrange( 1:m ), inter_bgrp_comm )
      IF ( gstart == 2 ) lagrange(1:m) = lagrange(1:m) - psi(1,1:m) * spsi(1)

      CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )

      psi_norm = lagrange(m)

      DO j = 1, m - 1
         psi(:,m)  = psi(:,m) - lagrange(j) * psi(:,j)
         psi_norm = psi_norm - lagrange(j)**2
      END DO

      psi_norm = SQRT( psi_norm )

      psi(:,m) = psi(:,m) / psi_norm

      IF ( gstart == 2 ) psi(1,m) = CMPLX( DBLE(psi(1,m)), 0.D0 ,kind=DP)

      END DO

   ELSE

      kdim = npwx * npol
      kdmx = npwx * npol

      CALL calbec ( npw, vkb, psi, becp, n_orbitals)

      DO m = 1, n_orbitals

         ! ... calculate S|psi>
      IF(okvan) THEN
         CALL s_1psi( npwx, npw, psi(1,m), spsi(1))
      ELSE
         spsi(:) = psi(:,m)
      END IF

         call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end

         products = ZERO

         if(m_start.le.m_end) &
         CALL ZGEMV( 'C', kdim, m_end-m_start+1, ONE, psi(1,m_start), kdmx, spsi(1), 1, ZERO, products(m_start), 1 )

      !   CALL ZGEMV( 'C', kdim, m, ONE, psi(1,1), kdmx, spsi(1), 1, ZERO, products(1), 1 )

         CALL mp_sum( products( 1:m ), inter_bgrp_comm )
         CALL mp_sum( products( 1:m ), intra_pool_comm )

         prod = DBLE( products(m) )

         DO j = 1, m - 1
            psi(:,m)  =  psi(:,m) - products(j) * psi(:,j)
            prod = prod - ( DBLE( products(j) )**2 + AIMAG( products(j) )**2 )
         END DO

         prod = SQRT( prod )
         psi(:,m) = psi(:,m) / prod

      END DO

   END IF    



END SUBROUTINE orthonormalize






SUBROUTINE orthonormalize2(npw)

  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE uspp,                 ONLY : nkb, vkb
  USE noncollin_module,     ONLY : npol
  USE wvfct,                ONLY : npwx, current_k
  USE wavefunctions,        ONLY : evc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npw
  INTEGER :: kdim, kdmx, j, m, l, ik, count, ortho_time
  COMPLEX(DP) :: prod, zdotc, s1psi(npol*npwx), products(nbnd)
  REAL(DP) :: rand
  external h_psi, s_psi, s_1psi



         kdim = npwx * npol
         kdmx = npwx * npol


        CALL calbec ( npw, vkb, evc, becp, nbnd)

        DO m = 1, nbnd

          ! ... calculate S|psi>
          !
        IF(okvan) THEN
          !
          CALL s_1psi( npwx, npw, evc(1,m), s1psi(1))
          !
        ELSE
          !
          s1psi(:) = evc(:,m)
          !
        END IF

       !   call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
          products = ZERO

       !   if(m_start.le.m_end) &
       !   CALL ZGEMV( 'C', kdim, m_end-m_start+1, ONE, evc(1,m_start), kdmx, &
       !                             s1psi(1,m), 1, ZERO, products(m_start), 1 )

           CALL ZGEMV( 'C', kdim, m, ONE, evc(1,1), kdmx, s1psi(1), 1, ZERO, products(1), 1 )

        !   CALL mp_sum( products( 1:m ), inter_bgrp_comm )
          !
          CALL mp_sum( products( 1:m ), intra_pool_comm )
          !
          prod = DBLE( products(m) )
          !
          DO j = 1, m - 1
             !
             evc(:,m)  =  evc(:,m) - products(j) * evc(:,j)
             ! 
             !
             prod = prod - ( DBLE( products(j) )**2 + AIMAG( products(j) )**2 )
             !
          END DO
          !
          prod = SQRT( prod )
          !
          evc(:,m) = evc(:,m) / prod
          !   
          !
       END DO


END SUBROUTINE orthonormalize2




SUBROUTINE orthonormalize_gpu( psi_d, npw, n_orbitals )

  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type 
  USE uspp,                 ONLY : nkb, vkb, okvan 
  USE noncollin_module,     ONLY : npol 
  USE wvfct,                ONLY : npwx, current_k 
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm 
  USE gvect,                ONLY : gstart 
  USE mp,                   ONLY : mp_sum 
  USE io_files,             ONLY : iunwfc, nwordwfc 
#if defined(__CUDA)
  USE cudafor 
  USE cublas 
#endif  
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN) :: psi_d(:,:)
  INTEGER, INTENT(IN) :: npw, n_orbitals
#if defined(__CUDA)  
  INTEGER :: kdim, kdmx, j, m, l, i, count, ortho_time, m_start, m_end, npw2, npwx2
  REAL(DP) :: rand, psi_norm
  REAL(DP) :: psi_norm_d
  REAL(DP),    ALLOCATABLE :: lagrange(:), lagrange_d(:)
  COMPLEX(DP), ALLOCATABLE :: products(:), products_d(:)
  COMPLEX(DP), ALLOCATABLE :: psi_aux(:)  
  COMPLEX(DP), ALLOCATABLE :: spsi_d(:)
  COMPLEX(DP) :: spsi1
  external h_psi, s_psi, s_1psi, lr_sm1_psi, s_1psi_gpu
  attributes(DEVICE) :: lagrange_d, spsi_d, psi_norm_d, psi_d, products_d

  ALLOCATE( spsi_d( npwx ) )
  
   ! code below taken from rcgdiagg_gpu.f90

   IF( gamma_only ) THEN
write(stdout,*) "in orthonormalize gpu gamma only"
   npw2 = 2 * npw
   npwx2 = 2 * npwx

   ALLOCATE( psi_aux( n_orbitals ) )
   ALLOCATE( lagrange( n_orbitals ) )  
   ALLOCATE( lagrange_d( n_orbitals ) )

   DO m = 1, n_orbitals
     ! ... calculate S|psi>
     !
     IF(okvan) THEN
       !
       CALL s_1psi_gpu( npwx, npw, psi_d(1,m), spsi_d )
       !
     ELSE
       !
       !$cuf kernel do(1) <<<*,*>>>
       DO i = 1, npw
         !
         spsi_d(i) = psi_d(i,m)
         !
       END DO       
       !
     END IF 
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     ! 
     call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
     lagrange = 0.d0
     if(m_start.le.m_end) &
     CALL DGEMV_gpu( 'T', npw2, m_end-m_start+1, 2.D0, psi_d(1,m_start), npwx2, spsi_d, 1, 0.D0, lagrange_d(m_start), 1 )
     if(m_start.le.m_end) lagrange( m_start:m_end ) = lagrange_d( m_start:m_end )
     !print *, 'lagrange ', lagrange(1:m)
     CALL mp_sum( lagrange( 1:m ), inter_bgrp_comm )
     !
     IF ( gstart == 2 ) THEN
        psi_aux(1:m) = psi_d(1,1:m)
        spsi1 = spsi_d(1)
        lagrange(1:m) = lagrange(1:m) - psi_aux(1:m) * spsi1
     END IF
     !     
     CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
     ! 
     psi_norm = lagrange(m)
     lagrange_d(1:m) = lagrange(1:m)
     !
     DO j = 1, m - 1
        !     
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, npwx
           psi_d(i,m)  = psi_d(i,m) - lagrange_d(j) * psi_d(i,j)
        END DO
        !
        !print *, 'psi_norm ', j, psi_norm
        psi_norm = psi_norm - lagrange(j)**2
        !
     END DO
     !
     psi_norm_d = SQRT( psi_norm )
     !print *, 'psi_norm 178', psi_norm
     !     
     !$cuf kernel do(1) <<<*,*>>>
     DO i = 1, npwx
        psi_d(i,m) = psi_d(i,m) / psi_norm_d
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF (i == 1) THEN
          IF ( gstart == 2 ) psi_d(1,m) = CMPLX( DBLE(psi_d(1,m)), 0.D0 ,kind=DP)
        END IF
     END DO
     !
   END DO

   DEALLOCATE( psi_aux )
   DEALLOCATE( lagrange_d )
   DEALLOCATE( lagrange )  

   ELSE
write(stdout,*) "in orthonormalize gpu and not in gamma only"

   IF ( npol == 1 ) THEN
      kdim = npw
      kdmx = npwx
   ELSE
      kdim = npwx*npol
      kdmx = npwx*npol
   END IF



   ALLOCATE( products( n_orbitals ) )  
   ALLOCATE( products_d( n_orbitals ) )

   DO m = 1, n_orbitals
     ! ... calculate S|psi>
     !
     IF(okvan) THEN
       !
   !    CALL s_1psi_gpu( npwx, npw, psi_d(1,m), spsi_d )
       !
     ELSE
       !
       !$cuf kernel do(1) <<<*,*>>>
       DO i = 1, npw
         !
         spsi_d(i) = psi_d(i,m)
         !
       END DO       
       !
     END IF 
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     ! 
     call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
     products = 0.d0


!psi_aux = spsi_d
!print *, 'print(sum(spsi_d))', sum(psi_aux)
!do l = m_start, m_end
!DO i = 1, npw
!psi_aux(i) = psi_d(i,l)
!END DO  
!print *, 'print(sum(psi_d(:,l)))', l, sum(psi_aux)
!end do

     if(m_start.le.m_end) &
     CALL ZGEMV_gpu( 'C', kdim, m_end-m_start+1, 1.D0, psi_d(1,m_start), kdmx, spsi_d, 1, 0.D0, products_d(m_start), 1 )
     if(m_start.le.m_end) products( m_start:m_end ) = products_d( m_start:m_end )
!     print *, 'products ', products(1:m)
     CALL mp_sum( products( 1:m ), inter_bgrp_comm )
     !     
     CALL mp_sum( products( 1:m ), intra_bgrp_comm )
     ! 
     psi_norm = products(m)
!     print *, 'products(m)', products(m)
     products_d(1:m) = products(1:m)
     !
     DO j = 1, m - 1
        !     
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, npwx
           psi_d(i,m)  = psi_d(i,m) - products_d(j) * psi_d(i,j)
        END DO
        !
!        print *, 'psi_norm ', j, psi_norm
        psi_norm = psi_norm - ( DBLE( products(j) )**2 + AIMAG( products(j) )**2 )
        !
     END DO
     !
     psi_norm_d = SQRT( psi_norm )
!     print *, 'psi_norm 178', psi_norm
     !     
     !$cuf kernel do(1) <<<*,*>>>
     DO i = 1, npwx
        psi_d(i,m) = psi_d(i,m) / psi_norm_d
     END DO
     !
   END DO

   DEALLOCATE( products_d )
   DEALLOCATE( products )  

   END IF

  DEALLOCATE( spsi_d )
  

#endif

END SUBROUTINE orthonormalize_gpu





SUBROUTINE compute_overlap( psi1, psi2 )
  ! Compute <psi_i|psi_j> = et
  USE noncollin_module,     ONLY : npol
  USE wvfct,                ONLY : nbnd, et, current_k
  USE mp_pools,             ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE io_global,            ONLY : stdout, ionode
  USE kinds,                ONLY : DP
  USE uspp,                 ONLY : nkb, vkb
  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_max, mp_bcast
  USE control_flags,        ONLY : gamma_only
  USE mp_bands_util,        ONLY : gstart
  IMPLICIT NONE
  INTEGER :: npw
  COMPLEX(DP), INTENT(IN) :: psi1(:,:), psi2(:,:) 
  INTEGER :: iband, jband, i, j, nbnd1, nbnd2
  REAL(DP) :: zdotc
  COMPLEX(DP), ALLOCATABLE :: overlap(:,:) !  <psi_i|psi_j>
  external s_psi
  COMPLEX(dp), ALLOCATABLE ::  spsi2(:,:)


   nbnd1 = size(psi1(1,:))
   nbnd2 = size(psi2(1,:))
   npw   = size(psi2(:,1))

   ALLOCATE(spsi2(npw,nbnd2))
!   CALL calbec ( npw, vkb, psi2, becp, nbnd2)
!   CALL s_psi ( npw, npw, nbnd2, psi2, spsi2 )

   ALLOCATE(overlap(nbnd1,nbnd2))
   spsi2 = psi2

      IF(gamma_only) THEN

      DO iband = 1, nbnd1

         IF(okvan) THEN
           CALL s_1psi( npwx, npw, psi2(1,iband), spsi2(1,iband))
         ELSE
           spsi2(:,iband) = psi2(:,iband)
         END IF

         DO jband = 1, nbnd2

            overlap(iband,jband) = 2 * SUM( CONJG( psi1(:, iband) ) * spsi2(:, jband) ) 

            IF( gstart == 2 ) overlap(jband,iband) = overlap(jband,iband) - CONJG( psi1(1, jband) ) * spsi2(1, iband)

         END DO

      END DO

      ELSE

      DO iband = 1, nbnd1

         IF(okvan) THEN
           CALL s_1psi( npwx, npw, psi2(1,iband), spsi2(1,iband))
         ELSE
           spsi2(:,iband) = psi2(:,iband)
         END IF

         DO jband = 1, nbnd2

            overlap(iband,jband) = SUM( CONJG( psi1(:, jband) ) * spsi2(:, iband) ) 

         END DO

      END DO

      END IF

    CALL mp_sum ( overlap, intra_bgrp_comm )

    WRITE(stdout,*) "Overlap Matrix:"


   ! WRITE(stdout,*) overlap

  do i=1,nbnd1
     write(stdout,'(5x,20f10.6)') (dreal(overlap(i,j)),j=1,nbnd2)
  end do

DEALLOCATE(spsi2)

END SUBROUTINE compute_overlap




  SUBROUTINE compute_band_energy_itddft(hpsi,hpsi_d)
  ! Compute <psi_i|H|psi_i> = et
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif  
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, et, current_k, et
  USE mp_pools,             ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE klist,                ONLY : nks, xk, ngk, igk_k    
  USE noncollin_module,     ONLY : npol
  USE uspp,                 ONLY : nkb, vkb
  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type, bec_type
  USE ldaU,                 ONLY : wfcU, Hubbard_projectors
  USE lsda_mod,             ONLY : current_spin, isk
  USE ener,                 ONLY : eband
  USE wavefunctions_gpum,   ONLY : evc_d, using_evc_d, using_evc, using_evc
  USE wvfct_gpum,           ONLY : using_et, using_wg, et_d, using_et, using_et_d, using_wg_d
  USE control_flags,        ONLY : use_gpu
  USE becmod_subs_gpum,     ONLY : using_becp_auto  
  IMPLICIT NONE
  INTEGER :: ibnd, count, band_energy_time, ik, npw
  COMPLEX(DP), EXTERNAL :: zdotc
  COMPLEX(DP), ALLOCATABLE :: products(:), products_d(:)
  COMPLEX(DP), INTENT(IN) :: hpsi(:,:), hpsi_d(:,:)
  REAL(DP) :: aux
  LOGICAL :: xdft_ham

#if defined(__CUDA)
   attributes(DEVICE) :: hpsi_d, products_d
#endif

      !CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )  

   IF( .NOT. use_gpu ) THEN 

      ALLOCATE(products(nbnd))
      CALL using_evc(0)
      CALL using_wg(0)
      CALL using_et(2)

      npw = ngk(current_k)

      IF(gamma_only) THEN

         DO ibnd = 1, nbnd

            products(ibnd) = 2 * zdotc (npw*npol, evc(1, ibnd), 1,  hpsi(1, ibnd), 1) 

            IF( gstart == 2 ) products(ibnd) = products(ibnd) - CONJG(evc(1, ibnd)*hpsi(1, ibnd))

         END DO

      ELSE

         DO ibnd = 1, nbnd

            products(ibnd) = zdotc (npw*npol, evc(1, ibnd), 1,  hpsi(1, ibnd), 1)

         END DO

      END IF
      
      CALL mp_sum ( products, intra_bgrp_comm )
      
      et(:,current_k) = products(:)

      DEALLOCATE(products)

   ELSE

      ALLOCATE(products_d(nbnd))
      CALL using_evc_d(0)
      CALL using_wg_d(0)
      CALL using_et_d(2)

      npw = ngk(current_k)

      products_d = 0.0_DP

      IF(gamma_only) THEN
         
         !$cuf kernel do(1) <<<*,*>>>
         DO ibnd = 1, nbnd
            DO i = 1, npw*npol
               products_d(ibnd) = products_d(ibnd) + 2.D0 * DCONJG(evc_d(i, ibnd)) * hpsi_d(i, ibnd)
            END DO
            IF (gstart == 2) products_d(ibnd) = products_d(ibnd) - DCONJG(evc_d(1, ibnd)) * hpsi_d(1, ibnd) 
         END DO

      ELSE

         !$cuf kernel do(1) <<<*,*>>>
         DO ibnd = 1, nbnd
            DO i = 1, npw*npol
               products_d(ibnd) = products_d(ibnd) + DCONJG(evc_d(i, ibnd)) * hpsi_d(i, ibnd)
            END DO
         END DO

      END IF

      CALL mp_sum ( products_d, intra_bgrp_comm )

      !$cuf kernel do(1) <<<*,*>>>
      DO ibnd = 1, nbnd      
         et_d(ibnd ,current_k) = products_d(ibnd )
      END DO

      DEALLOCATE(products_d)

   END IF

   !CALL deallocate_bec_type ( becp )

  END SUBROUTINE compute_band_energy_itddft

!! end modification


  !
END SUBROUTINE electrons_scf
!
!----------------------------------------------------------------------------
FUNCTION exxenergyace( )
  !--------------------------------------------------------------------------
  !! Compute exchange energy using ACE
  !
  USE kinds,           ONLY : DP
  USE buffers,         ONLY : get_buffer
  USE exx,             ONLY : vexxace_gamma, vexxace_k, domat
  USE klist,           ONLY : nks, ngk
  USE wvfct,           ONLY : nbnd, npwx, current_k
  USE lsda_mod,        ONLY : lsda, isk, current_spin
  USE io_files,        ONLY : iunwfc, nwordwfc
  USE mp_pools,        ONLY : inter_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  USE control_flags,   ONLY : gamma_only
  USE wavefunctions,   ONLY : evc
  !
  USE wavefunctions_gpum, ONLY : using_evc
  !
  IMPLICIT NONE
  !
  REAL(DP) :: exxenergyace
  !! computed energy
  !
  ! ... local variables
  !
  REAL(DP) :: ex
  INTEGER :: ik, npw
  !
  domat = .TRUE.
  exxenergyace=0.0_dp
  !
  CALL using_evc(0)
  !
  DO ik = 1, nks
     npw = ngk (ik)
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF (nks > 1) CALL using_evc(2)
     IF (gamma_only) THEN
        CALL vexxace_gamma( npw, nbnd, evc, ex )
     ELSE
        CALL vexxace_k( npw, nbnd, evc, ex )
     ENDIF
     exxenergyace = exxenergyace + ex
  ENDDO
  !
  CALL mp_sum( exxenergyace, inter_pool_comm )
  !
  domat = .FALSE.
  !
END FUNCTION exxenergyace


!! begin modification

!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_sm1_psi_2 (recalculate, ik, lda, n, m, psi, spsi)
  !----------------------------------------------------------------------------
  !
  ! This subroutine applies the S^{-1} matrix to m wavefunctions psi
  ! and puts the results in spsi.
  ! See Eq.(13) in B. Walker and R. Gebauer, JCP 127, 164106 (2007).
  ! Requires the products of psi with all beta functions
  ! in array becp(nkb,m) (calculated in h_psi or by ccalbec).
  !
  ! INPUT: recalculate   Decides if the overlap of beta 
  !                      functions is recalculated or not.
  !                      This is needed e.g. if ions are moved 
  !                      and the overlap changes accordingly.
  !        ik            k point under consideration
  !        lda           leading dimension of arrays psi, spsi
  !        n             true dimension of psi, spsi
  !        m             number of states psi
  !        psi           the wavefunction to which the S^{-1} 
  !                      matrix is applied
  ! OUTPUT: spsi = S^{-1}*psi
  !
  ! Original routine written by R. Gebauer
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE uspp,             ONLY : okvan, vkb, nkb, qq_nt
  USE uspp_param,       ONLY : nh, upf
  USE ions_base,        ONLY : ityp,nat,ntyp=>nsp
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : noncolin, npol
  USE matrix_inversion
  !
  ! mod !
  USE uspp_init,            ONLY : init_us_2  
  ! mod !

  IMPLICIT NONE
  LOGICAL, INTENT(in)      :: recalculate
  INTEGER, INTENT(in)      :: lda, n, m, ik
  COMPLEX(DP), INTENT(in)  :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(out) :: spsi(lda*npol,m)
  !
  LOGICAL :: recalc
  !


  REAL (DP),    ALLOCATABLE :: bbg(:,:)      ! nkb, nkb)
  ! for gamma_only     
  COMPLEX (DP), ALLOCATABLE :: bbk(:,:,:)    ! nkb, nkb, nks)
  ! for k points
  COMPLEX (DP), ALLOCATABLE :: bbnc(:,:,:,:) ! nkb, nkb, nspin_mag, nks)
  ! for the noncollinear case
  ! bbg = < beta^N_i | beta^P_j > 
  ! bbg/bbk/bbnc are the scalar products of beta functions 
  ! localized on atoms N and P.




  CALL start_clock( 'lr_sm1_psi_2' )
  !
  recalc = recalculate
  !
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSEIF (noncolin) THEN
     CALL errore( 'lr_sm1_psi_2', 'Noncollinear case is not implemented', 1 )
  ELSE
     CALL sm1_psi_k()
  ENDIF
  !
  CALL stop_clock( 'lr_sm1_psi_2' )
  !
  RETURN
  !
CONTAINS
  ! 
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    ! gamma_only version
    !
    ! Note: 1) vkb must be computed before calling this routine
    !       2) the array bbg must be deallocated somewhere
    !          outside of this routine.
    !
    USE becmod,   ONLY : bec_type,becp,calbec
    USE realus,   ONLY : real_space, invfft_orbital_gamma, initialisation_level, &
                         fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                         v_loc_psir, s_psir_gamma
 !   USE lrus,     ONLY : bbg
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    ! counters
    REAL(DP), ALLOCATABLE :: ps(:,:)
    LOGICAL, SAVE :: first_entry = .true.
    !
    ! Initialize spsi : spsi = psi



  !  REAL (DP),    ALLOCATABLE :: bbg(:,:)      ! nkb, nkb)
    ! for gamma_only   



    !
    CALL ZCOPY( lda * npol * m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B but
    ! use the ones which were already calculated and saved (if recalc=.false.).
    !
    IF (first_entry) THEN
       !
       IF (allocated(bbg)) DEALLOCATE(bbg)
       first_entry = .false.
       recalc = .true.
       !
    ENDIF
    !
    IF (.not.allocated(bbg)) recalc = .true.
    IF (recalc .and. allocated(bbg)) DEALLOCATE(bbg)
    !
    IF (recalc) THEN
       !
       ALLOCATE(bbg(nkb,nkb))
       bbg = 0.d0
       !
       ! The beta-functions vkb must have beed calculated elsewhere. 
       ! So there is no need to call the routine init_us_2.
       !
       ! Calculate the coefficients B_ij defined by Eq.(15).
       ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
       !
       CALL calbec (n, vkb, vkb, bbg(:,:), nkb)
       !
       ALLOCATE(ps(nkb,nkb))
       !
       ps(:,:) = 0.d0
       !
       ! Calculate the product of q_nm and B_ij of Eq.(16).
       !
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + &
                               & qq_nt(ih,jh,nt) * bbg(jkb,ii)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO
       !
       ! Add identity to q_nm * B_ij [see Eq.(16)].
       ! ps = (1 + q*B)
       !
       DO ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + 1.d0
       ENDDO
       !
       ! Invert matrix: (1 + q*B)^{-1}
       !
       CALL invmat( nkb, ps )
       !
       ! Use the array bbg as a workspace in order to save memory
       !
       bbg(:,:) = 0.d0
       !
       ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
       !
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            bbg(ii,jkb) = bbg(ii,jkb) &
                                    & - ps(ii,ikb) * qq_nt(ih,jh,nt)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO
       DEALLOCATE(ps)
    ENDIF
    !
    ! Compute the product of the beta-functions vkb with the functions psi,
    ! and put the result in becp.
    ! becp(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n,vkb,psi,becp,m)
    !
    ! Use the array ps as a workspace
    ALLOCATE(ps(nkb,m))
    ps(:,:) = 0.D0
    !
    ! Now let us apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : ps = lambda * <beta|psi>
    ! Here, lambda = bbg, and <beta|psi> = becp%r
    !
    call DGEMM( 'N','N',nkb,m,nkb,1.d0,bbg(:,:),nkb,becp%r,nkb,0.d0,ps,nkb)
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta>
    !
    call DGEMM('N','N',2*n,m,nkb,1.d0,vkb,2*lda,ps,nkb,1.d0,spsi,2*lda)
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_gamma
  
  SUBROUTINE sm1_psi_k()
    !-----------------------------------------------------------------------
    !
    ! k-points version
    ! Note: the array bbk must be deallocated somewhere
    ! outside of this routine.
    !
    USE becmod,   ONLY : bec_type,becp,calbec
    USE klist,    ONLY : nks, xk, ngk, igk_k
  !  USE lrus,     ONLY : bbk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ik1
    ! counters
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)



    !-mod-!
  !  COMPLEX (DP), ALLOCATABLE :: bbk(:,:,:)    ! nkb, nkb, nks)
    ! for k points
    !-mod-!



    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*npol*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! If this is the first entry, we calculate and save the coefficients B from Eq.(15)
    ! B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
    ! If this is not the first entry, we do not recalculate the coefficients B 
    ! (they are already allocated) but use the ones which were already calculated 
    ! and saved (if recalc=.false.).
    !
    IF (.NOT.ALLOCATED(bbk)) recalc = .true.
    IF (recalc .AND. ALLOCATED(bbk)) DEALLOCATE(bbk)
    !
    ! If recalc=.true. we (re)calculate the coefficients lambda_nm defined by Eq.(16).
    !
    IF (recalc) THEN
       !
       ALLOCATE(bbk(nkb,nkb,nks))
       bbk = (0.d0,0.d0)
       !
       ALLOCATE(ps(nkb,nkb))
       !
       DO ik1 = 1, nks
          !
          ! Calculate beta-functions vkb for a given k point.
          !
          CALL init_us_2(ngk(ik1),igk_k(:,ik1),xk(1,ik1),vkb)
          !
          ! Calculate the coefficients B_ij defined by Eq.(15).
          ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
          !
          CALL zgemm('C','N',nkb,nkb,ngk(ik1),(1.d0,0.d0),vkb, &
                 & lda,vkb,lda,(0.d0,0.d0),bbk(1,1,ik1),nkb)
          !
#if defined(__MPI)
          CALL mp_sum(bbk(:,:,ik1), intra_bgrp_comm)
#endif
          !
          ps(:,:) = (0.d0,0.d0)
          !
          ! Calculate the product of q_nm and B_ij of Eq.(16).
          !
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               ps(ikb,ii) = ps(ikb,ii) + &
                                & bbk(jkb,ii, ik1) * qq_nt(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO
          !
          ! Add an identity to q_nm * B_ij [see Eq.(16)].
          ! ps = (1 + q*B)
          !
          DO ii = 1, nkb
             ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
          ENDDO
          !
          ! Invert matrix: (1 + q*B)^{-1}
          !
          CALL invmat( nkb, ps )
          ! 
          ! Use the array bbk as a work space in order to save memory.
          !
          bbk(:,:,ik1) = (0.d0,0.d0)
          !
          ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
          !
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               bbk(ii,jkb,ik1) = bbk(ii,jkb,ik1) &
                                           & - ps(ii,ikb) * qq_nt(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO
          !
       ENDDO ! loop on k points
       DEALLOCATE(ps)
       !
    ENDIF
    !
    IF (n.NE.ngk(ik)) CALL errore( 'sm1_psi_k', &
                    & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k point 'ik'.
    !
    CALL init_us_2(n,igk_k(:,ik),xk(1,ik),vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n,vkb,psi,becp,m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb,m))
    ps(:,:) = (0.d0,0.d0)
    !
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's do this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi> 
    !
    DO ibnd=1,m
       DO jkb=1,nkb
          DO ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd) + bbk(jkb,ii,ik) * becp%k(ii,ibnd)
          ENDDO
       ENDDO
    ENDDO
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta>
    !
    CALL ZGEMM( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
               & lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_k

END SUBROUTINE lr_sm1_psi_2








!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
#define ONE  ( 1._DP, 0._DP )
#define MONE (-1._DP, 0._DP )
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_k_gpu( npwx, npw, nbnd, npol, psi_d, hpsi_d, spsi_d, e, &
                           uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
    USE io_global,            ONLY : stdout
  !
  USE util_param,     ONLY : DP, eps16
  USE mp,            ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands_util, ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id
  USE device_fbuff_m,         ONLY : buffer => dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi_d (npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi_d(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: eigen
  LOGICAL,     INTENT(IN)    :: reorder
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  LOGICAL                  :: eigen_
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  INTEGER                  :: jblock_start, jblock_end
  INTEGER                  :: ibnd_start, ibnd_end
  INTEGER                  :: jbnd_start, jbnd_end
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
  INTEGER                  :: buf_start, buf_end, buf_size
  !
  ! ... device variables 
  !
  INTEGER :: ii, jj, kk, info
  COMPLEX(DP), ALLOCATABLE :: phi_d(:,:), hphi_d(:,:), sphi_d(:,:)
#if defined (__CUDA)
  attributes(device) :: psi_d, hpsi_d, spsi_d
  attributes(device) :: phi_d, hphi_d, sphi_d
#endif
  !
  COMPLEX(DP), ALLOCATABLE :: sc_d(:), sc2_d(:,:)
#if defined (__CUDA)
  attributes(device) :: sc_d, sc2_d
#endif 
  !

write(stdout,*) "hello from gram schmidt k gpu"

  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  eigen_ = eigen
  !
  IF ( reorder ) THEN
     !
     eigen_ = .TRUE.
     !
  END IF
  !
  write(stdout,*) "gram schmidt 1", nbsize
  nblock = nbnd / nbsize
  write(stdout,*) "gram schmidt 2", nbsize  
  IF ( MOD( nbnd, nbsize ) /= 0 ) nblock = nblock + 1
  !
  CALL divide( inter_bgrp_comm, nblock, iblock_start, iblock_end )
  !
  IF ( my_bgrp_id >= nblock ) THEN
     !
     iblock_start = nblock + 1
     iblock_end   = nblock
     !
  END IF
  !
  ALLOCATE( phi_d ( kdmx, nbnd ) )
  IF ( eigen_ ) ALLOCATE( hphi_d( kdmx, nbnd ) )
  IF ( uspp )   ALLOCATE( sphi_d( kdmx, nbnd ) )
  !
  ALLOCATE( owner_bgrp_id( nblock ) )
  !
!$cuf kernel do(2)
  DO ii = 1, kdmx  
    DO jj = 1, nbnd
      phi_d(ii, jj) = ZERO 
    END DO 
  END DO   
  !
  IF ( eigen_ ) THEN 
!$cuf kernel do(2)
    DO ii = 1, kdmx  
      DO jj = 1, nbnd
        hphi_d(ii, jj) = ZERO 
      END DO 
    END DO   
  END IF  
  !
  IF ( uspp )   THEN 
!$cuf kernel do(2)
    DO ii = 1, kdmx  
      DO jj = 1, nbnd
        sphi_d(ii, jj) = ZERO 
      END DO 
    END DO   
  END IF  
  !
  ! ... Set owers of blocks
  !
  owner_bgrp_id = 0
  !
  DO iblock = 1, nblock
     !
     IF ( iblock_start <= iblock .AND. iblock <= iblock_end ) &
     owner_bgrp_id(iblock) = my_bgrp_id
     !
  END DO
  !
  CALL mp_max( owner_bgrp_id, inter_bgrp_comm )
  !
  ! ... Set initial : |phi_j> = |psi_j>
  !
  CALL ZCOPY_gpu( kdmx * nbnd, psi_d(1,1), 1, phi_d(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY_gpu( kdmx * nbnd, hpsi_d(1,1), 1, hphi_d(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY_gpu( kdmx * nbnd, spsi_d(1,1), 1, sphi_d(1,1), 1 )
  !
  !
  ! ... Allocate buffers 
  !
  ALLOCATE (sc_d(nbsize)) 
  IF ( nbnd .GT. nbsize)  ALLOCATE (sc2_d(nbsize, nbnd - nbsize)) 
  !
  !
  ! ... Blocking loop
  !
  DO iblock = 1, nblock
     !
     ! ... Orthogonalize diagonal block by standard Gram-Schmidt
     !
     ibnd_start = ( iblock - 1 ) * nbsize + 1
     ibnd_end   = MIN( iblock * nbsize, nbnd )
     !
     IF ( owner_bgrp_id(iblock) == my_bgrp_id ) &
         CALL gram_schmidt_diag( ibnd_start, ibnd_end )
     !
     !
     ! ... Bcast diagonal block
     !
     CALL mp_bcast( phi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( eigen_ ) &
     CALL mp_bcast( hphi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( uspp ) &
     CALL mp_bcast( sphi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     ! ... Project off-diagonal block outside of diagonal block
     !
     jblock_start = MAX( iblock_start, iblock + 1 )
     jblock_end   = iblock_end
     !
     jbnd_start = ( jblock_start - 1 ) * nbsize + 1
     jbnd_end   = MIN( jblock_end * nbsize, nbnd )
     !
     IF ( jblock_start <= jblock_end .AND. jbnd_start <= jbnd_end ) &
     !
     !
     CALL project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
     !
     !
  END DO
  !
  ! ... Buffer Realese
  !
  DEALLOCATE (sc_d) 
  IF (nbnd .GT. nbsize) DEALLOCATE (sc2_d) 
  !
  !
  ! ... Copy psi <- phi
  !
  CALL ZCOPY_gpu( kdmx * nbnd, phi_d(1,1), 1, psi_d(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY_gpu( kdmx * nbnd, hphi_d(1,1), 1, hpsi_d(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY_gpu( kdmx * nbnd, sphi_d(1,1), 1, spsi_d(1,1), 1 )
  !
  ! ... Calculate energy eigenvalues
  !
  IF ( eigen_ ) CALL energyeigen_gpu( )
  !
  ! ... Sort wave functions
  !
  IF ( reorder ) CALL sort_vectors_gpu( )
  !
  DEALLOCATE( phi_d )
  IF ( eigen_ ) DEALLOCATE( hphi_d )
  IF ( uspp )   DEALLOCATE( sphi_d )
  !
  DEALLOCATE( owner_bgrp_id )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE gram_schmidt_diag( ibnd_start, ibnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    !
    INTEGER                  :: ibnd
    REAL(DP)                 :: norm
    COMPLEX(DP), EXTERNAL    :: ZDOTC_gpu
    !
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL ZGEMV_gpu( 'C', kdim, ibnd - ibnd_start, ONE, phi_d(1,ibnd_start), kdmx, &
                         spsi_d(1,ibnd), 1, ZERO, sc_d(1), 1 )
             !
          ELSE
             !
             CALL ZGEMV_gpu( 'C', kdim, ibnd - ibnd_start, ONE, phi_d(1,ibnd_start), kdmx, &
                         psi_d(1,ibnd), 1, ZERO, sc_d(1), 1 )
             !
          END IF
          !
          !
          CALL mp_sum( sc_d, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, phi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, phi_d(1,ibnd), 1 )
          !
          !
          IF ( eigen_ ) &
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, hphi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, hphi_d(1,ibnd), 1 )
          !
          IF ( uspp ) &
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, sphi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, sphi_d(1,ibnd), 1 )
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          norm = DBLE( ZDOTC_gpu( kdim, phi_d(1,ibnd), 1, sphi_d(1,ibnd), 1 ) )
          !
       ELSE
          !
          norm = DBLE( ZDOTC_gpu( kdim, phi_d(1,ibnd), 1, phi_d(1,ibnd), 1 ) )
          !
       END IF
       !
       CALL mp_sum( norm, intra_bgrp_comm )
       !
       norm = SQRT( MAX( norm, 0._DP ) )
       !
       IF ( norm < eps16 ) &
       CALL errore( ' gram_schmidt_k ', ' vectors are linear dependent ', 1 )
       !
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, phi_d(1,ibnd), 1 )
       !
       IF ( eigen_ ) &
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, hphi_d(1,ibnd), 1 )
       !
       IF ( uspp ) &
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, sphi_d(1,ibnd), 1 )
       !
    END DO
    !
    !
    RETURN
    !
  END SUBROUTINE gram_schmidt_diag
  !
  !
  SUBROUTINE project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    INTEGER, INTENT(IN)  :: jbnd_start, jbnd_end
    !
    INTEGER                  :: ibnd_size
    INTEGER                  :: jbnd_size 
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1  
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL gpu_ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi_d(1,ibnd_start), kdmx, &
                   spsi_d(1,jbnd_start), kdmx, ZERO, sc2_d(1,1), nbsize )
       !
    ELSE
       !
       CALL gpu_ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi_d(1,ibnd_start), kdmx, &
                   psi_d(1,jbnd_start), kdmx, ZERO, sc2_d(1,1), nbsize )
       !
    END IF
    !
    CALL mp_sum( sc2_d, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, phi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, phi_d(1,jbnd_start), kdmx )
    !
    IF ( eigen_ ) &
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, hphi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, hphi_d(1,jbnd_start), kdmx )
    !
    IF ( uspp ) &
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, sphi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, sphi_d(1,jbnd_start), kdmx )
    !
    RETURN
    !
  END SUBROUTINE project_offdiag_gpu
  !
  !
  SUBROUTINE energyeigen_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd, ibnd_start, ibnd_end
    !
    COMPLEX(DP), EXTERNAL :: ZDOTC_gpu
    !
    ! ... <psi_i| H |psi_i>
    !
    e(:) = 0._DP
    !
    CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       e(ibnd) = DBLE( ZDOTC_gpu( kdim, psi_d(1,ibnd), 1, hpsi_d(1,ibnd), 1 ) )
       !
    END DO
    !
    CALL mp_sum( e(ibnd_start:ibnd_end), intra_bgrp_comm )
    CALL mp_sum( e, inter_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE energyeigen_gpu
  !
  !
  SUBROUTINE sort_vectors_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER  :: ibnd
    INTEGER  :: nswap
    REAL(DP) :: e0
    !
10  nswap = 0
    !
    DO ibnd = 2, nbnd
       !
       IF ( e(ibnd) < e(ibnd-1) ) THEN
          !
          nswap = nswap + 1
          !
          e0        = e(ibnd)
          e(ibnd)   = e(ibnd-1)
          e(ibnd-1) = e0
          !
          CALL ZSWAP_gpu( kdim, psi_d(1,ibnd), 1, psi_d(1,ibnd-1), 1 )
          !
          IF ( eigen_ ) &
          CALL ZSWAP_gpu( kdim, hpsi_d(1,ibnd), 1, hpsi_d(1,ibnd-1), 1 )
          !
          IF ( uspp ) &
          CALL ZSWAP_gpu( kdim, spsi_d(1,ibnd), 1, spsi_d(1,ibnd-1), 1 )
          !
       END IF
       !
    END DO
    !
    IF ( nswap > 0 ) GOTO 10
    !
    RETURN
    !
  END SUBROUTINE sort_vectors_gpu
  !
  !
END SUBROUTINE gram_schmidt_k_gpu


!! end modification
