! radiation_liquid_optics_socrates.F90 - SOCRATES method for parameterizing liquid droplet optics
!
! Copyright (C) 2014-2016 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_liquid_optics_nielsen

  implicit none

  ! SOCRATES (Edwards-Slingo) parameterizes info on the dependence of
  ! the scattering properties in each band on effective radius in
  ! terms of 16 coefficients
  integer, parameter :: NLiqOpticsCoeffsNielsenSW = 8

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties using a
  ! parameterization consisting of Pade approximants from the
  ! SOCRATES (Edwards-Slingo) code
  subroutine calc_liq_optics_nielsen(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! Effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um
    !real(jprb) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    re_um= re*1e6_jprb
    inv_re_um = 1.0_jprb / re_um

    od=lwp_gm_2 * coeff(1:nb,1)*(re_um**coeff(1:nb,2))
    scat_od=od*(coeff(1:nb,3)+coeff(1:nb,4)*re_um)
    g=coeff(1:nb,5)+coeff(1:nb,6)*re_um+(coeff(1:nb,7)*EXP(coeff(1:nb,8)*re_um))

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_nielsen

end module radiation_liquid_optics_nielsen
