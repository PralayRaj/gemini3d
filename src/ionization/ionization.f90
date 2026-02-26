module ionization

use phys_consts, only: elchrg, lsp, kb, mn, re, pi, wp, lwave, debug
use ionize_fang, only: fang2008, fang2010, fang2010_spectrum
!! we need the unperturbed msis temperatures to apply the simple chapman theory used by this module
use grid, only: lx1,lx2,lx3
use meshobj, only: curvmesh
use timeutils, only: ymd2doy

implicit none (type, external)
private
public :: ionrate_fang, ionrate_glow98, eheating, photoionization

interface
  module subroutine glow_run(W0,PhiWmWm2,date_doy,UTsec,xf107,xf107a,xlat,xlon,alt,nn,Tn,ns,Ts,&
    ionrate,eheating,iver)
    real(wp), dimension(:), intent(in) :: W0,PhiWmWm2,alt,Tn
    real(wp), dimension(:,:), intent(in) :: nn,ns,Ts
    real(wp), dimension(:,:), intent(inout) :: ionrate
    !! intent(out)
    real(wp), dimension(:), intent(inout) :: eheating, iver
    !! intent(out)
    real(wp), intent(in) :: UTsec, xlat, xlon, xf107, xf107a
    integer, intent(in) :: date_doy
  end subroutine glow_run
end interface

contains
  function photoionization(t,ymd,UTsec,x,nn,chi,f107,f107a,gavg,Tninf,Iinf)
    !------------------------------------------------------------
    !-------COMPUTE PHOTOIONIZATION RATES PER SOLOMON ET AL, 2005
    !------------------------------------------------------------
    real(wp), intent(in) :: t
    integer, intent(in), dimension(3) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    !real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), dimension(:,:,:), intent(in) :: chi
    real(wp), parameter :: chi0 = pi/2._wp
    real(wp), parameter :: dchi = 5._wp * pi / 180._wp   ! 5° smooth transition
    real(wp) :: w
    real(wp) :: Fday, Fnight
    real(wp), intent(in) :: f107,f107a
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), intent(in) :: Iinf
    integer, parameter :: ll=23     !number of wavelength bins
    integer :: il,isp,ix1,ix2,ix3
    !Line bin indices
    integer, parameter :: i_heii=9      ! 29.0–32.0 nm
    integer, parameter :: i_hei =11     ! 54.0–65.0 nm
    integer, parameter :: i_lyb =21     ! 98.7–102.7 nm
    integer, parameter :: i_lya =23     ! 105.0–121.6 nm
    real(wp), dimension(ll) :: lambda1,lambda2,sigmaO,sigmaN2,sigmaO2
    !Strobel line cross sections (convert 10^-18 cm^2 to m^2)
    real(wp), parameter :: cs = 1e-22_wp
    real(wp), parameter :: sigO_hei  = 10._wp*cs
    real(wp), parameter :: sigN2_hei = 23._wp*cs
    real(wp), parameter :: sigO2_hei = 23._wp*cs

    real(wp), parameter :: sigO_heii  = 13._wp*cs
    real(wp), parameter :: sigN2_heii = 24._wp*cs
    real(wp), parameter :: sigO2_heii = 11._wp*cs

    real(wp), parameter :: sigO_lyb  = 0._wp*cs
    real(wp), parameter :: sigN2_lyb = 0._wp*cs
    real(wp), parameter :: sigO2_lyb = 0.98_wp*cs

    real(wp), parameter :: sigO_lya  = 0._wp*cs
    real(wp), parameter :: sigN2_lya = 0._wp*cs
    real(wp), parameter :: sigO2_lya = 0._wp*cs
    
    real(wp) :: sigO_eff, sigN2_eff, sigO2_eff

    real(wp), dimension(ll) :: brN2i,brN2di,brO2i,brO2di,pepiO,pepiN2i,pepiN2di,pepiO2i,pepiO2di
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: bigX,y,Chfn
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol,nN2col,nO2col
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: phototmp
    real(wp) :: H
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),ll) :: Iflux
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),lsp-1) :: photoionization    !don't need a separate rate for electrons
    logical, dimension(size(nn,1),size(nn,2),size(nn,3)) :: is_night, is_day
    real(wp) :: alt_km, sza_deg, logF, sza
    real(wp) :: a, b, c
    ! --- DEBUG: choose one diagnostic point (set these to a known night location)
    integer, parameter :: ix1_dbg = 20, ix2_dbg = 1, ix3_dbg = 1
    logical :: do_dbg
    real(wp) :: Fn_fit, tau_night, att_night, prod_scale
    real(wp), dimension(size(chi,1),size(chi,2),size(chi,3)) :: chi_eff
    real(wp) :: tauO, tauN2, tauO2, tau_hei, tau_heii, tau_lyb, tau_lya, tau_vert, tau_slant
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol_slant,nN2col_slant,nO2col_slant
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol_vert ,nN2col_vert ,nO2col_vert
    
    !WAVELENGTH BIN BEGINNING AND END (THIS IDEALLY WOULD BE DATA STATEMENTS OR SOME KIND OF STRUCTURE THAT DOESN'T GET REASSIGNED AT EVERY CALL).  Actually all of these array assignments are static...
    lambda1=[0.05, 0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 65.0, &
        79.8, 79.8, 79.8, 91.3, 91.3, 91.3, 97.5, 98.7, 102.7, 105.0]*1e-9
    lambda2=[0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 79.8, 79.8, &
         91.3, 91.3, 91.3, 97.5, 97.5, 97.5, 98.7, 102.7, 105.0, 121.6]*1e-9
    
    
    !TOTAL ABSORPTION CROSS SECTIONS
    sigmaO=[0.0023, 0.0170, 0.1125, 0.1050, 0.3247, 1.3190, 3.7832, 6.0239, &
         7.7205, 10.7175, 13.1253, 8.5159, 4.7889, 3.0031, 4.1048, 3.7947, &
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]*1e-18*1e-4         !convert to m^2
    sigmaN2=[0.0025, 0.0201, 0.1409, 1.1370, 0.3459, 1.5273, 5.0859, 9.9375, &
        11.7383, 19.6514, 23.0931, 23.0346, 54.5252, 2.1434, 13.1062, 71.6931, &
        2.1775, 14.4390, 115.257, 2.5465, 0.0, 0.0, 0.0]*1e-18*1e-4
    sigmaO2=[0.0045, 0.034, 0.2251, 0.2101, 0.646, 2.6319, 7.6283, 13.2125, &
        16.8233, 20.3066, 27.0314, 23.5669, 24.9102, 10.4980, 10.9075, 13.3122, &
        13.3950, 14.4042, 32.5038, 18.7145, 1.6320, 1.15, 0.0]*1e-18*1e-4
    
    
    !BRANCHING RATIOS
    brN2i=[0.040,0.040,0.040,0.040, 0.717, 0.751, 0.747, 0.754, 0.908, 0.996, 1.0, 0.679,  &
        0.429, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    brN2di=[0.96, 0.96,0.96,0.96,0.282, 0.249, 0.253, 0.246, 0.093, 0.005, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    brO2i=[0.0, 0.0, 0.0, 0.0, 0.108, 0.347, 0.553, 0.624, 0.649, 0.759, 0.874, 0.672, 0.477, &
        0.549, 0.574, 0.534, 0.756, 0.786, 0.620, 0.830, 0.613, 0.0, 0.0]
    brO2di=[1.0, 1.0, 1.0, 1.0, 0.892, 0.653, 0.447, 0.376, 0.351, 0.240, 0.108, 0.001, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    
    !PHOTOELECTRON TO DIRECT PRODUCTION RATIOS
    pepiO=[217.12, 50.593, 23.562, 71.378, 4.995, 2.192, 1.092, 0.694, 0.418, &
        0.127, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    pepiN2i=[263.99, 62.57, 25.213, 8.54, 6.142, 2.288, 0.786, 0.324, 0.169, 0.031, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    pepiN2di=[78.674, 18.310, 6.948, 2.295, 1.647, 0.571, 0.146, 0.037, 0.008, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    pepiO2i=[134.69, 32.212, 13.309, 39.615, 2.834, 1.092, 0.416, 0.189, 0.090, 0.023, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    pepiO2di=[76.136, 17.944, 6.981, 20.338, 1.437, 0.521, 0.163, 0.052, 0.014, 0.001, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        
    !Column densities (O, N2, O2)
!     call compute_column_density(nn(:,:,:,1), chi, x, Tninf, gavg, mn(1), nOcol)
!     call compute_column_density(nn(:,:,:,2), chi, x, Tninf, gavg, mn(2), nN2col)
!     call compute_column_density(nn(:,:,:,3), chi, x, Tninf, gavg, mn(3), nO2col)

   !DEBUG    
    chi_eff = min(chi, chi0)   ! chi0 = pi/2
    call compute_column_density(nn(:,:,:,1), chi_eff, x, Tninf, gavg, mn(1), nOcol)
    call compute_column_density(nn(:,:,:,2), chi_eff, x, Tninf, gavg, mn(2), nN2col)
    call compute_column_density(nn(:,:,:,3), chi_eff, x, Tninf, gavg, mn(3), nO2col)
    
    call compute_column_density_vertical(nn(:,:,:,1), x, nOcol_vert)
    call compute_column_density_vertical(nn(:,:,:,2), x, nN2col_vert)
    call compute_column_density_vertical(nn(:,:,:,3), x, nO2col_vert)
    
    
    
!     !O COLUMN DENSITY
!     H=kB*Tninf/mn(1)/gavg                         !scalar scale height
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H          !a reduced altitude
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
!     !      Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1._wp-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nOcol=nn(:,:,:,1)*H*Chfn
!     
!     
!     !N2 COLUMN DENSITY
!     H=kB*Tninf/mn(2)/gavg     !all of these temp quantities need to be recomputed for eacb neutral species being ionized
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nN2col=nn(:,:,:,2)*H*Chfn
!     
!     
!     !O2 COLUMN DENSITY
!     H=kB*Tninf/mn(3)/gavg
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nO2col=nn(:,:,:,3)*H*Chfn
      

      is_night = (chi(1:lx1,1:lx2,1:lx3) >= pi/2._wp)
      is_day   = .not. is_night
      
      !DEBUG
!       write(*,*) 'DEBUG: lx1 =', lx1
!       write(*,*) 'DEBUG: lx2 =', lx2
!       write(*,*) 'DEBUG: lx3 =', lx3

      ! Flux computation
      
    do ix3 = 1, lx3
    do ix2 = 1, lx2
    do ix1 = 1, lx1

      alt_km = x%alt(ix1, ix2, ix3) / 1000._wp
      sza    = chi(ix1, ix2, ix3)

      ! Classify day/night using SZA
!       if (sza < pi/2._wp) then
!       ! Daytime: compute attenuated flux using column densities
!         !print *, 'Daytime SZA =', sza * 180.0_wp / pi
!         do il = 1, ll
!           Iflux(ix1,ix2,ix3,il) = Iinf(ix1,ix2,ix3,il) * exp( - &
!             ( sigmaO(il)  * nOcol(ix1,ix2,ix3)  + &
!               sigmaN2(il) * nN2col(ix1,ix2,ix3) + &
!               sigmaO2(il) * nO2col(ix1,ix2,ix3) ) )
!         end do
!       elseif (sza >= pi/2._wp .and. sza < pi) then
!       !Nighttime: use quadratic SZA fits from Strobel 1974 (for Ly-alpha, Ly-beta, He I, He II only)
! !       print *, 'Nighttime SZA =', sza * 180.0_wp / pi
!         Iflux(ix1,ix2,ix3,:) = 0._wp
!         ! He II
!         Iflux(ix1,ix2,ix3,i_heii) = get_nightflux("heii",  alt_km, sza) 
!         Iflux(ix1,ix2,ix3,i_hei) = get_nightflux("hei",    alt_km, sza) 
!         Iflux(ix1,ix2,ix3,i_lyb) = get_nightflux("lybeta", alt_km, sza)
!         Iflux(ix1,ix2,ix3,i_lya) = get_nightflux("lyalpha",alt_km, sza)
!         
!         do il = 1, ll
!           if (il /= i_hei .and. il /= i_heii .and. il /= i_lyb .and. il /= i_lya) then 
!           Iflux(ix1,ix2,ix3,il) = 0._wp
!           end if
!         end do
!       else
!       !--- Invalid SZA: zero everything
!       !print *, 'Invalid SZA at (',ix1,ix2,ix3,') SZA =', sza * 180.0_wp / pi
!         Iflux(ix1,ix2,ix3,:) = 0._wp
!       end if
    do il = 1, ll
    !day flux
   Fday = Iinf(ix1,ix2,ix3,il) * exp( - &
       ( sigmaO(il)  * nOcol(ix1,ix2,ix3)  + &
         sigmaN2(il) * nN2col(ix1,ix2,ix3) + &
         sigmaO2(il) * nO2col(ix1,ix2,ix3) ) )
         
    ! night flux
   Fnight = 0._wp
   !print *, 'Nighttime SZA =', sza * 180.0_wp / pi
   if (il == i_heii)  then
   tau_heii = sigO_heii  * nOcol_vert(ix1,ix2,ix3)  + sigN2_heii * nN2col_vert(ix1,ix2,ix3) + &
              sigO2_heii * nO2col_vert(ix1,ix2,ix3)
     
   Fnight = get_nightflux("heii",  alt_km, sza) * exp(-tau_heii)
   end if
   if (il == i_hei)  then
   tau_hei = sigO_hei  * nOcol_vert(ix1,ix2,ix3)  + sigN2_hei * nN2col_vert(ix1,ix2,ix3) + &
              sigO2_hei * nO2col_vert(ix1,ix2,ix3)
   Fnight = get_nightflux("hei",   alt_km, sza) * exp(-tau_hei) 
   end if
   if (il == i_lyb)   then
   tau_lyb = sigO_lyb  * nOcol_vert(ix1,ix2,ix3)  + sigN2_lyb * nN2col_vert(ix1,ix2,ix3) + &
              sigO2_lyb * nO2col_vert(ix1,ix2,ix3)
   Fnight = get_nightflux("lybeta",alt_km, sza) * exp(-tau_lyb)
   end if
   if (il == i_lya)   then
   tau_lya = sigO_lya  * nOcol_vert(ix1,ix2,ix3)  + sigN2_lya * nN2col_vert(ix1,ix2,ix3) + &
              sigO2_lya * nO2col_vert(ix1,ix2,ix3)
   Fnight = get_nightflux("lyalpha",alt_km, sza) * exp(-tau_lya) 
   end if
    
   !smooth transition from day to night
   !w = 1._wp / (1._wp + exp((sza - chi0)/dchi))
   w = 0.5_wp * (1._wp - tanh((sza - chi0)/(2._wp*dchi))) ! more stable option

   Iflux(ix1,ix2,ix3,il) = w * Fday + (1._wp - w) * Fnight
   
   tau_vert = sigO_hei  * nOcol_vert(ix1,ix2,ix3)  + sigN2_hei * nN2col_vert(ix1,ix2,ix3) + &
              sigO2_hei * nO2col_vert(ix1,ix2,ix3)
   tau_slant = sigO_hei  * nOcol(ix1,ix2,ix3)  + sigN2_hei * nN2col(ix1,ix2,ix3) + &
              sigO2_hei * nO2col(ix1,ix2,ix3)
!    write(*,*) 'tau_vert=',tau_vert,' tau_slant=',tau_slant
   
!    !DEBUG
!       do_dbg = (ix1==ix1_dbg .and. ix2==ix2_dbg .and. ix3==ix3_dbg)
! 
!     if (do_dbg) then
!     if (il == i_heii) then
!       Fn_fit = get_nightflux("heii", alt_km, sza)
! 
!       tau_night = sigO_heii  * nOcol(ix1,ix2,ix3)  + &
!                   sigN2_heii * nN2col(ix1,ix2,ix3) + &
!                   sigO2_heii * nO2col(ix1,ix2,ix3)
! 
!       att_night = exp(-tau_night)
! 
!       !prod_scale = nn(ix1,ix2,ix3,1) * (Fn_fit*att_night) * sigO_heii
!       !prod_scale = nn(ix1,ix2,ix3,1) * (Fn_fit) * sigO_heii
! 
! !       write(*,'(A,3I6,A,F8.2,A,F8.2)') 'DBG pt=',ix1,ix2,ix3, &
! !      '  alt_km=',alt_km,'  SZAdeg=',sza*180._wp/pi
! ! 
! !       write(*,'(A,ES12.3,A,ES12.3)') 'DBG HeI : Fday=',Fday, &
! !      '  Fnight_used=',Fnight
! ! 
! !       write(*,'(A,F10.6,A,ES12.3)') 'DBG HeI : w=',w, &
! !      '  Iflux=',Iflux(ix1,ix2,ix3,il)
! ! 
! !         ! Production scale using the actual flux driving the rate
! !       prod_scale = nn(ix1,ix2,ix3,1) * Iflux(ix1,ix2,ix3,il) * sigO_heii
! !       write(*,'(A,ES12.3)') 'DBG HeI : n*Iflux*sigO =', prod_scale
!     end if
! 
!     if (il == i_hei) then
!       Fn_fit = get_nightflux("hei", alt_km, sza)
! 
!       tau_night = sigO_hei  * nOcol(ix1,ix2,ix3)  + &
!                   sigN2_hei * nN2col(ix1,ix2,ix3) + &
!                   sigO2_hei * nO2col(ix1,ix2,ix3)
! 
!       att_night = exp(-tau_night)
! 
!       !prod_scale = nn(ix1,ix2,ix3,1) * (Fn_fit*att_night) * sigO_hei
!       !prod_scale = nn(ix1,ix2,ix3,1) * (Fn_fit) * sigO_hei
! 
! !       write(*,'(A,3I6,A,F8.2,A,F8.2)') 'DBG pt=',ix1,ix2,ix3, &
! !      '  alt_km=',alt_km,'  SZAdeg=',sza*180._wp/pi
! ! 
! !       write(*,'(A,ES12.3,A,ES12.3)') 'DBG HeII: Fday=',Fday, &
! !      '  Fnight_used=',Fnight
! ! 
! !       write(*,'(A,F10.6,A,ES12.3)') 'DBG HeII: w=',w, &
! !      '  Iflux=',Iflux(ix1,ix2,ix3,il)
! ! 
! !       prod_scale = nn(ix1,ix2,ix3,1) * Iflux(ix1,ix2,ix3,il) * sigO_hei
! !       write(*,'(A,ES12.3)') 'DBG HeII: n*Iflux*sigO =', prod_scale
!            
!     tauO  = sigO_hei  * nOcol(ix1,ix2,ix3)
!     tauN2 = sigN2_hei * nN2col(ix1,ix2,ix3)
!     tauO2 = sigO2_hei * nO2col(ix1,ix2,ix3)
! 
! !     write(*,'(A,ES12.3,A,ES12.3,A,ES12.3)') 'DBG cols: nOcol=', nOcol(ix1,ix2,ix3), &
! !      ' nN2col=', nN2col(ix1,ix2,ix3), ' nO2col=', nO2col(ix1,ix2,ix3)
! !     write(*,'(A,ES12.3,A,ES12.3,A,ES12.3,A,ES12.3)') 'DBG tau: tauO=', tauO, &
! !      ' tauN2=', tauN2, ' tauO2=', tauO2, ' tauSum=', tauO+tauN2+tauO2
!     end if
!     end if

    end do


    end do
    end do
    end do


    ! Ionization rate calculation
    photoionization = 0._wp

    do il = 1, ll

    where (is_day)
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
        photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
        photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
        photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
    end where

    ! Night time
    if (il == i_heii) then    
    where (is_night)
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux(:,:,:,il)*sigO_heii
        photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigN2_heii
        photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigN2_heii*brN2di(il)*(1 + pepiN2di(il))
        photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigO2_heii
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigO2_heii*brO2di(il)*(1 + pepiO2di(il))
    end where
    elseif (il == i_hei) then
    where (is_night)
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux(:,:,:,il)*sigO_hei
        photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigN2_hei
        photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigN2_hei*brN2di(il)*(1 + pepiN2di(il))
        photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigO2_hei
        photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigO2_hei*brO2di(il)*(1 + pepiO2di(il))
    end where
    elseif (il == i_lyb) then
    where (is_night)
        ! Ly-beta: only O2 has nonzero cross section
        photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigO2_lyb
    end where
    elseif (il == i_lya) then
    where (is_night)
        ! Ly-alpha: this contributes zero
    end where
    end if

    end do

    photoionization(:,:,:,2) = 0._wp
    photoionization(:,:,:,6) = 0._wp
!     
!     
!     do isp = 1, 6
!     phototmp = photoionization(:,:,:,isp)
!     where (x%nullpts)
!     phototmp = 0._wp
!     end where
!     photoionization(:,:,:,isp) = phototmp
!     end do 

    ! Ionization rate calculation (unified, smooth-flux compatible)

!     photoionization = 0._wp
! 
!     do il = 1, ll
! 
!   ! Effective cross sections (override for Strobel line bins)
!       sigO_eff  = sigmaO(il)
!       sigN2_eff = sigmaN2(il)
!       sigO2_eff = sigmaO2(il)
! 
!     if (il == i_heii) then
!       sigO_eff  = sigO_heii
!       sigN2_eff = sigN2_heii
!       sigO2_eff = sigO2_heii
!     elseif (il == i_hei) then
!       sigO_eff  = sigO_hei
!       sigN2_eff = sigN2_hei
!       sigO2_eff = sigO2_hei
!     elseif (il == i_lyb) then
!       sigO_eff  = sigO_lyb
!       sigN2_eff = sigN2_lyb
!       sigO2_eff = sigO2_lyb
!     elseif (il == i_lya) then
!       sigO_eff  = sigO_lya
!       sigN2_eff = sigN2_lya
!       sigO2_eff = sigO2_lya
!     end if
! 
!   !PRIMARY AND SECONDARY IONIZATION RATES
!   !direct O+ production
!     photoionization(:,:,:,1) = photoionization(:,:,:,1) + &
!        nn(:,:,:,1) * Iflux(:,:,:,il) * sigO_eff * (1._wp + pepiO(il))
! 
!   !direct N2+
!     photoionization(:,:,:,3) = photoionization(:,:,:,3) + &
!        nn(:,:,:,2) * Iflux(:,:,:,il) * sigN2_eff * brN2i(il) * (1._wp + pepiN2i(il))
! 
!   !dissociative ionization of N2 leading to N+
!     photoionization(:,:,:,5) = photoionization(:,:,:,5) + &
!        nn(:,:,:,2) * Iflux(:,:,:,il) * sigN2_eff * brN2di(il) * (1._wp + pepiN2di(il))
! 
!   !direct O2+
!     photoionization(:,:,:,4) = photoionization(:,:,:,4) + &
!        nn(:,:,:,3) * Iflux(:,:,:,il) * sigO2_eff * brO2i(il) * (1._wp + pepiO2i(il))
! 
!   !dissociative ionization of O2 leading to O+
!    photoionization(:,:,:,1) = photoionization(:,:,:,1) + &
!        nn(:,:,:,3) * Iflux(:,:,:,il) * sigO2_eff * brO2di(il) * (1._wp + pepiO2di(il))
! 
!     end do
!     !direct NO+
!     photoionization(:,:,:,2) = 0._wp
!    !H+ production
!     photoionization(:,:,:,6) = 0._wp
!    
   !THERE SHOULD BE SOME CODE HERE TO ZERO OUT THE BELOW-GROUND ALTITUDES.
   where (photoionization < 0._wp)
     photoionization = 0
   end where
   do isp=1,lsp-1
      phototmp = photoionization(:,:,:,isp)
     where (x%nullpts)
      phototmp = 0._wp
     end where
    photoionization(:,:,:,isp) = phototmp
    end do
    
    contains
    
      subroutine compute_column_density(nn_species, chi, x, Tninf, gavg, mass, n_col)
      real(wp), intent(in) :: nn_species(:,:,:), chi(:,:,:), Tninf, gavg, mass
      class(curvmesh), intent(in) :: x
      real(wp), intent(out) :: n_col(:,:,:)
      real(wp) :: H
      real(wp), dimension(size(nn_species,1), size(nn_species,2), size(nn_species,3)) :: bigX, y, Chfn
      
!       real(wp), parameter :: chi_tol = 1.0e-8_wp   ! rad: treat as "chi = 0" if |chi| < chi_tol
!       real(wp) :: dz
! 
!       integer :: i1, i2, i3
! 
!       real(wp), dimension(size(nn_species,1), size(nn_species,2), size(nn_species,3)) :: n_col_vert
!       logical,  dimension(size(nn_species,1), size(nn_species,2), size(nn_species,3)) :: is_vertical
      
          H = kB*Tninf/mass/gavg                      !scalar scale height
          bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H          !a reduced altitude
          y = sqrt(bigX/2._wp)*abs(cos(chi))
          
          where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
              !      Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1._wp-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
              Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
          elsewhere
              Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
          end where
              n_col = nn_species * H * Chfn
              
!         ! ---- Vertical override for chi ~ 0 (direct integral) ----
!           is_vertical = abs(chi) < chi_tol
! 
!   ! Only do the vertical integral if at least one point needs it
!     if (any(is_vertical)) then
!         n_col_vert = 0._wp
! 
!     do i3 = 1, size(nn_species,3)
!     do i2 = 1, size(nn_species,2)
! 
!         ! top boundary ~ 0 above model top
!         n_col_vert(size(nn_species,1), i2, i3) = 0._wp
! 
!     do i1 = size(nn_species,1)-1, 1, -1
!         dz = x%alt(i1+1,i2,i3) - x%alt(i1,i2,i3)   ! [m]
!         n_col_vert(i1,i2,i3) = n_col_vert(i1+1,i2,i3) + 0.5_wp * dz * &
!                                  ( nn_species(i1,i2,i3) + nn_species(i1+1,i2,i3) )
!     end do
! 
!     end do
!     end do
! 
!     ! Replace only where chi ~ 0
!     where (is_vertical)
!       n_col = n_col_vert
!     end where
!     end if
              
    end subroutine compute_column_density
      
      !DEBUG
      subroutine compute_column_density_vertical(nn_species, x, n_col)
      real(wp), intent(in) :: nn_species(:,:,:)
      class(curvmesh), intent(in) :: x
      real(wp), intent(out) :: n_col(:,:,:)

      integer :: i1,i2,i3
      real(wp) :: dz

      n_col = 0._wp

      do i3 = 1, size(nn_species,3)
        do i2 = 1, size(nn_species,2)
          n_col(size(nn_species,1), i2, i3) = 0._wp   ! top boundary ~0 above model top

        do i1 = size(nn_species,1)-1, 1, -1
          dz = x%alt(i1+1,i2,i3) - x%alt(i1,i2,i3)   ! meters
          n_col(i1,i2,i3) = n_col(i1+1,i2,i3) + 0.5_wp * dz * &
                          ( nn_species(i1,i2,i3) + nn_species(i1+1,i2,i3) )
        end do
        end do
      end do
    end subroutine compute_column_density_vertical

    function get_nightflux(line_type, alt_km, sza_rad) result(F)
    character(len=*), intent(in) :: line_type
    real(wp), intent(in) :: alt_km, sza_rad
    real(wp) :: F, sza_deg, a, b, c, logF
    sza_deg = sza_rad * 180.0_wp / pi
    sza_deg = min(max(sza_deg, 90._wp), 180._wp)
    !print *, 'Nighttime SZA =', sza_rad * 180.0_wp / pi
    select case (trim(adjustl(line_type)))
    case ('lyalpha')
    a = 9.7e-05_wp; b = -0.0377_wp; c = 12.7900_wp
    case ('lybeta')
    a = -1.3e-4_wp; b = -0.0500_wp; c = 11.5362_wp
!     print *, 'Case Ly-Beta'
    case ('hei')
    a = -2.4e-05_wp; b = -0.0596_wp; c = 13.5588_wp
!     print *, 'Case He I'
    case ('heii')
    a = -1.3e-04_wp; b = 0.0287_wp; c = 6.0416_wp
!     print *, 'Case He II'
    case default
    F = 0._wp; return
    end select


    logF = a * sza_deg**2 + b * sza_deg + c
    F = 10.0_wp**logF*1e4 !m-2s-1
!      print *, 'Flux for', trim(line_type), 'is', F
!      print *, 'SZA is:', sza_deg
!     print *, 'a =', a, ', b =', b, ', c =', c  
    end function get_nightflux

     
!     !Solar flux at each point on the grid, i.e. the attenuated vacuum flux
!     do il=1,ll
!       do ix3=1,x%lx3
!         do ix2=1,x%lx2
!           do ix1=1,x%lx1 
!           Iflux(ix1,ix2,ix3,il)=Iinf(ix1,ix2,ix3,il)*exp(-(sigmaO(il)*nOcol(ix1,ix2,ix3)+sigmaN2(il)*nN2col(ix1,ix2,ix3)+ &
!                     sigmaO2(il)*nO2col(ix1,ix2,ix3)))
!           end do
!         end do
!       end do
!     end do   
! 
! 
!     !PRIMARY AND SECONDARY IONIZATION RATES
!     photoionization=0
!     
!     !direct O+ production
!     do il=1,ll
!       photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
!     end do
!     
!     !direct NO+
!     photoionization(:,:,:,2) = 0
!     
!     !direct N2+
!     do il=1,ll
!       photoionization(:,:,:,3)=photoionization(:,:,:,3)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
!     end do
!     
!     !dissociative ionization of N2 leading to N+
!     do il=1,ll
!       photoionization(:,:,:,5)=photoionization(:,:,:,5)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
!     end do
!     
!     !direct O2+
!     do il=1,ll
!       photoionization(:,:,:,4)=photoionization(:,:,:,4)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
!     end do
!     
!     !dissociative ionization of O2 leading to O+
!     do il=1,ll
!       photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
!     end do
!     
!     !H+ production
!     photoionization(:,:,:,6) = 0
!     
!     
!     !THERE SHOULD BE SOME CODE HERE TO ZERO OUT THE BELOW-GROUND ALTITUDES.
!     where (photoionization < 0)
!       photoionization = 0
!     end where
!     do isp=1,lsp-1
!       phototmp=photoionization(:,:,:,isp)
!     !  where (x%nullpts>0.9 .and. x%nullpts<1.1)
!       where(x%nullpts)
!         phototmp=0
!       end where
!       photoionization(:,:,:,isp) = phototmp
!     end do
end function photoionization


  pure function ionrate_fang(W0, PhiWmWm2, alt, nn, Tn, g1, flag_fang, diff_num_flux, kappa, bimax_frac, W0_char)
  real(wp), dimension(:,:), intent(in) :: W0,PhiWmWm2
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: alt,Tn
  integer, intent(in) :: flag_fang, diff_num_flux
  real(wp), intent(in) :: kappa, bimax_frac, W0_char
  real(wp), dimension(:,:,:), intent(in) :: g1
  real(wp) :: W0keV, PhiW, W0_char_keV
  real(wp), dimension(1:size(nn,1)) :: massden,meanmass
  integer :: ix2,ix3,lx2,lx3
  real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: Ptot,PO,PN2,PO2
  real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1) :: ionrate_fang


  lx2=size(nn,2)
  lx3=size(nn,3)

  !IONIZATION RATES ARE COMPUTED ON A PER-PROFILE BASIS

  !zero flux should really be check per field line
    if ( maxval(PhiWmWm2) > 0) then   !only compute rates if nonzero flux given
      do ix3=1,lx3
        do ix2=1,lx2
          !CONVERSION TO DIFFERENTIAL NUMBER FLUX
          PhiW=PhiWmWm2(ix2,ix3)*1e-3_wp/elchrg    !from mW/m^2 to eV/m^2/s
          PhiW=PhiW/1e3_wp/1e4_wp    !to keV/cm^2/s
          W0keV=W0(ix2,ix3)/1e3_wp
          W0_char_keV=W0_char/1e3_wp
    
          massden=mn(1)*nn(:,ix2,ix3,1)+mn(2)*nn(:,ix2,ix3,2)+mn(3)*nn(:,ix2,ix3,3)
          !! mass densities are [kg m^-3] as per neutral/neutral.f90 "call meters(.true.)" for MSIS.
          meanmass=massden/(nn(:,ix2,ix3,1)+nn(:,ix2,ix3,2)+nn(:,ix2,ix3,3))
          !! mean mass per particle [kg]
    
          !> TOTAL IONIZATION RATE
          !! [cm^-3 s^-1] => [m^-3 s^-1]
          select case (flag_fang)
          case (8, 2008)
            Ptot(:,ix2,ix3) = fang2008(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3)) * 1e6_wp
          case (10, 2010)
            Ptot(:,ix2,ix3) = fang2010(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3)) * 1e6_wp
          case (0) ! composite spectrum
            Ptot(:,ix2,ix3) = fang2010_spectrum(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3), &
              diff_num_flux, kappa, bimax_frac, W0_char_keV) * 1e6_wp
          case default
            error stop 'ERROR:ionization:ionrate_fang: unknown flag_fang'
          end select
        end do
      end do
    
  
      !NOW THAT TOTAL IONIZATION RATE HAS BEEN CALCULATED BREAK IT INTO DIFFERENT ION PRODUCTION RATES
      PO = 0
      PN2 = 0
      PO2 = 0
    
      where (nn(:,:,:,1) + nn(:,:,:,2) + nn(:,:,:,3) > 1e-10_wp )
              PN2 = Ptot * 0.94_wp * nn(:,:,:,2) / &
                               (nn(:,:,:,3) + 0.94_wp*nn(:,:,:,2) + 0.55_wp * nn(:,:,:,1))
      endwhere
    
      where (nn(:,:,:,2) > 1e-10_wp)
        PO2 = PN2 * 1.07_wp * nn(:,:,:,3) / nn(:,:,:,2)
        PO = PN2 * 0.59_wp * nn(:,:,:,1) / nn(:,:,:,2)
      endwhere
    
    
      !SPLIT TOTAL IONIZATION RATE PER VALLANCE JONES, 1973
      ionrate_fang(:,:,:,1) = PO + 0.33_wp * PO2
      ionrate_fang(:,:,:,2) = 0
      ionrate_fang(:,:,:,3) = 0.76_wp * PN2
      ionrate_fang(:,:,:,4) = 0.67_wp * PO2
      ionrate_fang(:,:,:,5) = 0.24_wp * PN2
      ionrate_fang(:,:,:,6) = 0
    else
      ionrate_fang(:,:,:,:) = 0
    end if
  end function ionrate_fang
  
  
  pure function eheating(nn,ionrate,ns)
    !------------------------------------------------------------
    !-------COMPUTE SWARTZ AND NISBET, (1973) ELECTRON HEATING RATES.
    !-------ION ARRAYS (EXCEPT FOR RATES) ARE EXPECTED TO INCLUDE
    !-------GHOST CELLS.
    !------------------------------------------------------------
    
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1), intent(in) :: ionrate
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns    !includes ghost cells
    
    real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: totionrate,R,avgenergy
    integer :: lx1,lx2,lx3
    
    real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: eheating
    
    lx1=size(nn,1)
    lx2=size(nn,2)
    lx3=size(nn,3)
    
    R=log(ns(1:lx1,1:lx2,1:lx3,lsp)/(nn(:,:,:,2)+nn(:,:,:,3)+0.1_wp*nn(:,:,:,1)))
    avgenergy=exp(-(12.75_wp+6.941_wp*R+1.166_wp*R**2+0.08034_wp*R**3+0.001996_wp*R**4))
    totionrate=sum(ionrate,4)
    
    eheating=elchrg*avgenergy*totionrate
  end function eheating
  
  
  subroutine ionrate_glow98(W0,PhiWmWm2,ymd,UTsec,f107,f107a,glat,glon,alt,nn,Tn,ns,Ts, &
                                 eheating, iver, ionrate)
    !! COMPUTE IONIZATION RATES USING GLOW MODEL RUN AT EACH
    !! X,Y METHOD.
    
    real(wp), dimension(:,:,:), intent(in) :: W0,PhiWmWm2
    
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec, f107, f107a
    real(wp), dimension(:,:), intent(in) :: glat,glon
    
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
    real(wp), dimension(:,:,:), intent(in) :: alt,Tn
    
    real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)), intent(inout) :: eheating
    !! intent(out)
    real(wp), dimension(1:size(nn,2),1:size(nn,3),lwave), intent(inout) :: iver
    !! intent(out)
    real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1), intent(inout) :: ionrate
    !! intent(out)
    
    integer :: ix2,ix3,lx1,lx2,lx3,date_doy
    
    lx1=size(nn,1)
    lx2=size(nn,2)
    lx3=size(nn,3)
    
    !! zero flux should really be checked per field line
    if ( maxval(PhiWmWm2) > 0) then   !only compute rates if nonzero flux given
    
      date_doy = modulo(ymd(1), 100)*1000 + ymd2doy(ymd(1), ymd(2), ymd(3))
      !! date in format needed by GLOW (yyddd)
      do ix3=1,lx3
        do ix2=1,lx2
          !W0eV=W0(ix2,ix3) !Eo in eV at upper x,y locations (z,x,y) normally
    
          if ( maxval(PhiWmWm2(ix2,ix3,:)) <= 0) then    !only compute rates if nonzero flux given *here* (i.e. at this location)
            ionrate(:,ix2,ix3,:) = 0
            eheating(:,ix2,ix3) = 0
            iver(ix2,ix3,:) = 0
          else
            !Run GLOW here with the input parameters to obtain production rates
            !GLOW outputs ion production rates in [cm^-3 s^-1]
            call glow_run(W0(ix2,ix3,:), PhiWmWm2(ix2,ix3,:), &
              date_doy, UTsec, f107, f107a, glat(ix2,ix3), glon(ix2,ix3), alt(:,ix2,ix3), &
              nn(:,ix2,ix3,:),Tn(:,ix2,ix3), ns(1:lx1,ix2,ix3,:), Ts(1:lx1,ix2,ix3,:), &
              ionrate(:,ix2,ix3,:), eheating(:,ix2,ix3), iver(ix2,ix3,:))
    !        print*, 'glow called, max ionization rate: ', maxval(ionrate(:,ix2,ix3,:))
    !        print*, 'max iver:  ',maxval(iver(ix2,ix3,:))
    !        print*, 'max W0 and Phi:  ',maxval(W0(ix2,ix3,:)),maxval(PhiWmWm2(ix2,ix3,:))
          end if
        end do !Y coordinate loop
      end do !X coordinate loop
      eheating=eheating*elchrg
    else
      ionrate(:,:,:,:)=0 !No Q for incoming electrons, no electron impact
      eheating(:,:,:)=0
      iver(:,:,:)=0
    end if
  end subroutine ionrate_glow98
end module ionization
