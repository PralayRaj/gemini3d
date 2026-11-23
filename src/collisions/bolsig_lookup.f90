module bolsig_lookup

use phys_consts, only: wp, elchrg, ms, lsp
implicit none
private

public :: init_bolsig_table, bolsig_mu_lookup, bolsig_loaded

logical :: bolsig_loaded = .false.

integer :: nE_bolsig = 0
integer :: nz_bolsig = 0

real(wp), allocatable :: bolsig_Etd(:)
real(wp), allocatable :: bolsig_alt_km(:)
real(wp), allocatable :: bolsig_logmu(:,:)   ! [E/N index, altitude index]

contains


subroutine init_bolsig_table(filename)

  character(*), intent(in) :: filename
  logical :: ex
  integer :: iu, i

  print *, 'Trying to read BOLSIG table from: ', trim(filename)
  inquire(file=filename, exist=ex)

  if (.not. ex) then
    print *, 'BOLSIG table not found. Using default GEMINI mobilities.'
    print *, trim(filename)
    bolsig_loaded = .false.
    return
  end if

  open(newunit=iu, file=filename, status='old', action='read')

  read(iu,*) nE_bolsig, nz_bolsig

  allocate(bolsig_Etd(nE_bolsig))
  allocate(bolsig_alt_km(nz_bolsig))
  allocate(bolsig_logmu(nE_bolsig,nz_bolsig))

  read(iu,*) bolsig_Etd
  read(iu,*) bolsig_alt_km

  do i = 1,nE_bolsig
    read(iu,*) bolsig_logmu(i,:)
  end do

  close(iu)

  bolsig_loaded = .true.

  print *, 'Loaded BOLSIG mobility table:'
  print *, '  nE = ', nE_bolsig
  print *, '  nz = ', nz_bolsig
  print *, '  E/N range Td = ', bolsig_Etd(1), bolsig_Etd(nE_bolsig)
  print *, '  altitude range km = ', bolsig_alt_km(1), bolsig_alt_km(nz_bolsig)

end subroutine init_bolsig_table


real(wp) function bolsig_mu_lookup(EoverN_Td, alt_km) result(mu)

  real(wp), intent(in) :: EoverN_Td
  real(wp), intent(in) :: alt_km

  real(wp) :: x, z
  real(wp) :: x1, x2, z1, z2
  real(wp) :: tx, tz
  real(wp) :: f11, f12, f21, f22, flog
  integer :: iE, iz

  if (.not. bolsig_loaded) then
    mu = -1.0_wp
    return
  end if

  ! Clamp query to table range.
  x = log10(max(min(EoverN_Td, bolsig_Etd(nE_bolsig)), bolsig_Etd(1)))
  z = max(min(alt_km, bolsig_alt_km(nz_bolsig)), bolsig_alt_km(1))

  ! Find E/N bracket in log space.
  iE = 1
  do while (iE < nE_bolsig-1 .and. log10(bolsig_Etd(iE+1)) < x)
    iE = iE + 1
  end do

  ! Find altitude bracket.
  iz = 1
  do while (iz < nz_bolsig-1 .and. bolsig_alt_km(iz+1) < z)
    iz = iz + 1
  end do

  x1 = log10(bolsig_Etd(iE))
  x2 = log10(bolsig_Etd(iE+1))

  z1 = bolsig_alt_km(iz)
  z2 = bolsig_alt_km(iz+1)

  tx = (x - x1) / (x2 - x1)
  tz = (z - z1) / (z2 - z1)

  f11 = bolsig_logmu(iE,  iz)
  f12 = bolsig_logmu(iE,  iz+1)
  f21 = bolsig_logmu(iE+1,iz)
  f22 = bolsig_logmu(iE+1,iz+1)

  flog = (1.0_wp-tx)*(1.0_wp-tz)*f11 + &
         (1.0_wp-tx)*tz        *f12 + &
         tx        *(1.0_wp-tz)*f21 + &
         tx        *tz        *f22

  mu = 10.0_wp**flog

end function bolsig_mu_lookup


end module bolsig_lookup
