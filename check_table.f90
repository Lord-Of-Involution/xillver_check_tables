program check_tables
  implicit none

  integer, parameter  :: spec_dim = 2999
  integer             :: i, unit, rownum, rowmax
  real                :: spec(spec_dim)
  ! real, allocatable   :: spec_ref(:), spec(:)
  character (len=200) :: filename  
  filename = '/Users/gullo/Software/relxill/tables/xillver-comp.fits'
!Open the fits file and leave it opened 
  call open_fits(filename, unit)
!Read how many rows there are in the extension SPECTRA 
  call read_rownum(unit, rowmax)
  print *, 'rowmax', rowmax

  ! if (.not. allocated(spec_ref)) allocate(spec_ref(spec_dim))
  ! if (.not. allocated(spec)    ) allocate(spec    (spec_dim))

  !Read the spectrum
  rownum = 1 
  call read_spec(unit, rownum, spec_dim, spec)

  call close_fits(unit)
  
  do i = 1, spec_dim
     write(10, *) i, spec(i)
  enddo
  
end program check_tables



