
!This code checks all the spectra in the table and look if there are any NaN 


program check_tables_NaN
  implicit none

  integer, parameter :: spec_dim = 2999
  integer, parameter :: nmaxp = 10000
  double precision, parameter :: ergsev  = 1.602197d-12 ! Convert eV to ergs

  integer    :: i, j, unit, rownum, rowmax, start_ind, end_ind, index_scat, index_xil
  integer    :: par_dim, tot_spec_num, cont_bad_spec, cont, &
                tot_spec_analysed
  real       :: variation, norm, sum, sum_ref, threshold, E, dE, ratio, temp
  real       :: spec(spec_dim), spec_ref(spec_dim)
  real       :: ene_lo(spec_dim), ene_hi(spec_dim)
  real       :: earx(0 : spec_dim)
  ! real       :: en(spec_dim)
  real, allocatable :: par_val(:), save_par_val(:)
  logical    :: check      
  character (len=200) :: filename, hdu_name  
  integer           , allocatable :: num_of_val_par(:)
  character (len=12), allocatable :: par_names(:)

  double precision :: wp(nmaxp), spec_dble
  real             :: Icomp_scat(nmaxp)


!Simple check to find out if the IF statement recognize a NaN 
  ! spec(2) = 0.0
  ! spec(1) =  0.0/ spec(2)
  ! write(*, *) spec(1)
  
  ! if (spec(1) .ne. spec(1)) then
  !    write(*,*) 'It is a NaN'
  ! endif
  ! if (spec(1) .eq. spec(1)) then
  !    write(*,*) 'It is NOT a NaN'
  ! endif

  ! stop

  
  filename = '/Users/gullo/Software/relxill/tables/xillverCp_v3.6.fits' !xillver-comp.fits
  write(*,*) 'Reading table: ', trim(filename), ' (press Enter to continue)'
  read(*,*)
  
  ! Emax = 100.0
  ! Emin = 0.07177
  ! dloge = log10( Emax / Emin ) / float(spec_dim)
  ! do i = 0, spec_dim
  !    earx(i) = Emin * (Emax/Emin)**(float(i)/float(spec_dim))
  !    ! write(*,*) i, earx(i)
  ! end do

!Open the fits file and leave it opened 
  call open_fits(filename, unit)

  hdu_name = 'PARAMETERS'
  call get_rownum(unit, hdu_name, rowmax)
  par_dim = rowmax
  write(*,*) 'number of parameters in the table', par_dim

  allocate(par_names(par_dim))
  allocate(num_of_val_par(par_dim))
  call get_parameters(unit, par_dim, par_names, num_of_val_par)

  call get_energy(unit, spec_dim, ene_lo, ene_hi)

  do i = 1, spec_dim 
     earx(i - 1) =  ene_lo(i)
  enddo
  earx(spec_dim) = ene_hi(spec_dim)
  
  allocate(par_val(par_dim))
  allocate(save_par_val(par_dim))
  
  hdu_name = 'SPECTRA'
  call get_rownum(unit, hdu_name, rowmax)
  tot_spec_num = rowmax

!range between 20keV and 50keV (1713 - 1991)  
!range between 25keV and 35keV (1782 - 1883)  
  start_ind = 1713
  end_ind   = 1991
  save_par_val = 0.0
  cont_bad_spec = 0
  tot_spec_analysed = 0
  threshold = 1e2
  cont  = 0

  
  open (10, file='bad_xill_spectra.qdp')
  write(10,*) 'skip on'
  ! open (20, file='xill_spectra_parameters_ratio.txt')
  ! write(20,*) '  Spec number      Gamma        A_Fe             log_xi            kTe              Dens             Incl               ratio'

!loop through all the spectra in the xillver table   
  ! tot_spec_num = 1
  write(*,*) 
  write(*,*) 'I start reading the table...'
  do j = 1, tot_spec_num

     rownum = j
     ! rownum = 338382

     ! get the parameters (par_val) and the spectrum (spec) 
     call get_par_and_spec(unit, rownum, par_dim, par_val, spec_dim, spec)


     !Check if the spectrum has NaN
     do i = 1, spec_dim
        if(spec(i) .ne. spec(i)) then
           write(*,*) 'Found NaN in the spectrum: ', par_val(i)
           write(*,*) 'The value of the spectrum is: ', spec(i)
           write(*,*) 'Spectrum number on the table: ',  rownum
           cont_bad_spec = cont_bad_spec + 1
        end if
     enddo

     


     do i = 1, spec_dim                 
        dE = (ene_hi(i) - ene_lo(i))
        E  = (ene_hi(i) + ene_lo(i)) * 0.5
        write(10, *) E, spec(i) * E / dE
     enddo
     write(10,*) 'no no'

        !    cont_bad_spec = cont_bad_spec + 1 

           
        ! endif
        tot_spec_analysed = tot_spec_analysed + 1
     ! endif


     enddo
  write(*,*) "I've finished reading."
  write(*,*) 
     

666 continue
  
  write(10,*) 'scr white'
  write(10,*) 'log y x on'
  close(10)
  call close_fits(unit)

  write(*,*)  
  write(*,*)  
  write(*,*)  
  write(*,*) 'total number of spectra in the table', tot_spec_num
  write(*,*) 'total number of spectra analysed', tot_spec_analysed
  write (*,'(A, I6, A)') 'there are ', cont_bad_spec, ' bad spectra in the table'  

end program check_tables_NaN
