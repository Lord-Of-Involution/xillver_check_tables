program plot_single_spectrum
  implicit none

  integer, parameter :: spec_dim = 2999
  integer, parameter :: nmaxp = 5000
  double precision, parameter :: ergsev  = 1.602197d-12 ! Convert eV to ergs

  integer    :: i, unit, rownum, rowmax, start_ind, end_ind, index_scat, index_xil
  integer    :: par_dim
  real       :: threshold, E, dE, ratio, temp, spec(spec_dim)
  real       :: ene_lo(spec_dim), ene_hi(spec_dim)
  real       :: earx(0 : spec_dim)
  ! real       :: en(spec_dim)
  real, allocatable :: par_val(:), save_par_val(:)
  character (len=200) :: filename, hdu_name  
  integer           , allocatable :: num_of_val_par(:)
  character (len=12), allocatable :: par_names(:)

  double precision :: wp(nmaxp), spec_dble
  real             :: Icomp_scat(nmaxp)

  ! real, allocatable   :: spec_ref(:), spec(:)
  ! if (.not. allocated(spec_ref)) allocate(spec_ref(spec_dim))
  ! if (.not. allocated(spec)    ) allocate(spec    (spec_dim))


     filename = '/Users/gullo/Software/relxill/tables/xillverCp_v3.6.fits' !xillver-comp.fits

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

     !range between 20keV and 50keV (1713 - 1991)  
     !range between 25keV and 35keV (1782 - 1883)  
     start_ind = 1713
     end_ind   = 1991
     threshold = 1e2


     open (21, file='xillver_comp_scatt_spectrum.qdp')
     write(21,*) 'skip on'
     ! open (20, file='xill_spectra_parameters.txt')

     write(*,*) '   Insert the number of the spectrum'
     read(*,*) rownum
     ! rownum = 338382

     if (rownum .gt. rowmax) then
        write(*,'(A, I8, A)') '  there are only ', rowmax,' spectra in the table'
        stop
     else 
        ! get the parameters (par_val) and the spectrum (spec) 
        call get_par_and_spec(unit, rownum, par_dim, par_val, spec_dim, spec)
        write(*,*) '  Spec number      Gamma        A_Fe             log_xi            kTe              Dens             Incl               ratio'
        write(*,*) par_val

        !let's work on the predicted flux for the parameter that we extracted 
        call compton_hump_estimate(nmaxp, par_dim, par_val, wp, Icomp_scat)

        do i = 1, spec_dim - 1 
           dE = (ene_hi(i) - ene_lo(i))
           E  = (ene_hi(i) + ene_lo(i)) * 0.5
           spec_dble = dble(spec(i)) * dble(E/dE) 
           write(21, *) E, spec_dble * E   ! (E * 1e20 / dE) this is the factor to get the specifit intensity in keV
        enddo
        write(21,*) 'no no'

        do i = 1, nmaxp - 1 
          write(21,*) wp(i) / 1e3, Icomp_scat(i) * (wp(i) / 1e3) / 1e20
        enddo
        write(21,*) 'no no'

        !find energy index that peaks at the compton hump         
        temp = 0.0
        do i = start_ind,  end_ind
           index_xil = i
           ! E  = (ene_hi(index_xil) + ene_lo(index_xil)) * 0.5
           ! write(*,*) E,  spec(i)
           !this if statement checks when the spectrum starts to decrease 
           if (temp .gt. spec(i)) then
              exit
           endif
           temp = spec(i)
        enddo

        dE  = ene_hi(index_xil) - ene_lo(index_xil)
        E  = (ene_hi(index_xil) + ene_lo(index_xil)) * 0.5

        ! write(*,*) 'peak of the compton hump in keV (energy index and spectrum value)', E, index_xil, spec(index_xil)
        ! write(*,*)

        !find the compton hump energy at 
        do i = 1, nmaxp
           index_scat = i
           ! write(*,*) 'compare ', (wp(i) / 1e3), E
           if (abs((wp(i) / 1e3) - E) .lt. 0.1) then
              ! write(*,*) 'the corrisponding energy of the scattering spectrum is (index)', wp(i) / 1e3, index_scat
              exit
           endif
        enddo

        spec_dble = dble(spec(index_xil)) * dble(E/dE) * 1e20
        ratio = abs(dble(Icomp_scat(index_scat)) - spec_dble ) / dble(Icomp_scat(index_scat))
        ! write(12, *) j, ratio, Icomp_scat(index_scat), spec(index_xil) * E * 1e20 / dE

        write(21,*) 'scr white'
        write(21,*) 'log y x on'
        close(21)
        call close_fits(unit)

        write(*,*)
        write(*,*) 'Value of the scattering spectrum and the xillver spectrum', Icomp_scat(index_scat), spec_dble
        write(*,*) 'Value of the energies of the compton hump in the xillver spectrum and in the scattering spectrum ',&
             E, wp(index_scat) / 1e3
        write(*,*)
        write(*,*) 'relative ratio between the scattering spectrum and xillver at the compton hump energy', ratio

        write(*,*)  
        write(*,*)
     endif

   end program plot_single_spectrum
