program check_tables
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
  tot_spec_num = rowmax

!range between 20keV and 50keV (1713 - 1991)  
!range between 25keV and 35keV (1782 - 1883)  
  start_ind = 1713
  end_ind   = 1991
  save_par_val = 0.0
  cont_bad_spec = 0
  tot_spec_analysed = 0
  threshold = 1e5
  cont  = 0

  
  open (10, file='xill_spectra_ratio.qdp')
  write(10,*) 'skip on'
  open (20, file='xill_spectra_parameters_ratio.txt')
  write(20,*) '  Spec number      Gamma        A_Fe             log_xi            kTe              Dens             Incl               ratio'

!loop through all the spectra in the xillver table   
  ! tot_spec_num = 250000
  do j = 1, tot_spec_num

     rownum = j
     ! rownum = 338382

     ! get the parameters (par_val) and the spectrum (spec) 
     call get_par_and_spec(unit, rownum, par_dim, par_val, spec_dim, spec)

     ! if ((par_val(1) .gt. 1.4) .and. (par_val(4) .gt. 80) .and. (par_val(5) .gt. 20)) then
     ! write(*,*)  par_val(1), par_val(4)
     if ((par_val(1) .gt. 1.4) .and. (par_val(4) .gt. 80) ) then
        ! write(*,*) 'spectrum number ', j 

     
       !let's work on the predicted flux for the parameter that we extracted 
        call compton_hump_estimate(nmaxp, par_dim, par_val, wp, Icomp_scat)
        
        ! do i = 1, spec_dim
        !    dE = (ene_hi(i) - ene_lo(i))
        !    E  = (ene_hi(i) + ene_lo(i)) * 0.5
        !    spec_dble = dble(spec(i)) * dble(E/dE) * 1e20
        !    write(21, *) E, spec_dble   ! (E * 1e20 / dE) this is the factor to get the specifit intensity in keV
        ! enddo
        ! write(21,*) 'no no'

        ! do i = 1, nmaxp
        !    write(21,*) wp(i) / 1e3, Icomp_scat(i)
        ! enddo
        ! write(21,*) 'no no'

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

        ! if (i .eq. nmaxp) then
        !    do i = 1, nmaxp
        !       index_scat = i
        !       ! write(*,*) 'compare ', (wp(i) / 1e3), E
        !       if (abs((wp(i) / 1e3) - E) .lt. 0.5) then
        !          ! write(*,*) 'the corrisponding energy of the scattering spectrum is (index)', wp(i) / 1e3, index_scat
        !          exit
        !       endif
        !    enddo
        ! endif
           
        
        spec_dble = dble(spec(index_xil)) * dble(E/dE) * 1e20
        ! ratio = abs(dble(Icomp_scat(index_scat)) - spec_dble ) / dble(Icomp_scat(index_scat))
        ratio = dble(Icomp_scat(index_scat)) /  spec_dble
        ! write(12, *) j, ratio, Icomp_scat(index_scat), spec(index_xil) * E * 1e20 / dE
        
           ! write(*,*)
           ! write(*,*) 'value of the scattering spectrum and the xillver spectrum', Icomp_scat(index_scat), spec_dble, E, wp(index_scat) / 1e3, index_scat, index_xil
           ! write(*,*)
           ! write(*,*) 'relative ratio between the scattering spectrum and xillver at the compton hump energy', ratio

        if (ratio  .gt. threshold) then

           
           ! write(13,*)  rownum, E, wp(index_scat) / 1e3,  index_scat, index_xil, Icomp_scat(index_scat), spec(index_xil) * E * 1e20 / dE, ratio

           write(20,*)  rownum, par_val, ratio

           do i = 1, spec_dim                 
              dE = (ene_hi(i) - ene_lo(i))
              E  = (ene_hi(i) + ene_lo(i)) * 0.5
              write(10, *) E, spec(i) * E / dE
           enddo
           write(10,*) 'no no'

           cont_bad_spec = cont_bad_spec + 1 

           
        endif
        tot_spec_analysed = tot_spec_analysed + 1
     endif


  enddo

666 continue
  
  write(10,*) 'scr white'
  write(10,*) 'log y x on'
  close(10)
  close(20)
  call close_fits(unit)

  write(*,*)  
  write(*,*)  
  write(*,*)  
  write(*,*) 'total number of spectra in the table', tot_spec_num
  write(*,*) 'total number of spectra analysed', tot_spec_analysed
  write (*,'(A, I6, A)') 'there are ', cont_bad_spec, ' bad spectra in the table'  





end program check_tables





        ! write(*,*)
        ! write(*,*)

        
        ! write(20,*)  rownum, par_val

!renormalisation that happens inside xillver to match reflionx
      ! do j=2001,nmaxp   !just write from 0.1 keV - **
      !   coflux=flux(m,j)
      !   if(coflux.lt.0.d0)coflux=0.d0
      !   e1=wp(j-1)/1.d3
      !   e2=wp(j)/1.d3
      !   emis=coflux*(e2-e1)/wp(j)/1.d20*1.d3
      !   write(6,100)e1,e2,emis
      ! enddo
        
        ! if (tot_spec_analysed .gt. 1) then
        !    goto 666
        ! endif

        ! ! write(*,*) j, save_par_val
        ! ! write(*,*) j, par_val
        ! if ((par_val(1) .ne. save_par_val(1)) .or. (par_val(4) .ne. save_par_val(4))) then
        !    cont = cont + 1 
        !    spec_ref = spec
        !    sum_ref = 0.0
        !    do i = start_ind,  end_ind
        !       sum_ref = sum_ref + spec_ref(i)
        !    enddo
        !    sum_ref = sum_ref / real( end_ind - start_ind)
        !    save_par_val = par_val
           
        ! else
        !    !Bin 1836 corresponds to 30keV
        !    norm = spec_ref(1836) / spec(1836)

        !    sum = 0.0
        !    do i = start_ind,  end_ind
        !       sum = sum + (spec(i) * norm)
        !       ! variation = spec_ref(i) - (spec(i) * norm)
        !    enddo
        !    sum = sum / real( end_ind - start_ind)
        !    if (abs(sum_ref - sum)  .gt. threshold) then
        !       ! write(*,*) 'problem',  j, abs(sum_ref - sum)
        !       ! write(*,*) 'parameters',  par_val
        !       cont_bad_spec = cont_bad_spec + 1 
        !       write(20,*)  rownum, par_val
        !       do i = 1, spec_dim
                 
        !          write(10, *) (ene_hi(i) + ene_lo(i)) * 0.5, spec(i) / (ene_hi(i) - ene_lo(i))  
        !       enddo
        !       write(10,*) 'no no'
        !    endif
        ! endif




!Get the energies for the ENERGIES extension 
  ! call get_energy(unit, spec_dim, ene_lo, ene_hi)

  ! do i = 1, spec_dim
  !    ! write(*,*) ene_hi(i) , ene_lo(i)
  !    en(i) = (ene_hi(i) + ene_lo(i)) * 0.5
  !    ! write(*,*) en(i)
  ! enddo


  ! open (10, file='bad_spectra.qdp')
  ! write(10,*) 'skip on'

  ! do j = 1, tot_spec_num
  !    check = .FALSE.
 
  !    rownum = j
  !    call get_par_and_spec(unit, rownum, par_dim, par_val, spec_dim, spec)
  !    ! write(*,*) j, save_par_val
  !    ! write(*,*) j, par_val
  !    if ((par_val(1) .ne. save_par_val(1)) .or. (par_val(4) .ne. save_par_val(4))) then
  !       cont = cont + 1 
  !       spec_ref = spec
  !       sum_ref = 0.0
  !       do i = start_ind,  end_ind
  !          sum_ref = sum_ref + spec_ref(i)
  !       enddo
  !       save_par_val = par_val

  !    else
  ! !Bin 1836 corresponds to 30keV
  !       norm = spec_ref(1836) / spec(1836)

  !       sum = 0.0
  !       do i = start_ind,  end_ind
  !          sum = sum + (spec(i) * norm)
  !          ! variation = spec_ref(i) - (spec(i) * norm)
  !       enddo
  !       if (abs(sum_ref - sum)  .gt. threshold) then
  !          write(*,*) 'problem',  j, abs(sum_ref - sum)
  !          write(*,*) 'parameters',  par_val
  !          check  = .TRUE.
  !          cont_bad_spec = cont_bad_spec + 1 
  !       endif
  !       if (check) then 
  !          do i = 1, spec_dim
  !             write(10, *) i, spec(i)  
  !          enddo
  !          write(10,*) 'no no'
  !       endif
  !    endif
  ! enddo

  ! write(10,*) 'scr white'
  ! write(10,*) 'log y on'
  ! call close_fits(unit)

  ! close(10)
  ! write (*,'(A, I6, A)') 'there are ', cont_bad_spec, ' spectra in the table'  
  ! write(*,*) 'total numebr of spectra', tot_spec_num
  ! write(*,*) 'number of reference spectra', cont

!------------------------------------------------
  
  ! open (10, file='bad_spectra.qdp')
  ! write(10,*) 'skip on'
  ! do j = 1, tot_spec_num
  !    check = .FALSE.
 
  !    rownum = j
  !    call get_par_and_spec(unit, rownum, par_dim, par_val, spec_dim, spec)
  !    ! write(*,*) j, save_par_val
  !    ! write(*,*) j, par_val
  !    ! write(*,*) 'parameters',  par_val

  !    if (spec(spec_dim) .gt. 1e8) then
  !       ! write(*,*) spec(spec_dim)
  !       write(20,*) par_val
  !       cont_bad_spec = cont_bad_spec + 1 
  !       do i = 1, spec_dim
  !          write(10, *) (earx(i) + earx(i-1)) * 0.5, spec(i)  
  !       enddo
  !       write(10,*) 'no no'
  !    endif     
  ! enddo
  
  ! write(10,*) 'scr white'
  ! write(10,*) 'log y on'
  ! write(10,*) 'log x on'
  ! close(10)

  ! call close_fits(unit)
 
  ! write(*,*) 'total numebr of spectra', tot_spec_num
  ! write (*,'(A, I6, A)') 'there are ', cont_bad_spec, ' bad spectra in the table'  



!------------------------------------------------
