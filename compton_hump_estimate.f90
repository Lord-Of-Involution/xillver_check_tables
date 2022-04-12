  subroutine compton_hump_estimate(nmaxp, par_dim, par_val, wp, Icomp_scat)
    implicit none

    integer, INTENT(IN)  :: par_dim, nmaxp
    real   , INTENT(IN)  :: par_val(par_dim)
    real*8 , INTENT(OUT) :: wp(nmaxp)
    real   , INTENT(OUT) :: Icomp_scat(nmaxp)

    double precision, parameter :: pi = acos(-1.0) 
    double precision, parameter :: ergsev  = 1.602197d-12 ! Convert eV to ergs
    double precision     :: gamma, xi, kTe, nh, kTbb, mu, f_x
    integer              :: np, ifl
    double precision     :: par(5), ear(0:nmaxp - 1), photar(0:nmaxp - 1), prim(0 : nmaxp - 1)
    double precision     :: fionnr, tau_max, Iinc, fioniz(nmaxp)
    double precision     :: emin, emax, emax2, df(nmaxp)
    
    
! assign parameters (ionisation and density)
    gamma = dble(par_val(1))
    xi    = dble(10**par_val(3))
    kTe   = dble(par_val(4))
    nh    = dble(10**par_val(5))
    mu    = cos(par_val(6))
    kTbb  = 0.05d0

    tau_max = 10.d0

    emin  = 1.d0
    emax  = 9.9e5
    emax2 = 1e6
    
!    Case a.i (Assuming an nthComp continuum)

    ! write(*,*) 'par_val', par_val
    ! write(*,*) 'parameters:  gamma, xi, nh, kTe', gamma, xi, nh, kTe
    f_x = xi * nh/ (4.d0 * pi)

       ! write(*,*) '  Bild internal energy grid'       
       call enegrd(nmaxp, emin, emax, emax2, wp, df)
       ! write(*,*) ' Done'
       ! write(*,*) 
       
       do np = 1, nmaxp
          ear(np - 1) = wp(np) / 1.d3
       enddo
       
       par(1) = gamma
       par(2) = kTe
       par(3) = kTbb
       par(4) = 1.d0
       par(5) = 0.d0
       Ifl=0


       ! write(*,*) '  Call nthcomp'       
       call donthcomp(ear, nmaxp - 1, par, ifl, photar, prim)
       ! write(*,*) '  Done'       
       ! write(*,*) 

       ! write(11,*) 'skip on'
       ! do np = 1, nmaxp - 1
       !    write(11,*) (ear(np) + ear(np - 1)) * 0.5, photar(np)
       ! enddo
       ! write(11,*) 'no no'
       
       fionnr = 0.d0
       do np = 1, nmaxp
          fioniz(np) = photar(np - 1)
          if (wp(np) .ge. 100.d0) then  ! Fx between 0.1 keV - 1 MeV
             fionnr = fionnr + (fioniz(np) + fioniz(np - 1)) &
                  * (wp(np) - wp(np - 1)) * ergsev/2.d0
             ! fionnr = fionnr + (fioniz(np) + fioniz(np - 1)) &
             !      * (wp(np) - wp(np - 1)) * 0.5d0
          endif
       enddo

       ! write(*,*) 'fionnr, f_x', fionnr, f_x

       do np = 1, nmaxp
          fioniz(np) = fioniz(np)/fionnr * f_x          
       enddo

       fionnr = 0.d0
       do np = 1, nmaxp
          if (wp(np) .ge. 100.d0) then  ! Fx between 0.1 keV - 1 MeV
             fionnr = fionnr + (fioniz(np) + fioniz(np - 1)) &
                  * (wp(np) - wp(np - 1)) * 0.5d0
          endif
       enddo

       ! write(*,*) 'Integration after the re-normalisation ', fionnr, f_x
       
       ! do np = 1, nmaxp
       !    write(11,*) (ear(np) + ear(np - 1)) * 0.5, fioniz(np)
       ! enddo
       ! write(11,*) 'no no'

       ! write(11, *) 'scri white'
       ! write(11, *) 'log x y on '
       
! let's calculate the expected flux
       do np = 1, nmaxp
          Iinc  = fioniz(np) * 2.d0 / cos(45.0)
          Icomp_scat(np) = real(Iinc * (tau_max / (tau_max + 2.d0 * mu)))
       enddo

       close(11)

     end subroutine compton_hump_estimate

!energy grid routine
  subroutine enegrd(nmaxp, ppemin, ppemax, ppemax2,&  !inp
                    wp, df) !out
    implicit none
    integer nmaxp
    real*8 wp(nmaxp), ppemin, ppemax, ppemax2
    integer numcon, numcon2, numcon3, ll, ll2
    real*8 ebnd1, ebnd2, ebnd2o, dele, df(nmaxp)

    numcon = nmaxp
    if (numcon.lt.4) stop 'in ener: numcon error'
    numcon2=max(2,nmaxp/50)
    numcon3=numcon-numcon2
    ebnd1=ppemin
    ebnd2=ppemax
    ebnd2o=ebnd2
    dele=(ebnd2/ebnd1)**(1./dfloat(numcon3-1))
    wp(1)=ebnd1

    do ll=2,numcon3
       wp(ll)=wp(ll-1)*dele
    enddo

    ebnd2=ppemax2
    ebnd1=ebnd2o
    dele=(ebnd2/ebnd1)**(1./dfloat(numcon2-1))

    do ll2=1,numcon2
       ll=ll2+numcon3
       wp(ll)=wp(ll-1)*dele
    enddo

    df(1)=0.d0

    do ll=2,nmaxp
       df(ll) = wp(ll) - wp(ll-1)
    enddo

    return

  end subroutine enegrd
