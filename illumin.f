      subroutine illumin(nmaxp, itrans, wp, gioniz, ecut, ffrac,  !inp
     1                   ikbol, ergsev, tau, pi, nh, xi, kTbb,    !inp
     1                   f_x, fdisk, bbody, fioniz, aph, Teff)    !out
      implicit none
      integer nmaxp, itrans, np, iz, nb
      real*8 wp(nmaxp), gioniz, ecut, ikbol, ergsev, ffrac
      real*8 tau(itrans), f_x, fdisk, kTbb, ecut1, pi, nh, xi
      real*8 bbody(nmaxp), fioniz(nmaxp), aph(nmaxp,itrans)
      real*8 bbnorm, fionnr, pht, phtmax, amax, k1, ebr, Teff
      real*8 par(5), ear(0:nmaxp-1),Photar(0:nmaxp-1), prim(0:nmaxp-1)
      integer Ifl
c
      ecut1 = ecut*1.d3
      f_x = 0.d0
      fdisk = 0.d0
      bbnorm = 0.d0
      fionnr = 0.d0
      Teff = kTbb*1.d3*ikbol
c
      do np = 1, nmaxp
         fioniz(np) = 0.d0
         bbody(np) = 0.d0
      enddo
c
c     Possibilities for the illumination:
c
c       a) Only illumination from above, single:  i) PL;        (ffrac=0, gioniz > 0)
c                                                ii) BB;        (ffrac=0, gioniz <= 0)
c       b) Only illumination from above, dual:  i) BB+PL        (ffrac<0, |ffrac|=PL/BB)
c       c) Illumination from above and below: i) PL+BB;         (ffrac>0, ffrac=PL/BB)
c                                            ii) BB+BB          (ffrac>0, not sure yet!!)
c
c
c
c      Case a.i (Assuming an nthComp continuum)
       if (ffrac.eq.0.and.gioniz.gt.0) then
          print *, 'Running Case a.i: PL illumination at the top (nth)'
c
          f_x = xi*nh/(4.d0*pi)
c
          do np=1,nmaxp
             ear(np-1)=wp(np)/1.d3
          enddo
          par(1)=gioniz
          par(2)=ecut
          par(3)=kTbb
          par(4)=1.d0
          par(5)=0.d0
          Ifl=0
c
          call donthcomp(ear,nmaxp-1,par,Ifl,Photar,prim)
c
          do np=1,nmaxp
             fioniz(np)=Photar(np-1)
             if (wp(np).ge.100.d0) then  ! Fx between 0.1 keV - 1 MeV
                fionnr = fionnr + (fioniz(np)+fioniz(np-1))
     1                   * (wp(np)-wp(np-1))*ergsev/2.d0
             endif
          enddo
c
          do np = 1, nmaxp
             fioniz(np) = fioniz(np)/fionnr*f_x
             write(70,*)wp(np),fioniz(np)
          enddo
c
c      Case a.ii (assuming a Blackbody illumination from above)
       else if (ffrac.eq.0.and.gioniz.le.0) then
c
          fdisk = xi*nh/(4.d0*pi)
          f_x = fdisk
c
          print *, 'Running Case a.ii: BB illumination at the top'
          print *, 'Fdisk=',real(fdisk)
          print *, ' '
c
          pht = kTbb*1.d3
          do np = 1, nmaxp
             fioniz(np) = (wp(np)**3.) * dexp(-wp(np)/pht)/
     1                   (1.d0 - dexp(-wp(np)/pht))
             if (wp(np).ge.100.d0) then  ! Fx between 0.1 keV - 1 MeV
                fionnr = fionnr + (fioniz(np)+fioniz(np-1))
     1                   * (wp(np)-wp(np-1))*ergsev/2.d0
             endif
          enddo
c
          do np = 1, nmaxp
             fioniz(np) = fioniz(np)/fionnr*f_x
             if (fioniz(np).lt.1.d-30)fioniz(np)=1.d-30
!             write(70,*)wp(np),fioniz(np)
          enddo
c
c
c      Case b.i 
       else if (ffrac.lt.0) then
          stop 'Case not implemented yet!'
c      Case c.i 
       else if (ffrac.gt.0) then                    ! c.i
          print *, 'Illumination Case C.i'
!          f_x = xi*nh/(4.d0*pi)
!          fdisk = f_x/ffrac
          fdisk = 5.67d-5*Teff**4.
          f_x = ffrac*fdisk
          xi = 4.d0*pi*f_x/nh
c
c         Right now does PL here... need to include switch for nthComp case!!!

          pht = Teff/ikbol
          do np = 1, nmaxp
            if(wp(np).ge.pht*2.821439d0)then
               nb = np
               goto 10
            endif
          enddo
10        continue
c
          k1=2.43d-4*ikbol
          phtmax = wp(nb)/k1
          amax= dexp(-wp(nb)/ecut1)/(wp(nb)**(gioniz-1.d0))
     1         *(dexp(k1)-1.d0)/(phtmax*k1)**3.
c
          do np = 1, nmaxp
              bbody(np) = 1.d-34 + (wp(np)**3.) * dexp(-wp(np)/pht)/
     1                    (1.d0 - dexp(-wp(np)/pht))
c
!             if(np.lt.nb)then
!               fioniz(np) = amax*(wp(np)**3.) * dexp(-wp(np)/phtmax)/
!     1                      (1.d0 - dexp(-wp(np)/phtmax))
!             else
               fioniz(np) = (dexp(-wp(np)/ecut1))/
     1                      (wp(np)**(gioniz-1.d0))
!             endif
             if (np.gt.1) then
                bbnorm = bbnorm + (bbody(np)+bbody(np-1))
     1                   * (wp(np)-wp(np-1))*ergsev/2.d0
!               if (wp(np).ge.100.d0) then  ! Fx between 0.1 keV - 1 MeV
                  fionnr = fionnr + (fioniz(np)+fioniz(np-1))
     1                     * (wp(np)-wp(np-1))*ergsev/2.d0
!               endif
             endif
          enddo
          print *, '-------fionnr',fionnr,' fx',f_x
          do np = 1, nmaxp
             bbody(np) = bbody(np)/bbnorm*fdisk
             fioniz(np) = fioniz(np)/fionnr*f_x
          enddo
c
      endif

c     Initial photon distribution
      do iz = 1, itrans
         do np = 1, nmaxp
            if (ffrac.gt.0) then                    ! c.i
               aph(np,iz) = fioniz(np) + bbody(np)
            else
               aph(np,iz) = fioniz(np)
            endif
         enddo
      enddo
c
      return
      end subroutine
