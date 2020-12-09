!-----------------------------------------------------------------
      subroutine read_rownum(unit, rowmax)
      implicit none
      integer           :: unit, rowmax

      logical           :: anynul
      character (len=30):: hdu_name 
      integer           :: status,  chdu,  hdutype, colnum, felem
      integer           :: nrow
! Initialize status
      status    = 0
!Let's move in the correct HDU (with the name in hdu_name)       
      hdu_name = 'SPECTRA'
      call hdu_move(unit, hdu_name, status)
      if (status .gt. 0) call printerror(status)
      call ftghdn(unit, chdu)
      if (status .gt. 0) call printerror(status)
      
! Get the number of rows and columns in the CHDU
      call ftgnrw(unit, nrow, status)
      if (status .gt. 0) call printerror(status)
      rowmax = nrow
      ! call ftgncl(unit,ncol,status)

      ! ! Go to the extension 4
      ! chdu = 4
      ! call FTMAHD(unit,chdu,hdutype,status)


      return
    end subroutine read_rownum
!-----------------------------------------------------------------
      
!-----------------------------------------------------------------
      subroutine read_spec(unit, rownum, spec_dim, spec)
      implicit none
      integer           :: unit, rownum, spec_dim
      double precision  :: spec(spec_dim)

      logical           :: anynul
      character (len=30):: hdu_name 
      integer :: status,  chdu, hdutype, colnum, felem
      
! Initialize status
      status    = 0
!Let's move in the correct HDU (with the name in hdu_name)       
      hdu_name = 'SPECTRA'
      call hdu_move(unit,hdu_name,status)
      call ftghdn(unit,chdu)
      
! Read in the function 
      felem  = 1
      colnum = 2
      call ftgcve(unit, colnum, rownum, felem, spec_dim, 0.0,&
          spec, anynul, status)
      if (status .gt. 0) call printerror(status)

      return
    end subroutine read_spec
!-----------------------------------------------------------------
      
      
!------------------------------------------------------------------
    subroutine open_fits(filename, unit)
      implicit none
      integer :: status, unit
      character* (*) filename

      integer :: readwrite, blocksize
      
! Initialize status
      status    = 0
      readwrite = 0
      blocksize = 1
! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
! Open the FITS file
      call ftopen(unit,filename,readwrite,blocksize,status)

!! The FITS file must always be closed before exiting the program. 
!! Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
!      call ftclos(unit, status)
!      call ftfiou(unit, status)
! Check for any error, and if so print out error messages.
      if (status .gt. 0) call printerror(status)

      return
    end subroutine open_fits
!------------------------------------------------------------------
                                          
!------------------------------------------------------------------
    subroutine close_fits(unit)
      implicit none
      integer :: unit, status

! Initialize status
      status    = 0
! The FITS file must always be closed before exiting the program. 
! Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
! Check for any error, and if so print out error messages.
      if (status .gt. 0) call printerror(status)

      return
    end subroutine close_fits
!------------------------------------------------------------------

!---------------------------------------------------------------------!
    subroutine hdu_move(unit,hdu_name,status)
!PURPOSE:  Move to the HDU based on the name. 
      implicit none
      integer :: unit,hdutype,hdu_extver,status
      character (len=30) hdu_name,error_description

! The hdutype parameter may have a value of IMAGE HDU(0), ASCII TBL (1),
! BINARY TBL (2), or ANY HDU (-1) where ANY HDU means that
! only the extname and extver values will be used to locate the correct extension.
! If the input value of extver is 0 then the EXTVER keyword is ignored and the first HDU with a matching
! EXTNAME (or HDUNAME) keyword will be found. If no matching HDU is found in the file
! then the current HDU will remain unchanged and a status = BAD HDU NUM (301) will be returned.     
      hdutype = -1
      hdu_extver = 0
      call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      call ftgerr(status,error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
      
! if the extension is not RATE the code asks for a new name until it gets a right one      
      do while(status .eq. 301)
         write(*,*) "   No extension called: ",hdu_name
         write(*,*) "   Please specify the correct name or type 'no' to quit"
         read(*,*) hdu_name
         if (hdu_name .eq. 'no') stop
         hdutype = -1
         hdu_extver = 0
         status = 0
         call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      enddo
     
    end subroutine  hdu_move
!---------------------------------------------------------------------!


    
!------------------------------------------------------------------
    subroutine printerror(status)
!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      integer status
      character errtext*30,errmessage*80
!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return
!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext
!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
    end subroutine printerror
!----------------------------------------------------------------



