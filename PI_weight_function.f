      program calc_PI_weight_function

      implicit none

      integer pbead,PmP            ! P: number of replacs or beads
      integer niter,lw,lpp,nout,npo
      integer i,j,k,ir,iw

      REAL*16 dw,doutw,gx,x,dev,pmax,irmin,prcf,stout
      REAL*16 PI

      REAL*16 hx,hgx,htgx,suro,cro

      REAL*16, allocatable :: w(:),UA(:,:),OPM(:,:),ro(:)
      REAL*16, allocatable :: rotest(:),rto(:),irmx(:)
      REAL*16, allocatable :: fri(:)
      INTEGER, allocatable :: ri(:),shift(:),pox(:)

      PI=4.e0*ATAN(1.e0)

      write (*,821)
  821 format (/,' Number of beads:',
     &             '',$)
      read (*,*) pbead
      write (*,*) pbead

      niter=100

      dw=1.0e-2

      if (pbead.eq.1) then

          niter=2

          dw=1.0e0

      endif

      if (pbead.eq.2) then

          dw=4.0e0

      endif

      if (pbead.eq.3) then

          dw=6.75e0

      endif

      if (pbead.eq.4) then
 
          dw=8.0e0

      endif

      if (pbead.eq.6) then

          dw=9.0e0

      endif

      write (*,822)
  822 format (/,' Start writting of g(x) at:',
     &             '',$)
      read (*,*) stout
      write (*,*) stout

      write (*,823)
  823 format (/,' Interval of x-s:',
     &             '',$)
      read (*,*) doutw
      write (*,*) doutw

      write (*,824)
  824 format (/,' Number of datapoints:',
     &             '',$)
      read (*,*) nout
      write (*,*) nout

      prcf=1.0e-16

      write(*,*) "pbead,niter"

      write(*,*) pbead,niter


      PmP=pbead*pbead

      lpp=nint(PmP/dw)

      lw=nint(PmP*niter/dw)

      write(*,*) "lw,lpp"

      write(*,*) lw,lpp

      allocate (w(lw))

      allocate (ri(pbead))

      allocate (fri(pbead))

      allocate (shift(pbead))

      allocate (pox(pbead))

      allocate (ro(lw))

      allocate (rto(lw))

      allocate (irmx(lw))

      w=0.0e0

      pox=0.0e0

      npo=0

      pmax=0.0e0

      do k=0,pbead-1

        fri(k+1)=PmP*sin(k*PI/pbead)**2

        ri(k+1)=nint(sin(k*PI/pbead)**2*lpp+1)

      enddo

      do k=1,pbead

        shift(k)=ri(k)-1

      enddo

      i=1

      do k=1,pbead

        w(shift(k)+1)=w(shift(k)+1)+1.0e0/dble(pbead)

        if (pox(i).eq.0.0e0.and.shift(k)+1.gt.pmax) then

          pmax=shift(k)+1

          pox(i)=pmax

          i=i+1

        endif

      enddo

      npo=i-1

      write(*,*) "npo"

      write(*,*) npo

      write(*,*) "pox"

      write(*,*) pox


      w(1)=w(1)-1.0e0

      ro=0.0e0

      ro(1)=1.0e0

      write(*,*) "start iteration"

      do i=1,niter

        rto=0.0e0

        do j=1,lw

          cro=0.0e0

          k=1

          iw=pox(1)

          ir=iw+j-1

          do while (k.le.npo.and.ir.le.lw)

            rto(ir)=rto(ir)-ro(j)*w(iw)

            k=k+1

            iw=pox(k)

            ir=iw+j-1

          enddo


        enddo

        ro=rto

        ro(1)=1.0e0

        write(*,*) i

        suro=sum(ro)

        write(*,*) suro

      enddo

      irmin=ro(lw)

      do i=1,lw-1

        j=lw-i

        irmin=max(irmin,ro(j))

        irmx(j)=irmin/prcf

      enddo

      write(*,*) "omega,wP"

      do i=0,nout-1

        x=real(i,16)*doutw+stout

        write(*,*) x,gx(x,ro,irmx,lw,dw)

      enddo


      write(*,*) "omega,Cvv,Ta,error"

      do i=0,nout-1

        x=real(i,16)*doutw+stout

        hx=x/tanh(x)

        htgx=hgx(x,ro,irmx,lw,dw,fri,pbead)

        dev=htgx-hx

        write(*,*) x,hx,htgx,dev

      enddo
  




      deallocate (w)

      deallocate (ri)

      deallocate (fri)

      deallocate (shift)

      deallocate (pox)

      deallocate (ro)

      deallocate (rto)

      deallocate (irmx)

      return

      end

      FUNCTION gx(x,ro,irmx,lw,dw)
      INTEGER lw
      REAL*16 gx,dw,ro(lw),irmx(lw)
      INTEGER i
      REAL*16 x,y,xsx,sqy,hfy

      xsx=x*x

      y=xsx

      gx=0.0e0

      i=1

      do while (hfy.le.irmx(i))           

        sqy=sqrt(y)

        hfy=tanh(sqy)*sqy

        gx=gx+ro(i)/hfy

        y=y+dw

        i=i+1

      enddo 

      gx=gx*xsx

      RETURN
      END FUNCTION

      FUNCTION hgx(x,ro,irmx,lw,dw,fri,pbead)
      INTEGER lw,pbead
      REAL*16 hgx,gx,dw,ro(lw),irmx(lw),fri(pbead)
      INTEGER i
      REAL*16 x,y,xsx,sqy,xk,xsk

      xsx=x*x

      y=xsx

      hgx=0.0e0

      do k=1,pbead

        xsk=xsx+fri(k)

        xk=sqrt(xsk)

        hgx=hgx+gx(xk,ro,irmx,lw,dw)/xsk

      enddo

      hgx=hgx*xsx

      RETURN
      END FUNCTION

