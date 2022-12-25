      SUBROUTINE BSYNBplatt(NALLIN)
*
*-----------------------------------------------------------------------
*
* BSYNB is merely a continuation of BSYN, data transfered
* via scratch-file 14
*
* B takes line-absorption coefficients generated by A and
* calculates the spectra.
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      use cubint_module

      include 'spectrum.inc'
*
      character*50 mcode
      real ETAD(ndp),XC(ndp),maxetad
      real fluxme(lpoint),icenter(lpoint),prof(lpoint)
      real mum,prf,iprf(lpoint)
      doubleprecision XL1,XL2,eqwidth
      doubleprecision XLBEG,XLEND,DEL
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),
     &  nnNTAU
*
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON/CSURF/ HSURF,Y1(NRAYS)
      COMMON/CANGLE/ NMY,XMY(nrays),XMY2(nrays),WMY(nrays)
      COMMON/TAUC/ TAU(ndp),DTAULN(ndp),NTAU
      COMMON/PIECES/ XL1,XL2,DEL,EPS,NMX,NLBLDU,IINT,XMYC,IWEAK
      COMMON/ROSSC/ ROSS(NDP),cross(ndp)
!
      COMMON/SPACE3/ SOURCEplatt(ndp),ERROR(ndp),
     &               tomatch(7*nrays*ndp+8*ndp+1+3*nrays)
!
*
* extension for large number of wavelengths and lines (monster II)
      doubleprecision xlambda
      doubleprecision source_function
      common/large/ xlambda(lpoint),source_function(ndp,lpoint),
     & maxlam,ABSO(NDP,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
*
*
      real fcfc(lpoint),y1cy1c(lpoint),y1y1(nrays,lpoint),
     & xlm(lpoint)

* mu-points for intensity output.
* Use 12 Gauss-Radau points
      logical extrap
      integer nangles
      parameter (nangles=12)
      integer iout(nangles)
      real muout(nangles),yout(nangles),isurf(nangles,lpoint),
     &     icsurf(nangles,lpoint),uin(nrays)
      data muout /0.010018, 0.052035, 0.124619, 0.222841, 0.340008,
     &            0.468138, 0.598497, 0.722203, 0.830825, 0.916958,
     &            0.974726, 1.000000/

* special version NLTE
      logical nlte
      common /nlte_common/ nlte

      logical findtau1,hydrovelo
      real velocity
      common/velo/velocity(ndp),hydrovelo
      logical debug
      data debug /.false./ 
*
* NLBLDU is a dummy
*
      DATA eqwidth/0./,profold/0./

      PI=3.141593
*
* Initiate angle quadrature points 
*
      NMY=NMX
      CALL GAUSI(NMY,0.,1.,WMY,XMY)
      DO 1 I=1,NMY
        XMY2(I)=XMY(I)*XMY(I)
    1 CONTINUE
*
* Initiate mode of calculation
* IINT =1  : intensity at MY=XMYC
* IINT =0  : flux
* XL1      : wavelength where synthetic spectrum starts
* XL2      : wavelength where synthetic spectrum stops
* DEL      : wavelength step
* IWEAK =1 : weak line approximation for L/KAPPA le. EPS ! not implemented anymore.
* note  XL2.gt.XL1
*
      iweak=0
      if (iint.gt.0) then

! if intensity flag, add a mu-point. It should preferably be mu=1.0 so that
! interpolation at the Gauss-Radau points is properly made at mu=1.0
! xmyc is read in input.f, just after the intensity flag
! so for safety we impose here xmyc=1.0

        if (abs(xmyc-1.0).gt.0.0001) then
          xmyc=1.0
          print*,'WARNING! the additional mu-point was set to 1.0!'
        endif
        nmy=nmy+1
        xmy(nmy)=xmyc
        wmy(nmy)=0.
        WRITE(6,200) XMYC,XL1,XL2,DEL
      end iF
      if(iint.eq.0) write(6,201) xl1,xl2,del
      if(iweak.gt.0) write(6,202) eps
*
* Continuum calculations:
*
*
* 13/03-2019 BPz. We now compute continuum at all wavelengths
* a little more costly in time, but avoid interpolation in wavelength
* continuum flux.
*

      do j=1,maxlam
        xlsingle=xlambda(j)
        do k=1,ntau
          x(k)=absocont(k,j)
          s(k)=absoscont(k,j)
!          s(k)=0.     !test
*
* NLTE case not implemented for continuum
*
          if (nlte) then
! test            bplan(k)=source_function(k,j)
            bplan(k)=bpl(T(k),xlsingle)
          else
            bplan(k)=bpl(T(k),xlsingle)
          endif
*
        enddo
       
!        do k=1,ntau
!         write(59,*) xlsingle,k,bplan(k)
!        enddo

        call traneqplatt(0)

!        do k=1,ntau
!          if (abs(source_function(k,j)/sourceplatt(k)-1.0).gt.
!     &         0.01) then
!            print*,'source differ',xlsingle,k,source_function(k,j),
!     &          sourceplatt(k)
!          endif
!        enddo
*
* calculate emergent continuum intensity at prescribed mu-points
*
        if (iint.gt.0) then 
          call cubintp14(y1,xmy,uin,muout,yout,iout,nmy,
     &                     nangles,3,0,0,0,
     &                     extrap)
          do k=1,nangles
            icsurf(k,j)=yout(k)
          enddo
        endif

! intensity at special angle is now disk center intensity (mu forced to 1.0)
        Y1CY1C(j)=Y1(NMY)

        FCFC(j)=4.*HSURF*pi
      enddo
*
* cont + line flux
*
      numb=0
      do j=1,maxlam
        xlsingle=xlambda(j)
!        if (abs(xlsingle-5001.0).lt.1.e-4) then
!          print*,'k','61','kappa line',abso(61,j)*ross(61),
!     &       absos(61,j)*ross(61)
!        endif

        do k=1,ntau
! the continuum opacity is not included in abso

          x(k)=abso(k,j)+absocont(k,j)
! test            x(k)=abso(k,j)
! test 
          s(k)=absos(k,j)+absoscont(k,j)
!          s(k)=0. !test
*
* NLTE case implemented for lines
*
          if (nlte) then
! compute total source function, for lines and continuum
            bplan(k)=(source_function(k,j)*abso(k,j)+
     &                bpl(T(k),xlsingle)*absocont(k,j))/x(k)
          else
            bplan(k)=bpl(T(k),xlsingle)
          endif
*
        enddo

        idebug=0
       
!        do k=1,ntau
!         write(60,*) xlsingle,k,bplan(k)
!        enddo
cc        if (abs(xlsingle-5349.40).lt.0.005) then
cc          print*,'bsynbplatt final check', xlsingle,
cc     &       source_function(10,j),abso(10,j),bpl(T(10),xlsingle),
cc     &       absocont(10,j),x(10),s(10)
cc        endif
!          do k=1,ntau
!            if (xlsingle.gt.5889.7.and.xlsingle.lt.5889.93) then
!              print555,xlsingle,k,
!     &          source_function(k,j),abso(k,j),
!     &          x(k),
!     &          bpl(T(k),xlsingle)
!555           format(f9.3,1x,i3,4(1x,1pe9.3))
!            endif
!          enddo

        call traneqplatt(idebug)

!        if (xlsingle.gt.5889.7.and.xlsingle.lt.5889.93) then
!        do k=1,ntau
!          if (abs(source_function(k,j)/sourceplatt(k)-1.0).gt.
!     &         0.01) then
!            print556,xlsingle,k,source_function(k,j),
!     &          sourceplatt(k)
!556         format('source differ',f9.3,1x,i3,1x,2(1x,1pe9.3))
!          endif
!        enddo
!        endif

        if (iint.gt.0) then
*
* calculate emergent line intensity at prescribed mu-points
*
          call cubintp14(y1,xmy,uin,muout,yout,iout,nmy,
     &                     nangles,3,0,0,0,
     &                     extrap)
          do k=1,nangles
            isurf(k,j)=yout(k)
          enddo
        endif

* starting with version 12.1, flux is not divided by pi anymore.
* F_lambda integrated over lambda is sigma.Teff^4
        prf=4.*pi*hsurf/fcfc(j)
        fluxme(j)=hsurf*4.*pi
        if (iint.gt.0) then
* save intensity at the prescribed angle (now mu=1.0)
          iprf(j)=y1(nmy)/y1cy1c(j)
          icenter(j)=y1(nmy)
        endif
        prof(j)=1.-prf
*
* find depth where tau_lambda=1
        if (hydrovelo) then
          findtau1=.true.
        else
          findtau1=.false.
        endif
        if (findtau1) then
          taulambda=tau(1)*(x(1)+s(1))
          do k=2,ntau
            taulambda=taulambda+(tau(k)-tau(k-1))*(x(k)+s(k)+
     &               x(k-1)+s(k-1))*0.5
            if (taulambda.ge.1.) then
              print333,xlambda(j),k,tau(k),taulambda,t(k),ro(k)
333           format(f10.3,x,i3,x,1pe11.4,x,1pe11.4,x,0pf7.1,x,
     &            1pe11.4,0p)

              goto 1966
            endif
          enddo
1966      continue
        endif
*
* save intensities
cc        do nlx=1,nmy
cc          y1y1(nlx,j)=y1(nlx)
cc        enddo
*
* End loop over wavelengths
*
      enddo
*
* Write spectrum on file for convolution with instrument profile
*
      if (iint.gt.0) then
        write(46,1111) muout(1:nangles)
1111    format ('# mu-points ',12(1x,1pe13.6))
      endif
      do j=1,maxlam
        plezflux=1.-prof(j)
        if (iint.eq.0) then
* fluxme is sigma teff^4
          write(46,1964) xlambda(j),plezflux,fluxme(j)
1964      format(f11.3,1x,f10.5,1x,1pe12.5)
        else
* We add intensities at prescribed angles after spectrum output,
* absolute and normalised
          write(46,1965) xlambda(j),plezflux,fluxme(j),
     &                   (isurf(k,j),isurf(k,j)/icsurf(k,j),k=1,nangles)
cc     &                   icenter(j),iprf(j)
1965      format(f11.3,1x,f10.5,1x,1pe12.5,12(1x,1pe12.5,1x,0pf8.5))
        endif
      enddo

      call clock
      return
*
  100 FORMAT(4X,I1,6X,I3)
  207 FORMAT(' SPECTRUM CALCULATED FOR THE FOLLOWING ',I3,' LINES'
     &      /' ELEMENT LAMBDA XI(EV) LOG(GF) LOG(ABUND) LINE NO')
  208 FORMAT('   ',A2,I2,F9.2,F5.2,1X,F6.2,3X,F7.2,4X,I5)
  200 FORMAT(' INTENSITY SPECTRUM AT MU=',F6.2,' BETWEEN',F10.3,
     &      ' AND',F10.3,' A WITH STEP',F6.3,' A')
  201 FORMAT(' FLUX SPECTRUM BETWEEN',F10.3,' A AND ',F10.3,' A WITH',
     &       ' STEP',F6.3,' A')
  202 FORMAT(' ***WEAK LINE APPROX FOR L/KAPPA L.T.',F7.4)
  203 FORMAT(' MODEL ATMOSPHERE:',A)
  204 FORMAT(' CONTINUUM FLUX=',E12.4, ' AT ',f10.2,'AA')
  205 FORMAT(' CONTINUUM INTENSITY=',E12.4,' at ',f10.2,'AA')
  206 FORMAT(1X,10F8.4)
 2341 FORMAT(1X,'Block number: ',I8)
*
      END
