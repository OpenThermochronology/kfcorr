!	Program kfcorr (version 2.10)
!                 
!   Version 2.0 completed 10 July, 2012
!   Version 2.1 September, 2023

!	Version of Oscar Lovera's corrfft code for calculating the
!	correlation coefficient between age spectra and logR/Ro plots
    
! 	Changes made by Peter Zeitler in 2009 and (mostly) 2012
!	are largely about simpler methods of input and output,
!	and also addition of a plotting option if the user has the gmt package
!	installed.

!	USAGE:
!
!	    ./kfcorr  spectrumfilename  RRofilename  plotflag  samplename
!
!   The first file is floss-and-age and the second is floss-and-log10(R/Ro).
!	Both files must have the same floss values and number of steps; last floss should be 1.0
!	The file needs to have UNIX line endings.
!
!	plotflag: 0 - no plotting   1 - user has gmt installed and wants a plot
!
!	The samplename is used to tag the output filename. This can be omitted for a generic result.
!
!   You can place a single copy of the executable in any directory in your PATH (e.g. /usr/local/bin)
!   and the code will just run from your working directory like this (no'./'):
!
!   kfcorr  spectrumfilename  RRofilename  plotflag  samplename
!
!	COMPILATION:
!
!	There are issues with compiling the original legacy program in that gfortran will not compile
!	this code to produce a reliable working executable unless the relevant compiler flags are
!	used to allow flexibility in the use of variables (modern f90 implementations are more strict).
!	This code will now compile under freeform mode in gfortran, as follows:
!
!       gfortran kfcorr210.f90 -o kfcorrm2 -fno-automatic -O2  -fallow-argument-mismatch -w
!
!   4 January 2013: small update to permit executable to run from any working directory
!   provided directory that holds executable is in user's PATH. Also, moved output plot down
!   just a bit so plot doesn't crowd top margin.
!
!	September, 2023. Version 2.10 based on version 2.01: updated gmt plot calls to work with gmt 5 and gmt 6.
!

	parameter(ns=128,nc=100,ns2=256)
	dimension cros(ns2),cor1(ns2),cor2(ns2)
	dimension yy(ns)
	dimension f39(0:ns),age(0:ns),xlogr(0:ns),cint(nc),f39in(0:ns),agein(0:ns),xlogrin(0:ns)
	dimension agem(ns,0:ns),xlgrm(ns,0:ns)
	dimension c1(nc),c2(nc),cc(nc,nc)
	external favr
	
	integer gmt, inputs, agerange, logrange
	integer nspec,j
	
	double precision dummy, yagemax, yagemin, ylogmin

	character rfragment1*50, rfragment2*50, yfragment1*50, yfragment2*50, s1*8,s2*8
	character  tab1*9, filename*50, systemstring*250, stringfrag*200, morestring*200
	character samplename*20, command*50

    character(len=255) :: path
	
! ********************************************

! Before anything else, check that command line entries hold any hope of success

	inputs = command_argument_count()

	if (inputs.lt.3) then
		print *, ' '
		print *, ' kfcorr 2.10:'
		print *, ' USAGE:  kfcorr  spectrumfilename  RRofilename  plotflag  samplename'
		print *, ''
		print *, ' Input files are white-space delimited table of floss-and-age and floss-and-log10(R/Ro).'
		print *, ''
		stop
	end if
	
	if (inputs.gt.4) then
		print *, ' WARNING: extra command-line inputs. Program will continue, but be careful.'
		print *, ''
	end if	
	
	if (inputs.eq.3) then
		print *, ' NOTE. No samplename - output file will be kfcorr-generic.ps .'
		print *, ''
	end if

! Demonstrate program is running by greeting user

	write(*,*) ' '	
	write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'	    
	write(*,*) '                 PROGRAM kfcorr version 2.10'
	write(*,*) ' '
	write(*,*) '     Cross-correlation of 40-39 age spectra and logR/Ro data'
	write(*,*) ' '
	write(*,*) '           Original code corrfft by Oscar Lovera'
	write(*,*) '                 modified by Peter Zeitler'
	write(*,*) ' '    
	write(*,*) '                 Last updated September, 2023'
	write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	write(*,*) ' '

! make sure we are operating in user's working directory

	call getcwd(path)
	call chdir(TRIM(path))
	
	tab1=char(09)
	gmt = 0   ! give this a value; real value should supercede this from command line	

! Get range for modeling

	write (*,'("Enter fractional-loss range for modeling: ")')

	fmin = -100.
	do while ((fmin.lt.0.).or.(fmin.gt.1.))
		write (*,'("Lower limit: ")',advance='no')
		read *,fmin
	end do
	
	fmax = -100.
	do while ((fmax.le.fmin).or.(fmax.gt.1.))
		write (*,'("Upper limit: ")',advance='no')
		read *,fmax
	end do
	write (*,*) ' '	

! Parse the command-line entries (no mercy: they have to be in the correct order!)

	call get_command_argument(1, command)
	read(command,'(A50)') filename
	open(unit=22,file=TRIM(filename),status='old')
	
	call get_command_argument(2, command)
	read(command,'(A50)') filename
	open(unit=23,file=TRIM(filename),status='old')
	
	call get_command_argument(3, command)
	read(command,'(I1)') gmt
		
	if (inputs.gt.3) then
		call get_command_argument(4, command)
		read(command,'(A20)') samplename
	else
		samplename = 'generic'
	end if

! Read the input files
	
	do j=1,500
		read(22,*,END=20) f39(j),age(j)
		read(23,*,END=20) dummy,xlogr(j)
	end do
20	ni=j	! note that as coded we get one extra step on input because test for EOF is inside loop
	close(22)
	close(23)	

! ****** Here we need to recast the input data into the required forms
!
!        We want to plot the data as spectra,so we need arrays to hold that.
!        We have to feed Oscar's code the midpoint-smoothed data.
!        
!   When we're done f39(), age(), and xlogr() need to consistent with midpoint format, and
!	f39in(), agein(), and xlogrin() should be consistent with spectrum format for plotting.
!

! first deal with making spectra

	nspec = 2*(ni-1)
	
	f39in(1) = 0.0
	do j = 2,nspec,2
		f39in(j-1) = f39(j/2 - 1)
		f39in(j) = f39(j/2)
		agein(j-1) = age(j/2)
		agein(j) = age(j/2)
		xlogrin(j-1) = xlogr(j/2)
		xlogrin(j) = xlogr(j/2)
!			print *,j,f39(j),xlogr(j)
	end do	

! then deal with making midpoint data
	f39(ni) = 1.d0
	age(ni) = age(ni-1)
	xlogr(ni) = xlogr(ni-1)	

	do j = 2,nspec,2
		f39(j/2) = 0.5*(f39in(j) + f39in(j-1))
	end do

! ****** End of data diddling

	x0 =  0.
	mc = 15
	n = 1
   
!   MAKE MATRIX

      
	xl = fmax - fmin
	
	call chebft(x0,xl,c1,mc,favr,f39,age,ni,fmin)
	
	call chebft(x0,xl,c2,mc,favr,f39,xlogr,ni,fmin)
	
	call chint(x0,xl,c1,cint,mc)
	
	aget = chebev(x0,xl,cint,mc,xl)/xl
	
	call chint(x0,xl,c2,cint,mc)
	
	xlogt = chebev(x0,xl,cint,mc,xl)/xl     
	dx=xl/(ns-1)
	
	do j = 1,ns
		yy(j) = dx*(j-1)
		if (yy(j).gt.xl) yy(j) = xl
		
		agem(n,j) =  chebev(x0,xl,c1,mc,yy(j))-aget
		
		xlgrm(n,j) = chebev(x0,xl,c2,mc,yy(j))-xlogt
		
	end do
	agem(n,0) = fmin
	xlgrm(n,0) = fmax

!  CALCULATION OF CROSS AND PERMUTATIONS

	k1 = 1
	do j = 1,ns
		age(j) = agem(k1,j)
		xlogr(j) = xlgrm(k1,j)
	end do
  
	call correl(age,xlogr,ns,cros)
	
	call correl(age,age,ns,cor1)
	
	call correl(xlogr,xlogr,ns,cor2)
	
	cornor = sqrt(cor1(1)*cor2(1))
	cc(k1,k1) = cros(1)/cornor
	fmin = agem(k1,0)
	fmax = xlgrm(k1,0)
	
	write(*,*) ' '
	write(*, '("RESULTS FROM SAMPLE ",A20)') samplename
	write(*,'("---------------------------------------")')
	write(*,'(A,f7.4)')'             fmin: ',fmin          
	write(*,'(A,f7.4)')'             fmax: ',fmax    	
	write(*,'(A,f7.4)')'Cross-correlation: ',cc(k1,k1)
	write(*,'("---------------------------------------")')
	write(*,*) ' '
	
!	do j = -ns/2,-1
!		nj = ns+j+1
!		write(28,*)dx*(j),cros(nj)/cornor,cor1(nj)/cor1(1),cor2(nj)/cor2(1)
!	end do
	
!	do j = 1,ns/2
!		nj = j
!	      write(28,*)dx*(j-1),cros(nj)/cornor,cor1(nj)/cor1(1),cor2(nj)/cor2(1)
!	end do

! ***********************  plotting routines using gmt (if installed and user selected it) **********

	if (gmt.eq.1) then
	! set some gmt parameters
		call system("gmt gmtset LABEL_FONT_SIZE 14p") 
		call system("gmt gmtset ANNOT_FONT_SIZE_PRIMARY 12p")
		call system("gmt gmtset MAP_GRID_PEN_PRIMARY 0.25p,100/100/100,3_3:0p")
		call system("gmt gmtset CHAR_ENCODING ISOLatin1+")
		call system("gmt gmtset TICK_LENGTH 0.1i")
		call system("gmt gmtset TICK_PEN 0.75p")
		call system("gmt gmtset FRAME_PEN 1p")	
	
	! Get some values for Y-axis ranges
		yagemax = 0.
		yagemin = 5000.
		ylogmax = -10.
		ylogmin = -0.5  ! we'll just keep this fixed
		do  j=1,nspec
			if (agein(j).gt.yagemax.and.f39in(j).gt.0.20) yagemax = agein(j)
			if (xlogrin(j).gt.ylogmax) ylogmax = xlogrin(j)
			if (agein(j).lt.yagemin) yagemin = agein(j)
		end do

	! now deal with those ranges
	
		agerange = INT(yagemax - yagemin)
	
		select case (agerange)
			case ( : 5 )
				yfragment1 = "a1f0.5g1:'Age (Ma)':WS"
			case ( 6 : 10 )
				yfragment1 = "a2f1g2:'Age (Ma)':WS"
			case ( 11 : 50 )
				yfragment1 = "a5f2.5g5:'Age (Ma)':WS"
			case ( 51 : 100 )
				yfragment1 = "a10f5g10:'Age (Ma)':WS"
			case ( 101 : 250 )
				yfragment1 = "a20f10g20:'Age (Ma)':WS"
			case ( 251 : 500 )
				yfragment1 = "a50f25g50:'Age (Ma)':WS"
			case ( 501 : )
				yfragment1 = "a100f50g100:'Age (Ma)':WS"
		end select
	
		logrange = INT(ylogmax - ylogmin)
		
		select case (logrange)
			case ( : 1 )
				yfragment2 = "a0.5f0.1:'log10(R/R@-o@-)':E"
			case ( 2 )
				yfragment2 = "a0.5f0.25:'log10(R/R@-o@-)':E"
			case ( 3 )
				yfragment2 = "a0.5f0.5:'log10(R/R@-o@-)':E"
			case ( 4 :  )
				yfragment2 = "a0.5f0.5:'log10(R/R@-o@-)':E"
		end select

! Build -R option for age axis and plotting		
		y1min = 10.*FLOOR(yagemin/10.)
		write(s1,'(F6.0)') y1min
		y1max = 10.*CEILING(yagemax/10.)
		write(s2,'(F6.0)') y1max	
		rfragment1 = '-R/0/1.0/'//TRIM(ADJUSTL(s1))//'/'//TRIM(ADJUSTL(s2))
	
! Build -R option for log10RRo axis and plotting
		y2min = 0.5*CEILING(ylogmin/0.5)
		write(s1,'(F6.1)') y2min	
		y2max = 0.5*CEILING(ylogmax/0.5)

		write(s2,'(F6.1)') y2max
		rfragment2 = '-R/0/1.0/-0.5/'//TRIM(ADJUSTL(s2))

! Plot frame, plus left (age) axis		
		stringfrag = "gmt psbasemap -Xc -Y6.5i "//TRIM(rfragment1)//" -Bf0.05a0.2g0.1:'Fractional @+39@+Ar Loss':/"//TRIM(yfragment1)//&
			" -JX3.5i/2.75i -K -P > "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))

! plot right (log) axis	
		stringfrag = "gmt psbasemap "//TRIM(rfragment2)//" -B:/:/"//TRIM(yfragment2)//&
			" -JX3.5i/2.75i -O -K -P >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))
		
		stringfrag = "gmt psbasemap "//TRIM(rfragment2)//" -Bf0::/::n  -JX3.5i/2.75i -O -K -P >> "//&
			TRIM(samplename)//"-correlation.ps"	
		call system(TRIM(stringfrag))	

! plot age spectrum as observed
		open(unit=99,file='line.scr',status='UNKNOWN')
		do j = 1,nspec
			write (99,*) f39in(j),agein(j)
		end do
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i "//TRIM(rfragment1)//" -W3p,200/200/220 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))	

! plot logRRo spectrum as observed
		open(unit=99,file='line.scr',status='UNKNOWN')
		do j = 1,nspec
			write (99,*) f39in(j),xlogrin(j)
		end do
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i "//TRIM(rfragment2)//" -W3p,220/200/200 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))	

! plot age spectrum	over range modeled
		open(unit=99,file='line.scr',status='UNKNOWN')
		do j = 1,nspec
			if (f39in(j).ge.fmin.and.f39in(j).le.fmax) then
				write (99,*) f39in(j),agein(j)
			end if
		end do
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i "//TRIM(rfragment1)//" -W3p,100/100/255 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))	

! plot logRRo spectrum over range modeled
		open(unit=99,file='line.scr',status='UNKNOWN')
		do j = 1,nspec
			if (f39in(j).ge.fmin.and.f39in(j).le.fmax) then
				write (99,*) f39in(j),xlogrin(j)
			end if
		end do
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i "//TRIM(rfragment2)//" -W3p,255/100/100 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))	

! Some reporting text and ornaments - for the following calls to gmt we use a fictitious set of values for -R, spanning
!     y between 0 and 10. We need to do this because both axes will change scaling depending on the input data, so
!     the positioning of text labels would change from plot to plot if we used data coordinates from either axis.

! report cross-correlation and range	
		write(morestring,'(F6.3)') cc(k1,k1)
		stringfrag = "echo 0.615 1.2 8 0 0 LB Cross-correlation is "//trim(adjustl(morestring))//&
			" | gmt pstext -JX3.5i/2.75i -R0/1/0/10 "
		systemstring = TRIM(stringfrag)//" -W255/255/255 -O -K -P  >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(systemstring))
	
		write(s1,'(F6.3)') fmin
		write(s2,'(F6.3)') fmax
		morestring = "between "//trim(adjustl(s1))//" and "//trim(adjustl(s2))
		stringfrag = "echo 0.615 0.55 8 0 0 LB "//trim(adjustl(morestring))//&
			" | gmt pstext -JX3.5i/2.75i -R0/1/0/10 "
		systemstring = TRIM(stringfrag)//" -W255/255/255 -O -K -P  >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(systemstring))

! title plot		
		stringfrag = "echo 0.38 10.35 10 0 0 LB Sample "//trim(adjustl(samplename))//" | gmt pstext -JX3.5i/2.75i -R0/1/0/10"
		systemstring = TRIM(stringfrag)//" -W255/255/255 -N -O -K -P  >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(systemstring))	

! draw legend
	! age line
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'0.08 8.85'
		write (99,*)'0.21 8.85'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R0/1/0/10 -W3p,100/100/255 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))
		
		stringfrag = "echo 0.24 8.73 9 0 0 LB age | gmt pstext -JX3.5i/2.75i -R0/1/0/10"
		systemstring = TRIM(stringfrag)//" -W255/255/255 -N -O -K -P  >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(systemstring))			

	! logRRo line
		open(unit=99,file='line.scr',status='UNKNOWN')
		write (99,*)'0.08 8.1'
		write (99,*)'0.21 8.1'
		close(99)
		stringfrag = "gmt psxy line.scr -A -JX3.5i/2.75i -R0/1/0/10 -W3p,255/100/100 -O -K -P  >> "&
		    //TRIM(samplename)//"-correlation.ps"
		call system(TRIM(stringfrag))
		
		stringfrag = "echo 0.24 7.98 9 0 0 LB 'log10(R/R@-o@-)' | gmt pstext -JX3.5i/2.75i -R0/1/0/10"
		systemstring = TRIM(stringfrag)//" -W255/255/255 -N -O -P  >> "//TRIM(samplename)//"-correlation.ps"
		call system(TRIM(systemstring))	

! wipe scratch file	
		call system('rm line.scr')

! open plot
		stringfrag = TRIM(samplename)//"-correlation.ps"
		systemstring = "gmt psconvert "//stringfrag//"-Tf -Z "
		call system(TRIM(systemstring))
		stringfrag = "open "//TRIM(samplename)//"-correlation.pdf"
		call system(TRIM(stringfrag))
	
	end if
		
	stop
end

subroutine correl(data1,data2,n,ans)

	parameter(nmax = 8192)
	dimension data1(n),data2(n)
	complex fft(nmax),ans(n)
	
	call twofft(data1,data2,fft,ans,n)
	
	no2 = float(n)/2.0
	
	do i = 1,n/2+1
		ans(i) = fft(i)*conjg(ans(i))/no2
	end do

	ans(1) = cmplx(real(ans(1)),real(ans(n/2+1)))
	
	call realft(ans,n/2,-1)
	
	return
end

subroutine twofft(data1,data2,fft1,fft2,n)

	dimension data1(n),data2(n)
	complex fft1(n),fft2(n),h1,h2,c1,c2
	
	c1 = cmplx(0.5,0.0)
	c2 = cmplx(0.0,-0.5)
	
	do j = 1,n
		fft1(j) = cmplx(data1(j),data2(j))
	end do
	
      call four1(fft1,n,1)
      
      fft2(1) = cmplx(aimag(fft1(1)),0.0)
      fft1(1) = cmplx(real(fft1(1)),0.0)
      n2 = n+2
      
	do j = 2,n/2+1
		h1 = c1*(fft1(j)+conjg(fft1(n2-j)))
		h2 = c2*(fft1(j)-conjg(fft1(n2-j)))
		fft1(j) = h1
		fft1(n2-j) = conjg(h1)
		fft2(j) = h2
		fft2(n2-j) = conjg(h2)
	end do
	
      return
end

subroutine realft(data,n,isign)

	real*8 wr,wi,wpr,wpi,wtemp,theta
	dimension data(*)

	theta = 6.28318530717959d0/2.0d0/dble(n)
	c1 = 0.5
	
	if (isign.eq.1) then
		c2 = -0.5
		call four1(data,n,+1)
	else
		c2 = 0.5
		theta = -theta
	endif
	
	wpr = -2.0d0*dsin(0.5d0*theta)**2
	wpi = dsin(theta)
	wr = 1.0d0+wpr
	wi = wpi
	n2p3 = 2*n+3
	
	do i = 2,n/2+1
		i1 = 2*i-1
		i2 = i1+1
		i3 = n2p3-i2
		i4 = i3+1
		wrs = sngl(wr)
		wis = sngl(wi)
		h1r = c1*(data(i1)+data(i3))
		h1i = c1*(data(i2)-data(i4))
		h2r = -c2*(data(i2)+data(i4))
		h2i = c2*(data(i1)-data(i3))
		data(i1) = h1r+wrs*h2r-wis*h2i
		data(i2) = h1i+wrs*h2i+wis*h2r
		data(i3) = h1r-wrs*h2r+wis*h2i
		data(i4) = -h1i+wrs*h2i+wis*h2r
		wtemp = wr
		wr = wr*wpr-wi*wpi+wr
		wi = wi*wpr+wtemp*wpi+wi
	end do
	
	if (isign.eq.1) then
		h1r = data(1)
		data(1) = h1r+data(2)
		data(2) = h1r-data(2)
	else
		h1r = data(1)
		data(1) = c1*(h1r+data(2))
		data(2) = c1*(h1r-data(2))
		call four1(data,n,-1)
	endif
	
    return
end

subroutine four1(data,nn,isign)

	real*8 wr,wi,wpr,wpi,wtemp,theta
	dimension data(*)
	
	n = 2*nn
	j = 1
	
	do i = 1,n,2
		if (j.gt.i) then
			tempr = data(j)
			tempi = data(j+1)
			data(j) = data(i)
			data(j+1) = data(i+1)
			data(i) = tempr
			data(i+1) = tempi
		endif
		
		m = n/2

1		if ((m.ge.2).and.(j.gt.m)) then
			j = j-m
			m = m/2
			go to 1
		endif
	
		j = j+m
	end do
	
	mmax = 2
	
2	if (n.gt.mmax) then
		istep = 2*mmax
		theta = 6.28318530717959d0/(isign*mmax)
		wpr = -2.d0*dsin(0.5d0*theta)**2
		wpi = dsin(theta)
		wr = 1.d0
		wi = 0.d0
		
		do m = 1,mmax,2
			do i = m,n,istep
				j = i+mmax
				tempr = sngl(wr)*data(j)-sngl(wi)*data(j+1)
				tempi = sngl(wr)*data(j+1)+sngl(wi)*data(j)
				data(j) = data(i)-tempr
				data(j+1) = data(i+1)-tempi
				data(i) = data(i)+tempr
				data(i+1) = data(i+1)+tempi
			end do
			
			wtemp = wr
			wr = wr*wpr-wi*wpi+wr
			wi = wi*wpr+wtemp*wpi+wi
			
		end do
		
		mmax = istep
		go to 2
		
	endif
	return
end

subroutine fourn(data,nn,ndim,isign)

	real*8 wr,wi,wpr,wpi,wtemp,theta
	dimension nn(ndim),data(*)
	
	ntot = 1
	
	do idim = 1,ndim
		ntot = ntot*nn(idim)
	end do
	
	nprev = 1
	
	do idim = 1,ndim
        n = nn(idim)
        nrem = ntot/(n*nprev)
        ip1 = 2*nprev
        ip2 = ip1*n
        ip3 = ip2*nrem
        i2rev = 1
		do i2 = 1,ip2,ip1
			if (i2.lt.i2rev) then
				do i1 = i2,i2+ip1-2,2
					do i3 = i1,ip3,ip2
						i3rev = i2rev+i3-i2
						tempr = data(i3)
						tempi = data(i3+1)
						data(i3) = data(i3rev)
						data(i3+1) = data(i3rev+1)
						data(i3rev) = tempr
						data(i3rev+1) = tempi
					end do
				end do
			endif
			
			ibit = ip2/2
			
1			if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
				i2rev = i2rev-ibit
				ibit = ibit/2
				go to 1
			endif
			
			i2rev = i2rev+ibit
		end do
		
        ifp1 = ip1
        
2       if (ifp1.lt.ip2) then
			ifp2 = 2*ifp1
			theta = isign*6.28318530717959d0/(ifp2/ip1)
			wpr = -2.d0*dsin(0.5d0*theta)**2
			wpi = dsin(theta)
			wr = 1.d0
			wi = 0.d0
			
			do i3 = 1,ifp1,ip1
				do i1 = i3,i3+ip1-2,2
					do i2 = i1,ip3,ifp2
						k1 = i2
						k2 = k1+ifp1
						tempr = sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
						tempi = sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
						data(k2) = data(k1)-tempr
						data(k2+1) = data(k1+1)-tempi
						data(k1) = data(k1)+tempr
						data(k1+1) = data(k1+1)+tempi
					end do
				end do
				wtemp = wr
				wr = wr*wpr-wi*wpi+wr
				wi = wi*wpr+wtemp*wpi+wi
			end do
			
			ifp1 = ifp2
			go to 2
        endif
        
        nprev = n*nprev
	end do !idim
	
	return
end

subroutine chebft(a,b,c,n,func,r39,f39,ni,fmin)

	parameter (nmax = 100, pi = 3.141592653589793d0, ns = 100)
	dimension c(n),f(nmax),r39(0:ns),f39(0:ns)

	bma = 0.5d00*(b-a)
	bpa = 0.5d00*(b+a)
	
	do k = 1,n
		y = dcos(pi*(k-0.5d0)/n)
		f(k) = func(r39,f39,ni,fmin,y*bma+bpa)
	end do
	
	fac = 2.d0/n
      
	do j = 1,n
        sum = 0.d0
		do k = 1,n
			sum = sum+f(k)*cos((pi*(j-1))*((k-0.5d0)/n))
		end do
		c(j) = fac*sum
	end do

	return
end

function chebev(a,b,c,m,x)

	dimension c(m)
	
	if ((x-a)*(x-b).gt.0.) then
		stop 'ERROR(CHEBEV): x not in range .'
	endif
	
	d = 0.d0
	dd = 0.d0
	y = (2.d0*x-a-b)/(b-a)
	y2 = 2.d0*y
	
	do j = m,2,-1
		sv = d
		d = y2*d-dd+c(j)
		dd = sv
	end do

	chebev = y*d-dd+0.5d0*c(1)
	
	return
end

function favr(fx,fy,ni,fmin,xa)

	parameter (ns = 100)
	dimension fx(0:ns),fy(0:ns)
	
	do j = 1,ni
		x = xa+fmin
		if (fx(j).ge.x.and.fx(j-1).lt.x) then
			slope  =   (fy(j)-fy(j-1))/(fx(j)-fx(j-1))
			favr = fy(j-1)+ slope * (x-fx(j-1))
			return
		endif
	end do
end

subroutine chint(a,b,c,cint,n)
	dimension c(n),cint(n)
	
	con = 0.25*(b-a)
	sum = 0.
	fac = 1.
	
	do j = 2,n-1
		cint(j) = con*(c(j-1)-c(j+1))/(j-1)
		sum = sum+fac*cint(j)
		fac = -fac
	end do

	cint(n) = con*c(n-1)/(n-1)
	sum = sum+fac*cint(n)
	cint(1) = 2.*sum
	
	return
end
