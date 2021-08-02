!	Program to calculate 2D temperature profiles
!
!
!     Developer and contact address: E. Rivera-Calderón & O.M. Espinoza-Ojeda, Instituto de Investigaciones en Ciencias de la Tierra, Universidad Michoacana de San Nicolás de Hidalgo, Morelia, Michoacán 58060, México.
!	Telephone number and e-mail: +52 443 269 8138; omespinozaoj@conacyt.mx
!	Required hardware: Intel core i7; 16 GB RAM; or similar
!	Required software: Fortran compiler
!
**********************************************************************************************
!	Statement of 2D arrays
!	DECLARACIÓN DE ARREGLOS PARA PROGRAMA EN DOS DIMENSIONES
!
**********************************************************************************************

	allocatable :: tnew(:,:), tc(:,:)
	allocatable :: a(:,:), b(:,:), c(:,:)
	allocatable :: p(:,:), q(:,:)
	allocatable :: d(:,:), DELTATEM(:,:)
	allocatable :: k(:,:), alfax(:,:)
	allocatable :: alfay(:,:), twf(:,:)
	allocatable :: tsf(:,:), tef(:,:)
	allocatable :: tnf(:,:)

	real  ka, kb, kc, kd, ke, kf, kg, kh, ki, kj
	real  kl, km, kn, ko, kp, kr, ks, kt, ku


**********************************************************************************************
!	Western boundary arrays (Tw = Tc1)
!      ARREGLOS PARA FRONTERA OESTE  (Tw = Tc1)
!     (DEPENDERA DE LA CONFIGURACIÓN DEL POZO)
!
**********************************************************************************************

	allocatable :: aw(:), bw (:), cw (:), dw(:)
	allocatable :: Tc1(:), Pw(:)
      double precision Ti, kz
	allocatable :: T(:), Qw(:), keac1(:)

      real  Tn, Ts, nvc, dx, lonx, wa
	real m1, Tnvc, m, T0, f, kq
	
		
**********************************************************************************************
!	Southern boundary arrays (Ts = Tc2)
!	
!     ARREGLOS ARREGLOS PARA FRONTERA SUR Ts = Tc2
!   (SU CONFIGURACIÓN DEPENDE DE LOS LIMITES DE LOS POZOS W-E A LA DE PROFUNDIDAD)
!
**********************************************************************************************
       
	allocatable :: as(:), bs(:), cs(:), ds(:)
	allocatable :: Qs(:), k2(:), Tc2(:), Ps(:)
	 	 
	real  Tws, Tes, nvcs, dxs, lonxs, was, Tis


**********************************************************************************************
!	Eastern boundary arrays (Te = Tc3)
!
!      ARREGLOS PARA LA FRONTERA ESTE  Te = Tc3
!     (DEPENDERA DE LA CONFIGURACIÓN DEL POZO)
!
**********************************************************************************************

	allocatable :: ae(:), be(:), ce(:), de(:)
	allocatable :: Qe(:), keac2(:)
	allocatable :: Tc3(:), Pe(:)

      real  Tne, Tse, nvce, dxe, lonxe, wae
	real  fe, kqe, me, m1e, qf


**********************************************************************************************
!	Statement of output files
!     ARCHIVO QUE GENERA EL PROGRAMA 2D CON LOS DATOS  DE  TEMPERATURA  CALCULADOS  DURANTE SU
!     EJECUCIÓN, EN EL CUAL SE ESPECIFICA LA UNIDAD DONDE SE GENERA(PUEDE SER 1, 2, 3, ETC)  Y
!     EN: file= ' AQUI VA NOMBRE DEL ARCHIVO', EL ESTADO: status='NUEVO (new) o ANTIGUO (old)'
!     PARA NUESTRO CASO ES NUEVO CUANDO QUEREMOS QUE GENERE DATOS.
!
**********************************************************************************************

	open (unit=1,file='twf_1.txt',status='unknown')
	open (unit=2,file='tef_1.txt',status='unknown')
	open (unit=3,file='Name_2D_TC.txt',status='unknown')
	 
**********************************************************************************************
!	Statement of the thermal conductivity configuration
!
!     DEFINICIÓN DE LAS VARIABLES DE ENTRADA PARA ESPECIFICAR LA CONDUCTIVIDAD TÉRMICA DEL 2D, 
!     CONSIDERANDO CAMBIOS DE MATERIALES: ESTRATOS ROCOSOS O ESTRATIGRAFIA DE CADA POZO.
!     
!	LOS VALORES DE LAS CONDUCTIVIDADES K, DEPENDEN DEL TIPO DE ROCA.
!
**********************************************************************************************
!	Size of the 2D mesh
!     xlong= Distance between Tw and Te (m)
!     ylong= Total depth for calculation of the thermal profile(m)
!     nvc= Number of volume control in x and y
!
!     DIMENSIONES QUE TENDRA EL PERFIL DE TEMPERATURA
!     xlong= DISTANCIA ENTRE LOS POZOS Tw y Te (m)
!     ylong= PROFUNDIDAD A LA QUE SE DESEA CALCULAR EL PERFIL TERMICO (m)
!     nvc= NUMERO DE VOLUMENES DE CONTROL EN x Y EN y
!    
**********************************************************************************************

	xlong=560
	ylong=3000
	nvcx=560
	nvcy=3000

	allocate (k(nvcx,nvcy))

!      Thermal conductivity values for western boundary (Tw) (as example EAC1 wellbore)

      ka= 1.7
	kb= 2
	kc= 3.6
	kd= 2.5
	ke= 2.5
	kf= 2
	kg= 2.5
	kh= 3
	ki= 2
	kj= 3.05

!      Thermal conductivity values for eastern boundary (Te) (as example EAC2 wellbore)

	kl= 3.6
      km= 2.26
	kn= 2
	ko= 2.26
	kp= 2.8
	kr= 2.1
	ks= 3.05	

!      2D Thermal conductivity configuration for western boundary (Tw) (as example EAC1 wellbore) 


	do i = 1,230
		do j= 1,130
	k(i,j)= ka
		end do
	enddo


	do i = 1,230
		do j= 131,210
	  k(i,j)= kb
		end do
      enddo
    

      do i = 1,230
		do j= 211,240
	  k(i,j)= kc
		end do
      enddo

      do i = 1,230
		do j= 241,790
	  k(i,j)= kd
		end do
      enddo

	do i = 1,230
		do j= 791,1200
	  k(i,j)= ke
		end do
      enddo

	do i = 1,230
		do j= 1201,1290
	  k(i,j)= kf
		end do
      enddo

	do i = 1,230
		do j= 1291,1550
	  k(i,j)= kg
		end do
      enddo


	do i = 1,230
		do j= 1551,1600
	  k(i,j)= kh
		end do
      enddo

	do i = 1,230
		do j= 1601,1650
	  k(i,j)= ki
		end do
      enddo

      do i = 1,230
		do j= 1651,3000
	  k(i,j)= kj
		end do
      enddo

!      2D Thermal conductivity configuration for eastern boundary (Te) (as example EAC2 wellbore)


	do i = 231,560
		do j= 1,120
	  k(i,j)= kl
		end do
      enddo


       do i = 231,560
		do j= 121,160
	  k(i,j)= km
		end do
       enddo
       

        do i = 231,560
		do j= 161,290
	  k(i,j)= kn
		end do
       enddo


        do i = 231,560
		do j= 291,350
	  k(i,j)= ko
		end do
       enddo

       do i = 231,560
		do j= 351,450
	  k(i,j)= kp
		end do
       enddo

       do i = 231,560
		do j= 451,1560
	  k(i,j)= kr
		end do
       enddo

    
	do i = 231,560
		do j= 1561,3000
	  k(i,j)= ks
		end do
       enddo

	
**********************************************************************************************
!
!      Western boundary:  Tw = Tc1  (Dirichlet and/or Neumann boundary conditions; Lithology profile)
!
**********************************************************************************************
	
	allocate (keac1(nvcy))
	

!	 1D Numerical configuration of the lithology (thermal conductivity values)
    

	do i= 1,130

	keac1(i)= ka

	enddo

	do i= 131,210

	keac1(i)= kb

      enddo

	do i= 211,240
	  
	  keac1(i)= kc

	enddo

	do i= 241,790
	  
	  keac1(i)= kd

	enddo

	do i= 791,1200
	  
	  keac1(i)= ke

	enddo

	do i= 1201,1290
	  
	  keac1(i)= kf

	enddo

	do i= 1291,1550
	  
	  keac1(i)= kg

	enddo


	do i= 1551,1600
	  
	  keac1(i)= kh

	enddo


	do i= 1601,1650
	  
	  keac1(i)= ki

	enddo

	do i= 1651,3000
	  
	  keac1(i)= kj

	enddo


!     Average surface temperature

	Tn= 25

!	Bottomhole temperature (last measurement)

       Ts= 307.3

!     nvc: Number of volume control

       nvc= 1970

!	Total depth logged (or depth of bottomhole temperature) 

      lonx= 1970

!	Node size

	dx= lonx/nvc

!	Penultimate volume control
      wa= nvc-1

!	Number of volume control to extrapolate the western temperatures 

      m1= nvc+1
	m= nvc+1030

!	Average surface conductive heat flow (W/m2)

       f=0.4268

!     Average geothermal gradient (oC/m)

	 g=0.1435

!	Harmonic thermal conductivity (W/m K)

       kq=2.98

	
**********************************************************************************************
!	Calculation of the 1D thermal profile for the western boundary by using TDMA
!                    
!               CALCULO DE FLUJO DE CALOR EN EL POZO OESTE W EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
!
**********************************************************************************************
!     Parameter bw of TDMA
**********************************************************************************************

	allocate (aw(nvcy), bw(nvcy), cw(nvcy), dw(nvcy))

!	First to penultimate Node

      do i=1,wa

		bw(i)=(keac1(i+1)/dx)

      enddo

!	Last Node (boundary)

	bw(nvc)=0

**********************************************************************************************
!     Parameter cw of TDMA
**********************************************************************************************

!	First Node
   
   	cw(1) = 0

!	Second to last Node

      do i=2,nvc

		cw(i)=(keac1(i-1)/dx)

      enddo

**********************************************************************************************
!     Parameter aw of TDMA
**********************************************************************************************

!	First Node
   
   	aw(1) = (((keac1(2)/dx)) + (2*(keac1(1)/dx)))

!	Central Nodes (second to penultimate)

       do i=2,wa

		aw(i)=(((keac1(i+1)/dx)) + ((keac1(i-1)/dx)))

      enddo

!	Last Node

	aw(nvc)=((2*(keac1(nvc)/dx)) + ((keac1(nvc-1)/dx)))

**********************************************************************************************
!     Parameter dw of TDMA
**********************************************************************************************

!	First Node
   
   	dw(1)=(Tn*(2*(keac1(1)/dx)))

!	Central Nodes (second to penultimate)

       do i=2,wa

		dw(i)= 0

       enddo

!	Last Node

	dw(nvc)=(Ts*(2*(keac1(nvc)/dx)))

**********************************************************************************************
!     Parameter Pw of TDMA
**********************************************************************************************
    
	allocate (Tc1(nvcy), Pw(nvcy), T(nvcy), Qw(nvcy))

!	First Node

	Pw(1)=bw(1)/aw(1)

!	Central Nodes (second to penultimate)

	do i=2,wa

		Pw(i)= bw(i)/(aw(i)-cw(i)*Pw(i-1))

	enddo

!	Last Node

	Pw(nvc) = 0

**********************************************************************************************
!     Parameter Qw of TDMA
**********************************************************************************************
    	
!	First Node

	Qw(1)=dw(1)/aw(1)

!	Second to last Node
	
	do i=2,nvc

		Qw(i) = (dw(i) + cw(i)*Qw(i-1))/(aw(i)-cw(i)*Pw(i-1))
		
      enddo

**********************************************************************************************
!     Calculation to extrapolate the western boundary temperatures
**********************************************************************************************

	   Tc1(nvc)=Qw(nvc)

	do i=wa,1,-1


!	Boundary condition values
        
        Tc1(100)=45
	  Tc1(1970)=307.3
	 	  
	    Tc1(i)=(Pw(i)*Tc1(i+1)) + Qw(i)

	enddo

      do i=1,m

		If (i <= wa) then 
		
		Tc1(i)=(Pw(i)*Tc1(i+1)) + Qw(i)
		 
			else if (i > wa) then 
		
	Tc1(i) =  (f*1/kq) + Tc1(i-1)
		
		endif

	enddo

**********************************************************************************************
!	Calculation of 1D temperature profile in the southern boundary (Ts)
!                    
!               CALCULO DE FLUJO DE CALOR EN LA FRONTERA SUR  Ts = Tc2 EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
!
**********************************************************************************************

	allocate (k2(nvcx))
	allocate (as(nvcx), bs(nvcx), cs(nvcx), ds(nvcx))
	allocate (Qs(nvcx), Tc2(nvcx), Ps(nvcx))

!	Thermal conductivity value at final depth (e.g., basement: 3000 m)
    
      k2= 2.72

!     Extrapolated western temperature at final depht Tc1 (e.g., basement: 3000 m)

      Tws= 367.47

!     Extrapolated eastern temperature at final depht Tc3 (e.g., basement: 3000 m) 

      Tes= 292.84

!     Number of volume control (Distance between Tw and Te)

      nvcs= 560

!	Lenght (Distance between Tw and Te)

      lonxs= 560

!	Size of differential volume control

	dxs= lonxs/nvcs

      was= nvcs-1


**********************************************************************************************
!	Parameter bs of TDMA
**********************************************************************************************

!	First to penultimate Node

      do i=1,was

		bs(i)=(k2(i+1)/dxs)

      enddo

!	Last Node

	bs(nvcs)=0

**********************************************************************************************
!	Parameter cs of TDMA
**********************************************************************************************

!	First Node
 
   	cs(1) = 0

!	Second to last Node

      do i=2,nvcs

		cs(i)=(k2(i-1)/dxs)

      enddo

**********************************************************************************************
!	Parameter as of TDMA
**********************************************************************************************

!	First Node
   
   	as(1) = (((k2(2)/dxs)) + (2*(k2(1)/dxs)))

!	Central Nodes (second to penultimate)

      do i=2,was

		as(i)=(((k2(i+1)/dxs)) + ((k2(i-1)/dxs)))

      enddo

!	Last Node

	as(nvcs)=((2*(k2(nvcs)/dxs)) + ((k2(nvcs-1)/dxs)))

**********************************************************************************************
!	Parameter ds of TDMA
**********************************************************************************************

!	First Node
   
   	ds(1)=(Tws*(2*(k2(1)/dxs)))

!	Central Nodes (second to penultimate)

      do i=2,was

		ds(i)= 0

      enddo

!	Last Node

	ds(nvcs)=(Tes*(2*(k2(nvcs)/dxs)))

**********************************************************************************************
!	Parameter Ps of TDMA
**********************************************************************************************

!	First Node
    
	Ps(1)=bs(1)/as(1)

!	Central Nodes (second to penultimate)

	do i=2,was

		Ps(i)= bs(i)/(as(i)-cs(i)*Ps(i-1))

	enddo

!	Last Node

	Ps(nvcs) = 0

**********************************************************************************************
!	Parameter Qs of TDMA
**********************************************************************************************
    	
!	First Node

	Qs(1)=ds(1)/as(1)
	
!	Second to last Node

	do i=2,nvcs
	
		Qs(i) = (ds(i) + cs(i)*Qs(i-1))/(as(i)-cs(i)*Ps(i-1))
			
      enddo

**********************************************************************************************
!	Southern boundary temperature calculation
**********************************************************************************************

	Tc2(nvcs)=Qs(nvcs)

	do i=was,1,-1

	     Tc2(i)=(Ps(i)*Tc2(i+1)) + Qs(i)

	enddo

**********************************************************************************************
!	Calculation of the 1D thermal profile for the eastern boundary Te by using TDMA
!                    
!               CALCULO DEl PERFIL DE TEMPERATURA EN LA FRONTERA ESTE  Te = Tc3 EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
**********************************************************************************************
!
!	Input Data
!
**********************************************************************************************

	allocate (keac2(nvcy))
	allocate (ae(nvcy), be(nvcy), ce(nvcy), de(nvcy))
	allocate (Qe(nvcy), Pe(nvcy), Tc3(nvcy))


!	1D Numerical configuration of the lithology (thermal conductivity values)

	do i= 1,120
	  
	  keac2(i)= kl

	enddo

  
    	do i= 121,160
	  
	  keac2(i)= km

	enddo


	do i= 161,290
	  
	  keac2(i)= kn

	enddo


	do i= 291,350
	  
	  keac2(i)= ko

	enddo


	do i= 351,450
	  
	  keac2(i)= kp

	enddo


	do i= 451,1560
	  
	  keac2(i)= kr

	enddo


	do i= 1561,3000
	  
	  keac2(i)= ks

	enddo


!     Average surface temperature

      Tne= 25

!	Bottomhole temperature (last measurement)

      Tse=240

!	nvc: Number of volume control

      nvce= 1850

!	Total depth logged (or depth of bottomhole temperature)

      lonxe= 1850

!     Node size

	dxe= lonxe/nvce

!	Penultimate volume control

      wae= nvce-1

!	Number of volume control to extrapolate the western temperatures 

      m1e= nvce+1

	me= nvce+1150

!	Average surface conductive heat flow (W/m2)

      fe=0.2632

!	Average geothermal gradient (oC/m)

	ge=0.1072

!     Harmonic thermal conductivity (W/m K)

      kqe=2.46


**********************************************************************************************
!	Parameter be of TDMA
**********************************************************************************************

!	First to penultimate Node

      do i=1,wae

		be(i)=(keac2(i+1)/dxe)

      enddo

!	Last Node

	be(nvce)=0

!   *********************************************************************************************
!	Parameter ce of TDMA
!   *********************************************************************************************

!	First Node
   
   	ce(1) = 0

!	Second to last Node

      do i=2,nvce

		ce(i)=(keac2(i-1)/dxe)

      enddo

**********************************************************************************************
!	Parameter ae of TDMA
**********************************************************************************************

!	First Node
   
   	ae(1) = (((keac2(2)/dxe)) + (2*(keac2(1)/dxe)))

!	Central Nodes (second to penultimate)

      do i=2,wae

		ae(i)=(((keac2(i+1)/dxe)) + ((keac2(i-1)/dxe)))

      enddo

!	Last Node

	ae(nvce)=((2*(keac2(nvce)/dxe)) + ((keac2(nvce-1)/dxe)))

**********************************************************************************************
!	Parameter de of TDMA
**********************************************************************************************

!	First Node
 
   	de(1)=(Tne*(2*(keac2(1)/dxe)))

!	Central Nodes (second to penultimate)

      do i=2,wae

		de(i)= 0

      enddo

!	Last Node

	de(nvce)=(Tse*(2*(keac2(nvce)/dxe)))

**********************************************************************************************
!	Parameter Pe of TDMA
**********************************************************************************************
    
!	First Node

	Pe(1)=be(1)/ae(1)

!	Central Nodes (second to penultimate)

	do i=2,wae

		Pe(i)= be(i)/(ae(i)-ce(i)*Pe(i-1))

	enddo

!	Last Node

	Pe(nvce) = 0

**********************************************************************************************
!	Parameter Qe of TDMA
**********************************************************************************************
    	
!	First Node

	Qe(1)=de(1)/ae(1)
	
!	Second to last Node

	do i=2,nvce

		Qe(i) = (de(i) + ce(i)*Qe(i-1))/(ae(i)-ce(i)*Pe(i-1))
		
      enddo

**********************************************************************************************
!	Calculation to extrapolate the eastern boundary temperatures Te = Tc3
**********************************************************************************************

	 Tc3(nvce)=Qe(nvce)

!	Boundary condition values

	 Tc3(20)=40
	 Tc3(190)=263.77

	do i=wae,1,-1

	        Tc3(i)=(Pe(i)*Tc3(i+1)) + Qe(i)

	enddo

      do i=1,me

          If (i <= wae) then 
		
		    Tc3(i)=(Pe(i)*Tc3(i+1)) + Qe(i)
		 
		      else if (i > wae) then 		
		
	 	    Tc3(i) = (fe*1/kqe) + Tc3(i-1)

		endif
	
	enddo

**********************************************************************************************
!	Calculation of the 2D temperature profile
!
!						CALCULO DEL PERFIL DE TEMPERATURAS EN 2D 
**********************************************************************************************

	allocate (twf(nvcx,nvcy), tnf(nvcx,nvcy), tsf(nvcx,nvcy))
	allocate (tef(nvcx,nvcy))
	allocate (alfax(nvcx,nvcy), alfay(nvcx,nvcy), tc(nvcx,nvcy))
	allocate (a(nvcx,nvcy), b(nvcx,nvcy), c(nvcx,nvcy))
	allocate (d(nvcx,nvcy))
	allocate (p(nvcx,nvcy), q(nvcx,nvcy))
	allocate (tnew(nvcx,nvcy), DELTATEM(nvcx,nvcy))

!	Statement of boundary temperature values in the 2D Model

!	Western boundary	

      do j= 1,3000
		do i=1,1
	
			twf(i,j)=Tc1(j)

		enddo
	enddo

 !	Southern boundary	     
	
	do j= 1,1
		do i=1,560
	
			tsf(i,j)=25
      
		enddo
	enddo

!	Eastern boundary	

      do j= 1,3000
		do i=560,560

			tef(i,j)=Tc3(j)
      
		enddo
	enddo

!	Northern boundary	

      do j= 3000,3000
		do i=1,560

			tnf(i,j)=Tc2(i)

		enddo
	enddo
		
!     Calculation of the Node's area (e.g., 1x1 m, 10x10 m, 100x100 m)

	deltax=xlong/nvcx
	deltay=ylong/nvcy

	write(*,*)'deltax-deltay', deltax, deltay

	n=nvcx
	m=nvcy	
		
	write(*,*)n,m

!     Calculation of dimensional thermal conductivity

      do i=1,560
		do j=1,3000

			alfax(i,j)=k(i,j)*deltay/deltax
			alfay(i,j)=k(i,j)*deltax/deltay

		enddo
	enddo

**********************************************************************************************
!     Initial temperature field
**********************************************************************************************

	do i=1, n
		do j= 1, m

		 tc(i,j) =50.0

		enddo
	enddo

**********************************************************************************************
!	 2D Parameter a of TDMA
**********************************************************************************************

	cont=1

15	write(*,*) cont

      a(1,1)= (3*alfax(1,1))+(3*alfay(1,1)) 
   
	a(1,m)= (3*alfax(1,m))+(3*alfay(1,m)) 

	a(n,1)= (3*alfax(n,1))+(3*alfay(n,1))
	 
	a(n,m)= (3*alfax(n,m))+(3*alfay(n,m))  

	
	do j=2,m-1
	a(1,j)=(3*alfax(1,j))+(2*alfay(1,j))
	a(n,j)=(3*alfax(n,j))+(2*alfay(n,j))
	enddo

	do i=2,n-1
	a(i,1)=(2*alfax(i,1))+(3*alfay(i,1))
	a(i,m)=(2*alfax(i,m))+(3*alfay(i,m))	
	enddo

	do j=2, m-1
		do i=2, n-1
		a(i,j)= 2.0*(alfax(i,j)+alfay(i,j))
		enddo
	enddo


**********************************************************************************************
!	 2D Parameter b of TDMA
**********************************************************************************************
	do j=1, m
		do i=1, n-1
			b(i,j)= alfax(i,j)
		enddo
	enddo

	do j=1,m
	b(n,j)=0.0
	enddo

**********************************************************************************************
!	 2D Parameter c of TDMA
**********************************************************************************************

	do j=1, m
		do i=2, n
			c(i,j)= alfax(i,j)
		enddo
	enddo

	do j=1,m
	c(1,j)=0.0
	enddo

**********************************************************************************************
!	 2D Parameter d of TDMA
**********************************************************************************************

	d(1,1)=(2*alfax(1,1)*twf(1,1))+(2*alfay(1,1)*tsf(1,1))+(alfay(1,1)
	1*Tc(1,2))
	d(1,m)=(2*alfax(1,m)*twf(1,m))+(2*alfay(1,m)*tnf(1,m))+(alfay(1,m)
	1*Tc(1,m-1))
	d(n,1)=(2*alfax(n,1)*tef(n,1))+(2*alfay(n,1)*tsf(n,1))+(alfay(n,1)
	1*Tc(n,2))
	d(n,m)=(2*alfax(n,m)*tef(n,m))+(2*alfay(n,m)*tnf(n,m))+(alfay(n,m)
	1*Tc(n,m-1))
	
	do j=2,m-1
	d(1,j)= (2*alfax(1,j)*twf(1,j))+(alfay(1,j)*Tc(1,j+1))+(alfay(1,j)
	1*Tc(1,j-1))
	d(n,j)= (2*alfax(n,j)*tef(n,j))+(alfay(n,j)*Tc(n,j+1))+(alfay(n,j)
	1*Tc(n,j-1))
	enddo

	do i=2,n-1
	d(i,1)= (2*alfay(i,1)*tsf(i,1))+(alfay(i,1)*Tc(i,2))
	d(i,m)= (2*alfay(i,m)*tnf(i,m))+(alfay(i,m)*Tc(i,m-1))
	enddo

!     para los nodos centrales
	do j=2, m-1
		do i=2, n-1
			d(i,j)= (alfay(i,j)*tc(i,j+1))+(alfay(i,j)*tc(i,j-1))
!			write(*,*)i,j, d(i,j)
		enddo
	enddo

**********************************************************************************************
!	 2D Parameter P of TDMA
**********************************************************************************************

!     calculation of P(1,j)
	do j=1, m
	p(1,j)= b(1,j)/a(1,j)
	enddo

!     calculation of P(i,j) at central nodes for i 
	do i=2, n
	do j= 1, m
	p(i,j)= b(i,j)/(a(i,j)-(c(i,j)*p(i-1,j)))
	enddo
	enddo

!     calculation of Q(1,1)
	do j=1,m
	q(1,j)= d(1,j)/a(1,j)
	enddo

!     calculation of Q(i,j) at central nodes for i
	do i= 2, n
	do j=1,m
	q(i,j)= (d(i,j)+(c(i,j)*q(i-1,j)))/(a(i,j)-(c(i,j)*p(i-1,j)))
	enddo
	enddo

!	Update and calculation of the temperature profile

	do j=1,m
	Tnew(n,j)= q(n,j)
	enddo

	do j=1, m
	do i=n-1, 1, -1
	Tnew(i,j)= (p(i,j)*tnew(i+1,j))+q(i,j)
	enddo 
	enddo


**********************************************************************************************
!     Convergence criteria
**********************************************************************************************

	do j=1,m
		do i=1, n
		deltatem(i,j) = abs(tnew(i,j) - Tc(i,j))
		enddo
	enddo

	xmaxdeltate = 0.001

 	DO j=1,m
		do i=1, n
		if(deltatem(i,j) .gt. xmaxdeltate) then
			xmaxdeltate = DELTATEM(i,j)
		ENDIF
		enddo
	ENDDO
	

	IF(xMAXDELTATE .GT. 0.001) THEN
		do i=1,n
		        do j=1,m
			tc(i,j) = tnew(i,j)
			enddo
                enddo
		
		xmaxdeltate = 0.001
	cont= cont+1
	
	GO TO 15
	endif

**********************************************************************************************
!     Output files
**********************************************************************************************
	
	DO j=1,m
			do i=1,1
				write (1,105) i, j, twf (i,j)
			enddo
	ENDDO


	DO j=1,m
			do i=n,n
				write (2,105) i, j, tef (i,j)
			enddo
	ENDDO


	DO j=1,m
			do i=1,n
				write (3,105) i, j, k(i,j)
			enddo
	ENDDO


	write(*,*) cont

105	format (i4,1x,i4,1x,f10.5,1x,f10.5)



	deallocate (keac1,aw,bw,cw,dw,Tc1,Pw,T,Qw)
	deallocate (as,bs,cs,ds,Qs,k2,Tc2,Ps)
	deallocate (ae,be,ce,de,Qe,keac2,Tc3,Pe)
	deallocate (tnew,tc,a,b,c,d,p,q,DELTATEM)
	deallocate (k,alfax,alfay,twf,tnf,tsf,tef)



	end
