!	Program to calculate 2D temperature profiles
!
!
!      PROGRAMA PARA CALCULAR PERFILES DE TEMPERATURA EN ESTADO SEMITRANSITORIO
!           EN DOS DIMENSIONES (2D) PARA CAMPOS GEOTÉRMICOS CONDUCTIVOS
!      
!
!
!Observación: los datos de longitudes están a escala de 10, para agilizar el calculo de los 
!datos, por lo que una vez obtenidos se deben multiplicar por 10 para tenerlos a escala real. 
!
!   
!
**********************************************************************************************
!	Statement of 2D arrays
!	 DECLARACIÓN DE ARREGLOS PARA PROGRAMA EN DOS DIMENSIONES
!
**********************************************************************************************


!	double precision tnew(3300,3300), tc(3300,3300)
	allocatable :: tnew(:,:), tc(:,:)
!	double precision a(3000,3000), b(3000,3000), c(3000,3000)
	allocatable :: a(:,:), b(:,:), c(:,:)
!	double precision p(3000,3000), q(3000,3000)
	allocatable :: p(:,:), q(:,:)
!	double precision d(3000,3000), DELTATEM(3300,3300)
	allocatable :: d(:,:), DELTATEM(:,:)

!	double precision k(5000,5000), alfax(3000,3000)
	allocatable :: k(:,:), alfax(:,:)
!	double precision alfay (3000,3000), twf(3000,3000)
	allocatable :: alfay(:,:), twf(:,:)

 !     double precision tsf(3000,3000), tef(3000,3000)
	allocatable :: tsf(:,:), tef(:,:)

  !    double precision tnf(3000,3000)
	allocatable :: tnf(:,:)

	real  ka, kb, kc, kd, ke, kf, kg, kh, ki, kj
	real  kl, km, kn, ko, kp, kr, ks, kt, ku


**********************************************************************************************
!	Western boundary arrays
!      ARREGLOS PARA FRONTERA OESTE  TW = Tc1
!     (DEPENDERA DE LA CONFIGURACIÓN DEL POZO)
!
**********************************************************************************************


!	double precision aw(10000), bw(10000), cw(10000), dw(10000) 
	allocatable :: aw(:), bw (:), cw (:), dw(:)
!	double precision T1(10000), T2(10000), T3, T4(10000)
!	allocatable :: T1(:), T2(:)
	 
!	double precision GG2, Tc1(10000), Pw(10000)
	allocatable :: Tc1(:), Pw(:)
      double precision Ti, kz
	allocatable :: T(:), Qw(:), keac1(:)


      real  Tn, Ts, nvc, dx, lonx, wa
	real m1, Tnvc, m, T0, f, kq
	
		
**********************************************************************************************
!	Southern boundary arrays
!     ARREGLOS ARREGLOS PARA FRONTERA SUR TS = Tc2
!   (SU CONFIGURACIÓN DEPENDE DE LOS LIMITES DE LOS POZOS W-E A LOS 300 metros DE PROFUNDIDAD)
!
**********************************************************************************************
       

!	double precision as(10000), bs(10000), cs(10000), ds(10000)
	allocatable :: as(:), bs(:), cs(:), ds(:)
!	double precision Ts1(10000), Ts2(10000), Ts3(10000), Ts4(10000)
!	allocatable :: Ts1(:), Ts2(:), Ts3(:), Ts4(:)
!	double precision Qs(10000), k2(10000), Tc2(10000), Ps(10000)
	allocatable :: Qs(:), k2(:), Tc2(:), Ps(:)
	 
	 
	 real  Tws, Tes, nvcs, dxs, lonxs, was, Tis


**********************************************************************************************
!	Eastern boundary arrays
!      ARREGLOS PARA LA FRONTERA ESTE  TE = Tc3
!     (DEPENDERA DE LA CONFIGURACIÓN DEL POZO)
!
**********************************************************************************************


!	double precision ae(10000), be(10000), ce(10000), de(10000)
	allocatable :: ae(:), be(:), ce(:), de(:)
!     double precision Qe(10000), keac2(10000), rk (1000) 
	allocatable :: Qe(:), keac2(:)
!     double precision Tc3(10000), Pe(10000)
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
	open (unit=3,file='AC_2D_1.txt',status='unknown')
	 
**********************************************************************************************
!	Statement of the thermal conductivity
!     DEFINICIÓN DE LAS VARIABLES DE ENTRADA PARA ESPECIFICAR LA CONDUCTIVIDAD TÉRMICA DEL 2D, 
!     CONSIDERANDO CAMBIOS DE MATERIALES: ESTRATOS ROCOSOS O ESTRATIGRAFIA DE CADA POZO.
!     
!	LOS VALORES DE LAS CONDUCTIVIDADES K, DEPENDEN DEL TIPO DE ROCA.
!
**********************************************************************************************
!	2D size of the mesh
!    DIMENSIONES QUE TENDRA EL PERFIL DE TEMPERATURA
!
!    xlong= DISTANCIA ENTRE LOS POZOS EN METROS (m)
!    xlong= PROFUNDIDAD A LA QUE SE DESEA CALCULAR LA EXTRAPOLACION EN METROS (m)
!    nvc= NUMERO DE VOLUMENES DE CONTROL EN x Y EN y
!
!    NOTA: RECORDAR QUE LOS DATOS DE LONGITUDES PARA NUESTRO MODELO ESTÁN A ESCALA DE 10
!
!    
**********************************************************************************************


	xlong=560
	ylong=3000
	nvcx=560
	nvcy=3000


	allocate (k(nvcx,nvcy))

!      CONDUCTIVIDADES PARA POZO DE LA FRONTERA OESTE (W) EN EL 2D EAC1


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

!      CONDUCTIVIDADES PARA POZO DE LA FRONTERA ESTE (E) EN EL 2D EAC2

	kl= 3.6
      km= 2.26
	kn= 2
	ko= 2.26
	kp= 2.8
	kr= 2.1
	ks= 3.05


	

!      CONDUCTIVIDADES PARA POZO DE LA FRONTERA OESTE (W) EN EL 2D: VARIABLES 


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

!     CONDUCTIVIDADES PARA POZO DE LA FRONTERA ESTE (E) EN EL 2D:  VARIABLES


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
!      FRONTERA OESTE:  TW = Tc1   (DECLARACIÓN DE ARREGLOS PARA POZO OESTE)
!
**********************************************************************************************
	
	allocate (keac1(nvcy))
	

!	 CONDUCTIVIDADES TERMICAS 1D PARA EL POZO OESTE
    


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


!     TEMPERATURA SUPERFICIAL (DEPENDE DE CADA SITIO O DE LA TEMPERATURA AMBIENTE DEL CAMPO GEOTÉRMICO)


	Tn= 25

!	TEMPERATURA SUR REGISTRADA EN EN POZO (ÚLTIMA TEMP DEL REGISTRO CONDUCTIVO)

       Ts= 307.3

!     nvc: Número de Volumenes de Vontrol

       nvc= 1970

!	LONGITUD DEL POZO OESTE 

      lonx= 1970

!	LONGITUD DEL POZO/NUMERO DE VOLUMENES DE CONTROL

	dx= lonx/nvc


      wa= nvc-1


!	CALCULO DE TEMPERATURAS EXTRAPOLADAS 

      m1= nvc+1
	m= nvc+1030

!


!	FLUJO DE CALOR PROMEDIO PARA EL POZO W (DETERMINADO)


       f=0.4268

!     GRADIENTE GEOTERMICO PROMEDIO DEL POZO W (DETERMINADO)


	 g=0.1435

!	CONDUCTIVIDAD PROMEDIO PARA EL POZO  W (DETERMINADO)


       kq=2.98



!	ECUACIÓN PARA CALCULAR EL FLUJO DE CALOR EN CASO DE NO TENERLO CALCULADO PARA EL POZO
    
!	kz=0


!     do i=1,300
!	kz=kz+(1/k(i))
!	sum (1/k(i))



!      enddo

    
!	qf=(Tn-Ts)/(kz)
	
**********************************************************************************************
!	Calculation of the heat flow for the western boundary, using TDMA
!                    
!               CALCULO DE FLUJO DE CALOR EN EL POZO OESTE W EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
!
**********************************************************************************************
!     PARAMETRO bw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************

	allocate (aw(nvcy), bw(nvcy), cw(nvcy), dw(nvcy))

!	NODO (2, nvc-1)

      do i=1,wa

		bw(i)=(keac1(i+1)/dx)

      enddo


!	NODO (nvc)

	bw(nvc)=0



**********************************************************************************************
!     PARAMETRO cw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1
   
   	cw(1) = 0


!	NODO (2, nvc-1)

      do i=2,nvc

		cw(i)=(keac1(i-1)/dx)

      enddo



**********************************************************************************************
!     PARAMETRO aw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1
   
   	aw(1) = (((keac1(2)/dx)) + (2*(keac1(1)/dx)))



!	NODO (2, nvc-1)

       do i=2,wa

		aw(i)=(((keac1(i+1)/dx)) + ((keac1(i-1)/dx)))

      enddo



!	NODO (nvc)

	aw(nvc)=((2*(keac1(nvc)/dx)) + ((keac1(nvc-1)/dx)))



**********************************************************************************************
!     PARAMETRO dw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1
   
   	dw(1)=(Tn*(2*(keac1(1)/dx)))



!	NODO (2, nvc-1)

       do i=2,wa

		dw(i)= 0

       enddo



!	NODO (nvc)

	dw(nvc)=(Ts*(2*(keac1(nvc)/dx)))



**********************************************************************************************
!     PARAMETRO Pw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************
    
	allocate (Tc1(nvcy), Pw(nvcy), T(nvcy), Qw(nvcy))

	Pw(1)=bw(1)/aw(1)



	do i=2,wa

		Pw(i)= bw(i)/(aw(i)-cw(i)*Pw(i-1))

	enddo



	Pw(nvc) = 0




**********************************************************************************************
!     PARAMETRO Qw DE LA ECUACIÓN DE FLUJO DE CALOR
**********************************************************************************************
    	


	Qw(1)=dw(1)/aw(1)


	
	do i=2,nvc

		Qw(i) = (dw(i) + cw(i)*Qw(i-1))/(aw(i)-cw(i)*Pw(i-1))
		
      enddo



**********************************************************************************************
!     REALIZAR CALCULO DE TEMPERATURAS EXTRAPOLADAS EN EL POZO OESTE = FRONTERA "W" OESTE
**********************************************************************************************



    !  do i=1,300
	   
	
	!	T3 = Ts - ((kz)*(qf))


!	enddo



	   Tc1(nvc)=Qw(nvc)



	do i=wa,1,-1


!	VALORES FIJOS-CONDUCTIVOS DE TEMPERATURA DEL POZO
        
        Tc1(100)=45
    	  Tc1(200)=62
	  Tc1(400)=85
	  Tc1(500)=100
	  Tc1(600)=115
	  Tc1(800)=140

	  Tc1(1000)=168
    	  Tc1(1200)=192
	  Tc1(1400)=220
	  Tc1(1600)=245
	  Tc1(1650)=252
	  Tc1(1700)=261

	  Tc1(1750)=270
    	  Tc1(1800)=280
	  Tc1(1850)=285
	  Tc1(1900)=292
	  Tc1(1950)=300
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
!	Calculation of heat flow in the southern boundary
!                    
!               CALCULO DE FLUJO DE CALOR EN LA FRONTERA SUR  TS = Tc2 EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
!
**********************************************************************************************

	allocate (k2(nvcx))
	allocate (as(nvcx), bs(nvcx), cs(nvcx), ds(nvcx))
	allocate (Qs(nvcx), Tc2(nvcx), Ps(nvcx))

!	CONDUCTIVIDAD TERMICA A LOS 3000 metros (BASAMENTO DEL CAMPO GEOTERMICO)
    
   

      k2= 2.72



!     TEMPERATURA FINAL DEL POZO OESTE Tc1 (m) A LOS 3000 m



      Tws= 367.47



!     TEMPERATURA FINAL DEL POZO ESTE Tc3 (m) A LOS 3000 m 



      Tes= 292.84192



!     NUMERO DE VOLUMENES DE CONTROL = DISTANCIA ENTRE LOS POZOS



      nvcs= 560




!	LONGITUD = DISTANCIA ENTRE LOS POZOS



      lonxs= 560



!	LONGITUD / NUMERO DE VOLUMENES DE CONTROL



	dxs= lonxs/nvcs


      was= nvcs-1




**********************************************************************************************
!	PARAMETRO  bs DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************





!	NODO (2, nvc-1)


      do i=1,was


		bs(i)=(k2(i+1)/dxs)


      enddo



!	NODO (nvc)



	bs(nvcs)=0



**********************************************************************************************
!	PARAMETRO  cs DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1

   
   	cs(1) = 0


!	NODO (2, nvc-1)


      do i=2,nvcs

		cs(i)=(k2(i-1)/dxs)

      enddo




**********************************************************************************************
!	PARAMETRO as DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1
   

   	as(1) = (((k2(2)/dxs)) + (2*(k2(1)/dxs)))


!	NODO (2, nvc-1)


      do i=2,was

		as(i)=(((k2(i+1)/dxs)) + ((k2(i-1)/dxs)))

      enddo


!	NODO (nvc)


	as(nvcs)=((2*(k2(nvcs)/dxs)) + ((k2(nvcs-1)/dxs)))




**********************************************************************************************
!	PARAMETRO ds DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************



!	NODO 1
   

   	ds(1)=(Tws*(2*(k2(1)/dxs)))



!	NODO (2, nvc-1)


      do i=2,was

		ds(i)= 0

      enddo


!	NODO (nvc)


	ds(nvcs)=(Tes*(2*(k2(nvcs)/dxs)))


**********************************************************************************************
!	PARAMETRO Ps DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************


    
	Ps(1)=bs(1)/as(1)



	do i=2,was


		Ps(i)= bs(i)/(as(i)-cs(i)*Ps(i-1))


	enddo



	Ps(nvcs) = 0




**********************************************************************************************
!	PARAMETRO Qs DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************
    	


	Qs(1)=ds(1)/as(1)
	


	do i=2,nvcs

	
		Qs(i) = (ds(i) + cs(i)*Qs(i-1))/(as(i)-cs(i)*Ps(i-1))
	
		
      enddo




**********************************************************************************************
!	CALCULO DE LA TEMPERATURA DE LA FRONTERA SUR
**********************************************************************************************



	Tc2(nvcs)=Qs(nvcs)

	 

	do i=was,1,-1


	     Tc2(i)=(Ps(i)*Tc2(i+1)) + Qs(i)


	enddo




**********************************************************************************************
!	Calculation of heat flow in the eastern boundary, using TDMA
!                    
!               CALCULO DE FLUJO DE CALOR EN LA FRONTERA ESTE  TE = Tc3 EN 1D UTILIZANDO 
!                 ECUACIONES DE TRANSFERENCIA DE CALOR Y EL METODO DE TDMA
!
!
**********************************************************************************************
!
!	ENTRADA DE DATOS
!
**********************************************************************************************

	allocate (keac2(nvcy))
	allocate (ae(nvcy), be(nvcy), ce(nvcy), de(nvcy))
	allocate (Qe(nvcy), Pe(nvcy), Tc3(nvcy))


!	CONDUCTIVIDADES TÈRMICAS DEL POZO ESTE EN 1D


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






!     TEMPERATURA NORTE-SUPERFICIAL



      Tne= 25



!	TEMPERATURA SUR-REGISTRO DEL POZO ESTE



      Tse=240


!	NUMERO DE VOLUMENES DE CONTROL



      nvce= 1850



!	LONGITUD DEL POZO ESTE EN METROS



      lonxe= 1850



!     LONGITUD / NUMERO DE VOLUMENES DE CONTROL



	dxe= lonxe/nvce



      wae= nvce-1



!	CALCULOS DE TEMPERATURAS EXTRAPOLADAS EN EL POZO ESTE



      m1e= nvce+1

	me= nvce+1150



!	FLUJO DE CALOR EN EL POZO ESTE (DETERMINADO)


      fe=0.2632


!	GRADIENTE GEOTÉRMICO DEL POZO ESTE (DETERMINADO)

	ge=0.2632


!     CONDUCTIVIDAD PROMEDIO DEL POZO ESTE

 
      kqe=2.46




!     CALCULO DE FLUJO DE CALOR (EN CASO DE NO CONTAR CON VALOR DETERMINADO)
    
!	kz=0


!     do i=1,300

!	kz=kz+(1/rk(i))

!      enddo


    
!	qf=(Tse-Tne)/(kz)



**********************************************************************************************
!	PARAMETRO  be  DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************





!	NODO (2, nvc-1)


      do i=1,wae

		be(i)=(keac2(i+1)/dxe)

      enddo


!	NODO (nvc)



	be(nvce)=0




!   *********************************************************************************************
!	PARAMETRO  ce  DE LA ECUACION DE FLUJO DE CALOR
!   *********************************************************************************************



!	NODO 1


   
   	ce(1) = 0



!	NODO (2, nvc-1)



      do i=2,nvce


		ce(i)=(keac2(i-1)/dxe)


      enddo




**********************************************************************************************
!	PARAMETRO ae DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************




!	NODO 1
   


   	ae(1) = (((keac2(2)/dxe)) + (2*(keac2(1)/dxe)))



!	NODO (2, nvc-1)



      do i=2,wae

		ae(i)=(((keac2(i+1)/dxe)) + ((keac2(i-1)/dxe)))

      enddo



!	NODO (nvc)



	ae(nvce)=((2*(keac2(nvce)/dxe)) + ((keac2(nvce-1)/dxe)))




**********************************************************************************************
!	PARAMETRO de DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************




!	NODO 1


   
   	de(1)=(Tne*(2*(keac2(1)/dxe)))



!	NODO (2, nvc-1)



      do i=2,wae


		de(i)= 0


      enddo



!	NODO (nvc)



	de(nvce)=(Tse*(2*(keac2(nvce)/dxe)))



**********************************************************************************************
!	PARAMETRO  Pe DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************
    

	Pe(1)=be(1)/ae(1)



	do i=2,wae


		Pe(i)= be(i)/(ae(i)-ce(i)*Pe(i-1))


	enddo



	Pe(nvce) = 0



**********************************************************************************************
!	PARAMETRO Qe DE LA ECUACION DE FLUJO DE CALOR
**********************************************************************************************
    	


	Qe(1)=de(1)/ae(1)
	


	do i=2,nvce


		Qe(i) = (de(i) + ce(i)*Qe(i-1))/(ae(i)-ce(i)*Pe(i-1))
		

      enddo




**********************************************************************************************
!	CALCULAR TEMPERATURAS EXTRAPOLADAS EN EL POZO ESTE  FRONTERA TE = Tc3
**********************************************************************************************


	
!	do i=1,300
	   
	
!		T3= Ts - ((kz)*(qf))


!	enddo



	 Tc3(nvce)=Qe(nvce)



!	VALORES FIJOS-CONDUCTIVOS DE TEMPERATURA DEL POZO

!	 Tc3(20)=40
!	 Tc3(40)=80
!	 Tc3(60)=105
!	 Tc3(80)=125
!	 Tc3(100)=150
!	 Tc3(105)=155
!	 Tc3(110)=160
	 
!	 Tc3(115)=165
!	 Tc3(120)=170
!	 Tc3(125)=175
!	 Tc3(130)=180
!	 Tc3(135)=185
!	 Tc3(140)=190
!	 Tc3(145)=195

	
!	 Tc3(150)=200
!	 Tc3(155)=202
!	 Tc3(160)=210
!	 Tc3(165)=215
!	 Tc3(170)=220
!	 Tc3(175)=225
!	 Tc3(180)=232
!	 Tc3(185)=240
!	 Tc3(190)=263.77



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
!	Calculation of the 2D temperature
!						CALCULO DE TEMPERATURAS EN EL PERFIL 2D 
**********************************************************************************************


	allocate (twf(nvcx,nvcy), tnf(nvcx,nvcy), tsf(nvcx,nvcy))
	allocate (tef(nvcx,nvcy))
	allocate (alfax(nvcx,nvcy), alfay(nvcx,nvcy), tc(nvcx,nvcy))
	allocate (a(nvcx,nvcy), b(nvcx,nvcy), c(nvcx,nvcy))
	allocate (d(nvcx,nvcy))
	allocate (p(nvcx,nvcy), q(nvcx,nvcy))
	allocate (tnew(nvcx,nvcy), DELTATEM(nvcx,nvcy))



!	DECLARACION DE VALORES DE TEMPERATURA DE LAS FRONTERAS DEL MODELO 2D



!	FRONTERA OESTE DEL 2D	

      do j= 1,3000
		do i=1,1
	
			twf(i,j)=Tc1(j)

		enddo
	enddo


 !	FRONTERA SUR DEL 2D	     
	
	do j= 1,1
		do i=1,560
	
			tsf(i,j)=25
      
		enddo
	enddo


!	FRONTERA ESTE DEL 2D	

      do j= 1,3000
		do i=560,560

			tef(i,j)=Tc3(j)
      
		enddo
	enddo


!	FRONTERA NORTE DEL 2D	

      do j= 3000,3000
		do i=1,560

			tnf(i,j)=Tc2(i)

		enddo
	enddo
		


!     calculo de delta x y y


	deltax=xlong/nvcx
	deltay=ylong/nvcy

	write(*,*)'deltax-deltay', deltax, deltay

	n=nvcx
	m=nvcy	
		
	write(*,*)n,m





!     calculo de alfax y alfay

      do i=1,560
		do j=1,3000

			alfax(i,j)=k(i,j)*deltay/deltax
			alfay(i,j)=k(i,j)*deltax/deltay

		enddo
	enddo

**********************************************************************************************
!     Temperaturas iniciales y de frontera supuesta
**********************************************************************************************


	do i=1, n
		do j= 1, m

		 tc(i,j) =50.0

		enddo
	enddo


**********************************************************************************************
**********************************************************************************************


	cont=1

15	write(*,*) cont

!	write(*,*)'para los coeficientes a'

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
!	write(*,*)'para los coeficientes b'
	do j=1, m
		do i=1, n-1
			b(i,j)= alfax(i,j)
		enddo
	enddo

	do j=1,m
	b(n,j)=0.0
	enddo

**********************************************************************************************
!	write(*,*)'para los coeficientes c'


	do j=1, m
		do i=2, n
			c(i,j)= alfax(i,j)
		enddo
	enddo

	do j=1,m
	c(1,j)=0.0
	enddo


**********************************************************************************************
!	write(*,*)'para los coeficientes d'


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

!	calculo de la matriz tidiagonal TDMA 

!     calculo de p (1,j)
	do j=1, m
	p(1,j)= b(1,j)/a(1,j)
!	write(*,*)1,j,p(1,j)
	enddo

!     calculo de p(i,j) desde i=2 a n y j=1 a m 
	do i=2, n
	do j= 1, m
	p(i,j)= b(i,j)/(a(i,j)-(c(i,j)*p(i-1,j)))
!	write(*,*)i,j,p(i,j)
	enddo
	enddo

!     calculo de q(1,1)
!	write(*,*)'q'
	do j=1,m
	q(1,j)= d(1,j)/a(1,j)
!	write(*,*)1,j,q(1,j)
	enddo

!     calculo de q(i,j) desde 2 a m-1 y 2 a n-1
	do i= 2, n
	do j=1,m
	q(i,j)= (d(i,j)+(c(i,j)*q(i-1,j)))/(a(i,j)-(c(i,j)*p(i-1,j)))
!	write(*,*)i,j,q(i,j)
	enddo
	enddo






	do j=1,m
	Tnew(n,j)= q(n,j)
!	write(*,*)n,j,Tnew(n,j)
	enddo

!     Calculo de las T(i)
	do j=1, m
	do i=n-1, 1, -1
	Tnew(i,j)= (p(i,j)*tnew(i+1,j))+q(i,j)
!	write(*,*)i,j, tnew(i,j)
	enddo 
	enddo


**********************************************************************************************
!    Criterio de convergencia

	do j=1,m
		do i=1, n
		deltatem(i,j) = abs(tnew(i,j) - Tc(i,j))
!		write(*,*)'deltatem',i,j,deltatem(i,j)
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
	
!	write(*,*) 'xmaxdeltate', xMAXDELTATE


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
	
	DO j=1,m
			do i=1,1

!				print *,'i, j, tnew(i,j):',i,j,tnew(i,j)
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

!				print *,'i, j, tnew(i,j):',i,j,tnew(i,j)
				write (3,105) i, j, tnew (i,j)

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
