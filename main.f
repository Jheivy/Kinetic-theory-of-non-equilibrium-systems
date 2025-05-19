      program kinetic_simulation
      implicit none
      real*8 vx,vy,vz,v_mod,T,ta,saved
      integer N,NMC,i,pos,nbin,saved_v4,saved_H,nbucle
      integer gs,n_v,j
      parameter (N = 100000)  ! Número de velocidades
      parameter (T = 5.d-1)  ! Temperatura del sistema
      parameter (v_mod = 1d0) ! Radio de la esfera de velocidades
      parameter (NMC=200)
      parameter (nbin=1000)
      parameter (n_v=1) ! Número de dimensiones del espacio de velocidades
      parameter (saved_v4=1) ! Si no quiero guardar los datos saved=0 y si quiero guardar los datos saved=1
      parameter (saved_H=0)
      parameter (gs=1) ! 0 Distribución aleatoria con módulo 1; 1 Distribución aleatoria en cubo
      parameter (nbucle=10) ! Número de simulaciones
      dimension vx(N),vy(N),vz(N)
      do j=1,200
      do i=1,nbucle ! Este bucle es para hacer varias simulaciones e ir guardandolas
            call init_random_seed()
c           Generar velocidades iniciales
            call GENERATE_INITIAL_VELOCITIES(N,T,v_mod,vx,vy,vz,gs)
c            call SAVE_VALUES(N, vx, vy, vz, "Dataframes/v0_n.dat") ! Guardar velocidades iniciales
            call DSMC(N,j,nbin,i,vx,vy,vz,ta,saved_v4,saved_H,n_v,gs)
c           ta=ta/(NMC*N)*100
c           print*,"Tasa de aceptacion: ", ta, "%"
c            call SAVE_VALUES(N, vx, vy, vz, "Dataframes/v1_n.dat") ! Guardar velocidades finales
            call VELOCITIES_HISTOGRAM(N,i,vx,vy,vz)
            print*, "Generando datos para la simulacion",i,"/",nbucle
      end do
      end do
      print*,"--------------------------------------"
      print*,"-------Simulacion finalizada----------"
      print*,"--------------------------------------"

      end program kinetic_simulation

c----------------------------------------------------------------------------------------
c---------------------   Guardar valores en un archivo  --------------------------------- 
c----------------------------------------------------------------------------------------
      subroutine SAVE_VALUES(N, vx, vy, vz, namesaved)
c     Subroutina para guardar los archivos 
c     [INPUT]
c     N (Integer) Número de puntos
c     vx,vy,vz (Real) velocidades
c     namesaved (CHARACTER) Es el nombre del archivo que quieres que se guarde
c     [OUTPUT]
c     Archivo guardado
      implicit none
      integer N,i
      real*8 vx(N), vy(N), vz(N)
      character(*), intent(in) :: namesaved !Si lo pongo de otra forma me da error


      ! Abrir el archivo con el nombre especificado
      open(10, file=trim(namesaved), status="replace")

      ! Escribir datos en el archivo
      do i = 1, N
            write(10, *) vx(i), vy(i), vz(i)
      end do

      ! Cerrar el archivo
      close(10)

      ! Mensaje de confirmación
      print *, "Datos de las velocidades guardados en ", trim(namesaved)

      end subroutine SAVE_VALUES
c----------------------------------------------------------------------------------------
c---------------------   Generar velocidades iniciales  --------------------------------- 
c----------------------------------------------------------------------------------------
      subroutine GENERATE_INITIAL_VELOCITIES(N,T,v_mod,vx,vy,vz,gs)
c     Subroutina para generar las velocidades iniciales
c     [INPUT]
c     N (integer): Número de partículas
c     V_mod (real): Módulo de la velocidad
c     T (real): Temperatura
c     [OUTPUT]
c     vx (real): Velocidad en x
c     vy (real): Velocidad en y
c     vz (real): Velocidad en z
      implicit none
      real*8 vx,vy,vz,v_mod,T
      real*8 rx, ry, rz,PI,norm
      real*8 vcm_x,vcm_y,vcm_z
      real*8 theta,phi
      integer N,i,gs
      dimension vx(N),vy(N),vz(N)
      dimension rx(N),ry(N),rz(N),norm(N)
      PARAMETER (PI = ATAN(1.D0) * 4.D0)


      if (gs.eq.0)then
            do i=1,N
                  call random_number(theta)
                  theta = acos(2.0d0 * theta - 1.0d0)
                  call random_number(phi)
                  phi = 2.0d0 *PI* phi 

                  vx(i) = SIN(theta) * cos(phi)
                  vy(i) = SIN(theta) * sin(phi)
                  vz(i) = cos(theta)
            end do
      end if
      if (gs.eq.1) then
            call random_number(vx)
            call random_number(vy)
            call random_number(vz)

            ! Transformar para obtener valores en [-1, 1]
            vx = 2.0d0 * vx - 1.0d0
            vy = 2.0d0 * vy - 1.0d0
            vz = 2.0d0 * vz - 1.0d0
      end if

      ! Imponer momento igual a 0  
      call ONE_MOMENT(N,vx,vy,vz)

      ! Escalar las velocidades con la temperatura del sistema
      call SST(N,T,vx,vy,vz)
      end subroutine GENERATE_INITIAL_VELOCITIES


c----------------------------------------------------------------------------------------
c--------------------   Momento total del sistema igual a 0  ----------------------------
c----------------------------------------------------------------------------------------
      subroutine ONE_MOMENT(N,vx,vy,vz)
c     Subroutina para poner momento igual a 0
c     [INPUT]
c     N (integer): Número de partículas
c     vz,vy,vz (real): Velocidades
c     [OUTPUT]
c     vx,vy,vz (real) Velocidades
      implicit none
      real*8 vx,vy,vz
      real*8 vcm_x,vcm_y,vcm_z
      integer N,i
      dimension vx(N),vy(N),vz(N)      

      ! Imponer momento igual a 0  
      Vcm_x = sum(vx) / N
      Vcm_y = sum(vy) / N
      Vcm_z = sum(vz) / N

      do i = 1, N
            vx(i) = vx(i) - Vcm_x
            vy(i) = vy(i) - Vcm_y
            vz(i) = vz(i) - Vcm_z
      end do
c      print '(A, F6.2)', "Momento total del sistema", sum(vx) + sum(vy) + sum(vz)
      end subroutine ONE_MOMENT
c----------------------------------------------------------------------------------------
c--------------------   Escalar las velocidades con temperatura (SST)  ------------------
c----------------------------------------------------------------------------------------

      subroutine SST(N,T,vx,vy,vz)
c     Subroutina para escalar las velocidades en función de la temperatura del sistema
c     [INPUT]
c     N (integer): Número de partículas
c     T (real): Temperatura
c     vz,vy,vz (real): Velocidades
c     [OUTPUT]
c     vx,vy,vz (real) Velocidades

c     Utilizar k=1, T=0.5, E=3N/4 y v^2=3/2   
      implicit none
      real*8 vx,vy,vz,sum_v_2,lambda,T
      integer N,i
      dimension vx(N),vy(N),vz(N)


      sum_v_2=sum(vx**2+vy**2+vz**2)
      lambda=sqrt(3.d0*N*T/sum_v_2) !!!!Al imprimir esto al cuadrado me da 3/2
      vx=vx*lambda
      vy=vy*lambda
      vz=vz*lambda

      ! Imprimimos la nueva energía cinética para verificar
c      sum_v_2 = 0d0
c      do i = 1, N
c          sum_v_2 = sum_v_2 + vx(i)**2 + vy(i)**2 + vz(i)**2
c      end do
c      print*, "Energia cinetica despues de escalar:", sum_v_2/(2.d0*N) !!! Da 3/4 que debería dar
c      print '(A, F6.2)', "Momento total del sistema", sum(vx) + sum(vy) + sum(vz) !!! Da 0 como debe ser
      end subroutine SST


c----------------------------------------------------------------------------------------
c-----------------   Metodo Monte Carlo simulación directa (DSMC)  ---------------------- 
c----------------------------------------------------------------------------------------
      subroutine DSMC(N,NMC,nbin,pos,vx,vy,vz,ta,saved_v4,
     & saved_H,n_v,gs)
c     Subroutina para hacer las colisiones 
c     [INPUT]
c     N (integer): Número de partículas
c     NMC (integer): Número de pasos de Monte Carlo
c     pos (integer): Posición en la que se guardará el archivo
c     saved_v4 (integer): Guardar los valores de v4
c     saved_H (integer): Guardar los valores de H
c     nbin (integer): Número de bins del histograma
c     vz,vy,vz (real): Velocidades
c     [OUTPUT]
c     vx,vy,vz (real) Velocidades   
c     ta (real) Tasa de aceptación   
      implicit none
      real*8 vx,vy,vz,norm,vel,hist,Ht,n_v_r
      real*8 sigma_x,sigma_y,sigma_z,rnd1,rnd2,w_max,w,ta,mean_v4
      real*8 v_rel_x,v_rel_y,v_rel_z,rand_acep,value
      integer N,i,p1,NMC,p2,pos,j,counter,saved_v4,nbin,saved_H,val
      integer n_v,gs
      dimension vx(N),vy(N),vz(N),hist(nbin),vel(N)
      character(len=1000) filename,filename2
      character(len=1000) num_str,nv_str

c      n_v_r=dble(n_v)/2.d0
      if (saved_v4.eq.1) then
            write(num_str, '(I0)') pos  ! I0 evita espacios extra
            write(nv_str, '(I0)') n_v  

            ! Construir el nombre de la carpeta y el archivo dinámicamente
            if (gs.eq.0) then
                  filename = "Dataframes/random_dis/mean_"
     &//trim(nv_str)//"/vel"//trim(nv_str)//"_"//trim(num_str)
            else if (gs.eq.1) then
                  filename = "Dataframes/normal_dis/mean_"
     &//trim(nv_str)//"/vel"//trim(nv_str)//"_"//trim(num_str)
            end if
            ! Abrir archivo con un número de unidad de archivo distinto
            open(33, file=trim(filename)//".dat", status="unknown")    
      end if
      if (saved_H.eq.1) then
            write(num_str, '(I0)') pos
            ! Crear el nombre del archivo dinámico
            if (gs.eq.0)then
            filename2 = "Dataframes/random_dis/H/h_"//trim(num_str)
            else if (gs.eq.1)then
            filename2 = "Dataframes/normal_dis/H/h_"//trim(num_str)
            end if
            ! Abrir archivo con un número de unidad de archivo distinto
            open(34, file=trim(filename2)//".dat")    
      end if

      ta=0d0
      w_max=1.d-15
      counter=0
      val=0
      do j=1,NMC
            do i=1,N !Numero de pasos de Monte Carlo
                  norm = 2d0
                  do while (norm.GT.1.d0)
                  call random_number(sigma_x)
                  call random_number(sigma_y)
                  call random_number(sigma_z)

                        ! Transformamos a valores en [-1,1]
                  sigma_x = 2.0d0 * sigma_x - 1.0d0
                  sigma_y = 2.0d0 * sigma_y - 1.0d0
                  sigma_z = 2.0d0 * sigma_z - 1.0d0

                  norm = sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2)
                  end do                

                  sigma_x =  sigma_x / norm
                  sigma_y =  sigma_y / norm
                  sigma_z =  sigma_z / norm 


                  ! Elegir la primera velocidad aleatoria
                  call random_number(rnd1)
                  p1=int(rnd1*N)+1 

                  ! Elegir la segunda partícula asegurando que sea diferente de p1
                  call random_number(rnd2)
                  p2 = int(rnd2 * N) + 1
                  do while (p2.eq.p1) !Para evitar que se haga repetidos
                        call random_number(rnd2)
                        p2 = int(rnd2 * N) + 1
                  end do
                  
                  ! Velocidades relativas
                  v_rel_x=vx(p1)-vx(p2)
                  v_rel_y=vy(p1)-vy(p2)
                  v_rel_z=vz(p1)-vz(p2)


                  w=v_rel_x*sigma_x+v_rel_y*sigma_y+v_rel_z*sigma_z
                  value=w/w_max
                  call random_number(rand_acep)
                  if (rand_acep.LE.value)then
                        vx(p1)=vx(p1)-w*sigma_x
                        vy(p1)=vy(p1)-w*sigma_y
                        vz(p1)=vz(p1)-w*sigma_z

                        vx(p2)=vx(p2)+w*sigma_x
                        vy(p2)=vy(p2)+w*sigma_y
                        vz(p2)=vz(p2)+w*sigma_z
                        ta=ta+1.d0
                  end if

            if (abs(w).GT.w_max)then
                  w_max=abs(w)
            end if
c          Cálculo de <v^4>      
            counter = counter + 1  ! Incrementar contador
            if (saved_v4.eq.1)then ! Si no quiero guardar los datos saved=0 y si quiero guardar los datos saved=1
                  if (mod(counter, N) .eq. 0) then 
                  val=val+1
                  if (n_v.eq.1)then
                        mean_v4=sum(sqrt(vx**2+vy**2+vz**2))/N
                  elseif (n_v.eq.2)then
                        mean_v4=sum((vx**2+vy**2+vz**2))/N
                  elseif (n_v.eq.3)then
                        mean_v4=sum(sqrt(vx**2+vy**2+vz**2)*(
     &vx**2+vy**2+vz**2))/N
                  elseif (n_v.eq.4)then
                        mean_v4=sum((vx**2+vy**2+vz**2)**(2))/N
                  elseif (n_v.eq.5)then
                        mean_v4=sum(sqrt(vx**2+vy**2+vz**2
     &)*(vx**2+vy**2+vz**2)**(2))/N
                  elseif (n_v.eq.6)then
                        mean_v4=sum((vx**2+vy**2+vz**2)**(3))/N
                  end if
                  write(33,*) j,mean_v4 ! counter
                  end if
            end if
c          Cálculo de H     
            if (saved_H.eq.1)then
                  if (mod(counter, N) .eq. 0) then 
                        vel=sqrt(vx**2+vy**2+vz**2)
                        call THEOREM_H(vel,N,nbin,Ht)
                        write(34,*) j,Ht
                  end if
            end if

            end do
      end do
      close(33)
      close(34)
      end subroutine DSMC

c----------------------------------------------------------------------------------------
c----------------------------   histograma teorema H ------------------------------------ 
c----------------------------------------------------------------------------------------
      subroutine THEOREM_H(vel, N, nbin, Ht)
c     [INPUT]
c     vel (real): Módulo de la velocidad
c     N (int): Número de partículas
c     nbin (int): Número de bins del histograma
c     [OUTPUT]
c     Ht (real): Valor del teorema H calculado
c
      implicit none
      integer N, nbin, bin_index, i
      real*8 vel, hist, vi, PI
      real*8 min_vel, max_vel, h
      real*8 Ht
      PARAMETER (PI = ATAN(1.D0) * 4.D0)
      dimension vel(N), hist(nbin)

      ! Inicializar el histograma y el teorema H
      hist = 0.d0
      Ht = 0.d0

      ! Calcular los límites del histograma
      min_vel = minval(vel)
      max_vel = maxval(vel)

      h = (max_vel - min_vel)/nbin

      ! Construir el histograma
      do i = 1, N
            bin_index = int((vel(i) - min_vel) / h) + 1
            if (bin_index > nbin) then
                  bin_index = nbin
            end if
            hist(bin_index) = hist(bin_index) + 1.d0
      end do
      hist=hist/(N*h) ! Necesario para normalizar el histograma (de esta forma la suma de histogramas da N)
c      print*, sum(hist*h) ! Comprobar que la suma de los histogramas da N

      ! Calcular el teorema H
      do i = 1, nbin
            vi = min_vel+(i)*h !Para que esté en el centro del bin
            if (hist(i).GT.0.d0) then
            Ht = Ht+hist(i)*log(hist(i)/(4.d0*PI*vi**2))*h
            end if
      end do


      end subroutine THEOREM_H



c----------------------------------------------------------------------------------------
c---------------------------------   Histograma ---------------------------------------- 
c----------------------------------------------------------------------------------------
      subroutine VELOCITIES_HISTOGRAM(N, NMC, vx, vy, vz)
c     Subroutina para guardar el módulo de velocidades para el histograma
c     [INPUT]
c     N (integer): Número de partículas
c     NMC (integer): Número de pasos de Monte Carlo
c     vz,vy,vz (real): Velocidades
c     [OUTPUT]
c     vx,vy,vz (real) Velocidades            
      implicit none
      real*8 vx, vy, vz, v_mod
      integer N, i, NMC
      character(len=100) filename
      character(len=100) num_str
      dimension vx(N), vy(N), vz(N)

      ! Convertir NMC a string
      write(num_str, '(I5)') NMC
      ! Crear el nombre del archivo dinámico
      filename = "Dataframes/normal_dis/hist/v_hist_"//trim(num_str)
      print*, filename
      ! Abrir archivo con un número de unidad de archivo distinto
      open(11, file=trim(filename)//".dat")

      ! Guardar datos
      do i=1, N
            v_mod = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
            write(11, *) v_mod
      end do

      close(11)

      end subroutine VELOCITIES_HISTOGRAM

c---------------------------------------------------------------------------
c------------------- RANDOM SUBROUTINA ------------------------------------
c --------------------------------------------------------------------------
      subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms

      call random_seed(size = n)
      allocate(seed(n))
! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream",
     +  form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
        read(un) seed
        close(un)
      else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
        call system_clock(count)
        if (count /= 0) then
          t = transfer(count, t)
        else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 
     -         + dt(2) * 31_8 * 24 * 60 * 60 * 1000 
     -         + dt(3) * 24 * 60 * 60 * 60 * 1000 
     -         + dt(5) * 60 * 60 * 1000 
     -         + dt(6) * 60 * 1000 + dt(7) * 1000 
     -         + dt(8)
          t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n.ge.3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
        else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed