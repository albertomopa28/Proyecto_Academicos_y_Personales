program ondaplana
!Programa inicial para comenzar con el TLM
!En una cavidad solo uno de los indices puede ser cero. Los modos son 101,110,011,111,210,...
!La frecuencia de resonancia es f=(c/2)*SQRT((l/lx)**2+(m/ly)**2+(n/lz)**2) donde
!(l,m,n) es el indice del nodo y lx, ly y lz las dimensiones de la cavidad

    !VARIABLES DE ENTRADA
parameter (imax= 500, jmax= 1, kmax= 1)         !Dimensiones maximas de la red
Parameter (nits=2000)                            !Numero de iteraciones temporales
real ::  epsilon=8.8541878128E-12, mu=1.25663706144E-6, c=2.99792458E8  !Ctes electromagneticas
real :: mur1=1., epsilon_r1=1., mur2=9., epsilon_r2=1.                  !Epsilon y mu relativas de los dos medios
real :: Dl=0.1                                   !Dimensiones en las que dividimos el espacio y dt
real :: gammax=0, gammay, gammaz                !Coeficientes de reflexion alimenta V6 y se trata de modo TM
real :: e0=1.                                    !amplitud del pulso
integer :: Inx=200, Iny=1, Inz=1                 !Punto de excitacion
integer :: i_inter=250                          !Interseccion de ambos medios
!Puntos de salida:
integer :: Ox1=3, Oy1=1, Oz1=1               
integer :: Ox2=4, Oy2=1, Oz2=1               
          
            


    !RESTO DE VARIABLES 
real v2, Dt, Zo, Z1, Z2, Zx, Zy, Zz, Yx, Yy, Yz, Gx, Gy, Gz !unidad de tiempo y impedancias de los stubs
real hmin, fmax, g, g2, tm !Parametros del pulso gaussiano
real V(18,imax,jmax,kmax), line(18)
real E1, H1, E2, H2, E3, H3
integer :: modotetm,lineafuente

!	Funciones campo electrico y magnetico, medio 1:
Ex1(i,j,k)= 2.0*(V(1,i,j,k)+V(2,i,j,k)+V(9,i,j,k)+V(12,i,j,k)+Yx1*V(13,i,j,k))/(Dx*(4.0+Yx1+Gx1))
Ey1(i,j,k)=2.0*(V(3,i,j,k)+V(4,i,j,k)+V(8,i,j,k)+V(11,i,j,k)+Yy1*V(14,i,j,k))/(Dy*(4.0+Yy1+Gy1))
Ez1(i,j,k)=2.0*(V(5,i,j,k)+V(6,i,j,k)+V(7,i,j,k)+ V(10,i,j,k)+Yz1*V(15,i,j,k))/(Dz*(4.0+Yz1+Gz1))
Hx1(i,j,k)=2.0*(-V(4,i,j,k)+V(5,i,j,k)-V(7,i,j,k)+V(8 ,i,j,k)+V(16,i,j,k))/(Zo*Dx*(4.0+Zx1))
Hy1(i,j,k)=2.0*( V(2,i,j,k)-V(6,i,j,k)-V(9,i,j,k)+V(10,i,j,k)+V(17,i,j,k))/(Zo*Dy*(4.+Zy1))
Hz1(i,j,k)=2.0*(-V(1,i,j,k)+V(3,i,j,k)-V(11,i,j,k)+V(12,i,j,k)+V(18,i,j,k))/(Zo*Dz*(4.+Zz1))

!	Funciones campo electrico y magnetico, medio 2:
Ex2(i,j,k)= 2.0*(V(1,i,j,k)+V(2,i,j,k)+V(9,i,j,k)+V(12,i,j,k)+Yx2*V(13,i,j,k))/(Dx*(4.0+Yx2+Gx2))
Ey2(i,j,k)=2.0*(V(3,i,j,k)+V(4,i,j,k)+V(8,i,j,k)+V(11,i,j,k)+Yy2*V(14,i,j,k))/(Dy*(4.0+Yy2+Gy2))
Ez2(i,j,k)=2.0*(V(5,i,j,k)+V(6,i,j,k)+V(7,i,j,k)+ V(10,i,j,k)+Yz2*V(15,i,j,k))/(Dz*(4.0+Yz2+Gz2))
Hx2(i,j,k)=2.0*(-V(4,i,j,k)+V(5,i,j,k)-V(7,i,j,k)+V(8 ,i,j,k)+V(16,i,j,k))/(Zo*Dx*(4.0+Zx2))
Hy2(i,j,k)=2.0*( V(2,i,j,k)-V(6,i,j,k)-V(9,i,j,k)+V(10,i,j,k)+V(17,i,j,k))/(Zo*Dy*(4.+Zy2))
Hz2(i,j,k)=2.0*(-V(1,i,j,k)+V(3,i,j,k)-V(11,i,j,k)+V(12,i,j,k)+V(18,i,j,k))/(Zo*Dz*(4.+Zz2))

Dx=Dl
Dy=Dl
Dz=Dl
pi=4.*atan(1.0)
v2=1/sqrt(mur2*mu*epsilon_r2*epsilon)
Dt=Dx/(2*c)       
Zo=sqrt(mu/epsilon)
Z1=sqrt((mur1*mu)/(epsilon_r1*epsilon))
Z2=sqrt((mur2*mu)/(epsilon_r2*epsilon))

    !Bloque de lectura de datos
open(1,file='salidas_cambio_medio.txt')
open(3,file='datos_cambio_medio.txt')
open(11,file='video.txt')
21 write(*,*) '¿Desea excitar Modo TE o TM (1/2)?'     

read(*,*) modotetm
if (modotetm.eq.1) then !Modo TE
    gammay=-1
    gammaz=1
    lineafuente=3
    write(3,*) modotetm
    write(3,*) 'Excitación +1 sólo en línea 3: Modos TE'
else if (modotetm.eq.2) then !Modo TM
    gammay=1
    gammaz=-1
    lineafuente=6
    write(3,*) modotetm
    write(3,*) 'Excitación +1 sólo en línea 6: Modos TM'
else
    goto 21
end if



!Dimensoines de la cavidad y numero de iteraciones
write(3, *) "Numero de iteraciones: nits=", nits
write(3, *) "Dimensiones: imax=", imax,", jmax=", jmax, "kmax=", kmax
write(3, *) ""
write(3, *) "Constantes electricas:"
write(3, *) "epsilon=", epsilon, " mu=", mu, " c=", c
write(3, *)  " Mur1=", mur1, " epsilon_r1=", epsilon_r1, "Velocidad medio 1: v1=", c
write(3, *)  " Mur2=", mur2, " epsilon_r2=", epsilon_r2, "Velocidad medio 2: v2=", v2
write(3, *)  "Impedancia de la linea: Z0=", Zo
write(3, *) ""
write(3, *) "Parametros del sistema"
write(3, *) "Impedancia medio 1: Z1=", Z1, "Medio 2: Z2=", Z2
write(3, *) "Coeficiente de transmision medio 1 a 2: T=", 2*Z2/(Z1+Z2), "Coeficiente de reflexion: Gamma=", (Z2-Z1)/(Z1+Z2)
write(3, *) "Para H: T=", 2*Z1/(Z1+Z2), "Coeficiente de reflexion: Gamma=", (Z2-Z1)/(Z1+Z2)
write(3, *) ""
write(3, *) "Parametros del nudo:"
write(3, *) "Dx=", Dx, " Dy=", Dy, " Dz=", Dz, " Dt=", Dt
write(3, *) ""
write(3, *) "Coeficientes de reflexion"
write(3, *) "gammax=", gammax, " gammay=", gammay, " gammaz=", gammaz
write(3, *) ""
write(3, *) "Puntos de excitacion y lectura de campo electrico y/o magnetico"
write(3, *) "Inx=", Inx, " , Iny=", Iny, ", Inz=", Inz
write(3, *) "Ox1=", Ox1, ", Oy1=", Oy1, ", Oz1=", Oz1
write(3, *) "Ox2=", Ox2, ", Oy2=", Oy2, ", Oz2=", Oz2
write(3, *) "Interseccion de los medios, i_inter=", i_inter
write(3, *) ""


!Poner a cero la matriz de las tensiones
do k=1,kmax
do j=1,jmax
do i=1,imax
do l=1,18
V(l,i,j,k)= 0.0
end do
end do
end do
end do



!Calculo de los stubs
!Para el medio 1(el vacio):
Zx1=(2.*mur1*Dy*Dz/(c*Dt*Dx))-4.0
Zy1=(2.*mur1*Dx*Dz/(c*Dt*Dy))-4.0
Zz1=(2.*mur1*Dx*Dy/(c*Dt*Dz))-4.0

Yx1=(2.*epsilon_r1*Dy*Dz/(c*Dt*Dx))-4.0
Yy1=(2.*epsilon_r1*Dx*Dz/(c*Dt*Dy))-4.0
Yz1=(2.*epsilon_r1*Dx*Dy/(c*Dt*Dz))-4.0

write(3,*) "Impedancias de los stubs en medio 1 (vacío):"
write(3,*) "Zx1=", Zx1, ", Zy1=", Zy1, ", Zz1=", Zz1
write(3,*) "Yx1=", Yx1, ", Yy1=", Yy1, ", Yz1=", Yz1

Gx1=0.
Gy1=0.
Gz1=0.
write(3,*) "Gx1=", Gx1, ", Gy1=", Gy1, ", Gz1=", Gz1
write(3,*) ""

!Determinacion de los elementos de la matriz de dispersion
Saxy1=-(Gx1+Yx1)/(2.*(Gx1+Yx1+4.))+Zy1/(2.*(Zy1+4.))
Saxz1=-(Gx1+Yx1)/(2.*(Gx1+Yx1+4.))+Zz1/(2.*(Zz1+4.))
Sayx1=-(Gy1+Yy1)/(2.*(Gy1+Yy1+4.))+Zx1/(2.*(Zx1+4.))
Sayz1=-(Gy1+Yy1)/(2.*(Gy1+Yy1+4.))+Zz1/(2.*(Zz1+4.))
Sazx1=-(Gz1+Yz1)/(2.*(Gz1+Yz1+4.))+Zx1/(2.*(Zx1+4.))
Sazy1=-(Gz1+Yz1)/(2.*(Gz1+Yz1+4.))+Zy1/(2.*(Zy1+4.))

Sbx1=4./(2.*(Gx1+Yx1+4.))  !Son iguales a los e
Sby1=4./(2.*(Gy1+Yy1+4.))
Sbz1=4./(2.*(Gz1+Yz1+4.))

Scxy1=-(Gx1+Yx1)/(2.*(Gx1+Yx1+4.))-Zy1/(2.*(Zy1+4.))
Scxz1=-(Gx1+Yx1)/(2.*(Gx1+Yx1+4.))-Zz1/(2.*(Zz1+4.))
Scyx1=-(Gy1+Yy1)/(2.*(Gy1+Yy1+4.))-Zx1/(2.*(Zx1+4.))
Scyz1=-(Gy1+Yy1)/(2.*(Gy1+Yy1+4.))-Zz1/(2.*(Zz1+4.))
Sczx1=-(Gz1+Yz1)/(2.*(Gz1+Yz1+4.))-Zx1/(2.*(Zx1+4.))
Sczy1=-(Gz1+Yz1)/(2.*(Gz1+Yz1+4.))-Zy1/(2.*(Zy1+4.))

Sdx1=4./(2.*(Zx1+4.))   !Son iguales a los i
Sdy1=4./(2.*(Zy1+4.))
Sdz1=4./(2.*(Zz1+4.))

Shx1=-(Gx1-Yx1+4.)/(Gx1+Yx1+4.)
Shy1=-(Gy1-Yy1+4.)/(Gy1+Yy1+4.)
Shz1=-(Gz1-Yz1+4.)/(Gz1+Yz1+4.)

Sjx1=(4.-Zx1)/(4.+Zx1)
Sjy1=(4.-Zy1)/(4.+Zy1)
Sjz1=(4.-Zz1)/(4.+Zz1)

Sfx1=Zx1*Sdx1
Sfy1=Zy1*Sdy1
Sfz1=Zz1*Sdz1

Sgx1=Yx1*Sbx1
Sgy1=Yy1*Sby1
Sgz1=Yz1*Sbz1



!Para el medio 2:
Zx2=(2.*mur2*Dy*Dz/(c*Dt*Dx))-4.0
Zy2=(2.*mur2*Dx*Dz/(c*Dt*Dy))-4.0
Zz2=(2.*mur2*Dx*Dy/(c*Dt*Dz))-4.0

Yx2=(2.*epsilon_r2*Dy*Dz/(c*Dt*Dx))-4.0
Yy2=(2.*epsilon_r2*Dx*Dz/(c*Dt*Dy))-4.0
Yz2=(2.*epsilon_r2*Dx*Dy/(c*Dt*Dz))-4.0

write(3,*) "Impedancias de los stubs en medio 2 (vacío):"
write(3,*) "Zx2=", Zx2, ", Zy2=", Zy2, ", Zz2=", Zz2
write(3,*) "Yx2=", Yx2, ", Yy2=", Yy2, ", Yz2=", Yz2

Gx2=0.
Gy2=0.
Gz2=0.
write(3,*) "Gx2=", Gx2, ", Gy2=", Gy2, ", Gz2=", Gz2
write(3,*) ""

!Determinacion de los elementos de la matriz de dispersion
Saxy2=-(Gx2+Yx2)/(2.*(Gx2+Yx2+4.))+Zy2/(2.*(Zy2+4.))
Saxz2=-(Gx2+Yx2)/(2.*(Gx2+Yx2+4.))+Zz2/(2.*(Zz2+4.))
Sayx2=-(Gy2+Yy2)/(2.*(Gy2+Yy2+4.))+Zx2/(2.*(Zx2+4.))
Sayz2=-(Gy2+Yy2)/(2.*(Gy2+Yy2+4.))+Zz2/(2.*(Zz2+4.))
Sazx2=-(Gz2+Yz2)/(2.*(Gz2+Yz2+4.))+Zx2/(2.*(Zx2+4.))
Sazy2=-(Gz2+Yz2)/(2.*(Gz2+Yz2+4.))+Zy2/(2.*(Zy2+4.))

Sbx2=4./(2.*(Gx2+Yx2+4.))  !Son iguales a los e
Sby2=4./(2.*(Gy2+Yy2+4.))
Sbz2=4./(2.*(Gz2+Yz2+4.))

Scxy2=-(Gx2+Yx2)/(2.*(Gx2+Yx2+4.))-Zy2/(2.*(Zy2+4.))
Scxz2=-(Gx2+Yx2)/(2.*(Gx2+Yx2+4.))-Zz2/(2.*(Zz2+4.))
Scyx2=-(Gy2+Yy2)/(2.*(Gy2+Yy2+4.))-Zx2/(2.*(Zx2+4.))
Scyz2=-(Gy2+Yy2)/(2.*(Gy2+Yy2+4.))-Zz2/(2.*(Zz2+4.))
Sczx2=-(Gz2+Yz2)/(2.*(Gz2+Yz2+4.))-Zx2/(2.*(Zx2+4.))
Sczy2=-(Gz2+Yz2)/(2.*(Gz2+Yz2+4.))-Zy2/(2.*(Zy2+4.))

Sdx2=4./(2.*(Zx2+4.))   !Son iguales a los i
Sdy2=4./(2.*(Zy2+4.))
Sdz2=4./(2.*(Zz2+4.))

Shx2=-(Gx2-Yx2+4.)/(Gx2+Yx2+4.)
Shy2=-(Gy2-Yy2+4.)/(Gy2+Yy2+4.)
Shz2=-(Gz2-Yz2+4.)/(Gz2+Yz2+4.)

Sjx2=(4.-Zx2)/(4.+Zx2)
Sjy2=(4.-Zy2)/(4.+Zy2)
Sjz2=(4.-Zz2)/(4.+Zz2)

Sfx2=Zx2*Sdx2
Sfy2=Zy2*Sdy2
Sfz2=Zz2*Sdz2

Sgx2=Yx2*Sbx2
Sgy2=Yy2*Sby2
Sgz2=Yz2*Sbz2




! Parametros del pulso gaussiano

e0=1. !amplitud del pulso
hmin=10*Dx !lambda minimo que podemos describir con Dx segun lo estudiado en diferencias finitas

fmax=v2/hmin

g=fmax*2 !Un pulso gaussiano excita desde la frec 0 a f_exc=g/2 siendo g la anchura del pulso
tm=2.6/g+14.*Dt
g2=g**2

write(3,*) "Parametros del pulso Gaussiano:"
write(3,*) "Apmlitud pulso: e0=", e0, 'Lambda_min=',hmin, 'Fmax=',fmax
write(3,*) 'Un pulso gaussiano excita desde la frec 0 a f_exc=fmax=g/2: g=',g
write(3,*) 'Retardo que habra hasta que llegue el pico del pulso: tm=', tm,'=',tm/dt,'Dt'
!write(3,*) 'Primer maximo:', tm/dt+(Ox1-Inx)*2., 'Segundo maximo:', tm/dt+(Ox2-Inx)*2, 'Tercer maximo:', tm/dt+(Ox3-Inx)*2
write(3,*) ""


!Comienzo de las iteraciones temporales
do n=1,nits



    !BLOQUE DE RECOPILACION DE DATOS

!Para ver que el bucle se ha hecho bien
if (mod(N, 50) == 0) then
print *, N
end if


!Se guardan SALIDAS en ficheros externos.
E1=Ey1(Ox1, Oy1, Oz1)
H1=Hz1(Ox1, Oy1, Oz1)

E2=Ey2(Ox2, Oy2, Oz2)
H2=Hz2(Ox2, Oy2, Oz2)

write(1, '(1x,i6,1x,e14.6,1x,18(e14.6,1x))') n, n*Dt, E1, H1, E2, H2!, E1/H1, E2/H2 !Para no perder precision usar formato cientifico

do i=1, i_inter
    write(11, '(1x,i6,1x,e14.6,1x,18(e14.6,1x))') n, n*Dt,i*Dx, Ey1(i,1,1), Hz1(i,1,1)
end do

do i=i_inter+1, imax
    write(11, '(1x,i6,1x,e14.6,1x,18(e14.6,1x))') n, n*Dt,i*Dx, Ey2(i,1,1), Hz2(i,1,1)
end do

!fin almacenamiento de salidas


    !BLOQUE ALIMENTACION

t=Dt*n
tr=t-tm
coef=g2*tr**2
exc=e0*exp(-coef)  !pulso gaussiano: V=e_0*e^(-g^2*(t-tm)^2)
v3exc=exc*Dx
!Añadimos la alimentacion sumandola
V(lineafuente,Inx,Iny,Inz)=V(3,Inx,Iny,Inz)+v3exc



    !BLOQUE SCATTERING
!nVr=S*nVi

!Bucle para medio 1
do  k=1,kmax
do  j=1,jmax
do  i=1,i_inter

!Guardamos las tensiones incidentes en la matriz LINE
do l=1,18
LINE(l)=V(l,i,j,k)
end do

!Multiplicacion de la matris S1 por las tensiones incidentes

V(1,i,j,k)=Saxz1*LINE(1)+Sbx1*LINE(2)+Sdz1*LINE(3)+Sbx1*LINE(9)-Sdz1*LINE(11)+Scxz1*LINE(12)+Sgx1*LINE(13)+Sdz1*LINE(18)

V(2,i,j,k)=Sbx1*LINE(1)+Saxy1*LINE(2)+Sdy1*LINE(6)+Scxy1*LINE(9)-Sdy1*LINE(10)+Sbx1*LINE(12)+Sgx1*LINE(13)-Sdy1*LINE(17)

V(3,i,j,k)=Sdz1*LINE(1)+Sayz1*LINE(3)+Sby1*LINE(4)+Sby1*LINE(8)+Scyz1*LINE(11)-Sdz1*LINE(12)+Sgy1*LINE(14)-Sdz1*LINE(18)

V(4,i,j,k)=Sby1*LINE(3)+Sayx1*LINE(4)+Sdx1*LINE(5)-Sdx1*LINE(7)+Scyx1*LINE(8)+Sby1*LINE(11)+Sgy1*LINE(14)+Sdx1*LINE(16)

V(5,i,j,k)=Sdx1*LINE(4)+Sazx1*LINE(5)+Sbz1*LINE(6)+Sczx1*LINE(7)-Sdx1*LINE(8)+Sbz1*LINE(10)+Sgz1*LINE(15)-Sdx1*LINE(16)

V(6,i,j,k)=Sdy1*LINE(2)+Sbz1*LINE(5)+Sazy1*LINE(6)+Sbz1*LINE(7)-Sdy1*LINE(9)+Sczy1*LINE(10)+Sgz1*LINE(15)+Sdy1*LINE(17)

V(7,i,j,k)=-Sdx1*LINE(4)+Sczx1*LINE(5)+Sbz1*LINE(6)+Sazx1*LINE(7)+Sdx1*LINE(8)+Sbz1*LINE(10)+Sgz1*LINE(15)+Sdx1*LINE(16)

V(8,i,j,k)=Sby1*LINE(3)+Scyx1*LINE(4)-Sdx1*LINE(5)+Sdx1*LINE(7)+Sayx1*LINE(8)+Sby1*LINE(11)+Sgy1*LINE(14)-Sdx1*LINE(16)

V(9,i,j,k)=Sbx1*LINE(1)+Scxy1*LINE(2)-Sdy1*LINE(6)+Saxy1*LINE(9)+Sdy1*LINE(10)+Sbx1*LINE(12)+Sgx1*LINE(13)+Sdy1*LINE(17)

V(10,i,j,k)=-Sdy1*LINE(2)+Sbz1*LINE(5)+Sczy1*LINE(6)+Sbz1*LINE(7)+Sdy1*LINE(9)+Sazy1*LINE(10)+Sgz1*LINE(15)-Sdy1*LINE(17)

V(11,i,j,k)=-Sdz1*LINE(1)+Scyz1*LINE(3)+Sby1*LINE(4)+Sby1*LINE(8)+Sayz1*LINE(11)+Sdz1*LINE(12)+Sgy1*LINE(14)+Sdz1*LINE(18)

V(12,i,j,k)=Scxz1*LINE(1)+Sbx1*LINE(2)-Sdz1*LINE(3)+Sbx1*LINE(9)+Sdz1*LINE(11)+Saxz1*LINE(12)+Sgx1*LINE(13)-Sdz1*LINE(18)

V(13,i,j,k)=Sbx1*(LINE(1)+LINE(2)+LINE(9)+LINE(12))+Shx1*LINE(13)

V(14,i,j,k)=Sby1*(LINE(3)+LINE(4)+LINE(8)+LINE(11))+Shy1*LINE(14)

V(15,i,j,k)=Sbz1*(LINE(5)+LINE(6)+LINE(7)+LINE(10))+Shz1*LINE(15)

V(16,i,j,k)=Sfx1*(LINE(4)-LINE(5)+LINE(7)-LINE(8))+Sjx1*LINE(16)

V(17,i,j,k)=Sfy1*(-LINE(2)+LINE(6)+LINE(9)-LINE(10))+Sjy1*LINE(17)

V(18,i,j,k)=Sfz1*(LINE(1)-LINE(3)+LINE(11)-LINE(12))+Sjz1*LINE(18)

end do
end do
end do


!Bucle para medio 2
do  k=1,kmax
do  j=1,jmax
do  i=i_inter+1, imax

!Guardamos las tensiones incidentes en la matriz LINE
do l=1,18
LINE(l)=V(l,i,j,k)
end do

!Multiplicacion de la matris S2 por las tensiones incidentes

V(1,i,j,k)=Saxz2*LINE(1)+Sbx2*LINE(2)+Sdz2*LINE(3)+Sbx2*LINE(9)-Sdz2*LINE(11)+Scxz2*LINE(12)+Sgx2*LINE(13)+Sdz2*LINE(18)

V(2,i,j,k)=Sbx2*LINE(1)+Saxy2*LINE(2)+Sdy2*LINE(6)+Scxy2*LINE(9)-Sdy2*LINE(10)+Sbx2*LINE(12)+Sgx2*LINE(13)-Sdy2*LINE(17)

V(3,i,j,k)=Sdz2*LINE(1)+Sayz2*LINE(3)+Sby2*LINE(4)+Sby2*LINE(8)+Scyz2*LINE(11)-Sdz2*LINE(12)+Sgy2*LINE(14)-Sdz2*LINE(18)

V(4,i,j,k)=Sby2*LINE(3)+Sayx2*LINE(4)+Sdx2*LINE(5)-Sdx2*LINE(7)+Scyx2*LINE(8)+Sby2*LINE(11)+Sgy2*LINE(14)+Sdx2*LINE(16)

V(5,i,j,k)=Sdx2*LINE(4)+Sazx2*LINE(5)+Sbz2*LINE(6)+Sczx2*LINE(7)-Sdx2*LINE(8)+Sbz2*LINE(10)+Sgz2*LINE(15)-Sdx2*LINE(16)

V(6,i,j,k)=Sdy2*LINE(2)+Sbz2*LINE(5)+Sazy2*LINE(6)+Sbz2*LINE(7)-Sdy2*LINE(9)+Sczy2*LINE(10)+Sgz2*LINE(15)+Sdy2*LINE(17)

V(7,i,j,k)=-Sdx2*LINE(4)+Sczx2*LINE(5)+Sbz2*LINE(6)+Sazx2*LINE(7)+Sdx2*LINE(8)+Sbz2*LINE(10)+Sgz2*LINE(15)+Sdx2*LINE(16)

V(8,i,j,k)=Sby2*LINE(3)+Scyx2*LINE(4)-Sdx2*LINE(5)+Sdx2*LINE(7)+Sayx2*LINE(8)+Sby2*LINE(11)+Sgy2*LINE(14)-Sdx2*LINE(16)

V(9,i,j,k)=Sbx2*LINE(1)+Scxy2*LINE(2)-Sdy2*LINE(6)+Saxy2*LINE(9)+Sdy2*LINE(10)+Sbx2*LINE(12)+Sgx2*LINE(13)+Sdy2*LINE(17)

V(10,i,j,k)=-Sdy2*LINE(2)+Sbz2*LINE(5)+Sczy2*LINE(6)+Sbz2*LINE(7)+Sdy2*LINE(9)+Sazy2*LINE(10)+Sgz2*LINE(15)-Sdy2*LINE(17)

V(11,i,j,k)=-Sdz2*LINE(1)+Scyz2*LINE(3)+Sby2*LINE(4)+Sby2*LINE(8)+Sayz2*LINE(11)+Sdz2*LINE(12)+Sgy2*LINE(14)+Sdz2*LINE(18)

V(12,i,j,k)=Scxz2*LINE(1)+Sbx2*LINE(2)-Sdz2*LINE(3)+Sbx2*LINE(9)+Sdz2*LINE(11)+Saxz2*LINE(12)+Sgx2*LINE(13)-Sdz2*LINE(18)

V(13,i,j,k)=Sbx2*(LINE(1)+LINE(2)+LINE(9)+LINE(12))+Shx2*LINE(13)

V(14,i,j,k)=Sby2*(LINE(3)+LINE(4)+LINE(8)+LINE(11))+Shy2*LINE(14)

V(15,i,j,k)=Sbz2*(LINE(5)+LINE(6)+LINE(7)+LINE(10))+Shz2*LINE(15)

V(16,i,j,k)=Sfx2*(LINE(4)-LINE(5)+LINE(7)-LINE(8))+Sjx2*LINE(16)

V(17,i,j,k)=Sfy2*(-LINE(2)+LINE(6)+LINE(9)-LINE(10))+Sjy2*LINE(17)

V(18,i,j,k)=Sfz2*(LINE(1)-LINE(3)+LINE(11)-LINE(12))+Sjz2*LINE(18)


end do
end do
end do


    !BLOQUE CONEXIONES (PASO DEL TIEMPO)

!Las tensiones reflejadas se convierten en incidentes en el proximo Dt
!(n+1)Vi=T*nVr

!ConexiOn general direcci0n x. Se dejan sin tocar las lIneas extremas.
DO  k=1,kmax
DO  j=1,jmax
DO  i=1,imax-1
i1= i+1

TEMP1= V(11,i,j,k)
TEMP2= V(10,i,j,k)
V(11,i,j,k) = V(3,i1,j,k)
V(10,i,j,k) = V(6,i1,j,k)
V(3,i1,j,k) = TEMP1
V(6,i1,j,k) = TEMP2
end do
end do
end do

!Conexion general direccion y. Se dejan sin tocar las lineas extremas.
DO k=1,kmax
DO j=1,jmax-1
j1= j+1
DO i=1,imax

TEMP1 = V(12,i,j,k)
TEMP2 = V(7 ,i,j,k)
V(12,i,j,k) = V(1,i,j1,k)
V(7 ,i,j,k) = V(5,i,j1,k)
V(1,i,j1,k) = TEMP1
V(5,i,j1,k) = TEMP2
end do
end do
end do

!Conexion general direccion z. Se dejan sin tocar las lineas extremas.
DO k=1,kmax-1
k1= k+1
DO j=1,jmax
DO i=1,imax

TEMP1 = V(8,i,j,k)
TEMP2 = V(9,i,j,k)
V(8,i,j,k)  = V(4,i,j,k1)
V(9,i,j,k)  = V(2,i,j,k1)
V(4,i,j,k1) = TEMP1
V(2,i,j,k1) = TEMP2
end do
end do
end do


!Fronteras en xmin y xmax
do j=1,jmax
do k=1,kmax
V(3,1,j,k) = gammax*V(3,1,j,k)  !Conductor
V(6,1,j,k) = gammax*V(6,1,j,k)  !Conductor
V(10,imax,j,k) = gammax*V(10,imax,j,k)  !Conductor
V(11,imax,j,k) = gammax*V(11,imax,j,k)  !Conductor
end do
end do


!Fronteras en ymin e ymax. Condicion de simetria
do k = 1,kmax
do i = 1,imax
V(1,i,1,k) = gammay*V(1,i,1,k)
V(5,i,1,k) =gammay*V(5,i,1,k)
V(7,i,jmax,k) = gammay*V( 7,i,jmax,k)
V(12,i,jmax,k) = gammay*V(12,i,jmax,k)
end do
end do


!Fronteras en zmin e zmax. Condicion de simetria
do  j = 1,jmax
do  i = 1,imax
V(2,i,j,1) = gammaz*V(2,i,j,1)
V(4,i,j,1) = gammaz*V(4,i,j,1)
V(9,i,j,kmax) = gammaz*V(9,i,j,kmax)
V(8,i,j,kmax) = gammaz*V(8,i,j,kmax)
end do
end do

!Se convierten las tensiones reflejadas por los stubs en tensiones incidentes
!En los stubs abiertos (13,14,15) el pulso vuelve sin cambiar de signo


Do  k=1,kmax
do  j=1,jmax
do  i=1,imax
V(16,i,j,k) = -V(16,i,j,k)
V(17,i,j,k) = -V(17,i,j,k)
V(18,i,j,k) = -V(18,i,j,k)
end do
end do
end do


End do

close(1)


stop

end