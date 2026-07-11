program ondaplana
    !16-2-06
!Programa inicial para comenzar con el TLM
!En una cavidad s�lo uno de los �ndices puede ser cero. Los modos son 101,110,011,111,210,...
!La frecuencia de resonancia es f=(c/2)*SQRT((l/lx)**2+(m/ly)**2+(n/lz)**2) donde
!(l,m,n) es el �ndice del nodo y lx, ly y lz las dimensiones de la cavidad
!	use msimsl
parameter (imax= 100, jmax= 1, kmax= 1) !Dimensiones m�ximas de la red
Parameter (nits=10000) !N�mero de iteraciones temporales
real V(18,imax,jmax,kmax), line(18)
real mur !Permeabilidad magn�tica relativa
integer Ox1,Oy1,Oz1,Ox2,Oy2,Oz2,Ox3,Oy3,Oz3 !Puntos salida
real Dx,Dy,Dz !Dimensiones en las que dividimos el espacio
double precision :: epsilon, mu


!	Funciones campo el�ctrico y magn�tico
Ex(i,j,k)=2.0*(V(1,i,j,k)+V(2,i,j,k)+V(9,i,j,k)+V(12,i,j,k)+Yx*V(13,i,j,k))/(Dx*(4.0+Yx+Gx))
Ey(i,j,k)=2.0*(V(3,i,j,k)+V(4,i,j,k)+V(8,i,j,k)+V(11,i,j,k)+Yy*V(14,i,j,k))/(Dy*(4.0+Yy+Gy))
Ez(i,j,k)=2.0*(V(5,i,j,k)+V(6,i,j,k)+V(7,i,j,k)+ V(10,i,j,k)+Yz*V(15,i,j,k))/(Dz*(4.0+Yz+Gz))
Hx(i,j,k)=2.0*(-V(4,i,j,k)+V(5,i,j,k)-V(7,i,j,k)+V(8 ,i,j,k)+V(16,i,j,k))/(Zo*Dx*(4.0+Zx))
Hy(i,j,k)=2.0*( V(2,i,j,k)-V(6,i,j,k)-V(9,i,j,k)+V(10,i,j,k)+V(17,i,j,k))/(Zo*Dy*(4.+Zy))
Hz(i,j,k)=2.0*(-V(1,i,j,k)+V(3,i,j,k)-V(11,i,j,k)+V(12,i,j,k)+V(18,i,j,k))/(Zo*Dz*(4.+Zz))


!Constantes el�ctricas
pi=4.*atan(1.0)
epsilon = 8.854187817E-12 ! Permitividad del medio (por ejemplo, para vacío)
mu = 1.2566370614E-6      ! Permeabilidad del medio (por ejemplo, para vacío)
Zo=120.0*pi!=sqrt(mu/epsilon)
c=2.99792458E8
Mur=1.0

!Par�metros del nudo
Dx=1
Dy=Dx
Dz=Dx
Dt=Dx/(2*c)

!Coeficientes de reflexion
gammay=-1. !simetria
gammaz=1.  !simetria
gammax=-1. !Conductor

!Punto de excitacion
Inx=2
Iny=1
Inz=1

!Puntos de lectura de campo electrico y/o magnetico
nouts=3
Ox1=21
Oy1=1
Oz1=1

Ox2=11
Oy2=1
Oz2=1

Ox3=16
Oy3=1
Oz3=1


!Abrimos los ficheros
open(1,file='salidas.txt')


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


!Excitacion con delta de dirac

!do k=1,kmax
!do j=1,jmax
!V(3,Inx,j,k)=100
!end do
!end do


!Estamos en el vacio con lo que las impedancias, admitancias y conductancias son 0
Zx=0.
Zy=0.
Zz=0.
Yx=0.
Yy=0.
Yz=0.

Gx=0.
Gy=0.
Gz=0.

!Determinacion de los elementos de la matriz de dispersion
Saxy  =-(Gx+Yx)/(2.*(Gx+Yx+4.))+Zy/(2.*(Zy+4.))
Saxz  =-(Gx+Yx)/(2.*(Gx+Yx+4.))+Zz/(2.*(Zz+4.))
Sayx  =-(Gy+Yy)/(2.*(Gy+Yy+4.))+Zx/(2.*(Zx+4.))
Sayz  =-(Gy+Yy)/(2.*(Gy+Yy+4.))+Zz/(2.*(Zz+4.))
Sazx  =-(Gz+Yz)/(2.*(Gz+Yz+4.))+Zx/(2.*(Zx+4.))
Sazy  =-(Gz+Yz)/(2.*(Gz+Yz+4.))+Zy/(2.*(Zy+4.))

Sbx   = 4./(2.*(Gx+Yx+4.))  !Son iguales a los e
Sby   = 4./(2.*(Gy+Yy+4.))
Sbz   = 4./(2.*(Gz+Yz+4.))

Scxy  =-(Gx+Yx)/(2.*(Gx+Yx+4.))-Zy/(2.*(Zy+4.))
Scxz  =-(Gx+Yx)/(2.*(Gx+Yx+4.))-Zz/(2.*(Zz+4.))
Scyx  =-(Gy+Yy)/(2.*(Gy+Yy+4.))-Zx/(2.*(Zx+4.))
Scyz  =-(Gy+Yy)/(2.*(Gy+Yy+4.))-Zz/(2.*(Zz+4.))
Sczx  =-(Gz+Yz)/(2.*(Gz+Yz+4.))-Zx/(2.*(Zx+4.))
Sczy  =-(Gz+Yz)/(2.*(Gz+Yz+4.))-Zy/(2.*(Zy+4.))

Sdx   = 4./(2.*(Zx+4.))   !Son iguales a los i
Sdy   = 4./(2.*(Zy+4.))
Sdz   = 4./(2.*(Zz+4.))

Shx   = -(Gx-Yx+4.)/(Gx+Yx+4.)
Shy   = -(Gy-Yy+4.)/(Gy+Yy+4.)
Shz   = -(Gz-Yz+4.)/(Gz+Yz+4.)

Sjx   = (4.-Zx)/(4.+Zx)
Sjy   = (4.-Zy)/(4.+Zy)
Sjz   = (4.-Zz)/(4.+Zz)

Sfx   =Zx*Sdx
Sfy   =Zy*Sdy
Sfz   =Zz*Sdz

Sgx   =Yx*Sbx
Sgy   =Yy*Sby
Sgz   =Yz*Sbz


! Parametros del pulso gaussiano

e0=1. !amplitud del pulso
hmin=10*Dx !lambda minimo que podemos describir con Dx segun lo estudiado en diferencias finitas

fmax=c/hmin

g=fmax/2 !Un pulso gaussiano excita desde la frec 0 a f_exc=g/2 siendo g la anchura del pulso 
tm=2.6/g
g2=g**2

open(3,file='datos.txt')
write(3,*) 'Nits=',nits
write(3,*) 'imax, jmax, kmax=',imax,jmax,kmax
write(3,*) 'Dx=',dx,' dt=',dt
write(3,*) 'Lambda_min=',hmin, 'Fmax=',fmax
write(3,*) 'g=',g,' tm=', tm,'=',tm/dt,'Dt' !tm es el retardo que habra hasta que llegue el pico del pulso
write(3,*) 'Primer maximo:', tm/dt+(Ox1-Inx)*2., 'Segundo maximo:', tm/dt+(Ox2-Inx)*2, 'Tercer maximo:', tm/dt+(Ox3-Inx)*2



!Comienzo de las iteraciones temporales
do n=1,nits



    !BLOQUE DE RECOPILACION DE DATOS

!Para ver que el bucle se ha hecho bien
if (mod(N, 50) == 0) then
print *, N
end if

    
!Se guardan SALIDAS en ficheros externos. 
!En este caso Ey para tres puntos distintos
!Out1 = Ey(Ox1, Oy1, Oz1)
!Out2 = Ey(Ox2, Oy2, Oz2)
!Out3 = Ey(Ox3, Oy3, Oz3)
Out1 = sqrt(Hx(Ox1, Oy1, Oz1)**2+Hy(Ox1, Oy1, Oz1)**2+Hz(Ox1, Oy1, Oz1)**2)
Out2 = sqrt(Hx(Ox2, Oy2, Oz2)**2+Hy(Ox2, Oy2, Oz2)**2+Hz(Ox2, Oy2, Oz2)**2)
Out3 = sqrt(Hx(Ox3, Oy3, Oz3)**2+Hy(Ox3, Oy3, Oz3)**2+Hz(Ox3, Oy3, Oz3)**2)


write (1,fmt='(1x,i7,4(1x,e14.6))') n,n*Dt,Out1,Out2,Out3
!fin almacenamiento de salidas


    !BLOQUE ALIMENTACION

t=Dt*n
tr=t-tm
coef=g2*tr**2
exc=e0*exp(-coef)  !pulso gaussiano: V=e_0*e^(-g^2*(t-tm)^2)
v3exc=exc*Dx

!Añadimos la alimentacion sumandola 
do k=1,kmax
do j=1,jmax
V(3,Inx,j,k)=V(3,Inx,j,k)+v3exc
end do
end do


    !BLOQUE SCATTERING

!nVr=S*nVi

do  k=1,kmax
do  j=1,jmax
do  i=1,imax

!Guardamos las tensiones incidentes en la matriz LINE
do l=1,18
LINE(l)=V(l,i,j,k)
end do

!Multiplicacion de la matris S por las tensiones incidentes
V(1,i,j,k)=Saxz*LINE(1)+Sbx*LINE(2)+Sdz*LINE(3)+Sbx*LINE(9)-Sdz*LINE(11)+Scxz*LINE(12)+Sgx*LINE(13)+Sdz*LINE(18)

V(2,i,j,k)=Sbx*LINE(1)+Saxy*LINE(2)+Sdy*LINE(6)+Scxy*LINE(9)-Sdy*LINE(10)+Sbx*LINE(12)+Sgx*LINE(13)-Sdy*LINE(17)

V(3,i,j,k)=Sdz*LINE(1)+Sayz*LINE(3)+Sby*LINE(4)+Sby*LINE(8)+Scyz*LINE(11)-Sdz*LINE(12)+Sgy*LINE(14)-Sdz*LINE(18)

V(4,i,j,k)=Sby*LINE(3)+Sayx*LINE(4)+Sdx*LINE(5)-Sdx*LINE(7)+Scyx*LINE(8)+Sby*LINE(11)+Sgy*LINE(14)+Sdx*LINE(16)

V(5,i,j,k)=Sdx*LINE(4)+Sazx*LINE(5)+Sbz*LINE(6)+Sczx*LINE(7)-Sdx*LINE(8)+Sbz*LINE(10)+Sgz*LINE(15)-Sdx*LINE(16)

V(6,i,j,k)=Sdy*LINE(2)+Sbz*LINE(5)+Sazy*LINE(6)+Sbz*LINE(7)-Sdy*LINE(9)+Sczy*LINE(10)+Sgz*LINE(15)+Sdy*LINE(17)

V(7,i,j,k)=-Sdx*LINE(4)+Sczx*LINE(5)+Sbz*LINE(6)+Sazx*LINE(7)+Sdx*LINE(8)+Sbz*LINE(10)+Sgz*LINE(15)+Sdx*LINE(16)

V(8,i,j,k)=Sby*LINE(3)+Scyx*LINE(4)-Sdx*LINE(5)+Sdx*LINE(7)+Sayx*LINE(8)+Sby*LINE(11)+Sgy*LINE(14)-Sdx*LINE(16)

V(9,i,j,k)=Sbx*LINE(1)+Scxy*LINE(2)-Sdy*LINE(6)+Saxy*LINE(9)+Sdy*LINE(10)+Sbx*LINE(12)+Sgx*LINE(13)+Sdy*LINE(17)

V(10,i,j,k)=-Sdy*LINE(2)+Sbz*LINE(5)+Sczy*LINE(6)+Sbz*LINE(7)+Sdy*LINE(9)+Sazy*LINE(10)+Sgz*LINE(15)-Sdy*LINE(17)

V(11,i,j,k)=-Sdz*LINE(1)+Scyz*LINE(3)+Sby*LINE(4)+Sby*LINE(8)+Sayz*LINE(11)+Sdz*LINE(12)+Sgy*LINE(14)+Sdz*LINE(18)

V(12,i,j,k)=Scxz*LINE(1)+Sbx*LINE(2)-Sdz*LINE(3)+Sbx*LINE(9)+Sdz*LINE(11)+Saxz*LINE(12)+Sgx*LINE(13)-Sdz*LINE(18)

V(13,i,j,k)=Sbx*(LINE(1)+LINE(2)+LINE(9)+LINE(12))+Shx*LINE(13)

V(14,i,j,k)=Sby*(LINE(3)+LINE(4)+LINE(8)+LINE(11))+Shy*LINE(14)

V(15,i,j,k)=Sbz*(LINE(5)+LINE(6)+LINE(7)+LINE(10))+Shz*LINE(15)

V(16,i,j,k)=Sfx*(LINE(4)-LINE(5)+LINE(7)-LINE(8)) +Sjx*LINE(16)

V(17,i,j,k)=Sfy*(-LINE(2)+LINE(6)+LINE(9)-LINE(10))+Sjy*LINE(17)

V(18,i,j,k)=Sfz*(LINE(1)-LINE(3)+LINE(11)-LINE(12))+Sjz*LINE(18)

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

!Condiciones de contorno
!Conductor en xmax y xmin: Gammax=-1.
!Simetría en ejes y,z
!Para alimentación de Ey:   Gammay=-1., Gammaz=1.


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



