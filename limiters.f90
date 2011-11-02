module Limiters

 use Parameters, only : KND
 
 implicit none

  private
  public FluxLimiter, Limiter, limitertype, limparam
  
  integer,parameter:: minmodlim=1,extminmodlim=2,gammalim=3,vanalbadalim=4,vanleerlim=5,superbeelim=6
  integer,save   :: limitertype
  real(KND),save :: limparam

 contains

  elemental real(KND) function FluxLimiter(r)
    real(KND),intent(in):: r
    real(KND) K

    K=(1+2._KND*r)/3._KND  !assuming kappa=1/3 scheme (Koren, 1993; Hundsdorfer, 1994)
    if (limparam>0) then
     FluxLimiter=max(0._KND,min(2._KND*r,min(limparam,K)))
    elseif (limparam<0) then
     FluxLimiter=K
    else
     FluxLimiter=0
    endif
  endfunction FluxLimiter


  elemental real(KND) function Limiter(a,b,c)
    real(KND),intent(in):: a,b
    real(KND),intent(in),optional:: c
    real(KND) R
    real(KND),parameter:: epsil=0.00001_KND

    if (limitertype==minmodlim) then
     Limiter=MinMod(a,b)
    elseif (limitertype==extminmodlim) then
     if (.not.present(c)) then
       Limiter=MinMod(limparam*a,limparam*b,(a+b)/2._KND)
     else
       Limiter=MinMod(limparam*a,limparam*b,c)
     endif
    elseif (limitertype==gammalim) then
     if (abs(b)>epsil.and.abs(a)>epsil) then
      R=a/b
      Limiter=b*GAMMAL(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limitertype==vanalbadalim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*VanAlbada(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limitertype==vanleerlim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*VanLeer(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limitertype==superbeelim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*SuperBee(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limitertype==0) then
     Limiter=0
    else
     Limiter=(a+b)/2._KND
    endif
  endfunction Limiter

  elemental real(KND) function Heaviside(x)
    real(KND),intent(in):: x
    Heaviside=(1._KND+SIGN(1._KND,x))/2._KND
  endfunction Heaviside

  elemental real(KND) function GAMMAL(R)
    real(KND),intent(in):: R
    GAMMAL=((1-limparam)/limparam)*R*(Heaviside(R) - Heaviside(r-(limparam/(1-limparam)))) + Heaviside(r-(limparam/(1-limparam)))
  endfunction GAMMAL

  elemental real(KND) function VanAlbada(R)
    real(KND),intent(in):: R
    VanAlbada=(R+R*R)/(1+R*R)
  endfunction VanAlbada

  elemental real(KND) function VanLeer(R)
    real(KND),intent(in):: R
    VanLeer=(R+abs(R))/(1+abs(R))
  endfunction VanLeer

  elemental real(KND) function SuperBee(R)
    real(KND),intent(in):: R
    SuperBee=MAX(0._KND,MAX(MIN(2*R,1._KND),MIN(R,2._KND)))
  endfunction SuperBee

  elemental real(KND) function MinMod(a,b,c)
    real(KND),intent(in):: a,b
    real(KND),intent(in),optional:: c

    if (.not. present(c)) then

     MinMod=(SIGN(1._KND,a)+SIGN(1._KND,b))*MIN(ABS(a),ABS(b))/2._KND

    else

     if ((a>0).and.(b>0).and.(c>0)) then
         MinMod=MIN(a,b,c)
     elseif ((a<0).and.(b<0).and.(c<0)) then
         MinMod=MAX(a,b,c)
     else
         MinMod=0
     endif
    endif
  endfunction MinMod

end module Limiters