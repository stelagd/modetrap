c*************************************************************************
c 
c calculate mode trapping for the simple case of oscillations on a uniform
c string on which "beads" have been placed at x1, x2, ...
c
c*************************************************************************

c      subroutine modetrap(npmax,np,n1,n2,xpert,epsi,width,eig)
      subroutine modetrap(np,n1,n2,xpert,epsi,width,eig)
      implicit none
      integer i,j,k,n,nmax,imax,ndim,ieig,npmax
      parameter (nmax=1000,ndim=10000)
      integer, intent(in) :: np,n1,n2
      double precision, intent(in) :: xpert(np), epsi(np), width(np) 
      double precision, intent(out) :: eig(n2)
      double precision aamp(nmax),bamp(nmax),delta,lambda,
     1     amat(2,2,nmax),bmat(2,2,nmax),pi,bc,dlam,freq(nmax),
     2     discr(ndim),lamvec(ndim),lammin,lamp(nmax),lamm(nmax),lam1,
     3     lam2,tol,zbrent,phi(nmax),eig2(nmax)
      external zbrent, bcval

      pi=4.*atan(1.0)
      dlam=2.*pi/40.
      i=0
      n=0
      lammin=1.d-02
      lambda=lammin

c Calculate discriminant and count the number of modes, n

      do while (n.lt.n2) 
         i=i+1
         call bcval(lambda,xpert,epsi,width,amat,bmat,bc,np,nmax)
c         write(*,10) lambda,bc
 10   format('lambda bc = ',2f12.6)
         lamvec(i) = lambda
         discr(i) = bc
         if (i.ne.1) then
            if (discr(i)*discr(i-1).le.0.d+00) then
               n=n+1
c               write(*,*) 'Found a mode: ',n,i,lamvec(i)
c   store lambda values which bracket eigenmodes
c               lamm(n)=lamvec(i-1) -  0.01*dlam
c               lamp(n)=lamvec(i)   +  0.01*dlam
               lamm(n)=lamvec(i-1) 
               lamp(n)=lamvec(i)  
            endif
         endif
         if (i.eq.ndim) then
            write(*,40) n2
 40       format('Array boundary exceeded before frequency',i4,' found')
         endif
         lambda=lambda+dlam
      enddo
 30   imax=i

c converge to mode frequencies/periods, i.e., lambda values

      tol=1.d-13
      do n=n1,n2
         lam1=lamm(n)
         lam2=lamp(n)
         eig(n)=zbrent(lam1,lam2,tol,xpert,epsi,width,amat,bmat,
     1      np,nmax)
      enddo

      end

c*************************************************************************

      subroutine bcval(lambda,xpert,epsi,width,amat,bmat,bc,npert,nmax)
      implicit none
      integer nmax,npert,i,j,k,m
      double precision xpert(nmax),epsi(nmax),lambda,li,bc,cvec(2),
     1     amat(2,2,nmax),bmat(2,2,nmax),temp,width(nmax),enorm

      li=lambda*xpert(1)
      enorm=epsi(1)*exp(-(lambda*width(1)/2.)**2)
      amat(1,1,1)= 1.d0 - enorm*sin(li)*cos(li)
      amat(1,2,1)=      - enorm*(cos(li))**2
      amat(2,1,1)=        enorm*(sin(li))**2
      amat(2,2,1)= 1.d0 + enorm*sin(li)*cos(li)

      do i=2,npert
         li=lambda*xpert(i)
         enorm=epsi(i)*exp(-(lambda*width(i)/2.)**2)
         bmat(1,1,i)= 1.d0 - enorm*sin(li)*cos(li)
         bmat(1,2,i)=      - enorm*(cos(li))**2
         bmat(2,1,i)=        enorm*(sin(li))**2
         bmat(2,2,i)= 1.d0 + enorm*sin(li)*cos(li)

c do matrix multiplication: c = b . a
         do j=1,2
            do k=1,2
               temp=0.0
               do m=1,2
                  temp=temp + bmat(j,m,i)*amat(m,k,i-1)
               enddo
               amat(j,k,i)=temp
            enddo
         enddo
      enddo

      cvec(1)=amat(1,1,npert)
      cvec(2)=amat(2,1,npert)

      bc = sin(lambda)*cvec(1) + cos(lambda)*cvec(2)

      return
      end

c*************************************************************************

      FUNCTION zbrent(x1,x2,tol,xpert,epsi,width,amat,bmat,npert,nmax)
      INTEGER ITMAX 
      INTEGER nmax,npert 
      double precision zbrent,tol,x1,x2,EPS
      double precision xpert(nmax),epsi(nmax),amat(2,2,nmax),
     1     bmat(2,2,nmax),bc,width(nmax)
      EXTERNAL bcval
      PARAMETER (ITMAX=100,EPS=1.e-13)
      INTEGER iter
      double precision a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      call bcval(a,xpert,epsi,width,amat,bmat,bc,npert,nmax)
      fa=bc
      call bcval(b,xpert,epsi,width,amat,bmat,bc,npert,nmax)
      fb=bc
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
     *'root must be bracketed for zbrent'
c      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) 
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        call bcval(b,xpert,epsi,width,amat,bmat,bc,npert,nmax)
        fb=bc
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END

c*************************************************************************

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10000)
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

c*************************************************************************

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      double precision x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      double precision a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

c*************************************************************************
