FUNCTION  bessj0(x)
    REAL(8) bessj0,x
    REAL(8) ax,xx,z
    DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,           &
    s1,s2,s3,s4,s5,s6,y
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
    DATA p1,p2,p3,p4,p5/1.d0, -0.1098628627d-2, 0.2734510407d-4,                &
                        -0.2073370639d-5, 0.2093887211d-6/
    DATA q1,q2,q3,q4,q5/-0.1562499995d-1, 0.1430488765d-3, -0.6911147651d-5,    &
                        0.7621095161d-6, -0.934945152d-7/
    DATA r1,r2,r3,r4,r5,r6/57568490574.d0, -13362590354.d0, 651619640.7d0,      &
                            -11214424.18d0, 77392.33017d0, -184.9052456d0/
    DATA s1,s2,s3,s4,s5,s6/57568490411.d0, 1029532985.d0, 9494680.718d0,        &
                        59272.64853d0, 267.8532712d0,1.d0/
    if(abs(x).lt.8.0d0)then
        y=x**2
        bessj0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/   &
                    (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
    else
        ax=abs(x)
        z=8.0d0/ax
        y=z**2
        xx=ax - 0.785398164d0
        bessj0=sqrt(0.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))      &
            -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
!    write (*,*) 'p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6'
!    write (*,*) p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
!    read (*,*)
  return
  END
