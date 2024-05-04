         parameter(m=393,n=153)
         real x(m,n),y(m,n),cori(m,n),dep(m,n),dep_swan(m,n)

          m1=m
          n1=n
          m2=m-1
          n2=n-1

         
          open(1,file='xx_curv.txt')
           do j=1,n1
            read(1,*)(x(i,j),i=1,m1)
           enddo
          close(1)

          open(1,file='yy_curv.txt')
           do j=1,n1
            read(1,*)(y(i,j),i=1,m1)
           enddo
          close(1)

          open(1,file='dep_curv.txt')
           do j=1,n1
            read(1,*)(dep(i,j),i=1,m1)
           enddo
          close(1)

         do j=1,n1
           do i=1,m1
             dep(i,j)=-dep(i,j)
!             if(dep(i,j).gt.100)dep(i,j)=100.0
           enddo
         enddo


         depmax=100.
         depmin=0.1
         do j=1,n1
           do i=1,m1
             dep_swan(i,j)=dep(i,j)
           enddo
         enddo

         goto 111

! ideal
         depmax=10.
         depmin=0.1
         do j=1,n1
           do i=1,m1
             cori(i,j)=37.0
             dep_swan(i,j)=dep(i,j)
           enddo
         enddo
! end ideal

111     continue

! for swan
       do j=1,n1
        do i=1,m1
          if(dep_swan(i,j).lt.0.1) dep_swan(i,j)=0.1
        enddo
        enddo

        open(1,file='x_circ.txt')
         do j=1,n2
           write(1,100)(x(i,j),i=1,m2)
         enddo
        close(1)

        open(1,file='y_circ.txt')
         do j=1,n2
           write(1,100)(y(i,j),i=1,m2)
         enddo
        close(1)

        open(1,file='cori.txt')
         do j=1,n2
           write(1,100)(cori(i,j),i=1,m2)
         enddo
        close(1)

        open(1,file='dep_circ.txt')
         do j=1,n2
           write(1,100)(dep(i,j),i=1,m2)
         enddo
        close(1)

! for swan
         open(1,file='xxyy_swan.txt')
         do j=1,n2
         write(1,100)(x(i,j),i=1,m2)
         enddo
         do j=1,n2
         write(1,100)(y(i,j),i=1,m2)
         enddo
         close(1)

         open(1,file='dep_swan.txt')
         do j=1,n2
         write(1,100)(dep_swan(i,j),i=1,m2)
         enddo
         close(1)

100     format(3000f12.3)

        end
