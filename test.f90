      program test
      implicit double precision (a-h,o-z)
      integer iatmcv1, iatmcv2, ncv1grid, ncv2grid, ifmopt
      real*8 cv1grid, cv1f, cv1fpp, cv2grid, cv2f, cv2fpp, &
             cv1min, cv1max, cv1step, cv2min, cv2max, cv2step
      common /CVSPLI/ iatmcv1(2),iatmcv2(2), ncv1grid, ncv2grid, ifmopt
      common /CVSPLF/ cv1grid(500), cv1f(500), cv1fpp(500), &
                      cv2grid(500), cv2f(500),cv2fpp(500), &
                      cv1min, cv1max, cv1step, cv2min, cv2max, cv2step

      integer iatmcv3, iatmcv4, ncv3grid, ncv4grid
      real*8 cv3grid, cv3f, cv3fpp, cv4grid, cv4f, cv4fpp, &
             cv3min, cv3max, cv3step, cv4min, cv4max, cv4step
      common /CVSPLJ/ iatmcv3(2),iatmcv4(2), ncv3grid, ncv4grid
      common /CVSPLG/ cv3grid(500), cv3f(500), cv3fpp(500), &
                      cv4grid(500), cv4f(500),cv4fpp(500), &
                      cv3min, cv3max, cv3step, cv4min, cv4max, cv4step

!------------------------------------------------------------------------------------------------
! cutoff distances to turn on and off the dp-QM/MM interaction between CMS & MM charges; 0320BK20
      real*8 dmpon, cutof
!-------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
! local variables for dp-QM/MM; molecular polarizability on QM; 0312BK20;
!---------------------------------------------------------------------------------------
      integer ist, ied, j, imol
      real*8 alp,chg,epol,xcm,ycm,zcm,efx,efy,efz,mux,muy,muz
      real*8 r1,r2,r3,xx,yy,zz,rr3,rr3xx,rr3yy,rr3zz,dxcm,dycm,dzcm 
      real*8 txww(6),r1p,r3p,rxn
      real*8 tmas,gpfac,dalpdr,ef2

! local variables
      integer i,k
      real*8 yp1, ypn ,x(7500),y(7500),z(7500)
      real*8 xcv1(500),ycv1(500),zcv1(500),xcv2(500),ycv2(500),zcv2(500)
      real*8 dxcv1,dycv1,dzcv1,dxcv2,dycv2,dzcv2

! additional variables for spline under tension (splnut) 1223BK18 
! add common block 1228BK18
      real*8 z_cent,l_scal
      real(8),  parameter :: PI  = 4 * atan (1.0_8)

      common /CVSPBKZ/ z_cent,l_scal

! -- dpol_atm variables, 0711BK20
      integer help,ainnode,aininp,ainout,aitype
      real*8 aisum,aiiw,aib1,ailw,aib2,aiatm1,aiatm2,aipr,aip,aip0, &
             aiarg,aia,aians,aialp,aih,aih2,aiddx,dax,day,daz,dalpdx, &
             dalpdy,dalpdz,aiefx,aiefy,aiefz,aimux,aimuy,aimuz,aidxcm, &
             aidycm,aidzcm,aief2

      common /DPAI/ help,ainnode,aininp,ainout,aitype(512,5)
      common /DPAR/ aisum,aiiw(512,512),aib1(512),ailw(512,512), &
             aib2(512),aiatm1(3),aiatm2(3),aipr,aip(512),aip0(512), &
             aiarg(512),aia(512),aians(512),aialp(512),aih(512), &
             aih2(512),aiddx(512,512),dax,day,daz,dalpdx(512,512), &
             dalpdy(512,512),dalpdz(512,512),aiefx(512),aiefy(512), &
             aiefz(512),aimux(512),aimuy(512),aimuz(512),aidxcm(512), &
             aidycm(512),aidzcm(512),aief2(512)

! -- local variable variabler force, 
      integer natomx
      real*8 dx0(7500),dy0(7500),dz0(7500)
      real*8 dx(7500),dy(7500),dz(7500)

      write (*,'(A)') 'USER: modify forces on CVs based on ann.dat'
      open(unit=90,file='ann.dat',form='formatted',status='old', &
           err=990)
      write (*,'(A)') 'USER: read in CRD FILE'
      open(unit=91,file='crd.txt',form='formatted',status='old', &
           err=990)
      write (*,'(A)') 'USER: read in FORCE FILE'
      open(unit=92,file='force.txt',form='formatted',status='old', &
           err=990)
      go to 995
  990 continue
      write(*,'(A)') 'wrong in opening ann.dat stop.'
      stop
  995 continue

! readin option: ifmopt=1 (cubic spline grid on deltaF; 2-d CVs)
!                ifmopt=2 (FM grid on deltaF and deltaF''; 2-d CVs)
!                ifmopt=3 (cubic spline grid on deltaF; 4-d CVs)
!                ifmopt=4 (FM grid on deltaF and deltaF''; 4-d CVs)
!                ifmopt=14 (atm_ipol(dist), dp-QM/MM)
      read(90,*) ifmopt
      if (ifmopt .eq. 1) then
          write(6,*) 'ann.dat: read cubic spline grid on deltaF.'
      else if (ifmopt .eq. 2) then
          write(6,*) 'ann.dat: read FM grid on deltaF and deltaFPP.'
      else if (ifmopt .eq. 3) then
          write(6,*) 'ann.dat: read cubic spline grid on deltaF;4-d.'
      else if (ifmopt .eq. 4) then
          write(6,*) 'ann.dat: read FM grid on deltaF&deltaFPP; 4-d.'
! dp-QM/MM; 0312BK20
      else if (ifmopt .eq. 14) then
          write(6,*) 'ann.dat: read atm_pol(dist), dp-QM/MM.'
      else
          write(6,*) 'ann.dat: error in reading IFMOPT; stop.'
      endif

!#############! 
!# BLOCK 1/4 #! 0711BK20
!#############!
! read in parameters, 0330BK20
! ifmopt=14 (atm_ipol(dist), dp-QM/MM)
! ----------------------------------------

      if (ifmopt .eq. 14) then
         read(90,*) help
         read(90,*) dmpon,cutof,z_cent,l_scal

! --- set z_cent and l_scal to zero to skip spline under tension!
         if (ifmopt .eq. 14) then
            aisum = 0.0d0
            aisum = z_cent+l_scal
            if (aisum .gt. 0) then
               write(6,'(A)') 'z_cent + l_scal .gt. 0, stop.'
               stop
            endif
         endif

! --- spline under tension variables
! --- recommended setting for z_cent = 0 A; (lower bound, z=0.0)
! --- recommended setting for l_scal = 2 A;    (midpoint, z=0.5)
! --- no settings for upper bound;   (infinty by default, z=1.0)
! --- l_scal is uniformly applied as a relative measure
         if (ifmopt .eq. 15) then
            aisum = 0.0d0
            aisum = z_cent+l_scal
            if (aisum .eq. 0) then
               write(6,'(A)') 'z_cent + l_scal .eq. 0, stop.'
               stop
            endif
         endif

! --- # of neural nodes, inputs, and outputs
         read(90,*) ainnode,aininp,ainout

! --- print help (HH)
         write(6,'(A,I3)') 'HH: help:', help
         write(6,'(A,4F9.3)') 'HH: dmpon,cutof,z_cent,l_scal:', &
            dmpon,cutof,z_cent,l_scal
         if (ifmopt .eq. 14) then
            write(6,'(A)') 'HH: z_cent and l_scal are INACTIVE!'
         endif
         if (ifmopt .eq. 15) then
            write(6,'(A)') 'HH: z_cent and l_scal are ACTIVE!'
         endif
         write(6,'(A,3I3)') 'HH: ainnode,aininp,ainout:', &
            ainnode,aininp,ainout
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

! --- input features composed of internal coordinates (IC)
! --- IC option is listed in aitype(i,1)
! --- IC atoms are listed in aitype(i,2-4), use 0's for bond/angles
! --- only bond option (2) is available.
! --- angle (3) and dihedral (4) option are removed.
! ----------------------------------------
         do i = 1, aininp
            read(90,*) (aitype(i,j),j=1, 5)
               do j = 1, 5
               end do
         end do

! --- net.IW{1} - INITIAL WEIGHTS [nnode x ninp]
         do i = 1, ainnode
            read(90,*) (aiiw(i,j),j=1,aininp)
               do j = 1, aininp
               end do
         end do

! --- net.b{1} - INITIAL BIASES [nnode x 1]
         do i = 1, ainnode
            read(90,*) aib1(i)
         end do

! --- net.LW{2} - OUTPUT WEIGHTS [nout x nnode]
         do i = 1, ainout
            read(90,*) (ailw(i,j),j=1,ainnode)
         end do

! --- net.b{2} - OUTPUT BIASES [nout x 1]
         do i = 1, ainout
            read(90,*) aib2(i)
         end do
      endif

!#############! 
!# CRD/FORCE #! 
!#############!

      natomx = 6465
      do i = 1, natomx
         read(91,*) x(i), y(i), z(i)
      end do

      do i = 1, ainout+3
         write(*,704) 'crd', i, x(i), y(i), z(i)
      end do
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

      do i = 1, natomx
         read(92,*) dx0(i), dy0(i), dz0(i)
         dx(i) = dx0(i)
         dy(i) = dy0(i)
         dz(i) = dz0(i)
      end do

      do i = 1, ainout
         write(*,704) 'force', i, dx(i), dy(i), dz(i)
      end do
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

!#############!
!# BLOCK 2/4 #! 0711BK20
!#############!
! ann gradient
! ifmopt=14 (atm_ipol(dist), dp-QM/MM)
! ----------------------------------------

! ---- 0330BK19 ann implementation IC
      if (ifmopt .eq. 14) then

! generate ANN inputs for BONDS
      do i = 1, aininp
         if (aitype(i,1) .eq. 2) then
            aiatm1(1) = x(aitype(i,2))
            aiatm1(2) = y(aitype(i,2))
            aiatm1(3) = z(aitype(i,2))

            aiatm2(1) = x(aitype(i,3))
            aiatm2(2) = y(aitype(i,3))
            aiatm2(3) = z(aitype(i,3))

            call get_r12(aiatm1,aiatm2,aipr)

! aip0 is a legacy variable for mapminmax in MATLAB
            aip0(i) = aipr

! assign real distance to input feature; 0331BK20
            if (ifmopt .eq. 14) then
               aip(i) = aip0(i)
            endif

! convert real distance to dimensionless variable; 0331BK20
            if (ifmopt .eq. 15) then
               aip(i) = (2/PI)*(atan((aip0(i)-z_cent)/l_scal))
            endif
         endif
      end do

! -- generate argument for tansig
      do i = 1, ainnode
         aisum = 0.0d0
         do j = 1, aininp
            aisum = aisum+aiiw(i,j)*aip(j)
         end do
         aiarg(i) = aisum+aib1(i)
      end do

! -- compute activation function for hidden layer (tansig)
      do i = 1, ainnode
         aia(i) = (2/(1+exp(-2*aiarg(i))))-1
      end do

! -- compute output from output layer
      do i = 1, ainout
         aisum = 0.0d0
         do j = 1, ainnode
            aisum = aisum+ailw(i,j)*aia(j)
         end do
         aians(i) = aisum+aib2(i)
      end do

! Isotropic atomic polarizability correction; 
      do i = 1, ainout
         aialp(i) = aians(i)
         if (help .eq. 1) then
         write (*,702) 'HH: aialp', i, aialp(i)
         endif
      end do
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

! ---- gradient on input(s) 
! -- variable substitution in expression for neural nodes
      do j = 1, ainnode
         aisum = 0.0
         do i = 1, aininp
            aisum = aisum+aiiw(j,i)*aip(i)
         end do
         aih(j) = exp(-2*(aib1(j)+aisum))
         aih2(j) = 1/((aih(j)+1)*(aih(j)+1))
      end do

! -- sum of expressions for each input
      do i = 1, ainout
         do j = 1, aininp
            aisum = 0.0d0
            do k = 1, ainnode
               aisum = aisum + &
                  (aih(k)*aih2(k)*aiiw(k,j)*ailw(i,k))
            end do
            aiddx(i,j) = 4*aisum
         end do
      end do

      if (help .eq. 1) then
         do i = 1, ainout
            do j = 1, aininp
               write(*,603) 'HH: nout,ninp,ddx', i,j,aiddx(i,j)
            end do
         end do
      endif

      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

      do i = 1, ainout
         do j = 1, ainout
            dalpdx(i,j) = 0.0d0
            dalpdy(i,j) = 0.0d0
            dalpdz(i,j) = 0.0d0
         end do
      enddo

      do i = 1, ainout
         do j = 1, aininp
            if (aitype(j,1) .eq. 2) then
               ii = aitype(j,2)
               jj = aitype(j,3)

               dax = x(ii)-x(jj)
               day = y(ii)-y(jj)
               daz = z(ii)-z(jj)

               dalpdx(i,ii) = dalpdx(i,ii)+aiddx(i,j)*(dax/aip(j))
               dalpdy(i,ii) = dalpdy(i,ii)+aiddx(i,j)*(day/aip(j))
               dalpdz(i,ii) = dalpdz(i,ii)+aiddx(i,j)*(daz/aip(j))

               dalpdx(i,jj) = dalpdx(i,jj)+aiddx(i,j)*(-dax/aip(j))
               dalpdy(i,jj) = dalpdy(i,jj)+aiddx(i,j)*(-day/aip(j))
               dalpdz(i,jj) = dalpdz(i,jj)+aiddx(i,j)*(-daz/aip(j))

            endif
         end do
      end do

      if (help .eq. 1) then
! format check
      do i = 1, ainout
         do j = 1, ainout
            write(*,801) 'HH: nout,natm,dx,dy,dz', i,j, &
               dalpdx(i,j),dalpdy(i,j),dalpdz(i,j)
         end do
         write(*,*) ''
      end do
      endif
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'
      endif ! if (ifmopt .eq. 14)

!#############! 
!# BLOCK 3/4 #! 0711BK20
!#############!
! electric field gradient
! ifmopt=14 (atm_ipol(dist), dp-QM/MM)
! ----------------------------------------

      if (ifmopt .eq. 14) then

! atomic polarizability obtained - compute the induction energy as a result of solvent polarization on atomic centers 

         do i = 1, ainout
            aiefx(i) = 0.0d0
            aiefy(i) = 0.0d0
            aiefz(i) = 0.0d0

            imol = 0

! Coordinates (X,Y,Z)
            xcm = x(i)
            ycm = y(i)
            zcm = z(i)

! obtain the permanent electric field generate by solvent on each atomic polarizability center
! --- loop over solvent atoms; temporary hack assuming solvent starts from solvi for all water
! --- use r3; 0315BK20; use xx/yy/zz & r1p/r3p to avoid confusion; 0316BK20
! --- bug fixed: xx -> -(xi-xj)/rij^3 - rank-1 tensor PIPF Eq.(7a); 0316BK20 

            do j = ainout+1, natomx
               xx = - (xcm - x(j))
               yy = - (ycm - y(j))
               zz = - (zcm - z(j)) 
               r1p = sqrt(xx*xx + yy*yy + zz*zz)
               r3p = r1p*r1p*r1p
!
! --- use [dmpon, cutof] as the cutoff interval to turn on the short-range/long-range dp-QM/MM; 0320BK20
! --- residue/molecule-based enabled use; 0316BK20
! --- set molecule indicator (if dmpon (e.g. 2.0 A) < CMS-O < cutof (e.g., 12.0 A), whole molecule in); 0320BK20
! --- the lower-bound dmpon is added as the short-range efield should be damped out & no force to prevent 
!     water get too close to CMS so far (if distributed to atomic polarizability, may be easy to damp; 0320BK20)
!
               jj = j-ainout ! index reset, 0509BK20

               if (mod(jj,3) .eq. 1) then
                  if (r1p .ge. dmpon .and. r1p .le. cutof)then
                     imol = 1
                  else
                     imol = 0
                  endif
               endif

! --- whether include in the field, using molecule indicator; if O excluded, the associated H's out; 0316BK20
               if (imol .eq. 1) then

! bug fixed; remainder 1 as O; otherwise, H; 0315BK20
                  if (mod(jj,3) .eq. 1) then
                     chg = -0.834d0
                  else 
                     chg = 0.417d0
                  endif

! --- permanent field: negative*rank-1 tensor follows PIPF Eq.(3); 0316PJ20 
                  aiefx(i) = aiefx(i) - (xx/r3p)*chg
                  aiefy(i) = aiefy(i) - (yy/r3p)*chg
                  aiefz(i) = aiefz(i) - (zz/r3p)*chg

               end if ! if (imol .eq. 1) then
            end do ! do j = ainout+1, natomx
         end do ! do i = 1, ainout

! compute E_pol (no damping) due to the induced dipole, mu_ind = alp*E^0 and E_pol = -(1/2) * mu_ind * E^0 (in charge*charge/A -> kcal/mol by *332 as conversion factor)
         do i = 1, ainout
            aimux(i) = aialp(i)*aiefx(i)
            aimuy(i) = aialp(i)*aiefy(i)
            aimuz(i) = aialp(i)*aiefz(i)
         end do

         epol = 0.0d0

         do i = 1, ainout
            epol = epol + -0.5d0 * 332.0716D0 * &
               (aimux(i)*aiefx(i) + &
                aimuy(i)*aiefy(i) + &
                aimuz(i)*aiefz(i))
         end do

! --- debug
         if (help .eq. 1) then
         write(*,'(A)') 'HH: --- USERE is called --- '
         do i = 1, ainout
            write(*,'(A,1F9.3)') 'HH: Alpha correction = ', aialp(i)
         end do
         write(*,'(A,1F9.3)') 'HH: E_pol correction = ', epol
         endif
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

! pass E_pol to user energy
         eu = epol

! --- now computing gradients associated with epol; bug fixed; 0316BK20

! --- assign the prefactor (1.0 prefactor in gradient is consistent 
! with 0.5 in E_pol as E_pol = -0.5*alp*dEf^2/dx = -0.5*alp*2*Ef*dEf/dx = alp*Ef*dEf/dx; 0316PJ20

         do i = 1, ainout
            gpfac = 1.0d0

            aidxcm(i) = 0.0d0
            aidycm(i) = 0.0d0
            aidzcm(i) = 0.0d0

            imol = 0

! Coordinates (X,Y,Z)
            xcm = x(i)
            ycm = y(i)
            zcm = z(i)

            do j = ainout+1, natomx
! --- xx remains (xi-xj); not to form the rank-1 tensor; but to form rank-2 tensor following PIPF Eq. (9)
               xx = xcm - x(j)
               yy = ycm - y(j)
               zz = zcm - z(j)
               r1p = sqrt(xx*xx + yy*yy + zz*zz)
!
! --- use [dmpon, cutof] as the cutoff interval to turn on the short-range/long-range dp-QM/MM; 0320BK20
! --- residue/molecule-based enabled use; 0316BK20
! --- set molecule indicator (if dmpon (e.g. 2.0 A) < CMS-O < cutof (e.g., 12.0 A), whole molecule in); 0320BK20
! --- the lower-bound dmpon is added as the short-range efield should be damped out & no force to prevent
!     water get too close to CMS so far (if distributed to atomic polarizability, may be easy to damp; 0320BK20)
!
               jj = j-ainout ! index reset, 0509BK20

               if (mod(jj,3) .eq. 1) then
                  if (r1p .ge. dmpon .and. r1p .le. cutof) then
                     imol = 1
                  else 
                     imol = 0
                  endif
                endif

! --- whether include in the field, using molecule indicator; if O excluded, the associated H's out; 0316BK20
               if (imol .eq. 1) then

! bug fixed; remainder 1 as O; otherwise, H; 0315BK20
                  if (mod(jj,3) .eq. 1) then
                     chg = -0.834d0
                  else 
                     chg = 0.417d0
                  endif

! --- get the rank-2 tensor for gradient; follow PIPF eq. (9) --- no damping - code adapted from pipf/epipf.src
                  r2 = 1.0d0/(r1p*r1p)
                  r1 = sqrt(r2)
                  r3 = r1*r2

                  rr3 = 3.0d0*r3*r2
                  rr3xx = rr3*xx
                  rr3yy = rr3*yy
                  rr3zz = rr3*zz

                  txww(1) = rr3xx*xx-r3 
                  txww(2) = rr3xx*yy 
                  txww(3) = rr3yy*yy-r3 
                  txww(4) = rr3xx*zz 
                  txww(5) = rr3yy*zz 
                  txww(6) = rr3zz*zz-r3 

!
!      WE FOLLOW THE EQ. (30) IN CHEM. PHYS. LETT. 1990, 166,180. NOTE THAT
!      CHARMM ARRAY DX, DY, DZ ACTUALLY CONTAINS GRADIENT (DV/DX) INSTEAD OF
!      FORCE (-DV/DX)
!
! --- first portion of first term gone for doubly polarizable QM/MM --- non charge on center of mass of solute
!      FIRST, WE COMPUTE THE GRADIENT ASSOCIATED WITH PERMENENT CHARGE
!      (THE SECOND TERM IN EQ. 30)
!
!          Q^K * PI^K   = Q^K * sum (T^KA u^A)     (r denotes lamda)
!                  r             A    rb   b
!
!               DX(I) = DX(I) - CGI*(TXWW(1,NB)*DMUIND(1,J) &
!                      +TXWW(2,NB)*DMUIND(2,J)+TXWW(4,NB)*DMUIND(3,J))
!               DY(I) = DY(I) - CGI*(TXWW(2,NB)*DMUIND(1,J) &
!                       +TXWW(3,NB)*DMUIND(2,J)+TXWW(5,NB)*DMUIND(3,J))
!               DZ(I) = DZ(I) - CGI*(TXWW(4,NB)*DMUIND(1,J) &
!                       +TXWW(5,NB)*DMUIND(2,J)+TXWW(6,NB)*DMUIND(3,J))

                  dx(j) = dx(j) - gpfac * chg * (txww(1)*aimux(i) + &
                     txww(2)*aimuy(i) + txww(4)*aimuz(i)) * 332.0716D0
                  dy(j) = dy(j) - gpfac * chg * (txww(2)*aimux(i) + &
                     txww(3)*aimuy(i) + txww(5)*aimuz(i)) * 332.0716D0
                  dz(j) = dz(j) - gpfac * chg * (txww(4)*aimux(i) + &
                     txww(5)*aimuy(i) + txww(6)*aimuz(i)) * 332.0716D0

!
! --- 2nd portion of 2nd term gone for doubly polarizable QM/MM --- non dipole on solvent
! --- 1st portion of 2nd term is action/reaction force from 1st portion of 1st term
!      NEXT, WE COMPUTE THE GRADIENT ASSOCIATED WITH PERMENENT FIELD GRADIENT
!      (THE THIRD TERM IN EQ. 30)
!
!          u^K * E^K   = u^K * [- sum (T^KA Q^A)]    (r denotes lamda)
!           a     ar      a        A    ar

                  aidxcm(i) = aidxcm(i) + gpfac * chg * &
                     (txww(1)*aimux(i) + txww(2)*aimuy(i) + &
                      txww(4)*aimuz(i))  * 332.0716D0
                  aidycm(i) = aidycm(i) + gpfac * chg * &
                     (txww(2)*aimux(i) + txww(3)*aimuy(i) + &
                      txww(5)*aimuz(i))  * 332.0716D0
                  aidzcm(i) = aidzcm(i) + gpfac * chg * &
                     (txww(4)*aimux(i) + txww(5)*aimuy(i) + &
                      txww(6)*aimuz(i)) * 332.0716D0

!               DX(J) = DX(J) + CGI*(TXWW(1,NB)*DMUIND(1,J)   &
!                       +TXWW(2,NB)*DMUIND(2,J)+TXWW(4,NB)*DMUIND(3,J))
!               DY(J) = DY(J) + CGI*(TXWW(2,NB)*DMUIND(1,J) &
!                       +TXWW(3,NB)*DMUIND(2,J)+TXWW(5,NB)*DMUIND(3,J))
!               DZ(J) = DZ(J) + CGI*(TXWW(4,NB)*DMUIND(1,J) &
!                       +TXWW(5,NB)*DMUIND(2,J)+TXWW(6,NB)*DMUIND(3,J))

               end if ! if (imol .eq. 1) then
            end do ! do j = ainout+1, natomx
         end do ! do i = 1, ainout

! now redistribute gradient to solute atoms: dE/dx_j = dE/dX_CMS * dX_CMS/dx_j = dE/dX_CMS * (m_j/M)
         do i = 1, ainout
            dx(i) = dx(i) + aidxcm(i)
            dy(i) = dy(i) + aidycm(i)
            dz(i) = dz(i) + aidzcm(i)
         end do

!
! Finally, take d_alp/dx into account based on spline function; dE_pol/dx = -0.5*d_alpha/dx * Ef^2; 0316PJ20
! Note: although loop over all solute atoms, only CV atoms matter, other atoms have no dalpd[x-z] contributions.
! Note in regular PIPF - no gradient on alp using fixed alp, but now alp is geometry dependent, therefore alp has
! nuclear gradient (computed from cubic spline function - endpoints may need more careful treatment if drift ); 0316PJ20
!
         do i = 1, ainout
            aief2(i) = aiefx(i)*aiefx(i) + &
               aiefy(i)*aiefy(i) + aiefz(i)*aiefz(i)
            if (help .eq. 1) then
            write(*,702) 'HH: aief2:', i, aief2(i)
            endif
         end do
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'

         do i = 1, ainout
            do j = 1, ainout
               dx(i) = dx(i) - 0.5d0 * dalpdx(i,j)*aief2(j) * 332.0716D0
               dy(i) = dy(i) - 0.5d0 * dalpdy(i,j)*aief2(j) * 332.0716D0
               dz(i) = dz(i) - 0.5d0 * dalpdz(i,j)*aief2(j) * 332.0716D0
            end do
         end do

         if (help .eq. 1) then
         do i = 1, ainout
            write(*,704) 'HH: ftot:', i, dx(i), dy(i), dz(i)
         end do
         endif
      if (help .eq. 1) write(*,'(A)') 'oooooooooooooooooooooooooooooooo'
      endif ! if (ifmopt .eq. 14)

!#############! 
!# BLOCK 4/4 #! 0711BK20
!#############!
! end of dp-QM/MM
! ifmopt=14 (atm_ipol(dist), dp-QM/MM)
! ----------------------------------------
      603 format(A,I3,I3,1F9.3)
      702 format(A,I3,1F9.3)
      704 format(A,I3,3F9.3)
      801 format(A,I3,I3,3F9.3)
      end

!pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
!
! add spline rountines, taken from Numerical Recipe; but named it as -pj
! to avoid potential conflicts with existing spline routines. 0418PJ16
!
!******************************************************************
!* given arrays y and x of dimension n such that y = f(x) and
!* x(1)<x(2)<...<x(n) and given values yp1 and ypn for the first
!* derivatives
!* at 1 and n , this subroutine returns an array y2(n) which gives
!* the 2nd
!* derivative of the interpolating function.  if yp1 and/or ypn are
!* >= 1.E30 the natural spline will be calculated, that is, for
!*  boundary
!* conditions >= 1.E30 the 2nd derivative will be zero.
! n is the maxium dimension of x and y
! nr is the real dimension of x and y
!******************************************************************

      subroutine splnpj(x,y,n,nr,yp1,ypn,y2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer n,nr
      real*8 x,y,y2,u,sig,p,qn,un
!      dimension x(n),y(n),y2(n),u(n)
      dimension x(*),y(*),y2(n),u(n)

      write(*,'(A)') 'call spln'
      if(yp1 .gt. 0.99E30)then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      write(*,'(A)') 'y2(1)=',y2(1),'u(1)=',u(1)


      do 11 i=2,nr-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1) + 2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
      (x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1))/p
         write(*,'(A)') 'i=',i,'sig=',sig,'p=',p,&
         'y2(i)=',y2(i),'u(i)=',u(i)
   11 continue
      if (ypn .gt. 0.99E30) then
         qn=0.
         un=0.
      else
         qn=0.5
         un=(3./(x(nr)-x(nr-1)))*(ypn-(y(nr)-y(nr-1))/(x(nr)-x(nr-1)))
      endif
      write(*,'(A)') 'qn=',qn,'un=',un
      y2(nr)= (un-qn*u(nr-1))/(qn*y2(nr-1)+1.)
      do 12 k=nr-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      write(*,'(A)') 'y2(k)=',y2(k),'y2(k+1)=',y2(k+1),'u(k)=',u(k)
   12 continue
      write(*,'(A)') 'y2(nr):'
      write(*,'(A)') (y2(i),i=1,nr)
      write(*,'(A)') 'return from spln'
      return
      end subroutine splnpj

!******************************************************************
!* given xa(n) and ya(n) such that ya = f(xa) and xa is
!* monotonically increasing, and the array y2a(n) (output from
!* subroutine spln), and a given value x, this subroutine gives the
!* interpolated value y at the given x. 
!* --- modified to also return 1st derivative at dy=y'(x); 0316PJ20
!******************************************************************
!********1*********2*********3*********4*********5*********6*******

      subroutine splintpj2(xa,ya,y2a,n,nr,x,y,dy)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer n,klo,khi,k,nr
      real*8 xa,ya,y2a,x,y,h,a,b,dy
!      dimension xa(n),ya(n),y2a(n)
      dimension xa(*),ya(*),y2a(*)
!------this subroutine will make a point which is out of the range-----
!------defined in your IUMB input file, to lie between the two -----
!------points which are the two closest points to that point----
!------but the result is still ok------
!      write(*,'(A)') 'call splint'
      klo=1
      khi=nr
      write(*,'(A)') 'klo=',klo,'khi=',khi
    1 if (khi-klo.gt.1)then
         k=(khi+klo)/2
         write(*,'(A)') 'khi>klo, k=',k
         if(xa(k).gt.x)then
            khi=k
            write(*,'(A)') 'xa(',k,')=',xa(k),'>rc(',x,'),khi=k=',khi
         else
            klo=k
            write(*,'(A)') 'xa(',k,')=',xa(k),'<rc(',x,'),klo=k=',klo
         endif
         goto 1
      endif
      write(*,'(A)') 'finally,rc(',x,') lies between',xa(klo),'and',xa(khi)
      h=xa(khi)-xa(klo)
      write(*,'(A)') 'h=',h
      if (h.eq.0.) then
         y = 0.0
! added to return 0 for 1st derivative; 0316PJ20
         dy = 0.0
         write(*,'(A)') 'bad xa input in splint'
         STOP
      else
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo) + b*ya(khi) + &
       ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6
!
! added to return 1st derivative; bug fixed; 0316PJ20    
!
         dy= -ya(klo)/h + ya(khi)/h + &
            (-(3*a*a-1)*y2a(klo) + (3*b*b-1)*y2a(khi))*h/6

      endif
      write(*,'(A)') 'a=',a,'b=',b,'EUMB=',y,'dy=',dy
      write(*,'(A)') 'return from splint'
      return
      end subroutine splintpj2

! calculate distance between two 3-d cartesian coordinates
      subroutine get_r12(x1,x2,r12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer p
      real*8 x1(3),x2(3),r2,r12
      r2 = 0d0
      do p = 1, 3
         r2 = r2+((x1(p)-x2(p))*(x1(p)-x2(p)))
      end do
      r12 = sqrt(r2)
!      write(*,'(A)') r12
      return
      end subroutine get_r12

