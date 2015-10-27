cc===================================================================
c
c  This UMAT subroutine implements the following model in ABAQUS.
c
c-------------------------------------------------------------------
!#define DEBUG
#define DAMAGE
!#define WRITEALL
#define CALLPYTHON

      subroutine umat (stress,  statev,  ddsdde,  sse,     spd,
     &			scd,     rpl,     ddsddt,  drplde,  drpldt,
     &			strain,  dstrain, time,    dtime,   temp,
     &			dtemp,   predef,  dpred,   cmname,  ndi,
     &			nshr,    ntens,   nstatv,  props,   nprops,
     &			coords,  drot,    pnewdt,  celent,  dfgrd0,
     &			dfgrd1,  noel,    npt,     layer,   kspt,
     &			kstep,   kinc )

      include 'ABA_PARAM.INC'
c      implicit double precision (a-h,o-z)

      real*8 mesh_size !make sure this is a float and not an integer so crack length correct
      integer phase
      real*8 elastic_modulus
      real*8 bulk_modulus
      real*8 shear_modulus
      real*8 lame_parameter
      real*8 poisson_ratio
#include "Geom_Def.txt"

      parameter(num_grains_UMAT = 1) !NOTE McGinty's code loops over grains in UMAT, should always be set to 1
									 !this is different from num_grains, which depends on mesh/input file and is used by UEXTERNALDB
      character*8 cmname, ans
      logical Converged, Improved
	  
c-------------------------------------------------------------------
c  Dimension arrays passed into the UMAT sub
c-------------------------------------------------------------------

      dimension
     &	coords(3),	! Coordinates of Gauss pt. being evaluated
     &	ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &	ddsddt(ntens),	! Change in stress per change in temperature
     &	dfgrd0(3,3),	! Deformation gradient at beginning of step
     &	dfgrd1(3,3),	! Deformation gradient at end of step
     &	dpred(1),	! Change in predefined state variables
     &	drplde(ntens),	! Change in heat generation per change in strain
     &	drot(3,3),	! Rotation matrix
     &	dstrain(ntens),	! Strain increment tensor stored in vector form
     &	predef(1),	! Predefined state vars dependent on field variables
     &	props(nprops),	! Material properties passed in
     &	statev(nstatv),	! State Variables
     &	strain(ntens),	! Strain tensor stored in vector form
     &	stress(ntens),	! Cauchy stress tensor stored in vector form
     &	time(2)		! Step Time and Total Time
                                                                                
c-------------------------------------------------------------------            
c  Dimension other arrays used in this UMAT
c-------------------------------------------------------------------            

      dimension
     &	array1	(3,3),		! Dummy array
     &	array2	(3,3),		! Dummy array
     &	array3	(num_slip_sys,num_slip_sys), ! A_matrix of Newton Raphson
     &	array4	(6,6),		! Dummy array used in Voigt notation
     &	array5	(6,6),		! Inverse of array4().
     &	array6	(3,3,3,3),	! 4th rank dummy array
     &	a0	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at beginning of step
     &	a	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at end of step
     &	a0_1	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at beginning of step
     &	a_1	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at end of step
     &	a0_2	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at beginning of step
     &	a_2	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at end of step
     &	a0_3	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at beginning of step
     &	a_3	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at end of step	 
     &	a0_4	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at beginning of step
     &	a_4	(num_slip_sys,num_grains_UMAT), ! Kinematic stress at end of step	 	 
     &  DAM      (3,3,3,3),      ! Local  4th rank tensor from (I-D)
     &	C0	(3,3,3,3),	! Local  4th rank elastic stiffness tensor
     &	C	(3,3,3,3),	! Global 4th rank elastic stiffness tensor
     &	C_avg	(3,3,3,3),	! Average over all grains
     &	del	(3,3),		! Kronecker delta tensor
     &	ddpdsig	(3,3,3,3),	! deriv of D_p wrt sig * dt
     &	ddsdde_4th(3,3,3,3),	! 4th rank tangent stiffness tensor.
     &	daadga	(num_slip_sys),	! deriv of a_alpha wrt del_gamma_alpha
     &	dgadgb	(num_slip_sys,num_slip_sys),! deriv of g_alpha wrt del_gamma_beta
     &	d_gamma_dot(num_slip_sys),! Gamma_dot step from Newt-Raphson
     &	dir_cos	(3,3,num_grains_UMAT),! Direction cosines 
     &  R       (12,3,3),       ! Rotation matrix towards normal of slipl plane
     &	dTaudGd	(num_slip_sys,num_slip_sys), ! Deriv of Tau wrt Gamma_dot
     &	E_el	(3,3),		! Elastic Green strain tensor
     &	E_tot	(3,3),		! Green strain tensor E_el + E_p
     &	F0	(3,3),		! F at beginning of sub increment
     &	F1	(3,3),		! F at end of sub increment
     &	F_el	(3,3),		! Elastic part of F
     &	F_el_inv(3,3),		! Inverse of elastic part of F
     &	F_p_inv_0(3,3,num_grains_UMAT),! Inverse of plastic part of F at beginning
     &	F_p_inv(3,3,num_grains_UMAT),! Inverse of plastic part of F at end of step
     &	func	(num_slip_sys),	! Function that is solved to get gamma_dot(k)
     &	g0	(num_slip_sys,num_grains_UMAT), ! Ref stress at beginning of step
     &	g	(num_slip_sys,num_grains_UMAT), ! Ref stress at end of step
     &	gamma_dot(num_slip_sys,num_grains_UMAT),! shear strain rate on system
     &	gamma_try(num_slip_sys),! Candidate gamma_dots for line search
     &	grad(num_slip_sys),	! gradient of sum-sqares-error
     &	psi	(3,num_grains_UMAT),	! Euler angles
     &	sig	(3,3),		! Stress tensor
     &	sig_avg	(3,3),		! Rate of stress change for a grain
     &	Spk2	(3,3),		! 2nd Piola Kirkhoff stress
     &	tau	(num_slip_sys,num_grains_UMAT),! Resolved shear stress for slip dir.
     &  sigma   (num_slip_sys),! Normal resolved stress for slip dir.
     &	xL_p	(3,3),		! Plastic vel grad in current configuration
     &	xs0	(3,num_slip_sys),! Inter config slip directions in global coords
     &	xs	(3,num_slip_sys),! Current config slip directions in global coords
     &	xm0	(3,num_slip_sys),! Inter config plane normals in global coords
     &	xm	(3,num_slip_sys),! Current config plane normals in global coords
     &	y	(3,num_slip_sys),! Miller indices of slip plane normals 
     &	z	(3,num_slip_sys), ! Miller indices of slip directions
     &  yxz     (3,num_slip_sys),         ! Cross product between slip plane normals &slip directions 
     &  delta_E_p(3,3),         !increment of plastic strain tensor
     &  delta_gamma(num_slip_sys), !increment of shear strain on each slip sys
     & delta_gamma_cum(num_slip_sys), ! Acumulated increment of shear strain on each slip sys
     &	E_p(3,3),
     &  delta_gamma_cum_max(num_slip_sys),
     &  delta_gamma_cum_min(num_slip_sys)
	 
	 

       common/KUEXT_FIP/FIP_gl(num_grains,2*max_num_layers+1,num_slip_sys),
     & FIP_elem(num_elem,num_slip_sys),
     & FIP_elem_max(num_elem,3),! FIP_elem_max, plane #, normal stress
     & FIP_elem_avg_max(num_elem,2)
	 
#include "Common_block_Alv02.txt"
#include "Definitions.txt" !defines constants   

     
c************************************      
c  Constants assigned            
      Pi = 4. * atan(1.)
c************************************

c-----------------------------------------------------
c  Set the phase (Al matrix=1 vs particle=2)
c-----------------------------------------------------
      phase = 1
      if(nprops>9) then
         phase = props(10)
      end if


      if(phase.eq.1) then
c       print*
c       print*, 'UMAT CALLED for: ',kstep, kinc, noel, dtime
c       print*, 'UMAT CALLED for: ',kstep, kinc, noel
      if (kinc.gt.4000) then
       print*, 'UMAT CALLED for: ',kstep, kinc, noel
      end if
c-------------------------------------------------------------------
c Terminate simulation if crack has arrested
c-------------------------------------------------------------------

#ifdef DAMAGE
      if ((kstep.gt.1).and.(crack_arrested)) then
          print*
          print*, 'SIMULATION TERMINATED DUE TO CRACK ARREST'
          CALL XIT
      end if
#endif

c-------------------------------------------------------------------
c  Assign props() array to logical variable names
c-------------------------------------------------------------------
      
c Constants form Gustavo;s fit 31:	  
      psi_ang   = props(1)
      theta_ang = props(2)
      phi_ang   = props(3)

      C11	= 107300.      !107300.
      C12	= 60900.       !60900.
      C44	= 28300.       !28300.
	  
      gamma_dot_zero= 0.001 !0.001
      flow_exp	= 75.      !65.
	  
      g_zero	= 35.     !240.
      Hdir	= 0.0d0         !225.
      Hdyn	= 0.0d0          !1.
      xLatent	=1.4       !1.4
	  
      a_zero    = 0.      !35.
      Adir      = 75600.     !350.
      Adyn	= 720.           !1.
	  
      a_zero_1    =  0.0d0      !35.
      Adir_1      = 2000000.     !350.
      Adyn_1	= 20000.          !1.
      b_1 = Adir_1/Adyn_1 ! Saturation level
c      print*, "b_1: ", (b_1)
	  
      a_zero_2    = 0.0d0      !35.
      Adir_2      = 135000.     !350.
      Adyn_2	= 1421.05        !1.
      b_2 = Adir_2/Adyn_2 ! Saturation level
c      print*, "b_2: ", (b_2)	  
	  
      a_zero_3    = 0.0d0      !35.
      Adir_3      = 0.0d0     !350.
      Adyn_3	= 0.0d0          !1.
c      b_3 = Adir_3/Adyn_3 ! Saturation level
c      print*, "b_3: ", (b_3)	 	  

      a_zero_4    = 0.     !35.
      Adir_4     = 0.     !350.
      Adyn_4	= 0.          !1.	  
c      b_4 = Adir_4/Adyn_4 ! Saturation level	
c      print*, "b_4: ", (b_4)	 	  
c      print*, "SATURATION LEVEL: ", (b_1+b_2+b_3+b_4)
	  
      m_OW2 = 70 !ohno-wang exponent (use even integer)
	 
	  
   
c-------------------------------------------------------------------
c  Initialize Kronnecker delta tensor
c-------------------------------------------------------------------

      do i = 1,3
         do j = 1,3
            del(i,j) = 0.0
         end do
         del(i,i) = 1.0
      end do

c-------------------------------------------------------------------
c  Assign slip system normals and slip directions for an FCC.
c-------------------------------------------------------------------

c     !   plane      dir
c     /  1, 1, 1,  0, 1,-1 /  
c     /  1, 1, 1, -1, 0, 1 /  
c     /  1, 1, 1,  1,-1, 0 /  
c     /  1,-1,-1,  0,-1, 1 /  
c     /  1,-1,-1, -1, 0,-1 /  
c     /  1,-1,-1,  1, 1, 0 /  
c     / -1,-1, 1,  0,-1,-1 /
c     / -1,-1, 1,  1, 0, 1 /
c     / -1,-1, 1, -1, 1, 0 /
c     / -1, 1,-1,  0, 1, 1 /  
c     / -1, 1,-1,  1, 0,-1 /  
c     / -1, 1,-1, -1,-1, 0 / 


        y(1,1) =  1.0
        y(2,1) =  1.0
        y(3,1) =  1.0
     
        y(1,2) =  1.0
        y(2,2) =  1.0
        y(3,2) =  1.0
        
        y(1,3) =  1.0
        y(2,3) =  1.0
        y(3,3) =  1.0
        
        y(1,4) =  1.0
        y(2,4) = -1.0
        y(3,4) = -1.0
       
        y(1,5) =  1.0
        y(2,5) = -1.0
        y(3,5) = -1.0
        
        y(1,6) =  1.0
        y(2,6) = -1.0
        y(3,6) = -1.0
        
        y(1,7) = -1.0
        y(2,7) = -1.0
        y(3,7) =  1.0

        y(1,8) = -1.0
        y(2,8) = -1.0
        y(3,8) =  1.0
        
        y(1,9) = -1.0
        y(2,9) = -1.0
        y(3,9) =  1.0
        
        y(1,10)= -1.0
        y(2,10)=  1.0
        y(3,10)= -1.0
        
        y(1,11)= -1.0
        y(2,11)=  1.0
        y(3,11)= -1.0
        
        y(1,12)= -1.0
        y(2,12)=  1.0
        y(3,12)= -1.0
      

	z(1,1) =  0.0               
        z(2,1) =  1.0
	z(3,1) = -1.0
               
	z(1,2) = -1.0               
	z(2,2) =  0.0
        z(3,2) =  1.0              

        z(1,3) =  1.0
	z(2,3) = -1.0               
	z(3,3) =  0.0               

	z(1,4) =  0.0               
	z(2,4) = -1.0               
        z(3,4) =  1.0

	z(1,5) = -1.0               
	z(2,5) =  0.0               
	z(3,5) = -1.0               

        z(1,6) =  1.0
        z(2,6) =  1.0
	z(3,6) =  0.0               

	z(1,7) =  0.0               
        z(2,7) = -1.0
        z(3,7) = -1.0

        z(1,8) =  1.0
	z(2,8) =  0.0               
	z(3,8) =  1.0               

	z(1,9) = -1.0               
	z(2,9) =  1.0               
	z(3,9) =  0.0               

	z(1,10)=  0.0               
	z(2,10)=  1.0               
	z(3,10)=  1.0               

        z(1,11)=  1.0
	z(2,11)=  0.0               
        z(3,11)= -1.0

	z(1,12)= -1.0               
        z(2,12)= -1.0 
	z(3,12)=  0.0  


c-------------------------------------------------------------------
c  Define normal to y and z by doing the cross product y x z
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
      call cross_product(y(1,i),y(2,i),y(3,i),z(1,i),
     &z(2,i),z(3,i),yxz(1,i), yxz(2,i), yxz(3,i))
      end do

c-------------------------------------------------------------------
c  Normalise the Miller indices to length one.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
         call normalize_vector( y(1,i), y(2,i), y(3,i) )
         call normalize_vector( z(1,i), z(2,i), z(3,i) )
         call normalize_vector( yxz(1,i), yxz(2,i), yxz(3,i) )
      end do

        do j=1,3
        do i=1,num_slip_sys
        R(i,j,1)=y(j,i)
        R(i,j,2)=z(j,i)
        R(i,j,3)=yxz(j,i)
c        print*,i, R(i,j,1),R(i,j,2),R(i,j,3)
        end do
        end do

c-------------------------------------------------------------------
c  Initialize internal variables for initial time step
c-------------------------------------------------------------------

      if (time(2) .eq. 0.0) then
	  
c-------------------------------------------------------------------
c  Check for normality of Miller indices.
c-------------------------------------------------------------------

      do k = 1,num_slip_sys
        sum = 0.0
        do i = 1,3
          sum = sum + y(i,k) * z(i,k)
        end do
        if (abs(sum) .gt. tolerance) then
          print*,'The Miller indices are WRONG!!!'
          print*,'on slip system # ',k
          STOP
        end if
      end do

#ifdef DEBUG
      print*
      print*, "Check01 - Time = ", time(2) ,noel
#endif

c-------------------------------------------------------------------
c  Generate Euler angles  1-3
c-------------------------------------------------------------------

         do i = 1,num_grains_UMAT
           psi(1,i) = psi_ang
           psi(2,i) = theta_ang
           psi(3,i) = phi_ang
         end do

#ifdef DEBUG
      print*
      print*, "Check01a - Time = ", time(2) ,noel
#endif
c-------------------------------------------------------------------
c  Initialize reference shear stress for each slip system of
c  each grain.  13-24
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do n = 1,num_slip_sys
             g0(n,m) = g_zero
           end do
         end do
            
#ifdef DEBUG
      print*
      print*, "Check01b - Time = ", time(2) ,noel
#endif
c-------------------------------------------------------------------
c  Initialize kinematic stress for each slip system of
c  each grain.   25-36
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do n = 1,num_slip_sys
             a0(n,m) = a_zero
             a0_1(n,m) = a_zero_1
             a0_2(n,m) = a_zero_2
             a0_3(n,m) = a_zero_3	
             a0_4(n,m) = a_zero_4				 
           end do
         end do
            
#ifdef DEBUG
      print*
      print*, "Check01c - Time = ", time(2) ,noel
#endif
c-------------------------------------------------------------------
c  Initialize F_p_inv_0    4-12
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,3
             do j = 1,3
               F_p_inv_0(i,j,m) = 0.0d0
             end do
             F_p_inv_0(i,i,m) = 1.0
           end do
         end do
#ifdef DEBUG
      print*
      print*, "Check01d - Time = ", time(2) ,noel
#endif
c-------------------------------------------------------------------
c  Initialize E_p  37-45
c-------------------------------------------------------------------

        do i = 1,3
          do j = 1,3
           E_p(i,j) = 0.0d0
          end do
       end do  
	   
#ifdef DEBUG
      print*
      print*, "Check01e - Time = ", time(2) ,noel
#endif
c-------------------------------------------------------------------
c  Initialize E_eff   46
c-------------------------------------------------------------------

        E_eff = 0.0d0

c-------------------------------------------------------------------
c  Initialize E_p_eff    47
c-------------------------------------------------------------------

        E_p_eff = 0.0d0

c-------------------------------------------------------------------
c  Initialize E_p_eff_cum   48
c-------------------------------------------------------------------

        E_p_eff_cum = 0.0d0  
                
c-------------------------------------------------------------------
c  Initialize  delta_gamma_cum   49 - (49+num_slip_sys)
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
       delta_gamma_cum(i) = 0.0d0
      end do

c-------------------------------------------------------------------
c  Initialize global variables
c-------------------------------------------------------------------
       do j=1,4
       Elem_pos(noel,j)=int(props(3 +j)) !layer to which element belongs
       end do

       Elem_pos(noel,5)=int(props(8))  ! Grain to which the element belongs

       do i = 1, num_slip_sys
       gamma_cum_element(noel,i)=0.0
       delta_gamma_cum_max(i)=0.0
       delta_gamma_cum_min(i)=0.0
       end do

      d1=0.0
      statev(134)=0.0 ! Damage
        
      ! Arrays for stress-strain plotting
      do i = 1,num_elem
      do j = 1,3
      do k = 1,3			 
        sig_elem(i,j,k) = 0.0
        E_tot_elem(i,j,k) = 0.0	
        E_p_elem(i,j,k) = 0.0	  	  
      end do
      end do 
      end do 	  
	  
      ! Arrays for back stress plotting
      do i = 1,num_elem	  
      do j = 1,num_slip_sys  	  
         a_elem(i,j) = 0.0 
      end do 
      end do
		
c-------------------------------------------------------------------
c  End of initializations.  Read in internal variables.
c-------------------------------------------------------------------


#ifdef DEBUG
      print*
      print*, "End initialization - Total Time = ", time(2)
#endif

      else  ! time<>0
c If time is not 0 read in values from previous step instead of initialize

         n = 0
            
c-------------------------------------------------------------------
c  Read in Euler Angles 1-3
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,3
             n = n + 1
             psi(i,m) = statev(n)
           end do
         end do
		 
c        print*, "n-chk-1", n
            
c-------------------------------------------------------------------
c  Read inverse of the plastic part of F     4-12
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do j = 1,3
             do i = 1,3
               n = n + 1
               F_p_inv_0(i,j,m) = statev(n)
             end do
           end do
         end do
            
c-------------------------------------------------------------------
c  Read reference shear values       13-24
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             g0(i,m) = statev(n)
           end do
         end do

c        print*, "n-chk-2", n		 
		 
c-------------------------------------------------------------------
c  Read kinematic stress values     25-84
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             a0(i,m) = statev(n)
           end do
         end do
		 
         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             a0_1(i,m) = statev(n)
           end do
         end do

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             a0_2(i,m) = statev(n)
           end do
         end do	

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             a0_3(i,m) = statev(n)
           end do
         end do		

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1		 
             a0_4(i,m) = statev(n)
c             print*, "Read a_4: ", n, a0_4(i,m)				 
           end do
         end do
		 
c        print*, "n-chk-3", n		
		 
c-------------------------------------------------------------------
c  Read E_p        85-93
c-------------------------------------------------------------------

        do i = 1,3
         do j = 1,3
          n = n + 1
          E_p(i,j) = statev(n)
         end do
        end do
          
c-------------------------------------------------------------------
c  Read E_eff      94
c-------------------------------------------------------------------
        n=n+1
        E_eff = statev(n)

c-------------------------------------------------------------------
c  Read E_p_eff      95
c-------------------------------------------------------------------
        n=n+1
        E_p_eff = statev(n)

c-------------------------------------------------------------------
c  Read E_p_eff_cum   96
c-------------------------------------------------------------------
      n = n+1
      E_p_eff_cum = statev(n)  
		
c      print*, "n-chk-4", n		

c-------------------------------------------------------------------
c  Read delta_gamma_cum   97-108
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
       n = n+1
       delta_gamma_cum(i) = statev(n) 
      end do
   
c      print*, "n-chk-4.5", n
c-------------------------------------------------------------------
C   Read accumulated plastic strain values   n= 109-132
c ------------------------------------------------------------------
      do i = 1,num_slip_sys
            n = n + 1
         delta_gamma_cum_max(i) = statev(n)
            n = n + 1
         delta_gamma_cum_min(i) = statev(n)
      end do
			 
c      print*, "n-chk-4.6", n				 
c-------------------------------------------------------------------
c  Read FIP  133
c-------------------------------------------------------------------
      n=  n + 1
      FIP = statev(n)
c      print*, "n-chk-5", n	  
	  
c-------------------------------------------------------------------
c  Read Damage  134
c-------------------------------------------------------------------	  
      n=  n + 1
      d1= statev(n)	  
c      print*, "n-chk-6", n	  

c-------------------------------------------------------------------
c  Read Stress normal to the crack plane of the element  150
c-------------------------------------------------------------------

      sigma_crack_plane = statev(150)
      
c-------------------------------------------------------------------
c  End of initializations
c-------------------------------------------------------------------

      end if ! (time = 0)

c-------------------------------------------------------------------
c  Calculate direction cosines based on Euler angles Bunge notation.
c http://en.wikipedia.org/wiki/Euler_angles Z1X2Z3 form
c-------------------------------------------------------------------


      do i = 1,num_grains_UMAT


        s1 = sin(psi(1,i))
        c1 = cos(psi(1,i))
        s2 = sin(psi(2,i))
        c2 = cos(psi(2,i))
        s3 = sin(psi(3,i))
        c3 = cos(psi(3,i))
            
        
         dir_cos(1,1,i) = c1*c3-s1*s3*c2
         dir_cos(2,1,i) = s1*c3+c1*s3*c2
         dir_cos(3,1,i) = s3*s2
         dir_cos(1,2,i) = -c1*s3-s1*c3*c2
         dir_cos(2,2,i) = -s1*s3+c1*c3*c2
         dir_cos(3,2,i) = c3*s2
         dir_cos(1,3,i) = s1*s2
         dir_cos(2,3,i) = -c1*s2
         dir_cos(3,3,i) = c2   

      end do             
		  
c-------------------------------------------------------------------
c  Initialize ANISOTROPIC 4th rank elastic stiffness tensor
c-------------------------------------------------------------------

          do i = 1,3          
           do j = 1,3         
            do k = 1,3        
             do l = 1,3       
                C0(i,j,k,l) = (1-d1)*  C12 * del(i,j) * del(k,l) +               
     &           (1-d1)*  C44 * (del(i,k)*del(j,l)+del(i,l)*del(k,j))        
             end do           
            end do            
           end do             
          end do              
          C0(1,1,1,1) = (1-d1)* C11   
          C0(2,2,2,2) = (1-d1)* C11   
          C0(3,3,3,3) = (1-d1)* C11
		  

      if (kinc.gt.4000) then
       print*, 'A'
      end if
	 
c--------------------------------------------------------------------
c  Apply damage to elastic stiffness tensor TODO + FIXME
c--------------------------------------------------------------------
c          print*, 'CHECK_A'
#ifdef DAMAGE

      if (dam_elem(noel).eq.1.0) then !element is in crack
	  
          if (dtime.gt.H_cycle/100) then
          pnewdt=0.5
c          print*, "halved due to damage incr", noel         
          end if
		  
      if (sigma_crack_plane(noel).ge.0.1) then !failed slipsystem is in tension
           d1 = d1 + 3.*(dtime/H_cycle) !increase damage to stiffness matrix
		   
      else if (sigma_crack_plane(noel).lt.0.1) then !failed slipsystem is in compression			
           d1 = d1 - 3.*(dtime/H_cycle) !decrease damage to stiffness matrix
      end if


c Place bounds on d1:
        if (d1 .gt. 0.99) then
          d1=0.99   ! range: 0 - 1
        else if (d1 .lt. 0.0) then
          d1=0.0    ! range: 0 - 1		
        end if
		
c        print*, 'element damage is: ', d1
	 
      end if !element is in crack
	  
#endif

c          print*, 'CHECK_B'
c--------------------------------------------------------------------
c  Initialize arrays for averaging over grains
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        sig_avg(i,j) = 0.0
       end do
      end do
      if (kinc.gt.4000) then
       print*, 'B'
      end if
c====================================================================
c  Begin loop over grains 
c====================================================================

      do m = 1,num_grains_UMAT

c--------------------------------------------------------------------
c  Rotate local anisotropic elasticity tensor to
c  global coordinates.
c--------------------------------------------------------------------

      call rotate_4th(dir_cos(1,1,m),C0,C)

c--------------------------------------------------------------------
c  Convert Miller Indices to global coordinates
c--------------------------------------------------------------------

      do n = 1,num_slip_sys
        do i = 1,3
          xs0(i,n) = 0.0
          xm0(i,n) = 0.0
          do j = 1,3
            xs0(i,n) = xs0(i,n) + dir_cos(i,j,m) * z(j,n)
            xm0(i,n) = xm0(i,n) + dir_cos(i,j,m) * y(j,n)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Initialize number of subincrements.  Note that the
c  subincrement initially equals the total increment.  This
c  remains the case unless the process starts to diverge.
c--------------------------------------------------------------------

      N_incr       = 1
      N_incr_total = 1

      if (kinc.gt.4000) then
       print*, 'C'
      end if
c====================================================================
c  Top of Subincrement Time Step Integration Loop.
c  Calculate subincrement time step size.
c====================================================================

  100 dt_incr = dtime / N_incr_total
c      print*,N_incr_total

c-------------------------------------------------------------------
c  Initialize ref shear stress
c-------------------------------------------------------------------

      do n = 1,num_slip_sys
        g(n,m) = g0(n,m)
      end do

c-------------------------------------------------------------------
c  Initialize back stress
c-------------------------------------------------------------------

      do n = 1,num_slip_sys
        a(n,m) = a0(n,m)
        a_1(n,m) = a0_1(n,m)
        a_2(n,m) = a0_2(n,m)
        a_3(n,m) = a0_3(n,m)
        a_4(n,m) = a0_4(n,m)		
      end do

c--------------------------------------------------------------------
c  Initialize deformation gradients for beginning and
c  end of subincrement.
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          F0(i,j) = dfgrd0(i,j) + (dfgrd1(i,j) - dfgrd0(i,j)) *
     &					(N_incr - 1) / N_incr_total
          F1(i,j) = dfgrd0(i,j) + (dfgrd1(i,j) - dfgrd0(i,j)) *
     &						N_incr / N_incr_total
        end do
      end do

c--------------------------------------------------------------------
c  Multiply F() by F_p_inv() to get F_el()
c--------------------------------------------------------------------

      call aa_dot_bb(3,F0,F_p_inv_0(1,1,m),F_el)
      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm.
c--------------------------------------------------------------------

      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate elastic Green Strain
c--------------------------------------------------------------------
c      print*, 'green strain'
      call transpose(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
      do i = 1,3
        E_el(i,i) = E_el(i,i) - 1
        do j = 1,3
          E_el(i,j) = E_el(i,j) / 2
c          print*, E_el(i,j)
        end do
      end do 

c--------------------------------------------------------------------
c  Multiply the anisotropic stiffness tensor by the Green strain 
c  to get the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------

      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from PK2 stress to Cauchy stress
c--------------------------------------------------------------------

      det = determinant(F_el)
      call transpose(3,F_el,array2)
      call aa_dot_bb(3,F_el,Spk2,array1)
      call aa_dot_bb(3,array1,array2,sig)
      do i = 1,3
       do j = 1,3
        sig(i,j) = sig(i,j) / det
       end do
      end do

c--------------------------------------------------------------------
c  Calculate stresses on each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        tau(k,m) = 0.0
        sigma(k) = 0.0
        do j = 1,3
          do i = 1,3
            tau(k,m) = tau(k,m) + xs(i,k) * xm(j,k) * sig(i,j)
            sigma(k) = sigma(k) + xm(i,k) * xm(j,k) * sig(i,j)
          end do
        end do
c        print*, tau(k,1)
      end do
c      print*, " "
      if (kinc.gt.4000) then
       print*, 'D'
      end if
c--------------------------------------------------------------------
c  Store the stress normal to the slip plane
c XXX FIXME Should this be done here or after increment converges?
c I would expect after increment converges if we are interested in accuracy, right?
c--------------------------------------------------------------------	  
	  
c          print*, 'CHECK_C'
c          do k = 1,num_slip_sys
c          sigma_gl(noel,k)=sigma(k) !into COMBLK array for FIP calculation
c          if (kstep.eq.1) then
c          print*, 'S_G_1:', sigma(k)		  
c          end if
c
c          end do
	
c          print*, 'CHECK_D'

c--------------------------------------------------------------------
c  Calculate 1st estimate of gamma_dot for each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
       gamma_dot(k,m)=gamma_dot_zero*power((tau(k,m)-a(k,m))/
     & g(k,m),flow_exp)
      end do
      
c--------------------------------------------------------------------
c  Calculate d(Tau)/d(Gamma_dot)
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
       do ib = 1,num_slip_sys
        dTaudGd(ia,ib) = 0.0

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ib) * xm0(j,ib)
         end do
        end do
        call aaaa_dot_dot_bb(3,C0,array1,array2)

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ia) * xm0(j,ia)
         end do
        end do
        call aa_dot_dot_bb(3,array1,array2,dTaudGd(ia,ib))

        dTaudGd(ia,ib) = dTaudGd(ia,ib) * dt_incr
       end do !ib
      end do !ia

#ifdef DEBUG
      print*
      print*, "Start NR - Time = ", time(2)
#endif


c          print*, 'CHECK_D1'

c====================================================================
c  Begin Newton-Raphson iterative loops.
c====================================================================

      converged = .false.

      do while (.not.converged)

      converged = .true.

c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, sse.
c--------------------------------------------------------------------

#ifdef DEBUG
      print*
      print*, "Before Eval_func 1 - Time = ", time(2)
#endif

      if (kinc.gt.4000) then
       print*, 'E'
      end if

      call eval_func(xs0,	dt_incr,	gamma_dot(1,m),	
     &		     xm0,	F1,		num_slip_sys,
     &		     C,		F_p_inv(1,1,m),	F_p_inv_0(1,1,m),
     &		     g0(1,m),	Hdir,		Hdyn,
     &		     a0(1,m),	Adir,		Adyn,
     &		     a0_1(1,m),	Adir_1,		Adyn_1,
     &		     a0_2(1,m),	Adir_2,		Adyn_2,	
     &		     a0_3(1,m),	Adir_3,		Adyn_3,
     &		     a0_4(1,m),	Adir_4,		Adyn_4,	 
     &		     g_sat,	F_el,		flow_exp,
     &		     sig,	tau(1,m),	gamma_dot_zero,
     &		     g(1,m),	func,		xL_p,
     &		     xLatent,	sse,		a(1,m),
     &		     a_1(1,m),     a_2(1,m),
     &		     a_3(1,m),     a_4(1,m), m_OW2)	 
	 
      sse_ref = sse

#ifdef DEBUG
      print*
      print*, "After Eval_func 1 - Time = ", time(2)
#endif

      if (kinc.gt.4000) then
       print*, 'F'
      end if

c         print*, 'CHECK_D2'
c--------------------------------------------------------------------
c  Begin calculation of the partial derivatives needed for 
c  the Newton-Raphson step!!!
c  Calculate derivative of the hardening variable, g-alpha,
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------

      sum = 0
      do ia = 1,num_slip_sys
        sum = sum + abs(gamma_dot(ia,m))
      end do
      do ia = 1,num_slip_sys
       do ib = 1,num_slip_sys
         temp = xLatent
         if (ia .eq. ib) temp = 1.0
         dgadgb(ia,ib) = (Hdir * temp - Hdyn*g(ia,m)) * dt_incr / 
     &			(1 + Hdyn * sum * dt_incr)
         if(gamma_dot(ib,m).lt.0.0)dgadgb(ia,ib)=-dgadgb(ia,ib)
       end do
      end do

c          print*, 'CHECK_D2.1'
c--------------------------------------------------------------------
c  Calculate derivative of kinematic stress, a-alpha
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
c       daadga(ia) = (Adir_1 - Adyn_1 * A_1(ia,m) * sgn(gamma_dot(ia,m)))
c     &       * dt_incr / (1 + Adyn_1 * dt_incr * abs(gamma_dot(ia,m)))
c     &       +      (Adir_2 - Adyn_2 * A_2(ia,m) * sgn(gamma_dot(ia,m)))
c     &       * dt_incr / (1 + Adyn_2 * dt_incr * abs(gamma_dot(ia,m)))
c     &       +      (Adir_3 - Adyn_3 * A_3(ia,m) * sgn(gamma_dot(ia,m)))
c     &       * dt_incr / (1 + Adyn_3* dt_incr * abs(gamma_dot(ia,m)))	
c     &       +      (Adir_4 - Adyn_4 * A_4(ia,m) * sgn(gamma_dot(ia,m)))
c     &       * dt_incr / (1 + Adyn_4 * dt_incr * abs(gamma_dot(ia,m)))


       daadga(ia) =
     & (((Adir_1-(Adir_1/b_1)*
     & (((A_1(ia,m)/b_1)**m_OW2)*A_1(ia,m)*sgn(gamma_dot(ia,m))))*dt_incr) 
     &  / (1+dt_incr*abs(gamma_dot(ia,m))*(Adir_1/b_1)*(
     &  A_1(ia,m)*(2./b_1)*((A_1(ia,m)/b_1)**(m_OW2-1))
     &	 + (A_1(ia,m)/b_1)**m_OW2  )))
     &  + (((Adir_2-(Adir_2/b_2)*
     & (((A_2(ia,m)/b_2)**m_OW2)*A_2(ia,m)*sgn(gamma_dot(ia,m))))*dt_incr) 
     &  / (1+dt_incr*abs(gamma_dot(ia,m))*(Adir_2/b_2)*(
     &  A_2(ia,m)*(2./b_2)*((A_2(ia,m)/b_2)**(m_OW2-1))
     &	 + (A_2(ia,m)/b_2)**m_OW2  )))
!     &  + (((Adir_3-(Adir_3/b_3)*
!     & (((A_3(ia,m)/b_3)**m_OW2)*A_3(ia,m)*sgn(gamma_dot(ia,m))))*dt_incr) 
!     &  / (1+dt_incr*abs(gamma_dot(ia,m))*(Adir_3/b_3)*(
!     &  A_3(ia,m)*(2./b_3)*((A_3(ia,m)/b_3)**(m_OW2-1))
!     &	 + (A_3(ia,m)/b_3)**m_OW2  )))
!     &  + (((Adir_4-(Adir_4/b_4)*
!     & (((A_4(ia,m)/b_4)**m_OW2)*A_4(ia,m)*sgn(gamma_dot(ia,m))))*dt_incr) 
!     &  / (1+dt_incr*abs(gamma_dot(ia,m))*(Adir_4/b_4)*(
!     &  A_4(ia,m)*(2./b_4)*((A_4(ia,m)/b_4)**(m_OW2-1))
!     &	 + (A_4(ia,m)/b_4)**m_OW2  )))	 
	 
      end do

c          print*, 'CHECK_D2.2'
c--------------------------------------------------------------------
c  Form "A-matrix" of derivatives wrt d_gamma_beta
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys
          array3(ia,ib) = dgadgb(ia,ib) * (tau(ia,m)-a(ia,m))/g(ia,m)
        end do
        array3(ia,ia) = array3(ia,ia) + 
     &		g(ia,m) / (flow_exp * gamma_dot_zero) *
     &		abs(power((tau(ia,m)-a(ia,m))/g(ia,m),(1.-flow_exp)))
      end do

c          print*, 'CHECK_D2.3'
c--------------------------------------------------------------------
c  Add d(Tau)/d(Gamma_dot) to the A-matrix for the Newton-
c  Raphson iteration.
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys
          array3(ia,ib) = array3(ia,ib) + dTaudGd(ia,ib)
        end do
      end do

c          print*, 'CHECK_D2.4'
c--------------------------------------------------------------------
c  Add d(a-alpha)/d(Gamma_dot) to the A-matrix for the Newton-
c  Raphson iteration.
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
        array3(ia,ia) = array3(ia,ia) + daadga(ia)
      end do

c          print*, 'CHECK_D2.5'
c--------------------------------------------------------------------
c  Calculate the gradient of sse wrt gamma_dot().  Will be used
c  later to ensure that line search is in the correct direction.
c--------------------------------------------------------------------

      do j = 1,num_slip_sys
        grad(j) = 0.0
        do i = 1,num_slip_sys
          grad(j) = grad(j) + func(i) * array3(i,j)
        end do
        grad(j) = 2 * grad(j)
      end do

c          print*, 'CHECK_D2.6'
c--------------------------------------------------------------------
c  Solve for increment of gamma_dot.  Solution is returned 
c  in the func() array.
c--------------------------------------------------------------------

      call simeq(num_slip_sys,array3,func)

c          print*, 'CHECK_D2.7'
c--------------------------------------------------------------------
c  Store offset in d_gamma_dot(k) 
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
         d_gamma_dot(k) = func(k)
      end do

c          print*, 'CHECK_D2.8'
c--------------------------------------------------------------------
c  Check to make sure that N-R step leads 'down hill' the 
c  sse surface.
c--------------------------------------------------------------------

      sum = 0.0
      do i = 1,num_slip_sys
        sum = sum - grad(i) * d_gamma_dot(i)
      end do

      if (sum .gt. 0.0) then
        do i = 1,num_slip_sys
          d_gamma_dot(i) = -d_gamma_dot(i)
        end do
      end if

c          print*, 'CHECK_D2.9'
c--------------------------------------------------------------------
c  Multiply step size by two 'cause next loop will divide it by 2.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        d_gamma_dot(k) = d_gamma_dot(k) * 2
      end do
      
c====================================================================
c  Begin line search.
c====================================================================

      improved = .false.

      do N_ctr = 1,max_loops

      sse_old = sse

c--------------------------------------------------------------------
c  Divide step size by two.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        d_gamma_dot(k) = d_gamma_dot(k) / 2
        gamma_try(k)   = gamma_dot(k,m) + d_gamma_dot(k)
      end do
      
c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, and sse based
c  on gamma_try(k)
c--------------------------------------------------------------------


#ifdef DEBUG
      print*
      print*, "Before Eval_func 2 - Time = ", time(2)
#endif


c          print*, 'CHECK_D3'
      call eval_func(xs0,	dt_incr,	gamma_try,	
     &		     xm0,	F1,		num_slip_sys,
     &		     C,		F_p_inv(1,1,m),	F_p_inv_0(1,1,m),
     &		     g0(1,m),	Hdir,		Hdyn,
     &		     a0(1,m),	Adir,		Adyn,
     &		     a0_1(1,m),	Adir_1,		Adyn_1,
     &		     a0_2(1,m),	Adir_2,		Adyn_2,	
     &		     a0_3(1,m),	Adir_3,		Adyn_3,
     &		     a0_4(1,m),	Adir_4,		Adyn_4,	 
     &		     g_sat,	F_el,		flow_exp,
     &		     sig,	tau(1,m),	gamma_dot_zero,
     &		     g(1,m),	func,		xL_p,
     &		     xLatent,	sse,		a(1,m),
     &		     a_1(1,m),     a_2(1,m),
     &		     a_3(1,m),     a_4(1,m), m_OW2)	


#ifdef DEBUG
      print*
      print*, "After Eval_func 2 - Time = ", time(2)
#endif

      if (kinc.gt.4000) then
       print*, 'G'
      end if

c          print*, 'CHECK_D4'
c--------------------------------------------------------------------
c  Check for 'convergence' of the line search.  Note that the line
c  search doesn't waste time converging closely.  This is 'cause
c  the N-R step does so much better.
c--------------------------------------------------------------------

      if ((sse_old.le.sse_ref).and.(sse.ge.sse_old).and.
     &					(N_ctr.gt.1)) improved=.true.

      if (improved) go to 200

      end do ! Linear Search

c--------------------------------------------------------------------
c  Add "d_gamma_dot" to gamma_dot to get new values for
c  this iteration. 
c--------------------------------------------------------------------

  200 do k = 1,num_slip_sys
        gamma_dot(k,m) = gamma_dot(k,m) + d_gamma_dot(k) * 2.0
      end do
      if (kinc.gt.4000) then
       print*, 'H'
      end if
c--------------------------------------------------------------------
c  If (sse_old > tolerance) then this step has not converged.
c--------------------------------------------------------------------

      if (sse_old .gt. tolerance) converged = .false.

c--------------------------------------------------------------------
c  If (sse_old > sse_ref/2) then convergence is too slow and
c  increment is divided into two subincrements.
c--------------------------------------------------------------------
 
        if ((sse_old.gt.sse_ref/2.0).and.(.not.converged)) then
        N_incr = 2 * N_incr - 1
        N_incr_total = 2 * N_incr_total

c        if (N_incr_total.gt.1000)  go to 500
           if (N_incr_total.gt.5e6) then
           pnewdt=0.5
           print*, 'increment halved due to convegence'
           return           
            end if
        go to 100
        end if
 
c--------------------------------------------------------------------
c  End iterative loop.
c--------------------------------------------------------------------

      end do ! 'Newton Raphson Iterative Loop'

c--------------------------------------------------------------------
c  If another subincrement remains to be done, then reinitialize
c  F_p_inv_0 and g0.  F0 and F1 gets reinitialized back at the
c  top of this loop.
c--------------------------------------------------------------------

      if (N_incr .lt. N_incr_total) then
        if (N_incr .eq. (N_incr/2)*2) then	! N_incr is 'even'
          N_incr = N_incr / 2 + 1
          N_incr_total = N_incr_total / 2
        else					! N_incr is 'odd'
          N_incr = N_incr + 1
        end if
        do i = 1,3
          do j = 1,3
            F_p_inv_0(i,j,m) = F_p_inv(i,j,m)
          end do
        end do
        do i = 1,num_slip_sys
          g0(i,m) = g(i,m)
          a0(i,m) = a(i,m)
          a0_1(i,m) = a_1(i,m)
          a0_2(i,m) = a_2(i,m)
          a0_3(i,m) = a_3(i,m)
          a0_4(i,m) = a_4(i,m)		  
        end do
        go to 100
      end if
c--------------------------------------------------------------------
c  Calculate the average Cauchy stress.
c--------------------------------------------------------------------
c      print*, 'CHECK_E'
      do i = 1,3
        do j = 1,3
          sig_avg(i,j) = sig_avg(i,j) + sig(i,j)
        end do
      end do

c--------------------------------------------------------------------
c  Write out Euler Angles in Kocks Convention.
c--------------------------------------------------------------------

      call aa_dot_bb(3,F_el,dir_cos(1,1,m),array1)
      call kocks_angles(npt,m,time(2),array1)

c--------------------------------------------------------------------
c  Calculate the average elasticity tensor.
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          if (m .eq. 1) C_avg(i,j,k,l) = 0.0
          C_avg(i,j,k,l) = C_avg(i,j,k,l) + C(i,j,k,l)
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm.
c--------------------------------------------------------------------

      call inverse_3x3(F_el,F_el_inv)
      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate the derivative of the plastic part of the rate of
c  deformation tensor in the current configuration wrt sigma.
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          if (m .eq. 1) ddpdsig(i,j,k,l) = 0.0
          do n = 1,num_slip_sys
           ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) + (xs(i,n)*xm(j,n) +
     &       xm(i,n) * xs(j,n)) * (xs(k,n) * xm(l,n) + xm(k,n) *
     &       xs(l,n)) * abs(power((tau(n,m)-a(n,m))/g(n,m),
     &		flow_exp-1.0)) / g(n,m)
          end do ! num_slip_sys
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c      write(7,'(12f7.3)')(gamma_dot(i,m)/0.0004,i=1,12)

c--------------------------------------------------------------------
c  End loop over all the grains
c--------------------------------------------------------------------

      end do ! m  = 1,num_grains_UMAT
c--------------------------------------------------------------------
c  Calculate Green Strain
c--------------------------------------------------------------------
      call transpose(3,F1,array1)
	call aa_dot_bb(3,array1,F1,E_tot)

      do i = 1,3
        E_tot(i,i) = E_tot(i,i) - 1
        do j = 1,3
          E_tot(i,j) = E_tot(i,j) / 2
        end do 
      end do
	  
c--------------------------------------------------------------------
c  Store sig_avg() into sig()
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          sig(i,j) = sig_avg(i,j) / num_grains_UMAT
        end do
      end do
	  
c====================================================================
c  Begin calculation of the Jacobian (the tangent stiffness matrix).
c====================================================================

c--------------------------------------------------------------------
c  Divide C_avg by num_grains_UMAT to get correct average.
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          C_avg(i,j,k,l) = C_avg(i,j,k,l) / num_grains_UMAT
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Output stress, strain, etc. for post processing.
c--------------------------------------------------------------------

      sum = 0.0
      do i = 1,3
       do j = 1,3
        sum = sum + sig(i,j) * sig(i,j)
       end do
      end do
      trace = sig(1,1) + sig(2,2) + sig(3,3)
      sig_eff=sqrt(1.5 * sum - 0.5 * trace**2)

c      write(7,'(a3,i5,5f14.6)')'BM1',ninc,time(2),
c     &	abs(log(dfgrd1(3,3)))+dfgrd1(1,3)/dfgrd1(3,3)/1.7320508,
c     &	sig_eff,sig(3,3),sig(1,3)

c      write(7,'(a3,i5,5f14.6)')'BM1',ninc,time(2),sig(1,2),sig(2,3)
c--------------------------------------------------------------------
c  Calculate the inverse of F_el
c--------------------------------------------------------------------

      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Scale by appropriate constants and divide by num_grains_UMAT to
c  get the average value.  ALSO multiply by 'dtime' which is
c  d(sig)/d(sig_dot).
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) * dtime * flow_exp *
     &		gamma_dot_zero / 4. / num_grains_UMAT
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Multiply the 4th rank elastic stiffness tensor by the derivative
c  of the plastic part of the rate of deformation tensor wrt sig_dot.
c--------------------------------------------------------------------

      call aaaa_dot_dot_bbbb(3,C_avg,ddpdsig,array6)

c--------------------------------------------------------------------
c  Add 4th rank identity tensor to array6()
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          array6(i,j,k,l) = array6(i,j,k,l) + 0.5 * 
     &		(del(i,k) * del(j,l) + del(i,l) * del(j,k))
         end do
        end do
       end do
      end do
      if (kinc.gt.4000) then
       print*, 'I'
      end if
c--------------------------------------------------------------------
c  Need to take the inverse of Array4.  Since it relates two 2nd
c  rank tensors that are both symmetric, Array4 can be xformed to 
c  Voigt notation style in order to do the inverse, then xformed back.
c--------------------------------------------------------------------

      call forth_to_Voigt (Array6,Array4)
      call inverse        (6,Array4,Array5)
      call Voigt_to_forth (Array5,Array6)

c--------------------------------------------------------------------
c  Multiply Array6 by C, the elastic stiffness matrix to
c  finally get the Jacobian, but as a 4th rank tensor.
c--------------------------------------------------------------------

      call aaaa_dot_dot_bbbb (3,Array6,C_avg,ddsdde_4th)


c--------------------------------------------------------------------
c  Store the Jacobian in Voigt notation form.
c--------------------------------------------------------------------

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=i+j+1
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=k+l+1
          array4(ia,ib) = ddsdde_4th(i,j,k,l)
          IF(IB.GE.4) ARRAY4(IA,IB) = 2 * ARRAY4(IA,IB)
         end do
        end do
       end do
      end do

c      call print_array(6,array4)

      do i =1,6
        do j = 1,6
          ddsdde(i,j) = 0.0
        end do
      end do

      if (ndi .eq. 1) then			! 1-D
         ddsdde(1,1) = array4(1,1)
      else if (ndi .eq. 2) then			! 2-D plane stress & axi
         do i = 1,2
            do j = 1,2
               ddsdde(i,j) = array4(i,j)
            end do
         end do
         ddsdde(1,3) = array4(1,4)
         ddsdde(2,3) = array4(2,4)
         ddsdde(3,1) = array4(4,1)
         ddsdde(3,2) = array4(4,2)
         ddsdde(3,3) = array4(4,4)
      else if (ndi .eq. 3 .and. nshr .eq. 1) then ! plane strain
         do i = 1,4
            do j = 1,4
               ddsdde(i,j) = array4(i,j)
            end do
         end do
      else					! Full 3-D
         do i = 1,6
            do j = 1,6
               ddsdde(i,j) = array4(i,j)
            end do
         end do
      end if
	  
	  


c--------------------------------------------------------------------
c  Store the stress tensor in the ABAQUS stress 'vector'
c--------------------------------------------------------------------

      do i = 1,ndi
         stress(i) = sig(i,i)
      end do
      if (nshr .eq. 1) stress(ndi+1) = sig(1,2)
      if (nshr .eq. 3) then
         stress(4) = sig(1,2)
         stress(5) = sig(1,3)
         stress(6) = sig(2,3)
      end if

	  
#ifdef DEBUG
      print*
      print*, "Save statev - Time = ", time(2)
#endif


c      print*, 'CHECK_G'
c-------------------------------------------------------------------
c  Store the internal variables in the statev() array
c-------------------------------------------------------------------
         n = 0
          
c-------------------------------------------------------------------
c  Store the Euler Angles  1-3
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,3
             n = n + 1
             statev(n) = psi(i,m)
           end do
         end do
		 
c        print*, "n-chk-1b", n		 
            
c-------------------------------------------------------------------
c  Store the plastic part of F   4-12
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
            do j = 1,3
               do i = 1,3
                  n = n + 1
                  statev(n) = F_p_inv(i,j,m)
               end do
            end do
         end do

c-------------------------------------------------------------------
c  Store the reference shear stresses.     13-24
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = g(i,m)
           end do
         end do
		 
c        print*, "n-chk-2b", n			 

c-------------------------------------------------------------------
c  Store the back stress and its 4 components.     25-84
c-------------------------------------------------------------------

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = a(i,m)
           end do
         end do

        do i = 1,num_slip_sys
           a_elem(noel,i) = a(i,1)
        end do

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = a_1(i,m)
           end do
         end do

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = a_2(i,m)
           end do
         end do		 
		 
         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = a_3(i,m)
           end do
         end do

         do m = 1,num_grains_UMAT
           do i = 1,num_slip_sys
             n = n + 1
             statev(n) = a_4(i,m)
           end do
         end do		

c       print*, "n-chk-3b", n			 
      if (kinc.gt.4000) then
       print*, 'J'
      end if		 
c-------------------------------------------------------------------
c  Plastic Strain Calculations
c-------------------------------------------------------------------
c  Increment of plastic shear strain accumulated on each slip system
c  over this time step.
      do i = 1,num_slip_sys
       delta_gamma(i) = gamma_dot(i,1)*dtime
      end do

      do j = 1,3
       do k = 1,3
        delta_E_p(j,k) =  0.0
       end do
      end do

c  Increment of the plastic strain tensor
      do j = 1,3
       do k = 1,3
        do l = 1,num_slip_sys
	      delta_E_p(j,k) = delta_E_p(j,k) + 0.5*delta_gamma(l)*
     &     (xs0(j,l)*xm0(k,l)+xs0(k,l)*xm0(j,l))
          end do
       end do
      end do 

c  Plastic strain tensor
      do i = 1,3
      do j = 1,3
       E_p(i,j) = E_p(i,j)+delta_E_p(i,j)
      end do
      end do 

c Store as SDV 85-93:	  
      do i = 1,3
       do j = 1,3
        n = n+1
        statev(n) = E_p(i,j)
       end do
      end do		      
 
c  Effective total strain
      call aa_dot_dot_bb(3,E_tot,E_tot,sum)
      E_eff = sqrt(2./3.*sum)

c  Store as SDV 94:      
      n = n+1
      statev(n)= E_eff
c      print*, "n-chk-3.5b", n


c  Effective Plastic Strain
      call aa_dot_dot_bb(3,E_p,E_p,sum1)
      E_p_eff = sqrt(2./3.*sum1)

      if (E_p_eff.lt.0) E_p_eff = 0.
		
c  Store as SDV 95:    		
      n = n+1
      statev(n)= E_p_eff


c  Effective plastic strain increment
      sum2 = 0.
      call aa_dot_dot_bb(3,delta_E_p,delta_E_p,sum2)
      sum2=sqrt(2./3.*sum2)
      E_p_eff_cum = sum2 + E_p_eff_cum

c  Store as SDV 96:    
      n = n+1       
      statev(n) = E_p_eff_cum
	  
c      print*, "n-chk-4b", n		
		


c  Plastic shear strain accumulated on each slip system
c  over the entire time SDV 97-108

c       print*, 'number of slip systems', num_slip_sys  
       do i = 1,num_slip_sys
         delta_gamma_cum(i) = delta_gamma_cum(i)+delta_gamma(i)
         n=n+1
		  
         statev(n) = delta_gamma_cum(i)

       end do

c        print*, "n-chk-4.5b", n
	   


         do i = 1,num_slip_sys
          gamma_cum_element(noel,i) = delta_gamma_cum(i)!store in COMBLK for FIP calc in UEXDB
		  
          if (delta_gamma_cum(i).gt.delta_gamma_cum_max(i) ) then
           delta_gamma_cum_max(i)=delta_gamma_cum(i)
		   
          else if (delta_gamma_cum(i).lt.delta_gamma_cum_min(i)) then
           delta_gamma_cum_min(i)=delta_gamma_cum(i)
		   
          end if
	 
        end do


         do i = 1,num_slip_sys ! 109-132
         n=n+1
         statev(n)= delta_gamma_cum_max(i)
         n=n+1
         statev(n)=delta_gamma_cum_min(i)
         end do
		 
c        print*, "n-chk-4.6b", n		 

c remnants of Gustavo's old FIP calculation done in UMAT
c might be a useful check if it is working properly
c which i dont think it currently is		
         n=n+1
         do i = 1,num_slip_sys
         FIP=0.5*(delta_gamma_cum_max(i)-delta_gamma_cum_min(i))
     &        *(1.+0.5*sigma(i)/512.)
         FIP_elem(noel,i)=FIP
         if  (FIP .gt. statev(n)) then
         statev(n)=FIP   !133
         FIP_elem_max(noel,2)=i
         FIP_elem_max(noel,3)=sigma(i)
         end if
         end do
		 

c        print*, "n-chk-5b", n		

         FIP_elem_max(noel,1)=statev(n) !133

         n=n+1		
         statev(n)= d1 !134
		 
c        print*, "n-chk-6b", n			 


c--------------------------------------------------------------------
c  Re-calculate stresses normal to slip system at end of inc
c  and store in COMMONBLOCK for FIP calculation
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        sigma(k) = 0.0
        do j = 1,3
          do i = 1,3
            sigma(k) = sigma(k) + xm(i,k) * xm(j,k) * sig(i,j)
          end do
        end do
		
      sigma_gl(noel,k)=sigma(k) !into COMBLK array for FIP calculation	
	  
          if (kstep.eq.1) then
c          print*, 'S_G_2:', sigma(k)		  
          end if
		  
      end do

		   
c--------------------------------------------------------------------
c  Recalc stress normal to the crack plane and store as SDV 150
c--------------------------------------------------------------------
        ! Reset to zero
        sigma_crack_plane(noel) = 0.0   

        ! Calc normal component of stress
        do j = 1,3
         do i = 1,3
          sigma_crack_plane(noel) = sigma_crack_plane(noel) +      
     &    Crack_plane_norm(noel,i) * Crack_plane_norm(noel,j) * sig(i,j)
         end do
        end do
        
        ! Store as SDV 150         
        statev(150) = sigma_crack_plane(noel)
		   
c-------------------------------------------------------------------
c  Store the FIPs on each slip system and the total FIP 134-145 and 146
c-------------------------------------------------------------------
c          print*, "n-chk-6.5b", n			   
           do i = 1,num_slip_sys
           n = n + 1		   
           statev(n) = FIP_total(noel,i)
           end do
		   
           statev(n) = FIP_total(noel,13) !SDV146
		   
c          print*, "n-chk-7b", n	


c-------------------------------------------------------------------
c  Store the Cauchy Stress tensor into COMBLK array for use in UEXDB
c-------------------------------------------------------------------

      do j = 1,3
      do k = 1,3			 
      sig_elem(noel,j,k) = sig(j,k)
      end do
      end do 
 
c-------------------------------------------------------------------
c  Store the Total Strain tensor into COMBLK array for use in UEXDB
c-------------------------------------------------------------------	

      do j = 1,3
      do k = 1,3
         E_tot_elem(noel,j,k) = E_tot(j,k)
      end do
      end do 
	  
c-------------------------------------------------------------------
c  Store the Plastic Strain tensor into COMBLK array for use in UEXDB
c-------------------------------------------------------------------			 

      do j = 1,3
       do k = 1,3
         E_p_elem(noel,j,k) = E_p(j,k)
       end do
      end do 
			 	   
      if (kinc.gt.4000) then
       print*, 'ND'
      end if

      end if ! if phase is elastic-plastic
	 
      if(phase.eq.2) then
      d1 = 0
c-----------------------------------------------------
c Apply damage and store variables if damaged element and damage is active
c-----------------------------------------------------
#ifdef DAMAGE
      d1= statev(134)
      sigma_crack_plane(noel) = statev(150)
      if (dam_elem(noel).eq.1.0) then !element is in crack
	  
          if (dtime.gt.H_cycle/100) then
          pnewdt=0.5
          print*, "halved due to damage incr", noel         
          end if
		  
      if (sigma_crack_plane(noel).ge.0.1) then !failed slipsystem is in tension
           d1 = d1 + 3.*(dtime/H_cycle) !increase damage to stiffness matrix
		   
      else if (sigma_crack_plane(noel).lt.0.1) then !failed slipsystem is in compression			
           d1 = d1 - 3.*(dtime/H_cycle) !decrease damage to stiffness matrix
      end if

c Place bounds on d1:
        if (d1 .gt. 0.99) then
          d1=0.99   ! range: 0 - 1
        else if (d1 .lt. 0.0) then
          d1=0.0    ! range: 0 - 1		
        end if
		
c        print*, 'element damage is: ', d1
      statev(134) = d1
c      print*, 'damage'
c      print*, statev(134)
      end if !element is in crack
      
#endif
       elastic_modulus      = 161000.0*(1-d1)							! The Elastic Modulus
       poisson_ratio     = 0.3											! The Poisson's Ratio
       bulk_modulus      = elastic_modulus /
     &  (3.0d0 * (1.0d0 - 2.0d0 * poisson_ratio))						! The Bulk Modulus
       shear_modulus      = elastic_modulus
     &	/ (2.0d0 * (1.0d0 + poisson_ratio))								! Shear Modulus
       lame_parameter = (elastic_modulus * poisson_ratio) /
     &  ((1.0d0 + poisson_ratio) * (1.0d0 - 2.0d0 * poisson_ratio))		! Lame's 1st Parameter
	   
! The Elastic Stiffness
! If we consider the 6x6 Stiffness Matrix:
!               | L+2G   L     L   0   0   0 |
!               |  L    L+2G   L   0   0   0 |
! [Stress] =    |  L     L    L+2G 0   0   0 | [Strain]
!               |  0     0     0   G   0   0 |
!               |  0     0     0   0   G   0 |
!               |  0     0     0   0   0   G |
       Do i = 1, ndi
       Do j = 1, ndi
		     ddsdde(j, i) = lame_parameter
          Enddo
          ddsdde(i, i) = lame_parameter + (2.0d0 * shear_modulus)
       Enddo
       Do i = ndi + 1, ntens
          ddsdde(i, i) = shear_modulus
       Enddo
	   
! Calculate the Stress
! Stress(new) = Stress(old) + dStress/dStrain * dStrain
       Do j = 1, ntens
          stress(j) = 0
       end do
       Do i = 1, ntens
          Do j = 1, ntens
             stress(j) = stress(j) + 
     & (ddsdde(j, i) * (strain(i) + dstrain(i)))
          Enddo
       Enddo
      sigma_crack_plane(noel) = 0.0   
        sig(1,1) = stress(1)
        sig(2,2) = stress(2)
        sig(3,3) = stress(3)
        sig(1,2) = stress(4)
        sig(2,1) = stress(4)
        sig(1,3) = stress(5)
        sig(3,1) = stress(5)
        sig(2,3) = stress(6)
        sig(3,2) = stress(6)
        ! Calc normal component of stress
        do j = 1,3
         do i = 1,3
          sigma_crack_plane(noel) = sigma_crack_plane(noel) +      
     &    Crack_plane_norm(noel,i) * Crack_plane_norm(noel,j) * sig(i,j)
         end do
        end do
        
        ! Store as SDV 150         
        statev(150) = sigma_crack_plane(noel)
c-------------------------------------------------------------------
c  Store the Cauchy Stress tensor into COMBLK array for use in UEXDB
c-------------------------------------------------------------------		
        do j = 1,3
        do k = 1,3			 
        sig_elem(noel,j,k) = sig(j,k)
        end do
        end do 
		
         E_tot_elem(noel,1,1) = strain(1)
         E_tot_elem(noel,2,2) = strain(2)
         E_tot_elem(noel,3,3) = strain(3)
         E_tot_elem(noel,1,2) = strain(4)
         E_tot_elem(noel,2,1) = strain(4)
         E_tot_elem(noel,1,3) = strain(5)
         E_tot_elem(noel,3,1) = strain(5)
         E_tot_elem(noel,2,3) = strain(6)
         E_tot_elem(noel,3,2) = strain(6)
		
      end if
	  
      return
      end
	  
	  
c********************************************************************
c********************************************************************
c  Begin UVARM Subroutine
c********************************************************************
c********************************************************************
      
         SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(200)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(200),JARRAY(200),JMAC(*),JMATYP(*),COORD(*)

      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)

c	Plastic strain

      UVAR(1)  = ARRAY(1)
      UVAR(2)  = ARRAY(2)
      UVAR(3)  = ARRAY(3)
      RETURN
      END

	  
 
c********************************************************************
c********************************************************************
c  Begin UEXTERNALDB Subroutine
c********************************************************************
c********************************************************************
c WRITTEN BY: Gustavo Castelluccio - 2012
C MODIFIED BY: Conor Hennessey - 2014-2015
C MODIFIED BY: Paul Kern 2015

c Tags used to keep track of sections that need revision or work:
c FIXME to mark potential problematic code that requires special attention and/or review
c NOTE to document inner workings of code and indicate potential pitfalls
c TODO to indicate planned enhancements
c XXX to warn other programmers of problematic or misguiding code

       SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
	  
      logical eval_MSC_k
      real*8 mesh_size !make sure this is a float and not an integer so crack length correct
#include "Geom_Def.txt" !parameters such as number of grains, elements, etc.
	 
       DIMENSION TIME(2),
     &  FIP_cycle(num_elem,num_slip_sys), ! FIP values for each element (evaluated for each slip system)
     &  Min_dist_temp(num_lines_Min_dist*6), 
     &  Grains_info_temp(num_grains*4),
     &  gamma_ratch(num_elem,num_slip_sys),
     &  FIP_max(2),
     &  FIP_max_glss(10,4), ! Ten max FIPS for cycle (ordered highest to lowest), (FIP value, grain, layer, slip system)
     &  sig_elem_avg(3,3), !Plotting variable
     &  E_tot_elem_avg(3,3), !Plotting variable
     &  E_p_elem_avg(3,3), !Plotting variable
     &  list_cracked_elem(num_elem_fail,4) !list of all cracked elem


     
        character*512 OUTDIR
        character*511 dmkname
        character*250 FNAMER1,FNAMER2,FNAMER3,FNAMER4,FNAMER5
        character*250 CRCK_ELEM, PYMAIN, ARREST_FILE	
        character*250 FNAMEW1,FNAMEW2,FNAMEW3,FNAMEW4
        character*250 FNAMEW5,FNAMEW6,FNAMEW7,FNAMEW8
        character*250 CSSC_FILE !file to save the cyclic stress strain data		
        character*250 B_FILE !file to save the cyclic stress strain data
        character*20 filenames,filenames1
        CHARACTER(len=2)::filenum
        logical flag1, Eval_Nuc_logical, Eval_MSC_logical, pre_step
        logical reset_shear_strain, calculate_FIP, EXT
        logical isMax_E_p_in, isMin_E_p_in, isMin_E_p_in_num
        logical isMax_E_p_in_num
        real*8 C_int,d_grMSC,FIP_vec_band_norm,FIP_dam_temp,CGR_msc_vec_band_norm,CGR_msc_temp
        INTEGER*8 N_El_GB_temp,Nel_temp,N_msc_temp,N_msc_vec_band_norm
        INTEGER*8 t_start, t_end, line_counter	
        real*8 list_cracked_elem !list of all cracked elem
        INTEGER*8 n_steps_per_cycle
		 

#include "Common_block_Alv02.txt"	 
#include "Definitions.txt"

      print*
      print*, 'UEXTERNALDB Called'	  
      print*,'STEP:',kstep," INC:", kinc
      print*,'LOP:', LOP 	  

	  !Time how long the execution of UEXDB takes:
      call system_clock(t_start)	  
	  
!#define DEBUG_UEXTER
  
c===============================================================
c   DATA INPUT TASKS - initialize variables, read in files, etc.
c===============================================================

      n_steps_per_cycle = 3
      angle_thresh=20*3.14159/180 ! 20 degrees in radians	   
	   
       CALL GETOUTDIR(OUTDIR, LENOUTDIR ) ! Abaqus utility routine to get output directory of current job and length of string
	   
c---------------------------------------------------------------------------------
c Initializations of variables and arrays
c---------------------------------------------------------------------------------	   
        if ((time(2).eq.0.0).and.(kstep.eq.0)) then
	
         print*, 'INITALIZE UEXTERNALDB VARs'
		 
         

         UEXTERNALDB_time = 0
         call system_clock(t_start)	  
         START_time = t_start
c         print*, START_time

          crack_arrested=.false.		 
          Nfailed=0	!number of grains that have failed	

		  

          do k=1,num_grains
           do i=1,max_num_layers
            do j=1,4
             Cnum(k,i,j)= 0.0
            end do
           end do
          end do  		  
		  
          do k=1,num_grains
          do kk=1,num_grains
	     disAngle(k,kk)=100
	     disAngleMax(k,kk)=0.0
           end do
          end do

          do k=1,num_grains
          do kk=1,num_grains
          do i=1,max_num_layers
          do ii=1,max_num_layers
          do j=1,4
          do jj=1,4
         Min_dist(k,kk,i,ii,j,jj)=.FALSE.
          end do
          end do
          end do
          end do
          end do
          end do		  

      do i=1,num_lines_Min_dist
      do k = 1,6	  
      Min_dist_temp(i*k) = 0.0
      end do
      end do
		  

		
        do i=1,num_grain_fail+1

        Cnum_fail(i)=num_elem

        end do		
		
         do i=1,num_elem
           do j=1,num_slip_sys
             FIP_cycle(i,j)=0.0
             gamma_cum_element_max(i,j)=0.0
             gamma_cum_element_min(i,j)=0.0
             gamma_cum_init_even_step(i,j)=0.0 !added by CDH			 
           end do
           end do
		   
       do k=1,num_elem
        do j=1,13
          FIP_total(k,j)=0.0
        end do
       end do		  
		  
       do i=1,num_elem
        dam_elem(i)=0.
        sigma_crack_plane(i) = 0.0		
          do j=1,3
           Crack_plane_norm(i,j) = 0.0 !contains the crack normal within an element
          end do		
       end do
      
c---------------------------------------------------------------------------------
c Read in grain neighbour list
c---------------------------------------------------------------------------------
 	   
       FNAMER1=dmkname('Neighbor_grains.txt',OUTDIR(1:LENOUTDIR),' ')
       OPEN(60,FILE=FNAMER1,STATUS='OLD',ACTION='READ')
       READ(60,*,END=70) Neighbor1   !Jumps to 70 when the last line read 
70     CLOSE(60)	   

       jump=0
       do i=1,num_grains
       Neighbor_num(i)= Neighbor1(i+jump) ! Neighbor_num stores the quantity of Neighbouring grains around grain i. 
       jump=jump+Neighbor_num(i)
       end do

	   
	   !TODO: combine this loop and previous for efficiency
       jump=0
       do i=1,num_grains
       do j=1,Neighbor_num(i)
       Neighbor_grains(i,j)=Neighbor1(i+j+jump)  !Neighbor_grains(i,j) stores the number of Neighbouring grains (j) around grain i. 
       end do
       jump=jump+Neighbor_num(i)
       end do

       do k=1,num_grains	   
         do i=1,Neighbor_num(k)
           print*, Neighbor_grains(k,i)
         end do
       end do		

       call READCRACKEDELEMENTS(list_cracked_elem)	   
	   
       print*, 'Finished time = 0 initializations'
       end if ! end time=0.0
	   
	   

c-------------------------------------------------------	   
c Read in Min_dist and assign booleans to Min_dist array:
c-------------------------------------------------------

        if ((kinc.eq.2).and.(kstep.eq.1).and.(LOP.eq.1)) then
c	    NOTE: added LOP=1 because of Illegal Mem ref during lop=2

        FNAMER5=dmkname('Min_dist.txt',OUTDIR(1:LENOUTDIR),' ')
        OPEN(73,FILE=FNAMER5,STATUS='OLD',ACTION='READ')
        READ(73,*,END=71) Min_dist_temp !Jumps to 71 when the last line read
71        close(73)

         print*, 'chck1'

         do i=1,num_lines_Min_dist
c         print*, 'chck1.1', i		 
         Min_dist(Min_dist_temp(i*6-5),Min_dist_temp(i*6-4),Min_dist_temp(i*6-3),
     &   Min_dist_temp(i*6-2),Min_dist_temp(i*6-1),Min_dist_temp(i*6) )=.TRUE.
         end do
c         print*, 'chck1.5'
		 
          do k=1,num_grains
          do kk=1,num_grains
          do i=1,max_num_layers
          do ii=1,max_num_layers
          do j=1,4
          do jj=1,4
         Min_dist(kk,k,ii,i,jj,j)=Min_dist(k,kk,i,ii,j,jj)
          end do
          end do
          end do
          end do
          end do
          end do
		  
         print*, 'chck2'		  

c-------------------------------------------------------
c Read in information about grains, assign misorientation info to Grs() array:	
c-------------------------------------------------------
	  
        FNAMER1=dmkname('Grains.txt',OUTDIR(1:LENOUTDIR),' ')
        OPEN(74,FILE=FNAMER1,STATUS='OLD',ACTION='READ')
        READ(74,*,END=72) Grains_info_temp !Reads the next line, but to jump to statement 72 when the last line has already been read. 
72      CLOSE(74)

          do k=1,num_grains
         Grs(k,1)= Grains_info_temp(k*4-3)
         Grs(k,2)= Grains_info_temp(k*4-2)
         Grs(k,3)= Grains_info_temp(k*4-1)
         Grs(k,4)= Grains_info_temp(k*4)
          end do	

c         print*, 'chck3'			  

c---------------------------------------------------------------------------------
c Calculate the misorientation among grains using the original angle distributions
c---------------------------------------------------------------------------------

         do k=1,num_grains
         do lll=1,Neighbor_num(k)
          k_nd=Neighbor_grains(k,lll)   ! 2nd grain, neighboiring of the potential grain to fail that is neighbor of one that already failed
      disAngle_temp = 100.0
      disAngle_tempMax = 0.0
      call Gr_misor(Grs(k,2),Grs(k,3),Grs(k,4),Grs(k_nd,2),Grs(k_nd,3),Grs(k_nd,4),disAngle_temp)

       disAngle(k,k_nd)=disAngle_temp 
c       print*, 'disAngle:', disAngle(k,k_nd)	   
       disAngleMax(k,k_nd)=disAngle_tempMax
       end do
       end do

       FNAMEW6=dmkname('disAngle.txt',OUTDIR(1:LENOUTDIR),' ')
       OPEN(59,FILE=FNAMEW6, status="unknown")
        do k=1,num_grains
        do lll=1,Neighbor_num(k)
           k_nd=Neighbor_grains(k,lll) 
        if (disAngle(k,k_nd).lt.angle_thresh) then
        d_gr_nd_temp=(1-disAngle(k,k_nd)/angle_thresh)*Grs(k_nd,1)
        else
        d_gr_nd_temp=0
        end if 
        write(59,'(i6,i6,F,F,F)') k, k_nd, disAngle(k,k_nd),
     &   d_gr_nd_temp, disAngleMax(k,k_nd)
        end do
        end do
        close(59)
	  
       end if ! ((kinc.eq.2).and.(kstep.eq.1))

c---------------------------------------------------------------------------------
c Add up all the elements in a given band to determine Cnum:
c Must be done after 1st increment because Elem_pos is initialized in UMAT
c---------------------------------------------------------------------------------
	   
      if ((kinc.eq.1).and.(kstep.eq.1).and.(LOP.eq.2)) then
        do i=1,num_elem
c        print*, 'Cnum initialized for element:',i
        do j=1,4
        Cnum(Elem_pos(i,5),Elem_pos(i,j),j) =
     & Cnum(Elem_pos(i,5),Elem_pos(i,j),j)+1.0
        end do
c       print*,"Elem_pos",Elem_pos(i,1),Elem_pos(i,2),Elem_pos(i,3),Elem_pos(i,4),Elem_pos(i,5)		
        end do
      end if

	  
c---------------------------------------------------------------------------------
c Flow control: determine which loops should run this UEXTERNALDB call
c---------------------------------------------------------------------------------

c Criteria for call occurring before beginning the next loading step
      if ((kinc.eq.1).and.(LOP.eq.1))then 
                   pre_step = .true.
                   print*, "pre_step is true"
                   print*				   
      else
                   pre_step = .false.
      end if

c NOTE if an increment starts to diverge, LOP1 may occur multiple times for a given increment,
c and bad data may be read from COMMONBLOCK!

c Criteria to enter life evaluation loops:	 
 
      Eval_Nuc_logical = .false.
      Eval_MSC_logical = .false.

      if ((pre_step).and.(.not.crack_arrested)) then
      j = SIZE(eval_life_array)	  
      print*,"eval_life_array size:",j 
	  
          do i=1,j
          
          print*,i
          print*,eval_life_array(1)
          print*,eval_life_array(i)

              if (kstep.eq.eval_life_array(1)) then !Nucleation evaluation step
                 Eval_Nuc_logical = .true.
                 print*, "Eval_nuc is true"
                 exit
                 
               else if (kstep.eq.eval_life_array(i)) then !MSC evaluation step
                 Eval_MSC_logical = .true.
                 print*, 'Eval MSC life for kstep:' , kstep
                 exit
                 
               end if
               
          end do   
      
      end if

c Criteria to reset extreme values of shear strain following completion of a cycle:
      if ((kstep.gt.1).and.(mod(kstep,n_steps_per_cycle)
     &    .eq.0).and.pre_step) then
                   reset_shear_strain = .true.
      else
                   reset_shear_strain = .false.
      end if

c Criteria to enter FIP calculation:
c Must be performed every increment b/c we record the change in plastic strain
c along each slip system when this code block is evaluated
      if (((kstep.eq.1).and.(KINC.gt.2))
     & .or.((kstep.gt.1).and.(LOP.eq.2))
     & .or.(Eval_Nuc_logical).or.(Eval_MSC_logical)) then
                   calculate_FIP = .true.
      else
                   calculate_FIP = .false.
      end if
	 
c Criteria to record max/min plastic strain based on loading steps
        if ((kstep.eq.eval_life_array(1)).and.pre_step) then
                   isMax_E_p_in = .true.
                   isMin_E_p_in = .false.
      else if ((kstep.eq.(eval_life_array(1)-1)).and.pre_step) then
                   isMin_E_p_in = .true.
                   isMax_E_p_in = .false.
      else 
                   isMin_E_p_in = .false.
                   isMax_E_p_in = .false.
      end if
#ifdef WRITEALL
        if ((kstep.gt.1).and.(mod(kstep,n_steps_per_cycle)
     &    .eq.0).and.pre_step) then
                   isMax_E_p_in_num = .true.
                   isMin_E_p_in_num = .false.
      else if ((kstep.gt.1).and.(mod(kstep,n_steps_per_cycle)
     &    .eq.(n_steps_per_cycle-1)).and.pre_step) then
                   isMin_E_p_in_num = .true.
                   isMax_E_p_in_num = .false.
      else 
                   isMin_E_p_in_num = .false.
                   isMax_E_p_in_num = .false.
      end if
#endif

c=====================================================================
c DATA ANALYSIS TASKS - Calculate FIP values, evaluate life of bands
c=====================================================================

c BEGIN FIP CALCULATION

      if (calculate_FIP) then
      print*, 'Calculate FIP'	  
		  
c----------------------------------------------------------------
c Calculate the change of shear strain on each slip plane, for each element
c----------------------------------------------------------------

        do k=1,num_elem
        do j=1,num_slip_sys
		
c        if (k.eq.6) then
c        print*,	gamma_cum_element(k,j)
c        print*,	gamma_cum_element_max(k,j),gamma_cum_element_min(k,j)	
c        end if	

      if (gamma_cum_element(k,j).lt.gamma_cum_element_min(k,j)) then
         gamma_cum_element_min(k,j)=gamma_cum_element(k,j)
c        if (k.eq.6) print*, 'min updated'		 
      end if
	  
      if (gamma_cum_element(k,j).gt.gamma_cum_element_max(k,j)) then 
         gamma_cum_element_max(k,j)=gamma_cum_element(k,j)
c        if (k.eq.6) print*, 'max updated'		 
      end if 
	 
        end do
        end do			  
c----------------------------------------------------------------
c Calculate ratcheting shear strain on each slip plane, for each element
c----------------------------------------------------------------

      do k=1,num_elem
       do j=1,num_slip_sys
       gamma_ratch(k,j)=abs(gamma_cum_element(k,j)
     & - gamma_cum_init_even_step(k,j)) 
c       gamma_ratch(k,j) = 0.0 !test for FIP calculation  
       end do
      end do

c----------------------------------------------------------------
c Calculate FIP based on the change in a cycle on each slip plane, for each element
c if element is damaged, set FIP = 0.0 for visualization purposes
c----------------------------------------------------------------

       do k=1,num_elem    
       do j = 1,num_slip_sys
       if (dam_elem(k).ne.1.) then
       FIP_cycle(k,j)=
     & 0.5 *( ( gamma_cum_element_max(k,j)-gamma_cum_element_min(k,j) )-
     & ( gamma_ratch(k,j) )) *
     & (  1. + Cont_FS_trans*sigma_gl(k,j)/Sigma_y) 
	 
c       print*, 'Normal Stress Term:',Cont_FS_trans*sigma_gl(k,j)/Sigma_y	 
	 
       else
          FIP_cycle(k,j)= 0.0 ! if element is damaged
       end if
       end do
       end do	
	   
c----------------------------------------------------------------
c Mash all layers into a single FIP for each element --> FIP_total
c FIP_total is a global variable stored in COMBLK and written to 
c SDV for plotting purposes
c----------------------------------------------------------------

c       print*
c       print*, "FIP List at step:", kstep," inc:", kinc
c       print*, "LOP:", LOP 
	   
       do k=1,num_elem  
	   
       FIP_total(k,13) = 0.0	 
	   
       do j = 1,num_slip_sys
       FIP_total(k,j) = FIP_cycle(k,j)
       FIP_total(k,13)= FIP_total(k,13) + FIP_cycle(k,j)   
       end do
c       print*, "FIP total of element", k, "=",  FIP_total(k,13) 	   
       end do	


	 
c END FIP CALCULATION
       end if !calculate_FIP

c-----------------------------------------------------
c Record the appropriate data for the relative position in load cycle
c-----------------------------------------------------	  
      if (isMin_E_p_in) then
        print*, "Record minimum of plastic strains for all elements "
        print*
        FNAMEW1=dmkname('Ep_el_min.csv',OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')E_p_elem(k,1,1),
     &  E_p_elem(k,1,2),
     &  E_p_elem(k,1,3),
     &  E_p_elem(k,2,1),
     &  E_p_elem(k,2,2),
     &  E_p_elem(k,2,3),
     &  E_p_elem(k,3,1),
     &  E_p_elem(k,3,2),
     &  E_p_elem(k,3,3)
        end do
        close(50)
        FNAMEW1=dmkname('S_el_min.csv',OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')sig_elem(k,1,1),
     &  sig_elem(k,1,2),
     &  sig_elem(k,1,3),
     &  sig_elem(k,2,1),
     &  sig_elem(k,2,2),
     &  sig_elem(k,2,3),
     &  sig_elem(k,3,1),
     &  sig_elem(k,3,2),
     &  sig_elem(k,3,3)
        end do
        close(50)
      end if
	  
      if (isMax_E_p_in) then
        print*, "Record maximum of plastic strains for all elements "
        print*
        FNAMEW1=dmkname('Ep_el_max.csv',OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')E_p_elem(k,1,1),
     &  E_p_elem(k,1,2),
     &  E_p_elem(k,1,3),
     &  E_p_elem(k,2,1),
     &  E_p_elem(k,2,2),
     &  E_p_elem(k,2,3),
     &  E_p_elem(k,3,1),
     &  E_p_elem(k,3,2),
     &  E_p_elem(k,3,3)
        end do
        close(50)
        FNAMEW1=dmkname('S_el_max.csv',OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')sig_elem(k,1,1),
     &  sig_elem(k,1,2),
     &  sig_elem(k,1,3),
     &  sig_elem(k,2,1),
     &  sig_elem(k,2,2),
     &  sig_elem(k,2,3),
     &  sig_elem(k,3,1),
     &  sig_elem(k,3,2),
     &  sig_elem(k,3,3)
        end do
        close(50)
      end if
#ifdef WRITEALL	 
      if (isMin_E_p_in_num) then
        write(filenum,'(I2)'), (kstep/n_steps_per_cycle+1)
        print*, "Record minimum of plastic strains for all elements "
        print*
        filenum = adjustl(filenum)
        filenames= 'Ep_el_min' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')E_p_elem(k,1,1),
     &  E_p_elem(k,1,2),
     &  E_p_elem(k,1,3),
     &  E_p_elem(k,2,1),
     &  E_p_elem(k,2,2),
     &  E_p_elem(k,2,3),
     &  E_p_elem(k,3,1),
     &  E_p_elem(k,3,2),
     &  E_p_elem(k,3,3)
        end do
        close(50)
        filenames= 'S_el_min' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')sig_elem(k,1,1),
     &  sig_elem(k,1,2),
     &  sig_elem(k,1,3),
     &  sig_elem(k,2,1),
     &  sig_elem(k,2,2),
     &  sig_elem(k,2,3),
     &  sig_elem(k,3,1),
     &  sig_elem(k,3,2),
     &  sig_elem(k,3,3)
        end do
        close(50)
        filenames= 'E_el_min' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')E_tot_elem(k,1,1),
     &  E_tot_elem(k,1,2),
     &  E_tot_elem(k,1,3),
     &  E_tot_elem(k,2,1),
     &  E_tot_elem(k,2,2),
     &  E_tot_elem(k,2,3),
     &  E_tot_elem(k,3,1),
     &  E_tot_elem(k,3,2),
     &  E_tot_elem(k,3,3)
        end do
        close(50)
        filenames= 'SS_gam_el_min' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 12(",", F))')gamma_cum_element(k,1),
     &  gamma_cum_element(k,2),
     &  gamma_cum_element(k,3),
     &  gamma_cum_element(k,4),
     &  gamma_cum_element(k,5),
     &  gamma_cum_element(k,6),
     &  gamma_cum_element(k,7),
     &  gamma_cum_element(k,8),
     &  gamma_cum_element(k,9),
     &  gamma_cum_element(k,10),
     &  gamma_cum_element(k,11),
     &  gamma_cum_element(k,12)
        end do
        close(50)
      end if
	  
      if (isMax_E_p_in_num) then
        print*, "Record maximum of plastic strains for all elements "
        print*
        write(filenum,'(I2)'), (kstep/n_steps_per_cycle)
        filenum = adjustl(filenum)
        filenames= 'Ep_el_max' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')E_p_elem(k,1,1),
     &  E_p_elem(k,1,2),
     &  E_p_elem(k,1,3),
     &  E_p_elem(k,2,1),
     &  E_p_elem(k,2,2),
     &  E_p_elem(k,2,3),
     &  E_p_elem(k,3,1),
     &  E_p_elem(k,3,2),
     &  E_p_elem(k,3,3)
        end do
        close(50)
        filenames= 'S_el_max' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 8(",", F))')sig_elem(k,1,1),
     &  sig_elem(k,1,2),
     &  sig_elem(k,1,3),
     &  sig_elem(k,2,1),
     &  sig_elem(k,2,2),
     &  sig_elem(k,2,3),
     &  sig_elem(k,3,1),
     &  sig_elem(k,3,2),
     &  sig_elem(k,3,3)
        end do
        close(50)
        filenames= 'SS_gam_el_max' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 12(",", F))')gamma_cum_element(k,1),
     &  gamma_cum_element(k,2),
     &  gamma_cum_element(k,3),
     &  gamma_cum_element(k,4),
     &  gamma_cum_element(k,5),
     &  gamma_cum_element(k,6),
     &  gamma_cum_element(k,7),
     &  gamma_cum_element(k,8),
     &  gamma_cum_element(k,9),
     &  gamma_cum_element(k,10),
     &  gamma_cum_element(k,11),
     &  gamma_cum_element(k,12)
        end do
        close(50)
        filenames= 'FIP_components' // trim(filenum) // '.csv'
        FNAMEW1=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do k=1,num_elem
        write(50,'(1x, F, 60(",", F))')gamma_cum_init_even_step(k,1),
     &  gamma_cum_init_even_step(k,2),
     &  gamma_cum_init_even_step(k,3),
     &  gamma_cum_init_even_step(k,4),
     &  gamma_cum_init_even_step(k,5),
     &  gamma_cum_init_even_step(k,6),
     &  gamma_cum_init_even_step(k,7),
     &  gamma_cum_init_even_step(k,8),
     &  gamma_cum_init_even_step(k,9),
     &  gamma_cum_init_even_step(k,10),
     &  gamma_cum_init_even_step(k,11),
     &  gamma_cum_init_even_step(k,12),
     &  gamma_cum_element(k,1),
     &  gamma_cum_element(k,2),
     &  gamma_cum_element(k,3),
     &  gamma_cum_element(k,4),
     &  gamma_cum_element(k,5),
     &  gamma_cum_element(k,6),
     &  gamma_cum_element(k,7),
     &  gamma_cum_element(k,8),
     &  gamma_cum_element(k,9),
     &  gamma_cum_element(k,10),
     &  gamma_cum_element(k,11),
     &  gamma_cum_element(k,12),
     &  gamma_cum_element_max(k,1),
     &  gamma_cum_element_max(k,2),
     &  gamma_cum_element_max(k,3),
     &  gamma_cum_element_max(k,4),
     &  gamma_cum_element_max(k,5),
     &  gamma_cum_element_max(k,6),
     &  gamma_cum_element_max(k,7),
     &  gamma_cum_element_max(k,8),
     &  gamma_cum_element_max(k,9),
     &  gamma_cum_element_max(k,10),
     &  gamma_cum_element_max(k,11),
     &  gamma_cum_element_max(k,12),
     &  gamma_cum_element_min(k,1),
     &  gamma_cum_element_min(k,2),
     &  gamma_cum_element_min(k,3),
     &  gamma_cum_element_min(k,4),
     &  gamma_cum_element_min(k,5),
     &  gamma_cum_element_min(k,6),
     &  gamma_cum_element_min(k,7),
     &  gamma_cum_element_min(k,8),
     &  gamma_cum_element_min(k,9),
     &  gamma_cum_element_min(k,10),
     &  gamma_cum_element_min(k,11),
     &  gamma_cum_element_min(k,12),
     &  sigma_gl(k,1),
     &  sigma_gl(k,2),
     &  sigma_gl(k,3),
     &  sigma_gl(k,4),
     &  sigma_gl(k,5),
     &  sigma_gl(k,6),
     &  sigma_gl(k,7),
     &  sigma_gl(k,8),
     &  sigma_gl(k,9),
     &  sigma_gl(k,10),
     &  sigma_gl(k,11),
     &  sigma_gl(k,12)
        end do
        close(50)
      end if
#endif

c-----------------------------------------------------
c Reset the extreme values of shear strain after every cycle
c NOTE done before 1st increment of step following completion of cycle
c but after FIP and NUC or MSC lives have been calculated - hence why its
c at the end of the code
c-----------------------------------------------------

      if (reset_shear_strain) then
	  
        print*, "Extreme Shear strains reset on kstep", kstep
        print*	  
		
         do k=1,num_elem
         do j=1,num_slip_sys
		 
          gamma_cum_element_min(k,j)=gamma_cum_element(k,j)
          gamma_cum_element_max(k,j)=gamma_cum_element(k,j)
          
         end do
         end do
		 
      end if !reset_shear_strain
      if ((kstep.gt.1).and.(mod(kstep,n_steps_per_cycle)
     &    .eq.0).and.pre_step) then
c      if (reset_shear_strain) then
        print*, "Reset Ratcheting Shear strains on kstep", kstep
        print*	  
		
         do k=1,num_elem
         do j=1,num_slip_sys
          gamma_cum_init_even_step(k,j)=gamma_cum_element(k,j)
         end do
         end do
      end if ! needs to reset at maximum so that ratcheting works properly
	
c=============================================================	  
c DATA OUTPUT TASKS (DAT) - write to files, print to JOBNAME.log 
c=============================================================

c--------------------------------------------------------		
c Write files for Nucleation
c--------------------------------------------------------

        if (Eval_Nuc_logical) then
		
         print*, "WRITE NUCLEATION FILES"

c Update number of failed grains:
        Nfailed = 1	  
c Write FIP value for each element following nucleation		
        FNAMEW1=dmkname('FIP_Nuc_el.txt',OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW1,status="unknown")
        do i=1,num_elem
           do j=1,num_slip_sys
        write(50,'(i9,2x,i6,2x,i6,2x,F,2x,F)')i,j,Elem_pos(i,5),FIP_cycle(i,j), 
     &  Cnum(Elem_pos(i,5),Elem_pos(i,ceiling(j/3.)),ceiling(j/3.)) 
        end do
        end do
        close(50)
		

c Write file of neighbour grain influence:		
        FNAMEW3=dmkname('d_gr_nd.txt',OUTDIR(1:LENOUTDIR),' ')
        OPEN(77,FILE=FNAMEW3,STATUS='unknown')
         do j=1,12
         do i=1,max_num_layers
         do k=1,num_grains
        write(77,'(i4,i4,i4,F)'), k,i,j,d_gr_nd(k,i,j)
          end do
          end do
          end do
        close(77) 

c-------------------------------------------------------
c Call python program:	
c-------------------------------------------------------	  
  
        PYMAIN=dmkname('Main.py',OUTDIR(1:LENOUTDIR),' ')
c        print*,PYMAIN	
c        print*, ('python2.7 -E ' // PYMAIN)
#ifdef CALLPYTHON	
        print*,'Call Python'		

c          CALL SYSTEM ('python2.7 -E ' // PYMAIN // ' -SI' ) !moleculos
c           CALL SYSTEM ('python -E ' // PYMAIN // ' -SI' ) !granulous (make sure .pbs scrpit changes to python 2.7 before execution!)
           CALL SYSTEM ('python -E ' // PYMAIN // ' -eN' ) !granulous (make sure .pbs scrpit changes to python 2.7 before execution!)
c          CALL SYSTEM ('python2.7 -E ' // PYMAIN)

         print*,'Return from Python'
#endif
c Check for crack arrest:
		 
         print*, 'Check for Crack Arrest'

        ARREST_FILE=dmkname('CRACK_ARRESTED.txt',OUTDIR(1:LENOUTDIR),' ')		 
		 
         INQUIRE(FILE = ARREST_FILE, EXIST=EXT)	
         print*, 'File Exists?: ', EXT
#ifdef DAMAGE		 
      if (EXT) then ! crack has arrested
          crack_arrested=.true.
          print*
          print*, 'CRACK HAS ARRESTED'
          print*		  
      end if		 
#endif		 
      call READCRACKEDELEMENTS(list_cracked_elem)



        end if !write nucleation files
		
c--------------------------------------------------------		
c Write files for MSC life of bands
c--------------------------------------------------------	

      if (Eval_MSC_logical) then 
	  
         print*, "WRITE MSC FILES"

c Update number of failed grains:
       Nfailed = Nfailed + 1	  
c Write FIP values for each element to a file	   
        write(filenum,'(I2)'), Nfailed
        filenum = adjustl(filenum)
        filenames= 'FIP_MSC' // trim(filenum) // '_el.txt'
         FNAMEW2=dmkname(filenames,OUTDIR(1:LENOUTDIR),' ')
        OPEN(50,FILE=FNAMEW2,status="unknown")
        do i=1,num_elem
           do j=1,12
         write(50,'(i9,2x,i6,2x,i6,2x,F)') i,j, Elem_pos(i,5),FIP_cycle(i,j)
        end do
        end do
        close(50)	   

         print*, "MSC FILES WRITTEN"


c-------------------------------------------------------
c Call python program:	
c-------------------------------------------------------	  

        PYMAIN=dmkname('Main.py',OUTDIR(1:LENOUTDIR),' ')
#ifdef CALLPYTHON
         print*,'Call Python'			   
c          CALL SYSTEM ('python2.7 -E ' // PYMAIN // ' -SI' ) !Moleculos SI
c          CALL SYSTEM ('python -E ' // PYMAIN // ' -SI' ) !Granulous SI
          CALL SYSTEM ('python -E ' // PYMAIN // ' -eN') !Granulous SII
         print*,'Return from Python'
#endif
c Check for crack arrest:
		 
         print*, 'Check for Crack Arrest'

        ARREST_FILE=dmkname('CRACK_ARRESTED.txt',OUTDIR(1:LENOUTDIR),' ')		 
		 
         INQUIRE(FILE = ARREST_FILE, EXIST=EXT)	
         print*, 'File Exists?: ', EXT
		 
      if (EXT) then ! crack has arrested
          crack_arrested=.true.
          print*
          print*, 'CRACK HAS ARRESTED'
          print*		  
      end if		 
	  
      call READCRACKEDELEMENTS(list_cracked_elem)

      end if !MSC life evaluation call
	  
c--------------------------------------------------------		
c Print summary to the JOBNAME.log file at end of job
c--------------------------------------------------------		  
	  
      if (LOP.eq.3) then ! This loop actives at the end of the job
        print*
        print*,'Final Call of UEXTERNALDB'		
           
      print*
      call system_clock(t_end)
      print*, 'UEXBD_time (s)        : ', UEXTERNALDB_time/1000000.
      print*, 'Elapsed_time (s)      : ', (t_end-START_time)/1000000.
      print*

      end if	  



c--------------------------------------------------------		
c Print stress/total-strain/plastic-strain data to a file
c for post-processing uses 
c--------------------------------------------------------
       if (LOP.eq.2) then !comment out this loop if you want to see where the solution diverges on CSSC
       do i=1,3	
       do j=1,3		 
         sig_elem_avg(i,j) = 0.0
         E_tot_elem_avg(i,j) = 0.0
         E_p_elem_avg(i,j) = 0.0
       end do
       end do

      do i = 1,3
      do j = 1,3
      do k = 1,num_elem			 
      sig_elem_avg(i,j) = sig_elem_avg(i,j) + (sig_elem(k,i,j)/num_elem)
      E_tot_elem_avg(i,j) = E_tot_elem_avg(i,j) + E_tot_elem(k,i,j)/num_elem
      E_p_elem_avg(i,j) = E_p_elem_avg(i,j) + E_p_elem(k,i,j)/num_elem 	  
      end do
      end do 
      end do 	  

c      print*, sig_elem(1,1,1)
c      print*, sig_elem_avg(1,1),E_tot_elem_avg(1,1),E_p_elem_avg(1,1)
	  
      CSSC_FILE=dmkname('CSSC_Data_22.txt',OUTDIR(1:LENOUTDIR),' ')
      OPEN(69,FILE=CSSC_FILE, ACCESS = 'APPEND' , status="unknown")
         write(69,'(F,F,F)')
     &   E_tot_elem_avg(2,2),E_p_elem_avg(2,2),sig_elem_avg(2,2)
      close(69)
	  
      B_FILE=dmkname('B_Data.txt',OUTDIR(1:LENOUTDIR),' ')
      OPEN(42,FILE=B_FILE, ACCESS = 'APPEND' , status="unknown")
c      do k = 1,num_elem
         k = 65	  
         write(42,'(F,F,F,F,F,F,F,F,F,F,F,F)')
     & a_elem(k,1),a_elem(k,2),a_elem(k,3),a_elem(k,4),a_elem(k,5),a_elem(k,6),
     & a_elem(k,7),a_elem(k,8),a_elem(k,9),a_elem(k,10),a_elem(k,11),a_elem(k,12)	 
c      end do
      close(42)	  
	  
	  
      end if
      

  
c Time how long the execution of UEXBD takes (printed every call):
      call system_clock(t_end)

      i = (t_end-t_start)

c      print*, 'UEXDB Execution took: ',(i/1000000.),' seconds'

      UEXTERNALDB_time = UEXTERNALDB_time + i
c      print*, 'UEXBD_time (s)        : ', UEXTERNALDB_time/1000000.
c      print*, 'Elapsed_time (s)      : ', (t_end-START_time)/1000000.
 

       print*, "UXDB Fin"

     	  
	  
           RETURN
           END


c====================================================================
c====================================================================
c====================== S U B R O U T I N E S =======================
c====================================================================
c====================================================================

c  Read in cracked elements and return modified arrays
      SUBROUTINE READCRACKEDELEMENTS(list_cracked_elem)
	 
#include "Geom_Def.txt"
#include "Common_block_Alv02.txt"
	  
         dimension list_cracked_elem(num_elem_fail,4)
         real*8 list_cracked_elem
         logical EXT
         character*512 OUTDIR
         character*511 dmkname
         character*250 CRCK_ELEM
         CALL GETOUTDIR(OUTDIR, LENOUTDIR )
		 
c-------------------------------------------------------
c Read in list of cracked elements and store in arrays:	
c-------------------------------------------------------
         CRCK_ELEM=dmkname('cracked_elem.txt',OUTDIR(1:LENOUTDIR),' ')
		 
         INQUIRE(FILE = CRCK_ELEM, EXIST=EXT)	
         print*, 'Cracked Element File Exists?: ', EXT
         if (EXT) then ! cracked elements exist, add to damaged list
          print*, 'Read in list of cracked elements:'	

         do jj = 1,num_elem_fail
         do nn = 1,4
           list_cracked_elem(jj,nn) = 0.0
         end do
         end do		 
		
         max_col = 4		
         max_row = num_elem_fail
         line_counter = 0         

         print*, 'Num_elem_fail: ', num_elem_fail
        
        OPEN(71,FILE=CRCK_ELEM,STATUS='OLD',ACTION='READ')
        do n_row = 1,max_row
        print*, n_row
        line_counter=n_row        
       READ(71,*,END=89)(list_cracked_elem(n_row,n_col),n_col=1,max_col) !Reads the next line, but to jump to statement 89 when the last line has already been read. 
        end do		
89      CLOSE(71)

         !Store the normals of the cracked elements:
         print*
         print*, 'print list of cracked elements:'		 
         do jj = 1,(line_counter-1) !to account for new line at end of file



             print*, 'Elem: ', nint(list_cracked_elem(jj,1))		   
		     do mm = 1,3			 
          Crack_plane_norm(nint(list_cracked_elem(jj,1)),mm) = 
     &			 list_cracked_elem(jj,(mm+1))
          print*,'n:',mm,' ',Crack_plane_norm(
     &           nint(list_cracked_elem(jj,1)),mm)


             end do
         end do
         print*		 
		 
         print*, 'Set the damage in each element'
	
         print*, 'line counter:', line_counter

	
         !Set dam_elem of elements in crack to 1 (TRUE):
         do ii = 1,(line_counter-1) !to account for new line at end of file
           print*, 'ii:', ii


             print*, 'Damage Element: ' , nint(list_cracked_elem(ii,1))
             dam_elem(nint(list_cracked_elem(ii,1))) = 1.0

         end do			

         print*, 'Set the damage in each element complete'
		  
         end if
         
      end


c
c  Calculate a vector cross product.
c
c  c = a cross b
c
c--------------------------------------------------------------------

      subroutine cross_product(a1,a2,a3,b1,b2,b3,c1,c2,c3)
      
      implicit double precision (a-h,o-z)
      
      c1 = a2 * b3 - a3 * b2
      c2 = a3 * b1 - a1 * b3
      c3 = a1 * b2 - a2 * b1

      return
      end

c====================================================================
c====================================================================
c
c  Normalize the length of a vector to one.
c
c--------------------------------------------------------------------

      subroutine normalize_vector(x,y,z)

      implicit double precision (a-h,o-z)

      xlength = sqrt(x*x+y*y+z*z)
      x = x / xlength
      y = y / xlength
      z = z / xlength

      return
      end

c====================================================================
c====================================================================
c
c  Transpose an ( n x n ) tensor.
c
c--------------------------------------------------------------------

      subroutine transpose(n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      do i = 1,n
         do j = 1,n
            b(i,j) = a(j,i)
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine aa_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n), c(n,n)

      do i = 1,n
         do j = 1,n
            c(i,j) = 0
            do k = 1,n
               c(i,j) = c(i,j) + a(i,k) * b(k,j)
            end do
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 2nd rank tensors.
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bb(n,a,b,sum)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      sum = 0.0
      do i = 1,n
         do j = 1,n
            sum = sum + a(i,j) * b(i,j)
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 4th rank tensors.
c  Result is stored in c(i,j,k,l)
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n,n,n), b(n,n,n,n), c(n,n,n,n)

      do i = 1,n
       do j = 1,n
        do k = 1,n
         do l = 1,n
          c(i,j,k,l) = 0
          do m1 = 1,n
           do m2 = 1,n
            c(i,j,k,l) = c(i,j,k,l) + a(i,j,m1,m2) * b(m1,m2,k,l)
           end do !m2
          end do !m1
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 4th rank tensor and
c  a 2nd rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n,n,n), b(n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(i,j,k,l) * b(k,l)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 2nd rank tensor and
c  a 4th rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n,n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(k,l) * b(k,l,i,j)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Rotates any 3x3x3x3 tensor by a rotation matrix.
c
c  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
c
c--------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &		a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &		a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &		a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &		a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &		a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &		a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &		a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &		a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &		a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &		a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
      b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
      b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
      b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
      b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
      b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
      b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

      det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j) / det
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Solve simultaneous equations using LU decomposition (Crout's method)
c  Result is stored in b(i)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine simeq(n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n), index(n)

      call LU_Decomp(n,a,index)
      call LU_BackSub(n,a,index,b)

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a matrix using 
c  LU decomposition (Crout's method)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine inverse(n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n), c(n,n), index(n)

      do i = 1,n
         do j = 1,n
            c(i,j) = a(i,j)
         end do
      end do

      do i = 1,n
         do j = 1,n
            b(i,j) = 0.0
         end do
         b(i,i) = 1.0
      end do

      call LU_Decomp(n,c,index)
      do j = 1,n
         call LU_BackSub(n,c,index,b(1,j))
      end do

      return
      end

c====================================================================
c====================================================================
c
c  This sub performs an LU Decomposition (Crout's method) on the 
c  matrix "a". It uses partial pivoting for stability. The index()
c  vector is used for the partial pivoting.  The v() vector is 
c  a dummy work area.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_Decomp(n,a,index)

      implicit double precision (a-h,o-z)

      dimension a(n,n), index(n), v(n)

      tiny = 1.0e-20

c--------------------------------------------------------------------
c  Loop over the rows to get the implicit scaling info.
c--------------------------------------------------------------------

      do i = 1,n
         a_max = 0.0
         do j = 1,n
            a_max = max(a_max,abs(a(i,j)))
         end do !j
         v(i) = 1.0 / a_max
      end do !i

c--------------------------------------------------------------------
c  Begin big loop over all the columns.
c--------------------------------------------------------------------

      do j = 1,n

         do i = 1,j-1
            sum = a(i,j)
            do k = 1,i-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
         end do

         a_max = 0.0
         do i = j,n
            sum = a(i,j)
            do k = 1,j-1
               sum = sum - a(i,k) * a(k,j)
            end do
            a(i,j) = sum
            dummy = v(i) * abs(sum)
            if ( dummy .gt. a_max ) then
               imax = i
               a_max = dummy
            end if
         end do

c--------------------------------------------------------------------
c  Pivot rows if necessary.
c--------------------------------------------------------------------

         if ( j .ne. imax ) then
            do k = 1,n
               dummy = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dummy
            end do
            v(imax) = v(j)
         end if
         index(j) = imax

c--------------------------------------------------------------------
c  Divide by the pivot element.
c--------------------------------------------------------------------

         if ( a(j,j) .eq. 0.0 ) a(j,j) = tiny
         if ( j .ne. n ) then
            dummy = 1.0 / a(j,j)
            do i = j+1,n
               a(i,j) = a(i,j) * dummy
            end do
         end if

      end do !j

      return
      end

c====================================================================
c====================================================================
c
c  Solves a set of simultaneous equations by doing back substitution.
c  The answer in returned in the b() vector.  The a(,) matrix
c  must have already been "LU Decomposed" by the above subroutine.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_BackSub(n,a,index,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), index(n), b(n)

      ii = 0

c--------------------------------------------------------------------
c  Do the forward substitution.
c--------------------------------------------------------------------

      do i = 1,n
         m = index(i)
         sum = b(m)
         b(m) = b(i)
         if ( ii .ne. 0 ) then
            do j = ii,i-1
               sum = sum - a(i,j) * b(j)
            end do
         else if ( sum .ne. 0.0 ) then
            ii = i
         end if
         b(i) = sum
      end do

c--------------------------------------------------------------------
c  Do the back substitution.
c--------------------------------------------------------------------

      do i = n,1,-1
         sum = b(i)
         if ( i .lt. n ) then
            do j = i+1,n
               sum = sum - a(i,j) * b(j)
            end do
         end if
         b(i) = sum / a(i,i)
      end do

      return
      end
      
c====================================================================
c====================================================================
c
c  Restore a symmetric 4th rank tensor stored in Voigt notation 
c  back to its 4th rank form.
c
c--------------------------------------------------------------------

      subroutine Voigt_to_forth(b,a)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          a(i,j,k,l) = b(ia,ib)
          if (ia.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
          if (ib.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
         end do
        end do
       end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Store a SYMMETRIC 4th rank tensor in Voigt notation.
c
c--------------------------------------------------------------------

      subroutine forth_to_Voigt(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=9-k-l
          b(ia,ib) = a(i,j,k,l)
         end do
        end do
       end do
      end do

      return
      end



c====================================================================
c====================================================================
c
c  Perform x**y but while retaining the sign of x.
c
c--------------------------------------------------------------------

      function power(x,y)

      implicit double precision (a-h,o-z)

      if (x.eq.0.0) then
        if (y.gt.0.0) then
          power = 0.0
        else if (y .lt. 0.0) then
          power = 1.0d+300
        else
          power = 1.0
        end if
      else
         power = y * log10(abs(x))
         if (power .gt. 300.) then
           power = 1.d+300
         else
           power = 10.d0 ** power
         end if
         if (x .lt. 0.0) power = -power
      end if

      return
      end

c====================================================================
c====================================================================
c
c  Return the sign of a number.
c
c--------------------------------------------------------------------

      function sgn(a)

      implicit double precision (a-h,o-z)

      sgn = 1.0
      if (a .lt. 0.0) sgn = -1.0

      return
      end 

c====================================================================
c====================================================================
c
c  Calculate the determinant of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      function determinant(a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end

c===================================================================
c===================================================================
c
c  Print out an array.
c
c-------------------------------------------------------------------

      subroutine print_array(n,a)

      implicit double precision (a-h,o-z)
      dimension a(n,n)

      do i = 1,n
         write(6,'(10f12.3)')(a(i,j),j=1,n)
      end do
      print*,' '

      return
      end

c===================================================================
c===================================================================
c
c  Print out Euler angles in Kocks notation.
c
c-------------------------------------------------------------------

      subroutine kocks_angles(npt,m,time,array1)

      implicit double precision (a-h,o-z)

      dimension array1(3,3)

      pi = 4 * atan(1.0)

      if (abs(array1(3,3)) .gt. 0.99999) then
        psi   = atan2(array1(2,1),array1(1,1))
        theta = 0.0
        phi   = 0.0
      else
        psi   = atan2(array1(2,3),array1(1,3))
        theta = acos(array1(3,3))
        phi   = atan2(array1(3,2),-array1(3,1))
      end if

      psi   = 180 * psi   / pi
      theta = 180 * theta / pi
      phi   = 180 * phi   / pi


c      write(7,'(a3,2i5,4f10.3)')'BM2',ninc,m,time,psi,theta,phi

      return
      end

c=======================================================================
c=======================================================================
c
c  Evaluate function to be minimized to zero.  gamma_dot()'s are
c  input and several things are output.
c
c-----------------------------------------------------------------------


      subroutine eval_func( xs0,	dtime,		gamma_dot,	
     &			    xm0,	F1,		num_slip_sys,
     &			    C,		F_p_inv,	F_p_inv_0,
     &			    g0,		Hdir,		Hdyn,
     &			    a0,		Adir,		Adyn,
     &			    a0_1,		Adir_1,		Adyn_1,
     &			    a0_2,		Adir_2,		Adyn_2,
     &			    a0_3,		Adir_3,		Adyn_3,
     &			    a0_4,		Adir_4,		Adyn_4,	 
     &			    g_sat,	F_el,		flow_exp,
     &			    sig,	tau,		gamma_dot_zero,
     &			    g,		func,		xL_p,
     &			    xLatent,	sse,		a,
     &			    a_1,  a_2,   a_3,   a_4, m_OW2)

      implicit double precision (a-h,o-z)
      
      dimension
     &	xs0(3,num_slip_sys),	xm0(3,num_slip_sys),	F_p_inv_0(3,3),
     &	F1(3,3),		C(3,3,3,3),		F_p_inv(3,3),
     &	g0(num_slip_sys),	xL_p_inter(3,3),	F_el(3,3),
     &	E_el(3,3),		Spk2(3,3),		sig(3,3),
     &	tau(num_slip_sys),	g(num_slip_sys),	array1(3,3),
     &	func(num_slip_sys),	gamma_dot(num_slip_sys),array2(3,3),
     &	F_el_inv(3,3),		xL_p(3,3),
     &	xs(3,num_slip_sys),	xm(3,num_slip_sys),
     &	a0(num_slip_sys),	a(num_slip_sys),
     &	a0_1(num_slip_sys),	a_1(num_slip_sys),
     &	a0_2(num_slip_sys),	a_2(num_slip_sys),
     &	a0_3(num_slip_sys),	a_3(num_slip_sys),
     &	a0_4(num_slip_sys),	a_4(num_slip_sys)	 


c*** Note that xs0 and xm0 are in INTERMEDIATE configuration!!!

c--------------------------------------------------------------------
c  Calculate the plastic part of the
c  velocity gradient in the intermediate configuration.
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          xL_p_inter(i,j) = 0.0
          do k = 1,num_slip_sys
            xL_p_inter(i,j) = xL_p_inter(i,j) + 
     &			xs0(i,k) * xm0(j,k) * gamma_dot(k)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Begin calculation process of F_p_n+1 = exp(xL_p_inter*dt).F_p_n
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          array1(i,j) = xL_p_inter(i,j) * dtime
        end do
      end do

c--------------------------------------------------------------------
c  Calculate omega.
c--------------------------------------------------------------------

      sum = 0
      do i = 1,3
        do j = 1,3
          sum = sum + array1(i,j) * array1(i,j)
        end do
      end do
c      omega = sqrt(0.5 * sum)

#ifdef DEBUG
      print*
      print*, "before comp - Time = ", time(2)
#endif



c--------------------------------------------------------------------
c  Continue calculating intermediate stuff needed for F_p_n+1
c--------------------------------------------------------------------

      call aa_dot_bb(3,array1,array1,array2)

        do i = 1,3
         if (abs(sum) .gt. 1.E-12) then   ! if omega=0 then no need.

#ifdef DEBUG
      print*
      print*, "Before omega - Time=", time(2)
#endif


       omega = sqrt(0.5 * sum)

            do j = 1,3
               array1(i,j) = array1(i,j) * sin(omega) / omega +
     &            array2(i,j) * (1-cos(omega)) / omega**2
            end do
         end if

         array1(i,i) = 1 + array1(i,i)
       end do

#ifdef DEBUG
      print*
      print*, "after omega - Time=", time(2)
#endif


c--------------------------------------------------------------------
c   Finally multiply arrays to get F_p_inv at end of time step.
c--------------------------------------------------------------------

      call inverse_3x3(array1,array2)
      call aa_dot_bb(3,F_p_inv_0,array2,F_p_inv)

c--------------------------------------------------------------------
c  Multiply F() by F_p_inv() to get F_el()
c--------------------------------------------------------------------

      call aa_dot_bb(3,F1,F_p_inv,F_el)
      call inverse_3x3(F_el,F_el_inv)

c--------------------------------------------------------------------
c  Rotate director vectors from intermediate configuration to
c  the current configuration.
c--------------------------------------------------------------------

      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate elastic Green Strain
c--------------------------------------------------------------------

      call transpose(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
      do i = 1,3
        E_el(i,i) = E_el(i,i) - 1
        do j = 1,3
          E_el(i,j) = E_el(i,j) / 2
        end do
      end do 

c--------------------------------------------------------------------
c  Multiply the stiffness tensor by the Green strain to get
c  the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------

      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from PK2 stress to Cauchy stress
c--------------------------------------------------------------------

      det = determinant(F_el)
      call transpose(3,F_el,array2)
      call aa_dot_bb(3,F_el,Spk2,array1)
      call aa_dot_bb(3,array1,array2,sig)
      do i = 1,3
        do j = 1,3
          sig(i,j) = sig(i,j) / det
        end do
      end do

c--------------------------------------------------------------------
c  Calculate resolved shear stress for each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        tau(k) = 0.0
        do j = 1,3
          do i = 1,3
            tau(k) = tau(k) + xs(i,k) * xm(j,k) * sig(i,j)
          end do 
        end do
      end do
      
c--------------------------------------------------------------------
c  Calculate hardening law for ref shear stress
c--------------------------------------------------------------------

      sum = 0
      do i = 1,num_slip_sys
        sum = sum + abs(gamma_dot(i))
      end do
      do k = 1,num_slip_sys
        sum00 = xLatent * sum - (xLatent - 1) * abs(gamma_dot(k))
        g(k) = (g0(k) + Hdir*sum00*dtime) / (1 + Hdyn*sum*dtime)
      end do

c--------------------------------------------------------------------
c  Calculate hardening law for back stress
c--------------------------------------------------------------------

      b_1 = Adir_1/Adyn_1 ! Saturation level
      b_2 = Adir_2/Adyn_2 ! Saturation level
!      b_3 = Adir_3/Adyn_3 ! Saturation level
!      b_4 = Adir_4/Adyn_4 ! Saturation level	  

      do k = 1,num_slip_sys
        a_1(k) = (a0_1(k) + Adir_1 * dtime * gamma_dot(k))  /
     & (1 + (Adir_1/b_1)*((a0_1(k)/b_1)**m_OW2)*dtime*abs(gamma_dot(k)))
	 
        a_2(k) = (a0_2(k) + Adir_2 * dtime * gamma_dot(k))  /
     & (1 + (Adir_2/b_2)*((a0_2(k)/b_2)**m_OW2)*dtime*abs(gamma_dot(k)))
	 
!        a_3(k) = (a0_3(k) + Adir_3 * dtime * gamma_dot(k))  /
!     & (1 + (Adir_3/b_3)*((a0_3(k)/b_3)**m_OW2)*dtime*abs(gamma_dot(k)))
	 
!        a_4(k) = (a0_4(k) + Adir_4 * dtime * gamma_dot(k))  /
!     & (1 + (Adir_4/b_4)*((a0_4(k)/b_4)**m_OW2)*dtime*abs(gamma_dot(k)))
	  
c Bound the back stress values for stability in damaged elements:
        if (a_1(k).gt.b_1) then
c          print*, "backstress1>Sat"  
c          print*, b_1, a_1(k)        		
           a_1(k) = b_1
        end if		
		
        if (a_1(k).lt.(-b_1)) then
c          print*, "backstress1<Sat"
c          print*, -b_1, a_1(k)     
           a_1(k) = -b_1
        end if
		
        if (a_2(k).gt.b_2) then
c          print*, "backstress2>Sat"
c          print*, b_2, a_2(k)  		
           a_2(k) = b_2
        end if		
		
        if (a_2(k).lt.(-b_2)) then
c          print*, "backstress2<Sat"
c          print*, -b_2, a_2(k)  
           a_2(k) = -b_2
        end if		  				
		
		!Sum the 4 back stresses:
        a(k) = a_1(k) + a_2(k)! + a_3(k) + a_4(k)
	   
      end do

c--------------------------------------------------------------------
c  Calculate function values.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        func(k) = tau(k)- a(k) - g(k) * power( (gamma_dot(k)/
     &			gamma_dot_zero), (1./flow_exp) )
      end do

c--------------------------------------------------------------------
c  Calculate root mean square error that is to be minimized to zero.
c--------------------------------------------------------------------

      sse = 0.0
      gamma_dot_max = abs(gamma_dot(1))
      do k = 1,num_slip_sys
        sse = sse + abs(gamma_dot(k)) * func(k) ** 2
        gamma_dot_max = max(abs(gamma_dot(k)),gamma_dot_max)
      end do

      if (gamma_dot_max .gt. 0.0) sse = sqrt(sse / gamma_dot_max) / 
     &							num_slip_sys

      return
      end

c====================================================================
c====================================================================
c
c  Convert Euler angles to a rotation matrix.  Euler angles are
c  in Bunge notation.
c  1) phi_1   around z-axis
c  2) Phi around x'-axis
c  3) phi_2   around z'-axis
c
c The transpose of the direction cosigns as defined by Bunge pg17 is
c taken here to perform rotations correctly using McGintys original
c conventions.
c i.e. The convention here for a 2nd order tensor transformation
c          Crystal_Frame = g^T * Global_Frame * g
c
c--------------------------------------------------------------------

      subroutine euler_to_matrix(phi_1,phi,phi_2,a)

          include 'ABA_PARAM.INC' 

      dimension a(3,3)

      s1 = sin(phi_1)
      c1 = cos(phi_1)
      s2 = sin(phi)
      c2 = cos(phi)
      s3 = sin(phi_2)
      c3 = cos(phi_2)
            
      a(1,1) = c1*c3-s1*s3*c2
      a(2,1) = s1*c3+c1*s3*c2
      a(3,1) = s3*s2
      a(1,2) = -c1*s3-s1*c3*c2
      a(2,2) = -s1*s3+c1*c3*c2
      a(3,2) = c3*s2
      a(1,3) = s1*s2
      a(2,3) = -c1*s2
      a(3,3) = c2

      return
      end

	  
c====================================================================
c=================================================================== 
c     
c  Calculate misorientation between two grains
c     
c------------------------------------------------------------------- 
      subroutine Gr_misor(phi_1_old,theta_old,phi_2_old,phi_1_cur,theta_cur,phi_2_cur,disAngle_temp)

      include 'ABA_PARAM.INC'   

c       real*8 g1,g2,g3,dis,angle1,angle2,angle3,angle4,angle5,angle6,angle7,angle8,angle9,angle10,
c     &angle11,angle12,angle13,angle14,angle15,angle16,angle17,angle18,angle19,angle20,
c     &angle21,angle22,angle23,angle24
       dimension  dis(3,3),
     &  g1 (3,3), !dir cos for old angles 
     &  g2 (3,3), !dir cos for updates angles
     &  g3 (3,3) !dir cos other

      CALL euler_to_matrix(phi_1_old,theta_old,phi_2_old,g1)
      CALL euler_to_matrix(phi_1_cur,theta_cur,phi_2_cur,g2)
      CALL transpose(3,g1,g3)
      CALL aa_dot_bb(3,g2,g3,dis)
      
c      d11 = dis(1,1)
c      print*, dis(1,1)
c      d12 = dis(1,2)
cc      print*, d12
c      d13 = dis(1,3)
cc      print*, d13
c      d21 = dis(2,1)
cc      print*, d21
c      d22 = dis(2,2)
cc      print*, d22
c      d23 = dis(2,3)top

cc      print*, d23
c      d31 = dis(3,1)
c      d32 = dis(3,2)
c      d33 = dis(3,3)
       
c	print*,dis(1,1),dis(2,2),dis(3,3),dis(3,1),dis(1,2),dis(2,3)

        CALL calc_misorient(dis(1,1),dis(2,2),dis(3,3),angle1)
        if (angle1.lt.disAngle_temp) then
            disAngle_temp=angle1
        end if

        CALL calc_misorient(-dis(1,1),dis(2,2),-dis(3,3),angle2)
        if (angle2.lt.disAngle_temp) then
            disAngle_temp=angle2
        end if

        CALL calc_misorient(-dis(1,1),-dis(2,2),dis(3,3),angle3)
        if (angle3.lt.disAngle_temp) then
            disAngle_temp=angle3
        end if

        CALL calc_misorient(dis(1,1),-dis(2,2),-dis(3,3),angle4)
        if (angle4.lt.disAngle_temp) then
            disAngle_temp=angle4
        end if

        CALL calc_misorient(dis(2,1),dis(3,2),dis(1,3),angle5)
        if (angle5.lt.disAngle_temp) then
            disAngle_temp=angle5
        end if

        CALL calc_misorient(-dis(2,1),dis(3,2),-dis(1,3),angle6)
        if (angle6.lt.disAngle_temp) then
            disAngle_temp=angle6
        end if

        CALL calc_misorient(-dis(2,1),-dis(3,2),dis(1,3),angle7)
        if (angle7.lt.disAngle_temp) then
            disAngle_temp=angle7
        end if

        CALL calc_misorient(dis(2,1),-dis(3,2),-dis(1,3),angle8)
        if (angle8.lt.disAngle_temp) then
            disAngle_temp=angle8
        end if

        CALL calc_misorient(dis(3,1),dis(1,2),dis(2,3),angle9)
        if (angle9.lt.disAngle_temp) then
            disAngle_temp=angle9
        end if

        CALL calc_misorient(-dis(3,1),dis(1,2),-dis(2,3),angle10)
        if (angle10.lt.disAngle_temp) then
            disAngle_temp=angle10
        end if

        CALL calc_misorient(-dis(3,1),-dis(1,2),dis(2,3),angle11)
        if (angle11.lt.disAngle_temp) then
            disAngle_temp=angle11
        end if

        CALL calc_misorient(dis(3,1),-dis(1,2),-dis(2,3),angle12)
        if (angle12.lt.disAngle_temp) then
            disAngle_temp=angle12
        end if

        CALL calc_misorient(-dis(3,1),-dis(2,2),-dis(1,3),angle13)
        if (angle13.lt.disAngle_temp) then
            disAngle_temp=angle13
        end if

        CALL calc_misorient(dis(3,1),-dis(2,2),dis(1,3),angle14)
        if (angle14.lt.disAngle_temp) then
            disAngle_temp=angle14
        end if

        CALL calc_misorient(dis(3,1),dis(2,2),-dis(1,3),angle15)
        if (angle15.lt.disAngle_temp) then
            disAngle_temp=angle15
        end if

        CALL calc_misorient(-dis(3,1),dis(2,2),dis(1,3),angle16)
        if (angle16.lt.disAngle_temp) then
            disAngle_temp=angle16
        end if

        CALL calc_misorient(-dis(1,1),-dis(3,2),-dis(2,3),angle17)
        if (angle17.lt.disAngle_temp) then
            disAngle_temp=angle17
        end if

        CALL calc_misorient(dis(1,1),-dis(3,2),dis(2,3),angle18)
        if (angle18.lt.disAngle_temp) then
            disAngle_temp=angle18
        end if

        CALL calc_misorient(dis(1,1),dis(3,2),-dis(2,3),angle19)
        if (angle19.lt.disAngle_temp) then
            disAngle_temp=angle19
        end if

        CALL calc_misorient(-dis(1,1),dis(3,2),dis(2,3),angle20)
        if (angle20.lt.disAngle_temp) then
            disAngle_temp=angle20
        end if

        CALL calc_misorient(-dis(2,1),-dis(1,2),-dis(3,3),angle21)
        if (angle21.lt.disAngle_temp) then
            disAngle_temp=angle21
        end if

        CALL calc_misorient(dis(2,1),-dis(1,2),-dis(3,3),angle22)
        if (angle22.lt.disAngle_temp) then
            disAngle_temp=angle22
        end if

        CALL calc_misorient(dis(2,1),dis(1,2),-dis(3,3),angle23)
        if (angle23.lt.disAngle_temp) then
            disAngle_temp=angle23
        end if

        CALL calc_misorient(-dis(2,1),dis(1,2),dis(3,3),angle24)
        if (angle24.lt.disAngle_temp) then
            disAngle_temp=angle24
        end if

c      print*, 'disAngle_temp ' ,disAngle_temp


      return
      end	  

      subroutine calc_misorient(a,b,c,angle)
      include 'ABA_PARAM.INC'         
c      implicit double precision (a-h,o-z)

c       real*8

      
           angle = (a+b+c)/2-.5
           if (angle.ge.1) then
             angle = 1
           end if 
           if (angle.le.-1) then
            angle = -1
           end if
           angle = acos(angle)
c           print*,angle

      return
      end	  
	  
c====================================================================
c=================================================================== 
c     
c Function to concatenate path and file name
c     
c------------------------------------------------------------------- 
              
      character*(*) function dmkname(fname,dname,exten)
C
      character*(*) fname,dname,exten
C     fname  I   jobname
C     dname  I   directory
C     exten  I   extension
C     dmkname O directory/jobname.exten

      ltot = len(fname)
      lf = 0
      do k1 = ltot,2,-1
        if (lf.eq.0.and.fname(k1:k1).ne.' ')  lf = k1
      end do

      ltot = len(dname)
      ld = 0
      do k1 = ltot,2,-1
        if (ld.eq.0.and.dname(k1:k1).ne.' ')  ld = k1
      end do

      ltot = len(exten)
      le = 0
      do k1 = ltot,2,-1
        if (le.eq.0.and.exten(k1:k1).ne.' ')  le = k1
      end do

      if ((lf + ld + le) .le. len(dmkname)) then
        dmkname = dname(1:ld)//'/'//fname(1:lf)
        ltot = ld + lf + 1
        if ( le.gt.0) then
           dmkname = dmkname(1:ltot)//exten(1:le)
        end if
      end if

      return
      end
      
