	! me = methylene
	! be = berendsen
	! si = silicon
	! hy = hydrogen
	subroutine init_arrays(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, i, j

	real*8 rmass_hy
	real*8 rmass_carbon
	real*8 rmass_oxygen
	real*8 rmass_si
	real*8 rmass_me
	real*8 rmass_methyl
	
	real*8 sigma_alcohol_hy_hy
	real*8 epsilon_alcohol_hy_hy
	real*8 sigma_alcohol_oxygen_oxygen
	real*8 epsilon_alcohol_oxygen_oxygen
	real*8 sigma_opls_me_me
	real*8 epsilon_opls_me_me
	real*8 sigma_opls_oxygen_oxygen
	real*8 epsilon_opls_oxygen_oxygen
	real*8 sigma_trappe_me_me
	real*8 epsilon_trappe_me_me
	real*8 sigma_be_si_si
	real*8 epsilon_be_si_si
	real*8 sigma_be_oxygen_oxygen
	real*8 epsilon_be_oxygen_oxygen
	real*8 sigma_be_methyl_methyl
	real*8 epsilon_be_methyl_methyl
	
	real*8 rbond_alcohol_hy_oxygen
	real*8 rbond_alcohol_oxygen_me
	real*8 rbond_opls_me_me
	real*8 rbond_opls_oxygen_me
	real*8 rbond_opls_me_oxygen
	real*8 rbond_trappe_me_me
	real*8 rbond_be_si_oxygen
	real*8 rbond_be_oxygen_si
	real*8 rbond_be_si_methyl
	real*8 rbond_be_me_si
	
	real*8 rangle_alcohol_hy_oxygen_me
	real*8 rangle_opls_oxygen_me_me
	real*8 rangle_opls_me_oxygen_me
	real*8 rangle_opls_me_me_oxygen
	real*8 rangle_trappe_me_me_me
	real*8 rangle_be_tetrahedral
	real*8 rangle_be_si_oxygen_si
	
	real*8 rcoeff_alcohol_hy_oxygen_me
	real*8 rcoeff_opls_oxygen_me_me
	real*8 rcoeff_opls_me_oxygen_me
	real*8 rcoeff_opls_me_me_oxygen
	real*8 rcoeff_trappe_me_me_me
	real*8 rcoeff_be_si_oxygen_si
	real*8 rcoeff_be_oxygen_si_oxygen
	real*8 rcoeff_be_oxygen_si_methyl
	real*8 rcoeff_be_methyl_si_oxygen
	real*8 rcoeff_be_methyl_si_methyl
	
	real*8 dv_be
	real*8 dv_trappe_a
	real*8 dv_trappe_b
	real*8 dv_trappe_c
	real*8 dv_opls_meth_oxy_meth_meth_a
	real*8 dv_opls_meth_oxy_meth_meth_b
	real*8 dv_opls_meth_oxy_meth_meth_c
	real*8 dv_opls_meth_meth_oxy_meth_a
	real*8 dv_opls_meth_meth_oxy_meth_b
	real*8 dv_opls_meth_meth_oxy_meth_c
	real*8 dv_alcohol_a
	real*8 dv_alcohol_b
	real*8 dv_alcohol_c
	real*8 dv_alcohol_oxy_meth_meth_meth_a
	real*8 dv_alcohol_oxy_meth_meth_meth_b
	real*8 dv_alcohol_oxy_meth_meth_meth_c
	real*8 dv_opls_oxy_meth_meth_oxy_a
	real*8 dv_opls_oxy_meth_meth_oxy_b
	real*8 dv_opls_oxy_meth_meth_oxy_c
	
	real*8 pi180
	
	pi180 = 4.0D0 * atan(1.0D0) / 180.0D0

	! begin assignment of factorial
	rifact(1) = 1.0D0
	rifact(2) = 1.0D0
	rifact(3) = 2.0D0
	rifact(4) = 6.0D0
	rifact(5) = 24.0D0
	rifact(6) = 120.0D0
	rifact(7) = 720.0D0
	rifact(8) = 5040.0D0
	rifact(9) = 40320.0D0
	rifact(10) = 362880.0D0
	rifact(11) = 3628800.0D0
	rifact(12) = 39916800.0D0
	rifact(13) = 479001600.0D0
	rifact(14) = 6227020800.0D0
	rifact(15) = 87178291200.0D0
	rifact(16) = 1307674368000.0D0
	rifact(17) = 20922789888000.0D0
	rifact(18) = 355687428096000.0D0
	rifact(19) = 6402373705728000.0D0
	rifact(20) = 121645100408832000.0D0
	rifact(21) = 2432902008176640000.0D0
	rifact(22) = 51090942171709440000.0D0
	rifact(23) = 1124000727777607680000.0D0
	rifact(24) = 25852016738884976640000.0D0
	rifact(25) = 620448401733239439360000.0D0
	rifact(26) = 15511210043330985984000000.0D0
	rifact(27) = 403291461126605635584000000.0D0
	rifact(28) = 10888869450418352160768000000.0D0
	rifact(29) = 304888344611713860501504000000.0D0
	rifact(30) = 8841761993739701954543616000000.0D0
	rifact(31) = 265252859812191058636308480000000.0D0
	rifact(32) = 8222838654177922817725562880000000.0D0
	rifact(33) = 263130836933693530167218012160000000.0D0
	! end assignment of factorial (too precision for large rifact)
	
	! begin assignment of nonbonded intramolecular interaction pairs
	iself(1) = 16
	iself(2) = 16
	iself(3) = 16
	iself(4) = 16
	iself(5) = 16
	iself(6) = 16
	iself(7) = 16
	iself(8) = 16
	
	jself(1) = 21
	jself(2) = 22
	jself(3) = 23
	jself(4) = 24
	jself(5) = 26
	jself(6) = 27
	jself(7) = 28
	jself(8) = 29
	  
	iself(9)  = 17
	iself(10) = 17
	iself(11) = 17
	iself(12) = 17
	iself(13) = 17
	iself(14) = 17
	  
	jself(9)  = 22
	jself(10) = 23
	jself(11) = 24
	jself(12) = 27
	jself(13) = 28
	jself(14) = 29
	  
	iself(15) = 19
	iself(16) = 19
	iself(17) = 19
	iself(18) = 19
	iself(19) = 19
	iself(20) = 19
	  
	jself(15) = 22
	jself(16) = 23
	jself(17) = 24
	jself(18) = 27
	jself(19) = 28
	jself(20) = 29
	  
	iself(21) = 20
	iself(22) = 20
	iself(23) = 20
	  
	jself(21) = 27
	jself(22) = 28
	jself(23) = 29
	  
	iself(24) = 21
	iself(25) = 21
	iself(26) = 21
	iself(27) = 21
	  
	jself(24) = 26
	jself(25) = 27
	jself(26) = 28
	jself(27) = 29
	  
	iself(28) = 22
	iself(29) = 22
	iself(30) = 22
	iself(31) = 22
	iself(32) = 22
	  
	jself(28) = 25
	jself(29) = 26
	jself(30) = 27
	jself(31) = 28
	jself(32) = 29
	  
	iself(33) = 23
	iself(34) = 23
	iself(35) = 23
	iself(36) = 23
	iself(37) = 23
	
	jself(33) = 25
	jself(34) = 26
	jself(35) = 27
	jself(36) = 28
	jself(37) = 29
	
	iself(38) = 24
	iself(39) = 24
	iself(40) = 24
	iself(41) = 24
	iself(42) = 24
	
	jself(38) = 25
	jself(39) = 26
	jself(40) = 27
	jself(41) = 28
	jself(42) = 29
	! end assignment of nonbonded intramolecular interaction pairs
	
	! begin assignment of surfactant bonds
	do i = 1, 18
	
		ibond(i) = i
		jbond(i) = i + 1
	
	enddo
	
	ibond(19) = 18
	jbond(19) = 20
	
	ibond(20) = 20
	jbond(20) = 21
	
	ibond(21) = 21
	jbond(21) = 22
	
	ibond(22) = 21
	jbond(22) = 23
	
	ibond(23) = 21
	jbond(23) = 24
	
	ibond(24) = 18
	jbond(24) = 25
	
	ibond(25) = 25
	jbond(25) = 26
	
	ibond(26) = 26
	jbond(26) = 27
	
	ibond(27) = 26
	jbond(27) = 28
	
	ibond(28) = 26
	jbond(28) = 29
	! end assignment of surfactant bonds
	
	rmass_hy = 1.00794D0
	rmass_carbon = 12.011D0
	rmass_oxygen = 15.9994D0
	rmass_si = 28.0855D0
	rmass_me = rmass_carbon + 2.0D0 * rmass_hy
	rmass_methyl = rmass_carbon + 3.0D0 * rmass_hy
	
	! begin assignment of atomic mass
	ems(1) = rmass_hy
	ems(2) = rmass_oxygen
	ems(3) = rmass_me
	ems(4) = rmass_me
	ems(5) = rmass_oxygen
	ems(6) = rmass_me
	ems(7) = rmass_me
	ems(8) = rmass_oxygen
	ems(9) = rmass_me
	ems(10) = rmass_me
	ems(11) = rmass_oxygen
	ems(12) = rmass_me
	ems(13) = rmass_me
	ems(14) = rmass_oxygen
	ems(15) = rmass_me
	ems(16) = rmass_me
	ems(17) = rmass_me
	ems(18) = rmass_si
	ems(19) = rmass_methyl
	ems(20) = rmass_oxygen
	ems(21) = rmass_si
	ems(22) = rmass_methyl
	ems(23) = rmass_methyl
	ems(24) = rmass_methyl
	ems(25) = rmass_oxygen
	ems(26) = rmass_si
	ems(27) = rmass_methyl
	ems(28) = rmass_methyl
	ems(29) = rmass_methyl
	
	do i = 1, 29
	   
	   ems(i) = ems(i) / rmass_oxygen
	   
	enddo
	
	emw(1) = rmass_oxygen
	emw(2) = rmass_hy
	emw(3) = rmass_hy
	
	do i = 1, 3
	   
	   emw(i) = emw(i) / rmass_oxygen
	   
	enddo
	! end assignment of atomic mass
	
	! begin assignment of Coulomb and Lennard-Jones parameters
	sigma_alcohol_hy_hy = 2.6D0
	epsilon_alcohol_hy_hy = 33.5D0
	sigma_alcohol_oxygen_oxygen = 3.07D0
	epsilon_alcohol_oxygen_oxygen = 711.76D0
	sigma_opls_me_me = 3.983D0
	epsilon_opls_me_me = 478.12D0
	sigma_opls_oxygen_oxygen = 3.047D0
	epsilon_opls_oxygen_oxygen = 817.67D0
	sigma_trappe_me_me = 3.95D0
	epsilon_trappe_me_me = 382.47D0
	sigma_be_si_si = 3.385D0
	epsilon_be_si_si = 2448.0D0
	sigma_be_oxygen_oxygen = 2.955D0
	epsilon_be_oxygen_oxygen = 849.3D0
	sigma_be_methyl_methyl = 3.786D0
	epsilon_be_methyl_methyl = 753.2D0
	
	sgs(1) = sigma_alcohol_hy_hy
	eps(1) = epsilon_alcohol_hy_hy
	rqs(1) = 0.435D0
	
	sgs(2) = sigma_alcohol_oxygen_oxygen
	eps(2) = epsilon_alcohol_oxygen_oxygen
	rqs(2) = -0.700D0
	
	sgs(3) = sigma_opls_me_me
	eps(3) = epsilon_opls_me_me
	rqs(3) = 0.265D0
	
	sgs(4) = sigma_opls_me_me
	eps(4) = epsilon_opls_me_me
	rqs(4) = 0.29D0
	
	sgs(5) = sigma_opls_oxygen_oxygen
	eps(5) = epsilon_opls_oxygen_oxygen
	rqs(5) = -0.58D0
	
	sgs(6) = sigma_opls_me_me
	eps(6) = epsilon_opls_me_me
	rqs(6) = 0.29D0
	
	sgs(7) = sigma_opls_me_me
	eps(7) = epsilon_opls_me_me
	rqs(7) = 0.29D0
	
	sgs(8) = sigma_opls_oxygen_oxygen
	eps(8) = epsilon_opls_oxygen_oxygen
	rqs(8) = -0.58D0
	
	sgs(9) = sigma_opls_me_me
	eps(9) = epsilon_opls_me_me
	rqs(9) = 0.29D0
	
	sgs(10) = sigma_opls_me_me
	eps(10) = epsilon_opls_me_me
	rqs(10) = 0.29D0
	
	sgs(11) = sigma_opls_oxygen_oxygen
	eps(11) = epsilon_opls_oxygen_oxygen
	rqs(11) = -0.58D0
	
	sgs(12) = sigma_opls_me_me
	eps(12) = epsilon_opls_me_me
	rqs(12) = 0.29D0
	
	sgs(13) = sigma_opls_me_me
	eps(13) = epsilon_opls_me_me
	rqs(13) = 0.29D0
	
	sgs(14) = sigma_opls_oxygen_oxygen
	eps(14) = epsilon_opls_oxygen_oxygen
	rqs(14) = -0.58D0
	
	sgs(15) = sigma_opls_me_me
	eps(15) = epsilon_opls_me_me
	rqs(15) = 0.29D0
	
	sgs(16) = sigma_trappe_me_me
	eps(16) = epsilon_trappe_me_me
	rqs(16) = 0.0D0
	
	sgs(17) = sigma_trappe_me_me
	eps(17) = epsilon_trappe_me_me
	rqs(17) = 0.0D0
	
	sgs(18) = sigma_be_si_si
	eps(18) = epsilon_be_si_si
	rqs(18) = 0.3D0
	
	sgs(19) = sigma_be_methyl_methyl
	eps(19) = epsilon_be_methyl_methyl
	rqs(19) = 0.0D0
	
	sgs(20) = sigma_be_oxygen_oxygen
	eps(20) = epsilon_be_oxygen_oxygen
	rqs(20) = -0.3D0
	
	sgs(21) = sigma_be_si_si
	eps(21) = epsilon_be_si_si
	rqs(21) = 0.15D0
	
	sgs(22) = sigma_be_methyl_methyl
	eps(22) = epsilon_be_methyl_methyl
	rqs(22) = 0.0D0
	
	sgs(23) = sigma_be_methyl_methyl
	eps(23) = epsilon_be_methyl_methyl
	rqs(23) = 0.0D0
	
	sgs(24) = sigma_be_methyl_methyl
	eps(24) = epsilon_be_methyl_methyl
	rqs(24) = 0.0D0
	
	sgs(25) = sigma_be_oxygen_oxygen
	eps(25) = epsilon_be_oxygen_oxygen
	rqs(25) = -0.3D0
	
	sgs(26) = sigma_be_si_si
	eps(26) = epsilon_be_si_si
	rqs(26) = 0.15D0
	
	sgs(27) = sigma_be_methyl_methyl
	eps(27) = epsilon_be_methyl_methyl
	rqs(27) = 0.0D0
	
	sgs(28) = sigma_be_methyl_methyl
	eps(28) = epsilon_be_methyl_methyl
	rqs(28) = 0.0D0
	
	sgs(29) = sigma_be_methyl_methyl
	eps(29) = epsilon_be_methyl_methyl
	rqs(29) = 0.0D0
	
	do i = 1, 29
	
		sgs(i) = sgs(i) / sigma
		eps(i) = eps(i) / epsilon

	enddo

	sgw(1) = sigma
	epw(1) = epsilon
	rqw(1) = -0.8476D0
	
	sgw(2) = 0.0D0
	epw(2) = 0.0D0
	rqw(2) = 0.4238D0
	
	sgw(3) = 0.0D0
	epw(3) = 0.0D0
	rqw(3) = 0.4238D0
	
	do i = 1, 3
	
		sgw(i) = sgw(i) / sigma
		epw(i) = epw(i) / epsilon

	enddo
	! end assignment of Coulomb and Lennard-Jones parameters

	! begin assignment of pre-computed sigma and epsilon values
	do i = 1, 3
		do j = 1, 3

			sgww6(i, j) = ((sgw(i) + sgw(j)) / 2.0D0)**6
			epww(i, j) = sqrt(epw(i) * epw(j))

		enddo
	enddo

	do i = 1, 29
		do j = 1, 29

			sgss6(i, j) = ((sgs(i) + sgs(j)) / 2.0D0)**6
			epss(i, j) = sqrt(eps(i) * eps(j))

		enddo
	enddo

	do i = 1, 3
		do j = 1, 29

			sgws6(i, j) = ((sgw(i) + sgs(j)) / 2.0D0)**6
			epws(i, j) = sqrt(epw(i) * eps(j))

		enddo
	enddo

	do i = 1, 29

		sgcs6(i) = ((3.19D0 / sigma + sgs(i)) / 2.0D0)**6
		epcs(i) = sqrt((392.0D0 / epsilon) * eps(i))

	enddo
	! end assignment of pre-computed sigma and epsilon values
	
	! begin assignment of squared bond lengths
	rbond_alcohol_hy_oxygen = 0.945D0
	rbond_alcohol_oxygen_me = 1.430D0
	rbond_opls_me_me = 1.516D0
	rbond_opls_oxygen_me = 1.410D0
	rbond_opls_me_oxygen = 1.410D0
	rbond_trappe_me_me = 1.540D0
	rbond_be_si_oxygen = 1.60D0
	rbond_be_oxygen_si = 1.60D0
	rbond_be_si_methyl = 1.88D0
	rbond_be_me_si = 1.88D0
	
	dsqs(1) = rbond_alcohol_hy_oxygen
	dsqs(2) = rbond_alcohol_oxygen_me
	dsqs(3) = rbond_opls_me_me
	dsqs(4) = rbond_opls_me_oxygen
	dsqs(5) = rbond_opls_oxygen_me
	dsqs(6) = rbond_opls_me_me
	dsqs(7) = rbond_opls_me_oxygen
	dsqs(8) = rbond_opls_oxygen_me
	dsqs(9) = rbond_opls_me_me
	dsqs(10) = rbond_opls_me_oxygen
	dsqs(11) = rbond_opls_oxygen_me
	dsqs(12) = rbond_opls_me_me
	dsqs(13) = rbond_opls_me_oxygen
	dsqs(14) = rbond_opls_oxygen_me
	dsqs(15) = rbond_opls_me_me
	dsqs(16) = rbond_trappe_me_me
	dsqs(17) = rbond_be_me_si
	dsqs(18) = rbond_be_si_methyl
	dsqs(19) = rbond_be_si_oxygen
	dsqs(20) = rbond_be_oxygen_si
	dsqs(21) = rbond_be_si_methyl
	dsqs(22) = rbond_be_si_methyl
	dsqs(23) = rbond_be_si_methyl
	dsqs(24) = rbond_be_si_oxygen
	dsqs(25) = rbond_be_oxygen_si
	dsqs(26) = rbond_be_si_methyl
	dsqs(27) = rbond_be_si_methyl
	dsqs(28) = rbond_be_si_methyl
	
	do i = 1, 28
	   
	   dsqs(i) = dsqs(i)**2 / sigma**2
	   
	enddo
	
	dsqw(1) = 1.0D0
	dsqw(2) = 2.0D0 * sin(109.47D0 * pi180 / 2.0D0)
	dsqw(3) = 1.0D0
	
	do i = 1, 3
	   
	   dsqw(i) = dsqw(i)**2 / sigma**2
	 
	enddo
	! end assignment of squared bond lengths
	
	! begin assignment of valence angles and coefficients
	rangle_alcohol_hy_oxygen_me = 108.5D0
	rangle_opls_oxygen_me_me = 108.0D0
	rangle_opls_me_oxygen_me = 112.0D0
	rangle_opls_me_me_oxygen = 108.0D0
	rangle_trappe_me_me_me = 114.0D0
	rangle_be_tetrahedral = 109.5D0
	rangle_be_si_oxygen_si = 144.0D0
	
	th(1) = rangle_alcohol_hy_oxygen_me
	th(2) = rangle_opls_oxygen_me_me
	th(3) = rangle_opls_me_me_oxygen
	th(4) = rangle_opls_me_oxygen_me
	th(5) = rangle_opls_oxygen_me_me
	th(6) = rangle_opls_me_me_oxygen
	th(7) = rangle_opls_me_oxygen_me
	th(8) = rangle_opls_oxygen_me_me
	th(9) = rangle_opls_me_me_oxygen
	th(10) = rangle_opls_me_oxygen_me
	th(11) = rangle_opls_oxygen_me_me
	th(12) = rangle_opls_me_me_oxygen
	th(13) = rangle_opls_me_oxygen_me
	th(14) = rangle_opls_oxygen_me_me
	th(15) = rangle_trappe_me_me_me
	th(16) = rangle_be_tetrahedral
	th(17) = rangle_be_tetrahedral
	th(18) = rangle_be_tetrahedral
	th(19) = rangle_be_tetrahedral
	th(20) = rangle_be_tetrahedral
	th(21) = rangle_be_tetrahedral
	th(22) = rangle_be_si_oxygen_si
	th(23) = rangle_be_si_oxygen_si
	th(24) = rangle_be_tetrahedral
	th(25) = rangle_be_tetrahedral
	th(26) = rangle_be_tetrahedral
	th(27) = rangle_be_tetrahedral
	th(28) = rangle_be_tetrahedral
	th(29) = rangle_be_tetrahedral
	th(30) = rangle_be_tetrahedral
	th(31) = rangle_be_tetrahedral
	th(32) = rangle_be_tetrahedral
	th(33) = rangle_be_tetrahedral
	th(34) = rangle_be_tetrahedral
	th(35) = rangle_be_tetrahedral
	th(36) = rangle_be_tetrahedral
	
	do i = 1, 36
	   
	   th(i) = th(i) * pi180
	   
	enddo
	
	rcoeff_alcohol_hy_oxygen_me = 519700.0D0
	rcoeff_opls_oxygen_me_me = 519700.0D0
	rcoeff_opls_me_oxygen_me = 519700.0D0
	rcoeff_opls_me_me_oxygen = 519700.0D0
	rcoeff_trappe_me_me_me = 519700.0D0
	rcoeff_be_si_oxygen_si = 118400.0D0
	rcoeff_be_oxygen_si_oxygen = 791200.0D0
	rcoeff_be_oxygen_si_methyl = 418400.0D0
	rcoeff_be_methyl_si_oxygen = 418400.0D0
	rcoeff_be_methyl_si_methyl = 418400.0D0
	
	rkth(1) = rcoeff_alcohol_hy_oxygen_me
	rkth(2) = rcoeff_opls_oxygen_me_me
	rkth(3) = rcoeff_opls_me_me_oxygen
	rkth(4) = rcoeff_opls_me_oxygen_me
	rkth(5) = rcoeff_opls_oxygen_me_me
	rkth(6) = rcoeff_opls_me_me_oxygen
	rkth(7) = rcoeff_opls_me_oxygen_me
	rkth(8) = rcoeff_opls_oxygen_me_me
	rkth(9) = rcoeff_opls_me_me_oxygen
	rkth(10) = rcoeff_opls_me_oxygen_me
	rkth(11) = rcoeff_opls_oxygen_me_me
	rkth(12) = rcoeff_opls_me_me_oxygen
	rkth(13) = rcoeff_opls_me_oxygen_me
	rkth(14) = rcoeff_opls_oxygen_me_me
	rkth(15) = rcoeff_trappe_me_me_me
	rkth(16) = rcoeff_trappe_me_me_me
	rkth(17) = rcoeff_be_methyl_si_methyl
	rkth(18) = rcoeff_be_methyl_si_oxygen
	rkth(19) = rcoeff_be_methyl_si_oxygen
	rkth(20) = rcoeff_be_methyl_si_oxygen
	rkth(21) = rcoeff_be_methyl_si_oxygen
	rkth(22) = rcoeff_be_si_oxygen_si
	rkth(23) = rcoeff_be_si_oxygen_si
	rkth(24) = rcoeff_be_oxygen_si_methyl
	rkth(25) = rcoeff_be_oxygen_si_methyl
	rkth(26) = rcoeff_be_oxygen_si_methyl
	rkth(27) = rcoeff_be_methyl_si_methyl
	rkth(28) = rcoeff_be_methyl_si_methyl
	rkth(29) = rcoeff_be_methyl_si_methyl
	rkth(30) = rcoeff_be_oxygen_si_methyl
	rkth(31) = rcoeff_be_oxygen_si_methyl
	rkth(32) = rcoeff_be_oxygen_si_methyl
	rkth(33) = rcoeff_be_methyl_si_methyl
	rkth(34) = rcoeff_be_methyl_si_methyl
	rkth(35) = rcoeff_be_methyl_si_methyl
	rkth(36) = rcoeff_be_oxygen_si_oxygen
	
	do i = 1, 36
	   
	   rkth(i) = rkth(i) / epsilon
	
	enddo
	! end assignment of valence angles and coefficients
	
	! begin assignment of dihedral coefficients
	dv_be = 3770.0D0
	dv_trappe_a = 5903.8D0
	dv_trappe_b = -1134.0D0
	dv_trappe_c = 13159.0D0
	dv_opls_meth_oxy_meth_meth_a = 19860.0D0
	dv_opls_meth_oxy_meth_meth_b = -5853.0D0
	dv_opls_meth_oxy_meth_meth_c = 8956.0D0
	dv_opls_meth_meth_oxy_meth_a = 19860.0D0
	dv_opls_meth_meth_oxy_meth_b = -5853.0D0
	dv_opls_meth_meth_oxy_meth_c = 8956.0D0
	dv_alcohol_a = 3490.0D0
	dv_alcohol_b = -486.0D0
	dv_alcohol_c = 3130.0D0
	dv_alcohol_oxy_meth_meth_meth_a = 2940.0D0
	dv_alcohol_oxy_meth_meth_meth_b = -888.0D0
	dv_alcohol_oxy_meth_meth_meth_c = 12810.0D0
	dv_opls_oxy_meth_meth_oxy_a = 2940.0D0
	dv_opls_oxy_meth_meth_oxy_b = -888.0D0
	dv_opls_oxy_meth_meth_oxy_c = 12810.0D0
	
	dva(1) = dv_alcohol_a
	dvb(1) = dv_alcohol_b
	dvc(1) = dv_alcohol_c
	
	dva(2) = dv_opls_oxy_meth_meth_oxy_a
	dvb(2) = dv_opls_oxy_meth_meth_oxy_b
	dvc(2) = dv_opls_oxy_meth_meth_oxy_c
	
	dva(3) = dv_opls_meth_meth_oxy_meth_a
	dvb(3) = dv_opls_meth_meth_oxy_meth_b
	dvc(3) = dv_opls_meth_meth_oxy_meth_c
	
	dva(4) = dv_opls_meth_oxy_meth_meth_a
	dvb(4) = dv_opls_meth_oxy_meth_meth_b
	dvc(4) = dv_opls_meth_oxy_meth_meth_c
	
	dva(5) = dv_opls_oxy_meth_meth_oxy_a
	dvb(5) = dv_opls_oxy_meth_meth_oxy_b
	dvc(5) = dv_opls_oxy_meth_meth_oxy_c
	
	dva(6) = dv_opls_meth_meth_oxy_meth_a
	dvb(6) = dv_opls_meth_meth_oxy_meth_b
	dvc(6) = dv_opls_meth_meth_oxy_meth_c
	
	dva(7) = dv_opls_meth_oxy_meth_meth_a
	dvb(7) = dv_opls_meth_oxy_meth_meth_b
	dvc(7) = dv_opls_meth_oxy_meth_meth_c
	
	dva(8) = dv_opls_oxy_meth_meth_oxy_a
	dvb(8) = dv_opls_oxy_meth_meth_oxy_b
	dvc(8) = dv_opls_oxy_meth_meth_oxy_c
	
	dva(9) = dv_opls_meth_meth_oxy_meth_a
	dvb(9) = dv_opls_meth_meth_oxy_meth_b
	dvc(9) = dv_opls_meth_meth_oxy_meth_c
	
	dva(10) = dv_opls_meth_oxy_meth_meth_a
	dvb(10) = dv_opls_meth_oxy_meth_meth_b
	dvc(10) = dv_opls_meth_oxy_meth_meth_c
	
	dva(11) = dv_opls_oxy_meth_meth_oxy_a
	dvb(11) = dv_opls_oxy_meth_meth_oxy_b
	dvc(11) = dv_opls_oxy_meth_meth_oxy_c
	
	dva(12) = dv_opls_meth_meth_oxy_meth_a
	dvb(12) = dv_opls_meth_meth_oxy_meth_b
	dvc(12) = dv_opls_meth_meth_oxy_meth_c
	
	dva(13) = dv_opls_meth_oxy_meth_meth_a
	dvb(13) = dv_opls_meth_oxy_meth_meth_b
	dvc(13) = dv_opls_meth_oxy_meth_meth_c
	
	dva(14) = dv_alcohol_oxy_meth_meth_meth_a
	dvb(14) = dv_alcohol_oxy_meth_meth_meth_b
	dvc(14) = dv_alcohol_oxy_meth_meth_meth_c
	
	dva(15) = dv_trappe_a
	dvb(15) = dv_trappe_b
	dvc(15) = dv_trappe_c
	
	dva(16) = dv_trappe_a
	dvb(16) = dv_trappe_b
	dvc(16) = dv_trappe_c
	
	dva(17) = dv_trappe_a
	dvb(17) = dv_trappe_b
	dvc(17) = dv_trappe_c
	
	dva(18) = dv_trappe_a
	dvb(18) = dv_trappe_b
	dvc(18) = dv_trappe_c
	
	dva(19) = dv_be
	dvb(19) = 0.0D0
	dvc(19) = 0.0D0
	
	dva(20) = dv_be
	dvb(20) = 0.0D0
	dvc(20) = 0.0D0
	
	dva(21) = dv_be
	dvb(21) = 0.0D0
	dvc(21) = 0.0D0
	
	dva(22) = dv_be
	dvb(22) = 0.0D0
	dvc(22) = 0.0D0
	
	dva(23) = dv_be
	dvb(23) = 0.0D0
	dvc(23) = 0.0D0
	
	dva(24) = dv_be
	dvb(24) = 0.0D0
	dvc(24) = 0.0D0
	
	dva(25) = dv_be
	dvb(25) = 0.0D0
	dvc(25) = 0.0D0
	
	dva(26) = dv_be
	dvb(26) = 0.0D0
	dvc(26) = 0.0D0
	
	dva(27) = dv_be
	dvb(27) = 0.0D0
	dvc(27) = 0.0D0
	
	dva(28) = dv_be
	dvb(28) = 0.0D0
	dvc(28) = 0.0D0
	
	dva(29) = dv_be
	dvb(29) = 0.0D0
	dvc(29) = 0.0D0
	
	dva(30) = dv_be
	dvb(30) = 0.0D0
	dvc(30) = 0.0D0
	
	! the following is due to the assumed form of eval_dihedral
	do i = 1, 18
	   
	   dva(i) = dva(i) / 2.0D0
	   dvb(i) = dvb(i) / 2.0D0
	   dvc(i) = dvc(i) / 2.0D0
	   
	enddo
	
	do i = 1, 30
	   
	   dva(i) = dva(i) / epsilon
	   dvb(i) = dvb(i) / epsilon
	   dvc(i) = dvc(i) / epsilon
	   
	enddo
	! end assignment of dihedral coefficients
	
	! begin initialization of forces
	fxw1 = 0.0D0
	fyw1 = 0.0D0
	fzw1 = 0.0D0
	
	fxw2 = 0.0D0
	fyw2 = 0.0D0
	fzw2 = 0.0D0
	
	fxs = 0.0D0
	fys = 0.0D0
	fzs = 0.0D0
	! end initialization of forces
	
	! begin writing out input parameters
	open(10, file="output_input.dat")
	write(10,*) "SURFACTANT ATOMIC MASS (gm/mole)"
	do i = 1, 29
	   write(10,*) i, (ems(i) * rmass_oxygen)
	enddo
	write(10,*) "WATER ATOMIC MASS (gm/mole)"
	do i = 1, 3
	   write(10,*) i, (emw(i) * rmass_oxygen)
	enddo
	write(10,*) "SURFACTANT LENNARD-JONES PARAMETERS AND CHARGE (Angstrom, J/mole, e)"
	do i = 1, 29
	   write(10,*) i, (sgs(i) * sigma), (eps(i) * epsilon), (rqs(i))
	enddo
	write(10,*) "WATER LENNARD-JONES PARAMETERS AND CHARGE (Angstrom, J/mole, e)"
	do i = 1, 3
	   write(10,*) i, (sgw(i) * sigma), (epw(i) * epsilon), (rqw(i))
	enddo
	write(10,*) "SURFACTANT BOND LENGTHS (Angstrom)"
	do i = 1, 28
	   write(10,*) i, dsqrt((dsqs(i) * sigma**2))
	enddo
	write(10,*) "WATER BOND LENGTHS (Angstrom)"
	do i = 1, 3
	   write(10,*) i, dsqrt((dsqw(i) * sigma**2))
	enddo
	write(10,*) "SURFACTANT BOND ANGLES (degree)"
	do i = 1, 36
	   write(10,*) i, (th(i) / pi180)
	enddo
	write(10,*) "SURFACTANT BOND ANGLE COEFFICIENTS (J/mole/rad/rad)"
	do i = 1, 36
	   write(10,*) i, (rkth(i) * epsilon)
	enddo
	write(10,*) "SURFACTANT TORSIONAL COEFFICIENTS (J/mole)"
	do i = 1, 18
	   write(10,*) i, (dva(i) * epsilon * 2.0D0), (dvb(i) * epsilon * 2.0D0), (dvc(i) * epsilon * 2.0D0)
	enddo
	do i = 19, 30
	   write(10,*) i, (dva(i) * epsilon), (dvb(i) * epsilon), (dvc(i) * epsilon)
	enddo
	close(10)
	! end writing out input parameters

	return
	end subroutine
