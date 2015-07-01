	program mybox
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 ierr, mype, nprs
	integer*4 i, j, k, m, kb
	integer*4 ifar, iwrite, ibox64, ibox2
	integer*4 ib4, ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000
	integer*4 np, nc, df
	real*8 rlxtime, starttime
	
	! topology variables
	integer*4 ndims, idims(3), iposit(3)
	integer*4 my3D, comm3D, ieast, iwest, inorth, isouth, iup, idown
	logical isperiodic(3), reorder
	
	call mpi_init(ierr)
	call mpi_comm_rank(mpi_comm_world, mype, ierr)
	call mpi_comm_size(mpi_comm_world, nprs, ierr)
	call init_box_dimensions()

	ifar = 10
	iwrite = 100

	! setup 3D topology
	ndims = 3
	idims(1) = ibdz2
	idims(2) = ibdy2
	idims(3) = ibdx2
	isperiodic(1) = .false.
	isperiodic(2) = .false.
	isperiodic(3) = .false.
	reorder = .true.

	call mpi_cart_create(mpi_comm_world, ndims, idims, isperiodic, reorder, comm3D, ierr)
	call mpi_cart_get(comm3D, ndims, idims, isperiodic, iposit, ierr)
	call mpi_cart_rank(comm3D, iposit, my3D, ierr)
	call mpi_cart_shift(comm3D, 2, 1, iwest, ieast, ierr)
	call mpi_cart_shift(comm3D, 1, 1, isouth, inorth, ierr)
	call mpi_cart_shift(comm3D, 0, 1, idown, iup, ierr)
	
	np = 3 * iwater + 29 * isurfactant
	nc = 3 * iwater + 28 * isurfactant
	df = 3 * np - nc
	
	tr = 3.81087D0
	rlxtime = 7.5D0 * dt
	pnh1 = dt / df / tr / rlxtime**2
	pnh2 = -2.0D0 * df * tr
	eta = 0.0D0
	
	call get_box_center(comm3D, my3D, iposit)
	call init_arrays(my3D)
	
	call read_config(my3D)
	call get_molecule_count(my3D)
	call read_graphite(my3D, iposit(1))
	call init_temp(my3D)
	call compute_temperature(my3D)
	if(my3D.ne.0) call config_to_root(my3D)
	if(my3D.eq.0) call config_from_root(my3D, nprs, 0)
	
	call eastSR(comm3D, my3D, ieast, iwest)
	call westSR(comm3D, my3D, iwest, ieast)
	call northSR(comm3D, my3D, inorth, isouth)
	call southSR(comm3D, my3D, isouth, inorth)
	call up_downSR_isr_update(comm3D, my3D, iup, idown, iposit(1))
	
	call buildSnnn(my3D)
	call buildSnn(my3D)
	call buildSn(my3D)
	call buildW2nn(my3D)
	call buildW2n(my3D)
	
	call computeOnm4(my3D)
	call east_west_Onm4(comm3D, my3D, ieast, iwest)
	call north_south_Onm4(comm3D, my3D, inorth, isouth)
	call up_down_Onm4_isr(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown, iposit)
	call computeLjk4(my3D)
	
	call evalLjk4(my3D)
	
	call evalw1w1n(my3D)
	call evalw1w2n(my3D)
	call evalw1w2nn(my3D)
	call evalw1sn(my3D)
	call evalw1snn(my3D)
	call evalw1snnn(my3D)
	call evalw2w1n(my3D)
	call evalw2w2n(my3D)
	call evalw2w2nn(my3D)
	call evalw2sn(my3D)
	call evalw2snn(my3D)
	call evalw2snnn(my3D)
	call evalw2nw1n(my3D)
	call evalw2nw2n(my3D)
	call evalw2nw2nn(my3D)
	call evalw2nS(my3D)
	call evalw2nSn(my3D)
	call evalw2nSnn(my3D)
	call evalw2nSnnn(my3D)
	call evalSw1n(my3D)
	call evalSw2nn(my3D)
	call evalSSn(my3D)
	call evalSSnn(my3D)
	call evalSSnnn(my3D)
	call evalSnw1n(my3D)
	call evalSnw2nn(my3D)
	call evalSnSn(my3D)
	call evalSnSnn(my3D)
	call evalSnSnnn(my3D)
	call evalSnnw1n(my3D)
	call evalSnnw2nn(my3D)
	call evalSnnSnn(my3D)
	call evalSnnSnnn(my3D)
	
	call buildW2force(my3D)
	call buildSforce(my3D)
	call ewnsud_W2force(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
	call ewnsud_Sforce(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
	call assign_split_forces(my3D)
	
	!if(my3D.ne.0) call force_to_root(my3D)
	!if(my3D.eq.0) call force_from_root(my3D, nprs, 0)
	
	call init_accel(my3D)
	
	do kb = 99101, 150000
	
		if(my3D.eq.0) write(*,*) "*********"
		if(my3D.eq.0) write(*,*) "kb = ", kb
		starttime = mpi_wtime()
		
		call move_surfactant_a(my3D)
		call move_w1_a(my3D)
		call move_w2_a(comm3D, my3D)
		
		if(mod(kb, ifar).eq.0) then
		
			call delete_s(my3D)
			call ewnsud_deleted_s(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
			call add_s(my3D)
			
			call delete_w1(my3D)
			call ewnsud_deleted_w1(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
			call add_w1(my3D)
			
			call delete_w2(my3D)
			call ewnsud_deleted_w2(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
			call add_w2(my3D)
			
		endif
		
		call eastSR(comm3D, my3D, ieast, iwest)
		call westSR(comm3D, my3D, iwest, ieast)
		call northSR(comm3D, my3D, inorth, isouth)
		call southSR(comm3D, my3D, isouth, inorth)
		
		if(mod(kb, ifar).eq.0) then
		
			call up_downSR_isr_update(comm3D, my3D, iup, idown, iposit(1))
			
			call buildSnnn(my3D)
			call buildSnn(my3D)
			call buildSn(my3D)
			call buildW2nn(my3D)
			call buildW2n(my3D)
			
			call computeOnm4(my3D)
			call east_west_Onm4(comm3D, my3D, ieast, iwest)
			call north_south_Onm4(comm3D, my3D, inorth, isouth)
			call up_down_Onm4_isr(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown, iposit)
			call computeLjk4(my3D)
		
		else
		
			call up_downSR_isr(comm3D, my3D, iup, idown, iposit(1))
		
		endif
		
		call evalLjk4(my3D)
		
		call evalw1w1n(my3D)
		call evalw1w2n(my3D)
		call evalw1w2nn(my3D)
		call evalw1sn(my3D)
		call evalw1snn(my3D)
		call evalw1snnn(my3D)
		call evalw2w1n(my3D)
		call evalw2w2n(my3D)
		call evalw2w2nn(my3D)
		call evalw2sn(my3D)
		call evalw2snn(my3D)
		call evalw2snnn(my3D)
		call evalw2nw1n(my3D)
		call evalw2nw2n(my3D)
		call evalw2nw2nn(my3D)
		call evalw2nS(my3D)
		call evalw2nSn(my3D)
		call evalw2nSnn(my3D)
		call evalw2nSnnn(my3D)
		call evalSw1n(my3D)
		call evalSw2nn(my3D)
		call evalSSn(my3D)
		call evalSSnn(my3D)
		call evalSSnnn(my3D)
		call evalSnw1n(my3D)
		call evalSnw2nn(my3D)
		call evalSnSn(my3D)
		call evalSnSnn(my3D)
		call evalSnSnnn(my3D)
		call evalSnnw1n(my3D)
		call evalSnnw2nn(my3D)
		call evalSnnSnn(my3D)
		call evalSnnSnnn(my3D)
		
		call buildW2force(my3D)
		call buildSforce(my3D)
		call ewnsud_W2force(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
		call ewnsud_Sforce(comm3D, my3D, ieast, iwest, inorth, isouth, iup, idown)
		call assign_split_forces(my3D)
		
		call move_surfactant_b(my3D)
		call move_w1_b(my3D)
		call move_w2_b(my3D)
		
		if(mod(kb, iwrite).eq.0) then
		if(my3D.ne.0) call config_to_root(mype)
		if(my3D.eq.0) call config_from_root(mype, nprs, kb)
		
		!if(my3D.ne.0) call force_to_root(my3D)
		!if(my3D.eq.0) call force_from_root(my3D, nprs, kb)
		endif
		
		call get_molecule_count(mype)
		call compute_temperature(mype)
		
		if(mype.eq.4) write(*,*) "time (s) = ", mpi_wtime() - starttime
	enddo

	call mpi_finalize(ierr)

	end program
