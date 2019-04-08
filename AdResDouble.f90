program AdRes
	
use randgen
implicit none

integer, parameter :: long = selected_real_kind(15, 50)
double precision, parameter :: PI=3.141592653589793238462
double precision, parameter :: Dp = 2e-17
double precision, parameter :: rI = 0.5e-6
double precision, parameter :: db = 0.01e-6
double precision, parameter :: dinf = 1e-6
double precision, parameter :: d0= 1e-6
integer, parameter :: R = 100
double precision, parameter :: dt = 0.8e-6
double precision, parameter :: b = 60e-9
double precision, parameter :: sigma = 1.2e-9
double precision, parameter :: kB = 1.4e-23
double precision, parameter :: T = 300
double precision, parameter :: vs = 0.001
integer, parameter :: s = 3
integer, parameter :: Z = 2

double precision, dimension(3) :: sphere, p, midb
double precision, dimension(2) :: seed
double precision, dimension(R+1, 10) :: q
double precision, dimension(s**2, 3) :: newc, newsp
double precision, allocatable, dimension(:, :) :: beads, springs, buffer, spbuffer
integer, dimension(R) :: resolutions = spread(0, 1, R)
double precision :: min_distance, avrad, drag, k, suml=0
logical :: sts

integer :: i, j, x, tot, jx, jx_1, minj, time, j1, jN, ebcount, m, a, rA, rB, w

! INITIALIZATION
! q holds information on boundary beads
! r holds information on all beads
! springs holds information on the springs and associated timescales + radii

sts = .TRUE.

call FJC([real(0,8), real(0, 8), real(0, 8)], R, s*b, q(:, 1:3))

q(:, 4:6) = q(:, 1:3)
q(:, 7:9) = q(:, 1:3)
q(:, 10) = spread(1, 1, R+1)

allocate(beads(R+1, 10))
allocate(springs(R, 3))
allocate(buffer(R+1, 10))
allocate(spbuffer(R, 3))

time = 0
beads = q
springs = spread([double precision :: s**4, b*s, sigma*s**2], 1, R)

midb = q(51, 1:3)
call random_number(seed)
seed(1) = 2*PI*seed(1)
seed(2) = acos(2*seed(2) - 1)
sphere = [d0*cos(seed(1))*sin(seed(2)), d0*sin(seed(1))*sin(seed(2)), d0*cos(seed(2))]
p = midb + sphere

do while (sts)
	time = time +1
	do x=1, size(springs, 1)-1
		jx = int(springs(x, 1))
		jx_1 = int(springs(x+1, 1))
		minj = min(jx, jx_1)
		if (mod(time, minj)==0) then
			avrad = 0.5*(springs(x,3) + springs(x+1, 3))
			drag = 6*vs*avrad*PI
			beads(x+1, 1:3) = beads(x+1, 4:6) + sqrt(2*kB*T*minj*dt/drag)*&
			[random_normal(), random_normal(), random_normal()]
			if (mod(time, jx)==0) then
				k = 3*kB*T/(springs(x,2)**2)
				if (jx == s**4) then
					beads(x+1, 1:3) = beads(x+1, 1:3) + (k/drag)*(beads(x,7:9) - beads(x+1, 7:9))*jx*dt
				else
					beads(x+1, 1:3) = beads(x+1, 1:3) + (k/drag)*(beads(x,4:6) - beads(x+1, 4:6))*jx*dt
				end if
			end if
			if (mod(time, jx_1)==0) then
				k = 3*kB*T/(springs(x+1,2)**2)
				if (jx_1 == s**4) then
					beads(x+1, 1:3) = beads(x+1, 1:3) + (k/drag)*(beads(x+2,7:9) - beads(x+1, 7:9))*jx_1*dt
				else
					beads(x+1, 1:3) = beads(x+1, 1:3) + (k/drag)*(beads(x+2,4:6) - beads(x+1, 4:6))*jx_1*dt
				end if
			end if
		else
			beads(x+1, 1:3) = beads(x+1, 4:6)
		end if
	end do
	
	j1 = int(springs(1,1))
	jN = int(springs(size(springs, 1), 1))
	
	if (mod(time, j1)==0) then
		avrad = 0.5*(sigma + springs(1, 3))
		drag = 6*vs*PI*avrad
		k = 3*kB*T/springs(1, 2)**2
		if (j1 == s**4) then
			beads(1, 1:3) = beads(1, 4:6) + (k/drag)*(beads(2, 7:9) - beads(1, 7:9))*j1*dt + &
			sqrt(2*kB*T*j1*dt/drag)*[random_normal(), random_normal(), random_normal()]	
		else
			beads(1, 1:3) = beads(1, 4:6) + (k/drag)*(beads(2, 4:6) - beads(1, 4:6))*j1*dt + &
			sqrt(2*kB*T*j1*dt/drag)*[random_normal(), random_normal(), random_normal()]	
		end if
	else
		beads(1, 1:3) = beads(1, 4:6)
	end if
	
	if (mod(time, jN)==0) then
		avrad = 0.5*(sigma + springs(size(springs,1), 3))
		drag = 6*vs*PI*avrad
		k = 3*kB*T/springs(size(springs,1), 2)**2
		if (jN == s**4) then
			beads(size(springs,1) + 1, 1:3) = beads(size(springs,1) + 1, 4:6) + (k/drag)*(beads(size(springs,1), 7:9)&
			- beads(size(springs,1)+1, 7:9))*jN*dt + sqrt(2*kB*T*jN*dt/drag)*[random_normal(), random_normal(), random_normal()]	
		else
			beads(size(springs,1) + 1, 1:3) = beads(size(springs,1) + 1, 4:6) + (k/drag)*(beads(size(springs,1), 4:6)&
			- beads(size(springs,1)+1, 4:6))*jN*dt + sqrt(2*kB*T*jN*dt/drag)*[random_normal(), random_normal(), random_normal()]	
		end if
	else
		beads(size(springs,1) + 1, 1:3) = beads(size(springs,1) + 1, 4:6)
	end if
	
	! Updating array before next timestep
	
	beads(:, 4:6) = beads(:, 1:3)
	if (mod(time,int(s**4))==0) then
		beads(:, 7:9) = beads(:, 1:3)
	end if

	! Updating end beads array
	
	ebcount=1
	do m=1, size(beads, 1)
		if (beads(m, 10)==1) then
			q(ebcount, :) = beads(m, :)
			ebcount = ebcount +1
		end if
	end do
	
	! Updating position of protein
	
	p = p + sqrt(2*Dp*dt)*random_normal()
	
	! Now checking for successful binding and resolution increase/decrease
	
	if (min_distance(beads, p, size(beads,1)) < db) then
		Print*, "Binding successful"
		sts = .FALSE.
	else if (min_distance(q, p, size(q, 1)) > dinf) then
		Print*, "Binding failed"
		sts = .FALSE.
	end if
	
	do a=1,R
		if(resolutions(a)==0) then
		! spring is currently in low resolution - check to see if we should zoom in
			if (norm2(q(a, 1:3) - p) < rI .OR. norm2(q(a+1,1:3) - p) < rI) then
				! Introduce new beads using MH algorithm 
				print*, "zooming in", a, time
				call findloc(beads(:, 1), q(a, 1), size(beads(:, 1), 1), rA)
				call findloc(beads(:, 1), q(a+1, 1), size(beads(:, 1), 1), rB)
				call mh_chain(q(a, 1:3), q(a+1, 1:3), s, newc)
				deallocate(buffer)
				allocate(buffer(size(beads, 1), 10))
				buffer = beads
				deallocate(beads)
				allocate(beads(size(buffer, 1) + 8, 10))
				beads(1:rA, :) = buffer(1:rA, :)
				beads(rA+1:rA+8, 1:3) = newc(2:, :)
				beads(rA+1:rA+8, 4:6) = newc(2:, :)
				beads(rA+1:rA+8, 7:9) = newc(2:, :)
				beads(rA+1:rA+8, 10) = spread(0, 1, 8) ! not end beads!
				beads(rA+9:, :) = buffer(rB:, :)
				newsp = spread([double precision :: 1.0, b, sigma], 1, 9)
				deallocate(spbuffer)
				allocate(spbuffer(size(springs, 1), 3))
				spbuffer = springs
				deallocate(springs)
				allocate(springs(size(spbuffer, 1) + 8, 3))
				springs(1:rA-1, :) = spbuffer(1:rA-1, :)
				springs(rA:rA+8, :) = newsp
				springs(rA+9:, :) = spbuffer(rA+1:, :)
				resolutions(a) = 1
			end if
			
		! spring is currently in high resolution - check to see if we should zoom out
		else if (resolutions(a)==1) then
			if (norm2(q(a, 1:3) - p) > Z*rI .AND. norm2(q(a+1, 1:3) - p) > Z*rI) then
				! zoom out and delete beads as necessary
				print*, "zooming out", a, time
				call findloc(beads(:, 1), q(a, 1), size(beads(:, 1), 1), rA)
				call findloc(beads(:, 1), q(a+1, 1), size(beads(:, 1), 1), rB)
				deallocate(buffer)
				allocate(buffer(size(beads, 1), 10))
				buffer = beads
				deallocate(beads)
				allocate(beads(size(buffer, 1) - 8, 10))
				beads(1:rA, :) = buffer(1:rA, :)
				beads(rA+1:, :) = buffer(rB:, :)
				deallocate(spbuffer)
				allocate(spbuffer(size(springs, 1), 3))
				spbuffer = springs
				deallocate(springs)
				allocate(springs(size(spbuffer, 1) - 8, 3))
				springs(1:rA-1, :) = spbuffer(1:rA-1, :)
				springs(rA, :) = [double precision :: s**4, b*s, sigma*s**2]
				springs(rA+1:, :) = spbuffer(rB:, :)
				resolutions(a) = 0
			end if
		end if
	
	end do
	!print*, min_distance(beads, p, size(beads,1))
end do
end program

function min_distance(beads, p, nbeads)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	integer, intent(in) :: nbeads
	double precision, dimension(nbeads, 3), intent(in) :: beads
	double precision, dimension(3), intent(in) :: p
	double precision, dimension(nbeads) :: distances
	double precision :: min_distance
	integer :: i
	
	do i=1, nbeads
		distances(i) = dot_product(p-beads(i,:), p-beads(i,:))
	end do
	
	min_distance = sqrt(minval(distances))
	
end function

subroutine FJC(start, n, r, chain)
	
	implicit none
	double precision, parameter :: PI=3.141592653589793238462
	integer, parameter :: long = selected_real_kind(15, 50)
	double precision, dimension(3), intent(in) :: start
	integer, intent(in) :: n
	double precision, intent(in) :: r
	double precision, dimension(n+1, 3), intent(out) :: chain
	double precision :: lambda, phi
	double precision, dimension(3) :: sphere
	double precision, dimension(2) :: seed
	integer :: i
	
	chain(1, :) = start
	
	do i=2, n+1
		call random_number(seed)
		lambda = 2*PI*seed(1)
		phi = acos(2*seed(2) - 1)
		sphere = [r*cos(lambda)*cos(phi), r*cos(phi)*sin(lambda), r*sin(phi)]
		chain(i, :) = chain(i-1, :) + sphere
	end do
	
end subroutine FJC

function ddelta(x)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)

	double precision, intent(in) :: x
	integer :: ddelta
	
    if(x == 0) then
        ddelta = 1
    else
        ddelta = 0
	end if
	
end function
	
function phi(dd, F, rB, length)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	double precision, intent(in) :: dd
	integer, intent(in) :: length
	double precision, dimension(3), intent(in) :: rB
	double precision, dimension(length, 3), intent(in) :: F
	double precision :: lambda, prod, phi
	integer :: i, ddelta
	
	lambda = norm2(rB - F(length, :))
	
	do i=2,length
		prod = prod*ddelta(norm2(F(i, :) - F(i-1,:)) - dd)
	end do
	if (lambda > 2*dd) then
		phi = 0
	else
		phi = sqrt(dd**2 - 0.25*lambda**2)*prod
	end if

end function

subroutine mh_chain(rA, rB, s, F)

! For the moment, I'm going to stick with this initial chain. I will add the M-H aspect in later

	implicit none
	double precision, parameter :: PI=3.141592653589793238462
	integer, parameter :: long = selected_real_kind(15, 50)
	double precision, dimension(3), intent(in) :: rA, rB
	integer, intent(in) :: s
	double precision, dimension(s**2, 3), intent(out) :: F
	
	double precision :: d, dd, l, theta
	double precision, dimension(3) :: diff, midpoint, axis1, axis2
	logical :: redraw=.TRUE.
	
	d = norm2(rA-rB)
	dd = d/s
	
	do while(redraw)
		call FJC(rA, s**2-2, dd, F(1:s**2-1, :))
		if (norm2(rB - F(s**2-1, :)) < 2*dd) then
			redraw = .FALSE.
			diff = rB - F(s**2-1, :)
			midpoint = F(s**2-1, :) + 0.5*diff
			axis1 = [1, 1, 1]
			if (diff(3) /= 0) then
				axis1(3) = -(diff(1) + diff(2))/diff(3)
			else if (diff(2) /= 0) then
				axis1(2) = -(diff(1) + diff(3))/diff(2)
			else if (diff(1) /= 0) then
				axis1(1) = -(diff(2) + diff(3))/diff(1)
			end if
			axis1 = axis1/norm2(axis1)
			
			axis2 = [1, 2, 3]
			if (diff(3) /= 0) then
				axis2(3) = -(diff(1) + 2*diff(2))/diff(3)
			else if (diff(2) /= 0) then
				axis2(2) = -(diff(1) + 3*diff(3))/diff(2)
			else if (diff(1) /= 0) then
				axis2(1) = -(2*diff(2) + 3*diff(3))/diff(1)
			end if
			axis2 = axis2/norm2(axis2)
			
			! Now use Gram Schmidt
			
			axis2 = axis2 - (dot_product(axis1, axis2)/norm2(axis1))*axis1
			axis2 = axis2/norm2(axis2)
			
			call random_number(theta)
			theta = theta*2*PI
			l = sqrt(dd**2 - 0.25*norm2(rB - F(s**2-1, :))**2)
			F(s**2, :) = midpoint + cos(theta)*l*axis1 + sin(theta)*l*axis2
		end if
	end do
	
end subroutine

subroutine findloc(array1, target_value, num_elements, loc)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	integer, intent(in) :: num_elements
	double precision, intent(in) :: target_value
	double precision, dimension(num_elements), intent(in) :: array1
	integer, intent(out) :: loc
	integer :: i
	
	do i = 1, num_elements
		if (array1(i) .eq. target_value) then
			loc = i
			exit
		endif
	end do

end subroutine
