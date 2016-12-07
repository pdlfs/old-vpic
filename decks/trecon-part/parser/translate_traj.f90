program tranlate
implicit none
logical check,write_select
integer nsp,n,myid
integer (kind=4) itype,tindex,tindex0,tstart,ncpu,nc,i,dt,ipart,pid
integer (kind=4) record_part,iselect,iselect_e,iselect_i,nselect,itotal
real(kind=4) r1,r2,ekmax_i,ekmax_e,ek1,ethe,ethi,ek1_i,ek2_i,ek1_e,ek2_e,mime
real(kind=8) d1, d2
integer (kind=4) npartmax,nsteps
parameter (nsteps=20000,npartmax=1000000)
integer (kind=4)tag_part,numpar_i(npartmax),numpar_e(npartmax),tag_select(npartmax)
real(kind=4) buffer_part(13),buffer_select(13,nsteps)
character(40) fname, fname_select

type:: p0header
  integer(kind=4) :: itype
  real(kind=4) ::  r1
  real(kind=8) :: d1
  real(kind=4) :: r2
  real(kind=8) :: d2
end type p0header

type(p0header) :: p0

inquire(iolength=record_part)real(buffer_part,kind=4)
write_select = .true.
tstart = 0
ncpu = 128
pid = 21

tindex = 12000
tindex0 = 0
mime = 25.0
ethe = 0.01
ethi = 0.01/sqrt(mime)
ek1_e = 10.50
ek2_e = 10.55
ek1_i = 5.8
ek2_i = 6.00
iselect = 0
iselect_e = 0
iselect_i = 0

!open(unit=11,file='selected-e.dat',status='unknown',form='binary')
!open(unit=12,file='selected-i.dat',status='unknown',form='binary')

!do nc = 0, ncpu-1

  write(fname,"(A5,A16,I0)") "traj/","electron_tracer.",pid
  inquire(file=trim(fname),exist=check)
  ekmax_e = 0.0
  if (check) then
    open(unit=13,file=trim(fname),status='unknown',form='binary')
    !read(13) ipart,r1,d1,r2,d2
    do i=1,20
      read(13) buffer_part
      ek1 = sqrt(1.0+buffer_part(5)**2+buffer_part(6)**2+buffer_part(7)**2)-1.0
      print*, buffer_part(1),ek1
      !ek1 = ek1/ethe
      !if (ek1>ekmax_e) ekmax_e = ek1
      !if (ek1 > ek1_e .and. ek1 < ek2_e) then
      !  write(11) tag_part,ek1
      !  iselect_e = iselect_e + 1
      !  if (iselect_e>npartmax) iselect_e = npartmax
      !  numpar_e(iselect_e) = tag_part
      !endif
    enddo
    close(13)
    print*, fname, ekmax_e, iselect_e
  endif

if (write_select) then
  nselect = iselect_e
  itotal = 0
  do nc = 0, ncpu-1
    write(fname,"(A9,I0,A17,I0)") "tracer/T.",tindex0,"/electron_tracer.",nc
    write(fname_select,"(A16,I0,A17,I0)") "tracer_select/T.",tindex0,"/electron_tracer.",nc

    ipart = 0
    inquire(file=trim(fname),exist=check)
    if (check) then
      open(unit=10,file=trim(fname),status='unknown',form='binary')
      read(10) itype,r1,d1,r2,d2
      do i=1, itype
        read(10) tag_part,buffer_part
        do iselect=1,nselect
          if (tag_part==numpar_e(iselect)) then
            !print*, "got",tag_part, numpar_e(iselect), iselect
            ipart = ipart + 1
            itotal = itotal + 1
            tag_select(ipart) = tag_part
            buffer_select(:,ipart) = buffer_part(:)
          endif
        enddo
      enddo
      close(10)

      open(unit=11,file=trim(fname_select),status='unknown',form='binary')
      write(11) ipart,r1,d1,r2,d2
      if (ipart>0) then
        do i=1,ipart
          write(11) tag_select(i),buffer_select(:,i)
        enddo
      endif
      close(11)
    endif
  enddo
  print*, itotal

endif


end program tranlate
