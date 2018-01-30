subroutine real2char(aaa,ndecim,char)
character*10 char
if(aaa.gt.10000) then
	write(*,*)' aaa=',aaa,' is larger than max value: 10,000'
	stop
end if
if(ndecim.gt.4) then
	write(*,*)' ndecim=',ndecim,' is larger than max value: 4'
	stop
end if

if(aaa.ge.0) then

	if(ndecim.eq.0) then
		kkk=aaa
		write(char,'(i4)')kkk
		if(kkk.lt.1000) write(char,'(i3)')kkk
		if(kkk.lt.100) write(char,'(i2)')kkk
		if(kkk.lt.10) write(char,'(i1)')kkk
	end if
	if(ndecim.eq.1) then
		write(char,'(f6.1)')aaa
		if(aaa.lt.1000) write(char,'(f5.1)')aaa
		if(aaa.lt.100) write(char,'(f4.1)')aaa
		if(aaa.lt.10) write(char,'(f3.1)')aaa
	end if
	if(ndecim.eq.2) then
		write(char,'(f7.2)')aaa
		if(aaa.lt.1000) write(char,'(f6.2)')aaa
		if(aaa.lt.100) write(char,'(f5.2)')aaa
		if(aaa.lt.10) write(char,'(f4.2)')aaa
	end if
	if(ndecim.eq.3) then
		write(char,'(f8.3)')aaa
		if(aaa.lt.1000) write(char,'(f7.3)')aaa
		if(aaa.lt.100) write(char,'(f6.3)')aaa
		if(aaa.lt.10) write(char,'(f5.3)')aaa
	end if
	if(ndecim.eq.4) then
		write(char,'(f9.4)')aaa
		if(aaa.lt.1000) write(char,'(f8.4)')aaa
		if(aaa.lt.100) write(char,'(f7.4)')aaa
		if(aaa.lt.10) write(char,'(f6.4)')aaa
	end if
else
	if(ndecim.eq.0) then
		kkk=aaa
		write(char,'(i5)')kkk
		if(abs(kkk).lt.1000) write(char,'(i4)')kkk
		if(abs(kkk).lt.100) write(char,'(i3)')kkk
		if(abs(kkk).lt.10) write(char,'(i2)')kkk
	end if
	if(ndecim.eq.1) then
		write(char,'(f7.1)')aaa
		if(abs(aaa).lt.1000) write(char,'(f6.1)')aaa
		if(abs(aaa).lt.100) write(char,'(f5.1)')aaa
		if(abs(aaa).lt.10) write(char,'(f4.1)')aaa
	end if
	if(ndecim.eq.2) then
		write(char,'(f8.2)')aaa
		if(abs(aaa).lt.1000) write(char,'(f7.2)')aaa
		if(abs(aaa).lt.100) write(char,'(f6.2)')aaa
		if(abs(aaa).lt.10) write(char,'(f5.2)')aaa
	end if
	if(ndecim.eq.3) then
		write(char,'(f9.3)')aaa
		if(abs(aaa).lt.1000) write(char,'(f8.3)')aaa
		if(abs(aaa).lt.100) write(char,'(f7.3)')aaa
		if(abs(aaa).lt.10) write(char,'(f6.3)')aaa
	end if
	if(ndecim.eq.4) then
		write(char,'(f10.4)')aaa
		!write(*,'(f10.4,1x,a10)')aaa,char
		if(abs(aaa).lt.1000) write(char,'(f9.4)')aaa
		if(abs(aaa).lt.100) write(char,'(f8.4)')aaa
		if(abs(aaa).lt.10) write(char,'(f7.4)')aaa
	end if

end if


return
end