module POTwrapper


!~**********************************************************************
!~* Purpose: temporary-before-major-cleanup wrapper to old phd subs for POT sampling
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 03/09/2009
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1.
!~*		2.
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none

Private

public:: POTindices_daily


contains

subroutine POTIndices_daily(xin,&
				YearStart,YearEnd, & !define season of interest (IN)
				T, &! threshold
				C1,C2,& !independence constraints
				POT,&!peaks
				IED,& !inter-event duration
				mvflag,& !any mv during flood event?
				YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				err,mess)

!^**********************************************************************
!^* Purpose: extract POT indices from daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:06/02/2009
!^**********************************************************************
!^* Comments: change yearlist from real to integer pointer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. xin, input time series;
!^*		2.YearStart
!^*		3.YearEnd
!^*		4.T, threshold
!^*		5.C1,C2 independence constraints
!^* OUT
!^*		1. POT
!^*		2. N
!^*		3. mvflag
!^*		4. YearList
!^*		5. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6. mess, error message
!^**********************************************************************
use TimeSeries_tools, only: OneDTimeSeriesType
use Dates_tools

type(OneDTimeSeriesType), intent(in)::xin
type(DateType), intent(in),optional::YearStart, YearEnd
real(mrk),intent(in)::T,C1,C2
real(mrk), pointer::POT(:),IED(:)
integer(mik),pointer::yearlist(:),mvflag(:)
integer(mik), intent(out)::err
 character(*),intent(out)::mess
! locals
integer(mik)::i,n
real(mrk),allocatable::jour(:),annee(:),debit(:)
real(mrk)::sdeb,sfin
real, pointer::resultat(:, :)
integer, pointer::indx(:)

if(allocated(jour)) deallocate(jour)
allocate(jour(xin%n))
if(allocated(annee)) deallocate(annee)
allocate(annee(xin%n))
if(allocated(debit)) deallocate(debit)
allocate(debit(xin%n))

do i=1,xin%n
	jour(i)=date2Days(xin%ts(i)%date - datetype(xin%ts(i)%date%year,1,1,0,0,0))+1
	annee(i)=xin%ts(i)%date%year
	debit(i)=xin%ts(i)%q
enddo
where(debit<0) debit=-9999.
sdeb=date2Days(YearStart .Yminus. datetype(1,1,1,0,0,0))+1
sfin=date2Days(YearEnd .Yminus. datetype(1,1,1,0,0,0))+1

call SUPSEUILC3S(real(jour), real(annee), real(debit), real(T), real(C1), real(C2), &
				resultat, real(sdeb), real(sfin))

if(.not. associated(resultat)) then
    err=1;
    mess='POTIndices_daily:Fatal:problem extracting threshold exceedances; threshold might never be exceeded'
    return
endif
!sort resultat
n=size(resultat, dim=1)
if(n==0) then
    err=1;
    mess='POTIndices_daily:Fatal:problem extracting threshold exceedances; threshold is never exceeded'
    return
endif

call indexs(resultat(:,6), indx)
! compute POT, IED, Yearlist and mvflag
if (associated(IED)) deallocate(IED)
allocate(IED(n))
IED(1)=resultat(indx(1),6)
do i=2,n
	IED(i)=(resultat(indx(i),6)-resultat(indx(i-1),6))*1._mrk
enddo

if (associated(POT)) deallocate(POT)
allocate(POT(n))
POT=resultat(indx,4)*1000._mrk

if (associated(mvflag)) deallocate(mvflag)
allocate(mvflag(n))
mvflag=int(resultat(indx,13))

if (associated(YearList)) deallocate(YearList)
allocate(YearList(n))
YearList=int(resultat(indx,2))

end subroutine POTIndices_daily

!!!!!!!!!
!PRIVATE!
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SUPSEUILC3S(jour, annee, debit, S, C1, C2, resultat, sdeb, sfin)
! extraction de variables sup-seuil pour une saison comprise entre sdeb et sfin
! Par convention, si sdeb>sfin, on s'intéresse à la saison 2 (commençant en sdeb)
! (cf Renard 2004 pour notations et definitions)
! la matrice résultat contient les colonnes (dans cet ordre):
! rang de l'evt/annee/jour/q/teta/temps cumulé corrigé/delta/m/r/sigma/gamma/
! indice dans Xt/code lacune/complexité (=nb de maximums locaux)
! sdeb est le premier jour de l'année hydro

integer, PARAMETER:: npt=int(1.E3)
real, intent(in)::jour(:), annee(:), debit(:), S, C1, C2, sdeb, sfin
real, pointer::resultat(:, :), pics(:, :), lacunes(:), time(:), episode(:), debis(:, :)
integer n, d1, f1, p, k, m, tpic, jo, bar,compt,i
integer, pointer::tepisode(:), chsigne(:), chsigne0(:)
real prems, maxi, toto(npt, 15), coeffdir, intercept, joey, starr, som, cpt

prems=1.*int(sdeb)
n=size(jour)
if(associated(time)) deallocate(time)
if(associated(lacunes)) deallocate(lacunes)
allocate(time(n), lacunes(n))
time=(/(1.*i, i=1, n)/)
call Saisonnalise(jour, debit, min(sdeb, sfin), max(sdeb, sfin), debis)

if(sdeb<sfin) then
	call CumulLacunes(debis(:, 1), lacunes)
	Call ExtraitPicsC3(time, debis(:, 1), S, C1, C2, pics)
else
	call CumulLacunes(debis(:, 2), lacunes)
	Call ExtraitPicsC3(time, debis(:, 2), S, C1, C2, pics)
endif

if(.not.associated(pics)) return
p=size(pics(:, 1))

do compt=1, p
	toto(compt, 13)=0.
	toto(compt, 12)=pics(compt, 1)
	k=int(pics(compt, 1))
	maxi=pics(compt, 2)
	!Extraction des variables
	!rang
	toto(compt, 1)=compt
	!Annee
	toto(compt, 2)=annee(time(k))
	!Jour
	toto(compt, 3)=jour(time(k))
	!valeur du maximum
	toto(compt, 4)=0.001*maxi
	!date angulaire
	toto(compt, 5)=DateAngulaire(jour(time(k)), prems)
	!cumul corrigé
	toto(compt, 6)=time(k)-lacunes(k)

	! extraction de l'épisode
	m=k+1
	do while (debit(m)>maxi/2.)
		m=m+1
		if(debit(m+1)==-9999.) toto(compt, 13)=1.
	enddo
	f1=m-1
	m=max(k-1,1)
	do while (debit(m)>maxi/2.)
		m=m-1
		if(m<=0) exit
		if(m>=2) then
		    if(debit(m-1)==-9999.) toto(compt, 13)=1.
		endif
	enddo
	d1=m+1

	!calcul de la durée sup- 0.5*pic
	if(d1>1) then
	    coeffdir=(debit(d1)-debit(d1-1))/(time(d1)-time(d1-1))
	    intercept=debit(d1)-coeffdir*time(d1)
	    starr=(0.5*maxi-intercept)/coeffdir
	else
	    starr=time(d1)
	endif
	coeffdir=(debit(f1)-debit(f1+1))/(time(f1)-time(f1+1))
	intercept=debit(f1)-coeffdir*time(f1)
	joey=(0.5*maxi-intercept)/coeffdir
	toto(compt, 7)=joey-starr
	!Temps de montée
	toto(compt, 8)=time(k)-starr
	!Temps de redescente
	toto(compt, 9)=joey-time(k)
	!Coeffs de forme
	if(associated(episode)) deallocate(episode)
	if(associated(tepisode)) deallocate(tepisode)
	if(associated(chsigne0)) deallocate(chsigne0)
	if(associated(chsigne)) deallocate(chsigne)
	allocate(episode(f1-d1+1), tepisode(f1-d1+1), chsigne0(f1-d1+2), chsigne(f1-d1+2))
	episode=debit(d1:f1)
	tepisode=int(time(d1:f1))
	chsigne0(1)=1
	chsigne0(f1-d1+2)=-1
	if(f1/=d1) then
		do i=2,f1-d1+1
			chsigne0(i)=signe(episode(i)-episode(i-1))
		enddo
	endif
	jo=0
	bar=0
	do i=1, f1-d1+2
		if(chsigne0(i)/=0) then
			jo=jo+1
			chsigne(jo)=chsigne0(i)
		else
			bar=bar+1
			chsigne(f1-d1+2-bar+1)=-1
		endif
	enddo

	!Ecart-type
	!toto(compt, 10)=maxi/(joey-starr)
	tpic=k-d1+1
	som=0.
	do i=1,f1-d1+1
		som=som+(1.*(tepisode(i)-tepisode(tpic))**2)*episode(i)
	enddo
	som=som+((starr-1.*tepisode(tpic))**2)*0.5*maxi+((joey-1.*tepisode(tpic))**2)*0.5*maxi
	toto(compt, 10)=sqrt(som/(sum(episode)+maxi))
	!Asymetrie
	som=0.
	do i=1,f1-d1+1
		som=som+(4.*(tepisode(i)-tepisode(tpic))**3)*episode(i)
	enddo
	som=som+((starr-1.*tepisode(tpic))**3)*0.5*maxi*4.+((joey-1.*tepisode(tpic))**3)*0.5*maxi*4.
toto(compt, 11)=(joey-time(k))/(joey-starr)
	!Complexité
	cpt=0.
	do i=2, f1-d1+2
		if(abs(chsigne(i)-chsigne(i-1))==2) then
			cpt=cpt+1.
		endif
	enddo
	toto(compt, 14)=1.*cpt
	!Volume
	som=0.
	do i=1,f1-d1
		som=som+0.5*0.001*(episode(i)+episode(i+1))*(tepisode(i+1)-tepisode(i))
	enddo
	som=som+(0.5*0.001*(0.5*maxi+episode(1))*(tepisode(1)-starr))&
			+(0.5*0.001*(0.5*maxi+episode(f1-d1+1))*(joey-tepisode(f1-d1+1)))
	toto(compt, 15)=0.000001*60.*60.*24.*som
	deallocate(episode, tepisode, chsigne, chsigne0)
enddo

if(associated(resultat)) deallocate(resultat)
allocate(resultat(p, 15))
resultat=toto(1:p, :)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Saisonnalise(jour, debit, sdeb, sfin, debis)

! A partir de la chronique débit, crée 2 chroniques saisonnalisées, contenues dans 'debis'
! sdeb<sfin sont les limites des saisons (jours)

real, intent(in)::debit(:), jour(:), sdeb, sfin
real, pointer::debis(:, :)
integer nd,debi


nd=size(debit)
if(associated(debis)) deallocate(debis)
allocate(debis(nd, 2))
	do debi=1, nd
		if(jour(debi)>=int(sdeb)) then
			if(jour(debi)<=int(sfin)) then
				debis(debi, 1)=debit(debi)
				debis(debi, 2)=-9999.
			else
				debis(debi, 2)=debit(debi)
				debis(debi, 1)=-9999.
			endif
		else
			debis(debi, 2)=debit(debi)
			debis(debi, 1)=-9999.
		endif
	enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CumulLacunes(debit, lacunes)
! Calcule, à chaque pas de temps, la durée cumulée des lacunes

real, intent(in)::debit(:)
real, pointer::lacunes(:)
integer n, k,i

n=size(debit)
if(associated(lacunes)) deallocate(lacunes)
allocate(lacunes(n))

k=0
do i=1, n
	if (debit(i)==-9999.) k=k+1
	lacunes(i)=k
enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! La variante C3 consiste à appliquer un troisième critère d'indépendance:
! sur l'évènement (dépassement de la moitié du pic), le débit ne doit pas repasser
! au dessus du pic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtraitPicsC3(t, Xt, S, C1, C2, resultat)

! Extrait tous les pics supérieurs à un seuil S, avec les contraintes C1 et C2
! le tableau result est de la forme (en colonne): (date (pic), valeur(pic))

integer, PARAMETER::npt=int(1.E3)
real, intent(in)::t(:), Xt(:), S, C1, C2
real, pointer:: resultat(:, :), provisoire(:, :), provisoire2(:, :), tempo(:), debito(:),&
				joe(:), provisoire3(:, :), provisoire4(:, :)
integer, pointer:: indexo(:)
integer n, avirer(npt), compt, decalage, j, pb, m, nc3, debut, fin,k,ii,l
real mini, toto(npt, 2), seuil

call ExtraitPicsSansContraintes(t, Xt, S, provisoire)

if(.not.associated(provisoire)) return
! Application de la contrainte de redescente
call indexs(provisoire(:, 1), indexo)
n=size(indexo)
if(associated(tempo)) deallocate(tempo)
if(associated(debito)) deallocate(debito)
Allocate(tempo(n), debito(n))
tempo=provisoire(indexo, 1)
debito=provisoire(indexo, 2)

compt=0
decalage=1
j=1
do while (j+decalage<=n)
	seuil=C1*debito(j)
	mini=minimum(Xt(int(tempo(j)):int(tempo(j+decalage))))
	if(mini<seuil) then
		pb=0
	else
		pb=1
	endif
	if (pb==0) then
		j=j+decalage
		decalage=1
	else
		compt=compt+1
		if(debito(j)<debito(j+decalage)) then
			avirer(compt)=j
			j=j+decalage
			decalage=1
		else
			avirer(compt)=j+decalage
			decalage=decalage+1
		endif
	endif
enddo

! Nouvel échantillon après application de la première contrainte
if(associated(provisoire2)) deallocate(provisoire2)
if(associated(joe)) deallocate(joe)
Allocate(provisoire2(n-compt, 2))
Allocate(joe(compt))
call tri(1.*avirer(1:compt), joe)

if(compt>0) then
	j=1
	m=1
	do k=1, n
		if(k==int(joe(j))) then
			if(j+1<=size(joe)) j=j+1
		else
			provisoire2(m, 1)=tempo(k)
			provisoire2(m, 2)=debito(k)
			m=m+1
		endif
	enddo
else
	provisoire2=provisoire
endif

! Application de la contrainte d'espacement temporel
deallocate(indexo, tempo, debito)
call indexs(provisoire2(:, 2), indexo)
n=size(indexo)

if(associated(tempo)) deallocate(tempo)
if(associated(debito)) deallocate(debito)
Allocate(tempo(n), debito(n))
tempo=provisoire2(indexo, 1)
debito=provisoire2(indexo, 2)
tempo=tempo(n:1:-1)
debito=debito(n:1:-1)

toto(1, 1)=tempo(1)
toto(1, 2)=debito(1)
compt=1
if (n>1) then
	do k=2, n
		pb=0
		do l=1, compt
			if(abs(tempo(k)-toto(l, 1))<C2) pb=1
		enddo
		if (pb==0) then
			compt=compt+1
			toto(compt, 1)=tempo(k)
			toto(compt, 2)=debito(k)
		endif
	enddo
endif
! Nouvel échantillon après C2
if(associated(provisoire3)) deallocate(provisoire3)
allocate(provisoire3(compt, 2))
provisoire3(:, :)=toto(1:compt, :)


! Application de la contrainte C3
deallocate(indexo, tempo, debito)
n=size(provisoire3(:, 1))
if(associated(tempo)) deallocate(tempo)
if(associated(debito)) deallocate(debito)
if(associated(provisoire4)) deallocate(provisoire4)
allocate(tempo(n), debito(n), provisoire4(n, 2))
tempo=provisoire3(:, 1)
debito=provisoire3(:, 2)

nc3=0
do ii=1, n
	!extraction d'un épisode
	k=int(tempo(ii))
	m=k+1
	do while (Xt(m)>=0.5*Xt(k))
		m=m+1
	enddo
	fin=m-1
	m=max(k-1,1)
	do while (Xt(m)>=0.5*Xt(k))
		m=m-1
		if(m<=0) exit
	enddo
	debut=m+1
	if(maxval(Xt(debut:fin))==Xt(k)) then
		nc3=nc3+1
		provisoire4(nc3, :)=provisoire3(ii, :)
	endif
enddo
if(associated(resultat)) deallocate(resultat)
allocate(resultat(nc3, 2))
resultat=provisoire4(1:nc3, :)
deallocate(tempo, debito, provisoire4)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tri(tin, tout)
!Renvoie le tableau tin trié dans l'ordre croissant!

real, intent(in)::tin(:)
real, pointer::tout(:), bret(:)
integer N, ind,k
real temp

N=size(tin)
if(associated(tout)) deallocate(tout)
if(associated(bret)) deallocate(bret)
Allocate(tout(N), bret(N))
bret=tin
do k=1, N
	ind=minval(minloc(bret(k:N)))
	tout(k)=bret(ind+k-1)
	temp=bret(ind+k-1)
	bret(ind+k-1)=bret(k)
	bret(k)=temp
end do

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine ExtraitPicsSansContraintes(t, Xt, S, resultat)

! Extrait tous les pics supérieurs à un seuil S, sans contraintes d'indépendance)
! le tableau result est de la forme (en colonne): (date (pic), valeur(pic))

integer, PARAMETER::npt=int(1.E3)
real, intent(in)::t(:), Xt(:), S
real, pointer::resultat(:, :), time(:), dis(:), titi(:)
integer n, rang(1), k, compt, m, debut, fin
real maxi, toto(npt, 2)

n=size(t)

if(associated(time)) deallocate(time)
if(associated(dis)) deallocate(dis)
allocate(time(n), dis(n))
time=t
dis=Xt
rang=maxloc(dis)
k=rang(1)
maxi=dis(k)
compt=0

do while (maxi>=S)
compt=compt+1
!extraction d'un épisode
m=k+1
do while (dis(m)>=S)
	m=m+1
enddo
fin=min(m-1,n)
m=max(k-1,1)
do while (dis(m)>=S)
	m=m-1
	if(m<=0) exit
enddo
debut=m+1

toto(compt, 1)=time(k)
toto(compt, 2)=maxi

!retrait de l'épisode de la chronique
if(associated(titi)) deallocate(titi)
allocate(titi(size(time)))
titi=time
deallocate(time)
allocate(time(size(titi)-(fin-debut+1)))
time=(/titi(1:debut-1), titi(fin+1:)/)
deallocate(titi)
allocate(titi(size(dis)))
titi=dis
deallocate(dis)
allocate(dis(size(titi)-(fin-debut+1)))
dis=(/titi(1:debut-1), titi(fin+1:)/)
deallocate(titi)

rang=maxloc(dis)
k=rang(1)
maxi=dis(k)

enddo

if (compt==0) return
if(associated(resultat)) deallocate(resultat)
allocate(resultat(compt, 2))
resultat=toto(1:compt, :)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function minimum(X)
! renvoie le minimum non nul d'une chronique avec lacunes

real, intent(in)::X(:)
integer n,i
real mini

mini=1.E10
n=size(X)
do i=1, n
	if(X(i)>0) then
		if(X(i)<mini) mini=X(i)
	endif
enddo
minimum=mini
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function DateAngulaire(jour, debut)
! convertit une date en date angulaire, à partir du début
! de l'année hydrologique

real, parameter::pi=3.14159265
real, intent(in)::jour, debut
real day

if(jour>=debut) then
	day=jour-debut+1
else
	day=365.25-debut+1.+jour
endif
DateAngulaire=(2.*pi*day)/365.25

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function signe(x)

real, intent(in)::x
integer signe

if(x==0) then
	signe=0
else if(x>0) then
	signe=1
else if (x<0) then
	signe=-1
end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine indexs(tin, tout)

!Renvoie le tableau des index correpondant au tableau tin!

real, intent(in)::tin(:)
integer, pointer::tout(:)
real, pointer:: chemoul(:)
integer N, ind, temp2,k,l
real temp

N=size(tin)
if(associated(tout)) deallocate(tout)
if(associated(chemoul)) deallocate(chemoul)
Allocate(tout(N), chemoul(N))
tout=(/(l, l=1, N)/)
chemoul=tin
do k=1, N
	ind=minval(minloc(chemoul(k:N)))
	temp=chemoul(ind+k-1)
	chemoul(ind+k-1)=chemoul(k)
	chemoul(k)=temp
	temp2=tout(ind+k-1)
	tout(ind+k-1)=tout(k)
	tout(k)=temp2
end do

end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module POTwrapper
