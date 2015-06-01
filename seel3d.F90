!*************************************************************************
!*************************************************************************
!*****
!***** SEEL3D.f90
!***** Simulation des Equations d'Euler Linearisees 3D
!*****
!*****
!***** Version: 1.3
!*****
!*************************************************************************
!*************************************************************************
!
!
!*************************************************************************
!****  LISTE DES MODULES
!*************************************************************************
!
module mod_grille
  implicit none
  real, dimension(:), allocatable :: xg,dxg
  real, dimension(:), allocatable :: yg,dyg
  real, dimension(:), allocatable :: zg,dzg
  real :: deltax,deltay,deltaz,deltat,CFL
  integer :: nymin,nymax,nxmax,nxray,nyray,nzray
  integer :: nt0,nx,ny,nz,ntfin,record_step
  real :: time
end module mod_grille
!
module mod_scheme
  implicit none
  real a(-3:3),a24(7),a15(7),a06(7)
  real rk(5),ck(5)
  real dfc2(1:3),dfc20(1:3),dfc02(1:3),dfc4(1:5),dfc13(1:5),dfc31(1:5)
  real dfc04(1:5),dfc40(1:5),dfc6(-3:3)
  real dfilt8(-4:4),dfilt6(-3:3),dfilt4(-2:2),dfilt2(-1:1)
  integer :: numrk,nrk
  real :: pi=acos(-1.)
end module mod_scheme
!
module mod_filtrage
  implicit none
  real :: cfx,cfx6,cfx4,cfx2,cfx1
  integer :: xfmin_x,xfmax_x,yfmin_x,yfmax_x,zfmin_x,zfmax_x
  real :: cfy,cfy6,cfy4,cfy2,cfy1
  integer ::xfmin_y,xfmax_y,yfmin_y,yfmax_y,zfmin_y,zfmax_y
  real :: cfz,cfz6,cfz4,cfz2,cfz1
  integer ::xfmin_z,xfmax_z,yfmin_z,yfmax_z,zfmin_z,zfmax_z
end module mod_filtrage
!
module mod_condlim
  implicit none
  integer :: ibc_right
end module mod_condlim
!
module mod_options
  implicit none
  integer :: icas=100
  logical :: o_reduction_order_bords
  logical :: o_damping,o_damping_sup
  logical :: o_restart=.false.
  character (len=40) :: filename_load_restart="bin_restart_run01.bin"
  integer :: irun=1         !! numero du run, valeur par defaut ici, modifiee si restart 
end module mod_options
!
module mod_scales
  implicit none
  real :: Lref
end module mod_scales
!
module mod_onde
  implicit none
  real :: omega,amp,alpha,mo
end module mod_onde
!
module mod_vectors
  use mod_grille
  implicit none
  real, dimension(:,:,:,:), allocatable :: U,Un,Ut,E,F,H,S,G
  real, dimension(:,:,:), allocatable :: rhoo,uo,vo,wo,po,coo
  real, dimension(:,:,:), allocatable :: duox,duoy,duoz,dvox,dvoy,dvoz,dwox,dwoy,dwoz,dpox,dpoy,dpoz
  real :: gamma,rgp,theta,phi
  real, dimension(:,:,:,:), allocatable :: VORT
  real, dimension(:), allocatable :: residu
  real, dimension(:,:,:,:,:),allocatable :: rect
  integer, dimension(:),allocatable :: xrect,yrect,zrect
end module mod_vectors
!
module mod_record
  implicit none
  real :: record_dummy=123.321
  real :: read_dummy
  logical :: o_record_vort=.false.
  logical :: o_rect=.true.
  character (len=40) :: record_filename_champs="bin_champs2d"
  character (len=40) :: record_filename_vort="bin_vort2d"
  character (len=40) :: record_filename_ecoulmoy="bin_ecoulmoy"
  integer :: irecord=1   !!!valeur par defaut, modifiee si redemarrage
  integer :: nxrect,nyrect,nzrect
end module mod_record

!
!
!
!
!*************************************************************************
!****  DIFFERENTS CAS
!*************************************************************************
!  icas=100    Pulse de pression gaussien avec ou sans ecoulement
!  icas=120    Monopole sans ecoulement
!*************************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
program seel2d
  !
  !***** Programme principal
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_condlim
  use mod_options
  use mod_scales
  use mod_onde
  use mod_vectors
  implicit none
integer :: Time_1,clock_rate,Time_2
  !
  call setcas
  call inivar
  call coeff_schemas
  call maillage
  call ecoulmoyen
  call pastemps
  call affichage
  !
  write(6,*) '////////////////////////// DEBUT INTEGRATION ///////////////////////'
  call system_clock(Time_1,clock_rate)
  call integ
  CALL system_clock(Time_2)
  print*,'TEMPS DE CALCUL ',(Time_2-Time_1)*1./clock_rate
  call record_for_restart   
  call closevar
  write(6,*) '////////////////////////// FIN PROGRAMME ///////////////////////////'
  !
  stop
  !
end program seel2d
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine setcas
  !
  ! Regle les parametres des cas pre-existants
  !******************************************************************
  !******************************************************************
  use mod_grille
  use mod_scheme
  use mod_filtrage
  use mod_condlim
  use mod_options
  use mod_scales
  use mod_onde
  use mod_record
  !
  implicit none
  !
  !
  !!///// CONDITIONS PAR DEFAUT
  !
  ! Maillage et CL et RK4
  deltax=1
  deltay=deltax
  deltaz=deltax
  nrk=4
  CFL=1
  ntfin=500
  record_step=1
  ! Source 
  omega=2*pi/10
  amp=0.1
  alpha=log(2.)/4.
  ! Ecoulement
  mo= 0
  Lref=10
  ! booleen
  o_reduction_order_bords=.false.
  o_damping=.true.
  o_damping_sup=.true.
  o_rect=.true.
  ! Filtrage
  cfx=0.3
  cfx6=cfx
  cfx4=cfx
  cfx2=cfx
  cfx1=0.
  !
  cfy=cfx
  cfy6=cfx6
  cfy4=cfx4
  cfy2=cfx2
  cfy1=cfx1
  !
  cfz=cfy
  cfz6=cfx6
  cfz4=cfx4
  cfz2=cfx2
  cfz1=cfy1
  ! Conditions aux limites
  ibc_right=1
  !
  !
  !///// CONDITIONS MODIFIEES EN FONCTION DU CAS ETUDIE
  !
  !
  !!///// Pulse de pression avec ou sans ecoulement
  if (icas==100) then
     nx=26!101
     ny=26!101
     nz=26!101
     ntfin=200
     record_step=500
     mo=0.5
     amp=0.001
     alpha=log(2.)/9.
     o_damping=.true.
     cfx=0.3;cfy=cfx;cfz=cfx
     o_damping_sup=.true.
     cfx6=0.6 ;cfx4=0.6 ;cfx2=0.6 ;cfx1=0.
     cfy6=cfx6;cfy4=cfx4;cfy2=cfx2;cfy1=cfx1
     cfz6=cfx6;cfz4=cfx4;cfz2=cfx2;cfz1=cfx1
     o_rect=.true.
     o_reduction_order_bords=.true.
     ibc_right=1    ! 2 marche aussi
     !
     !
     !!///// Monopole oscillant sans ecoulement
  elseif (icas==120) then
     nx=101
     ny=55
     nz=75
     ntfin= 500 !!5000 pour voir assez loin
     record_step=25
     mo=0.5
     omega=2*pi/10
     amp=0.1
     alpha=log(2.)/4.
     o_damping = .true.
     o_damping_sup=.true.
     o_rect=.true.
  end if
  !
  !
  !///// LIMITES FILTRAGE
  !
  !
  xfmin_x=2
  xfmax_x=nx-1
  yfmin_x=-2
  yfmax_x=ny+3
  zfmin_x=-2
  zfmax_x=nz+3
  !
  xfmin_y=-2
  xfmax_y=nx+3
  yfmin_y=2
  yfmax_y=ny-1
  zfmin_y=-2
  zfmax_y=nz+3
  !
  xfmin_z=-2
  xfmax_z=nx+3
  yfmin_z=-2
  yfmax_z=ny+3
  zfmin_z=2
  zfmax_z=nz-1
  !
  !
  !
end subroutine setcas
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine affichage
  !
  ! Affichage de diverses donnees avant le demarrage du calcul
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_scales
  use mod_onde
  use mod_vectors
  use mod_grille
  use mod_scheme
  use mod_options
  use mod_condlim
  use mod_record
  !
  implicit none
  !
27 format(A39,I4)
28 format(A39,G11.4)
  !
  !
  !! ///// INFOS SUR LES PARMETRES GENERAUX
  !
  !
  write(6,*) '/////'
  write(6,*) '////////////////////////// PARAMETRES DU CALCUL /////////////////////////////////////////'
  write(6,*) '/////'
  write(6,27) ' Nombre de points du maillage     nx = ',nx
  write(6,27) ' Nombre de points du maillage     ny = ',ny
  write(6,27) ' Nombre de points du maillage     nz = ',nz
  write(6,27) ' Numero pas temporel debut       nt0 = ',nt0
  write(6,27) ' Numero pas temporel fin       ntfin = ',ntfin
  write(6,28) ' CFL                             CFL = ',CFL
  write(6,28) ' Pas de temps                 deltat = ',deltat
  write(6,28) ' Nombre de Mach de reference      Mo = ',mo
  write(6,27) ' Frequence enregistr     record_step = ',record_step
  write(6,28) ' gamma                               = ',gamma
  write(6,28) ' constante des gaz parfaits          = ',rgp
  write(6,28) ' Longueur d''onde                     = ',2.*pi/omega
  write(6,28) ' Amplitude onde acoustique         A = ',amp
  write(6,28) ' Produit k.dx                 k.dx   = ',omega
  write(6,28) ' Demi-largeur de la source        bs = ',  sqrt(log(2.)/alpha)
  write(6,*)
  write(6,*) 'Filtrage              : ',o_damping
  write(6,*) 'Filtrage sup bords    : ',o_damping_sup
  write(6,*) 'Reduction ordre bords : ',o_reduction_order_bords 
  write(6,*)
  write(6,27) ' Schema Runge-Kutta en temps :   nrk = ',nrk
  write(6,*)
  write(6,27) ' Point reference rayonnement   nxray = ',nxray
  write(6,27) '                               nyray = ',nyray
  write(6,27) '                               nzray = ',nzray
  write(6,*) 'Conditions aux limites :'
  write(6,*) '/////'
  !
  !
end subroutine affichage
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ecoulmoyen
  !
  !***** Definition de l'ecoulement moyen
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_options
  use mod_condlim
  use mod_scales
  use mod_onde
  use mod_vectors
  implicit none
  integer x,y,z,i
  logical :: gradient

  gamma=1.4
  rgp=286.8875
  gradient=.true.
  !
  !///// ECOULEMENT MOYEN SI NECESSAIRE PAR DEFAUT
  !
  rhoo=1./1.
  uo=mo
  vo=0.
  wo=0.
  po=1./gamma
  coo=1.
  !
  duox=0.
  duoy=0.
  duoz=0.
  !
  dvox=0.
  dvoy=0.
  dvoz=0.
  !
  dwox=0.
  dwoy=0.
  dwoz=0.
  !
  dpox=0.
  dpoy=0.
  dpoz=0.
  !
  gradient=.false.
  uo=mo
  !
  !
  !///// CALCUL NUMERIQUE DES GRADIENTS MOYENS
  !             (en chantier...)
  !
  if (gradient) then
     write(6,*) ' .. calcul numerique des gradients de vitesse '
     do x=1,nx
        do y=1,ny
           do z=1,nz
              do i=-3,3
                 duox(x,y,z) = duox(x,y,z) + dxg(x)*a(i)*uo(x+i,y,z)
                 duoy(x,y,z) = duoy(x,y,z) + dyg(y)*a(i)*uo(x,y+i,z)
                 duoz(x,y,z) = duoz(x,y,z) + dzg(y)*a(i)*uo(x,y,z+i)

                 dvox(x,y,z) = dvox(x,y,z) + dxg(x)*a(i)*vo(x+i,y,z)
                 dvoy(x,y,z) = dvoy(x,y,z) + dyg(y)*a(i)*vo(x,y+i,z)
                 dvoz(x,y,z) = dvoz(x,y,z) + dzg(y)*a(i)*vo(x,y,z+i)

                 dwox(x,y,z) = dwox(x,y,z) + dxg(x)*a(i)*wo(x+i,y,z)
                 dwoy(x,y,z) = dwoy(x,y,z) + dyg(y)*a(i)*wo(x,y+i,z)
                 dwoz(x,y,z) = dwoz(x,y,z) + dzg(y)*a(i)*wo(x,y,z+i)
              end do
           end do
        end do
     end do
     !
     do x=-2,0
        do y=1,ny
           do z=1,nz
              duox(x,y,z) =duox(1,y,z)
              duoy(x,y,z) =duoy(1,y,z)
              duoz(x,y,z) =duoz(1,y,z)

              dvox(x,y,z) =dvox(1,y,z)
              dvoy(x,y,z) =dvoy(1,y,z)
              dvoz(x,y,z) =dvoz(1,y,z)

              dwox(x,y,z) =dwox(1,y,z)
              dwoy(x,y,z) =dwoy(1,y,z)
              dwoz(x,y,z) =dwoz(1,y,z)
           end do
        end do
     end do
     !
     do x=nx+1,nx+3
        do y=1,ny
           do z=1,nz
              duox(x,y,z) =duox(nx,y,z)
              duoy(x,y,z) =duoy(nx,y,z)
              duoz(x,y,z) =duoz(nx,y,z)

              dvox(x,y,z) =dvox(nx,y,z)
              dvoy(x,y,z) =dvoy(nx,y,z)
              dvoz(x,y,z) =dvoz(nx,y,z)

              dwox(x,y,z) =dwox(nx,y,z)
              dwoy(x,y,z) =dwoy(nx,y,z)
              dwoz(x,y,z) =dwoz(nx,y,z)
           end do
        end do
     end do

  end if
  !
  !
end subroutine ecoulmoyen
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
 subroutine pastemps
  !
  ! Calcul du pas de temps (critere cfl)
  !
  !******************************************************************
  !******************************************************************
  use mod_vectors
  use mod_grille
  use mod_options
  implicit none
  real :: uopmax,vopmax,dxmin,dxmax,dymin,dymax
  real :: wopmax,dzmin,dzmax
  integer :: x,y,z
  !
  !
  uopmax = 0.
  vopmax = 0.
  wopmax = 0.
  dxmin = xg(nx) - xg(1)
  dymin = yg(ny) - yg(1)
  dzmin = zg(nz) - zg(1)
  dxmax = 0.
  dymax = 0.
  dzmax = 0.
  do x=-2,nx+2
     dxmin = min( dxmin,xg(x+1)-xg(x) )
     dxmax = max( dxmax,xg(x+1)-xg(x) )
     do y=-2,ny+2
        dymin = min( dymin,yg(y+1)-yg(y) )
        dymax = max( dymax,yg(y+1)-yg(y) )
        do z=-2,nz+2
           dzmin = min( dzmin,zg(z+1)-zg(z) )
           dzmax = max( dzmax,zg(z+1)-zg(z) )
           uopmax = max( uopmax,uo(x,y,z)+sqrt(gamma*po(x,y,z)*rhoo(x,y,z)) )
           vopmax = max( vopmax,vo(x,y,z)+sqrt(gamma*po(x,y,z)*rhoo(x,y,z)) )
           wopmax = max( wopmax,wo(x,y,z)+sqrt(gamma*po(x,y,z)*rhoo(x,y,z)) )
        end do
     end do
  end do
  !
  deltat=CFL*min(dxmin/uopmax,dymin/vopmax,dzmin/wopmax)

end subroutine pastemps
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine integ
  !
  !***** Integration temporelle
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_onde
  use mod_options
  use mod_condlim
  use mod_vectors
  use mod_record
  use mod_filtrage
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  integer :: itime,irk,x,y,z,j,nt
  !
  !
  call sauveparametre
  !
!$OMP PARALLEL 


!/// NON NECESSAIRE
#ifdef _OPENMP
nt=OMP_get_num_threads()
!$OMP single
print*, "Nombre de Thread : ",nt
!$OMP end single
#endif
!/// NON NECESSAIRE


!$OMP DO  
  do z=-2,nz+3
    Un(:,:,:,z)=U(:,:,:,z)
enddo
!$OMP end DO nowait
!$OMP BARRIER
  !
  do itime=nt0,ntfin

     call set_champ(itime) 
     call record(itime)           
     do irk=1,nrk

        call ts(irk)                ! calcul de S
        call fluxes                 ! Un -> E,F,G,H

        call ptsint(irk)            ! U +      E,F,G,H,S -> Ut
        !// FACES                   ! U + Un + E,F,G,H,S -> Ut
        call ptsright(irk,ibc_right)
        call ptsleft(irk)
        call ptsbot(irk)
        call ptstop(irk)
        call ptsfront(irk)
        call ptsback(irk)
        !// ARRETES
        call ptsbottomfront(irk)
        call ptsbottomleft(irk)
        call ptsleftfront(irk)
        call ptsbottomright(irk)
        call ptstopleft(irk)
        call ptstopfront(irk)
        call ptsrightfront(irk)
        call ptsbottomback(irk)
        call ptsleftback(irk)
        call ptstopback(irk)
        call ptsrightback(irk)
        call ptstopright(irk)
        !// COINS
        call pts_c_toprightfront(irk)
        call pts_c_toprightback(irk)
        call pts_c_bottomrightfront(irk)
        call pts_c_bottomrightback(irk)
        call pts_c_topleftfront(irk)
        call pts_c_topleftback(irk)
        call pts_c_bottomleftfront(irk)
        call pts_c_bottomleftback(irk)
!
!$OMP BARRIER
!$OMP DO  
  do z=-2,nz+3
    Un(:,:,:,z)=Ut(:,:,:,z)
enddo
!$OMP end DO nowait
!$OMP BARRIER
     end do

     if (o_damping) then
         call filtrage8x          !!rem: entree Ut --> sortie: Un=Ut
         call filtrage8y          !!rem: entree Ut --> sortie: Un=Ut
         call filtrage8z          !!rem: entree Ut --> sortie: Un=Ut
        !
        if(o_damping_sup) then
           call filtragex_sup !!rem: entree Ut --> sortie: Un=Ut
           call filtragey_sup !!rem: entree Ut --> sortie: Un=Ut
           call filtragez_sup !!rem: entree Ut --> sortie: Un
        end if
     end if
!$OMP DO 
  do z=-2,nz+3
    U(:,:,:,z)=Un(:,:,:,z)
enddo
!$OMP end DO nowait
!$OMP BARRIER
     time=time+deltat
  end do
!$OMP END PARALLEL
  print*,U(:,5,5,5)
  !
  !
end subroutine integ
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine inivar
  !
  !***** Initialisation des tableaux de travail et champs initiaux
  !******************************************************************
  !******************************************************************
  use mod_condlim
  use mod_options
  use mod_scales
  use mod_onde
  use mod_vectors
  use mod_record
  use mod_grille
  implicit none
  integer :: x,y,z
  integer :: nt0_save, ntfin_save, irun_save, irecord_save
  real :: time_save
  !
  !
  !///// ALLOCATIONS MEMOIRE
  !
  !
  allocate(xg(-2:nx+3))
  allocate(dxg(-2:nx+3))
  allocate(yg(-2:ny+3))
  allocate(dyg(-2:ny+3))
  allocate(zg(-2:nz+3))
  allocate(dzg(-2:nz+3))
  !
  allocate(rhoo(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(uo(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(vo(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(wo(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(po(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(coo(-2:nx+3,-2:ny+3,-2:nz+3))
  !
  allocate(duox(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(duoy(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(duoz(-2:nx+3,-2:ny+3,-2:nz+3))
  !
  allocate(dvox(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dvoy(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dvoz(-2:nx+3,-2:ny+3,-2:nz+3))
  !
  allocate(dwox(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dwoy(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dwoz(-2:nx+3,-2:ny+3,-2:nz+3))
  !
  allocate(dpox(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dpoy(-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(dpoz(-2:nx+3,-2:ny+3,-2:nz+3))
  !
  allocate( U(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(Un(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate(Ut(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate( E(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate( F(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate( G(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate( H(5,-2:nx+3,-2:ny+3,-2:nz+3))
  allocate( S(5,-2:nx+3,-2:ny+3,-2:nz+3))
  ! 
  if(o_record_vort) then
     allocate(VORT(3,-2:nx+3,-2:ny+3,-2:nz+3))
     VORT=0.
  end if
  !
  !
  !///// INITIALISATION ECOULEMENT MOYEN
  !
  !
  duox=0.
  duoy=0.
  duoz=0.
  !
  dvox=0.
  dvoy=0.
  dvoz=0.
  !
  dwox=0.
  dwoy=0.
  dwoz=0.
  !
  dpox=0.
  dpoy=0.
  dpoz=0.
  !
  rhoo=HUGE(1.)
  uo=0.
  vo=0.
  wo=0.
  po=0.
  !
  !
  !///// INITIALISATION DES TABLEAUX DE VARIABLES
  !
  !
  if(o_restart) then
     write(6,*) '***'
     write(6,*) '*** LECTURE D''UN FICHIER RESTART ***'
     write(6,*) '***'
     write(6,*) 'Fichier lu: ', filename_load_restart
     open(507,file=filename_load_restart,form='unformatted',status='unknown')
     read(507) read_dummy
     write(6,*) '   check: read_dummy=', read_dummy
     read(507) nt0_save
     read(507) ntfin_save
     read(507) irun_save
     read(507) irecord_save
     read(507) time_save
     read(507) read_dummy
     write(6,*) '   check: read_dummy=', read_dummy
     read(507) (((U(1,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     read(507) (((U(2,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     read(507) (((U(3,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     read(507) (((U(4,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     read(507) (((U(5,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     read(507) read_dummy
     write(6,*) '   check: read_dummy=', read_dummy
     close(507)
     !
     irun=irun_save+1
     nt0=ntfin_save+1
     ntfin=ntfin+ntfin_save+1
     irecord=irecord_save
     time=time_save
     Un=U
     Ut=U
     E=0
     F=0
     G=0
     H=0
     S=0
27   format(A39,I4)
     write(6,27) ' Numero pas temporel debut       nt0 = ',nt0
     write(6,27) ' Numero pas temporel fin       ntfin = ',ntfin
     write(6,27) ' Numero pas temporel fin        irun = ',irun
     write(6,*) '***'
     write(6,*) '***'

  else
     U=0
     Un=0
     Ut=0
     E=0
     F=0
     G=0
     H=0
     S=0
     nt0=0
     time=0
  endif
  !
  !
  !
  !
  !///// ALLOCATIONS MEMOIRE TABLEAUX POUR ENREGISTREMENT
  !
  !
  nxrect=2
  nyrect=2
  nzrect=2

  allocate(rect(nt0:ntfin,1:nxrect,1:nyrect,1:nzrect,1:5))      !! temps / x / y / z / numero_champ
  rect=0

  allocate(xrect(1:nxrect))
  allocate(yrect(1:nyrect))
  allocate(zrect(1:nzrect))

  xrect(1)=(nx+1)/2+10
  yrect(1)=(ny+1)/2+10
  zrect(1)=(nz+1)/2+10
  xrect(2)=nx
  yrect(2)=ny
  zrect(2)=nz
  !
  !
end subroutine inivar
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine closevar
  !
  !***** Deallocations, fermetures fichiers,...
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_vectors
  use mod_record
  implicit none
  !
  !///// DEALLOCATIONS MEMOIRE
  !
  deallocate(xg)
  deallocate(dxg)
  deallocate(yg)
  deallocate(dyg)
  deallocate(zg)
  deallocate(dzg)
  !
  deallocate(rhoo)
  deallocate(uo)
  deallocate(vo)
  deallocate(wo)
  deallocate(po)
  deallocate(coo)
  !
  deallocate(duox)
  deallocate(duoy)
  deallocate(duoz)
  !
  deallocate(dvox)
  deallocate(dvoy)
  deallocate(dvoz)
  !
  deallocate(dwox)
  deallocate(dwoy)
  deallocate(dwoz)
  !
  deallocate(dpox)
  deallocate(dpoy)
  deallocate(dpoz)
  !
  deallocate(U)
  deallocate(Un)
  deallocate(Ut)
  deallocate(E)
  deallocate(F)
  deallocate(G)
  deallocate(H)
  deallocate(S)
  !
  !!/// Pulse pression ou tourbillon dans ecoulement uniforme
  !deallocate(residu)
  !
  if(o_rect) then
     deallocate(rect)
     deallocate(xrect)
     deallocate(yrect)
     deallocate(zrect)
     if(o_record_vort) then
        deallocate(VORT)
     end if
  end if
  !
end subroutine closevar
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine set_champ(itime)
  !
  !***** Champ initial (pulse de pression,...)
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_onde
  use mod_vectors
  use mod_grille
  implicit none
  integer :: itime
  integer :: x,y,z
  real    :: arg
  !
  !
  !///// CHAMP INITIAL
  !
  !
  !!///// Pulse de pression
  if (itime.eq.0) then
     if (icas==100) then
!$OMP DO  
        do z=1,nz
           do y=1,ny
!$OMP SIMD private(arg)
              do x=1,nx
                 arg = (xg(x)+0.)**2 + yg(y)**2 + zg(z)**2
                 U(1,x,y,z) = amp * exp( - alpha*arg )
                 U(5,x,y,z) = U(1,x,y,z)
              end do
           end do
        end do
!$OMP end DO nowait
!$OMP BARRIER 
!$OMP DO 
  do z=-2,nz+3
        Un(:,:,:,z)=U(:,:,:,z)
enddo
!$OMP end DO nowait
!$OMP BARRIER
     end if
  end if
  !
  !
end subroutine set_champ
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsint(irk)
  !
  !***** Calcul des points centraux (interieurs)
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_vectors
  implicit none
  real :: CC(5),DDx(5),DDy(5),DDz(5)
  integer :: irk,x,y,z
  !
!$OMP DO  collapse(2)
     do z=1,nz
        do y=1,ny
!$OMP SIMD private(ddx,ddy,ddz,cc)
           do x=1,nx
              DDx = a(3)*(E(:,x+3,y,z) - E(:,x-3,y,z))    &
                   +a(2)*(E(:,x+2,y,z) - E(:,x-2,y,z))    &
                   +a(1)*(E(:,x+1,y,z) - E(:,x-1,y,z))
              ! 
              DDy = a(3)*(F(:,x,y+3,z) - F(:,x,y-3,z))    &
                   +a(2)*(F(:,x,y+2,z) - F(:,x,y-2,z))    &
                   +a(1)*(F(:,x,y+1,z) - F(:,x,y-1,z))
              ! 
              DDz = a(3)*(G(:,x,y,z+3) - G(:,x,y,z-3))    &
                   +a(2)*(G(:,x,y,z+2) - G(:,x,y,z-2))    &
                   +a(1)*(G(:,x,y,z+1) - G(:,x,y,z-1))
              !
              CC = - dxg(x)*DDx - dyg(y)*DDy - dzg(z)*DDz - H(:,x,y,z) + S(:,x,y,z)
              Ut(:,x,y,z) = U(:,x,y,z) + deltat*CC*rk(irk)
           end do
        end do
     end do
!$OMP end DO nowait
  !
end subroutine ptsint
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsright(irk,ibc)
!
!***** Traitement des points de la face droite (x>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_vectors
   implicit none
   integer :: irk,ibc,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=0
dirx=1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=z
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
            if(ibc==1) then
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
            elseif(ibc==2) then
               DD(5) = vray * (xc*Dx(5) + yc*Dy(5)  + zc*Dz(5)  + Un(5,x,y,z)*r)
               DD(2) = uo(x,y,z) * Dx(2) + vo(x,y,z)*Dy(2) + wo(x,y,z)*Dz(2) + Dx(5)*rhoo(x,y,z)
               DD(3) = uo(x,y,z) * Dx(3) + vo(x,y,z)*Dy(3) + wo(x,y,z)*Dz(3) + Dy(5)*rhoo(x,y,z)
               DD(4) = uo(x,y,z) * Dx(4) + vo(x,y,z)*Dy(4) + wo(x,y,z)*Dz(4) + Dz(5)*rhoo(x,y,z)
               DD(1) = uo(x,y,z) * Dx(1) + vo(x,y,z)*Dy(1) + wo(x,y,z)*Dz(1) &
                     -(uo(x,y,z) * Dx(5) + vo(x,y,z)*Dy(5) + wo(x,y,z)*Dz(5) - DD(5) ) / coo(x,y,z)**2
            end if
!
             Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsright
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsleft(irk)
!
!***** Traitement des points de la face gauche
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_onde
   use mod_vectors
   implicit none
   integer :: irk,j,x,y,z,i
   real :: DD(5),r,r2d,xc,yc,zc,vray,co,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=0
dirx=-1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=z
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsleft
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptstop(irk)
!
!***** Traitement des points de la face haute
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=1
dirx=0

!$OMP DO  collapse(2)
       do z=1,nz
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=z
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
! 
end subroutine ptstop
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsbot(irk)
!
!***** Points de la face du bas
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=-1
dirx=0

!$OMP DO  collapse(2)
       do z=1,nz
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=z
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
!
end subroutine ptsbot
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsback(irk)
!
!***** Points de la face arriere
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=0
dirx=0

!$OMP DO  collapse(2)
       do z=-2,0
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=1
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
!
end subroutine ptsback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsfront(irk)
!
!***** Points de la face avant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=0
dirx=0

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=nz
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsbottomleft(irk)
!
!***** Points de l'arrete en bas a gauche
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=-1
dirx=-1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=z
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsbottomleft
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsbottomback(irk)
!
!***** Points de l'arrete bas derriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=-1
dirx=0

!$OMP DO  collapse(2)
       do z=-2,0
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=1
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsbottomback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsleftback(irk)
!
!***** Points de l'arrete gauche arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real ::DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=0
dirx=-1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=1
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsleftback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptstopright(irk)
!
!***** Points de l'arrete haut droit
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=1
dirx=1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=z
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptstopright
!******************************************************************
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsbottomright(irk)
!
!***** Points de l'arrete en bas a droite
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=-1
dirx=1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=z
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsbottomright
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptstopleft(irk)
!
!***** Points de l'arrete haut gauche
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=0
diry=1
dirx=-1

!$OMP DO  collapse(2)
       do z=1,nz
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=z
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
!
end subroutine ptstopleft
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptstopback(irk)
!
!***** Points de l'arrete haut arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=1
dirx=0

!$OMP DO  collapse(2)
       do z=-2,0
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=1
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
!
end subroutine ptstopback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsrightback(irk)
!
!***** Points de l'arrete droite arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=0
dirx=1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=1
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsrightback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsrightfront(irk)
!
!***** Points de l'arrete droite devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=0
dirx=1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=nz
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsrightfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsleftfront(irk)
!
!***** Points de l'arrete droite devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=0
dirx=-1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=1,ny
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=nz
         y3=y ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsleftfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptsbottomfront(irk)
!
!***** Points de l'arrete bas devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=-1
dirx=0

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=nz
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!      
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptsbottomfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ptstopfront(irk)
!
!***** Points de l'arrete haut devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=1
dirx=0

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=1,nx

         z3=nz
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=x

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine ptstopfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_bottomleftback(irk)
!
!***** Points du coin bas gauche arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=-1
dirx=-1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=1
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
	 end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_bottomleftback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_bottomrightback(irk)
!
!***** Points du coin bas droit arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=-1
dirx=1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=1
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_bottomrightback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_toprightback(irk)
!
!***** Points du coin haut droits arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=1
dirx=1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=1
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_toprightback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_topleftback(irk)
!
!***** Points du coin haut gauche arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=-1
diry=1
dirx=-1

!$OMP DO  collapse(2)
       do z=-2,0
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=1
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
	    end do
          end do
        end do
!$OMP end DO nowait
!
end subroutine pts_c_topleftback
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_bottomleftfront(irk)
!
!***** Points du coin gauche bas devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=-1
dirx=-1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=nz
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_bottomleftfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_bottomrightfront(irk)
!
!***** Points du coin bas droit avant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=-1
dirx=1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=-2,0
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd)
         do x=nx+1,nx+3

         z3=nz
         y3=1 ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_bottomrightfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_toprightfront(irk)
!
!***** Points du coin haut droit devant (z>0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=1
dirx=1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=nx+1,nx+3

         z3=nz
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=nx

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_toprightfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine pts_c_topleftfront(irk)
!
!***** Points du coin gauche haut arriere (z<0)
!******************************************************************
!******************************************************************
   use mod_scheme
   use mod_condlim
   use mod_vectors
   implicit none
   integer :: irk,i,j,x,y,z
   real :: DD(5),r,r2d,xc,yc,zc,vray,costheta,sintheta,sinphi,cosphi
   real, dimension(5) :: Dx,Dy,Dz
   real :: uer, uetheta, uephi
  real,dimension(-3:3,-3:3)      :: ax,ay,az
  integer :: dirx,diry,dirz,x1,x2,y1,y2,z1,z2,x3,y3,z3
!
!
  ax(-3:-1   , 0)=-a(3:1:-1) ; ax(  0     , 0)= 0.  ; ax(1:3    , 0)= a(1:3)
  ay(-3:-1   , 0)=-a(3:1:-1) ; ay(  0     , 0)= 0.  ; ay(1:3    , 0)= a(1:3) ! Rearangement des coeffs
  az(-3:-1   , 0)=-a(3:1:-1) ; az(  0     , 0)= 0.  ; az(1:3    , 0)= a(1:3)
  ax(  :     , 1)= a24       ; ax(  :     , 2)= a15 ; ax( :     , 3)= a06
  ay(  :     , 1)= a24       ; ay(  :     , 2)= a15 ; ay( :     , 3)= a06
  az(  :     , 1)= a24       ; az(  :     , 2)= a15 ; az( :     , 3)= a06
  ax( 3:-3:-1,-1)=-a24       ; ax( 3:-3:-1,-2)=-a15 ; ax(3:-3:-1,-3)=-a06
  ay( 3:-3:-1,-1)=-a24       ; ay( 3:-3:-1,-2)=-a15 ; ay(3:-3:-1,-3)=-a06
  az( 3:-3:-1,-1)=-a24       ; az( 3:-3:-1,-2)=-a15 ; az(3:-3:-1,-3)=-a06
!

dirz=1
diry=1
dirx=-1

!$OMP DO  collapse(2)
       do z=nz+1,nz+3
       do y=ny+1,ny+3
!$OMP SIMD private(x3,y3,z3,Dx,dy,dz,xc,yc,zc,r,r2d,costheta,sintheta,cosphi,sinphi,uer,uetheta,uephi,vray,dd) 
         do x=-2,0

         z3=nz
         y3=ny ! x3,y3 et z3 permettent de gerer le decentrement
         x3=1

       Dx=0.
       Dy=0.
       Dz=0.

         do j=-3,3
           Dx = Dx - ax(j,x-x3)*Un(:,x3-j,y   ,z   )
           Dy = Dy - ay(j,y-y3)*Un(:,x   ,y3-j,z   )
           Dz = Dz - az(j,z-z3)*Un(:,x   ,y   ,z3-j)
         enddo

           Dx = Dx* dxg(x)
           Dy = Dy* dyg(y)
           Dz = Dz* dxg(z)
!
	    xc = xg(x) - xg(nxray)
	    yc = yg(y) - yg(nyray)
	    zc = zg(z) - zg(nzray)
!
             r  =1./sqrt(xc**2 + yc**2 + zc**2)
	    r2d= sqrt(xc**2 + yc**2) 
             costheta =  zc*r
             sintheta = r2d*r
            if(r2d.ne.0) then
               cosphi = xc/r2d
               sinphi = yc/r2d
            else                            !Dans ce cas d/dr=d/dz
               cosphi=1.                     !car alors r2d=0 et sintheta=0.
               sinphi=1.                     !Donc peu importe la valeur (finie) de cosphi et sinphi.
            end if
!
	    xc = sintheta*cosphi
	    yc = sintheta*sinphi
	    zc = costheta
!
            uer     =  uo(x,y,z)*xc + vo(x,y,z)*yc + wo(x,y,z)*zc
             uetheta =  (uo(x,y,z)*cosphi + vo(x,y,z)*sinphi)*costheta - wo(x,y,z)*sintheta
             uephi   =   vo(x,y,z)*cosphi - uo(x,y,z)*sinphi
            vray=uer+ sqrt(coo(x,y,z)**2 -uetheta**2 -uephi**2)
!
                DD = vray * (xc*Dx + yc*Dy  + zc*Dz  + Un(:,x,y,z)*r)
!
	      Ut(:,x,y,z) = U(:,x,y,z) - DD*deltat*rk(irk)
!
         end do
      end do
   end do
!$OMP end DO nowait
!
end subroutine pts_c_topleftfront
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtrage8x
!
!***** Filtrage des champs dans la direction x
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!///// FILTRAGE EN X DES POINTS INTERIEURS
!
!
!$OMP DO  collapse(3)
      do z=zfmin_x,zfmax_x
         do y=yfmin_x,yfmax_x
            do x=xfmin_x,xfmax_x
               Un(:,x,y,z) = Ut(:,x,y,z) - cfx*                   &
                           ( dfilt8(4)*(Ut(:,x+4,y,z)+Ut(:,x-4,y,z))  &
                            +dfilt8(3)*(Ut(:,x+3,y,z)+Ut(:,x-3,y,z))  &
                            +dfilt8(2)*(Ut(:,x+2,y,z)+Ut(:,x-2,y,z))  &
                            +dfilt8(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                            +dfilt8(0)*Ut(:,x,y,z) )
         end do
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO   collapse(2)
  do z=zfmin_x,zfmax_x
     do y=yfmin_x,yfmax_x
!$OMP SIMD
        do x=xfmin_x,xfmax_x
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtrage8x
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtrage8y
!
!***** Filtrage des champs dans la direction y
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!///// FILTRAGE EN Y DES POINTS INTERIEURS
!
!
!$OMP DO  collapse(3)
      do z=zfmin_y,zfmax_y
         do y=yfmin_y,yfmax_y
            do x=xfmin_y,xfmax_y
               Un(:,x,y,z) = Ut(:,x,y,z) - cfy*                  &
                          ( dfilt8(4)*(Ut(:,x,y+4,z)+Ut(:,x,y-4,z))  &
                           +dfilt8(3)*(Ut(:,x,y+3,z)+Ut(:,x,y-3,z))  &
                           +dfilt8(2)*(Ut(:,x,y+2,z)+Ut(:,x,y-2,z))  &
                           +dfilt8(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                           +dfilt8(0)*Ut(:,x,y,z) )
         end do
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_y,zfmax_y
     do y=yfmin_y,yfmax_y
!$OMP SIMD 
        do x=xfmin_y,xfmax_y
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtrage8y
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtrage8z
!
!***** Filtrage des champs dans la direction z
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!///// FILTRAGE EN Z DES POINTS INTERIEURS
!
!
!$OMP DO  collapse(3)
      do z=zfmin_z,zfmax_z
         do y=yfmin_z,yfmax_z
            do x=xfmin_z,xfmax_z
               Un(:,x,y,z) = Ut(:,x,y,z) - cfz*                  &
                          ( dfilt8(4)*(Ut(:,x,y,z+4)+Ut(:,x,y,z-4))  &
                           +dfilt8(3)*(Ut(:,x,y,z+3)+Ut(:,x,y,z-3))  &
                           +dfilt8(2)*(Ut(:,x,y,z+2)+Ut(:,x,y,z-2))  &
                           +dfilt8(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                           +dfilt8(0)*Ut(:,x,y,z) )
         end do
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_z,zfmax_z
     do y=yfmin_z,yfmax_z
!$OMP SIMD 
        do x=xfmin_z,xfmax_z
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtrage8z
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtragex_sup
!
!***** Filtrage additionel en x sur les bords
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!// Points droits (x=nx ; x=nx+1 ; x=nx+2)
!
!
   x=nx
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx6*                         &
                          ( dfilt6(3)*(Ut(:,x+3,y,z)+Ut(:,x-3,y,z))  &
                           +dfilt6(2)*(Ut(:,x+2,y,z)+Ut(:,x-2,y,z))  &
                           +dfilt6(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   x=nx+1
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx4*                       &
                         ( dfilt4(2)*(Ut(:,x+2,y,z)+Ut(:,x-2,y,z))  &
                          +dfilt4(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   x=nx+2
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx2*                        &
                         ( dfilt2(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   x=nx+3
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx1*0.5*(Ut(:,x,y,z)-Ut(:,x-1,y,z))
      end do
   end do
!$OMP end DO nowait
!
!
!// Points gauches  (x=1 ; x=0 ; x=-1)
!
!
   x=1
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx6*                         &
                          ( dfilt6(3)*(Ut(:,x+3,y,z)+Ut(:,x-3,y,z))  &
                           +dfilt6(2)*(Ut(:,x+2,y,z)+Ut(:,x-2,y,z))  &
                           +dfilt6(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   x=0
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx4*                       &
                         ( dfilt4(2)*(Ut(:,x+2,y,z)+Ut(:,x-2,y,z))  &
                          +dfilt4(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   x=-1
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx2*                        &
                         ( dfilt2(1)*(Ut(:,x+1,y,z)+Ut(:,x-1,y,z))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
!
   x=-2
!$OMP DO  !collapse(2)
      do z=zfmin_x,zfmax_x
!$OMP SIMD 
         do y=yfmin_x,yfmax_x
            Un(:,x,y,z) = Ut(:,x,y,z) - cfx1*0.5*(Ut(:,x,y,z)-Ut(:,x+1,y,z))
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_x,zfmax_x
     do y=yfmin_x,yfmax_x
!$OMP SIMD 
        do x=-2,1
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_x,zfmax_x
     do y=yfmin_x,yfmax_x
!$OMP SIMD 
        do x=nx,nx+3
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtragex_sup
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtragey_sup
!
!***** Filtrage additionel en y sur les bords
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!// Points hauts (y=ny; y=ny+1; y=ny+2)
!
!
   y=ny
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy6*                         &
                          ( dfilt6(3)*(Ut(:,x,y+3,z)+Ut(:,x,y-3,z))  &
                           +dfilt6(2)*(Ut(:,x,y+2,z)+Ut(:,x,y-2,z))  &
                           +dfilt6(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=ny+1
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy4*                        &
                         ( dfilt4(2)*(Ut(:,x,y+2,z)+Ut(:,x,y-2,z))  &
                          +dfilt4(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=ny+2
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy2*                        &
                         ( dfilt2(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=ny+3
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy1*0.5*(Ut(:,x,y,z)-Ut(:,x,y-1,z))
      end do
   end do
!$OMP end DO nowait
!
!
!// Points bas (y=1; y=0; y=-1)
!
!
   y=1
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy6*                         &
                          ( dfilt6(3)*(Ut(:,x,y+3,z)+Ut(:,x,y-3,z))  &
                           +dfilt6(2)*(Ut(:,x,y+2,z)+Ut(:,x,y-2,z))  &
                           +dfilt6(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=0
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy4*                        &
                         ( dfilt4(2)*(Ut(:,x,y+2,z)+Ut(:,x,y-2,z))  &
                          +dfilt4(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=-1
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy2*                        &
                         ( dfilt2(1)*(Ut(:,x,y+1,z)+Ut(:,x,y-1,z))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   y=-2
!$OMP DO  !collapse(2)
      do z=zfmin_y,zfmax_y
!$OMP SIMD 
         do x=xfmin_y,xfmax_y
            Un(:,x,y,z) = Ut(:,x,y,z) - cfy1*0.5*(Ut(:,x,y,z)-Ut(:,x,y+1,z))
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_y,zfmax_y
     do y=-2,1
!$OMP SIMD 
        do x=xfmin_y,xfmax_y
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!$OMP DO  collapse(2)
  do z=zfmin_y,zfmax_y
     do y=ny,ny+3
!$OMP SIMD 
        do x=xfmin_y,xfmax_y
        Ut(:,x,y,z)=Un(:,x,y,z)
enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtragey_sup
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine filtragez_sup
!
!***** Filtrage additionel en z sur les bords
!******************************************************************
!******************************************************************
   use mod_filtrage
   use mod_vectors
   use mod_scheme
   implicit none
   integer :: x,y,z,i,j
!
!
!// Points avant (z=nz; z=nz+1; z=nz+2)
!
!
   z=nz
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz6*                         &
                          ( dfilt6(3)*(Ut(:,x,y,z+3)+Ut(:,x,y,z-3))  &
                           +dfilt6(2)*(Ut(:,x,y,z+2)+Ut(:,x,y,z-2))  &
                           +dfilt6(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=nz+1
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz4*                        &
                         ( dfilt4(2)*(Ut(:,x,y,z+2)+Ut(:,x,y,z-2))  &
                          +dfilt4(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=nz+2
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz2*                        &
                         ( dfilt2(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=nz+3
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz1*0.5*(Ut(:,x,y,z)-Ut(:,x,y,z-1))
         end do
   end do
!$OMP end DO nowait
!
!// Points arrieres (z=1; z=0; z=-1)
!
!
   z=1
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz6*                         &
                          ( dfilt6(3)*(Ut(:,x,y,z+3)+Ut(:,x,y,z-3))  &
                           +dfilt6(2)*(Ut(:,x,y,z+2)+Ut(:,x,y,z-2))  &
                           +dfilt6(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                           +dfilt6(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=0
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz4*                        &
                         ( dfilt4(2)*(Ut(:,x,y,z+2)+Ut(:,x,y,z-2))  &
                          +dfilt4(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                          +dfilt4(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=-1
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz2*                        &
                         ( dfilt2(1)*(Ut(:,x,y,z+1)+Ut(:,x,y,z-1))  &
                          +dfilt2(0)*Ut(:,x,y,z) )
      end do
   end do
!$OMP end DO nowait
!
   z=-2
!$OMP DO  !collapse(2)
      do y=yfmin_z,yfmax_z
!$OMP SIMD 
         do x=xfmin_z,xfmax_z
            Un(:,x,y,z) = Ut(:,x,y,z) - cfz1*0.5*(Ut(:,x,y,z)-Ut(:,x,y,z+1))
      end do
   end do
!$OMP end DO nowait
!$OMP BARRIER
!
!
end subroutine filtragez_sup
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine fluxes
  !
  !***** Calcul des flux
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_vectors
  implicit none
  integer :: x,y,z
  !
  !
!$OMP DO  collapse(3)
  do z=-2,nz+3
     do y=-2,ny+3
!!$OMP SIMD 
        do x=-2,nx+3
  E(1,x,y,z) = uo(x,y,z)*Un(1,x,y,z)+Un(2,x,y,z)
  E(2,x,y,z) = uo(x,y,z)*Un(2,x,y,z)+Un(5,x,y,z)
  E(3,x,y,z) = uo(x,y,z)*Un(3,x,y,z)
  E(4,x,y,z) = uo(x,y,z)*Un(4,x,y,z)
  E(5,x,y,z) = uo(x,y,z)*Un(5,x,y,z) + gamma*po(x,y,z)*rhoo(x,y,z)*Un(2,x,y,z)
  !
  F(1,x,y,z) = vo(x,y,z)*Un(1,x,y,z)+Un(3,x,y,z)
  F(2,x,y,z) = vo(x,y,z)*Un(2,x,y,z)
  F(3,x,y,z) = vo(x,y,z)*Un(3,x,y,z)+Un(5,x,y,z)
  F(4,x,y,z) = vo(x,y,z)*Un(4,x,y,z)
  F(5,x,y,z) = vo(x,y,z)*Un(5,x,y,z) + gamma*po(x,y,z)*rhoo(x,y,z)*Un(3,x,y,z)

  G(1,x,y,z) = wo(x,y,z)*Un(1,x,y,z)+Un(4,x,y,z)
  G(2,x,y,z) = wo(x,y,z)*Un(2,x,y,z)
  G(3,x,y,z) = wo(x,y,z)*Un(3,x,y,z)
  G(4,x,y,z) = wo(x,y,z)*Un(4,x,y,z)+Un(5,x,y,z)
  G(5,x,y,z) = wo(x,y,z)*Un(5,x,y,z) + gamma*po(x,y,z)*rhoo(x,y,z)*Un(4,x,y,z)
  !
  H(1,x,y,z) = 0.
  H(2,x,y,z) = duox(x,y,z) * E(1,x,y,z)      &
             + duoy(x,y,z) * F(1,x,y,z)      &
             + duoz(x,y,z) * G(1,x,y,z)
  H(3,x,y,z) = dvox(x,y,z) * E(1,x,y,z)      &
             + dvoy(x,y,z) * F(1,x,y,z)      &
             + dvoz(x,y,z) * G(1,x,y,z)
  H(4,x,y,z) = dwox(x,y,z) * E(1,x,y,z)      &
             + dwoy(x,y,z) * F(1,x,y,z)      &
             + dwoz(x,y,z) * G(1,x,y,z)

  H(5,x,y,z) = (gamma-1.) * ( (duox(x,y,z)+dvoy(x,y,z)+dwoz(x,y,z))*Un(5,x,y,z)       &
                             -(dpox(x,y,z)*Un(2,x,y,z)                    &
                             + dpoy(x,y,z)*Un(3,x,y,z)                    &
                             + dpoz(x,y,z)*Un(4,x,y,z))*rhoo(x,y,z) )

enddo
enddo
enddo
!$OMP end DO nowait
!$OMP BARRIER
end subroutine fluxes
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine ts(irk)
  !
  !***** Termes sources
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_options
  use mod_onde
  use mod_vectors
  use mod_scales
  implicit none
  integer,intent(in) :: irk
  integer :: x,y,z
  !
  !
  !!///// Monopole avec ou sans ecoulement
  if (icas==120) then
!$OMP DO  collapse(2)
     do z=-2,nz+3
        do y=-2,ny+3
!$OMP SIMD 
           do x=-2,nx+3   
              S(1,x,y,z) = amp * sin(omega*(time+ck(irk)*deltat)) &
                               * exp(-alpha*(xg(x)**2 + yg(y)**2 + zg(z)**2))                
              S(5,x,y,z) = S(1,x,y,z)
           end do
        end do
     end do
!$OMP end DO nowait
  end if
  !
  !
end subroutine ts
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine maillage
  !
  !***** Maillage
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_condlim
  use mod_options
  use mod_scales
  use mod_grille
  implicit none
  integer :: i,j,x,y,z
  !
  !
  !//// MAILLAGE PAR DEFAUT
  !
  !
  xg((nx+1)/2)=0
  do i=(nx+1)/2+1,nx+3
     xg(i)=xg(i-1)+deltax
  enddo
  do i=(nx+1)/2-1,-2,-1
     xg(i)=xg(i+1)-deltax
  enddo
  !
  yg((ny+1)/2)=0
  do i=(ny+1)/2+1,ny+3
     yg(i)=yg(i-1)+deltay
  enddo
  do i=(ny+1)/2-1,-2,-1
     yg(i)=yg(i+1)-deltay
  enddo
  !
  zg((nz+1)/2)=0
  do i=(nz+1)/2+1,nz+3
     zg(i)=zg(i-1)+deltaz
  enddo
  do i=(nz+1)/2-1,-2,-1
     zg(i)=zg(i+1)-deltaz
  enddo
  !
  !
  !//// ORIGINE DU RAYONNEMENT ACOUSTIQUE
  !
  !
  nxray = (nx+1)/2
  nyray = (ny+1)/2
  nzray = (nz+1)/2
  !
  !
  !///// DERIVEES DU MAILLAGE
  !
  !
  do x=1,nx
     dxg(x)=0.
     do j=-3,3
        dxg(x) = dxg(x) + a(j)*xg(x+j)
     end do
  end do
  !
  x=0
  dxg(x)=0.
  do j=-2,4
     dxg(x) = dxg(x) + a24(3+j)*xg(x+j)
  end do

  x=-1
  dxg(x)=0.
  do j=-1,5
     dxg(x) = dxg(x) + a15(2+j)*xg(x+j)
  end do

  x=-2
  dxg(x)=0.
  do j=0,6
     dxg(x) = dxg(x) + a06(1+j)*xg(x+j)
  end do
  !
  x=nx+1
  dxg(x)=0.
  do j=-4,2
     dxg(x) = dxg(x) - a24(3-j)*xg(x+j)
  end do

  x=nx+2
  dxg(x)=0.
  do j=-5,1
     dxg(x) = dxg(x) - a15(2-j)*xg(x+j)
  end do

  x=nx+3
  dxg(x)=0.
  do j=-6,0
     dxg(x) = dxg(x) - a06(1-j)*xg(x+j)
  end do
  !
  do y=1,ny
     dyg(y)=0.
     do j=-3,3
        dyg(y) = dyg(y) + a(j)*yg(y+j)
     end do
  end do
  !
  y=0
  dyg(y)=0.
  do j=-2,4
     dyg(y) = dyg(y) + a24(3+j)*yg(y+j)
  end do
  y=-1
  dyg(y)=0.
  do j=-1,5
     dyg(y) = dyg(y) + a15(2+j)*yg(y+j)
  end do
  y=-2
  dyg(y)=0.
  do j=0,6
     dyg(y) = dyg(y) + a06(1+j)*yg(y+j)
  end do
  !
  y=ny+1
  dyg(y)=0.
  do j=-4,2
     dyg(y) = dyg(y) - a24(3-j)*yg(y+j)
  end do
  y=ny+2
  dyg(y)=0.
  do j=-5,1
     dyg(y) = dyg(y) - a15(2-j)*yg(y+j)
  end do
  y=ny+3
  dyg(y)=0.
  do j=-6,0
     dyg(y) = dyg(y) - a06(1-j)*yg(y+j)
  end do
  !
  do z=1,nz
     dzg(z)=0.
     do j=-3,3
        dzg(z) = dzg(z) + a(j)*zg(z+j)
     end do
  end do
  !
  z=0
  dzg(z)=0.
  do j=-2,4
     dzg(z) = dzg(z) + a24(3+j)*zg(z+j)
  end do
  z=-1
  dzg(z)=0.
  do j=-1,5
     dzg(z) = dzg(z) + a15(2+j)*zg(z+j)
  end do
  z=-2
  dzg(z)=0.
  do j=0,6
     dzg(z) = dzg(z) + a06(1+j)*zg(z+j)
  end do
  !
  z=nz+1
  dzg(z)=0.
  do j=-4,2
     dzg(z) = dzg(z) - a24(3-j)*zg(z+j)
  end do
  z=nz+2
  dzg(z)=0.
  do j=-5,1
     dzg(z) = dzg(z) - a15(2-j)*zg(z+j)
  end do
  z=nz+3
  dzg(z)=0.
  do j=-6,0
     dzg(z) = dzg(z) - a06(1-j)*zg(z+j)
  end do
  !
  do x=-2,nx+3
     dxg(x) = 1./dxg(x)
  end do
  !
  do y=-2,ny+3
     dyg(y) = 1./dyg(y)
  end do
  !
  do z=-2,nz+3
     dzg(z) = 1./dzg(z)
  end do
  !
  !
end subroutine maillage
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine record(itime)
  !
  !***** Sauvegarde de l'ecoulement moyen
  !***** Sauvegarde des champs complets a intervalles reguliers
  !***** Sauvegarde de certaines donnees a tous las pas de temps
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_scales
  use mod_onde
  use mod_vectors
  use mod_grille
  use mod_scheme
  use mod_options
  use mod_record
  !
  implicit none
  integer :: itime,x,y,z,i,j
  character (len=60) :: fichier_a_ecrire
  character (len=6) :: char_de_irecord
  character (len=1) :: sep="_"
  character (len=2) :: char_de_irun   
  character (len=4) :: extension=".bin"
  !
  !
  !!///// SAUVEGARDE DE L'ECOULEMENT MOYEN AU DEBUT DU CALCUL

  if (mod(itime,record_step).eq.0) then
     call dm_num2char_6(char_de_irecord,irecord)
     if(o_record_vort) then
        call calculvort
     endif
  endif

!$OMP SECTIONS
!$OMP SECTION
  if (itime.eq.0) then
     fichier_a_ecrire=trim(record_filename_ecoulmoy)//extension
     write(6,*) 'Enregistrement champ moyen:',  '     t/T = ',time/(2.*pi/omega)
     open(501,file=fichier_a_ecrire,form='unformatted',status='unknown')
     write(501) record_dummy
     write(501) (((1./rhoo(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((uo(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((vo(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((wo(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((po(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) record_dummy
     write(501) (((duox(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((duoy(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((duoz(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) record_dummy
     write(501) (((dvox(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dvoy(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dvoz(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) record_dummy
     write(501) (((dwox(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dwoy(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dwoz(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) record_dummy
     write(501) (((dpox(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dpoy(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) (((dpoz(x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(501) record_dummy
     close(501)
  endif
  !
  !
  !!///// SAUVEGARDE DES CHAMPS DE PRESSION TOUS LES record_step ITERATIONS
!$OMP SECTION
  if (mod(itime,record_step).eq.0) then
     fichier_a_ecrire=trim(record_filename_champs)//sep//char_de_irecord//extension
     write(6,*) 'Enregistrement champs complets:',irecord,'/',1+int((ntfin-nt0)/record_step),  '     t/T = ',time/(2.*pi/omega)
     open(502,file=fichier_a_ecrire,form='unformatted',status='unknown')
     write(502) record_dummy
     write(502) itime
     write(502) record_dummy
     write(502) (((Un(1,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(502) record_dummy
     write(502) (((Un(2,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(502) record_dummy
     write(502) (((Un(3,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(502) record_dummy
     write(502) (((Un(4,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(502) record_dummy
     write(502) (((Un(5,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
     write(502) record_dummy       
     close(502)
  endif
     !
     !!// Enregistrement vorticite
!$OMP SECTION
  if (mod(itime,record_step).eq.0) then
     if(o_record_vort) then
        fichier_a_ecrire=trim(record_filename_vort)//sep//char_de_irecord//extension
        write(6,*) '   + Enregistrement vorticite'
        open(503,file=fichier_a_ecrire,form='unformatted',status='unknown')
        write(503) record_dummy
        write(503) itime
        write(503) record_dummy
        write(503) ((((VORT(i,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3),i=1,3)
        write(503) record_dummy
        close(503)
     end if
     !
  end if

  !
  !
  !!///// SAUVEGARDE DU RESIDU
  !!      do x=1,nx
  !         do y=1,ny
  !            residu(itime) = residu(itime) + U(x,y,4)**2
  !         end do
  !      end do
  !      residu(itime) = sqrt(residu(itime)/real(nx*ny))
  !      if(itime.eq.ntfin) then
  !         fichier_a_ecrire=trim("bin_residu_p")//extension
  !         write(6,*) 'Enregistrement du residu de pression.' 
  !        open(504,file=fichier_a_ecrire,form='unformatted',status='unknown')
  !         write(504) record_dummy
  !         write(504) (residu(i),i=0,ntfin)
  !         write(504) record_dummy
  !         close(504)
  !   endif
  !
  !!///// SAUVEGARDE A TOUS LES PAS DE TEMPS
  !
  !
!$OMP SECTION
  if (o_rect) then
     if(itime.eq.ntfin) then
       do z=1,nzrect
          do y=1,nyrect
             do x=1,nxrect
                rect(itime,x,y,z,:)= U(:,xrect(x),yrect(y),zrect(z))
             end do
          end do
       end do
        call dm_num2char_2(char_de_irun,irun)
        fichier_a_ecrire=trim("bin_rect")//sep//trim("run")//char_de_irun//extension
        write(6,*) 'Enregistrement du fichier rect a tous les pas de temps.' 
        open(505,file=fichier_a_ecrire,form='unformatted',status='unknown')
        write(505) record_dummy
        write(505) nxrect
        write(505) nyrect
        write(505) nzrect
        write(505) record_dummy
        write(505) xrect
        write(505) yrect
        write(505) zrect
        write(505) record_dummy
        write(505) (((((rect(i,x,y,z,j),i=nt0,ntfin),x=1,nxrect),y=1,nyrect),z=1,nzrect),j=1,5)
        write(505) record_dummy
        close(505)
     endif
  end if
  !
  !
!$OMP END SECTIONS nowait
  if (mod(itime,record_step).eq.0) irecord=irecord+1
end subroutine record
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine record_for_restart
  !
  !***** Sauvegarde les informations a la fin d'un run
  !***** pour pouvoir faire un restart
  !******************************************************************
  !******************************************************************
  use mod_options
  use mod_scales
  use mod_onde
  use mod_vectors
  use mod_grille
  use mod_scheme
  use mod_options
  use mod_record
  implicit none
  character (len=2) :: char_de_irun
  character (len=60) :: fichier_a_ecrire
  character (len=1) :: sep="_"
  character (len=4) :: extension=".bin"
  integer :: x,y,z
  !
  !
  call dm_num2char_2(char_de_irun,irun)
  fichier_a_ecrire=trim("bin_restart")//sep//trim("run")//char_de_irun//extension
  write(6,*) 'Enregistrement du fichier restart (fin du run ', char_de_irun,').' 
  open(506,file=fichier_a_ecrire,form='unformatted',status='unknown')
  write(506) record_dummy
  write(506) nt0 
  write(506) ntfin
  write(506) irun
  write(506) irecord
  write(506) time
  write(506) record_dummy
  write(506) (((Un(1,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
  write(506) (((Un(2,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
  write(506) (((Un(3,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
  write(506) (((Un(4,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
  write(506) (((Un(5,x,y,z),x=-2,nx+3),y=-2,ny+3),z=-2,nz+3)
  write(506) record_dummy
  close(506)
  !
end subroutine record_for_restart
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine sauveparametre
  !
  !***** Sauvegarde des parametres essentiels de la simulation
  !******************************************************************
  !******************************************************************
  use mod_grille
  use mod_filtrage
  use mod_onde
  use mod_options
  use mod_record
  use mod_condlim
  implicit none
  integer :: int_option
  !
  open(500,file='bin_parametre.bin',form='unformatted',status='unknown')
  write(500) record_dummy
  write(500) deltax
  write(500) deltay
  write(500) deltaz
  write(500) deltat
  write(500) record_dummy
  write(500) nx
  write(500) ny
  write(500) nz
  write(500) nymin
  write(500) nymax
  write(500) nt0
  write(500) ntfin
  write(500) record_step
  write(500) record_dummy
  write(500) xg
  write(500) yg
  write(500) zg
  write(500) record_dummy
  write(500) mo
  write(500) amp
  write(500) omega 
  write(500) alpha
  write(500) record_dummy
  if(o_damping) then
     int_option=1
  else
     int_option=0
  end if
  write(500) int_option
  write(500) record_dummy
  if(o_damping_sup) then
     int_option=1
  else
     int_option=0
  end if
  write(500) int_option
  write(500) record_dummy
  write(500) cfx
  write(500) cfy
  write(500) cfz
  write(500) icas
  write(500) irun
  write(500) record_dummy
  close(500)
  !
  !
end subroutine sauveparametre
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine dm_num2char_6(char,num)
  !
  !***** Transforme un entier num avec 6 digits en un caractere char
  !***** ex: num=5 devient char="000005"
  !******************************************************************
  !******************************************************************
  implicit none
  integer :: num
  character(len=6) :: char
  !
  if (num.ge.100000) then
     write (char,1) num
1    format (i6)
  else
     if (num.ge.10000) then
        write (char,2) num
2       format ('0',i5)
     else
        if(num.ge.1000) then
           write (char,3) num
3          format ('00',i4)
        else
           if (num.ge.100) then
              write (char,4) num
4             format ('000',i3)
           else
              if (num.ge.10) then
                 write (char,5) num
5                format ('0000',i2)
              else
                 write (char,6) num
6                format ('00000',i1)
              endif
           endif
        endif
     endif
  endif
  !
end subroutine dm_num2char_6
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine dm_num2char_2(char,num)
  !
  !***** Transforme un entier num avec 2 digits en un caractere char
  !***** ex: num=5 devient char="05"
  !******************************************************************
  !******************************************************************
  implicit none
  integer :: num
  character(len=2) :: char
  !
  if (num>=10) then
     write (char,1) num
1    format (i2)
  else
     write (char,2) num
2    format ('0',i1)
  endif
  !
end subroutine dm_num2char_2
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine calculvort
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_vectors
  implicit none
  integer :: x,y,z,j
  !
  !
  VORT=0.
!$OMP DO  collapse(2)
  do z=1,nz
     do y=1,ny
!$OMP SIMD 
        do x=1,nx
           do j=-3,3
              VORT(1,x,y,z) = VORT(1,x,y,z) + dyg(y)*a(j)*U(4,x,y+j,z)          &
                                            - dzg(z)*a(j)*U(3,x,y,z+j)
              VORT(2,x,y,z) = VORT(2,x,y,z) + dzg(z)*a(j)*U(2,x,y,z+j)          &
                                            - dxg(x)*a(j)*U(4,x+j,y,z)
              VORT(3,x,y,z) = VORT(3,x,y,z) + dxg(x)*a(j)*U(3,x+j,y,z)          &
                                            - dyg(y)*a(j)*U(2,x,y+j,z)
           end do
        end do
     end do
  end do
!$OMP end DO nowait
!$OMP BARRIER
  !
  !
end subroutine calculvort
!******************************************************************
!
!
!
!
!******************************************************************
!******************************************************************
subroutine coeff_schemas
  !
  !***** Coefficients des schemas DRP, des filtres, de Runge-Kutta
  !***** Refs: - Tam & Shen, AIAA Paper 93-4325
  !*****       - Tam, AIAA J 33, 1995, 1788-1796
  !*****       - Bogey et Bailly, JCP, 2004, 194-214
  !******************************************************************
  !******************************************************************
  use mod_scheme
  use mod_options
  implicit none
  !
  !///// COEFF DRP (Tam&Shen,1993)
  !
  a(1) =  0.770882380518
  a(2) = -0.166705904415
  a(3) =  0.020843142770
  a(0) =  0.
  a(-1) = -a(1)
  a(-2) = -a(2)
  a(-3) = -a(3)
  !
  !///// COEFF DRP DECENTRES (Tam,1995)
  !
  a24(1) =  0.049041958
  a24(2) = -0.468840357
  a24(3) = -0.474760914
  a24(4) =  1.273274737
  a24(5) = -0.518484526
  a24(6) =  0.166138533
  a24(7) = -0.026369431
  !
  a15(1) = -0.209337622
  a15(2) = -1.084875676
  a15(3) =  2.147776050
  a15(4) = -1.388928322
  a15(5) =  0.768949766
  a15(6) = -0.281814650
  a15(7) =  0.048230454
  !
  a06(1) = -2.192280339
  a06(2) =  4.748611401
  a06(3) = -5.108851915
  a06(4) =  4.461567104
  a06(5) = -2.833498741
  a06(6) =  1.128328861
  a06(7) = -0.203876371
  !
  !.. DFC ordre 6 / attention a l'ordre des coef. a(:) = dfc6(:)
  !
  dfc6(1) =  3./4.
  dfc6(2) = -3./20.
  dfc6(3) =  1./60.
  dfc6(0) =  0.
  dfc6(-1) = -dfc6(1)
  dfc6(-2) = -dfc6(2)
  dfc6(-3) = -dfc6(3)
  !
  !.. DFC ordre 4
  !
  dfc4(3) = 0.
  dfc4(4) = 2./3.
  dfc4(5) = -1./12.
  dfc4(2) = -dfc4(4)
  dfc4(1) = -dfc4(5)
  !
  dfc13(1) = -3./12.
  dfc13(2) = -10./12.
  dfc13(3) = 18./12.
  dfc13(4) = -6./12.
  dfc13(5) = 1./12.
  !
  dfc31(1:5) = -dfc13(5:1:-1)
  !
  dfc04(1) = -25./12.
  dfc04(2) = 48./12.
  dfc04(3) = -36./12.
  dfc04(4) = 16./12.
  dfc04(5) = -3./12.
  !
  dfc40(1:5) = -dfc04(5:1:-1)
  !
  !.. DFC ordre 2
  !
  dfc2(2) = 0.
  dfc2(3) = 1./2.
  dfc2(1) = -dfc2(3)
  !
  dfc02(1) = -3./2.
  dfc02(2) = 2.
  dfc02(3) = -1./2.
  !
  dfc20(1) = 1./2.
  dfc20(2) = -2.
  dfc20(3) = 3./2.
  !
  !///// FILTRAGE ORDRE 8 (Bogey et Bailly, 2004)
  !
  dfilt8(0) = 35./128.
  dfilt8(1) = -7./32.
  dfilt8(2) = 7./64.
  dfilt8(3) = -1./32.
  dfilt8(4) = 1./256.
  dfilt8(-1) = dfilt8(1)
  dfilt8(-2) = dfilt8(2)
  dfilt8(-3) = dfilt8(3)
  dfilt8(-4) = dfilt8(4)
  !
  !///// FILTRAGE ORDRE 6 STANDARD
  !
  dfilt6(0) =  5./16.
  dfilt6(1) = -15./64.
  dfilt6(2) =  3./32.
  dfilt6(3) = -1./64.
  dfilt6(-1) = dfilt6(1)
  dfilt6(-2) = dfilt6(2)
  dfilt6(-3) = dfilt6(3)
  ! 
  !///// FILTRAGE ORDRE 4 STANDARD (Tam,1995)
  ! 
  dfilt4(0) = 0.375
  dfilt4(1) = -0.25
  dfilt4(2) = 0.0625
  dfilt4(-1) = dfilt4(1)
  dfilt4(-2) = dfilt4(2)
  !
  !///// FILTRAGE ORDRE 2 STANDARD (Tam,1995)
  !
  dfilt2(0) = 0.5
  dfilt2(1) = -0.25
  dfilt2(-1) = dfilt2(1)
  !
  !
  !///// REDUCTION DE L'ORDRE DUR LES BORDS
  !
  !
  if(o_reduction_order_bords) then
     write(6,*) '*****    REDUCTION ORDRE SUR LES BORDS'
     write(6,*) '*****    --> ATTENTION : coeffs a15 et a24 MODIFIES!!!'
     !
!!!a06=0
!!!a06(1:3)=dfc02(1:3)
!!!a06(1:5)=dfc04(1:5)
     !
     a15=0
     a15(1:3)=dfc2(1:3)
     !
     a24=0
     a24(1:5)=dfc4(1:5)
     !
  end if
  !
  !
  !// SCHEMAS DE RUNGE-KUTTA
  !
  !    - numrk : choix du schema de Runge-Kutta
  !            = 1 : RK ordre 4 (lineaire)
  !            = 2 : RK ordre 2 Hu et al.
  !            = 3 : RK ordre 2 Bogey
  rk = 0.
  ck = 0.
  if (nrk.eq.4) then
     rk(1) = 1./4.
     rk(2) = 1./3.
     rk(3) = 1./2.
     rk(4) = 1.
     !
     ck(1) = 0.
     ck(2) = 1./4.
     ck(3) = 1./3.
     ck(4) = 1./2.
  elseif (nrk.eq.5) then
     rk(1) = 0.169193539
     rk(2) = 0.23717924
     rk(3) = 0.333116
     rk(4) = 1./2. 
     rk(5) = 1.
     !
     ck(1) = 0.
     ck(2) = 0.169193539
     ck(3) = 0.23717924
     ck(4) = 0.333116
     ck(5) = 1./2.
  end if
  !
  !
end subroutine coeff_schemas
!******************************************************************


