c*********** No-Core Slater Determinant shell model code (NCSD) with hashing search
c*********** P. Navratil
c*********** Version with linear hash table and three shift key generation
c*********** Basis generation for protons and neutrons separate
c*********** Global configurations not constructed
c*********** memsave option to avoid allocation of the third real(8) nsd-size array
c*********** BLAS library function calls included
c*********** TRIUMF, 11/02/2011
c*********** Compile by: mpiifort -O3  -lblas
c*********** Input and output the same as in the MFD Arizona version
c*********** by J. Vary, D.C. Zheng, P. Navratil                 ***********
c*********** Name of the input file: mfdp.dat
c*********** Two-body interaction file name given in mfdp.dat
c ********** In order to run MPI, comment out or remove dummy subroutines **
c ********** at the end of this file ***************************************
c ********** In a sequential run, the dummy subroutines must be present ****
c ********** In addition, a phony mpif.h file needs to be supplied or ******
c ********** the statements " include 'mpif.h' " should be commented out *** 
c
c ********** Number of processors: Any *************************************
c
c ********** Restart options:                                   ************
c ********** irest=0 ... standard calculation                   ************
c ********** irest=4 ... generate more Lanczos iterations       ************
c ********** irest=6 ... IT iloc available, start Lanczos iterations *******      
c ********** All restart options require the same number of     ************
c ********** processors                                         ************
c ********** Three-body interaction read from the file          ************
c ********** named in mfdp.dat                                  ************
c ********** In mfdp.dat set:                                   ************
c ********** nbit=3  ... Three-body interaction run             ************
c ********** nbit/= 3 ...standard run from TBME file            ************
c
c----------------------------------------------------------------------------
      module parameters
c     ----- NOTE: nbit, nbit1 have to be changed according to the
c     -----            computer used: 
c**************************************************************************
c******************** Current version:
c     -----            nbit=64 for a 64-bit machine
c**************************************************************************
      integer nbit,nbit1
      parameter (nbit=64,nbit1=nbit-1)
c------------------------------------------------------------------------

c     ***** ONE-MAJOR-SHELL SPACE *****
c     ----- For the 0p shell: jmx=3, nobt=2, ntbme=15
c     ----- For the sd shell: jmx=5, nobt=3, ntbme=63
c     ----- For the fp shell: jmx=7, nobt=4, ntbme=195
c     ----- mxnsd   = Maximum number of Slater Determinants allowed.
c     ----- mxnwd   = Maximum number of words allowed. This defines the 
c     -----           Maximum number of s.p. states (mxnwd*nbit) allowed for
c     -----           each class of nucleons.
c     -----           nwd=1 for 32 or less s.p. states;
c     -----           nwd=2 for more than 32 but less than 65 s.p. states.
c     ----- jmx     = maximum value of J for all the two-body states allowed.
c     ----- mxnobt  = Maximum number of s.p. orbitals in the model space.

c------------------------------------------------------------------------

c     ***** MULTI-MAJOR-SHELL SPACE *****
c     ----- The N=3 & N1+N2=6 space: nwd= 3, ntbme= 2308, jmx=5, nobt=10
c     ----- The N=4 & N1+N2=4 space: nwd= 3, ntbme=  994, jmx=5, nobt=15
c     ----- The N=5 & N1+N2=5 space: nwd= 4, ntbme= 3018, jmx=6, nobt=21
c     ----- The N=5 & N1+N2=6 space: nwd= 4, ntbme= 7557, jmx=7, nobt=21
c     ----- The N=6 & N1+N2=6 space: nwd= 6, ntbme= 8260, jmx=7, nobt=28
c     ----- The N=6 & N1+N2=7 space: nwd= 6, ntbme=19170, jmx=8, nobt=28
c     ----- The N=7 & N1+N2=7 space: nwd= 8, ntbme=20340, jmx=9, nobt=36
c     ----- The N=7 & N1+N2=8 space: nwd= 8, ntbme=44556, jmx=9, nobt=36
c     ----- The N=8 & N1+N2=8 space: nwd=11, ntbme=46408, jmx=9, nobt=45
c***** nwd/2(+1) is enough when nbit=64
      integer jmx
c------------------------------------------------------------------------

c     ----- PERMANENT: No change needed for the following line unless there
c     ----- are more than two classes of particles
c     ----- mxsps = nbit*nwd, mxsps2 = 2*mxsps
      integer mxsps,mxsps2

      end module parameters
c------------------------------------------------------------------------

      module config
c     mnpi(N): minimum number of protons  in the N-th shell
c     mnnu(N): minimum number of neutrons in the N-th shell
c     mxpi(N): maximum number of protons  in the N-th shell
c     mxnu(N): maximum number of neutrons in the N-th shell
      integer,allocatable,dimension(:) :: mnpi,mxpi,mnnu,mxnu
     +     ,mnpn,mxpn,ntempi,ntemnu
c     nsps0(N): number of s.p. states in the N-th shell
c     nspsb(N): number of s.p. states seen before the N-th shell
c     npib(N):  number of protons  seen before the N-th shell
c     nnub(N):  number of neutrons seen before the N-th shell
c     nops(N):  number of proton  many-body states in the N-th shell
c     nons(N):  number of neutron many-body states in the N-th shell
      integer,allocatable,dimension(:) :: nsps0,nspsb,npib,nnub,
     +     nops,nons,neng
c     nprtn(iconf,N): number of protons  in the major shell N in the iconf-th
c                     configuration
c     nneut(iconf,N): number of neutrons in the major shell N in the iconf-th
c                     configuration
      integer,allocatable,dimension(:,:) :: nprtn,nneut 
      end module config

      module ibas
c**      integer(4),allocatable::ibasis(:,:,:)
      integer(8),allocatable::ibasis(:,:,:)
c      integer(2),allocatable::iloc(:,:)
      integer(2),allocatable,target::ilocp(:,:),ilocn(:,:)
      integer(4),allocatable,target::iloc(:,:)
      integer :: nsdp,nsdn,nminho
      integer(4),allocatable:: nsd_hw(:)
      real(kind(0.0)),allocatable :: delta_IT(:)
      real(kind(0.d0)) :: deltaIT
      end module ibas

      module nuconf
      integer(2),allocatable :: negy(:) !,mconf(:),ndiff(:,:)
      integer(4),allocatable :: iendconf(:)
      integer(4) nconf 
      end module nuconf

      module lanczvec
      real(8),allocatable :: amp(:,:)
      end module lanczvec

      module nodeinfo
      integer nproc,iproc,icomm
      end module nodeinfo

      module multig
      integer nsetm1,nespmin,iset1,ippint,innint
      end module multig
      
      module ctrl
      integer ki,kf,nf
      real egs
      end module ctrl

      module jsgna
      integer,parameter :: jsgn(-127:128)=(/(-1,1,i=-127,127,2)/)
      end module jsgna

      module spodata
      integer,allocatable :: nnl(:),lx(:),j2x(:)
      end module spodata

      module nutbme
      integer ntbme
      end module nutbme

      module idxb
      integer idxobt(66)
      integer, allocatable :: iorder(:)
      end module idxb

      module consts
      real sqrt2
      real xxxxx
      end module consts

      module iters
      integer :: nsd,nite,major,ibuf
c     MKGK
      integer :: iter_old
      integer :: memsave
c      integer(8) :: nhme
      end module iters

      module spb
      integer :: nucleons,nprotons,nneutrns,mnop,igt,irest=0
      integer :: mjtotal,mttotal,nshll,nwd,nasps,nhw,iparity,nhwp,nhwn
      integer :: mnn2l,mxn2l
      real hbomeg
c     variables added by mk
      integer(2) :: nhw_boot, nhw0, nhw_min, nhw_max, nhw_restart=-1
      integer(2) :: nhw_now,nhw_hold 
      end module spb

      module v3b
c      integer(8),allocatable :: pointmtpi(:,:,:)
c      integer(4),allocatable :: index3(:),timerevsp(:)
c      integer(8) :: ih3ind
c      real(kind=kind(0.0)),allocatable :: H3(:)
      real(kind(0.0)),allocatable :: v3b_cJ(:)   ! must be defined in v3b module
      integer,allocatable :: index_abc(:,:,:),index_abcdef(:,:),
     $     start_abcdef(:),nlj_orb(:) ! must be defined in v3b module
      real(kind=kind(0.0)),allocatable :: cgj12(:,:,:,:,:),
     $     cgj123(:,:,:,:,:)
      real(kind=kind(0.0)) :: cgt12(0:1,0:1,0:1),
     $     cgt123(0:1,0:1,-1:1,0:1)
      integer :: E1max,E12max,E123max
      end module v3b

      module paramdef
      type jsp
      integer,allocatable :: j(:)
      end type jsp
      type lspjsp
      type(jsp),allocatable :: l(:)
      end type lspjsp
      type(lspjsp),allocatable :: nlj_st_TUD(:)
      end module paramdef

      module TUD_tbme
      integer :: l_Max,n_Max,ist1dim_TUD
      real(kind(0.d0)),allocatable :: V_NN_TUD(:),
     $     Trel_TUD(:),HOrel_TUD(:),V_3N_no1b_TUD(:,:)
      real(kind(0.d0)) :: V_3N_no0b
      real(kind(0.0)),allocatable :: H_NN_TUD(:)
      integer :: tot_num_of_NNmatel_TUD
      integer :: N1_max_tud,N12_max_tud
      integer,allocatable :: ist1_TUD(:,:),
     $     index_ab(:,:),index_abcd(:,:)
      integer :: dim_ab_TUD,dim_abcd_TUD
      logical :: Trelsave=.true.,tbmeTUD=.false.,no2bv3n=.false.
      real(kind=kind(0.0)) :: cgt12_tud(0:1,-1:1,0:1,0:1)
      character(len=120) :: v3intfile
      end module TUD_tbme

      module pn_tbme
      integer :: N1_max_pn,N12_max_pn
      integer :: num_of_2bst_pn(-1:1),num_of_2bel_pn(-1:1)
      integer,allocatable :: isp2npoi_pn(:,:,:),
     $     isp2ndim_pn(:,:,:),i2belnpoi_pn(:,:,:)
      logical :: tbmepn=.false.
      type ist_pn
      integer,allocatable :: ist(:,:)
      end type ist_pn
      type(ist_pn) :: isp2n_pn(-1:1)
      type el_pn
      real(kind(0.d0)),allocatable :: el(:)
      end type el_pn      
      type(el_pn) :: V_NN_pn(-1:1),Trel_pn(-1:1),HOrel_pn(-1:1),
     $     H_NN_pn(-1:1),rful_pn(-1:1)
      type sp_2
      integer :: sp2min,sp2max
      integer,allocatable :: sp2(:)
      end type sp_2
      type sp_1
      integer :: sp1min,sp1max
      type(sp_2),allocatable :: sp1(:)
      end type sp_1
      type Tz_tbd
      type(sp_1) :: Tz(-1:1)
      end type Tz_tbd
      type(Tz_tbd),allocatable :: pntbdst(:,:)
      end module pn_tbme

      module spbasis
c      type J3T3_block
c      integer :: begin,end,dim
c      end type J3T3_block
c      type(J3T3_block),allocatable :: block_J3T3(:)
c      integer :: block_J3T3_dim
      integer :: N1_max,isp1ntot,isp3ntot,dim_abc,dim_abcdef
      integer,allocatable :: isp1n_cJ(:,:),isp3n_cJ(:,:),sp3ntot(:,:,:)
      end module spbasis

      module spsdata
      integer,allocatable :: n2l_sp(:),n_sp(:),l_sp(:),
     +     j2_sp(:),m2_sp(:),mt2_sp(:),nobt_sp(:)
      integer :: N12_max
      end module spsdata

      module effop
      real epi,enu,glp,gln,gsp,gsn
      end module effop

      module tpsb
      integer,allocatable :: nmebjtp(:,:,:),ntps(:,:,:),
     +     itps(:,:,:,:) 
      end module tpsb

      module tbmes
      real,allocatable :: gful(:,:),cful(:),rful(:),Esp(:,:)
      real :: E_core
      end module tbmes

      module tbmepar
      integer :: intform,nobt,iclmb,itrel,ihrel,ispe,nskip
      real(kind(0.0)) :: strcm
      character(len=120) :: intfile
      end module tbmepar

      module facs
      real(8) fac(0:64),sqfac(0:64) ! fac(k) = k!, sqfac(k) = sqrt(k!)
      end module facs

      module the3js
      real the3t(0:1,-1:1,-1:1)
      real,allocatable :: the3j(:,:,:)
      end module the3js

      module bits
      integer,allocatable :: locan(:),locn(:),locz(:)
      end module bits

      module tempor
      integer nsdx,nmfss2x
      end module tempor

      module albe
      real(8),allocatable :: alpha(:),beta(:),energy(:)
c     MKGK
      real(8), allocatable :: energy_con(:)
      real(8) :: conset
      integer :: conflag
      end module albe

      module vect
      real,allocatable :: v(:,:)
      end module vect

c      module hamc
c      integer ii,jj,kk,ll
c      end module hamc

c      module i1i3
c      integer i1,i3
c      end module i1i3

      module hmes
      real,allocatable :: hsaved(:)
      integer,allocatable :: i3saved(:)
      integer :: i1maxsaved=0,bufferdimi1saved=1000000
      integer(8) :: bufferhsaved=200000000
      integer(8),allocatable :: dimi1saved(:)
      integer(8) :: memoryalloc=0,memoryallocva=0,memavail=0
      end module hmes

      module nnmaxmod
      integer nnmax
      end module nnmaxmod

      module obmes
      real,allocatable :: xlp(:,:),xln(:,:),
     +     xsp(:,:),xsn(:,:),
     +     xqp(:,:),xqn(:,:)
      end module obmes

      module r2HOme
      real(8),allocatable :: r2mea(:,:,:,:)
      end module r2HOme

      module gamalog
      integer nrmax,lrelm,maxgam
      double precision,allocatable,dimension(:,:) :: dsq
      double precision,allocatable,dimension(:) :: gamal
      end module gamalog 

c     MKGK added modules
      module ITvar
      real(8) :: kappa, cmin
      logical :: piv_saved = .false.
      integer,allocatable :: nsd_kappa(:)
      integer :: kappa_counter,kappa_points,kappa_restart=-1
      real(8),allocatable :: kappa_store(:)
      end module ITvar
c     end of MKGK added modules

c     MKGK
c     module occ and hash_tables is now in a seperate file
c     named bhashed.f

c     module is part of the NCSD code
      module occ
      type occ_bas
      integer,pointer :: nucl,protn,neutn,dim_st,dim_st_p,dim_st_n
      integer(4),pointer :: occ(:,:)
      integer(2),pointer :: occp(:,:),occn(:,:)
      end type occ_bas
      type(occ_bas) :: ty
      end module occ

      module hash_tables
      use occ
      private
      public construct_hash_table,get_state_index,init_hash_table_only,
     $     insert_hash_table_entry,insert_hash_table_entry_word,
     $     dim_constr,occ_constr,max_el_constr,
     $     insert_hash_table_entry_N,
c* PN
     $     dim_constr_p, dim_constr_n,dim_aloc_p,dim_aloc_n  !,
c     $     get_state_index_loc
      integer :: dim,maxdisp,dimp,dimn,maxdispp,maxdispn
c     added from eigv_Jplus.f (one line below)
      integer :: mxnw,dim_constr,dim_aloc,max_el_constr
c* PN
      integer :: dim_constr_p, dim_constr_n,dim_aloc_p,dim_aloc_n
      integer,allocatable :: linhash(:),linhashp(:),linhashn(:)
      integer(2),allocatable :: occ_tmp(:)
      integer(8),allocatable :: intbas_tmp(:)
      integer(8),allocatable :: intbas(:),intbasN(:)
      integer(4),pointer :: occ_constr(:,:)

      contains

      integer function pos_in_iloc(occ,nucx,typ)
      use ibas, only : nsdp, nsdn, ilocp, ilocn
      implicit none
      integer, intent(IN) :: typ, nucx
      integer(2), intent(IN) :: occ(nucx)
      integer :: sd, pos, location

      location = -1
c* PN
c     search a proton iloc
      if (typ==1) then
         do sd=1,dim_constr_p !nsdp
            do pos=1,nucx
               if (occ(pos).ne.ilocp(pos,sd)) exit !cycle
               if (pos==nucx) then 
                  location = sd
                  exit
               end if
            end do
         end do
c     search the neutrons
      else if (typ==2) then
         do sd=1,dim_constr_n !nsdn
            do pos=1,nucx
               if (occ(pos).ne.ilocn(pos,sd)) exit !cycle
               if (pos==nucx) then
                  location = sd
                  exit
               end if
            end do
         end do
      end if
c      if (location==-1) print*,'location not found in pos_in_iloc'
      pos_in_iloc=location
      end function pos_in_iloc

      subroutine construct_hash_table(dimfac,type)
      use spb, only: nwd
      use nodeinfo
      use hmes, only: memoryallocva
      implicit none
c* PN
      integer :: dimfac
      integer,intent(IN) :: type
      integer,allocatable :: temp_ind(:)
      integer :: i,key,dim_st,nuc,max_el,temp,wd,cr_1,ibit,ip,in,
     $     nbit,dim_st_p,dim_st_n
c**      integer(8),allocatable :: intbas(:)

      nbit=64
      dim_st=ty%dim_st
      nuc=ty%nucl

c     MKGK
      dim_constr = dim_st
      dim_aloc = dim_st
      dim_aloc_p=ty%dim_st_p
      dim_aloc_n=ty%dim_st_n


c      if (dim_st<1 431 655 765) then
c* PN
      if (dim_st*dimfac<750 000 000) then
c         dim=dim_st*3/2
         i=nint(log(real(dim_st*dimfac,
     $        kind(0.d0)))/log(2.d0))
         if (i<=29) then
            dim=2**(i+1)
         else
            dim=2100000000
         endif
c         dim=2*(2**(nint(log(real(dim_st*dimfac,
c     $        kind(0.d0)))/log(2.d0))))
c         if (dim<0) dim=2 147 483 648
      else
c         dim=dim_st
         dim=2100000000
cc         dim=2147483648
      endif

      if (iproc==0) print *,' dim,hash table length=',dim_st,dim

c     MKGK
      memoryallocva=memoryallocva+dim*4
      if (allocated(linhash)) deallocate(linhash)
      allocate(linhash(0:dim-1))
      linhash=-1
      maxdisp=0

c     MKGK
      if (allocated(intbas)) deallocate(intbas)
      if (allocated(occ_tmp)) deallocate(occ_tmp)
      allocate(intbas(2*nwd))
      allocate(occ_tmp(nuc))

      do i=1,dim_st

         if (i==10000000*(i/10000000)) then
            if (iproc==0) print *,' construct_hash_table: i=',i
         endif

         ip=ty%occ(1,i)
         in=ty%occ(2,i)
         occ_tmp(1:ty%protn)=ty%occp(:,ip)
         occ_tmp(ty%protn+1:nuc)=ty%occn(:,in)
c         nwd=(maxval(occ_tmp)-1)/nbit+1
c         if (nwd>mxnwd) mxnwd=nwd
         intbas=0
         do cr_1=1,ty%nucl
            wd=(occ_tmp(cr_1)-1)/nbit+1
            ibit=mod(occ_tmp(cr_1)-1,nbit)
            intbas(wd)=ibset(intbas(wd),ibit)
         end do
c         print *,' occ_tmp=',occ_tmp

c         call get_key(nuc,occ_tmp,key,dim)
         call get_key(2*nwd,intbas,key,dim)

         if (key<0.or.key>dim-1) then
            print *,'*** error in construct_hash_table: key,dim=',
     +           key,dim
            stop
         endif
c         deallocate(intbas)
c         print *,' key=',key

         max_el=0
         do
            if (linhash(key)==-1) then
               linhash(key)=i
               exit
            endif
            key=mod(key+1,dim)
            max_el=max_el+1
         end do
         if (max_el>maxdisp) maxdisp=max_el
      end do
c*PN
cc      deallocate(intbas)
c     MKGK
      if (allocated(intbas_tmp)) deallocate(intbas_tmp)
      allocate(intbas_tmp(2*nwd))
      if (iproc==0) print *,' maximal displacement on proc0=',maxdisp

      if (type>0) then

         mxnw=nwd
         if (allocated(intbasN)) deallocate(intbasN)
         allocate(intbasN(2*mxnw))

         dimfac=max(3,dimfac/4)

         dim_st_p=ty%dim_st_p
         nuc=ty%protn
         dim_constr_p = dim_st_p

         if (dim_st_p*dimfac<750 000 000) then
c         dim=dim_st*3/2
            dimp=2*(2**(nint(log(real(dim_st_p*dimfac,
     $           kind(0.d0)))/log(2.d0))))
         else
c         dim=dim_st
            dimp=2147483648
         endif

         if (iproc==0) print *,' ilocp: dim,hash table length=',
     $        dim_st_p,dimp

         memoryallocva=memoryallocva+dimp*4
         if (allocated(linhashp)) deallocate(linhashp)
         allocate(linhashp(0:dimp-1))
         linhashp=-1
         maxdispp=0

c         if (allocated(occ_tmpN)) deallocate(occ_tmpN)
c         allocate(occ_tmpN(nuc))

         do i=1,dim_st_p

            if (i==1000000*(i/1000000)) then
               if (iproc==0) print *,' construct_hash_table: i=',i
            endif

            occ_tmp(1:nuc)=ty%occp(:,i)

            intbasN=0
            do cr_1=1,nuc
               wd=(occ_tmp(cr_1)-1)/nbit+1
               ibit=mod(occ_tmp(cr_1)-1,nbit)
               intbasN(wd)=ibset(intbasN(wd),ibit)
            end do

            call get_key(mxnw,intbasN(1:mxnw),key,dimp)
            if (key<0.or.key>dimp-1) then
               print *,'*** error in construct_hash_table: key,dimp=',
     +              key,dimp
               stop
            endif

            max_el=0
            do
               if (linhashp(key)==-1) then
                  linhashp(key)=i
                  exit
               endif
               key=mod(key+1,dimp)
               max_el=max_el+1
            end do
            if (max_el>maxdispp) maxdispp=max_el
         end do

         if (iproc==0) print *,' ilocp: maximal displacement on proc0=',
     $        maxdispp

         dim_st_n=ty%dim_st_n
         nuc=ty%neutn
         dim_constr_n = dim_st_n

         if (dim_st_n*dimfac<750 000 000) then
c         dim=dim_st*3/2
            dimn=2*(2**(nint(log(real(dim_st_n*dimfac,
     $           kind(0.d0)))/log(2.d0))))
         else
c         dim=dim_st
            dimn=2147483648
         endif

         if (iproc==0) print *,' ilocn: dim,hash table length=',
     $        dim_st_n,dimn

         memoryallocva=memoryallocva+dimn*4
         if (allocated(linhashn)) deallocate(linhashn)
         allocate(linhashn(0:dimn-1))
         linhashn=-1
         maxdispn=0

c         if (allocated(occ_tmpN)) deallocate(occ_tmpN)
c         allocate(occ_tmpN(nuc))

         do i=1,dim_st_n

            if (i==1000000*(i/1000000)) then
               if (iproc==0) print *,' construct_hash_table: i=',i
            endif

            occ_tmp(1:nuc)=ty%occn(:,i)

            intbasN=0
            do cr_1=1,nuc
               wd=(occ_tmp(cr_1)-1)/nbit+1
               ibit=mod(occ_tmp(cr_1)-1,nbit)
               intbasN(wd)=ibset(intbasN(wd),ibit)
            end do

            call get_key(mxnw,intbasN(mxnw+1:2*mxnw),key,dimn)
            if (key<0.or.key>dimn-1) then
               print *,'*** error in construct_hash_table: key,dimn=',
     +              key,dimn
               stop
            endif

            max_el=0
            do
               if (linhashn(key)==-1) then
                  linhashn(key)=i
                  exit
               endif
               key=mod(key+1,dimn)
               max_el=max_el+1
            end do
            if (max_el>maxdispn) maxdispn=max_el
         end do

         if (iproc==0) print *,' ilocn: maximal displacement on proc0=',
     $        maxdispn

      endif

      end subroutine construct_hash_table

      subroutine get_key(nwd,intbas,key,M)
      implicit none
      integer,intent(IN) :: nwd,M
      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: key
c      integer(8) :: intbas_tmp
      integer :: wd,ibit
      integer(8) :: isum,mul
      isum=intbas(1)
      mul=isum
c! ieor,ishft appears to work the best
      isum=ieor(isum,ishft(isum,13))
      isum=ieor(isum,ishft(isum,-17))
      isum=ieor(isum,ishft(isum,5))
cc      isum=ieor(isum,ishftc(isum,13))
cc      isum=ieor(isum,ishftc(isum,-17))
cc      isum=ieor(isum,ishftc(isum,5))
      do wd=2,nwd
         isum=isum+intbas(wd)
         isum=ieor(isum,ishft(isum,13))
         isum=ieor(isum,ishft(isum,-17))
         isum=ieor(isum,ishft(isum,5))
cc         isum=ieor(isum,ishftc(isum,13))
cc         isum=ieor(isum,ishftc(isum,-17))
cc         isum=ieor(isum,ishftc(isum,5))
         mul=mul*intbas(wd)
      end do
      isum=ieor(isum,mul)
      key=abs(mod(isum,int(M,kind=8)))
      end subroutine get_key

c     IN : give the integer basis representation
c     OUT : give the SD index, e.g. i3 in amp(i3,na) where
c     amp is the Lanczos vector. na refers to the initial Lanc Vec.
      subroutine get_state_index(nuc,nwd,intbas,index)
      implicit none
      integer,intent(IN) :: nuc,nwd
      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: index
      integer :: key,i,j,wd,ibit,ip,in
      logical :: OK

      call get_key(nwd,intbas,key,dim)

      index=linhash(key)
      if (index==-1) return
      ip=ty%occ(1,index)
      in=ty%occ(2,index)
      occ_tmp(1:ty%protn)=ty%occp(:,ip)
      occ_tmp(ty%protn+1:nuc)=ty%occn(:,in)

      intbas_tmp=0
      do i=1,nuc
         wd=(occ_tmp(i)-1)/64+1
         ibit=mod(occ_tmp(i)-1,64)
c     print *,' i,wd,ibit=',i,wd,ibit
         intbas_tmp(wd)=ibset(intbas_tmp(wd),ibit)
      end do

      OK=.true.
      do wd=1,nwd
         if (intbas(wd)/=intbas_tmp(wd)) then
            OK=.false.
            exit
         endif
      end do

      if (OK) then
         return
      endif

      do j=1,maxdisp
         key=mod(key+1,dim)
         index=linhash(key)
         if (index==-1) return
         ip=ty%occ(1,index)
         in=ty%occ(2,index)
         occ_tmp(1:ty%protn)=ty%occp(:,ip)
         occ_tmp(ty%protn+1:nuc)=ty%occn(:,in)

         intbas_tmp=0
         do i=1,nuc
            wd=(occ_tmp(i)-1)/64+1
            ibit=mod(occ_tmp(i)-1,64)
c     print *,' i,wd,ibit=',i,wd,ibit
            intbas_tmp(wd)=ibset(intbas_tmp(wd),ibit)
         end do

         OK=.true.
         do wd=1,nwd
            if (intbas(wd)/=intbas_tmp(wd)) then
               OK=.false.
               goto 1000
            endif
         end do

         if (OK) then
            return
         endif

 1000    continue
      end do
      index=-1
      end subroutine get_state_index

c      subroutine get_state_index_loc(nuc,nwd,intbas,index)
c      implicit none
c      integer,intent(IN) :: nuc,nwd
c      integer(8),intent(IN) :: intbas(nwd)
c      integer,intent(OUT) :: index
c      integer :: key,i,j,wd,ibit,ip,in
c      logical :: OK
c      integer(2),allocatable :: occ_tmp_loc(:)
c      integer(8),allocatable :: intbas_tmp_loc(:)

c      call get_key(nwd,intbas,key,dim)

c      index=linhash(key)
c      if (index==-1) return
c      if (allocated(occ_tmp_loc)) deallocate(occ_tmp_loc)
c      allocate(occ_tmp_loc(nuc))
c      if (allocated(intbas_tmp_loc)) deallocate(intbas_tmp_loc)
c      allocate(intbas_tmp_loc(nwd))
c      ip=ty%occ(1,index)
c      in=ty%occ(2,index)
c      occ_tmp_loc(1:ty%protn)=ty%occp(:,ip)
c      occ_tmp_loc(ty%protn+1:nuc)=ty%occn(:,in)

c      intbas_tmp_loc=0
c      do i=1,nuc
c         wd=(occ_tmp_loc(i)-1)/64+1
c         ibit=mod(occ_tmp_loc(i)-1,64)
cc     print *,' i,wd,ibit=',i,wd,ibit                                                                                                                                                                                                                      
c         intbas_tmp_loc(wd)=ibset(intbas_tmp_loc(wd),ibit)
c      end do

c      OK=.true.
c      do wd=1,nwd
c         if (intbas(wd)/=intbas_tmp_loc(wd)) then
c            OK=.false.
c            exit
c         endif
c      end do

c      if (OK) then
c         return
c      endif

c      do j=1,maxdisp
c         key=mod(key+1,dim)
c         index=linhash(key)
c         if (index==-1) return
c         ip=ty%occ(1,index)
c         in=ty%occ(2,index)
c         occ_tmp_loc(1:ty%protn)=ty%occp(:,ip)
c         occ_tmp_loc(ty%protn+1:nuc)=ty%occn(:,in)

c         intbas_tmp_loc=0
c         do i=1,nuc
c            wd=(occ_tmp_loc(i)-1)/64+1
c            ibit=mod(occ_tmp_loc(i)-1,64)
cc     print *,' i,wd,ibit=',i,wd,ibit                                                                                                                                                                                                                      
c            intbas_tmp_loc(wd)=ibset(intbas_tmp_loc(wd),ibit)
c         end do

c         OK=.true.
c         do wd=1,nwd
c            if (intbas(wd)/=intbas_tmp_loc(wd)) then
c               OK=.false.
c               goto 1000
c            endif
c         end do

c         if (OK) then
c            return
c         endif

c 1000    continue
c      end do
c      index=-1
c      end subroutine get_state_index_loc

      subroutine init_hash_table_only(dimfac)
      use nodeinfo
      use spb, only : nwd
      implicit none
      integer,intent(IN) :: dimfac  
      integer,allocatable :: temp_ind(:)
      integer :: i,key,dim_st,nuc,max_el,temp,wd,cr_1,ibit,ip,in,
     $     nbit

      nbit=64
      dim_st=ty%dim_st
      nuc=ty%nucl
c      mxnw=2*ty%nwd
c     MKGK
      mxnw = 2*nwd

c      if (dim_st<1 431 655 765) then
      if (dim_st*dimfac<750 000 000) then
c         dim=dim_st*3/2
         dim=2*(2**(nint(log(real(dim_st*dimfac,kind(0.d0)))
     $        /log(2.d0))))
      else
c         dim=dim_st
         dim=2147483648
      endif

      if (iproc==0) print *,' dim,hash table length=',dim_st,dim

      if (allocated(linhash)) deallocate(linhash) 
      allocate(linhash(0:dim-1))
      linhash=-1
      maxdisp=0
      max_el_constr=0
      dim_constr=0

      if (allocated(intbas)) deallocate(intbas)
      allocate(intbas(mxnw))
      if (allocated(occ_tmp)) deallocate(occ_tmp)
      allocate(occ_tmp(nuc))
      if (allocated(intbas_tmp)) deallocate(intbas_tmp)
      allocate(intbas_tmp(mxnw))
      dim_aloc=dim_st
      if (associated(occ_constr)) deallocate(occ_constr)
      allocate(occ_constr(nuc,dim_aloc))
      ty%occ=>occ_constr

      end subroutine init_hash_table_only


c     MKGK added routine
c     inserts hash table entry but expects intbas as input
c     get both the initial word and final word so you can determine
c     if the added state is a proton or neutron state

      subroutine insert_hash_table_entry_word(nuc,intbasi,intbasf,
     $     nwd2,index)
      use nodeinfo
      use parameters
      use spb, only : nwd, nprotons, nneutrns, nucleons,nhw_min,nhw_boot
      use ibas
      implicit none
      integer,intent(IN) :: nuc,nwd2
      integer(8), intent(IN) :: intbasi(nwd2), intbasf(nwd2)
      integer,intent(OUT) :: index
      integer :: cr_1,ibit,max_el,key,i,i1,wd
c      character(len=100) :: filenmtmp
      character(len=8) :: line
c     MKGK added variables
      integer, parameter :: buffer = 2000000 ! 2E6
c     for intbas -> occ
      integer :: k, pos, ind,kp,kn
c     these occ variables are not the same as in module occ
      integer(2) :: occx(nuc)
      integer(2) :: occp(nprotons), occn(nneutrns)
      logical :: prot_sd, neut_sd
      integer :: dim_aloc_p_read,dim_aloc_n_read
      integer,allocatable :: temp(:,:)
      real(kind(0.0)),allocatable :: temp_delta(:)
      real(kind(0.d0)) :: del

c     NOTES:
c     dim_constr initialized to dim_st = nsd
c     dim_constr_p = nsdp
c     dim_constr_n = nsdn

      call get_key(nwd2,intbasf,key,dim)

      if (key<0.or.key>dim-1) then
         print *,'*** error in construct_hash_table: key,dim=',
     +        key,dim
         stop
      endif
c         print *,' key=',key

c     MKGK : called with nwd2 instead of nwd (as in Jplus)
      call get_state_index(nuc,nwd2,intbasf,index)

      if (index>-1) then
c         print *,' inside insert_hash: index=',index
         if (nhw_boot>=nhw_min) then
            del=delta_IT(index)
            delta_IT(index)=max(deltaIT,del)
         endif
         return
      endif

c     convert intbasf to iloc
      k=0
      do wd=1,2*nwd
         if (intbasf(wd) == 0) cycle
         do pos=0,nbit1
            if (btest(intbasf(wd),pos)) then
               k=k+1
               ind = (pos+1)+nbit*(wd-1)
               occx(k)=ind
            end if
         end do                 !pos
c         if (wd==nwd.and.k/=nprotons) then
c            print*,'Warning: intbasf -> occ had problems'
c            print*, ' proton part of occx = ',occx
c            print*, 'intbasf = ', intbasf
c         endif
      end do                    ! wd

c     safety check
      if (k.ne.nucleons) then 
         print*,'Warning: intbasf -> occ had problems'
         print*, 'occx = ',occx
         print*, 'intbasf = ', intbasf
      end if

      occp(:)=occx(1:nprotons)
      occn(:)=occx(nprotons+1:nucleons)
c     end of conversion

c     now check if the state goes with a proton, neutron or
c     mixed type SD
c      prot_sd = .false.
c      neut_sd = .false.

c      do wd=1,nwd2
c         if (intbasi(wd) .ne. intbasf(wd)) then
c            if (wd > nwd) then 
c               neut_sd = .true.
c* PN
c               dim_constr_n = dim_constr_n + 1
c            else if (wd <= nwd) then
c               prot_sd = .true.
c* PN
c               dim_constr_p = dim_constr_p + 1
c            end if
c         end if
c      end do
c* PN
      call insert_hash_table_entry_N(nprotons,occp,1,kp)
c      kp=pos_in_iloc(occp,nprotons,1)
c      if (kp<1) then
c         prot_sd=.true.
c         dim_constr_p = dim_constr_p + 1
c      endif
      call insert_hash_table_entry_N(nneutrns,occn,2,kn)
c      kn=pos_in_iloc(occn,nneutrns,2)
c      if (kn<1) then
c         neut_sd=.true.
c         dim_constr_n = dim_constr_n + 1
c      endif

c     safety
c      if (.not.(prot_sd.or.neut_sd)) then
c         print*,'WARNING: unable to determine if P or N SD in bhashed.f'
c         print*, 'intbasi = ',intbasi
c         print*, 'intbasf = ',occp,occn
c         print *,' kp=',kp,'   kn=',kn
c         print *,' ilocp,ilocn=',ilocp(:,kp),ilocn(:,kn)
cc         stop
c      end if
c     end of proton/neutron check

c     obviously the state is not found otherwise we would not add it
      if (index == -1) then
         dim_constr = dim_constr+1
         if (dim_constr > dim) then
            print *,'***error: dim_constr,dim=',dim_constr,dim
            print *,' hash table is too small '
            stop
         endif
c     make space in iloc for new states; dump to file and read back
         if (dim_constr > dim_aloc) !.or.dim_constr_p>dim_aloc_p
c     $        .or.dim_constr_n>dim_aloc_n)
     $        then
c            write(line,'(i8)') iproc
c            line=adjustl(line)
c            filenmtmp = 'iloc_'//trim(line)//'.tmp'
c            open(47,file=trim(filenmtmp)//'.realoc',
c     $           form='unformatted',
c     $           status='unknown',action='readwrite')
cc            write(47) ((occ_constr(i,i1),i=1,nuc),i1=1,dim_aloc)

c     save iloc
c            write(47) (iloc(1,i1),i1=1,dim_aloc)
c            write(47) (iloc(2,i1),i1=1,dim_aloc)
cc            write(47) ((ilocp(i,i1),i=1,nprotons),i1=1,dim_aloc_p)
cc            write(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,dim_aloc_n)

c            deallocate(iloc)
cc            deallocate(ilocp)
cc            deallocate(ilocn)
            dim_aloc = dim_aloc + buffer
c* PN
cc            dim_aloc_p_read=dim_aloc_p
cc            dim_aloc_n_read=dim_aloc_n
cc            if (dim_constr_p>dim_aloc_p) then
cc               dim_aloc_p = dim_aloc_p + buffer/2
cc            endif
cc            if (dim_constr_n>dim_aloc_n) then
cc               dim_aloc_n = dim_aloc_n + buffer/2
cc            endif
            allocate(temp(2,dim_aloc))
cc            allocate(ilocp(nprotons,dim_aloc_p))
cc            allocate(ilocn(nneutrns,dim_aloc_n))

            temp(1:2,1:dim_aloc-buffer)=iloc
            call move_alloc(temp,iloc)

            if (nhw_boot>=nhw_min) then
               allocate(temp_delta(dim_aloc))
               temp_delta(1:dim_aloc-buffer)=delta_IT
               call move_alloc(temp_delta,delta_IT)
            endif

c            rewind(47)
c            read(47) (iloc(1,i1),i1=1,dim_aloc-buffer)
c            read(47) (iloc(2,i1),i1=1,dim_aloc-buffer)
cc            read(47) ((ilocp(i,i1),i=1,nprotons),i1=1,dim_aloc_p_read)
cc            read(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,dim_aloc_n_read)
c            close(47)

            ty%occ=>iloc
cc            ty%occp=>ilocp
cc            ty%occn=>ilocn
         endif

         max_el=0
         do
            if (linhash(key)==-1) then
               linhash(key)=dim_constr

c     only proton SD is "new"
c               if (.not.(prot_sd.or.neut_sd)) then
                  iloc(1,dim_constr)=kp
                  iloc(2,dim_constr)=kn

                  if (nhw_boot>=nhw_min) then
                     delta_IT(dim_constr)=deltaIT
                  endif
                  
c               elseif (prot_sd.and..not.neut_sd) then 
c                  ilocp(:,dim_constr_p)=occp(:)
c* PN
c                  iloc(1,dim_constr)=dim_constr_p
c     the neutron sd should already exist
c                  k=pos_in_iloc(occn,nneutrns,2)
c                  if (k < 1) then
c                     print*,'warning: could not locate neutrons'
c                     print*,'check in bhashed.f ~ line 530'
c                  end if
c                  iloc(2,dim_constr)=kn
c     only neutron "SD" is new
c               else if (neut_sd.and..not.prot_sd) then
c                  ilocn(:,dim_constr_n)=occn(:) !dim_constr is not right here
c* PN
c                  iloc(2,dim_constr)=dim_constr_n
c                  k=pos_in_iloc(occp,nprotons,1)
c                  if (k < 1) then
c                     print*,'warning: could not locate protons'
c                     print*,'check in bhashed.f ~ line 550'
c                  end if
c                  iloc(1,dim_constr)=kp
c     this combination of SD's is new for both P and N
c               else if (prot_sd.and.neut_sd) then
c* PN
c                  iloc(1,dim_constr)=dim_constr_p
c                  iloc(2,dim_constr)=dim_constr_n
c                  ilocp(:,dim_constr_p)=occp(:)
c                  ilocn(:,dim_constr_n)=occn(:)
c               end if
               exit
            endif
            key=mod(key+1,dim)
            max_el=max_el+1
         end do
c* PN         if (max_el>max_el_constr) then
         if (max_el>maxdisp) then
c            max_el_constr=max_el
c            maxdisp=max_el_constr
            maxdisp=max_el
            if (iproc==0) print *,' maximal displacement on proc0=',
     $           maxdisp
         endif
         index=dim_constr
c         print *,' dim_constr_p,dim_constr_n,index=',
c     $        dim_constr_p,dim_constr_n,index
      endif
      end subroutine insert_hash_table_entry_word

      subroutine insert_hash_table_entry(nuc,occ,type,index)
      use nodeinfo
      use parameters
      use spb, only : nwd, nprotons, nneutrns, nucleons,nhw_min,nhw_boot
      use ibas
      implicit none
      integer,intent(IN) :: nuc,occ(nuc),type
      integer,intent(OUT) :: index
      integer :: cr_1,ibit,max_el,key,i,i1,wd
c      character(len=100) :: filenmtmp
      character(len=8) :: line
c     MKGK added variables
      integer, parameter :: buffer = 2000000 ! 2E6
c     for intbas -> occ
      integer :: k, pos, ind,kp,kn
c     these occ variables are not the same as in module occ
c      integer(2) :: occx(nuc)
      integer(2) :: occp(nprotons), occn(nneutrns)
      logical :: prot_sd, neut_sd
      integer :: dim_aloc_p_read,dim_aloc_n_read
      integer,allocatable :: temp(:,:)
      real(kind(0.0)),allocatable :: temp_delta(:) 
      real(kind(0.d0)) :: del

c     NOTES:
c     dim_constr initialized to dim_st = nsd
c     dim_constr_p = nsdp
c     dim_constr_n = nsdn

c     convert occ to intbas
      intbas=0
      do cr_1=1,nuc
         wd=(occ(cr_1)-1)/nbit+1
         ibit=mod(occ(cr_1)-1,nbit)
         intbas(wd)=ibset(intbas(wd),ibit)
      end do

      call get_key(2*nwd,intbas,key,dim)

      if (key<0.or.key>dim-1) then
         print *,'*** error in construct_hash_table: key,dim=',
     +        key,dim
         stop
      endif
c         print *,' key=',key

c     MKGK : called with nwd2 instead of nwd (as in Jplus)
      call get_state_index(nuc,2*nwd,intbas,index)

      if (index/=-1) then
cc         print *,' inside insert_hash: index=',index
         if (nhw_boot>=nhw_min) then
            del=delta_IT(index)
            delta_IT(index)=max(deltaIT,del)
         endif
         return
      endif

      occp(:)=occ(1:nprotons)
      occn(:)=occ(nprotons+1:nucleons)
c     end of conversion

c     now check if the state goes with a proton, neutron or
c     mixed type SD
c      prot_sd = .false.
c      neut_sd = .false.

c      do wd=1,nwd2
c         if (intbasi(wd) .ne. intbasf(wd)) then
c            if (wd > nwd) then 
c               neut_sd = .true.
c* PN
c               dim_constr_n = dim_constr_n + 1
c            else if (wd <= nwd) then
c               prot_sd = .true.
c* PN
c               dim_constr_p = dim_constr_p + 1
c            end if
c         end if
c      end do
c* PN
      call insert_hash_table_entry_N(nprotons,occp,1,kp)
c      kp=pos_in_iloc(occp,nprotons,1)
c      if (kp<1) then
c         prot_sd=.true.
c         dim_constr_p = dim_constr_p + 1
c      endif
      call insert_hash_table_entry_N(nneutrns,occn,2,kn)
c      kn=pos_in_iloc(occn,nneutrns,2)
c      if (kn<1) then
c         neut_sd=.true.
c         dim_constr_n = dim_constr_n + 1
c      endif

c     safety
c      if (.not.(prot_sd.or.neut_sd)) then
c         print*,'WARNING: unable to determine if P or N SD in bhashed.f'
c         print*, 'intbasi = ',intbasi
c         print*, 'intbasf = ',occp,occn
c         print *,' kp=',kp,'   kn=',kn
c         print *,' ilocp,ilocn=',ilocp(:,kp),ilocn(:,kn)
cc         stop
c      end if
c     end of proton/neutron check

c     obviously the state is not found otherwise we would not add it
      if (index == -1) then
         dim_constr = dim_constr+1
         if (dim_constr > dim) then
            print *,'***error: dim_constr,dim=',dim_constr,dim
            print *,' hash table is too small '
            stop
         endif
c     make space in iloc for new states; dump to file and read back
         if (dim_constr > dim_aloc) !.or.dim_constr_p>dim_aloc_p
c     $        .or.dim_constr_n>dim_aloc_n)
     $        then
c            write(line,'(i8)') iproc
c            line=adjustl(line)
c            filenmtmp = 'iloc_'//trim(line)//'.tmp'
c            open(47,file=trim(filenmtmp)//'.realoc',
c     $           form='unformatted',
c     $           status='unknown',action='readwrite')
cc            write(47) ((occ_constr(i,i1),i=1,nuc),i1=1,dim_aloc)

c     save iloc
c            write(47) (iloc(1,i1),i1=1,dim_aloc)
c            write(47) (iloc(2,i1),i1=1,dim_aloc)
cc            write(47) ((ilocp(i,i1),i=1,nprotons),i1=1,dim_aloc_p)
cc            write(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,dim_aloc_n)

c            deallocate(iloc)
cc            deallocate(ilocp)
cc            deallocate(ilocn)
            dim_aloc = dim_aloc + buffer
c* PN
cc            dim_aloc_p_read=dim_aloc_p
cc            dim_aloc_n_read=dim_aloc_n
cc            if (dim_constr_p>dim_aloc_p) then
cc               dim_aloc_p = dim_aloc_p + buffer/2
cc            endif
cc            if (dim_constr_n>dim_aloc_n) then
cc               dim_aloc_n = dim_aloc_n + buffer/2
cc            endif
c            allocate(iloc(2,dim_aloc))
cc            allocate(ilocp(nprotons,dim_aloc_p))
cc            allocate(ilocn(nneutrns,dim_aloc_n))

c            rewind(47)
c            read(47) (iloc(1,i1),i1=1,dim_aloc-buffer)
c            read(47) (iloc(2,i1),i1=1,dim_aloc-buffer)
cc            read(47) ((ilocp(i,i1),i=1,nprotons),i1=1,dim_aloc_p_read)
cc            read(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,dim_aloc_n_read)
c            close(47)

            allocate(temp(2,dim_aloc))
            temp(1:2,1:dim_aloc-buffer)=iloc
            call move_alloc(temp,iloc)

            if (nhw_boot>=nhw_min) then
               allocate(temp_delta(dim_aloc))
               temp_delta(1:dim_aloc-buffer)=delta_IT
               call move_alloc(temp_delta,delta_IT)
            endif
            
            ty%occ=>iloc
c            ty%occp=>ilocp
c            ty%occn=>ilocn
         endif

         max_el=0
         do
            if (linhash(key)==-1) then
               linhash(key)=dim_constr

c     only proton SD is "new"
c               if (.not.(prot_sd.or.neut_sd)) then
                  iloc(1,dim_constr)=kp
                  iloc(2,dim_constr)=kn

                  if (nhw_boot>=nhw_min) then
                     delta_IT(dim_constr)=deltaIT
                  endif
                  
c               elseif (prot_sd.and..not.neut_sd) then 
c                  ilocp(:,dim_constr_p)=occp(:)
c* PN
c                  iloc(1,dim_constr)=dim_constr_p
c     the neutron sd should already exist
c                  k=pos_in_iloc(occn,nneutrns,2)
c                  if (k < 1) then
c                     print*,'warning: could not locate neutrons'
c                     print*,'check in bhashed.f ~ line 530'
c                  end if
c                  iloc(2,dim_constr)=kn
c     only neutron "SD" is new
c               else if (neut_sd.and..not.prot_sd) then
c                  ilocn(:,dim_constr_n)=occn(:) !dim_constr is not right here
c* PN
c                  iloc(2,dim_constr)=dim_constr_n
c                  k=pos_in_iloc(occp,nprotons,1)
c                  if (k < 1) then
c                     print*,'warning: could not locate protons'
c                     print*,'check in bhashed.f ~ line 550'
c                  end if
c                  iloc(1,dim_constr)=kp
c     this combination of SD's is new for both P and N
c               else if (prot_sd.and.neut_sd) then
c* PN
c                  iloc(1,dim_constr)=dim_constr_p
c                  iloc(2,dim_constr)=dim_constr_n
c                  ilocp(:,dim_constr_p)=occp(:)
c                  ilocn(:,dim_constr_n)=occn(:)
c               end if
               exit
            endif
            key=mod(key+1,dim)
            max_el=max_el+1
         end do
c* PN         if (max_el>max_el_constr) then
         if (max_el>maxdisp) then
c            max_el_constr=max_el
c            maxdisp=max_el_constr
            maxdisp=max_el
            if (iproc==0) print *,' maximal displacement on proc0=',
     $           maxdisp
         endif
         index=dim_constr
c         print *,' dim_constr_p,dim_constr_n,index=',
c     $        dim_constr_p,dim_constr_n,index
      endif
      end subroutine insert_hash_table_entry

      subroutine insert_hash_table_entry_N(nuc,occ,type,index)
      use ibas, only: ilocp,ilocn
      use nodeinfo
      implicit none
      integer,intent(IN) :: nuc,type
      integer(2),intent(IN) :: occ(nuc)
      integer,intent(OUT) :: index
      integer :: cr_1,ibit,max_el,key,i,i1,wd,nbit
c      character(len=100) :: filenmtmp
      character(len=8) :: line
      integer(2),allocatable :: temp(:,:)

      nbit=64

      intbasN=0
      do cr_1=1,nuc
         wd=(occ(cr_1)-1)/nbit+1
         ibit=mod(occ(cr_1)-1,nbit)
         intbasN(wd)=ibset(intbasN(wd),ibit)
      end do

      if (type==1) then
         call get_key(mxnw,intbasN(1:mxnw),key,dimp)
         if (key<0.or.key>dimp-1) then
            print *,'*** error in construct_hash_table: key,dimp=',
     +           key,dimp
            stop
         endif
      else
         call get_key(mxnw,intbasN(mxnw+1:2*mxnw),key,dimn)
         if (key<0.or.key>dimn-1) then
            print *,'*** error in construct_hash_table: key,dimn=',
     +           key,dimn
            stop
         endif
      endif
c         print *,' key=',key

      call get_state_index_N(nuc,occ,type,index)

      if (index==-1) then
         if (type==1) then
            dim_constr_p=dim_constr_p+1
            if (dim_constr_p>dimp) then
               print *,'***error: dim_constr_p,dimp=',dim_constr_p,dimp
               stop
            endif
            if (dim_constr_p>dim_aloc_p) then
c               write(line,'(i8)') iproc
c               line=adjustl(line)
c               filenmtmp='ilocp_'//trim(line)//'.tmp'
c               open(47,file=trim(filenmtmp)//'.realoc',
c     $              form='unformatted',
c     $              status='unknown',action='readwrite')
c               write(47) ((ilocp(i,i1),i=1,nuc),i1=1,dim_aloc_p)
c               deallocate(ilocp)
               dim_aloc_p=dim_aloc_p+100000
c               allocate(ilocp(nuc,dim_aloc_p))
c               rewind(47)
c               read(47) ((ilocp(i,i1),i=1,nuc),i1=1,dim_aloc_p-100000)
c               close(47)

               allocate(temp(nuc,dim_aloc_p))
               temp(1:nuc,1:dim_aloc_p-100000)=ilocp
               call move_alloc(temp,ilocp)

               ty%occp=>ilocp
            endif
            max_el=0
            do
               if (linhashp(key)==-1) then
                  linhashp(key)=dim_constr_p
                  ilocp(:,dim_constr_p)=occ(:)
                  exit
               endif
               key=mod(key+1,dimp)
               max_el=max_el+1
            end do
            if (max_el>maxdispp) then
               maxdispp=max_el
               if (iproc==0) print *,
     $              ' ilocp: maximal displacement on proc0=',
     $              maxdispp
            endif
            index=dim_constr_p
         else
            dim_constr_n=dim_constr_n+1
            if (dim_constr_n>dimn) then
               print *,'***error: dim_constr_n,dimn=',dim_constr_n,dimn
               stop
            endif
            if (dim_constr_n>dim_aloc_n) then
c               write(line,'(i8)') iproc
c               line=adjustl(line)
c               filenmtmp='ilocn_'//trim(line)//'.tmp'
c               open(47,file=trim(filenmtmp)//'.realoc',
c     $              form='unformatted',
c     $              status='unknown',action='readwrite')
c               write(47) ((ilocn(i,i1),i=1,nuc),i1=1,dim_aloc_n)
c               deallocate(ilocn)
               dim_aloc_n=dim_aloc_n+100000
c               allocate(ilocn(nuc,dim_aloc_n))
c               rewind(47)
c               read(47) ((ilocn(i,i1),i=1,nuc),i1=1,dim_aloc_n-100000)
c               close(47)

               allocate(temp(nuc,dim_aloc_n))
               temp(1:nuc,1:dim_aloc_n-100000)=ilocn
               call move_alloc(temp,ilocn)

               ty%occn=>ilocn
            endif
            max_el=0
            do
               if (linhashn(key)==-1) then
                  linhashn(key)=dim_constr_n
                  ilocn(:,dim_constr_n)=occ(:)
                  exit
               endif
               key=mod(key+1,dimn)
               max_el=max_el+1
            end do
            if (max_el>maxdispn) then
               maxdispn=max_el
               if (iproc==0) print *,
     $              ' ilocn: maximal displacement on proc0=',
     $              maxdispn
            endif
            index=dim_constr_n
         endif
      endif
      end subroutine insert_hash_table_entry_N

      subroutine get_state_index_N(nuc,occ,type,index)
c      subroutine get_state_index(nuc,nwd,intbas,index)
      implicit none
      integer,intent(IN) :: nuc,type
      integer(2),intent(IN) :: occ(nuc)
c      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: index
      integer :: key,i,j,wd,ibit,ip,in
      logical :: OK
      integer :: maxdis

      intbasN=0
      do i=1,nuc
         wd=(occ(i)-1)/64+1
         ibit=mod(occ(i)-1,64)
         intbasN(wd)=ibset(intbasN(wd),ibit)
      end do

      if (type==1) then
         call get_key(mxnw,intbasN(1:mxnw),key,dimp)
         index=linhashp(key)
      else
         call get_key(mxnw,intbasN(mxnw+1:2*mxnw),key,dimn)
         index=linhashn(key)
      endif
      if (index==-1) return
      if (type==1) then
         occ_tmp(1:nuc)=ty%occp(:,index)
         maxdis=maxdispp
      else
         occ_tmp(1:nuc)=ty%occn(:,index)
         maxdis=maxdispn
      endif
      OK=.true.
      do i=1,nuc
         if (occ_tmp(i)/=occ(i)) then
            OK=.false.
            exit
         endif
      end do

      if (OK) then
         return
      endif

      do j=1,maxdis
         if (type==1) then
            key=mod(key+1,dimp)
            index=linhashp(key)
         else
            key=mod(key+1,dimn)
            index=linhashn(key)
         endif
         if (index==-1) return
         if (type==1) then
            occ_tmp(1:nuc)=ty%occp(:,index)
         else
            occ_tmp(1:nuc)=ty%occn(:,index)
         endif

         OK=.true.
         do i=1,nuc
            if (occ_tmp(i)/=occ(i)) then
               OK=.false.
               goto 1000
            endif
         end do

         if (OK) then
            return
         endif

 1000    continue
      end do
      index=-1
      end subroutine get_state_index_N

      end module hash_tables

      program ncsd
      use parameters
      use nodeinfo
      use iters, only: nite,memsave
      use ITvar
      use spb
      use hmes, only: memoryalloc,memoryallocva
      implicit none
      include 'mpif.h'
c     MKGK
      integer :: ierr
      integer :: nhw_begin
      real(8) cputime,walltime,timeused
      integer omp_get_num_threads,omp_get_thread_num
      character(len=80) :: line,line1
      integer :: Nmin_HO,nminho
      integer :: ii
      logical :: kappadepwr
      
      call MPI_INIT(ierr)
      icomm = MPI_COMM_WORLD
      call MPI_COMM_RANK(icomm,iproc,ierr)
      call MPI_COMM_SIZE(icomm,nproc,ierr)
c*      if (iproc==0) then
c*      call clock(cputime,walltime,0,0)
c*      timeused = walltime
c*      endif

c----------------------------------------------------------------------------
c     ----- Input files:
c     Unit # 1 (mfdp.dat)  --- Input parameters.
c     Unit # 2 (TBME*.int) --- s.p. energies and TBME's for Trel, Hrel Clmb, G.
c     Search for "open" for more details.

c     Unit # 13 (mfdp.egv)
c     ----- Storage of true eigenvectors and other information needed for
c     ----- density or spectroscopic factor calculations

c     ----- If irest=4 (to restart after the last iteration and to
c     ----- evaluate expectations and transitions), input files also include:
c     ----- Unit # 14 (lzegv.tmp) which stores the eigenstates in the 
c     ----- Lanczos basis.

c     ----- For irest=0, units # 13 are output files.

c-----------------------------------------------------------------------
c     ----- Output files:
      if(iproc.eq.0)then
c     ----- detailed output
         open(8,file='mfd.log',status='unknown')  
      endif
c-----------------------------------------------------------------------
      call Input
      if(iproc.eq.0)write(6,*)' --- Input  called ---' 

c**      print *,' nite,iproc=',nite,iproc

      if (iproc==0) inquire(unit=675,opened=kappadepwr)

      if (iproc<=nite) then
c     ----- Temporary files ------
         write(line,'(i8)') iproc
         line=adjustl(line)
         line='lancz'//trim(line)//'.tmp'
         open(12,file=trim(line),access='sequential',form='unformatted',
     $        status='unknown')
      endif
c      if (nproc>1.and.iproc==nproc) then
c         write(line,'(i8)') iproc
c         line=adjustl(line)
c         line='ampb'//trim(line)//'.tmp'
c         open(15,file=trim(line),access='sequential',form='unformatted',
c     $        status='unknown')
c      endif

      call MPI_Barrier(icomm,ierr)

c     ======== Bootstrap loop ================
c     Do a loop over nmax's:
c     start at nmax=0 and work up to nmax wanted.
      
      if (nhw_restart == -1) then
         nhw_begin = nhw0
      else
         nhw_begin=nhw_restart
      endif

      do nhw_boot = nhw_begin, nhw_max, 2
         nhw = nhw_boot
         if (iproc == 0) then
            print*
            print*, ' ** nhw_boot = ', nhw_boot
         end if         
c* PN 
         if (iproc==0) then
            nminho=Nmin_HO(nprotons)+Nmin_HO(nneutrns)
            write(9,"(/,' Nhw=',i4,'  Nmax=',i4 )") nhw_boot,
     $           nhw_boot-nminho
            if (kappadepwr) write(675,"('# ',i4)") nhw_boot
         endif
c**
         memoryallocva=0
c     no importance truncation (full spaces)
         if (nhw_boot < nhw_min) then
            call Setup
            if(iproc.eq.0)write(6,*)' --- Setup  called ---' 

            call make_hash_iloc(.true.)
            if (iproc==0) print *, ' --- make_hash_iloc called ---'

            call MPI_Barrier(icomm,ierr)

            if(iproc.eq.0)write(6,*)' --- Doing  Lanczos ---'
            call Lanczos_hash
            call MPI_Barrier(icomm,ierr)
         
            if (iproc==0) print *, '--saving pivot--'
            call save_pivot
      
            if(iproc.eq.0)write(6,*)' --- Doing  Transitions ---' 
            if (iproc==0) then
c     ----- true eigenvectors
               write(line,'(i8)') nhw
               line=adjustl(line)
               if (memsave==1) then
                  open(13,file='ncsd_light_pn_'//trim(line)//'.egv',
     +                 form='unformatted',status='unknown')
               else
                  open(13,file='mfdp_'//trim(line)//'.egv',
     +                 form='unformatted',status='unknown')
               endif
            endif
            call Transitions
            call cleaner
            call file_open_close
c**
            if (iproc==0) then
               close(13) 
            endif

*C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nprocomp,iprocomp)
*      nprocomp=omp_get_num_threads()
*      iprocomp=omp_get_thread_num()
*      print *,' MPI process #',iproc,'  MPI nproc=',nproc,
*     +     ' OMP process #',iprocomp,'  OMP nproc=',nprocomp
*C$OMP END PARALLEL

c     end of if boot < nhw_min
c** PN
         else if (nhw_boot .ge. nhw_min) then
            nhw_hold=nhw-2
            if (irest==0) then
               kappa=kappa_store(1)
               call importance_truncate
               if (iproc==0)
     $              print *, ' --- importance_truncate called ---'
               call MPI_Barrier(icomm,ierr)
            elseif (irest==6.and.nhw_restart==nhw_boot
     $              .and.kappa_restart==-1) then
               kappa=kappa_store(1)
               call iloc_read(nhw_boot-2)
               if (iproc==0) print *, ' --- iloc_read called ---'
               call MPI_Barrier(icomm,ierr)
               call make_hash_iloc(.false.)
               if (iproc==0)
     $              print *, ' --- make_hash_iloc called ---'
               call importance_truncate
               if (iproc==0)
     $              print *, ' --- importance_truncate called ---'
               call MPI_Barrier(icomm,ierr)
               kappa_restart=kappa_points
            endif
            do ii=kappa_points,1,-1
               kappa=kappa_store(ii)
               if (iproc==0) then
                  print *,' kappa_min=',kappa
                  write(9,"(/,' kappa_min=',d10.3)") kappa
               endif
               if ((irest==4.or.irest==6).and.nhw_boot==nhw_restart
     $              .and.kappa_restart<ii) then
                  cycle
               elseif (irest==0.or.irest==6.or.
     $                 (irest==4.and.nhw_boot==nhw_restart
     $                 .and.kappa_restart==ii)) then
                  memoryallocva=0
                  call iloc_read(nhw_boot)
                  if (iproc==0) print *, ' --- iloc_read called ---'
                  call MPI_Barrier(icomm,ierr)
                  call make_hash_iloc(.false.)
                  if (iproc==0)
     $                 print *, ' --- make_hash_iloc called ---'
                  call MPI_Barrier(icomm,ierr)
c               else
c                  call importance_truncate
c                  if (iproc==0)
c     $                 print *, ' --- importance_truncate called ---'
c                  call MPI_Barrier(icomm,ierr)
               endif
               if (iproc.eq.0) write(6,*) ' --- Doing  Lanczos ---'
               call Lanczos_hash
               call MPI_Barrier(icomm,ierr)
         
               if (ii==1) then
                  if (iproc==0) print *, '--saving pivot--'
                  call save_pivot
               endif

               if(iproc.eq.0)write(6,*)' --- Doing  Transitions ---' 
               if (iproc==0) then
c     ----- true eigenvectors
                  write(line,'(i8)') nhw
                  line=adjustl(line)
                  write(line1,'(d10.3)') kappa
                  line1=adjustl(line1)
                  if (memsave==1) then
                     open(13,file=
     $  'ncsd_light_pn_'//trim(line)//'_'//trim(line1)//'.egv',
     +                    form='unformatted',status='unknown')
                  else
                     open(13,file=
     $                    'mfdp_'//trim(line)//'_'//trim(line1)//'.egv',
     +                    form='unformatted',status='unknown')
                  endif
               endif
               call Transitions
               call cleaner
               call file_open_close
c**
               if (iproc==0) then
                  close(13) 
               endif
               nhw_hold=nhw
            end do
         end if
         if (iproc==0.and.kappadepwr) write(675,"('&')")
      end do    ! end of bootstrap loop
      
c     HERE

      if(iproc==0)then
         call clock(cputime,walltime,9,1)
c*      timeused = walltime - timeused ! End timer

         write(8,*)'Calculations completed without fatal error.'
c*         if (timeused>0.d0) then
c*            write(6,*)'Total time used in this calculation (sec):',
c*     &           timeused
c*            write(8,*)'Total time used in this calculation (sec):',
c*     &           timeused
c*         write(9,*)'Total time used in this calculation (sec):',
c*     &        timeused
c*         endif   
         close(8) 
         close(9)
         inquire(unit=675,opened=kappadepwr)
         if (kappadepwr) close(675)
      endif   
c      close(10)
c      close(11)
      close(12)
      call MPI_Barrier(icomm,ierr)
      call MPI_Finalize(ierr)
      stop
      end


c     Read input from the file "mfdp.dat" and do some initializations:
      subroutine Input
      use parameters
      use config
      use nodeinfo
      use multig
      use ctrl
      use jsgna
      use spodata
      use nutbme
      use idxb
      use consts
      use iters
      use spb
      use spsdata
      use effop
      use ITvar
      use albe, only : conset
      use tbmepar
      use v3b, only: E1max,E12max,E123max
      use TUD_tbme
      use pn_tbme
      use hmes, only: memavail
      include 'mpif.h'
      character(len=120) :: outfile
c     mnpi(N): minimum number of protons  in the N-th shell
c     mnnu(N): minimum number of neutrons in the N-th shell
c     mxpi(N): maximum number of protons  in the N-th shell
c     mxnu(N): maximum number of neutrons in the N-th shell
c     nsps0(N): number of s.p. states in the N-th shell
c     nspsb(N): number of s.p. states seen before the N-th shell
c     npib(N):  number of protons  seen before the N-th shell
c     nnub(N):  number of neutrons seen before the N-th shell
c     nops(N):  number of proton  many-body states in the N-th shell
c     nons(N):  number of neutron many-body states in the N-th shell
c----------------------------------------------------------------------------
      integer, allocatable :: idxord(:),idxtmp(:)
c----------------------------------------------------------------------------
      real(8) cputime,walltime
      logical :: tbmebin
      real(kind(0.d0)) :: hbo_in,lambda_in,rmemavail

      sqrt2 = sqrt(2.0)
c
c     ----- Read input from unit 1:
      open(1, file='mfdp.dat',status='old')
      read(1,*)                 ! Title line
      read(1,*)
      read(1,'(a)') intfile      ! name of the interaction file
      read(1,*)intform
      read(1,'(a)') outfile      ! name of the main output file
      if (iproc==0) then
         open(9,file=trim(outfile),status='unknown')
         call clock(cputime,walltime,9,1)
         write(9,"(' Number of MPI tasks =',i5,/)") nproc
      endif   
      if (abs(intform)==3) then
         tbmeTUD=.true.
         if (intform==-3) no2bv3n=.true.
      elseif (intform==4) then
         tbmepn=.true.
      endif
      read(1,*)nprotons,nneutrns,hbomeg ! Z, N and hbar*Omega:
      read(1,*)nhw,iparity,mjtotal
c     ----- Read in min(N1 or N2), max(N1 or N2) and max(N1+N2),  here N=2n+l:

      if (tbmeTUD) then
         read(1,*)mnn2l,mxn2l,mxnn,n_Max,l_Max
         N1_max_tud=mxn2l
         N12_max_tud=mxnn
      else
         read(1,*)mnn2l,mxn2l,mxnn ! Restrictions on many-body states
         if (tbmepn) then
            N1_max_pn=mxn2l
            N12_max_pn=mxnn
         endif
      endif

      N12_max=mxnn

      jmx=mxnn+1                      ! Maximal possible two-nucleon J

      read(1,*)iham,iclmb,strcm
      if (iproc==0) then
         write(9,91)nprotons,nneutrns,hbomeg
 91      format(3x,'Z =',i3,'   N =',i3,'   hbar*Omega =',f8.5)
         write(9,92)nhw,iparity,mjtotal
 92      format(3x,'Nsum =',i2,'   Parity =',i2,'   2Msum =',i3)
         write(9,93)mnn2l,mxn2l,mxnn
 93      format(2x,' NN interaction triangle: min(N) =',i2,
     $        '   max(N) =',i2,'   max(N1+N2) =',i3)
         if (tbmeTUD) write(9,
     $        "('   TUD TBME ordering used; n_Max,l_Max=',2i4)")
     $        n_Max,l_Max
         if (tbmepn) write(9,"('   Proton-neutron TBME used')")
         if(iham.eq.1) write(9,*) '  Use H = Trel + G.'
         if(iham.eq.2) write(9,*) '  Use H = Hrel + G.'
         if(iham.eq.3) write(9,*) '  Use H = SPE + G.'
c         if(iclmb.eq.0)write(9,*) '  Without the Coulomb interaction.'
c         if(iclmb.eq.1)write(9,*) '  With the Coulomb interaction.'
         if(strcm.eq.0.)write(9,*)'  Without the c.o.m. projection.'
         if(strcm.gt.0.)write(9,3867)  strcm
 3867    format('   With the c.o.m. projection: beta= ',f8.4)
         write(9,*) '  Interaction file ',trim(intfile)
      endif    
c***      read(1,*)nbits            ! 32 bits or 64 bits (for machine)
c************ nbits not used **** pn **
      read(1,*) mnop
      if (mnop/=2.and.abs(mnop)/=3) then
         if (iproc==0) print *,'*** warning: mnop set to 2!!!****'
         mnop=2
      endif
      read(1,*) major
      read(1,*) nshll0

      nucleons = nprotons + nneutrns
      mttotal = nprotons - nneutrns
c     
      itrel = 0
      ihrel = 0
      ispe  = 0
      if(iham.eq.1) itrel=1
      if(iham.eq.2) ihrel=1
      if(iham.eq.3) ispe =1
      if(mnn2l.eq.mxn2l) strcm = 0.0
c
c     ----- If major = 0, then nshll0 = number of subshells in the model space
c     ----- If major > 0, nshll0.le.0 to order j's in the increasing order
c     -----               nshll0  > 0 to order j's in the decreasing order

      if (major.eq.0)then
         nshll = nshll0         ! Number of subshells 
         if(iproc.eq.0)write(9,*)'  With subshell occ. restriction.'
      else
         jorder = nshll0        ! Tells how the subshells in a major shell
                                ! should be ordered.
         if(iproc.eq.0)write(9,*)
     &        '  With major shell occupation restriction.'
         nshll = mxn2l - mnn2l + 1 ! Number of major shells
      endif
c     

      allocate(mnpi(nshll),mxpi(nshll),mnnu(nshll),mxnu(nshll))
      allocate(mnpn(nshll),mxpn(nshll),ntempi(nshll),ntemnu(nshll))
      allocate(nsps0(nshll),nspsb(nshll),npib(nshll),nnub(nshll))
      allocate(nops(nshll),nons(nshll),neng(nshll))

      if(iproc.eq.0)then
         write(9,*)'  Number of shells = ',nshll
         write(9,*)'  Occupation restrictions:'
         write(9,*)'  (Shell could be either a major or a sub shell)'
         write(9,*)
     &    '                   min#p max#p min#n max#n  min#pn max#pn'
      endif
c     Read in occupation restrictions:
      read(1,*) 
c     ----------------------------------------------------------------------
      if(major.eq.0)then        ! Read in subshell occupation restrictions
         if (allocated(iorder)) deallocate(iorder)
         allocate(iorder(nshll))
         allocate(idxtmp(nshll))
         allocate(idxord(nshll))
         allocate(nnl(nshll),lx(nshll),j2x(nshll))
         nasps=0
         do ishl = 1, nshll
            read(1,*)idx0,nnl(ishl),lx(ishl),j2x(ishl),
     &           mnpi(ishl),mxpi(ishl),mnnu(ishl),mxnu(ishl),
     &           mnpn(ishl),mxpn(ishl)
            nasps=nasps+j2x(ishl)+1
            if (iproc==0) then
               write(9,95)nnl(ishl),lx(ishl),j2x(ishl),
     &              mnpi(ishl),mxpi(ishl),mnnu(ishl),mxnu(ishl),
     &              mnpn(ishl),mxpn(ishl)
 95            format(2x,'(N,l,2j)=',3i2,4i6,2i7)
            endif
c     ----- idx0: index used in the TBME.int file to denote s.p. orbits
            if(idx0.gt.66) stop 'Increase the dimension of idxobt().'
            idxobt(idx0) = ishl
            idxord(ishl) = idx0
         enddo
         
c     ----- Work out the proper order of the single-particle energies:
         do i = 1, nshll
            idxtmp(i) = idxord(i)
         enddo
         do i = 1, nshll
            ismall = idxtmp(i)
            do j = i+1, nshll
               if(idxtmp(j).lt.ismall)then
                  ismall = idxtmp(j)
                  idxtmp(j) = idxtmp(i)
                  idxtmp(i) = ismall
               endif
            enddo
         enddo
         do i = 1, nshll
            iorder(i) = i
            itmp = idxtmp(i)
            do j = 1, nshll
               if(idxord(j).eq.itmp) iorder(i) = j
            enddo
         enddo
         if(iproc.eq.0)then
            write(9,*)'  Order of the single-particle orbits:'
            write(9,*)'    Input              Output:'
            write(9,104)(idxord(i),iorder(i),i=1,nshll)
 104        format(6x,i3,16x,i3)
         endif
         nwd = (nasps-1)/nbit + 1
         mxsps=nwd*nbit
         mxsps2=2*mxsps
         allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +        nobt_sp(mxsps2))
         allocate(n2l_sp(mxsps2),m2_sp(mxsps2),mt2_sp(mxsps2))
      endif
c     ----------------------------------------------------------------------
      if(major>0)then        ! Read in major shell occupation restrictions
         do ishl = 1, nshll
            read(1,*)mnpi(ishl),mxpi(ishl),mnnu(ishl),mxnu(ishl),
     &           mnpn(ishl),mxpn(ishl)
            if(iproc.eq.0)then
               write(9,102)ishl,mnpi(ishl),mxpi(ishl),
     &              mnnu(ishl),mxnu(ishl),mnpn(ishl),mxpn(ishl)
 102           format(1x,'   Shell #',i2,4x,6i6)
            endif
         enddo
      endif
c     
c     ----------------------------------------------------------
c     ----- Establish quantum numbers for proton s.p. states:
      if(major>0)then
         isps = 0               ! Running sum of the number of s.p. states
         iobt = 0               ! Running sum of the number of s.p. orbitals
         ishl = 0                ! shell index, from 1 to nshll
         do nmj = mnn2l, mxn2l
            ishl = ishl + 1
            lmin = nmj - 2*(nmj/2)
            if(jorder.le.0)then
               lstart = lmin
               lend = nmj
               lstep = 2
            else
               lstart = nmj
               lend = lmin
               lstep = -2
            endif
            do lrun = lstart, lend, lstep
               nspx = (nmj-lrun)/2
               llr = lrun + lrun
               j2min = max(1,llr-1)
               if(jorder.le.0)then
                  jstart = j2min
                  jend = llr+1
                  jstep = 2
               else
                  jstart = llr+1
                  jend = j2min
                  jstep = -2
               endif
               do j2run = jstart, jend, jstep
                  iobt = iobt + 1
                  do m2run = j2run, -j2run, -2
                     isps = isps + 1
                  enddo
               enddo
            enddo
         enddo
         nobt = iobt
         allocate(nnl(nobt),lx(nobt),j2x(nobt))
         nasps=isps
         nwd = (nasps-1)/nbit + 1
         mxsps=nwd*nbit
         mxsps2=2*mxsps
         allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +        nobt_sp(mxsps2))
         allocate(n2l_sp(mxsps2),m2_sp(mxsps2),mt2_sp(mxsps2))

         isps = 0               ! Running sum of the number of s.p. states
         iobt = 0               ! Running sum of the number of s.p. orbitals
         ishl = 0                ! shell index, from 1 to nshll
         do nmj = mnn2l, mxn2l
            ishl = ishl + 1
            nspsb(ishl) = min(isps,mxsps)
            neng(ishl) = nmj
            lmin = nmj - 2*(nmj/2)
            if(jorder.le.0)then
               lstart = lmin
               lend = nmj
               lstep = 2
            else
               lstart = nmj
               lend = lmin
               lstep = -2
            endif
            do lrun = lstart, lend, lstep
               nspx = (nmj-lrun)/2
               llr = lrun + lrun
               j2min = max(1,llr-1)
               if(jorder.le.0)then
                  jstart = j2min
                  jend = llr+1
                  jstep = 2
               else
                  jstart = llr+1
                  jend = j2min
                  jstep = -2
               endif
               do j2run = jstart, jend, jstep
                  iobt = iobt + 1
                  nnl(iobt) = nmj
                  lx(iobt) = lrun
                  j2x(iobt) = j2run
c***                  do m2run = j2run, -j2run, -2
                  do m2run = -j2run, j2run, 2
                     isps = isps + 1
                     if(isps.gt.mxsps)goto 97
                     nobt_sp(isps) = iobt
                     n2l_sp(isps) = nmj
                     n_sp(isps) = nspx
                     l_sp(isps) = lrun
                     j2_sp(isps) = j2run
                     m2_sp(isps) = m2run
                     mt2_sp(isps) = 1
 97                  continue
                  enddo
               enddo
            enddo
         enddo

      endif
      if(major.eq.0)then
         nobt = nshll
         isps = 0
         do ishl = 1, nshll
            nmj = nnl(ishl)
            lrun = lx(ishl)
            j2run = j2x(ishl)
            nspx = (nmj-lrun)/2
            nspsb(ishl) = min(isps,mxsps)
            neng(ishl) = nmj    ! Unperturbed one-body energy of this shell.
                                ! Could be different from (2n+l) if necessary
                                ! (e.g., when imposing subshell restrictions).
            do m2run = j2run, -j2run, -2
               isps = isps + 1
               if(isps.gt.mxsps)goto 98
               nobt_sp(isps) = ishl
               n2l_sp(isps) = nmj
               n_sp(isps) = nspx
               l_sp(isps) = lrun
               j2_sp(isps) = j2run
               m2_sp(isps) = m2run
               mt2_sp(isps) = 1
 98            continue
            enddo
         enddo
      endif
c
c     ----- nasps = number of active s.p. states.
      if (nasps/=isps)then
         if(iproc.eq.0)then
            write(6,*)'  nasps =',nasps,'  isps =',isps
            write(8,*)'  nasps =',nasps,'  isps =',isps
            call MPI_Abort(icomm,101,ierr)
            stop
         endif
      endif
c     
c     ----- Establish quantum numbers for neutron s.p. states:
      do isps = 1, nasps
         ispsnu = isps + mxsps
         nobt_sp(ispsnu) = nobt_sp(isps)
         n2l_sp(ispsnu) = n2l_sp(isps)
         n_sp(ispsnu) = n_sp(isps)
         l_sp(ispsnu) = l_sp(isps)
         j2_sp(ispsnu) = j2_sp(isps)
         m2_sp(ispsnu) = m2_sp(isps)
         mt2_sp(ispsnu) = -1
      enddo
c
      nsps0(nshll) = nasps - nspsb(nshll)
      do i = 1, nshll-1
         nsps0(i) = nspsb(i+1)-nspsb(i)
      enddo
c****************************************************************************
c     ----- Read number of G matrices, minimum unperturbed energies of the 
c     ----- spectators (0 for spectators in the 0s1/2 shell), number of G's 
c     ----- to skip, and iset1 (search for "iset1" for its usage)
      read(1,*)
      read(1,*)nsets,nespmin,nskip,iset1
      read(1,*)
      read(1,*)ki,kf,nf
      read(1,*)egs              ! the g.s. energy
      read(1,*)
      read(1,*)nite             ! Number of iterations
      read(1,*)igt              ! igt=1 to compute GT+/-, 0 for not to.
      read(1,*)
      read(1,*)irest            ! irest=4: Restarting; irest=0: Initial run; irest=6: IT iloc available
      if (irest/=0.and.irest/=4.and.irest/=6) irest=0
      read(1,*) memsave         ! =1: save memory during the second reorthogonalization

c     MKGK
c**   Read importance truncation variables
      read(1,*)
      read(1,*) nhw0, nhw_min, nhw_max, nhw_restart
c* PN
c      read(1,*) kappa, cmin, kappa_restart
      read(1,*) kappa_points, cmin, kappa_restart
      allocate(kappa_store(kappa_points))
      read(1,*) (kappa_store(i),i=1,kappa_points)
c     set convergence check in keV
      read(1,*) conset
c* PN
c      kappa = kappa*1E-5 
      kappa_store = kappa_store*1.d-5 
      cmin = cmin*kappa
      conset=conset/1000.d0

c* PN
      nhw_min=max(nhw_min,nhw0+2)  ! Required
c**

c     output information
      if (iproc == 0.and.nhw_min<=nhw_max) then
      write(9,*) 
      write(9,*) 'Importance truncation variables:'
      write(9,*)
      write(9,*) 'nhw0, nhw_min, nhw_max, nhw_restart:',nhw0,nhw_min,
     $     nhw_max, nhw_restart
      write(9,*) 'number of kappa points=',kappa_points
      write(9,*) 'kappa:'
      write(9,'(5d10.3)') kappa_store
c      write(9,"(' convergence check in keV:',d11.4)") conset
      write(9,*)

      write(6,*) 
      write(6,*) 'Importance truncation variables:'
      write(6,*)
      write(6,*) 'nhw0, nhw_min, nhw_max, nhw_restart:',nhw0,nhw_min,
     $     nhw_max, nhw_restart
      write(6,*) 'number of kappa points=',kappa_points
      write(6,*) 'kappa:'
      write(6,'(5d10.3)') kappa_store
c      write(6,"(' convergence check in keV:',d11.4)") conset
      
      open(675,file=trim(outfile)//'_kappa_dep.agr',
     $     status='unknown',form='formatted',position='append')

      end if
C**   end of importance truncation variables


      if (abs(mnop)==3) then
         if (iproc==0) write(9,*) 
     +                '  *** Three-body interaction calculation ***'
         read(1,'(a)') v3intfile   ! name of the interaction file
         if (iproc==0) write(9,*) '  Three-body interaction file ', 
     +        trim(v3intfile)
         read(1,*) E1max,E12max,E123max
         if (iproc==0)
     $        write(9,"('  Three-body interaction file cuts:',3i4)")
     $        E1max,E12max,E123max
         if (E1max/=mxn2l.or.E12max/=N12_max) then
            print *,' NN and 3N file do not match:',
     $           mxn2l,N_12max,E1max,E12max
            stop
         endif
      else
         read(1,'(a)') v3intfile
         read(1,*)
      endif
      read(1,*)
      read(1,*)epi,enu          ! Effective charges
      read(1,*)glp,gln,gsp,gsn  ! Effective g-factors

      read(1,*,end=7777,err=7777) piv_saved
!!! input in GB in real     
      read(1,*,end=7777,err=7777) rmemavail
      memavail=nint(1024.d0*rmemavail,kind(8))
      if (iproc==0) then
         print *,' rmemavail,memavail=',rmemavail,memavail
         write(9,"(/,'   Total memory available per MPI task=',
     $i8,' MB')") memavail
      endif
      
      close(1)                  ! End of reading input
 7777 continue

      if (no2bv3n) then
         if (iproc==0) then
            write(9,"(
     $'  *** Normal-ordered 3N interaction calculation
     $in NO2B approximation ***')")
            write(9,"('  NO2B interaction file ',a)") trim(v3intfile)
         endif
      endif
         
      if(iproc.eq.0)then
         if (nespmin>0.or.nskip>0.or.iset1>0) then
            write(9,99)nsets,nespmin,nskip,iset1
 99         format(3x,'Number of interaction matrices =',i2,/,
     &           '   min(Nsum) for spectators =',i2,/,
     &           '   Skip first',i2,' sets of G.   iset1 =',i2)
         else
            write(9,"('   Number of NN interaction sets (3=pn,pp,nn):'
     $,i2)")
     $           nsets
         endif

         write(9,101)ki,kf,nf
 101     format(3x,'Initial state ki =',i3,
     &        '   First final state kf =',i3,/,
     &        '   Expectations and transitions will be calculated for'
     &        ,i3,' final states.')
c     ----- Will calculate expectations for the kf-th, (kf+1)-th, through
c     ----- (kf+nf-1)-th states and transitions from the ki-th state to them.
         if(egs.ne.0.0)write(9,*)'  The g.s. energy is:',egs
         write(9,*)'  Number of Lanczos iterations =',nite
      endif
      
      nsetm1 = nsets - 1
      
      ippint=nsets/3
      innint=2*ippint 


c----------------------------------------------------------------------------
c     ----- Call twops to construct all possible two-particle states 
c     ----- in the modell space for fast retrieval of the TBME's:

      call MPI_Barrier(icomm,ierr)

      call twops(nobt,mnn2l,mxn2l,mxnn,ntbme0)

      ntbme = ntbme0

c----------------------------------------------------------------------------
c     ----- Call readtbme to read various TBMEs from TBME.int:

      call MPI_Barrier(icomm,ierr)

      if (tbmeTUD) then
         if (iproc==0) then
            print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            print *,' n_Max,l_Max=',n_Max,l_Max 
         endif
         call TUD_tbme_setup
c***  MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         if (iproc==0) print *,' TUD_tbme_setup called'

         inquire(file=trim(intfile)//'_VNNbin',exist=tbmebin)
         if (tbmebin) then
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.true.)
            if (iproc==0) then
               print *,' read_TUD_tbme_bin called: hbo,lambda=',
     $              hbo_in,lambda_in
               print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            endif
            if (abs(hbo_in-hbomeg)>1.d-4) then
               if (iproc==0) print *,'***error: hbomeg,hbo_in=',
     $              hbomeg,hbo_in
               stop
            endif
            if (N12_max_tud/=mxnn) then
               if (iproc==0) print *,
     $              '*** warning: N12_max,mxnn=',
     $              N12_max_tud,mxnn
            endif
         else
            call read_TUD_tbme(intfile)
            if (iproc==0) print *,' read_TUD_tbme called'
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.false.)
            if (iproc==0) print *,
     $           ' read_TUD_tbme_bin called for Trel, HOrel'
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call set_TUD_H(strcm)
         if (iproc==0) print *,' set_TUD_H called' 
c*** MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
         if (abs(mnop)/=3) then
            call cg_init(N1_max_TUD,N12_max_tud,0)
            if (iproc==0) print *,' cg_init called'
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
      elseif (tbmepn) then
         if (iproc==0) then
            print *,' N1_max,N12_max=',N1_max_pn,N12_max_pn
         endif
         call pn_tbme_setup
c***  MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         if (iproc==0) print *,' pn_tbme_setup called'
         call read_pntbme(intfile)
c***  MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         if (iproc==0) print *,' read_pntbme called'
         call set_pn_H(strcm)
c***  MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI         
         if (iproc==0) print *,' set_pn_H called'
         if (abs(mnop)/=3) then
            call cg_init(N1_max_pn,N12_max_pn,0)
            if (iproc==0) print *,' cg_init called'
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
      else
         call readtbme(intfile,intform,nobt,strcm,iclmb,
     &        itrel,ihrel,ispe,nskip)
      endif
c----------------------------------------------------------------------------
c     ----- Call threejs to set up 3j coefficients for (sp1,sp2)^J:
      call threejs(mxnn)
c----------------------------------------------------------------------------

      if (abs(mnop)==3) !.and.(irest==0.or.irest==4)) 
c     +     call threebodysetup(itrel,mxnn,v3intfile)
     +     call threebodysetup_cJ(itrel,mxnn,v3intfile)
      return
      end


c     This subroutine sets up (J,T,Parity) coupled two-particle states to 
c     be later used for fast retrieval of the two-body G matrix elements. 
c
      subroutine twops(nobt,mnn2l,mxn2l,mxnn,ngs)
      use parameters
      use nodeinfo
      use jsgna
      use spodata
      use tpsb
      allocate(nmebjtp(0:jmx,0:1,0:1),ntps(0:jmx,0:1,0:1))
c
c     ----- Initialization:
      do jt = 0,jmx
         do it = 0,1
            do iparity=0,1
               nmebjtp(jt,it,iparity) = 0
               ntps(jt,it,iparity) = 0
            enddo
         enddo
      enddo

      allocate(itps(nobt,nobt,0:jmx,0:1))
      
      do jt = 0,jmx
         do it = 0,1
            do ia = 1,nobt
               do ib = 1,nobt
                  itps(ia,ib,jt,it) = 0
               enddo
            enddo
         enddo
      enddo
c
      ntps_so_far=0
      ngs_so_far=0
      do jt = 0,jmx
         do it = 0,1
            do iparity=0,1
               njtp = 0         ! Initialize njtp to zero. 
               do 70 ia = 1,nobt
                  na = nnl(ia)
                  la = lx(ia)
                  if(na.lt.mnn2l)goto 70
                  if(na.gt.mxn2l)goto 70
                  do 60 ib = ia,nobt
                     nb = nnl(ib)
                     lb = lx(ib)
                     if(jsgn(iparity+la+lb).eq.-1)goto 60
                     if(nb.lt.mnn2l)goto 60
                     if(nb.gt.mxn2l)goto 60
                     if((na+nb).gt.mxnn)goto 60
                     mnj = (iabs(j2x(ia)-j2x(ib)))/2
                     mxj = (j2x(ia)+j2x(ib))/2
                     if(jt.lt.mnj.or.jt.gt.mxj)goto 60
                     if((ia.eq.ib).and.jsgn(jt+it).eq.1)goto 60
c     ----- This state, which survives all the tests, is to be recorded:
                     njtp = njtp+1
                     itps(ia,ib,jt,it) = njtp
 60               continue
 70            continue
               if(njtp.eq.0)goto 80
c     ----- The number of two-particle states that this (J,T,Parity) has is:
               ntps(jt,it,iparity) = njtp
c     ----- The number of TBMEs appeared before this (J,T,Parity) block is:
               nmebjtp(jt,it,iparity) = ngs_so_far
c     ----- The number of two-particle states has been seen so far:
               ntps_so_far = ntps_so_far+njtp
c     ----- The number of TBMEs has been seen so far:
               ngs_so_far = ngs_so_far+njtp*(njtp+1)/2
 80            continue
            enddo
         enddo
      enddo
      if(iproc.eq.0)then
         write(6,*)'Total number of two-body states:',ntps_so_far
         write(8,*)'Total number of two-body states:',ntps_so_far
         write(9,*)'Total number of two-body states:',ntps_so_far
         write(6,*)'Total number of two-body matrix elements:',
     &        ngs_so_far
         write(8,*)'Total number of two-body matrix elements:',
     &        ngs_so_far
         write(9,*)'Total number of two-body matrix elements:',
     &        ngs_so_far
      endif
      ngs = ngs_so_far
      return
      end


c
c     This subroutine reads all TBME's from TBME.int
c
      subroutine readtbme(intfile,intform,nobt,strcm,iclmb,
     &     itrel,ihrel,ispe,nskip)
      use parameters
      use nodeinfo
      use multig
      use spodata
      use nutbme
      use idxb
      use iters
      use spb
      use tbmes
      use hmes, only: memoryalloc
      include 'mpif.h'
      character(len=120),intent(IN) :: intfile
      dimension skip(5) !,spe(nobt)
      real,allocatable :: gin(:)
      integer :: nobt_in,izz
c
      memoryalloc=memoryalloc+ntbme*4*3
      if (allocated(gful)) deallocate(gful)
      if (allocated(cful)) deallocate(cful)
      if (allocated(rful)) deallocate(rful)
      allocate(gful(ntbme,0:nsetm1))
      allocate(cful(ntbme),rful(ntbme))
      gful=0.d0; cful=0.d0; rful=0.d0
      if (ispe==1) then
         if (allocated(Esp)) deallocate(Esp)
         allocate(Esp(nobt,-1:1))
      endif

      if(iproc>0) goto 100
      allocate(gin(0:nsetm1))
      comfr = 4.0*(197.327)**2/938.9185/((real(nucleons))**2*hbomeg)
      comft = 2.*hbomeg/real(nucleons) 
      comfspe = 0.0 !real(ispe)/real(nucleons-1)
      comfcm1 = strcm*hbomeg/real(nucleons-1)
      comfcm2 = strcm*2.*hbomeg/real(nucleons) 
c     ----- The Coulomb TBMEs are calculated for b = 1 fm. They must be
c     ----- scaled by the actual (1/b):
      comfcl = real(iclmb)*sqrt(938.9185*hbomeg)/197.327

      open(2,file=trim(intfile),status='old')

c      if(ispe.eq.1)then
c         if(major.eq.0)read(2,*)ntbme0,(spe(iorder(ispo)),ispo=1,nobt)
c         if(major>0)read(2,*)ntbme0,(spe(ispo),ispo=1,nobt)
c      else            
      read(2,*)ntbme0
c      do izz = 1,nobt
c         spe(izz) = 0.0
c      enddo
c      endif

      if (ispe/=1.and.ntbme0/=ntbme) goto 200

      if (ispe==1) then
c         if (allocated(Esp)) deallocate(Esp)
c         allocate(Esp(nobt,-1:1))
         Esp=0.0
         read(2,*) nobt_in
         if (nobt_in/=nobt) then
            print *,'***error in Esp: nobt_in,nobt=',nobt_in,nobt
            stop
         endif
         do izz=1,nobt
            read(2,*) ia,Esp(izz,1),Esp(izz,-1)
            if (ia/=izz) then
               print *,'***error in Esp: izz,ia=',izz,ia
               stop
            endif
         end do
         read(2,*) E_core
      endif
         
      do i = 1, ntbme0
         if(nskip.eq.0)then
            if(intform.eq.2)then
               read(2,*)ia,ib,ic,id,jt,it,trel,hrel,coul
     &              ,(gin(k),k=0,nsetm1,1)
            else
               trel=0.0
               hrel=0.0
               coul=0.0
               read(2,*)ia,ib,ic,id,jt,it,(gin(k),k=0,nsetm1,1)
            endif
         else
            if(intform.eq.2)then
               read(2,*)ia,ib,ic,id,jt,it,trel,hrel,coul
     &              ,(skip(j),j=1,nskip),(gin(k),k=0,nsetm1,1)
            else
               trel=0.0
               hrel=0.0
               coul=0.0
               read(2,*)ia,ib,ic,id,jt,it
     &              ,(skip(j),j=1,nskip),(gin(k),k=0,nsetm1,1)
            endif
         endif
         if(major.eq.0) then
            ia = idxobt(ia)
            ib = idxobt(ib)
            ic = idxobt(ic)
            id = idxobt(id)
         endif
         call findindx(ia,ib,ic,id,jt,it,idx,phse)
         trel = trel * phse
         hrel = hrel * phse
         coul = coul * phse
         do jset = 0, nsetm1, 1
            gin(jset) = gin(jset)*phse
         enddo
c     
c     ----- Rrel = Hrel - Trel:
         rful(idx) = comfr*(hrel-trel)
c**********************************************
         if (abs(mnop)==3.and.ihrel==1) gin=0.0
c**********************************************
         do jset = 0, nsetm1, 1
            gful(idx,jset) = gin(jset) 
     &           + comft * (itrel*trel + ihrel*hrel)
         enddo
c----------------------------------------------------------------------------
c     ----- Construct H + strcm*(Hcm-(N0+1.5)hw) 
c     -----         = H + strcm*(H0-Hrel-(N0+1.5)hw)
c     ----- (The second term is non-vanishing only for diagonal ME's.)
c     ----- ( NOTE:  H0 - (N0+1.5)hw = sum_{i<j} [hi+hj-((2N0+3)/A)hw] )
c
c     ----- First construct H + strcm*(H0-(N0+1.5)hw):
         if(ia.eq.ic.and.ib.eq.id.or.ia.eq.id.and.ib.eq.ic)then
            hcmspe = comfcm1*(nnl(ia)+nnl(ib)+3.0
     &           -(2.0*mnn2l+3.0)/real(nucleons))
!     &           + comfspe*(spe(ia)+spe(ib)) ! Single-particle energies
            do jset = 0, nsetm1, 1
               gful(idx,jset) = gful(idx,jset) + hcmspe
            enddo
         endif
c
c     ----- Now construct H + strcm*(H0-(N0+1.5)hw) - strcm*Hrel
c     -----             = H + strcm*(Hcm-(N0+1.5)hw)
c     ----- where (N0+1.5)hw is the lowest c.o.m. energy: N0=mnn2l

         do jset = 0, nsetm1, 1
            gful(idx,jset) = gful(idx,jset) - comfcm2*hrel
         enddo
c----------------------------------------------------------------------------
c     ----- Coulomb:
c     ----- Scale the Coulomb matrix elements to this value of hbomeg
         cful(idx) = comfcl*coul
      enddo
      close(2)
      deallocate(gin)
 100  continue
      if (nproc>1) then
         call MPI_Barrier(icomm,ierr)
         if (ispe==1) then
            call MPI_Bcast(E_core,1,MPI_REAL,0,icomm,ierr)
            call MPI_Bcast(Esp(1,-1),nobt*3,MPI_REAL,0,icomm,ierr)
         endif
         call MPI_Bcast(gful(1,0),(nsetm1+1)*ntbme,
     &        MPI_REAL,0,icomm,ierr)
         call MPI_Bcast(cful,ntbme,MPI_REAL,0,icomm,ierr)
         call MPI_Bcast(rful,ntbme,MPI_REAL,0,icomm,ierr)
      endif
      return
 200  if(ntbme0.ne.ntbme)then
         write(6,*)'  ntbme =',ntbme,'  ntbme0 =',ntbme0
         write(8,*)'  ntbme =',ntbme,'  ntbme0 =',ntbme0
         stop 'ntbme0.ne.ntbme'
      endif
      end


      subroutine findindx(ia,ib,ic,id,jt,it,idx,phse)
      use parameters
      use nodeinfo
      use jsgna
      use spodata
      use nutbme
      use tpsb
      include 'mpif.h'
      iphs = 0
      if(ia.le.ib)then
         i1 = ia
         i2 = ib
      else
         i1 = ib
         i2 = ia
         iphs = iphs + (j2x(i1)+j2x(i2))/2+jt+it
      endif
      if(ic.le.id)then
         i3 = ic
         i4 = id
      else
         i3 = id
         i4 = ic
         iphs = iphs + (j2x(i3)+j2x(i4))/2+jt+it
      endif
      if(i1.gt.i3.or.(i1.eq.i3.and.i2.gt.i4))then
c     ----- Interchange "12" and "34":
         ix = i1
         i1 = i3
         i3 = ix
         ix = i2
         i2 = i4
         i4 = ix
      endif
c     
      iprty = (1-jsgn(lx(i1)+lx(i2)))/2
      ks = ntps(jt,it,iprty)
      k1 = itps(i1,i2,jt,it)
      k2 = itps(i3,i4,jt,it)
      if(k1.le.0.or.k2.le.0) then
         write(6,*)'  Node #',iproc,'  k1 =',k1,'  k2 =',k2
         write(6,*)'  (abcd JT) =', ia,ib,ic,id,jt,it
         call MPI_Abort(icomm,101,ierr)
         stop 'k1 or k2 out of range.'
      endif
      idx = nmebjtp(jt,it,iprty)+(k1-1)*ks-k1*(k1-1)/2+k2
      if(idx.le.0.or.idx.gt.ntbme)then
         write(6,*)'  Node #',iproc,'  idx = ',idx
         write(6,*)'  (abcd JT) =', ia,ib,ic,id,jt,it
         call MPI_Abort(icomm,102,ierr)
         stop 'idx out of range.'
      endif
      phse = real(jsgn(iphs))
      return
      end
     

      subroutine threejs(mxnn)
      use parameters
      use nodeinfo
      use spodata
      use spb
      use spsdata
      use facs
      use the3js

      real(8) fk
c     
      fac(0) = 1.d0
      do k = 1,64
         fac(k)=k*fac(k-1)
      enddo
      sqfac(0) = 1.d0
      do k = 1,64
         fk = real(k)
         sqfac(k)=dsqrt(fk)*sqfac(k-1)
      enddo
      allocate(the3j(mxsps,mxsps,0:jmx))
      do ia = 1, mxsps
         do ib = 1, mxsps
            do jt = 0, jmx
               the3j(ia,ib,jt) = 0.0
            enddo
         enddo
      enddo
c
      do ia = 1, nasps
         nnla = n2l_sp(ia)
         ja = j2_sp(ia)
         ma = m2_sp(ia)
         do 10 ib = 1, nasps
            if((nnla + n2l_sp(ib)).gt.mxnn)goto 10
            jb = j2_sp(ib)
            mb = m2_sp(ib)
            do jab = iabs(ja-jb), min((ja+jb),2*jmx), 2
               if(iabs(ma+mb).le.jab)
     &              the3j(ia,ib,jab/2) = threej(ja,jb,jab,ma,mb)
            enddo
 10      continue
      enddo
c
      do it = 0, 1
         iit = it+it
         do itza = -1,1,2
            do itzb = -1,1,2
               if(iabs(itza+itzb).le.iit)then
                  the3t(it,itza,itzb) = threej(1,1,iit,itza,itzb)
               else
                  the3t(it,itza,itzb) = 0.0
               endif
            enddo
         enddo
      enddo
c     
      return
      end


      function cgcf(j1d,m1d,j2d,m2d,j3d)
      use jsgna
c     Clebsch-Gordan coefficients
c     j's and m's are twice the physical values.
      cgcf=threej(j1d,j2d,j3d,m1d,m2d)*sqrt(float(j3d+1))
      cgcf=jsgn((j2d-j1d-m1d-m2d)/2)*cgcf
      return
      end


      function threej(j1d,j2d,j3d,m1d,m2d)
      use facs
c     Wigner three-j symbols
c     j's and m's are twice the physical values.

      real(8) a
      if(j1d+j2d+j3d.gt.128) stop '2(j1+j2+j3) > 128.'
      threej = 0.0
      if(itri(j1d,j2d,j3d).eq.0) return
      if(iabs(m1d).gt.j1d)return
      if(iabs(m2d).gt.j2d)return
      m3d = -(m1d+m2d)
      if(iabs(m3d).gt.j3d)return
c     
      jjj  = (j1d+j2d+j3d)/2
      jjj3 = jjj - j3d
      j2j3m1 = (j2d-j3d-m1d)/2
      j1j3m2 = (j1d-j3d+m2d)/2
c
      j1pm1 = (j1d+m1d)/2
      j2pm2 = (j2d+m2d)/2
      j3pm3 = (j3d+m3d)/2
      j1mm1 = j1pm1 - m1d
c
      kmin=max0(0,j2j3m1,j1j3m2)
      kmax=min0(jjj3,j1mm1,j2pm2)
      if(kmax.lt.kmin) return
      a = 0.d0
      ib = jjj3 - kmin
      ic = j1mm1 - kmin
      id = j2pm2 - kmin
      ie = kmin - j2j3m1
      ig = kmin - j1j3m2
      do k = kmin, kmax
         a = -a + 1.d0/(fac(k)*fac(ib)*fac(ic)
     &        *fac(id)*fac(ie)*fac(ig))
         ib = ib-1
         ic = ic-1
         id = id-1
         ie = ie+1
         ig = ig+1
      enddo
c
      k = kmax + (j1d-j2d-m3d)/2
      if(k.ne.2*(k/2)) a = -a
c
      a = a*sqfac(jjj-j1d)*sqfac(jjj-j2d)*sqfac(jjj3)
     &     *sqfac(j1pm1)*sqfac(j1mm1)*sqfac(j2pm2)*sqfac(j2pm2-m2d)
     &     *sqfac(j3pm3)*sqfac(j3pm3-m3d)/sqfac(jjj+1)
      threej = sngl(a)
      return
      end


      function itri(na,nb,nc)
      itri=0
      if(((na+nb).ge.nc).and.(iabs(na-nb).le.nc)) itri=1
      return
      end

      subroutine threebodysetup_cJ(itrel,nhom2sp,v3intfile)
c** Setup for a three-body interaction ***
c** Petr Navratil, 30 June 2011, TRIUMF *****
      use parameters
      use nodeinfo
      use spb
      use spsdata
      use spbasis
      use v3b
      use hmes, only: memoryalloc
c      use hamc
      implicit none
      include 'mpif.h'
      integer,intent(IN) :: itrel,nhom2sp
      character(len=120),intent(IN) :: v3intfile
      integer :: J_3,T_3,tot_num_of_3bmatel,ini,fin
      integer :: minz,nspsmin,nhom,ierr
c      integer(8) :: index_v3b_cJ
c      integer(4) :: index_v3b_cJ
c      real(kind(0.d0)),allocatable :: tempspfi(:)
      integer :: ntot_a,n_a,l_a,j2_a,isp1n_a
      integer :: ntot_b,n_b,l_b,j2_b,isp1n_b
      integer :: ntot_c,n_c,l_c,j2_c,isp1n_c
      integer :: ntot_d,n_d,l_d,j2_d,isp1n_d
      integer :: ntot_e,n_e,l_e,j2_e,isp1n_e,isp1n_e_max
      integer :: ntot_f,n_f,l_f,j2_f,isp1n_f,isp1n_f_max
      integer :: ii,j_ab,t_ab,j_de,t_de,i_abc,i_def,i_abcdef
c     $     ,order_abc,order_def
      integer(4) :: ih3,ih3oibuf,iix
      integer :: i,ibuf,pi_i,pi_f
      integer :: Nmin_HO
      logical :: compfile

c      minz=max(0,nneutrns-2)+max(0,nprotons-2)
      minz=Nmin_HO(nprotons)+Nmin_HO(nneutrns)
cc      nspsmin=max(minz-3,0)
      if (iproc==0) print *,' Ntotmin=',minz
      nspsmin=min(Nmin_HO(nprotons-3)+Nmin_HO(nneutrns),
     $     Nmin_HO(nprotons-2)+Nmin_HO(nneutrns-1),
     $     Nmin_HO(nprotons-1)+Nmin_HO(nneutrns-2),
     $     Nmin_HO(nprotons)+Nmin_HO(nneutrns-3))
      if (iproc==0) print *,' nspsmin=',nspsmin
      nhom=nhw-nspsmin   
      if (iproc==0) print *,' Needed N123_max=',nhom
      if (nhom>E123max) then
         if (iproc==0) print *,' E123max too small:',E123max
         stop
      endif
      nhom=E123max
      if (iproc==0) print *,' Used N123_max=',nhom
c*** test
c      nhom=min(nhom,3*(2*n_sp(nasps)+l_sp(nasps)))
c      nhom=max(nhom,nhom2sp) ! good only for truncating 6Li,4He inter. use rather nhom=max(nhom,N123max)
c*** test

      N1_max=nshll-1
      if (E1max/=N1_max.or.E12max/=nhom2sp) then
         print *,' basis size and 3N file size do not match:',
     $        N1_max,nhom2sp,E1max,E12max
         stop
      endif
      call sp1nbas_cJ(N1_max,9)

      memoryalloc=memoryalloc+4*isp1ntot**3
      allocate(index_abc(isp1ntot,isp1ntot,isp1ntot))
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               ii=ii+1
               index_abc(isp1n_a,isp1n_b,isp1n_c)=ii
            end do
         end do
      end do
      dim_abc=ii
c     print *,' dim_abc=',dim_abc
      memoryalloc=memoryalloc+4*dim_abc**2
      allocate(index_abcdef(dim_abc,dim_abc))
      index_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
c               order_abc=(isp1n_a-1)*((isp1ntot-1)*(isp1ntot+1)+1)
c     $              +(isp1n_b-1)*isp1ntot+isp1n_c-1

c               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
c                  do isp1n_e=1,isp1n_d
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
c                     do isp1n_f=1,isp1n_e
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
c                        order_def=(isp1n_d-1)
c     $                       *((isp1ntot-1)*(isp1ntot+1)+1)
c     $                       +(isp1n_e-1)*isp1ntot+isp1n_f-1
c                        if (order_def>order_abc) then
                        if (i_def>i_abc) then
                           print *,'a,b,c=',isp1n_a,isp1n_b,isp1n_c
                           print *,'e,d,f=',isp1n_e,isp1n_d,isp1n_f
c                           print *,' order_abc=',order_abc
c                           print *,' order_def=',order_def
                           print *,' i_abc=',i_abc
                           print *,' i_def=',i_def
                           stop
                        endif
cc                        if (i_def<i_abc) cycle
                        ii=ii+1
                        index_abcdef(i_abc,i_def)=ii
                     end do
                  end do
               end do
            end do
         end do
      end do
      dim_abcdef=ii
c     print *,' dim_abcdef=',dim_abcdef
      memoryalloc=memoryalloc+4*dim_abcdef
      allocate(start_abcdef(dim_abcdef))
      start_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
c               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
c                  do isp1n_e=1,isp1n_d
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
c                     do isp1n_f=1,isp1n_e
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
cc                        if (i_def<i_abc) cycle
                        i_abcdef=index_abcdef(i_abc,i_def)
                        start_abcdef(i_abcdef)=ii+1
                        do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                           do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                              do J_3=max(abs(2*j_ab-j2_c),
     $                             abs(2*j_de-j2_f)),
     $                             min(2*j_ab+j2_c,2*j_de+j2_f),2
                                 do t_ab=0,1
c                                    if (isp1n_a==isp1n_b.and.
c     $                                   mod(j_ab+t_ab,2)==0) cycle
                                    do t_de=0,1
c                                       if (isp1n_d==isp1n_e.and.
c     $                                      mod(j_de+t_de,2)==0) cycle

                                       do T_3=max(abs(2*t_ab-1),
     $                                      abs(2*t_de-1)),
     $                                      min(2*t_ab+1,2*t_de+1),2
                                          ii=ii+1
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      tot_num_of_3bmatel=ii
      if (iproc==0) then
         print *,' tot_num_of_3bmatel=',tot_num_of_3bmatel
c         write(2,1000) tot_num_of_3bmatel
c 1000    format(/,' Number of 3N matrix elements in threebodysetup_cJ:',
c     $        i8)
      endif
      memoryalloc=memoryalloc+tot_num_of_3bmatel*4
      allocate(v3b_cJ(tot_num_of_3bmatel))
      if (iproc==0) then

         inquire(file=trim(v3intfile)//'_comp',exist=compfile)
         if (compfile) then

            open(27,file=trim(v3intfile)//'_comp',status='old',
     $           form='unformatted',action='read')
            write(9,*)
            write(9,'(a)')
     $           ' 3N interaction file '//trim(v3intfile)//'_comp used'

            read(27) (v3b_cJ(ii),ii=1,tot_num_of_3bmatel)
            close(27)

         else

            open(27,file=trim(v3intfile),status='old',
     $           form='unformatted',action='read')
            ii=0
         do isp1n_a=1,isp1ntot
            n_a=isp1n_cJ(1,isp1n_a)
            l_a=isp1n_cJ(2,isp1n_a)
            j2_a=isp1n_cJ(3,isp1n_a)
            ntot_a=2*n_a+l_a
            if (ntot_a>nhom) cycle
            do isp1n_b=1,isp1n_a
               n_b=isp1n_cJ(1,isp1n_b)
               l_b=isp1n_cJ(2,isp1n_b)
               j2_b=isp1n_cJ(3,isp1n_b)
               ntot_b=2*n_b+l_b
               if (ntot_a+ntot_b>nhom2sp) cycle
               do isp1n_c=1,isp1n_b
                  n_c=isp1n_cJ(1,isp1n_c)
                  l_c=isp1n_cJ(2,isp1n_c)
                  j2_c=isp1n_cJ(3,isp1n_c)
                  ntot_c=2*n_c+l_c
                  if (ntot_a+ntot_c>nhom2sp) cycle
                  if (ntot_b+ntot_c>nhom2sp) cycle
                  if (ntot_a+ntot_b+ntot_c>nhom) cycle
                  i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
                  pi_i=(-1)**(l_a+l_b+l_c)
c                  do isp1n_d=1,isp1ntot
                  do isp1n_d=1,isp1n_a
                     n_d=isp1n_cJ(1,isp1n_d)
                     l_d=isp1n_cJ(2,isp1n_d)
                     j2_d=isp1n_cJ(3,isp1n_d)
                     ntot_d=2*n_d+l_d
                     if (ntot_d>nhom) cycle
                     if (isp1n_d==isp1n_a) then
                        isp1n_e_max=isp1n_b
                     else
                        isp1n_e_max=isp1n_d
                     endif
c                     do isp1n_e=1,isp1n_d
                     do isp1n_e=1,isp1n_e_max
                        n_e=isp1n_cJ(1,isp1n_e)
                        l_e=isp1n_cJ(2,isp1n_e)
                        j2_e=isp1n_cJ(3,isp1n_e)
                        ntot_e=2*n_e+l_e
                        if (ntot_d+ntot_e>nhom2sp) cycle
                        if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                           isp1n_f_max=isp1n_c
                        else
                           isp1n_f_max=isp1n_e
                        endif
c                     do isp1n_f=1,isp1n_e
                        do isp1n_f=1,isp1n_f_max
                           n_f=isp1n_cJ(1,isp1n_f)
                           l_f=isp1n_cJ(2,isp1n_f)
                           j2_f=isp1n_cJ(3,isp1n_f)
                           ntot_f=2*n_f+l_f
                           if (ntot_d+ntot_f>nhom2sp) cycle
                           if (ntot_e+ntot_f>nhom2sp) cycle
                           if (ntot_d+ntot_e+ntot_f>nhom) cycle
                           pi_f=(-1)**(l_d+l_e+l_f)
                           if (pi_i/=pi_f) cycle
                           i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
cc                           if (i_def<i_abc) cycle
c     i_abcdef=index_abcdef(i_abc,i_def)
c                        start_abcdef(i_abcdef)=ii+1
                           do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                              do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                                 do J_3=max(abs(2*j_ab-j2_c),
     $                                abs(2*j_de-j2_f)),
     $                                min(2*j_ab+j2_c,2*j_de+j2_f),2
                                    do t_ab=0,1
c                                       if (isp1n_a==isp1n_b.and.
c     $                                      mod(j_ab+t_ab,2)==0) cycle
                                       do t_de=0,1
c                                          if (isp1n_d==isp1n_e.and.
c     $                                      mod(j_de+t_de,2)==0) cycle
                                          
                                          do T_3=max(abs(2*t_ab-1),
     $                                         abs(2*t_de-1)),
     $                                         min(2*t_ab+1,2*t_de+1),2
                                             ii=ii+1
                                             read(27) v3b_cJ(ii)
c     print *, isp1n_a,isp1n_b,
c     $                                         isp1n_c,isp1n_d,isp1n_e,
c     $                                         isp1n_f,j_ab,t_ab,
c     $                                         j_de,t_de,J_3,T_3,
c     $                                         v3b_cJ(ii)
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
            close(27)
         endif
      endif
      if (iproc==0) print *,' v3b_cJ read'
      call MPI_Barrier(icomm,ierr)
      ibuf=3000000
      if (tot_num_of_3bmatel<=ibuf) then
         i=tot_num_of_3bmatel
         call MPI_Bcast(v3b_cJ(1),i,MPI_REAL,0,icomm,ierr) 
      else
         ih3oibuf=tot_num_of_3bmatel/ibuf
         iix=1
         do i=1,ih3oibuf
            call MPI_Bcast(v3b_cJ(iix),ibuf,MPI_REAL,0,icomm,ierr) 
            if (iproc==0) print *,' Broadcasted iix=',iix
            iix=iix+ibuf
         end do
         i=mod(tot_num_of_3bmatel,int(ibuf,kind(4)))
         if (i/=0) then
            call MPI_Bcast(v3b_cJ(iix),i,MPI_REAL,0,icomm,ierr)
            if (iproc==0) print *,' Broadcasted iix=',iix
         endif
      endif
      if (iproc==0) print *,' H3 broadcasted'
      call cg_init(N1_max,nhom2sp,nhom)
      if (iproc==0) print *,' cg_init called'
      end

      subroutine sp1nbas_cJ(nhom,iunitout)
      use parameters
      use spsdata
      use spbasis
      use v3b
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer,intent(IN) :: nhom,iunitout
      integer :: n,l,j2
      integer :: ii,ntot
      ii=0
      do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
            end do
         end do
      end do
      isp1ntot=ii
      if (allocated(isp1n_cJ)) deallocate(isp1n_cJ)
      allocate(isp1n_cJ(3,isp1ntot))
      if (iproc==0) then
         write(iunitout,1000) isp1ntot
c      write(6,1000) isp1ntot
 1000    format(' Number of nlj states:',i5)
      endif
      ii=0
      do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
               isp1n_cJ(1,ii)=n
               isp1n_cJ(2,ii)=l
               isp1n_cJ(3,ii)=j2
c               write(2,1100) ii,n,l,j2
c               write(6,1100) ii,n,l,j2
c 1100          format(' #',i5,'  n=',i3,'  l=',i3,
c     +              '  j=',i3,'/2')
            end do
         end do
      end do
      allocate(nlj_orb(mxsps2))
      nlj_orb=0
      do ii=1,mxsps2
         do l=1,isp1ntot
            if (n_sp(ii)==isp1n_cJ(1,l).and.l_sp(ii)==isp1n_cJ(2,l)
     $           .and.j2_sp(ii)==isp1n_cJ(3,l)) then
               nlj_orb(ii)=l
               exit
            endif
         end do
      end do
      end

      subroutine cg_init(N1max,N12max,N123max)
      use v3b
      use spbasis
      use TUD_tbme, only: cgt12_tud
      implicit none
      integer,intent(IN) :: N1max,N12max,N123max
      integer :: mt1,mt2,t12,mt12,mt3,T3
      integer :: j1,m1,j2,m2,j12,m12,j3,m3,j123,j1max,j12max,j123max
      real(kind(0.d0)) :: clebd,cl
c      allocate(cgj12())
      cgt12=0.d0
      cgt123=0.d0
      cgt12_tud=0.d0
      do mt1=-1,1,2
cc         print *, ' mt1,mt1/2=',mt1,mt1/2
         do mt2=-1,1,2
            do t12=abs(mt1+mt2),2,2
               cl=clebd(1,mt1,1,mt2,t12,mt1+mt2)
               cgt12(t12/2,(mt1+1)/2,(mt2+1)/2)=cl
               cgt12_tud(t12/2,(mt1+mt2)/2,(mt1+1)/2,(mt2+1)/2)=cl
            end do
         end do
      end do
      do t12=0,2,2
         do mt12=-t12,t12,2
            do mt3=-1,1,2
               do T3=max(abs(t12-1),abs(mt12+mt3)),t12+1,2
                  cgt123(t12/2,T3/2,mt12/2,(mt3+1)/2)=
     $                 clebd(t12,mt12,1,mt3,T3,mt12+mt3)
               end do
            end do
         end do
      end do

      j1max=2*N1max+1
      j12max=2*(N12max+1)
      j123max=2*N123max+3

      allocate(cgj12(0:j12max/2,0:j1max/2,0:j1max,0:j1max/2,0:j1max))
      allocate(cgj123(0:j123max/2,0:j12max/2,-j12max/2:j12max/2,
     $     0:j1max/2,0:j1max))
      cgj12=0.d0
      cgj123=0.d0
      do j1=1,j1max,2
         do m1=-j1,j1,2
            do j2=1,j1max,2
               do m2=-j2,j2,2
                  do j12=max(abs(j1-j2),abs(m1+m2)),min(j1+j2,j12max),2
                     cgj12(j12/2,j1/2,(j1+m1)/2,j2/2,(j2+m2)/2)=
     $                    clebd(j1,m1,j2,m2,j12,m1+m2)
                  end do
               end do
            end do
         end do
      end do

      do j12=0,j12max,2
         do m12=-j12,j12,2
            do j3=1,j1max,2
               do m3=-j3,j3,2
                  do j123=max(abs(j12-j3),abs(m12+m3)),
     $                 min(j12+j3,j123max),2
                     cgj123(j123/2,j12/2,m12/2,j3/2,(j3+m3)/2)=
     $                    clebd(j12,m12,j3,m3,j123,m12+m3)
                  end do
               end do
            end do
         end do
      end do
      end

      subroutine TUD_tbme_setup
      use paramdef, only: nlj_st_TUD
      use TUD_tbme
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer :: n,l,j2,N1,ii,ia,ib,ic,id,na,la,ja,nb,lb,jb
      integer :: nc,lc,jc,nd,ld,jd,i_ab,i_cd,i_abcd,id_max,J12,T12,MT12
      integer :: N1_max,N12_max
      N1_max=N1_max_tud
      N12_max=N12_max_tud
      if (allocated(nlj_st_TUD)) deallocate(nlj_st_TUD)
      allocate(nlj_st_TUD(0:N1_max/2))
      do n=0,N1_max/2
         allocate(nlj_st_TUD(n)%l(0:N1_max-2*n))
         do l=0,N1_max-2*n
            allocate(nlj_st_TUD(n)%l(l)%j((2*l-1)/2:(2*l+1)/2))
            nlj_st_TUD(n)%l(l)%j(:)=-1
         end do
      end do
      ii=0
      do N1=0,N1_max
         do l=mod(N1,2),N1,2
            if (l>l_Max) cycle
            n=(N1-l)/2
            if (n>n_Max) cycle
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
            end do   
         end do
      end do
      ist1dim_TUD=ii
      if (iproc==0) print *,' ist1dim_TUD=',ist1dim_TUD 
      if (allocated(ist1_TUD)) deallocate(ist1_TUD)
      allocate(ist1_TUD(3,ist1dim_TUD))
      ii=0
      do N1=0,N1_max
         do l=mod(N1,2),N1,2
            if (l>l_Max) cycle
            n=(N1-l)/2
            if (n>n_Max) cycle
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
               ist1_TUD(1,ii)=n
               ist1_TUD(2,ii)=l
               ist1_TUD(3,ii)=j2
               nlj_st_TUD(n)%l(l)%j(j2/2)=ii
c               if (iproc==0) print *,' nlj_st_TUD=',ii,' n,l,j2:',n,l,j2
c               write(iunitout,3519) ii,2*n+l,n,l,j2
c 3519          format(' level#',i4,'    N=',i2,'   nr=',i2,
c     +              '   l=',i2,'  j=',i2,'/2')
            end do   
         end do
      end do
      
      if (allocated(index_ab)) deallocate(index_ab)
      allocate(index_ab(ist1dim_TUD,ist1dim_TUD))
      index_ab=-1
      ii=0
      do ia=1,ist1dim_TUD
         na=ist1_TUD(1,ia)
         la=ist1_TUD(2,ia)
         do ib=1,ia 
            nb=ist1_TUD(1,ib)
            lb=ist1_TUD(2,ib)
            if (2*na+la+2*nb+lb>N12_max) cycle
            ii=ii+1
            index_ab(ia,ib)=ii
         end do
      end do
      dim_ab_TUD=ii
      if (iproc==0) print *,' dim_ab_TUD=',dim_ab_TUD
      if (allocated(index_abcd)) deallocate(index_abcd)
      allocate(index_abcd(dim_ab_TUD,dim_ab_TUD))
      index_abcd=-1
      tot_num_of_NNmatel_TUD=0
      ii=0
      do ia=1,ist1dim_TUD
         na=ist1_TUD(1,ia)
         la=ist1_TUD(2,ia)
         ja=ist1_TUD(3,ia)
         do ib=1,ia 
            nb=ist1_TUD(1,ib)
            lb=ist1_TUD(2,ib)
            jb=ist1_TUD(3,ib)
            if (2*na+la+2*nb+lb>N12_max) cycle
            i_ab=index_ab(ia,ib)
            do ic=1,ia
               nc=ist1_TUD(1,ic)
               lc=ist1_TUD(2,ic)
               jc=ist1_TUD(3,ic)
               if (ic==ia) then
                  id_max=ib
               else
                  id_max=ic
               endif
               do id=1,id_max 
                  nd=ist1_TUD(1,id)
                  ld=ist1_TUD(2,id)
                  jd=ist1_TUD(3,id)
                  if (2*nc+lc+2*nd+ld>N12_max) cycle
                  if (mod(la+lb+lc+ld,2)==1) cycle
                  i_cd=index_ab(ic,id)
                  ii=ii+1
                  index_abcd(i_ab,i_cd)=tot_num_of_NNmatel_TUD+1
                  do J12=max(abs(ja-jb),abs(jc-jd))/2,min(ja+jb,jc+jd)/2
                     do T12=0,1
                        do MT12=-T12,T12
                           tot_num_of_NNmatel_TUD=
     $                          tot_num_of_NNmatel_TUD+1
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      dim_abcd_TUD=ii
      if (iproc==0) then
         print *,' dim_abcd_TUD=',dim_abcd_TUD
         print *,' tot_num_of_NNmatel_TUD=',tot_num_of_NNmatel_TUD
      endif
      end 

      subroutine read_TUD_tbme(filename)
      use paramdef, only: nlj_st_TUD
      use TUD_tbme
      use nodeinfo
      use hmes, only: memoryalloc
      implicit none
      include 'mpif.h'
      character(len=120), intent(IN) :: filename
      integer :: tot_num_of_el,tot_num_of_lines
      character(len=200) :: line
      integer :: i,i1,irec,j
      real(kind(0.0)) :: v2,vx(10)
      integer :: ibuf,iix,ih3oibuf,ierr

      memoryalloc=memoryalloc+tot_num_of_NNmatel_TUD*8
      if (allocated(V_NN_TUD)) deallocate(V_NN_TUD)
      allocate(V_NN_TUD(tot_num_of_NNmatel_TUD))
      V_NN_TUD=0.d0
      if (no2bv3n) then
         if (allocated(V_3N_no1b_TUD)) deallocate(V_3N_no1b_TUD)
         allocate(V_3N_no1b_TUD((ist1dim_TUD*(ist1dim_TUD+1))/2,0:1))
         V_3N_no1b_TUD=0.d0
         memoryalloc=memoryalloc+(ist1dim_TUD*(ist1dim_TUD+1)+1)*8
      endif
      
      if (iproc==0) then
         open(11,file=trim(filename),status='old',form ='formatted',
     $        action='read')

         read(11,*,end=11111)
         i=0
         do
            read(11,*,end=11111) v2
            i=i+1
         end do
11111    continue
         tot_num_of_lines=i
         tot_num_of_el=10*(i-1)

         if (tot_num_of_el>tot_num_of_NNmatel_TUD) then
            print *,'*** error: tot_num_of_el,tot_num_of_NNmatel_TUD=',
     $           tot_num_of_el,tot_num_of_NNmatel_TUD
            stop
         endif

         rewind(11)
         read(11,*)
         do i=1,tot_num_of_lines-1
            read(11,*) (V_NN_TUD(10*(i-1)+j),j=1,10)
         end do
         read(11,'(a)') line
         irec=0
         line=adjustl(line)
         do
            if (len_trim(line)==0) exit 
            irec=irec+1
            i1=index(line,' ')
            line=adjustl(line(i1:))
         end do
         backspace(11)
         read(11,*,err=3353,end=3353) (vx(j),j=1,irec)
         close(11)

         if (tot_num_of_el+irec/=tot_num_of_NNmatel_TUD) then
            print *,
     $'*** error: tot_num_of_el,irec,tot_num_of_NNmatel_TUD='
     $           ,tot_num_of_el,irec,tot_num_of_NNmatel_TUD
            stop
         endif      
         V_NN_TUD(tot_num_of_el+1:tot_num_of_NNmatel_TUD)=vx(1:irec)
         if (no2bv3n) then
            call read_TUD_nov3n
         endif
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=3000000
      if (tot_num_of_NNmatel_TUD<=ibuf) then
         i=tot_num_of_NNmatel_TUD
         call MPI_Bcast(V_NN_TUD(1),i,MPI_REAL8,0,icomm,ierr) 
      else
         ih3oibuf=tot_num_of_NNmatel_TUD/ibuf
         iix=1
         do i=1,ih3oibuf
            call MPI_Bcast(V_NN_TUD(iix),ibuf,MPI_REAL8,0,icomm,ierr) 
            if (iproc==0) print *,' Broadcasted iix=',iix
            iix=iix+ibuf
         end do
         i=mod(tot_num_of_NNmatel_TUD,int(ibuf,kind(4)))
         if (i/=0) then
            call MPI_Bcast(V_NN_TUD(iix),i,MPI_REAL8,0,icomm,ierr)
            if (iproc==0) print *,' Broadcasted iix=',iix
         endif
      endif
      if (iproc==0) print *,' V_NN_TUD broadcasted' 
      if (no2bv3n) then
         call MPI_Bcast(V_3N_no0b,1,MPI_REAL8,0,icomm,ierr) 
         i=ist1dim_TUD*(ist1dim_TUD+1)
         call MPI_Bcast(V_3N_no1b_TUD(1,0),i,MPI_REAL8,0,icomm,ierr) 
         if (iproc==0) print *,' V_3N_no1b_TUD broadcasted' 
      endif
      return
 3353 continue
      print *,' error in reading in read_TUD_tbme'
      stop
      end subroutine read_TUD_tbme

      subroutine read_TUD_tbme_bin(filename,hbo,lambda,NNbin)
      use TUD_tbme
      use hmes, only: memoryalloc
      use nodeinfo
      implicit none
      include 'mpif.h'
      character(len=120), intent(IN) :: filename
      integer :: N1_max,N12_max
      real(kind(0.d0)),intent(OUT) :: hbo,lambda
      logical,intent(IN) :: NNbin
      integer :: tot_num_of_NNmatel_TUD_in,N1_max_in,N12_max_in
      integer :: i,ibuf,iix,ih3oibuf,ierr
      memoryalloc=memoryalloc+tot_num_of_NNmatel_TUD*8*3
      hbo=-1.d0
      lambda=-1.d0
      if (NNbin) then
         if (iproc==0) then
            open(22,file=trim(filename)//'_VNNbin',status='old',
     $           form ='unformatted',action='read')

            write(9,"(/,' NN interaction file used: ',a)")
     $           trim(filename)//'_VNNbin'
            
            read(22) tot_num_of_NNmatel_TUD,N1_max,N12_max,hbo,lambda
            if (N1_max/=N1_max_TUD.or.N12_max/=N12_max_TUD) then
               print *,
     $'***error in read_TUD_tbme_bin:
     $N1_max,N1_max_TUD,N12_max,N12_max_TUD=',
     $              N1_max,N1_max_TUD,N12_max,N12_max_TUD
               stop
            endif
         endif
         call MPI_Bcast(tot_num_of_NNmatel_TUD,1,MPI_INTEGER,0,icomm,
     $        ierr) 
         call MPI_Bcast(N1_max,1,MPI_INTEGER,0,icomm,ierr) 
         call MPI_Bcast(N12_max,1,MPI_INTEGER,0,icomm,ierr) 
         call MPI_Bcast(hbo,1,MPI_REAL8,0,icomm,ierr) 
         call MPI_Bcast(lambda,1,MPI_REAL8,0,icomm,ierr)
         if (allocated(V_NN_TUD)) deallocate(V_NN_TUD)
         allocate(V_NN_TUD(tot_num_of_NNmatel_TUD))
         V_NN_TUD=0.d0
         if (no2bv3n) then
            if (allocated(V_3N_no1b_TUD)) deallocate(V_3N_no1b_TUD)
            allocate(V_3N_no1b_TUD((ist1dim_TUD*(ist1dim_TUD+1))/2,0:1))
            V_3N_no1b_TUD=0.d0
            memoryalloc=memoryalloc+((ist1dim_TUD*(ist1dim_TUD+1))+1)
     $           *8
         endif
         if (iproc==0) then
            read(22) (V_NN_TUD(i),i=1,tot_num_of_NNmatel_TUD)
            close(22)
            if (no2bv3n) then
               call read_TUD_nov3n
            endif
         endif
         call MPI_Barrier(icomm,ierr)
         ibuf=3000000
         if (tot_num_of_NNmatel_TUD<=ibuf) then
            i=tot_num_of_NNmatel_TUD
            call MPI_Bcast(V_NN_TUD(1),i,MPI_REAL8,0,icomm,ierr) 
         else
            ih3oibuf=tot_num_of_NNmatel_TUD/ibuf
            iix=1
            do i=1,ih3oibuf
               call MPI_Bcast(V_NN_TUD(iix),ibuf,MPI_REAL8,0,icomm,ierr) 
               if (iproc==0) print *,' Broadcasted iix=',iix
               iix=iix+ibuf
            end do
            i=mod(tot_num_of_NNmatel_TUD,int(ibuf,kind(4)))
            if (i/=0) then
               call MPI_Bcast(V_NN_TUD(iix),i,MPI_REAL8,0,icomm,ierr)
               if (iproc==0) print *,' Broadcasted iix=',iix
            endif
         endif
         if (iproc==0) print *,' V_NN_TUD broadcasted'
         if (no2bv3n) then
            call MPI_Bcast(V_3N_no0b,1,MPI_REAL8,0,icomm,ierr) 
            i=ist1dim_TUD*(ist1dim_TUD+1)
            call MPI_Bcast(V_3N_no1b_TUD(1,0),i,MPI_REAL8,0,icomm,ierr) 
            if (iproc==0) print *,' V_3N_no1b_TUD broadcasted' 
         endif
      endif
      if (Trelsave) then
         if (iproc==0) then
            open(22,file=trim(filename)//'_Trelbin',status='old',
     $           form ='unformatted',action='read')

            write(9,"(' Trel file used: ',a)")
     $           trim(filename)//'_Trelbin'
            
            read(22) tot_num_of_NNmatel_TUD_in,N1_max_in,N12_max_in
            if (tot_num_of_NNmatel_TUD_in/=tot_num_of_NNmatel_TUD
     $           .or.N1_max_in/=N1_max_TUD.or.N12_max_in/=N12_max_TUD)
     $           then
               print *,
     $'*** error: tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD='
     $              ,tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD
               print *,
     $'*** error: Trel N1(2)_max_in,N1(2)_max=',N1_max_in,N1_max_TUD,
     $              N12_max_in,N12_max_TUD        
               stop
            endif     
         endif
         if (allocated(Trel_TUD)) deallocate(Trel_TUD)
         allocate(Trel_TUD(tot_num_of_NNmatel_TUD))
         Trel_TUD=0.d0
         if (iproc==0) then
            read(22) (Trel_TUD(i),i=1,tot_num_of_NNmatel_TUD)
            close(22)
         endif
         call MPI_Barrier(icomm,ierr)
         ibuf=3000000
         if (tot_num_of_NNmatel_TUD<=ibuf) then
            i=tot_num_of_NNmatel_TUD
            call MPI_Bcast(Trel_TUD(1),i,MPI_REAL8,0,icomm,ierr) 
         else
            ih3oibuf=tot_num_of_NNmatel_TUD/ibuf
            iix=1
            do i=1,ih3oibuf
               call MPI_Bcast(Trel_TUD(iix),ibuf,MPI_REAL8,0,icomm,ierr) 
               if (iproc==0) print *,' Broadcasted iix=',iix
               iix=iix+ibuf
            end do
            i=mod(tot_num_of_NNmatel_TUD,int(ibuf,kind(4)))
            if (i/=0) then
               call MPI_Bcast(Trel_TUD(iix),i,MPI_REAL8,0,icomm,ierr)
               if (iproc==0) print *,' Broadcasted iix=',iix
            endif
         endif
         if (iproc==0) print *,' Trel_TUD broadcasted'    
         if (iproc==0) then
            open(22,file=trim(filename)//'_HOrelbin',status='old',
     $           form ='unformatted',action='read')
            
            write(9,"(' HOrel file used: ',a)")
     $           trim(filename)//'_HOrelbin'
            
            read(22) tot_num_of_NNmatel_TUD_in,N1_max_in,N12_max_in
            if (tot_num_of_NNmatel_TUD_in/=tot_num_of_NNmatel_TUD
     $           .or.N1_max_in/=N1_max_TUD.or.N12_max_in/=N12_max_TUD)
     $           then
               print *,
     $'*** error: tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD='
     $              ,tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD
               print *,
     $'*** error: HOrel N1(2)_max_in,N1(2)_max=',N1_max_in,N1_max_TUD,
     $              N12_max_in,N12_max_TUD        
               stop
            endif     
         endif
         if (allocated(HOrel_TUD)) deallocate(HOrel_TUD)
         allocate(HOrel_TUD(tot_num_of_NNmatel_TUD))
         HOrel_TUD=0.d0
         if (iproc==0) then
            read(22) (HOrel_TUD(i),i=1,tot_num_of_NNmatel_TUD)
            close(22)
         endif
         call MPI_Barrier(icomm,ierr)
         ibuf=3000000
         if (tot_num_of_NNmatel_TUD<=ibuf) then
            i=tot_num_of_NNmatel_TUD
            call MPI_Bcast(HOrel_TUD(1),i,MPI_REAL8,0,icomm,ierr) 
         else
            ih3oibuf=tot_num_of_NNmatel_TUD/ibuf
            iix=1
            do i=1,ih3oibuf
               call MPI_Bcast(HOrel_TUD(iix),ibuf,MPI_REAL8,0,
     $              icomm,ierr) 
               if (iproc==0) print *,' Broadcasted iix=',iix
               iix=iix+ibuf
            end do
            i=mod(tot_num_of_NNmatel_TUD,int(ibuf,kind(4)))
            if (i/=0) then
               call MPI_Bcast(HOrel_TUD(iix),i,MPI_REAL8,0,icomm,ierr)
               if (iproc==0) print *,' Broadcasted iix=',iix
            endif
         endif
         if (iproc==0) print *,' HOrel_TUD broadcasted'    
      endif
      end subroutine read_TUD_tbme_bin

      subroutine read_TUD_nov3n
      use paramdef, only: nlj_st_TUD
      use TUD_tbme
      implicit none
      integer :: tot_num_of_el,tot_num_of_lines
      character(len=200) :: line
      integer :: i,i1,irec,j
      integer :: tot_num_of_NNmatel_TUD_in,N1_max_in,N12_max_in,
     $     tot_num_of_1bmatel_in,nljt,nnljt,ia,ib,mta,mtb,i_ab
      real(kind(0.0)) :: v2,vx(10)
      logical :: binfile
      real(kind(0.d0)),allocatable :: temp_TUD(:)

      inquire(file=trim(v3intfile)//'_no2bV3Nbin',exist=binfile)

      if (binfile) then
         open(22,file=trim(v3intfile)//'_no2bV3Nbin',status='old',
     $        form ='unformatted',action='read')

          write(9,"(' NO2B 3N file used (2-body): ',a)")
     $           trim(v3intfile)//'_no2bV3Nbin'
         
         read(22) tot_num_of_NNmatel_TUD_in,N1_max_in,N12_max_in
         if (tot_num_of_NNmatel_TUD_in/=tot_num_of_NNmatel_TUD
     $        .or.N1_max_in/=N1_max_TUD.or.N12_max_in/=N12_max_TUD) then
            print *,
     $'*** error: tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD='
     $,tot_num_of_NNmatel_TUD_in,tot_num_of_NNmatel_TUD
            print *,
     $'*** error in no2bv3n: N1(2)_max_in,N1(2)_max=',
     $           N1_max_in,N1_max_TUD,
     $           N12_max_in,N12_max_TUD       
            stop
         endif     
         if (allocated(temp_TUD)) deallocate(temp_TUD)
         allocate(temp_TUD(tot_num_of_NNmatel_TUD))
         temp_TUD=0.d0
         read(22) (temp_TUD(i),i=1,tot_num_of_NNmatel_TUD)
         close(22)
         do j=1,tot_num_of_NNmatel_TUD
             V_NN_TUD(j)=V_NN_TUD(j)+temp_TUD(j)
         end do
         deallocate(temp_TUD)

         open(22,file=trim(v3intfile)//'_no1bV3Nbin',status='old',
     $        form ='unformatted',action='read')

         write(9,"(' NO2B 3N file used (1-body): ',a)")
     $        trim(v3intfile)//'_no1bV3Nbin'
         
         read(22) N1_max_in
         read(22) V_3N_no0b
         read(22) tot_num_of_1bmatel_in
         if (tot_num_of_1bmatel_in/=ist1dim_TUD*(ist1dim_TUD+1)/2
     $        .or.N1_max_in/=N1_max_TUD) then
            print *,
     $'*** error: tot_num_of_1bmatel_in,ist1dim_TUD*(ist1dim_TUD+1)/2='
     $,tot_num_of_1bmatel_in,ist1dim_TUD*(ist1dim_TUD+1)
            print *,
     $'*** error in no2bv3n: N1_max_in,N1_max=',N1_max_in,N1_max_TUD
            stop
         endif
         read(22) ((V_3N_no1b_TUD(i,j),i=1,tot_num_of_1bmatel_in),
     $        j=0,1)
cc         read(22) (V_3N_no1b_TUD(i),i=1,tot_num_of_1bmatel_in)
         close(22)

      else
         open(11,file=trim(v3intfile)//'_no2bV3N',status='old',
     $        form ='formatted',action='read')

         write(9,"(' NO2B 3N file used (2-body): ',a)")
     $        trim(v3intfile)//'_no2bV3N'
         
         read(11,*,end=11111)
         i=0
         do
            read(11,*,end=11111) v2
            i=i+1
         end do
11111    continue
         tot_num_of_lines=i
         tot_num_of_el=10*(i-1)
         if (tot_num_of_el>tot_num_of_NNmatel_TUD) then
            print *,
     $'*** error in no2bv3n: tot_num_of_el,tot_num_of_NNmatel_TUD=',
     $           tot_num_of_el,tot_num_of_NNmatel_TUD
            stop
         endif
         rewind(11)
         read(11,*)
         do i=1,tot_num_of_lines-1
            read(11,*) (vx(j),j=1,10)
            do j=1,10
               V_NN_TUD(10*(i-1)+j)=V_NN_TUD(10*(i-1)+j)+vx(j)
            end do
         end do
         read(11,'(a)') line
         irec=0
         line=adjustl(line)
         do
            if (len_trim(line)==0) exit 
            irec=irec+1
            i1=index(line,' ')
            line=adjustl(line(i1:))
         end do
         backspace(11)
         read(11,*,err=3353,end=3353) (vx(j),j=1,irec)
         close(11)

         if (tot_num_of_el+irec/=tot_num_of_NNmatel_TUD) then
            print *,
     $'*** error in no2bv3n: tot_num_of_el,irec,tot_num_of_NNmatel_TUD='
     $           ,tot_num_of_el,irec,tot_num_of_NNmatel_TUD
            stop
         endif 
         do j=1,irec
            V_NN_TUD(tot_num_of_el+j)=V_NN_TUD(tot_num_of_el+j)
     $           +vx(j)     
         end do

         open(11,file=trim(v3intfile)//'_no1bV3N',status='old',
     $        form ='formatted',action='read')

         write(9,"(' NO2B 3N file used (1-body): ',a)")
     $        trim(v3intfile)//'_no1bV3N'
         
         read(11,*,end=22222)
         read(11,*,end=22222)
         i=0
         do
            read(11,*,end=22222) v2
            i=i+1
         end do
22222    continue
         tot_num_of_lines=i
         tot_num_of_el=10*(i-1)
         print *,' tot_num_of_el+10=',tot_num_of_el+10
         if (tot_num_of_el>2*ist1dim_TUD*(ist1dim_TUD+1)) then
            print *,
     $'*** error in no2bv3n: tot_num_of_el,tot_num_of_no1bv3nme=',
     $           tot_num_of_el,2*ist1dim_TUD*(ist1dim_TUD+1)
            stop
         endif
         if (allocated(temp_TUD)) deallocate(temp_TUD)
         allocate(temp_TUD(tot_num_of_el+10))
         temp_TUD=0.d0
         rewind(11)
         read(11,*)
         read(11,*) V_3N_no0b
         do i=1,tot_num_of_lines-1
            read(11,*) (temp_TUD(10*(i-1)+j),j=1,10)
cc            read(11,*) (V_3N_no1b_TUD(10*(i-1)+j),j=1,10)
         end do
         read(11,'(a)') line
         irec=0
         line=adjustl(line)
         do
            if (len_trim(line)==0) exit 
            irec=irec+1
            i1=index(line,' ')
            line=adjustl(line(i1:))
         end do
         backspace(11)
         read(11,*,err=3353,end=3353) (vx(j),j=1,irec)
         close(11)

         if (tot_num_of_el+irec>2*ist1dim_TUD*(ist1dim_TUD+1)) then
            print *,
     $'*** error in no2bv3n: tot_num_of_el,irec,tot_num_of_no1bv3nme='
     $           ,tot_num_of_el,irec,2*ist1dim_TUD*(ist1dim_TUD+1)
            stop
         endif 
         do j=1,irec
            temp_TUD(tot_num_of_el+j)=vx(j)     
cc            V_3N_no1b_TUD(tot_num_of_el+j)=vx(j)     
         end do
         nljt=0
         do ia=1,ist1dim_TUD
            do mta=0,1
               nljt=nljt+1
               nnljt=0
               ibcyc: do ib=1,ist1dim_TUD
                  do mtb=0,1
                     nnljt=nnljt+1
                     if (nnljt>nljt) exit ibcyc
                     if (mta==mtb) then
                        if (ib>ia) then
                           print *,
     $                          '*** error: ia,ib,mta,mtb,nljt,nnljt=',
     $                          ia,ib,mta,mtb,nljt,nnljt
                           stop
                        endif
                        V_3N_no1b_TUD(ib+ia*(ia-1)/2,mta)=
     $                       temp_TUD(nnljt+nljt*(nljt-1)/2)
                     endif
                  end do
               end do ibcyc
            end do
         end do
      endif
      return
 3353 continue
      print *,' error in reading in read_TUD_nov3n'
      stop
      end subroutine read_TUD_nov3n

      subroutine set_TUD_H(strcm)
      use spb
      use TUD_tbme
      use tbmes, only: rful
      use tbmepar, only: itrel,ihrel
      use nodeinfo
      implicit none
      real(kind(0.0)),intent(IN) :: strcm
      real(kind(0.d0)) :: comfr,comft,comfcm1,comfcm2,mass
      integer :: ii
      real(kind(0.d0)) :: sumH,hcmspe
      integer :: ia,ib,ic,id,na,la,ja,nb,lb,jb
      integer :: nc,lc,jc,nd,ld,jd,i_ab,i_cd,i_abcd,id_max,J12,T12,MT12
      integer :: iTUD,jj
!      real(kind(0.d0)) :: fac_ab,fac_cd
 
      mass=2.d0*938.27231d0*939.56563d0/(938.27231d0+939.56563d0)
      comfr = 4.d0*(197.327053d0)**2/mass
     $     /((real(nucleons,kind(0.d0)))**2*hbomeg)
      comft = 2.d0*hbomeg/real(nucleons,kind(0.d0)) 
      comfcm1 = strcm*hbomeg/real(nucleons-1,kind(0.d0))
      comfcm2 = strcm*2.d0*hbomeg/real(nucleons,kind(0.d0))

      if (allocated(H_NN_TUD)) deallocate(H_NN_TUD)
      allocate(H_NN_TUD(tot_num_of_NNmatel_TUD))
      H_NN_TUD=0.0
      if (allocated(rful)) deallocate(rful)
      allocate(rful(tot_num_of_NNmatel_TUD))
      rful=0.0

      iTUD=0
      do ia=1,ist1dim_TUD
         na=ist1_TUD(1,ia)
         la=ist1_TUD(2,ia)
         ja=ist1_TUD(3,ia)
         do ib=1,ia 
            nb=ist1_TUD(1,ib)
            lb=ist1_TUD(2,ib)
            jb=ist1_TUD(3,ib)
            if (2*na+la+2*nb+lb>N12_max_tud) cycle
!            if (ia==ib) then
!               fac_ab=sqrt(2.d0)
!            else
!               fac_ab=1.d0
!            endif
            i_ab=index_ab(ia,ib)
            do ic=1,ia
               nc=ist1_TUD(1,ic)
               lc=ist1_TUD(2,ic)
               jc=ist1_TUD(3,ic)
               if (ic==ia) then
                  id_max=ib
               else
                  id_max=ic
               endif
               do id=1,id_max 
                  nd=ist1_TUD(1,id)
                  ld=ist1_TUD(2,id)
                  jd=ist1_TUD(3,id)
                  if (2*nc+lc+2*nd+ld>N12_max_tud) cycle
                  if (mod(la+lb+lc+ld,2)==1) cycle
!                  if (ic==id) then
!                     fac_cd=sqrt(2.d0)
!                  else
!                     fac_cd=1.d0
!                  endif
                  i_cd=index_ab(ic,id)
                  i_abcd=index_abcd(i_ab,i_cd)
                  if (ia==ic.and.ib==id.or.ia==id.and.ib==ic) then
                     hcmspe=comfcm1*(2*na+la+2*nb+lb+3.d0
     $                    -3.d0/real(nucleons,kind(0.d0)))
                  else
                     hcmspe=0.d0
                  endif
                  ii=0
                  do J12=max(abs(ja-jb),abs(jc-jd))/2,min(ja+jb,jc+jd)/2
                     do T12=0,1
                        do MT12=-T12,T12
                           iTUD=iTUD+1
                           jj=i_abcd+ii
                           if (ihrel==1) then
                              sumH=(comft-comfcm2)*HOrel_TUD(jj)
     $                             +V_NN_TUD(jj)+hcmspe
                           else
                              sumH=comft*Trel_TUD(jj)+V_NN_TUD(jj)
     $                             -comfcm2*HOrel_TUD(jj)+hcmspe
                           endif
                           H_NN_TUD(jj)=real(sumH,kind(0.0))
                           rful(jj)=comfr*(HOrel_TUD(jj)-Trel_TUD(jj))
                           ii=ii+1
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate(Trel_TUD)
      deallocate(HOrel_TUD)
      deallocate(V_NN_TUD)
      if (iTUD/=tot_num_of_NNmatel_TUD) then
         if (iproc==0) print *,'***error: iTUD,tot_num_of_NNmatel_TUD=',
     $        iTUD,tot_num_of_NNmatel_TUD
      endif      
      end subroutine set_TUD_H

      subroutine pn_tbme_setup
      use tbmepar, only: nobt
      use spodata
      use hmes, only: memoryalloc
      use pn_tbme
      use nodeinfo
      implicit none
      integer :: iunitout=9
      integer :: na,la,j2a,nb,lb,j2b,jt,tz,ipar,ip,nnla,nnlb
      integer :: ii,ia,ib,ipoicount,idimcount,ibmin,isum2,
     $     sp2_min,sp2_max
      allocate(isp2npoi_pn(0:N12_max_pn+1,0:1,-1:1),
     $     isp2ndim_pn(0:N12_max_pn+1,0:1,-1:1),
     $     i2belnpoi_pn(0:N12_max_pn+1,0:1,-1:1))

      allocate(pntbdst(0:N12_max_pn+1,0:1))

      tzloop: do tz=1,-1,-1
      ii=0
      isum2=0
      ipoicount=0
      do jt=0,N12_max_pn+1
         do ip=0,1
            ipar=(-1)**ip
            isp2npoi_pn(jt,ip,tz)=ipoicount
            i2belnpoi_pn(jt,ip,tz)=isum2

            pntbdst(jt,ip)%Tz(tz)%sp1max=nobt
            allocate(pntbdst(jt,ip)%Tz(tz)%sp1(nobt))
            memoryalloc=memoryalloc+nobt*4
            
            idimcount=0
            do ia=1,nobt
               nnla=nnl(ia)
               la=lx(ia)
               j2a=j2x(ia)
               if (nnla>N1_max_pn) exit
               sp2_min=nobt
               sp2_max=1
               if (tz==0) then
                  ibmin=1
               else
                  ibmin=ia
               endif
               do ib=ibmin,nobt
                  nnlb=nnl(ib)
                  lb=lx(ib)
                  j2b=j2x(ib)
                  if (nnlb>N1_max_pn) exit
                  if (nnla+nnlb>N12_max_pn) cycle
                  if (iabs(j2a-j2b)>2*jt.or.(j2a+j2b)<2*jt) cycle
                  if ((-1)**(la+lb)/=ipar) cycle
                  if ((tz==1.or.tz==-1)
     $                 .and.ia==ib.and.(-1)**(jt+1)/=-1) cycle
                  ii=ii+1
                  idimcount=idimcount+1
                  sp2_min=min(sp2_min,ib)
                  sp2_max=max(sp2_max,ib)
               end do
               pntbdst(jt,ip)%Tz(tz)%sp1(ia)%sp2max=sp2_max
               pntbdst(jt,ip)%Tz(tz)%sp1(ia)%sp2min=sp2_min
               allocate(pntbdst(jt,ip)%Tz(tz)%sp1(ia)
     $              %sp2(sp2_min:sp2_max))
               memoryalloc=memoryalloc+(sp2_max-sp2_min+1)*4
               pntbdst(jt,ip)%Tz(tz)%sp1(ia)%sp2(:)=-1
            end do
            ipoicount=ipoicount+idimcount
            isp2ndim_pn(jt,ip,tz)=idimcount
            isum2=isum2+idimcount*(idimcount+1)/2
         end do
      end do
      num_of_2bst_pn(tz)=ii
      num_of_2bel_pn(tz)=isum2
      allocate(isp2n_pn(tz)%ist(4,num_of_2bst_pn(tz)))
      memoryalloc=memoryalloc+4*num_of_2bst_pn(tz)*4
      if (iproc==0) then
         write(iunitout,1000) num_of_2bst_pn(tz),tz
 1000 format(' Number of two-nucleon sp states:',i5,'   for tz=',i2)
         write(iunitout,3456) num_of_2bel_pn(tz),tz
 3456 format(' Number of Hamiltonian 2-b matrix elements:',i8,
     $        '   for tz=',i2)
      endif
      ii=0
      do jt=0,N12_max_pn+1
         do ip=0,1
            ipar=(-1)**ip
            do ia=1,nobt
               nnla=nnl(ia)
               la=lx(ia)
               j2a=j2x(ia)
               if (nnla>N1_max_pn) exit
               if (tz==0) then
                  ibmin=1
               else
                  ibmin=ia
               endif
               do ib=ibmin,nobt
                  nnlb=nnl(ib)
                  lb=lx(ib)
                  j2b=j2x(ib)
                  if (nnlb>N1_max_pn) exit
                  if (nnla+nnlb>N12_max_pn) cycle
                  if (iabs(j2a-j2b)>2*jt.or.(j2a+j2b)<2*jt) cycle
                  if ((-1)**(la+lb)/=ipar) cycle
                  if ((tz==1.or.tz==-1)
     $                 .and.ia==ib.and.(-1)**(jt+1)/=-1) cycle
                  ii=ii+1
                  isp2n_pn(tz)%ist(1,ii)=ia
                  isp2n_pn(tz)%ist(2,ii)=ib
                  isp2n_pn(tz)%ist(3,ii)=jt
                  isp2n_pn(tz)%ist(4,ii)=ip
                  pntbdst(jt,ip)%Tz(tz)%sp1(ia)%sp2(ib)=ii
     $                 -isp2npoi_pn(jt,ip,tz)
               end do
            end do
         end do
      end do
      end do tzloop
      end subroutine pn_tbme_setup

      subroutine read_pntbme(filename)
      use pn_tbme
      use hmes, only: memoryalloc
      use nodeinfo
      implicit none
      include 'mpif.h'
      character(len=120),intent(IN) :: filename
      integer :: tz,i,N1_max_in,N12_max_in,num_of_el_in,
     $     ibuf,iix,ih3oibuf,ierr
      do tz=1,-1,-1
         if (allocated(V_NN_pn(tz)%el))
     $        deallocate(V_NN_pn(tz)%el)
         allocate(V_NN_pn(tz)%el(num_of_2bel_pn(tz)))
         V_NN_pn(tz)%el=0.d0
         if (allocated(Trel_pn(tz)%el))
     $        deallocate(Trel_pn(tz)%el)
         allocate(Trel_pn(tz)%el(num_of_2bel_pn(tz)))
         Trel_pn(tz)%el=0.d0
         if (allocated(HOrel_pn(tz)%el))
     $        deallocate(HOrel_pn(tz)%el)
         allocate(HOrel_pn(tz)%el(num_of_2bel_pn(tz)))
         HOrel_pn(tz)%el=0.d0
         memoryalloc=memoryalloc+num_of_2bel_pn(tz)*8*2 ! eventually only two arrays
      end do
      if (iproc==0) then
         open(22,file=trim(filename)//'_pnVNNbin',status='old',
     $        form ='unformatted',action='read')
         write(9,"(/,' NN interaction file used: ',a)")
     $        trim(filename)//'_pnVNNbin'
         read(22) N1_max_in,N12_max_in
         if (N1_max_in/=N1_max_pn.or.N12_max_in/=N12_max_pn) then
            print *,
     $'***error in read_pntbme file pnVNNbin:
     $N1_max_in,N1_max_pn,N12_max_in,N12_max_pn=',
     $           N1_max_in,N1_max_pn,N12_max_in,N12_max_pn
            stop
         endif
         tzloop: do tz=1,-1,-1
         read(22) num_of_el_in
         if (num_of_el_in/=num_of_2bel_pn(tz)) then
            print *,'***error in read_pntbme file pnVNNbin:
     $tz,num_of_el_in,num_of_2bel_pn=',
     $          tz,num_of_el_in,num_of_2bel_pn(tz)
            stop
         endif
         read(22) (V_NN_pn(tz)%el(i),i=1,num_of_2bel_pn(tz))
         end do tzloop
         close(22)
      
         open(22,file=trim(filename)//'_pnTrelbin',status='old',
     $        form ='unformatted',action='read')
         write(9,"(' Trel file used: ',a)")
     $        trim(filename)//'_pnTrelbin'
         read(22) N1_max_in,N12_max_in
         if (N1_max_in/=N1_max_pn.or.N12_max_in/=N12_max_pn) then
            print *,
     $'***error in read_pntbme file pnTrelbin:
     $N1_max_in,N1_max_pn,N12_max_in,N12_max_pn=',
     $           N1_max_in,N1_max_pn,N12_max_in,N12_max_pn
            stop
         endif
         tzloopTrel: do tz=1,-1,-1
         read(22) num_of_el_in
         if (num_of_el_in/=num_of_2bel_pn(tz)) then
            print *,'***error in read_pntbme file pnTrelbin:
     $tz,num_of_el_in,num_of_2bel_pn=',
     $          tz,num_of_el_in,num_of_2bel_pn(tz)
            stop
         endif
         read(22) (Trel_pn(tz)%el(i),i=1,num_of_2bel_pn(tz))
         end do tzloopTrel
         close(22)
           
         open(22,file=trim(filename)//'_pnHOrelbin',status='old',
     $        form ='unformatted',action='read')
         write(9,"(' HOrel file used: ',a)")
     $        trim(filename)//'_pnHOrelbin'
         read(22) N1_max_in,N12_max_in
         if (N1_max_in/=N1_max_pn.or.N12_max_in/=N12_max_pn) then
            print *,
     $'***error in read_pntbme file pnHOrelbin:
     $N1_max_in,N1_max_pn,N12_max_in,N12_max_pn=',
     $           N1_max_in,N1_max_pn,N12_max_in,N12_max_pn
            stop
         endif
         tzloopHOrel: do tz=1,-1,-1
         read(22) num_of_el_in
         if (num_of_el_in/=num_of_2bel_pn(tz)) then
            print *,'***error in read_pntbme file pnHOrelbin:
     $tz,num_of_el_in,num_of_2bel_pn=',
     $          tz,num_of_el_in,num_of_2bel_pn(tz)
            stop
         endif
         read(22) (HOrel_pn(tz)%el(i),i=1,num_of_2bel_pn(tz))
         end do tzloopHOrel
         close(22)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=3000000
      do tz=1,-1,-1
         if (num_of_2bel_pn(tz)<=ibuf) then
            i=num_of_2bel_pn(tz)
            call MPI_Bcast(V_NN_pn(tz)%el(1),i,MPI_REAL8,0,icomm,ierr)
            call MPI_Bcast(Trel_pn(tz)%el(1),i,MPI_REAL8,0,icomm,ierr) 
            call MPI_Bcast(HOrel_pn(tz)%el(1),i,MPI_REAL8,0,icomm,ierr) 
         else
            ih3oibuf=num_of_2bel_pn(tz)/ibuf
            iix=1
            do i=1,ih3oibuf
               call MPI_Bcast(V_NN_pn(tz)%el(iix),ibuf,MPI_REAL8,0,
     $              icomm,ierr) 
               call MPI_Bcast(Trel_pn(tz)%el(iix),ibuf,MPI_REAL8,0,
     $              icomm,ierr) 
               call MPI_Bcast(HOrel_pn(tz)%el(iix),ibuf,MPI_REAL8,0,
     $              icomm,ierr) 
               if (iproc==0) print *,' Broadcasted iix=',iix
               iix=iix+ibuf
            end do
            i=mod(num_of_2bel_pn(tz),int(ibuf,kind(4)))
            if (i/=0) then
               call MPI_Bcast(V_NN_pn(tz)%el(iix),i,MPI_REAL8,0,icomm,
     $              ierr)
               call MPI_Bcast(Trel_pn(tz)%el(iix),i,MPI_REAL8,0,icomm,
     $              ierr)
               call MPI_Bcast(HOrel_pn(tz)%el(iix),i,MPI_REAL8,0,icomm,
     $              ierr)
               if (iproc==0) print *,' Broadcasted iix=',iix
            endif
         endif
         if (iproc==0) print *,' V_NN_pn, Trel_pn HOrel_pn broadcasted'
      end do
      end subroutine read_pntbme

      subroutine set_pn_H(strcm)
      use spb
      use pn_tbme
      use spodata
      use tbmepar, only: itrel,ihrel,nobt
      use nodeinfo
      implicit none
      real(kind(0.0)),intent(IN) :: strcm
      real(kind(0.d0)) :: comfr,comft,comfcm1,comfcm2,mass
      real(kind(0.d0)) :: sumH,hcmspe
      integer :: ia,ib,ic,id,na,la,ja,nb,lb,jb,nnlc,nnld
      integer :: nc,lc,jc,nd,ld,jd,J12,tz,ip,J12f,ipf
      integer :: ii,ff,sff,sii,index
      real(kind(0.d0)) :: fac_ab,fac_cd
 
      mass=2.d0*938.27231d0*939.56563d0/(938.27231d0+939.56563d0)
      comfr = 4.d0*(197.327053d0)**2/mass
     $     /((real(nucleons,kind(0.d0)))**2*hbomeg)
      comft = 2.d0*hbomeg/real(nucleons,kind(0.d0)) 
      comfcm1 = strcm*hbomeg/real(nucleons-1,kind(0.d0))
      comfcm2 = strcm*2.d0*hbomeg/real(nucleons,kind(0.d0))

      tzloop: do tz=1,-1,-1
         if (allocated(H_NN_pn(tz)%el))
     $        deallocate(H_NN_pn(tz)%el)
         allocate(H_NN_pn(tz)%el(num_of_2bel_pn(tz)))
         H_NN_pn(tz)%el=0.d0
         if (allocated(rful_pn(tz)%el))
     $        deallocate(rful_pn(tz)%el)
         allocate(rful_pn(tz)%el(num_of_2bel_pn(tz)))
         rful_pn(tz)%el=0.d0

         do ff=1,num_of_2bst_pn(tz)
            ic=isp2n_pn(tz)%ist(1,ff)
            id=isp2n_pn(tz)%ist(2,ff)
            J12f=isp2n_pn(tz)%ist(3,ff)
            ipf=isp2n_pn(tz)%ist(4,ff)
            sff=ff-isp2npoi_pn(J12f,ipf,tz)
            if (ic==id.and.tz/=0) then
               fac_cd=sqrt(2.d0)
            else
               fac_cd=1.d0
            endif

            nnlc=nnl(ic)
            nnld=nnl(id)
      
            do sii=sff,isp2ndim_pn(J12f,ipf,tz)
               ii=sii+isp2npoi_pn(J12f,ipf,tz)
               ia=isp2n_pn(tz)%ist(1,ii)
               ib=isp2n_pn(tz)%ist(2,ii)
               J12=isp2n_pn(tz)%ist(3,ii)
               ip=isp2n_pn(tz)%ist(4,ii)
               if (J12/=J12f.or.ip/=ipf) then
                  if (iproc==0)
     $                 print *,'*** error in set_pn_H: 
     $J12,J12f,ip,ipf=',J12,J12f,ip,ipf
                  stop
               endif
               if (ia==ib.and.tz/=0) then
                  fac_ab=sqrt(2.d0)
               else
                  fac_ab=1.d0
               endif

               if (abs(tz)==1.and.
     $              (ia==ic.and.ib==id.or.ia==id.and.ib==ic)) then
                  hcmspe=comfcm1*(nnlc+nnld+3.d0
     $                 -3.d0/real(nucleons,kind(0.d0)))
               elseif (tz==0.and.ia==ic.and.ib==id) then
                  hcmspe=comfcm1*(nnlc+nnld+3.d0
     $                 -3.d0/real(nucleons,kind(0.d0)))
               else
                  hcmspe=0.d0
               endif

               index=i2belnpoi_pn(J12,ip,tz)+sff+sii*(sii-1)/2

               H_NN_pn(tz)%el(index)=(comft*Trel_pn(tz)%el(index)
     $              +V_NN_pn(tz)%el(index)
     $              -comfcm2*HOrel_pn(tz)%el(index)+hcmspe)
     $              *fac_ab*fac_cd

               rful_pn(tz)%el(index)=comfr*(HOrel_pn(tz)%el(index)
     $              -Trel_pn(tz)%el(index))*fac_ab*fac_cd
               
            end do
         end do
         deallocate(V_NN_pn(tz)%el)
         deallocate(Trel_pn(tz)%el)
         deallocate(HOrel_pn(tz)%el)
      end do tzloop
      end subroutine set_pn_H

c     This subroutine sets up all valid m-scheme many-body states
c
      subroutine Setup
      use parameters
      use config
      use ibas
      use nuconf
      use nodeinfo
      use jsgna
      use spodata
      use iters
      use spb
      use spsdata
      use bits
      use tempor
      include 'mpif.h'

c****pn****
c**      integer(4) iik
      integer(8) iik
c*      integer(4),allocatable nsds(:)
c*************************

      real(kind(0.d0)) :: cputime,walltime,timeused,amfss1

c     nprtn(iconf,N): number of protons  in the major shell N in the iconf-th
c                     configuration
c     nneut(iconf,N): number of neutrons in the major shell N in the iconf-th
c                     configuration
      integer(4) nmfss1,nmfss2
c***      integer popcntt    ! comment out for integer(8)
c
c***pn***
      integer :: class,ip,in,prot,neut,mj_p,neng_p,ipar_p,
     $     mj_n,neng_n,ipar_n,sps,nengx,exminho
      integer(2),allocatable :: iltmp(:,:)

      if (iproc==0) then
         iik=bit_size(iik)
         if (iik.ne.nbit) then
            write(6,*)' bitsize=',iik,'   nbit=',nbit
            write(8,*)' bitsize=',iik,'   nbit=',nbit
            call MPI_Abort(icomm,1000,ierr)
            stop
         endif
         iik=0
c***         iik=count((/(btest(ibset(iik,nbit1),ibi),ibi=0,nbit1)/))
c***         iik=popcntt(ibset(iik,nbit1))    ! integer(8) comment out
         iik=popcnt(ibset(iik,nbit1))
         if (iik.ne.1) then
            write(6,*)' # 1bits=',iik,'   for bit=',nbit1
            write(8,*)' # 1bits=',iik,'   for bit=',nbit1
            call MPI_Abort(icomm,1001,ierr)
            stop
         endif
      endif

      if (allocated(locan)) deallocate(locan)
      if (allocated(locn)) deallocate(locn)
      if (allocated(locz)) deallocate(locz)
      allocate(locan(nucleons),locn(mxsps2),locz(mxsps2))
c***************************************************              
c--------------------------------------------------------------------
c     ----- Obtain all valid configurations:
      call clock(cputime,walltime,0,0)
      timeused = cputime
      if(iproc.eq.0)write(8,*)'Valid configurations are:'

      nhwp=nhw-Nmin_HO(nneutrns)
      nhwn=nhw-Nmin_HO(nprotons)
      if (iproc==0) print *,' nhw,nhwp,nhwn=',nhw,nhwp,nhwn

      do class=1,2

      nconf = 0
      call loopconf(nshll,class,0)

*      print *,'#',iproc,' nconf=',nconf

      select case(class)
      case(1)
         if (allocated(nprtn)) then
            deallocate(nprtn)
         endif
         allocate(nprtn(nconf,nshll))
      case(2)
         if (allocated(nneut)) then
            deallocate(nneut)
         end if
         allocate(nneut(nconf,nshll))
      end select
c      allocate(ndiff(nconf,nconf),negy(nconf))
      if (allocated(negy)) deallocate(negy)
      allocate(negy(nconf))
c      ndiff=0
      negy=0
      nconf = 0
      call loopconf(nshll,class,1)
c
      if(iproc.eq.0)then
         write(8,*)
         write(8,*)'Maximum number of s.p. states allowed: ',mxsps
         write(8,*)'Number of s.p. states actually used:   ',nasps
         write(9,*)
         write(9,*)'Maximum number of s.p. states allowed: ',mxsps
         write(9,*)'Number of s.p. states actually used:   ',nasps
      endif
c     nmfss1 = nfac(nasps,nprotons)*nfac(nasps,nneutrns)
c      select case(class)
c      case(1)
c         amfss1 = real(nfac(nasps,nprotons),kind(0.d0))
c      case(2)
c         amfss1 = real(nfac(nasps,nneutrns),kind(0.d0))
c      end select
      if (iproc==0) then
c         write(6,*)'Number of unrestricted many-body states:   ',
c     &        amfss1
c         write(8,*)'Number of unrestricted many-body states:   ',
c     &        amfss1
c         write(9,*)'Number of unrestricted many-body states:   ',
c     &        amfss1
         write(6,*)'Number of valid configurations:        ',nconf
         write(8,*)'Number of valid configurations:        ',nconf
         write(9,*)'Number of valid configurations:        ',nconf
      endif

c      if(nconf.eq.0)stop 'nconf = 0.'

c     ----- Now that all the valid configurations have been found, proceed
c     ----- to finding out valid m-scheme states for each configuration.

      ipass=0
      nsdx = 0
      nmfss2x = 0
      do  iconf = 1, nconf, 1

         select case(class)
         case(1)
            npib(1) = 0
            do i = 2, nshll, 1
               npib(i) = npib(i-1) + nprtn(iconf,i-1)
            enddo
c     
c     ----- Upgrade the bit pattern & check if the new pattern is valid:
            do mj = 1, nshll, 1
               nsps0c = nsps0(mj)
               nops(mj) = nfac(nsps0c,nprtn(iconf,mj)) - 1
            enddo

         case(2)
            nnub(1) = 0
            do i = 2, nshll, 1
               nnub(i) = nnub(i-1) + nneut(iconf,i-1)
            enddo
c     
c     ----- Upgrade the bit pattern & check if the new pattern is valid:
            do mj = 1, nshll, 1
               nsps0c = nsps0(mj)
               nons(mj) = nfac(nsps0c,nneut(iconf,mj)) - 1
            enddo
         end select
c     
         do ibynd = nasps+1, mxsps, 1
            locn(ibynd) = 0
            locn(ibynd+mxsps) = 0
         enddo

         call loopbit(nshll,class,iconf,ipass)

      end do
c
c*      print *,'#',iproc,' nsdx=',nsdx

cc      if (allocated(mconf)) deallocate(mconf)
cc      allocate(mconf(nsdx))
c      if (allocated(iendconf)) deallocate(iendconf)
c      allocate(iendconf(nconf))
cc      mconf=0
      istartconf=0
      
      select case(class)
      case(1)
         nsdp=nsdx
         if (major==2) then
            if (allocated(ilocp)) deallocate(ilocp)
            allocate(ilocp(nprotons,nsdp))
            ilocp=0
         else   
            if (allocated(ibasis)) deallocate(ibasis)
            allocate(ibasis(nwd,2,nsdx))
            ibasis=0
         endif
      case(2)
         nsdn=nsdx
         if (major==2) then
            if (allocated(ilocn)) deallocate(ilocn)
            allocate(ilocn(nneutrns,nsdn))
            ilocn=0
         else   
            if (allocated(ibasis)) deallocate(ibasis)
            allocate(ibasis(nwd,2,nsdn))
            ibasis=0
         endif
      end select
      ipass=1
      nsdx = 0
      nmfss2x = 0
      do  iconf = 1, nconf, 1

         select case(class)
         case(1)
            npib(1) = 0
            do i = 2, nshll, 1
               npib(i) = npib(i-1) + nprtn(iconf,i-1)
            enddo
c     
c     ----- Upgrade the bit pattern & check if the new pattern is valid:
            do mj = 1, nshll, 1
               nsps0c = nsps0(mj)
               nops(mj) = nfac(nsps0c,nprtn(iconf,mj)) - 1
            enddo

         case(2)
            nnub(1) = 0
            do i = 2, nshll, 1
               nnub(i) = nnub(i-1) + nneut(iconf,i-1)
            enddo
c     
c     ----- Upgrade the bit pattern & check if the new pattern is valid:
            do mj = 1, nshll, 1
               nsps0c = nsps0(mj)
               nons(mj) = nfac(nsps0c,nneut(iconf,mj)) - 1
            enddo
         end select
c     
         do ibynd = nasps+1, mxsps, 1
            locn(ibynd) = 0
            locn(ibynd+mxsps) = 0
         enddo

         call loopbit(nshll,class,iconf,ipass)

c         iendconf(iconf)=nsdx
*         if (iproc==0) 
*     +        print *,' iconf=',iconf,'  iendconf=',iendconf(iconf)

      end do
c

      nmfss2 = nmfss2x
      
      if(iproc.eq.0)then
         select case(class)
         case(1)
            write(6,*)'No.of proton SDs after Energy test: ',nmfss2
            write(8,*)'No.of proton SDs after Energy test: ',nmfss2
            write(9,*)'No.of proton SDs after Energy test: ',nmfss2
         case(2)
            write(6,*)'No.of neutron SDs after Energy test: ',nmfss2
            write(8,*)'No.of neutron SDs after Energy test: ',nmfss2
            write(9,*)'No.of neutron SDs after Energy test: ',nmfss2
         end select
      endif
c     
      end do  ! class

      deallocate(nprtn,nneut)
      deallocate(negy)
c      deallocate(iendconf)

      if (allocated(iltmp)) deallocate(iltmp)
      allocate(iltmp(nprotons,nsdp))
      iltmp=ilocp
      nminho=Nmin_HO(nprotons)
      nsdx=0
      do nengx=nminho,nhwp
         do ip=1,nsdp
            neng_p=0
            do prot=1,nprotons
               sps=iltmp(prot,ip)
               neng_p=neng_p+2*n_sp(sps)+l_sp(sps)
            end do
            if (nengx==neng_p) then
               nsdx=nsdx+1
               ilocp(:,nsdx)=iltmp(:,ip)
            endif
         end do
      end do
c      print *,' nsdx,nsdp=',nsdx,nsdp
      deallocate(iltmp)
      allocate(iltmp(nneutrns,nsdn))
      iltmp=ilocn
      nminho=Nmin_HO(nneutrns)
      nsdx=0
      do nengx=nminho,nhwn
         do ip=1,nsdn
            neng_p=0
            do prot=1,nneutrns
               sps=iltmp(prot,ip)
               neng_p=neng_p+2*n_sp(sps)+l_sp(sps)
            end do
            if (nengx==neng_p) then
               nsdx=nsdx+1
               ilocn(:,nsdx)=iltmp(:,ip)
            endif
         end do
      end do
c      print *,' nsdx,nsdp=',nsdx,nsdp
      deallocate(iltmp)

      nminho=Nmin_HO(nprotons)+Nmin_HO(nneutrns)
      exminho=mod(nhw-nminho,2)
      if (allocated(nsd_hw)) deallocate(nsd_hw)
      allocate(nsd_hw(exminho:nhw-nminho))
      nsdx=0
      do nengx=nminho+exminho,nhw,2
         do ip=1,nsdp
            mj_p=0
            neng_p=0
            ipar_p=0
            do prot=1,nprotons
               sps=ilocp(prot,ip)
               mj_p=mj_p+m2_sp(sps)
               neng_p=neng_p+2*n_sp(sps)+l_sp(sps)
               ipar_p=ipar_p+l_sp(sps)
            end do

            if (neng_p>nengx) exit

c         print *,' ip,neng_p=',ip,neng_p

            do in=1,nsdn
               mj_n=0
               neng_n=0
               ipar_n=0
               do neut=1,nneutrns
                  sps=ilocn(neut,in)
                  mj_n=mj_n+m2_sp(sps)
                  neng_n=neng_n+2*n_sp(sps)+l_sp(sps)
                  ipar_n=ipar_n+l_sp(sps)
               end do

               if (neng_p+neng_n>nengx) exit

               if (neng_p+neng_n==nengx) then
                  if (mj_n+mj_p/=mjtotal) cycle
                  if (mod(ipar_p+ipar_n,2)/=iparity) cycle
                  nsdx=nsdx+1
c                  if (nengx==nminho) nsd_hw0=nsdx
               endif
         
            end do
         end do
         nsd_hw(nengx-nminho)=nsdx
         if (iproc==0) print *,' nengx,Nmax,nsd_hw=',nengx,nengx-nminho,
     $        nsd_hw(nengx-nminho)
      end do

      nsd = nsdx
c     
      if(iproc.eq.0)then
         write(6,*)'Number of valid SDs after the M test: ',nsd
         write(8,*)'Number of valid SDs after the M test: ',nsd
         write(9,*)'Number of valid SDs after the M test: ',nsd
      endif

      if(nsd.eq.0)then
         if(iproc.eq.0)write(8,*)'  Error: nsd = 0.'
         stop 'nsd = 0.'
      endif

      if (allocated(iloc)) deallocate(iloc)
      allocate(iloc(2,nsd))
      nsdx=0
cc      do nengx=nminho,nhw,2
      do nengx=nminho+exminho,nhw,2   ! a bug was here
         do ip=1,nsdp
            mj_p=0
            neng_p=0
            ipar_p=0
            do prot=1,nprotons
               sps=ilocp(prot,ip)
               mj_p=mj_p+m2_sp(sps)
               neng_p=neng_p+2*n_sp(sps)+l_sp(sps)
               ipar_p=ipar_p+l_sp(sps)
            end do

            if (neng_p>nengx) exit

c         print *,' ip,neng_p=',ip,neng_p

            do in=1,nsdn
               mj_n=0
               neng_n=0
               ipar_n=0
               do neut=1,nneutrns
                  sps=ilocn(neut,in)
                  mj_n=mj_n+m2_sp(sps)
                  neng_n=neng_n+2*n_sp(sps)+l_sp(sps)
                  ipar_n=ipar_n+l_sp(sps)
               end do

               if (neng_p+neng_n>nengx) exit

               if (neng_p+neng_n==nengx) then
                  if (mj_n+mj_p/=mjtotal) cycle
                  if (mod(ipar_p+ipar_n,2)/=iparity) cycle
                  nsdx=nsdx+1
                  iloc(1,nsdx)=ip
                  iloc(2,nsdx)=in
               endif
         
            end do
         end do
      end do
      if (nsdx/=nsd) then
         print *,' error: nsd,nsdx=',nsd,nsdx
         stop
      endif

      call clock(cputime,walltime,0,0)
      timeused = cputime - timeused
      if(iproc.eq.0)write(8,*)
     &     'CPU time spent in SETUP (seconds): ', timeused
      return
      end subroutine Setup

      integer function Nmin_HO(Z)
      implicit none
      integer,intent(IN) :: Z
      integer :: Nmin,N,Zrem
      Zrem=Z
      Nmin=0
      N=0
      do
         Nmin=Nmin+N*min((N+1)*(N+2),Zrem)
         Zrem=Zrem-(N+1)*(N+2)
         if (Zrem<=0) exit
         N=N+1
      end do
      Nmin_HO=Nmin
      end function Nmin_HO

      recursive subroutine loopconf(Nin,iclass,ipass)
      use parameters
      use config
      use nuconf
      use nodeinfo
      use jsgna
      use spb

      integer, intent(IN) :: Nin,iclass,ipass
      integer np,nn,iloop

***      print *,' nshll, Nin, iloop:',nshll,Nin,iloop

      iloop=nshll-Nin+1
      select case (iclass)
      case (1)
         do np = mxpi(iloop),mnpi(iloop),-1
            if (iloop==1) then
               npib(1)=np
            else
               if (npib(iloop-1)+np>nprotons) then
                  cycle
               else
                  npib(iloop)=npib(iloop-1)+np
               endif
            endif
            ntempi(iloop)=np
***            print *,' #',iproc,' Nin,np,nn:',Nin,np,nn
***            print *,' #',iproc,' iloop=',iloop
            if (Nin>1) then
               call loopconf(Nin-1,iclass,ipass)
            elseif (Nin==1) then

***            write(8,*)' #',iproc,' ntempi,ntempn:'
***            write(8,*)(ntempi(i),ntemnu(i),i=1,nshll)

               npi = 0
               nengy = 0
               do i = 1, nshll
                  npi = npi + ntempi(i)
                  nengy = nengy + neng(i)*ntempi(i)
               enddo
               if (npi/=nprotons.or.nengy>nhwp) cycle
c**               if (nengy>nhw) cycle
c     
c               nodd = 0
c               do i = 1,nshll
c                  if (jsgn(neng(i))==-1) then
c                     nodd = nodd + ntempi(i) + ntemnu(i)
c                  endif
c               enddo
c               if (jsgn(nodd+iparity)==-1) cycle
c     This configuration survives the Nhw test and the parity test.
               nconf = nconf + 1
               if (ipass==0) cycle
               negy(nconf) = nengy
               do i = 1, nshll
                  nprtn(nconf,i) = ntempi(i)
               enddo
               if(iproc==0)then
                  write(8,3)nconf,(nprtn(nconf,i),i=1,nshll)
 3                format(1x,'Conf #',i5,16(i3))
               endif
               
            else
               print *,'*** Error in loopcnf ***'
               stop
            endif    
         end do   

      case (2)
         do nn = mxnu(iloop),mnnu(iloop),-1
c         do nn = min(mxnu(iloop),mxpn(iloop)-np),
c     +           max(mnnu(iloop),mnpn(iloop)-np),-1
c**         do nn = mxnu(iloop),mnnu(iloop),-1
c**            if ((np+nn)<mnpn(iloop).or.(np+nn)>mxpn(iloop)) cycle 
            if (iloop==1) then
               nnub(1)=nn
            else
               if (nnub(iloop-1)+nn>nneutrns) then
                  cycle
               else
                  nnub(iloop)=nnub(iloop-1)+nn
               endif
            endif
            ntemnu(iloop)=nn
***            print *,' #',iproc,' Nin,np,nn:',Nin,np,nn
***            print *,' #',iproc,' iloop=',iloop
            if (Nin>1) then
               call loopconf(Nin-1,iclass,ipass)
            elseif (Nin==1) then

***            write(8,*)' #',iproc,' ntempi,ntempn:'
***            write(8,*)(ntempi(i),ntemnu(i),i=1,nshll)

               nnu = 0
               nengy = 0
               do i = 1, nshll
                  nnu = nnu + ntemnu(i)
                  nengy = nengy + neng(i)*ntemnu(i)
               enddo
               if (nnu/=nneutrns.or.nengy>nhwn) cycle
c**               if (nengy>nhw) cycle
c     
c               nodd = 0
c               do i = 1,nshll
c                  if (jsgn(neng(i))==-1) then
c                     nodd = nodd + ntempi(i) + ntemnu(i)
c                  endif
c               enddo
c               if (jsgn(nodd+iparity)==-1) cycle
c     This configuration survives the Nhw test and the parity test.
               nconf = nconf + 1
               if (ipass==0) cycle
               negy(nconf) = nengy
               do i = 1, nshll
                  nneut(nconf,i) = ntemnu(i)
               enddo
               if(iproc==0)then
                  write(8,3)nconf,(nneut(nconf,i),i=1,nshll)
c 3                format(1x,'Conf #',i5,16(i3))
               endif
               
            else
               print *,'*** Error in loopcnf ***'
               stop
            endif    

         end do 
      end select
      end 


      recursive subroutine loopbit(Nin,iclass,iconf,ipass)
      use parameters
      use config
      use ibas
      use nuconf
      use iters
      use spb
      use spsdata
      use bits
      use tempor

      integer, intent(IN) :: Nin,iclass,iconf,ipass
      integer is,loopmax

      if (iclass==2) then 
	 loopmax = nons(Nin)	   
      elseif (iclass==1) then
	 loopmax = nops(Nin)
      else
         print *,'*** Error in loopbit ***'
         stop	
      endif    	 	   	 
      call initbitptn(iconf,iclass,Nin)
      do is = 0, loopmax, 1
         if (is>0) call bitgen(iclass,Nin)
         if (Nin>1) then
            call loopbit(Nin-1,iclass,iconf,ipass)
c         elseif (Nin==1.and.ncls>1) then
c            call loopbit(nshll,iclass,iconf,ncls-1,ipass)
c         elseif (Nin==1.and.ncls==1) then
         elseif (Nin==1) then
c
            nmfss2x = nmfss2x + 1
c     
c     ----- Inner most loop: Check if the sum of mz is equal to mjtotal:
c            m2sum = 0
c            do nucleon = 1, nucleons, 1
c               m2sum = m2sum + m2_sp(locan(nucleon))
c            enddo
c            if (m2sum/=mjtotal) cycle
            nsdx = nsdx + 1

            if (ipass==0) cycle

cc            mconf(nsdx) = iconf
            
c            if(major/=2)then
c               do ix = 1, nucleons
c                  locx = locan(ix) - 1
c                  icls = locx/mxsps 
c                  locx = locx - icls*mxsps
c                  icls = icls + 1
c                  iwd = locx/nbit
c                  ibit = locx - iwd*nbit
c                  iwd = iwd + 1
c                  ibasis(iwd,icls,nsdx) = 
c     +                 ibset(ibasis(iwd,icls,nsdx),ibit)
c               enddo
c            else
               select case(iclass)
               case(1)
                  do ix = 1, nprotons
                     ilocp(ix,nsdx) = locan(ix)
                  enddo
               case(2)
                  do ix = 1, nneutrns
                     ilocn(ix,nsdx) = locan(ix)
                  enddo
               end select
c            endif
            
         else
            print *,'*** Error in loopbit ***'
            stop
         endif    
      end do

      end

c
c     This subroutine initializes the bit pattern for class "icls" 
c     within the major shell (N=mj). Called by SETUP.
c
      subroutine initbitptn(iconf,icls,mj)
      use parameters
      use config
      use spb
      use bits
      integer,intent(IN) :: iconf,icls,mj

      if(mj.gt.nshll)return
      goto (1,2),icls
 1    nspsr = nspsb(mj)
      npibr = npib(mj)
      nprtnr = nprtn(iconf,mj)
      do i = 1, nprtnr, 1
         nspsr = nspsr + 1
         locn(nspsr) = 1
         locan(npibr+i) = nspsr
      enddo
      do i = nprtnr+1, nsps0(mj), 1
         nspsr = nspsr + 1
         locn(nspsr) = 0
      enddo
      return
 2    nspsr = mxsps + nspsb(mj)
c      nnubr = nprotons + nnub(mj)
      nnubr = nnub(mj)
      nneutr = nneut(iconf,mj)
      do i = 1, nneutr, 1
         nspsr = nspsr + 1
         locn(nspsr) = 1
         locan(nnubr+i) = nspsr
      enddo
      do i = nneutr+1, nsps0(mj), 1
         nspsr = nspsr + 1
         locn(nspsr) = 0
      enddo
      return
      end


c
c     This subroutine generates a new bit pattern based on the old bit pattern.
c     Called by SETUP.
c     
      subroutine bitgen(icls,mj)
      use parameters
      use config
      use nodeinfo
      use spb
      use bits
      include 'mpif.h'
c     Take locn and generate the next higher bit pattern if possible, 
c     update locan accordingly.
      ibloc = 0
      if(icls.eq.1)then
         ioff = nspsb(mj)
         iofn = npib(mj)
      else
         ioff = nspsb(mj)+mxsps
c         iofn = nnub(mj)+nprotons
         iofn = nnub(mj)
      endif
      do 100 i = ioff+1, ioff+nsps0(mj)-1, 1
         if(locn(i).eq.0)goto 100
         if(locn(i+1).ne.0)then
            ibloc = ibloc + 1
            goto 100
         endif
c     Promote, downshift if needed and return
         locn(i) = 0
         locn(i+1) = 1
         locan(iofn+ibloc+1) = i + 1
c     If lowest occuppied orbital has been promoted on the previous call
c     and cannot be promoted again, then downshifting is also required.
c     Check for this case:
         if(ibloc*locn(ioff+1).gt.0)return
c     Downshifting required of lowest surviving continuous string of 1's:
         do j = ioff+1, ioff+ibloc, 1
            locn(j) = 1
            locan(iofn+j-ioff) = j
         enddo
         do j = ioff+ibloc+1, i-1, 1
            locn(j) = 0
         enddo
         return
 100  continue
c     Arriving here means no promotion of sps possible
      if(iproc.eq.0)write(8,*)'  Error:  mj = ',mj,'  icls =',icls
      call MPI_Abort(icomm,106,ierr)
      stop 'Error: Failed to promote the bit pattern.'
      end


c     The binary polynomial function (n,m) = n!/[(n-m)!m!]:
      function nfac(n,m)
      use nodeinfo
      include 'mpif.h'
      if(m.gt.n) then
         if(iproc.eq.0)write(8,*)'  Error: m > n in nfac: ',m,n
         stop 'm > n.'
      endif
      nfac = 1
      do i = 1, m
         nfac = (nfac*(n-i+1))/i
      enddo
      return
      end


c     MKGK routine
c     taken from Lanczos hash
      subroutine make_hash_iloc(ilocsave)
      use hash_tables
      use occ
      use ibas
      use nodeinfo
      use spb
      use iters
      use parameters
      use ITvar
      implicit none
      include 'mpif.h'
      logical,intent(IN) :: ilocsave
      integer :: i1, cr_1,wd, ibit, ionb
      integer(8), allocatable :: intbasi(:)
      integer, allocatable :: occi(:)
      integer :: prot_st,neut_st
c     TEST
      integer :: k,i
c* PN
      integer :: dimfac,type
      character(len=8) :: line
      character(len=100) :: filenmtmp

      if (associated(ty%nucl)) deallocate(ty%nucl)
      if (associated(ty%protn)) deallocate(ty%protn)
      if (associated(ty%neutn)) deallocate(ty%neutn)
      if (associated(ty%dim_st)) deallocate(ty%dim_st)
      if (associated(ty%dim_st_p)) deallocate(ty%dim_st_p)
      if (associated(ty%dim_st_n)) deallocate(ty%dim_st_n)

      allocate(ty%nucl)
      allocate(ty%protn)
      allocate(ty%neutn)
      allocate(ty%dim_st)
      allocate(ty%dim_st_p)
      allocate(ty%dim_st_n)
      ty%nucl=nucleons
      ty%protn=nprotons
      ty%neutn=nneutrns
      ty%dim_st=nsd
      ty%dim_st_p=nsdp
      ty%dim_st_n=nsdn
      ty%occ=>iloc
      ty%occp=>ilocp
      ty%occn=>ilocn

c     TEST code
c      open(59,file='iloc.tmp',status='unknown')
c      open(61,file='intbas.tmp',status='unknown')

c* PN
      if (nhw_boot==nhw_min-2.or.
     $     (nhw_boot==nhw_min.and.irest==6.and.nhw_restart==nhw_boot
     $     .and.kappa_restart==-1)) then
         if (nsd>50000000) then
            dimfac=4
         elseif (nsd>22000000) then
            dimfac=6
         elseif (nsd>9000000) then
            dimfac=10
         elseif (nsd>3000000) then
            dimfac=20
         else
            dimfac=100
         endif
         type=1
      elseif (nhw_boot>=nhw_min) then
         if (nsd>50000000) then
            dimfac=3
         elseif (nsd>22000000) then
            dimfac=5
         elseif (nsd>9000000) then
            dimfac=8
         elseif (nsd>3000000) then
            dimfac=12
         else
            dimfac=50
         endif
         type=1
      else
         dimfac=1
         type=0
      endif

      call construct_hash_table(dimfac,type)
      if (iproc==0) print *,' construct_hash_table called'
      allocate(intbasi(2*nwd))
      allocate(occi(nucleons))
      do i1=1,nsd
         prot_st=iloc(1,i1)
         neut_st=iloc(2,i1)
         occi(1:nprotons)=ilocp(:,prot_st)
         occi(nprotons+1:nucleons)=ilocn(:,neut_st)
         intbasi=0
         do cr_1=1,nucleons
            wd=(occi(cr_1)-1)/nbit+1
            ibit=mod(occi(cr_1)-1,nbit)
c     print *,' i,wd,ibit=',i,wd,ibit
            intbasi(wd)=ibset(intbasi(wd),ibit)
         end do

c     TEST code
c         write(59,*) occi(:)
c         write(61,199) (intbasi(k),k=1,4)
c 199     format(4(I20,2x))
c     end TEST

         call get_state_index(nucleons,2*nwd,intbasi,ionb)
         if (i1/=ionb) then
            print *,'*** error:'
            if (iproc==0) then
               print *,' i=',i1,'    index=',ionb
               print *, iloc(:,i1)
               print *,' occi=',occi
            endif
            stop
         endif
      end do

c     TEST 
c      close(59)
c      close(61)

      if (iproc==0.and.ilocsave) then
         write(line,'(i8)') nhw_boot
         line=adjustl(line)
         filenmtmp = 'iloc_'//trim(line)//'.tmp'
         open(47,file=trim(filenmtmp),form='unformatted',
     $        status='unknown',action='write')
         write(47) nsd,nsdp,nsdn
!         write(47) (iloc(1,i1),i1=1,nsd)
!         write(47) (iloc(2,i1),i1=1,nsd)
         write(47) ((iloc(i,i1),i=1,2),i1=1,nsd)
         write(47) ((ilocp(i,i1),i=1,nprotons),i1=1,nsdp)
         write(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,nsdn)
         if (nhw_boot>=nhw_min.and.allocated(delta_IT)) then
            write(47) (delta_IT(i1),i1=1,nsd)
            print *,' delta_IT saved'
         endif
         close(47)
      endif

      end subroutine make_hash_iloc

      subroutine iloc_read(Nmax)
      use ibas
      use spb
      use iters, only: nsd
      use ITvar
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer(2),intent(IN) :: Nmax
      character(len=8) :: line
      integer :: i1,i,temp(3)
      integer(2),allocatable :: ilocN_tmp_buf1(:)
      integer(4),allocatable :: iloc_tmp_buf1(:),iloc_tmp(:,:) !,
!     $     iloc_tmp1(:,:)
      integer :: ibuf,iallredbuf,iiy,nucl,iallred,ierr,nsdkappa
      character(len=100) :: filenmtmp
      real(8) :: cputime,walltime,timeused

      if (iproc==0) then
         call clock(cputime,walltime,0,0)
         timeused = cputime
         write(line,'(i8)') Nmax
         line=adjustl(line)
         filenmtmp = 'iloc_'//trim(line)//'.tmp'
         open(47,file=trim(filenmtmp),form='unformatted',
     $        status='old',action='read')
         read(47) nsd,nsdp,nsdn
         temp(1)=nsd
         temp(2)=nsdp
         temp(3)=nsdn
         print *,' read nsd,nsdp,nsdn=',nsd,nsdp,nsdn
         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
         print *,' CPU time used=',timeused
         timeused = cputime
         if (allocated(iloc)) deallocate(iloc)
         allocate(iloc(2,nsd))
!         read(47) (iloc(1,i1),i1=1,nsd)
!         print *,' iloc(1) read'
         read(47) ((iloc(i,i1),i=1,2),i1=1,nsd)
         print *,' iloc read'
         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
         print *,' CPU time used=',timeused
         timeused = cputime
!         read(47) (iloc(2,i1),i1=1,nsd)
!         print *,' iloc(2) read'
!         call clock(cputime,walltime,0,0)
!         timeused = cputime - timeused
!         print *,' CPU time used=',timeused
!         timeused = cputime
         if (allocated(ilocp)) deallocate(ilocp)
         allocate(ilocp(nprotons,nsdp))
         read(47) ((ilocp(i,i1),i=1,nprotons),i1=1,nsdp)
         print *,' ilocp read'
         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
         print *,' CPU time used=',timeused
         timeused = cputime
         if (allocated(ilocn)) deallocate(ilocn)
         allocate(ilocn(nneutrns,nsdn))
         read(47) ((ilocn(i,i1),i=1,nneutrns),i1=1,nsdn)
         print *,' ilocn read'
         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
         print *,' CPU time used=',timeused
         timeused = cputime
         if (Nmax>=nhw_min.and.kappa>kappa_store(1)) then
            if (allocated(delta_IT)) deallocate(delta_IT)
            allocate(delta_IT(nsd))
            read(47) (delta_IT(i1),i1=1,nsd)
            print *,' delta read'
            call clock(cputime,walltime,0,0)
            timeused = cputime - timeused
            print *,' CPU time used=',timeused
            timeused = cputime
            allocate(iloc_tmp(2,nsd))
            nsdkappa=0
            do i1=1,nsd
               if (abs(delta_IT(i1))>=kappa) then
                  nsdkappa=nsdkappa+1
                  iloc_tmp(:,nsdkappa)=iloc(:,i1)
               endif
            end do
            call clock(cputime,walltime,0,0)
            timeused = cputime - timeused
            print *,' CPU time used=',timeused
            timeused = cputime
            print *,' kappa,nsd,nsdkappa=',kappa,nsd,nsdkappa
            print *,' nsdp,nsdn=',nsdp,nsdn
            deallocate(delta_IT)
            nsd=nsdkappa
            temp(1)=nsd
            deallocate(iloc)
            allocate(iloc(2,nsd))
            do i1=1,nsd
               iloc(:,i1)=iloc_tmp(:,i1)
            end do            
            deallocate(iloc_tmp)
            call clock(cputime,walltime,0,0)
            timeused = cputime - timeused
            print *,' CPU time used=',timeused
            timeused = cputime
            print *,' iloc reset'
         endif
         if (Nmax>=nhw_min.and.kappa==kappa_store(1)) then
            print *,' kappa,nsd=',kappa,nsd
            print *,' nsdp,nsdn=',nsdp,nsdn
         endif
         close(47)
      endif
      call MPI_Barrier(icomm,ierr)
      call MPI_Bcast(temp,3,MPI_INTEGER,0,icomm,ierr)
      if (iproc/=0) then
         nsd=temp(1)
         nsdp=temp(2)
         nsdn=temp(3)
         if (allocated(iloc)) deallocate(iloc)
         if (allocated(ilocp)) deallocate(ilocp)
         if (allocated(ilocn)) deallocate(ilocn)
         allocate(iloc(2,nsd))
         allocate(ilocp(nprotons,nsdp))
         allocate(ilocn(nneutrns,nsdn))
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsd*2,3000000)
      if (nsd*2<=ibuf) then
         call MPI_Bcast(iloc(1,1),nsd*2,MPI_INTEGER,0,
     $        icomm,ierr)
      else
         if (allocated(iloc_tmp_buf1)) deallocate(iloc_tmp_buf1)
         allocate(iloc_tmp_buf1(ibuf))
         do nucl=1,2
            iallredbuf=nsd/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  iloc_tmp_buf1(:)=iloc(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(iloc_tmp_buf1,ibuf,MPI_INTEGER,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  iloc(nucl,iiy:iiy+ibuf-1)=iloc_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsd,ibuf)
            if (iallred/=0) then
               if (iallred/=nsd-iiy+1) then
                  print *,'***error:iallred,nsd,iiy:',iallred,nsd,iiy
               endif
               if (iproc==0) then
                  iloc_tmp_buf1(1:iallred)=iloc(nucl,iiy:nsd)
               endif
               call MPI_Bcast(iloc_tmp_buf1,iallred,MPI_INTEGER,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  iloc(nucl,iiy:nsd)=iloc_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(iloc_tmp_buf1)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdp*nprotons,3000000)
      if (nsdp*nprotons<=ibuf) then
         call MPI_Bcast(ilocp(1,1),nsdp*nprotons,MPI_INTEGER2,0,
     $        icomm,ierr)
      else
         if (allocated(ilocN_tmp_buf1)) deallocate(ilocN_tmp_buf1)
         allocate(ilocN_tmp_buf1(ibuf))
         do nucl=1,nprotons
            iallredbuf=nsdp/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  ilocN_tmp_buf1(:)=ilocp(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,ibuf,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocp(nucl,iiy:iiy+ibuf-1)=ilocN_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsdp,ibuf)
            if (iallred/=0) then
               if (iallred/=nsdp-iiy+1) then
                  print *,'***error:iallred,nsdp,iiy:',iallred,nsdp,iiy
               endif
               if (iproc==0) then
                  ilocN_tmp_buf1(1:iallred)=ilocp(nucl,iiy:nsdp)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,iallred,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocp(nucl,iiy:nsdp)=ilocN_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(ilocN_tmp_buf1)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdn*nneutrns,3000000)
      if (nsdn*nneutrns<=ibuf) then
         call MPI_Bcast(ilocn(1,1),nsdn*nneutrns,MPI_INTEGER2,0,
     $        icomm,ierr)
      else
         if (allocated(ilocN_tmp_buf1)) deallocate(ilocN_tmp_buf1)
         allocate(ilocN_tmp_buf1(ibuf))
         do nucl=1,nneutrns
            iallredbuf=nsdn/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  ilocN_tmp_buf1(:)=ilocn(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,ibuf,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocn(nucl,iiy:iiy+ibuf-1)=ilocN_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsdn,ibuf)
            if (iallred/=0) then
               if (iallred/=nsdn-iiy+1) then
                  print *,'***error:iallred,nsdp,iiy:',iallred,nsdn,iiy
               endif
               if (iproc==0) then
                  ilocN_tmp_buf1(1:iallred)=ilocn(nucl,iiy:nsdn)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,iallred,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocn(nucl,iiy:nsdn)=ilocN_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(ilocN_tmp_buf1)
      endif
      if (iproc==0) write(9,"(/,' iloc read: nsd,nsdp,nsdn=',i12,2i9)")
     $     nsd,nsdp,nsdn
      end subroutine iloc_read

      subroutine importance_truncate
      use parameters
      use ibas
      use nuconf
      use lanczvec
      use nodeinfo
      use multig
      use ctrl
      use iters
      use spb
      use bits
      use albe
c      use hamc
c      use i1i3
c      use hmes
      use nnmaxmod
      use tbmes
      use tbmepar
      use v3b
      use spsdata
      use occ
      use hash_tables
      use ITvar
      use TUD_tbme              !, only: tbmeTUD
      use pn_tbme, only: tbmepn
      implicit none
      include 'mpif.h'

      integer :: i1,i3
      integer :: ii,jj,kk,ll

c---------------------------------------------------------------------------
c**      integer(4) itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)
      integer(8) itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)

      real(kind(0.d0)) xamp,xxi,ovlp,ovlpmx,alph,betsq,betapre,x00

      real(8) cputime,walltime,timeused
      real(8),allocatable:: bmp0(:)
c      real(8),allocatable:: bmp1(:)
      real(kind=kind(0.d0)) denom
      integer ncut,mionb,ionb,ierr
c**      integer(4),allocatable :: ibascomp(:)
      integer(8),allocatable :: ibascomp(:)
c**      integer(4) ibascomi
      integer(8) ibascomi
c**      integer popcntt   ! integer(8) comment out
      integer stat(MPI_STATUS_SIZE)
c      integer(8) :: nhmer,nhmesum

      real :: ham,ham3b_cJ,ham2b_TUD,ham2b_pn
      real(kind(0.d0)) :: xcoef
      integer :: mysave,iix,iterest,iter,na,nb,nesd,iro,iphase
      integer,allocatable :: occi(:),occf(:),occim2(:),occim3(:)
      integer :: N_HO_i,wd,ibit,interm_energy,N_HO_f
      integer(8),allocatable :: intbasi(:),intbasf(:),intbasim3(:) !,
c     $     intbas(:)
      integer :: ia_1,an_st_1,cr_1,cr_st_1
      integer :: ia_2,an_st_2,cr_st_2,sp_st_an_2,sp_st_cr_2,
     $     cr_st_min_1,cr_st_max_1,cr_st_min_2,cr_st_max_2
      integer :: ia_3,an_st_3,cr_st_3,cr_st_min_3,cr_st_max_3
      type multiple_type
      integer :: dim
      integer,allocatable :: multiple(:,:)
      end type multiple_type
      type MJ_type
      type(multiple_type),allocatable :: MJ(:)
      end type MJ_type
      type(MJ_type),allocatable :: parity_MT(:,:)
      integer :: mshift(0:1),mj,pi,mt,num_of_mult,minz,nspsmin,
     $     N_1,N_2,N_3,mj_ind,mt_ind,N123,N123_max
      integer :: mult_nucl_dim,ij
      integer,allocatable :: mult_nucl(:,:)
      integer,allocatable :: multiple_dim(:,:,:),multiple_point(:,:,:),
     $     multiple_ist(:,:)
      integer :: iallred,iallredbuf,iiy,northg
      character(len=80) :: line,line1
      logical :: save47
      integer :: prot_st,neut_st,exminho
      real(kind(0.d0)) :: dnrm2, ddot
      integer :: Nmin_HO
c* PN
      integer :: nsdmax,nsd_hold,nsdp_hold,nsdn_hold,k1
      logical :: IT_test
      integer :: mnop_save
      real(kind(0.d0)) :: kappa_max

c     MKGK variables added
c      integer :: occx(nucleons)
c      integer :: wd,pos,k, index
c     do loop optimization
      integer :: N12_max_hold
      logical :: tbmebin
      real(kind(0.d0)) :: hbo_in,lambda_in

c     test variables
c      integer :: state_adder, state_adder_mpi
c      integer :: state_deleter, state_deleter_mpi

      ibuf=min(nsd,3000000)

c     ------------------------------------------------------------------------
c     ----- bmp() is used in TRANSITIONS to store the final eigenstates. 
c     ------------------------------------------------------------------------
c

      if (associated(ty%nucl)) deallocate(ty%nucl)
      if (associated(ty%protn)) deallocate(ty%protn)
      if (associated(ty%neutn)) deallocate(ty%neutn)
      if (associated(ty%dim_st)) deallocate(ty%dim_st)
      if (associated(ty%dim_st_p)) deallocate(ty%dim_st_p)
      if (associated(ty%dim_st_n)) deallocate(ty%dim_st_n)


      allocate(ty%nucl)
      allocate(ty%protn)
      allocate(ty%neutn)
      allocate(ty%dim_st)
      allocate(ty%dim_st_p)
      allocate(ty%dim_st_n)
      ty%nucl=nucleons
      ty%protn=nprotons
      ty%neutn=nneutrns
      ty%dim_st=nsd
      ty%dim_st_p=nsdp
      ty%dim_st_n=nsdn
      ty%occ=>iloc
      ty%occp=>ilocp
      ty%occn=>ilocn

c     Setup hold variables for do loop optimization
c      N12_max_hold = N12_max
c     N12_max is now correct for current nhw space (below)
c     No further mods need to be made in Lanczos_Hash
c      N12_max = nhw-(nhw_max-N12_max) ! CHECK! 

c     new parameters set for importance truncation
c* PN
c      nhw_hold = nhw-2
c      nhw = nhw + 2
      N12_max_hold = N12_max
      nspsmin=min(Nmin_HO(nprotons-2)+Nmin_HO(nneutrns),
     $     Nmin_HO(nprotons-1)+Nmin_HO(nneutrns-1),
     $     Nmin_HO(nprotons)+Nmin_HO(nneutrns-2))
      if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
      N12_max=nhw-nspsmin
      if (iproc==0) print *,' N12_max=',N12_max
c      N12_max = nhw - nhw0
      if (iproc==0) print*, 'in importance trunc nhw = ', nhw
c      print*, ' in " " N12_max = ',N12_max

      mnop_save=mnop
      if (mnop==-3) then
         mnop=2
      endif
      if (iproc==0) print *,' mnop,mnop_save=',mnop,mnop_save
      
c     set hash table info for adding states
      dim_constr_p=nsdp
      dim_constr_n=nsdn

      nsdmax=nsd

      if (allocated(delta_IT)) deallocate(delta_IT)
      allocate(delta_IT(nsd))
      delta_IT=999.99
      
      kappa_max=kappa_store(kappa_points)
      if (iproc==0) print *,' kappa,kappa_max=',kappa,kappa_max

c* PN: set beta_cm=0 for IT testing
      call MPI_Barrier(icomm,ierr)
      if (tbmeTUD) then
         if (iproc==0) then
            print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            print *,' n_Max,l_Max=',n_Max,l_Max 
         endif
         inquire(file=trim(intfile)//'_VNNbin',exist=tbmebin)
         if (tbmebin) then
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.true.)
            if (iproc==0) then
               print *,' read_TUD_tbme_bin called: hbo,lambda=',
     $              hbo_in,lambda_in
               print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            endif
            if (abs(hbo_in-hbomeg)>1.d-4) then
               if (iproc==0) print *,'***error: hbomeg,hbo_in=',
     $              hbomeg,hbo_in
               stop
            endif
c            if (N12_max_tud/=mxnn) then
c               if (iproc==0) print *,
c     $              '*** warning: N12_max,mxnn=',
c     $              N12_max_tud,mxnn
c            endif
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call TUD_tbme_setup
            if (iproc==0) print *,' TUD_tbme_setup called'
         else
            call TUD_tbme_setup
            if (iproc==0) print *,' TUD_tbme_setup called'
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call read_TUD_tbme(intfile)
            if (iproc==0) print *,' read_TUD_tbme called'
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.false.)
            if (iproc==0) print *,
     $           ' read_TUD_tbme_bin called for Trel, HOrel'
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call set_TUD_H(0.0)
         if (iproc==0) print *,' set_TUD_H called with beta_cm=0' 
c*** MPI
c         call MPI_Barrier(icomm,ierr)
c***  MPI
c         if (abs(mnop)/=3) then
c            call cg_init(N1_max_TUD,N12_max_tud,0)
c            if (iproc==0) print *,' cg_init called'
c         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
      else
         call readtbme(intfile,intform,nobt,0.0,iclmb,
     &        itrel,ihrel,ispe,nskip)
         call MPI_Barrier(icomm,ierr)      
         if (iproc==0) print *,' readtbme called with beta_cm=0'
      endif

c     how many states added
c      state_adder = 0

c     test code
c      if (iproc==0) write(6,*) 'N12_max=',N12_max


c     MKGK edit - this is called in make_hash_iloc
c      call construct_hash_table

      if (allocated(intbasi)) deallocate(intbasi)
      if (allocated(occi)) deallocate(occi)

      allocate(intbasi(2*nwd))
      allocate(occi(nucleons))
c      do i1=1,nsd
c         prot_st=iloc(1,i1)
c         neut_st=iloc(2,i1)
c         occi(1:nprotons)=ilocp(:,prot_st)
c         occi(nprotons+1:nucleons)=ilocn(:,neut_st)
c         intbasi=0
c         do cr_1=1,nucleons
c            wd=(occi(cr_1)-1)/nbit+1
c            ibit=mod(occi(cr_1)-1,nbit)
cc     print *,' i,wd,ibit=',i,wd,ibit
c            intbasi(wd)=ibset(intbasi(wd),ibit)
c         end do
c         call get_state_index(nucleons,2*nwd,intbasi,ionb)
c         if (i1/=ionb) then
c            print *,'*** error:'
c            if (iproc==0) then
c               print *,' i=',i1,'    index=',ionb
c               print *, iloc(:,i1)
c               print *,' occi=',occi
c            endif
c            stop
c         endif
c      end do

      if (mnop==3) then
         minz=max(0,nneutrns-2)+max(0,nprotons-2)
         nspsmin=max(minz-3,0)
         if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
         nspsmin=min(Nmin_HO(nprotons-3)+Nmin_HO(nneutrns),
     $        Nmin_HO(nprotons-2)+Nmin_HO(nneutrns-1),
     $        Nmin_HO(nprotons-1)+Nmin_HO(nneutrns-2),
     $        Nmin_HO(nprotons)+Nmin_HO(nneutrns-3))
         if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
         N123_max=nhw-nspsmin
         if (iproc==0) print *,' N123_max=',N123_max
         mult_nucl_dim=0
         do ia_1=1,nucleons-2
            do ia_2=ia_1+1,nucleons-1
               do ia_3=ia_2+1,nucleons
                  mult_nucl_dim=mult_nucl_dim+1
               end do
            end do
         end do
         if (iproc==0) print *,' mult_nucl_dim=',mult_nucl_dim
         allocate(mult_nucl(3,mult_nucl_dim))
         ij=0
         do ia_1=1,nucleons-2
            do ia_2=ia_1+1,nucleons-1
               do ia_3=ia_2+1,nucleons
                  ij=ij+1
                  mult_nucl(1,ij)=ia_1
                  mult_nucl(2,ij)=ia_2
                  mult_nucl(3,ij)=ia_3
               end do
            end do
         end do
         if (ij/=mult_nucl_dim) then
            print *,'*** error: ii,mult_nucl_dim=',
     $           ij,mult_nucl_dim
            stop
         endif

c         allocate(parity_MT(0:1,0:mnop))
         mshift=1
         do cr_st_1=1,nasps-2
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N123_max) exit
            do cr_st_2=cr_st_1+1,nasps-1
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N123_max) exit
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               do cr_st_3=cr_st_2+1,nasps
                  N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                  if (N_2+N_3>N12_max) exit
                  if (N_1+N_3>N12_max) exit
                  if (N_1+N_2+N_3>nhw-nspsmin) exit
c                  if (N_2+N_3>N12_max) cycle
c                  if (N_1+N_3>N12_max) cycle
c                  if (N_1+N_2+N_3>nhw-nspsmin) cycle
                  mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)+m2_sp(cr_st_3)
                  if (mj>mshift(1)) mshift(1)=mj
               end do
            end do
         end do
         do cr_st_1=1,nasps-1
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N123_max) exit
            do cr_st_2=cr_st_1+1,nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N123_max) exit
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               do cr_st_3=mxsps+1,mxsps+nasps
                  N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                  if (N_2+N_3>N12_max) exit
                  if (N_1+N_3>N12_max) exit
                  if (N_1+N_2+N_3>nhw-nspsmin) exit
c                  if (N_2+N_3>N12_max) cycle
c                  if (N_1+N_3>N12_max) cycle
c                  if (N_1+N_2+N_3>nhw-nspsmin) cycle
                  mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)+m2_sp(cr_st_3)
                  if (mj>mshift(0)) mshift(0)=mj
               end do
            end do
         end do
         if (iproc==0) print *,' mshift=',mshift

         if (iproc==0) print *,' maxval(mshift)=',maxval(mshift)
         allocate(multiple_dim(0:1,0:3,0:maxval(mshift)))
         multiple_dim=0
         allocate(multiple_point(0:1,0:3,0:maxval(mshift)))
         multiple_point=0

         ij=0
         do pi=0,1
            do mt=-3,3,2

               select case(mt)
               case(3)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_min_3=3
                  cr_st_max_1=nasps-2
                  cr_st_max_2=nasps-1
                  cr_st_max_3=nasps
               case(1)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
                  cr_st_min_3=mxsps+1
                  cr_st_max_3=mxsps+nasps
               case(-1)
                  cr_st_min_1=1
                  cr_st_min_2=mxsps+1
                  cr_st_min_3=mxsps+2
                  cr_st_max_1=nasps
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case(-3)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_min_3=mxsps+3
                  cr_st_max_1=mxsps+nasps-2
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
c               allocate(parity_MT(pi,(mt+mnop)/2)%MJ(
c     $              0:mshift(abs(mt)/2)))
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  multiple_point(pi,mt_ind,mj_ind)=ij
                  num_of_mult=0
                  do N123=0,N123_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123_max) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
                        if (N_1+N_2>N12_max) exit
c                        if (N_1+N_2>N12_max) cycle
                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
     $                       cr_st_max_3
                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                           if (N_1+N_2+N_3>N123) exit
                           if (N_2+N_3>N12_max) exit
                           if (N_1+N_3>N12_max) exit
c                           if (N_1+N_2+N_3>N123_max) exit
                           if (N_1+N_2+N_3/=N123) cycle
c                           if (N_2+N_3>N12_max) cycle
c                           if (N_1+N_3>N12_max) cycle
c                           if (N_1+N_2+N_3>N123_max) cycle
                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                          +mt2_sp(cr_st_3)/=mt) cycle
                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                          +m2_sp(cr_st_3)/=mj) cycle
                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                          +l_sp(cr_st_3),2)/=pi) cycle
                           num_of_mult=num_of_mult+1
                        end do
                     end do
                  end do
                  end do
c                  if (iproc==0) print *,' pi,mt,mj,num_of_mult=',
c     $                 pi,mt,mj,num_of_mult

c                  parity_MT(pi,mt_ind)%MJ(mj_ind)
c     $                 %dim=num_of_mult

                  ij=ij+num_of_mult
c                  if (num_of_mult>ij) ij=num_of_mult
                  multiple_dim(pi,mt_ind,mj_ind)=num_of_mult


                  if (num_of_mult==0) cycle
c                  allocate(parity_MT(pi,mt_ind)
c     $                 %MJ(mj_ind)
c     $                 %multiple(3,num_of_mult))
c                  num_of_mult=0
c                  do N123=0,N123_max
c                  do cr_st_1=cr_st_min_1,cr_st_max_1
c                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
c                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
c                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
c                        if (N_1+N_2>N12_max) cycle
c                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
c     $                       cr_st_max_3
c                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
c                           if (N_2+N_3>N12_max) cycle
c                           if (N_1+N_3>N12_max) cycle
c                           if (N_1+N_2+N_3>N123_max) cycle
c                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
c     $                          +mt2_sp(cr_st_3)/=mt) cycle
c                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
c     $                          +m2_sp(cr_st_3)/=mj) cycle
c                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
c     $                          +l_sp(cr_st_3),2)/=pi) cycle
c                           if (N_1+N_2+N_3/=N123) cycle
c                           num_of_mult=num_of_mult+1
cc                           if (iproc==0) print *,
cc     $                      ' cr_st_1,cr_st_2,cr_st_3,num_of_mult=',
cc     $                          cr_st_1,cr_st_2,cr_st_3,num_of_mult

cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(1,num_of_mult)
cc     $                          =cr_st_1
cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(2,num_of_mult)
cc     $                          =cr_st_2
cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(3,num_of_mult)
cc     $                          =cr_st_3
c                        end do
c                     end do
c                  end do
c                  end do
               end do
            end do
         end do
         
         if (iproc==0) print *,' total num_of_mult=',ij
         allocate(multiple_ist(3,ij))
         multiple_ist=0
         ij=0
         do pi=0,1
            do mt=-3,3,2

               select case(mt)
               case(3)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_min_3=3
                  cr_st_max_1=nasps-2
                  cr_st_max_2=nasps-1
                  cr_st_max_3=nasps
               case(1)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
                  cr_st_min_3=mxsps+1
                  cr_st_max_3=mxsps+nasps
               case(-1)
                  cr_st_min_1=1
                  cr_st_min_2=mxsps+1
                  cr_st_min_3=mxsps+2
                  cr_st_max_1=nasps
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case(-3)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_min_3=mxsps+3
                  cr_st_max_1=mxsps+nasps-2
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
c                  num_of_mult=0
                  do N123=0,N123_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
                        if (N_1+N_2>N12_max) exit
c                        if (N_1+N_2>N12_max) cycle
                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
     $                       cr_st_max_3
                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                           if (N_1+N_2+N_3>N123) exit
                           if (N_2+N_3>N12_max) exit
                           if (N_1+N_3>N12_max) exit
c                           if (N_1+N_2+N_3>N123_max) exit
                           if (N_1+N_2+N_3/=N123) cycle
c                           if (N_2+N_3>N12_max) cycle
c                           if (N_1+N_3>N12_max) cycle
c                           if (N_1+N_2+N_3>N123_max) cycle
                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                          +mt2_sp(cr_st_3)/=mt) cycle
                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                          +m2_sp(cr_st_3)/=mj) cycle
                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                          +l_sp(cr_st_3),2)/=pi) cycle
c                           num_of_mult=num_of_mult+1
                           ij=ij+1
                           multiple_ist(1,ij)
     $                          =cr_st_1
                           multiple_ist(2,ij)
     $                          =cr_st_2
                           multiple_ist(3,ij)
     $                          =cr_st_3
                        end do
                     end do
                  end do
                  end do
               end do
            end do
         end do

      elseif (mnop==2) then

         mult_nucl_dim=0
         do ia_1=1,nucleons-1
            do ia_2=ia_1+1,nucleons
               mult_nucl_dim=mult_nucl_dim+1
            end do
         end do
         if (iproc==0) print *,' mult_nucl_dim=',mult_nucl_dim
         allocate(mult_nucl(2,mult_nucl_dim))
         ij=0
         do ia_1=1,nucleons-1
            do ia_2=ia_1+1,nucleons
               ij=ij+1
               mult_nucl(1,ij)=ia_1
               mult_nucl(2,ij)=ia_2
            end do
         end do
         if (ij/=mult_nucl_dim) then
            print *,'*** error: ii,mult_nucl_dim=',
     $           ij,mult_nucl_dim
            stop
         endif

c         allocate(parity_MT(0:1,0:mnop))
         mshift=0
         do cr_st_1=1,nasps-1
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N12_max) exit
            do cr_st_2=cr_st_1+1,nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)
               if (mj>mshift(1)) mshift(1)=mj
            end do
         end do
         do cr_st_1=1,nasps
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N12_max) exit
            do cr_st_2=mxsps+1,mxsps+nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)
               if (mj>mshift(0)) mshift(0)=mj
            end do
         end do

         if (iproc==0) print *,' mshift=',mshift

         if (iproc==0) print *,' maxval(mshift)=',maxval(mshift)
         allocate(multiple_dim(0:1,0:2,0:maxval(mshift)))
         multiple_dim=0
         allocate(multiple_point(0:1,0:2,0:maxval(mshift)))
         multiple_point=0

         ij=0
         do pi=0,1
            do mt=-2,2,2

               select case(mt)
               case(2)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
               case(0)
                  cr_st_min_1=1
                  cr_st_max_1=nasps
                  cr_st_min_2=mxsps+1
                  cr_st_max_2=mxsps+nasps
               case(-2)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_max_1=mxsps+nasps-1
                  cr_st_max_2=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2

c               allocate(parity_MT(pi,(mt+mnop)/2)%MJ(
c     $              0:mshift(abs(mt)/2)))
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  multiple_point(pi,mt_ind,mj_ind)=ij
                  num_of_mult=0
                  do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
c                        if (N_1+N_2>N12_max) exit
                        if (N_1+N_2/=N123) cycle
c                        if (N_1+N_2>N12_max) cycle
                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                       /=mt) cycle
                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                       ,2)/=pi) cycle
                        num_of_mult=num_of_mult+1
                     end do
                  end do
                  end do
c                  parity_MT(pi,mt_ind)%MJ(mj_ind)
c     $                 %dim=num_of_mult

                  ij=ij+num_of_mult
c                  if (num_of_mult>ij) ij=num_of_mult
                  multiple_dim(pi,mt_ind,mj_ind)=num_of_mult

                  if (num_of_mult==0) cycle
c                  allocate(parity_MT(pi,mt_ind)
c     $                 %MJ(mj_ind)
c     $                 %multiple(2,num_of_mult))
c                  num_of_mult=0
c                  do N123=0,N12_max
c                  do cr_st_1=cr_st_min_1,cr_st_max_1
c                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
c                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
c                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
c                        if (N_1+N_2>N12_max) cycle
c                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
c     $                       /=mt) cycle
c                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
c     $                       /=mj) cycle
c                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
c     $                       ,2)/=pi) cycle
c                        if (N_1+N_2/=N123) cycle
c                        num_of_mult=num_of_mult+1
cc                        parity_MT(pi,mt_ind)
cc     $                       %MJ(mj_ind)
cc     $                       %multiple(1,num_of_mult)
cc     $                       =cr_st_1
cc                        parity_MT(pi,mt_ind)
cc     $                       %MJ(mj_ind)
cc     $                       %multiple(2,num_of_mult)
cc     $                       =cr_st_2
c                     end do
c                  end do
c                  end do
               end do
            end do
         end do
         
         if (iproc==0) print *,' total num_of_mult=',ij
         allocate(multiple_ist(2,ij))
         multiple_ist=0
         ij=0

         do pi=0,1
            do mt=-2,2,2
               select case(mt)
               case(2)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
               case(0)
                  cr_st_min_1=1
                  cr_st_max_1=nasps
                  cr_st_min_2=mxsps+1
                  cr_st_max_2=mxsps+nasps
               case(-2)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_max_1=mxsps+nasps-1
                  cr_st_max_2=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
c                        if (N_1+N_2>N12_max) exit
                        if (N_1+N_2/=N123) cycle
c                        if (N_1+N_2>N12_max) cycle
                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                       /=mt) cycle
                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                       ,2)/=pi) cycle
                        ij=ij+1
                        multiple_ist(1,ij)
     $                       =cr_st_1
                        multiple_ist(2,ij)
     $                       =cr_st_2
                     end do
                  end do
                  end do
               end do
            end do
         end do

      endif

      xxi=0.d0     
c     ----- Initialize amp:

      if (allocated(amp)) deallocate(amp)
c* PN
c*      allocate(amp(nsd,2))
      allocate(amp(nsd,nf))
      amp = 0.d0           ! New Lanczos vector
c      call pivot(amp(:,1))  ! Initial Lanczos vector

c     make the pivot call
      if (piv_saved) then
         if (iproc==0) print *,
     $        ' calling read_pivot from importance truncate'
         do i1=1,nf
            call read_pivot(nsd,amp(:,i1),i1)  ! PN: error for nf>2
         end do
      else 
      
         exminho=mod(nhw-nminho,2)
         if (nhw-nminho>=4) then
            x00=sqrt(0.65d0/real(nsd_hw(exminho),kind(0.d0)))
            amp(1:nsd_hw(exminho),1) = x00 
            x00=sqrt(0.25d0/real(nsd_hw(exminho+2)-nsd_hw(exminho),
     $           kind(0.d0)))
            do ij=nsd_hw(exminho)+1,nsd_hw(exminho+2)
               amp(ij,1) = x00 
            end do
            x00=sqrt(0.1d0/real(nsd_hw(exminho+4)-nsd_hw(exminho+2),
     $           kind(0.d0)))
            do ij=nsd_hw(exminho+2)+1,nsd_hw(exminho+4)
               amp(ij,1) = x00 
            end do
         elseif (nhw-nminho>=2) then
            x00=sqrt(0.65d0/real(nsd_hw(exminho),kind(0.d0)))
            amp(1:nsd_hw(exminho),1) = x00 
            x00=sqrt(0.35d0/real(nsd_hw(exminho+2)-nsd_hw(exminho),
     $           kind(0.d0)))
            do ij=nsd_hw(exminho)+1,nsd_hw(exminho+2)
               amp(ij,1) = x00 
            end do
         else
            x00=sqrt(1.d0/real(nsd,kind(0.d0)))
            do ij=1,nsd
               amp(ij,1) = x00*(-1)**ij
            end do
         endif
      end if ! piv_saved or not

c     check norm of new pivot
      x00=0.d0
c      print *,' x00=',x00
      do ij=1,nsd
         x00=x00+amp(ij,1)*amp(ij,1)
c         print *,' i1,amp(i1,1)=',i1,amp(i1,1),x00
      end do
      if (iproc==0) print *,' pivot norm=',x00

c     MKGK
c     nnmax = min(nite,nsd)     ! Number of Lanczos iterations
c     for Importance truncation 1 Lanczos iter is done
      nnmax = 1
      mysave=0

c     MKGK
      if (allocated(alpha)) deallocate(alpha)
      if (allocated(beta)) deallocate(beta)
      allocate(alpha(nnmax),beta(nnmax))
      alpha=0.d0
      beta=0.d0
      
c-------------------------------------------------------------------------
c     ----- THE MAIN LOOP: LANCZOS ITERATIONS
c     ----- The matrix elements of H are EITHER read from the files 
c     ----- mfdp.tmp1 and mfdp.tmp2 OR calculated in the first iteration.


      if (allocated(intbasf)) deallocate(intbasf)
cc      if (allocated(intbas)) deallocate(intbas)
      if (allocated(intbasim3)) deallocate(intbasim3)
      if (allocated(occim2)) deallocate(occim2)
      if (allocated(occim3)) deallocate(occim3)
      if (allocated(occf)) deallocate(occf)
      if (allocated(energy)) deallocate(energy)

      allocate(intbasf(2*nwd))
cc      allocate(intbas(2*nwd))
      allocate(intbasim3(2*nwd))
ccc   allocate(locz(2*mxsps))
      allocate(occim2(nucleons-2))
      allocate(occim3(nucleons-3))
      allocate(occf(nucleons))

      if (.not.allocated(locz)) allocate(locz(2*mxsps))
ccc      if (.not.allocated(intbas)) allocate(intbas(2*nwd))

c     MKGK
c     allocate(energy(nnmax))

c     set parameters
      na = 1
      nb = 2
      iter = 1

c***  Standard calculation ***********************************

      do 5000 i1 = iproc+1, nsd, nproc
cc      do 5000 i1 = 1, nsd
         prot_st=iloc(1,i1)
         neut_st=iloc(2,i1)
         occi(1:nprotons)=ilocp(:,prot_st)
         occi(nprotons+1:nucleons)=ilocn(:,neut_st)
c     occi(:)=iloc(:,i1)
         intbasi=0
         do cr_1=1,nucleons
            wd=(occi(cr_1)-1)/nbit+1
            ibit=mod(occi(cr_1)-1,nbit)
c     print *,' i,wd,ibit=',i,wd,ibit
            intbasi(wd)=ibset(intbasi(wd),ibit)
         end do
         locz(1:occi(1))=0
         do cr_1=2,nucleons
            locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
         end do
         locz(occi(nucleons)+1:2*mxsps)=nucleons


         N_HO_i=0
         do ia_1=1,nucleons
            an_st_1=occi(ia_1)
            N_HO_i=N_HO_i+2*n_sp(an_st_1)+l_sp(an_st_1)
         end do
         nesd = N_HO_i - nespmin
         xamp = amp(i1,na)

         if (i1.eq.10000*(i1/10000).and.nproc==1) then 
            print *, ' doing i1=',i1 
         endif
c     *pn* 
         if (mnop==3) then

c     print *,' i1',i1
c     print *,' occi=',occi

c     occf=occi ! test
c     locz=0 !test

            do ij=1,mult_nucl_dim

               ia_1=mult_nucl(1,ij)
               ia_2=mult_nucl(2,ij)
               ia_3=mult_nucl(3,ij)

c     do ia_1=1,nucleons-2
               an_st_1=occi(ia_1)
c     do ia_2=ia_1+1,nucleons-1
               an_st_2=occi(ia_2)
c     do ia_3=ia_2+1,nucleons
               an_st_3=occi(ia_3)

               intbasim3=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)

               mt=mt2_sp(an_st_1)
     $              +mt2_sp(an_st_2)+mt2_sp(an_st_3)
               mj=m2_sp(an_st_1)
     $              +m2_sp(an_st_2)+m2_sp(an_st_3)
               pi=mod(l_sp(an_st_1)+l_sp(an_st_2)
     $              +l_sp(an_st_3),2)

c     occim3(1:ia_1-1)=occi(1:ia_1-1)
c     occim3(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
c     occim3(ia_2-1:ia_3-3)=occi(ia_2+1:ia_3-1)
c     occim3(ia_3-2:nucleons-3)=occi(ia_3+1:nucleons)


c     print *,' ia_1,ia_2,ia_3=',ia_1,ia_2,ia_3
c     print *,' an_st_1,an_st_2,an_st_3=',
c     $                       an_st_1,an_st_2,an_st_3
c     print *,' occim3=',occim3

c     locz(1:occim3(1))=0
c     do cr_1=2,nucleons-3
c     locz(occim3(cr_1-1)+1:occim3(cr_1))=cr_1-1
c     end do
c     locz(occim3(nucleons-3)+1:2*mxsps)=nucleons-3

c     intbas=0
c     do cr_1=1,nucleons-3
c     wd=(occim3(cr_1)-1)/nbit+1
c     ibit=mod(occim3(cr_1)-1,nbit)
c     c     print *,' i,wd,ibit=',i,wd,ibit
c     intbas(wd)=ibset(intbas(wd),ibit)
c     end do
               interm_energy=-2*n_sp(an_st_1)-l_sp(an_st_1)
     $              -2*n_sp(an_st_2)-l_sp(an_st_2)
     $              -2*n_sp(an_st_3)-l_sp(an_st_3)+N_HO_i


               mt_ind=(mt+mnop)/2
               mj_ind=(mshift(abs(mt)/2)+mj)/2

c     print *,' interm_energy=',interm_energy
c     print *,' pi,mt,mj=',pi,mt,mj
c     print *,' mt_ind,mj_ind=',mt_ind,mj_ind
c     print *,' parity_MT(pi,mt_ind)%MJ(mj_ind)%dim=',
c     $                       parity_MT(pi,mt_ind)%MJ(mj_ind)%dim

c     do num_of_mult=1,parity_MT(pi,mt_ind)
c     $                       %MJ(mj_ind)%dim
c     cr_st_1=parity_MT(pi,mt_ind)
c     $                          %MJ(mj_ind)
c     $                          %multiple(1,num_of_mult)
c     cr_st_2=parity_MT(pi,mt_ind)
c     $                          %MJ(mj_ind)
c     $                          %multiple(2,num_of_mult)
c     cr_st_3=parity_MT(pi,mt_ind)
c     $                          %MJ(mj_ind)
c     $                          %multiple(3,num_of_mult)

               cr_1=multiple_point(pi,mt_ind,mj_ind)

               do num_of_mult=1,multiple_dim(pi,mt_ind,mj_ind)

                  cr_1=cr_1+1
c     cr_1=num_of_mult
c     $                          +multiple_point(pi,mt_ind,mj_ind)
                  cr_st_1=multiple_ist(1,cr_1)
                  cr_st_2=multiple_ist(2,cr_1)
                  cr_st_3=multiple_ist(3,cr_1)

c     if (2*n_sp(cr_st_1)+l_sp(cr_st_1)
c     +                          +interm_energy>nhw) exit

                  N_HO_f=2*n_sp(cr_st_1)+l_sp(cr_st_1)
     $                 +2*n_sp(cr_st_2)+l_sp(cr_st_2)
     $                 +2*n_sp(cr_st_3)+l_sp(cr_st_3)
     +                 +interm_energy

c     MKGK uncommented line below, commented next exit line
                  if (N_HO_f>nhw) cycle !exit
c                  if (N_HO_f>N_HO_i) exit
                  
                  intbasf=intbasim3
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbasf(wd),ibit)) cycle
                  intbasf(wd)=ibset(intbasf(wd),ibit)
                  wd=(cr_st_2-1)/nbit+1
                  ibit=mod(cr_st_2-1,nbit)
                  if (btest(intbasf(wd),ibit)) cycle
                  intbasf(wd)=ibset(intbasf(wd),ibit)
                  wd=(cr_st_3-1)/nbit+1
                  ibit=mod(cr_st_3-1,nbit)
                  if (btest(intbasf(wd),ibit)) cycle
                  intbasf(wd)=ibset(intbasf(wd),ibit)

c     occf(1:locz(cr_st_1))=
c     $                          occim3(1:locz(cr_st_1))
c     occf(locz(cr_st_1)+1)=cr_st_1
c     occf(locz(cr_st_1)+2:locz(cr_st_2)+1)=
c     $                          occim3(locz(cr_st_1)+1
c     $                          :locz(cr_st_2))
c     occf(locz(cr_st_2)+2)=cr_st_2
c     occf(locz(cr_st_2)+3:locz(cr_st_3)+2)=
c     $                          occim3(locz(cr_st_2)+1
c     $                          :locz(cr_st_3))
c     occf(locz(cr_st_3)+3)=cr_st_3
c     occf(locz(cr_st_3)+4:nucleons)=
c     +                          occim3(locz(cr_st_3)+1:nucleons-3)
                  

c     print *,' cr_st_1,cr_st_2,cr_st_3=',
c     $                                cr_st_1,cr_st_2,cr_st_3
c     print *,' occf=',occf

                  call get_state_index(nucleons,2*nwd,intbasf,
     $                 i3)
c     i3=1 !test

c     print *,' i3=',i3

c* PN
c                  if (i3==-1) cycle
                  if (i3 > 0) then
                     if (delta_IT(i3)>kappa_max) cycle
                  endif
c                  state_adder_mpi = state_adder_mpi+1

c                  iphase=ia_1+ia_2+ia_3
c     $                 +locz(cr_st_1)+locz(cr_st_2)
c     $                 +locz(cr_st_3)
c                  if (cr_st_1>an_st_1) iphase=iphase+1
c                  if (cr_st_1>an_st_2) iphase=iphase+1
c                  if (cr_st_1>an_st_3) iphase=iphase+1
c                  if (cr_st_2>an_st_1) iphase=iphase+1
c                  if (cr_st_2>an_st_2) iphase=iphase+1
c                  if (cr_st_2>an_st_3) iphase=iphase+1
c                  if (cr_st_3>an_st_1) iphase=iphase+1
c                  if (cr_st_3>an_st_2) iphase=iphase+1
c                  if (cr_st_3>an_st_3) iphase=iphase+1

c                  xcoef=real((-1)**iphase,kind(0.d0))
c     $                 *ham3b_cJ(cr_st_1,cr_st_2,cr_st_3,
c     $                 an_st_1,an_st_2,an_st_3)
                  xcoef=ham3b_cJ(cr_st_1,cr_st_2,cr_st_3,
     $                 an_st_1,an_st_2,an_st_3)

                  IT_test=.false.
                  deltaIT=0.d0
                  do k1=1,nf
                     deltaIT=max(deltaIT,
     $                    abs(xcoef*amp(i1,k1))/(2.d0*hbomeg))
c                     if (abs(xcoef*amp(i1,k1))/(2.d0*hbomeg)>kappa)then
c                        exit
                  end do
                  if (deltaIT>kappa)then
                     IT_test=.true.
                  endif

c     add state to hash table if kappa satisfied
                  if (IT_test) then
c                     state_adder = state_adder+1

c     test code
c                     print*,'intbasi=',intbasi
c                     print*,'intbasf=',intbasf

                     call insert_hash_table_entry_word(nucleons,intbasi,
     $                    intbasf,2*nwd,i3)
                     if (i3==-1) print*, 'i3 bummer in insert_hash;3820'
                     nsdmax=max(nsdmax,i3)
c                     print *,' #1:N_HO_i,N_HO_f,i3,nsdmax=',
c     $                    N_HO_i,N_HO_f,i3,nsdmax
                  end if

c     print *,' xcoef'
c     MKGK - don't uncomment line below (xamp)
c     xamp = amp(i1,na)
c                  amp(i3,nb) = amp(i3,nb) + xcoef * xamp

c                  if (N_HO_f<N_HO_i) amp(i1,nb) = amp(i1,nb)
c     $                 +xcoef*amp(i3,na)

c     ipht=locz(icrstate)+locz(icrstatep)+locz(icrstatepp)
c     +              +locz(ianstatepp)
c     +              +locz(ianstate)+locz(ianstatep)-1
c     if (icrstate.gt.ianstate) ipht=ipht-1
c     if (icrstatep.gt.ianstate) ipht=ipht-1
c     if (icrstatepp.gt.ianstate) ipht=ipht-1
c     if (icrstate.gt.ianstatep) ipht=ipht-1
c     if (icrstatep.gt.ianstatep) ipht=ipht-1
c     if (icrstatepp.gt.ianstatep) ipht=ipht-1
c     if (icrstate.gt.ianstatepp) ipht=ipht-1
c     if (icrstatep.gt.ianstatepp) ipht=ipht-1
c     if (icrstatepp.gt.ianstatepp) ipht=ipht-1
c     iphase=(-1)**ipht
c     xcoef=real(iphase)*ham3b(icrstate,icrstatep,
c     +              icrstatepp,ianstate,ianstatep,ianstatepp)
c**   xcoef=0.0
*************************************************************
               end do
c     end do
c     end do
            end do

c     if mnop==2
         else

            do ij=1,mult_nucl_dim

               ia_1=mult_nucl(1,ij)
               ia_2=mult_nucl(2,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)

               intbasim3=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)

               mt=mt2_sp(an_st_1)
     $              +mt2_sp(an_st_2)
               mj=m2_sp(an_st_1)
     $              +m2_sp(an_st_2)
               pi=mod(l_sp(an_st_1)+l_sp(an_st_2),2)

               interm_energy=-2*n_sp(an_st_1)-l_sp(an_st_1)
     $              -2*n_sp(an_st_2)-l_sp(an_st_2)+N_HO_i
               
               mt_ind=(mt+mnop)/2
               mj_ind=(mshift(abs(mt)/2)+mj)/2

               cr_1=multiple_point(pi,mt_ind,mj_ind)

               do num_of_mult=1,multiple_dim(pi,mt_ind,mj_ind)
                  
                  cr_1=cr_1+1
                  
                  cr_st_1=multiple_ist(1,cr_1)
                  cr_st_2=multiple_ist(2,cr_1)

                  N_HO_f=2*n_sp(cr_st_1)+l_sp(cr_st_1)
     $                 +2*n_sp(cr_st_2)+l_sp(cr_st_2)
     +                 +interm_energy

c     MKGK uncommented line below, commented next exit line
                  if (N_HO_f>nhw) cycle !exit
c                  if (N_HO_f>N_HO_i) exit

                  intbasf=intbasim3
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbasf(wd),ibit)) cycle
                  intbasf(wd)=ibset(intbasf(wd),ibit)
                  wd=(cr_st_2-1)/nbit+1
                  ibit=mod(cr_st_2-1,nbit)
                  if (btest(intbasf(wd),ibit)) cycle
                  intbasf(wd)=ibset(intbasf(wd),ibit)

                  call get_state_index(nucleons,2*nwd,intbasf,
     $                 i3)

c                  print *,' first call: i3=',i3

c     DO IMPORTANCE TRUNCATION HERE
c     cycle if S.D. already exists in hash table
c                  state_adder_mpi = state_adder_mpi+1

                  if (i3 > 0) then
                     if (delta_IT(i3)>kappa_max) cycle
                  endif
c                  if (i3 > 0) cycle
c                  state_adder_mpi = state_adder_mpi+1

c* PN
c                  iphase=ia_1+ia_2
c     $                 +locz(cr_st_1)+locz(cr_st_2)+1
c                  if (cr_st_1>an_st_1) iphase=iphase+1
c                  if (cr_st_1>an_st_2) iphase=iphase+1
c                  if (cr_st_2>an_st_1) iphase=iphase+1
c                  if (cr_st_2>an_st_2) iphase=iphase+1

c     iphase=(-1)**(ia_1+ia_2+locz(cr_st_1)
c     $                          +locz(cr_st_2)+1)
                  kk=an_st_1
                  ll=an_st_2
                  ii=cr_st_1
                  jj=cr_st_2
************************************************************
c* PN                  
c                  xcoef=real((-1)**iphase)*ham(nesd)

                  if (tbmepn) then 
                     xcoef=ham2b_pn(ii,jj,kk,ll)
                  elseif (tbmeTUD) then 
                     xcoef=ham2b_TUD(ii,jj,kk,ll)
                  else
                     xcoef=ham(ii,jj,kk,ll,nesd)
                  endif
***********************************************************
c* PN
                  IT_test=.false.
                  deltaIT=0.d0
                  do k1=1,nf
                     deltaIT=max(deltaIT,
     $                    abs(xcoef*amp(i1,k1))/(2.d0*hbomeg))
c                     if (dabs(xcoef*amp(i1,k1))/(2.d0*hbomeg)>kappa)then
c                        exit
                  end do
                  if (deltaIT>kappa)then
                     IT_test=.true.
                  endif

c     add state to hash table if kappa satisfied
                  if (IT_test) then
c                     state_adder = state_adder+1

c     test code
c                     print*,'intbasi=',intbasi
c                     print*,'intbasf=',intbasf

                     call insert_hash_table_entry_word(nucleons,intbasi,
     $                    intbasf,2*nwd,i3)
                     if (i3==-1) print*, 'i3 bummer in insert_hash;3820'
                     nsdmax=max(nsdmax,i3)
c                     print *,' #1:N_HO_i,N_HO_f,i3,nsdmax=',
c     $                    N_HO_i,N_HO_f,i3,nsdmax
                  end if

c     safety check for next if statement N_HO_f < ...
c* PN: Why the following block??
c                  if (i3 > nsd) cycle
c
c                  if ((N_HO_f<N_HO_i).and.
c     $                 dabs(xcoef*amp(i3,na)/(2.d0*hbomeg))>kappa) then
c                     state_adder = state_adder+1
c                     call insert_hash_table_entry_word(nucleons,intbasi,
c     $                    intbasf,2*nwd,i3)
c                     if (i3==-1) print*, 'i3 bummer in insert_hash;3820'
c                     nsdmax=max(nsd,i3)
cc                     print *,' #2:N_HO_i,N_HO_f,i3,nsdmax=',
cc     $                    N_HO_i,N_HO_f,i3,nsdmax
c                  end if

c     MKGK : Remember both statements below (old code)
c                  amp(i3,nb) = amp(i3,nb) + xcoef * xamp
c                  if (N_HO_f<N_HO_i) amp(i1,nb) = amp(i1,nb)
c     $                 +xcoef*amp(i3,na)

               end do
            end do
         endif
 5000 continue

      deallocate(amp)
      call MPI_Barrier(icomm,ierr)
      if (nproc>1) then
         call reconcile_iloc(nsdmax,dim_constr_p,dim_constr_n)
      endif

c      print*, 'state_adder = ',state_adder
c      print*, 'state_adder_mpi = ',state_adder_mpi
c* PN
      if (iproc==0) then
         write(6,"(' Importance Trunaction: kappa_min=',d10.3)") kappa
         write(6,*) ' from nhw=',nhw_hold
         print *,' original nsd=',nsd
         print *,' original nsdp=',nsdp
         print *,' original nsdn=',nsdn
      endif
      nsd_hold=nsd
      nsdp_hold=nsdp
      nsdn_hold=nsdn
      nsd=nsdmax
      nsdp=dim_constr_p
      nsdn=dim_constr_n
      if (iproc==0) then
         write(6,*) ' nhw=',nhw
         print *,' constructed nsd=',nsd
         print *,' constructed nsdp=',nsdp
         print *,' constructed nsdn=',nsdn
      endif

      if (iproc==0) then
         write(9,*)
         write(9,"(' Importance Trunaction: kappa_min=',d10.3)") kappa
         write(9,*) ' from nhw=',nhw_hold
         write(9,*) ' original nsd=',nsd_hold
         write(9,*) ' original nsdp=',nsdp_hold
         write(9,*) ' original nsdn=',nsdn_hold
         write(9,*)
         write(9,*) ' nhw=',nhw
         write(9,*) ' constructed nsd=',nsd
         write(9,*) ' constructed nsdp=',nsdp
         write(9,*)' constructed nsdn=',nsdn
      endif

      call MPI_Barrier(icomm,ierr)

      call basis_N_HO_reorder(nsd_hold)

c* PN: restore beta_cm for Lanczos
      call MPI_Barrier(icomm,ierr)
      if (tbmeTUD) then
         if (iproc==0) then
            print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            print *,' n_Max,l_Max=',n_Max,l_Max 
         endif
         inquire(file=trim(intfile)//'_VNNbin',exist=tbmebin)
         if (tbmebin) then
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.true.)
            if (iproc==0) then
               print *,' read_TUD_tbme_bin called: hbo,lambda=',
     $              hbo_in,lambda_in
               print *,' N1_max,N12_max=',N1_max_tud,N12_max_tud
            endif
            if (abs(hbo_in-hbomeg)>1.d-4) then
               if (iproc==0) print *,'***error: hbomeg,hbo_in=',
     $              hbomeg,hbo_in
               stop
            endif
c            if (N12_max_tud/=mxnn) then
c               if (iproc==0) print *,
c     $              '*** warning: N12_max,mxnn=',
c     $              N12_max_tud,mxnn
c            endif
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call TUD_tbme_setup
            if (iproc==0) print *,' TUD_tbme_setup called'
         else
            call TUD_tbme_setup
            if (iproc==0) print *,' TUD_tbme_setup called'
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call read_TUD_tbme(intfile)
            if (iproc==0) print *,' read_TUD_tbme called'
c***  MPI
            call MPI_Barrier(icomm,ierr)
c*** MPI
            call read_TUD_tbme_bin(intfile,hbo_in,lambda_in,.false.)
            if (iproc==0) print *,
     $           ' read_TUD_tbme_bin called for Trel, HOrel'
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call set_TUD_H(strcm)
         if (iproc==0) print *,' set_TUD_H called with original beta_cm' 
c*** MPI
c         call MPI_Barrier(icomm,ierr)
c***  MPI
c         if (abs(mnop)/=3) then
c            call cg_init(N1_max_TUD,N12_max_tud,0)
c            if (iproc==0) print *,' cg_init called'
c         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
      else
         call readtbme(intfile,intform,nobt,strcm,iclmb,
     &        itrel,ihrel,ispe,nskip)
         call MPI_Barrier(icomm,ierr)      
         if (iproc==0) print *,' readtbme called with original beta_cm'
      endif

      mnop=mnop_save
      
c     MKGK
c     reset values for Lanczos_hash
c      nhw = nhw_hold
c      N12_max = N12_max_hold
      return
      end subroutine importance_truncate

      subroutine reconcile_iloc(nsdf,nsdpf,nsdnf)
      use spb
      use ibas
      use occ
      use hash_tables
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer,intent(INOUT) :: nsdf,nsdpf,nsdnf
      integer :: i,nsdsend(3),ierr,tag1,tag2,j,i2,nsdmax(3),
     $     nsdmaxtmp(3),stat,ise,ip,in,tag3,tag4,tag5
      integer(2),allocatable :: ilocp_tmp(:,:),ilocn_tmp(:,:),
     $     ilocN_tmp_buf1(:)
      integer(4),allocatable :: iloc_tmp(:,:),iloc_tmp_buf1(:),occ0(:)
      integer :: ibuf,iallredbuf,iiy,nucl,iallred,k
      integer :: status(MPI_STATUS_SIZE) 
      integer,parameter :: iunitvout=6
      logical :: isent
      real(kind(0.0)),allocatable :: delta_IT_tmp(:),delta_IT_tmp_buf(:)

      nsdmax(1)=nsdf
      nsdmax(2)=nsdpf
      nsdmax(3)=nsdnf
      if (iproc==0) then
         write(iunitvout,
     $        "(/,' iloc reconciliation starting from nsdf=',3i10)")
     $        nsdmax
      endif
      allocate(occ0(nucleons))
      k=log(real(nproc,kind(0.d0)))/log(2.d0)
      if (2**k<nproc) k=k+1
      if (iproc==0) then
         print *,' nproc,k,2**k=',nproc,k,2**k
      endif
      isent=.false.
      call MPI_Barrier(icomm,ierr)
      do i=0,k-1
c         if (mod(iproc,2**i)/=0) cycle
         if (mod(iproc/(2**i),2)==1.and.mod(iproc,2**i)==0) then
            nsdsend=nsdmax
            ise=iproc-2**i
            tag1=(iproc+1)  !*(2**i)
c            print *,' i,iproc,tag1=',i,iproc,tag1
            call MPI_Send(nsdsend,3,MPI_INTEGER,
     +                    ise,tag1,icomm,ierr)
c            print *,' sent: i,iproc,ise,tag1,nsdsend=',
c     $          i,iproc,ise,tag1,nsdsend

            tag2=(nproc+iproc+2)  !*(2**i)
c            print *,' i,iproc,tag2,num=',i,iproc,tag2,nsdf*2
            call MPI_Send(iloc(1,1),nsdsend(1)*2,MPI_INTEGER,
     +                    ise,tag2,icomm,ierr)
c            print *,' sent iloc: i,iproc,ise,tag2=',i,iproc,ise,tag2

            tag3=(2*nproc+iproc+3)  !*(2**i)
c            print *,' i,iproc,tag3,num=',i,iproc,tag3,nsdpf*nprotons
            call MPI_Send(ilocp(1,1),nsdsend(2)*nprotons,MPI_INTEGER2,
     +                    ise,tag3,icomm,ierr)
c            print *,' sent ilocp: i,iproc,ise,tag3=',i,iproc,ise,tag3

            tag4=(3*nproc+iproc+4)  !*(2**i)
c            print *,' i,iproc,tag4,num=',i,iproc,tag4,nsdnf*nneutrns
            call MPI_Send(ilocn(1,1),nsdsend(3)*nneutrns,MPI_INTEGER2,
     +                    ise,tag4,icomm,ierr)
c            print *,' sent ilocn: i,iproc,ise,tag4=',i,iproc,ise,tag4

            tag5=(4*nproc+iproc+5)
            call MPI_Send(delta_IT(1),nsdsend(1),MPI_REAL,
     +                    ise,tag5,icomm,ierr)

            isent=.true.
         endif

         if (mod(iproc/(2**i),2)==0.and.iproc+2**i<=nproc-1
     $        .and.mod(iproc,2**i)==0) then
            ise=iproc+2**i
            tag1=(ise+1)  !*(2**i) !MPI_ANY_TAG 
c            print *,' i,iproc=',i,iproc
            call MPI_Recv(nsdsend,3,MPI_INTEGER,
     $           ise,tag1,icomm,stat,ierr)
c            print *,' received: i,iproc,ise,tag1,nsdsend=',
c     $           i,iproc,ise,tag1,nsdsend
            allocate(iloc_tmp(2,nsdsend(1)))
c            print *,' iloc_tmp allocated: i,iproc,num=',i,iproc,
c     $           2*nsdsend(1)
            tag2=(nproc+ise+2)  !*(2**i) !MPI_ANY_TAG 
            call MPI_Recv(iloc_tmp(1,1),nsdsend(1)*2,MPI_INTEGER,
     +                    ise,tag2,icomm,stat,ierr)
c            print *, ' received iloc_tmp: i,iproc,ise,tag2=',
c     $           i,iproc,ise,tag2

            allocate(ilocp_tmp(nprotons,nsdsend(2)))
c            print *,' ilocp_tmp allocated: i,iproc,num=',i,iproc,
c     $           nprotons*nsdsend(2)
            tag3=(2*nproc+ise+3)  !*(2**i) !MPI_ANY_TAG 
            call MPI_Recv(ilocp_tmp(1,1),nsdsend(2)*nprotons,
     $           MPI_INTEGER2,ise,tag3,icomm,stat,ierr)
c            print *, ' received ilocp_tmp: i,iproc,ise,tag3=',
c     $           i,iproc,ise,tag3

            allocate(ilocn_tmp(nneutrns,nsdsend(3)))
c            print *,' ilocn_tmp allocated: i,iproc,num=',i,iproc,
c     $           nneutrns*nsdsend(3)
            tag4=(3*nproc+ise+4)  !*(2**i) !MPI_ANY_TAG 
            call MPI_Recv(ilocn_tmp(1,1),nsdsend(3)*nneutrns,
     $           MPI_INTEGER2,ise,tag4,icomm,stat,ierr)
c            print *, ' received ilocn_tmp: i,iproc,ise,tag4=',
c     $           i,iproc,ise,tag4

            allocate(delta_IT_tmp(nsdsend(1)))
            tag5=(4*nproc+ise+5)  !*(2**i) !MPI_ANY_TAG 
            call MPI_Recv(delta_IT_tmp(1),nsdsend(1),MPI_REAL,
     +                    ise,tag5,icomm,stat,ierr)

            do j=1,nsdsend(1)
               ip=iloc_tmp(1,j)
               in=iloc_tmp(2,j)
               occ0(1:nprotons)=ilocp_tmp(1:nprotons,ip)
               occ0(nprotons+1:nucleons)=ilocn_tmp(1:nneutrns,in)
               deltaIT=delta_IT_tmp(j)
               call insert_hash_table_entry(nucleons,
     $              occ0,1,i2)
               nsdmax(1)=max(nsdmax(1),i2)
               nsdmax(2)=max(nsdmax(2),dim_constr_p)
               nsdmax(3)=max(nsdmax(3),dim_constr_n)
            end do

c            nsdf=nsdmax(1)
c            nsdpf=nsdmax(2)
c            nsdnf=nsdmax(3)

         endif

         if (allocated(iloc_tmp)) deallocate(iloc_tmp)
         if (allocated(ilocp_tmp)) deallocate(ilocp_tmp)
         if (allocated(ilocn_tmp)) deallocate(ilocn_tmp)
         if (allocated(delta_IT_tmp)) deallocate(delta_IT_tmp)
         call MPI_Barrier(icomm,ierr)
         if (iproc==0) then
            write(iunitvout,"(' iproc=0 after i=',i4,'   nsdf=',3i10)")
     $           i,nsdmax
         endif
      end do

      deallocate(occ0)

      if (.not.isent) then
         print *,' nothing sent from iproc=',iproc
      endif

      if (iproc==0) then
         nsdmaxtmp=nsdmax
      endif
      call MPI_Barrier(icomm,ierr)
      call MPI_Bcast(nsdmaxtmp,3,MPI_INTEGER,0,icomm,ierr)
      nsdf=nsdmaxtmp(1)
      nsdpf=nsdmaxtmp(2)
      nsdnf=nsdmaxtmp(3)
      if (iproc/=0) then
         deallocate(iloc)
         deallocate(ilocp)
         deallocate(ilocn)
         deallocate(delta_IT)
         allocate(iloc(2,nsdf))
         allocate(ilocp(nprotons,nsdpf))
         allocate(ilocn(nneutrns,nsdnf))
         allocate(delta_IT(nsdf))
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdf*2,3000000)
      if (nsdf*2<=ibuf) then
         call MPI_Bcast(iloc(1,1),nsdf*2,MPI_INTEGER,0,
     $        icomm,ierr)
      else
         if (allocated(iloc_tmp_buf1)) deallocate(iloc_tmp_buf1)
         allocate(iloc_tmp_buf1(ibuf))
         do nucl=1,2
            iallredbuf=nsdf/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  iloc_tmp_buf1(:)=iloc(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(iloc_tmp_buf1,ibuf,MPI_INTEGER,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  iloc(nucl,iiy:iiy+ibuf-1)=iloc_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsdf,ibuf)
            if (iallred/=0) then
               if (iallred/=nsdf-iiy+1) then
                  print *,'***error:iallred,nsd,iiy:',iallred,nsdf,iiy
               endif
               if (iproc==0) then
                  iloc_tmp_buf1(1:iallred)=iloc(nucl,iiy:nsdf)
               endif
               call MPI_Bcast(iloc_tmp_buf1,iallred,MPI_INTEGER,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  iloc(nucl,iiy:nsdf)=iloc_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(iloc_tmp_buf1)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdpf*nprotons,3000000)
      if (nsdpf*nprotons<=ibuf) then
         call MPI_Bcast(ilocp(1,1),nsdpf*nprotons,MPI_INTEGER2,0,
     $        icomm,ierr)
      else
         if (allocated(ilocN_tmp_buf1)) deallocate(ilocN_tmp_buf1)
         allocate(ilocN_tmp_buf1(ibuf))
         do nucl=1,nprotons
            iallredbuf=nsdpf/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  ilocN_tmp_buf1(:)=ilocp(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,ibuf,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocp(nucl,iiy:iiy+ibuf-1)=ilocN_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsdpf,ibuf)
            if (iallred/=0) then
               if (iallred/=nsdpf-iiy+1) then
                  print *,'***error:iallred,nsdp,iiy:',iallred,nsdpf,iiy
               endif
               if (iproc==0) then
                  ilocN_tmp_buf1(1:iallred)=ilocp(nucl,iiy:nsdpf)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,iallred,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocp(nucl,iiy:nsdpf)=ilocN_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(ilocN_tmp_buf1)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdnf*nneutrns,3000000)
      if (nsdnf*nneutrns<=ibuf) then
         call MPI_Bcast(ilocn(1,1),nsdnf*nneutrns,MPI_INTEGER2,0,
     $        icomm,ierr)
      else
         if (allocated(ilocN_tmp_buf1)) deallocate(ilocN_tmp_buf1)
         allocate(ilocN_tmp_buf1(ibuf))
         do nucl=1,nneutrns
            iallredbuf=nsdnf/ibuf
            iiy=1
            do iallred=1,iallredbuf
               if (iproc==0) then
                  ilocN_tmp_buf1(:)=ilocn(nucl,iiy:iiy+ibuf-1)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,ibuf,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocn(nucl,iiy:iiy+ibuf-1)=ilocN_tmp_buf1(:)
               endif
               iiy=iiy+ibuf
            end do
            iallred=mod(nsdnf,ibuf)
            if (iallred/=0) then
               if (iallred/=nsdnf-iiy+1) then
                  print *,'***error:iallred,nsdp,iiy:',iallred,nsdnf,iiy
               endif
               if (iproc==0) then
                  ilocN_tmp_buf1(1:iallred)=ilocn(nucl,iiy:nsdnf)
               endif
               call MPI_Bcast(ilocN_tmp_buf1,iallred,MPI_INTEGER2,0,
     $              icomm,ierr)
               if (iproc/=0) then
                  ilocn(nucl,iiy:nsdnf)=ilocN_tmp_buf1(1:iallred)
               endif
            endif
         end do
         deallocate(ilocN_tmp_buf1)
      endif

      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsdf,3000000)
      if (nsdf<=ibuf) then
         call MPI_Bcast(delta_IT(1),nsdf,MPI_REAL,0,
     $        icomm,ierr)
      else
         if (allocated(delta_IT_tmp_buf)) deallocate(delta_IT_tmp_buf)
         allocate(delta_IT_tmp_buf(ibuf))
         iallredbuf=nsdf/ibuf
         iiy=1
         do iallred=1,iallredbuf
            if (iproc==0) then
               delta_IT_tmp_buf(:)=delta_IT(iiy:iiy+ibuf-1)
            endif
            call MPI_Bcast(delta_IT_tmp_buf,ibuf,MPI_REAL,0,
     $              icomm,ierr)
            if (iproc/=0) then
               delta_IT(iiy:iiy+ibuf-1)=delta_IT_tmp_buf(:)
            endif
            iiy=iiy+ibuf
         end do
         iallred=mod(nsdf,ibuf)
         if (iallred/=0) then
            if (iallred/=nsdf-iiy+1) then
               print *,'***error:iallred,nsd,iiy:',iallred,nsdf,iiy
            endif
            if (iproc==0) then
               delta_IT_tmp_buf(1:iallred)=delta_IT(iiy:nsdf)
            endif
            call MPI_Bcast(delta_IT_tmp_buf,iallred,MPI_REAL,0,
     $           icomm,ierr)
            if (iproc/=0) then
               delta_IT(iiy:nsdf)=delta_IT_tmp_buf(1:iallred)
            endif
         endif
         deallocate(delta_IT_tmp_buf)
      endif

      end

      subroutine basis_N_HO_reorder(nsd_orig)
      use nodeinfo
      use iters, only: nsd
      use spb
      use ibas
      use spsdata
      use ctrl,only: nf
      use ITvar,only: piv_saved
      implicit none
      integer,intent(IN) :: nsd_orig
      integer :: i1,ia,i2,N,i_count,NHO_i1,ip,in,k
      integer,allocatable :: iloc_tmp(:,:),imap(:)
      real(kind(0.d0)),allocatable :: bmp_orig(:),bmp_new(:)
      integer :: exminho,Nmin_HO,ierr,nsave
      character(len=80) :: line
      real(kind(0.0)),allocatable :: delta_IT_tmp(:) 

      nminho=Nmin_HO(nprotons)+Nmin_HO(nneutrns)
      exminho=mod(nhw-nminho,2)
      if (allocated(nsd_hw)) deallocate(nsd_hw)
      allocate(nsd_hw(exminho:nhw-nminho))
      nsd_hw=0

      if (iproc==0) then
         print *, '--reordering iloc and pivot--'
      endif
      allocate(imap(nsd_orig))
      allocate(iloc_tmp(2,nsd))
      allocate(delta_IT_tmp(nsd))
      i_count=0
ccc      do N=nhw0,nhw,2 !!! error!!!
ccc      do N=nminho,nhw,2 ! error!
      do N=nminho+exminho,nhw,2
         do i1=1,nsd
            ip=iloc(1,i1)
            in=iloc(2,i1)
            NHO_i1=0
            do ia=1,nprotons
               k=ilocp(ia,ip)
               NHO_i1=NHO_i1+2*n_sp(k)+l_sp(k)
            end do
            do ia=1,nneutrns
               k=ilocn(ia,in)
               NHO_i1=NHO_i1+2*n_sp(k)+l_sp(k)
            end do
            if (N==NHO_i1) then
               i_count=i_count+1
               iloc_tmp(:,i_count)=iloc(:,i1)
               delta_IT_tmp(i_count)=delta_IT(i1)
               if (i1<=nsd_orig) then
c                  if (i1/=i_count) then
c                     if (iproc==0) print *,' N,i1,i_count=',N,i1,i_count
c                  endif
                  imap(i1)=i_count
               endif
            endif
         end do
         nsd_hw(N-nminho)=i_count
         if (iproc==0) then
            print *,' N,Nmax,nsd_hw=',N,N-nminho,nsd_hw(N-nminho)
         endif
      end do
      if (i_count/=nsd) then
         print *,'*** error in basis_N_HO_reorder: i_count,nsd=',
     $        i_count,nsd
      endif
      do i1=1,nsd
         iloc(:,i1)=iloc_tmp(:,i1)
         delta_IT(i1)=delta_IT_tmp(i1)
      end do
      deallocate(iloc_tmp)
      deallocate(delta_IT_tmp)
      if (iproc==0)
     $     print *, '--iloc reordered: nsd_orig,nsd,i_count=',
     $     nsd_orig,nsd,i_count
      if (piv_saved) then
!!! restriction on number of saved pivot vectors
         nsave=min(0,nf)
!!!
         allocate(bmp_orig(nsd_orig))
         allocate(bmp_new(nsd))
         do i2=0,nsave
            bmp_new=0.d0
            if (iproc==0) print *,
     $           ' calling read_pivot from basis_N_HO_reorder'
            call read_pivot(nsd_orig,bmp_orig(:),i2) 
            do i1=1,nsd_orig
               i_count=imap(i1)
               bmp_new(i_count)=bmp_orig(i1)
            end do
            if (iproc==0) then
               if (i2==0) then
                  open(30,file='piv.tmp',            
     +                 access='sequential',
     +                 form='unformatted',status='unknown')
                  write(30) nsd
                  write(30) (bmp_new(i1),i1=1,nsd)
                  close(30)
                  print *,' piv.tmp saved'
               else
                  write(line,'(i8)') i2-1
                  line=adjustl(line)
                  open(30,file='piv'//trim(line)//'.tmp',
     +                 access='sequential',form='unformatted',
     $                 status='unknown')
                  write(30) nsd
                  write(30) (bmp_new(i1),i1=1,nsd)
                  close(30)
                  print *,' piv'//trim(line)//'.tmp saved'
               endif
            endif
         end do
         deallocate(bmp_orig)
         deallocate(bmp_new)
      endif
      if (iproc==0)
     $     print *, '--pivot reordered: nsd_orig,nsd,i_count=',
     $     nsd_orig,nsd,i_count
      call MPI_Barrier(icomm,ierr)
      call make_hash_iloc(.true.)
      if (iproc==0) print *,' make_hash_iloc called'
      end subroutine basis_N_HO_reorder

      subroutine Lanczos_hash
      use parameters
      use ibas
      use nuconf
      use lanczvec
      use nodeinfo
      use multig
      use ctrl
      use iters
      use spb
      use bits
      use albe
c      use hamc
c      use i1i3
      use hmes
      use nnmaxmod
      use tbmes
      use v3b
      use spsdata
      use occ
      use hash_tables
      use ITvar
      use tbmepar, only: ispe
      use TUD_tbme, only: tbmeTUD,no2bv3n,V_3N_no0b
      use pn_tbme, only: tbmepn
      implicit none
      include 'mpif.h'

      integer :: i1,i3
      integer :: ii,jj,kk,ll

c---------------------------------------------------------------------------
c**      integer(4) itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)
      integer(8) :: itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)

      real(kind(0.d0)) :: xamp,xxi,ovlp,ovlpmx,alph,betsq,betapre,x00

      real(8) :: cputime,walltime,timeused
      real(8),allocatable:: bmp0(:)
c      real(8),allocatable:: bmp1(:)
      real(kind=kind(0.d0)) :: denom
      integer ncut,mionb,ionb,ierr
c**      integer(4),allocatable :: ibascomp(:)
      integer(8),allocatable :: ibascomp(:)
c**      integer(4) ibascomi
      integer(8) :: ibascomi
c**      integer popcntt   ! integer(8) comment out
      integer :: stat(MPI_STATUS_SIZE)
c      integer(8) :: nhmer,nhmesum

      real :: ham,ham3b_cJ,ham2b_TUD,ham2b_pn
      real(kind(0.d0)) :: xcoef,espe
      integer :: mysave,iix,iterest=0,iter,na,nb,nesd,iro,iphase
      integer,allocatable :: occi(:),occf(:),occim2(:),occim3(:)
      integer :: N_HO_i,wd,ibit,interm_energy,N_HO_f
      integer(8),allocatable :: intbasi(:),intbasf(:),intbasim3(:) !,
c     $     intbas(:)
      integer :: ia_1,an_st_1,cr_1,cr_st_1
      integer :: ia_2,an_st_2,cr_st_2,sp_st_an_2,sp_st_cr_2,
     $     cr_st_min_1,cr_st_max_1,cr_st_min_2,cr_st_max_2
      integer :: ia_3,an_st_3,cr_st_3,cr_st_min_3,cr_st_max_3
      type multiple_type
      integer :: dim
      integer,allocatable :: multiple(:,:)
      end type multiple_type
      type MJ_type
      type(multiple_type),allocatable :: MJ(:)
      end type MJ_type
      type(MJ_type),allocatable :: parity_MT(:,:)
      integer :: mshift(0:1),mj,pi,mt,num_of_mult,minz,nspsmin,
     $     N_1,N_2,N_3,mj_ind,mt_ind,N123,N123_max
      integer :: mult_nucl_dim,ij
      integer,allocatable :: mult_nucl(:,:)
      integer,allocatable :: multiple_dim(:,:,:),multiple_point(:,:,:),
     $     multiple_ist(:,:)
      integer :: iallred,iallredbuf,iiy,northg
      character(len=80) :: line,line1
      logical :: save47,op12
      integer :: prot_st,neut_st,exminho
      real(kind(0.d0)) :: dnrm2, ddot
      integer :: Nmin_HO

c     do loop optimization
      integer :: N12_max_hold
      integer :: mnop_save
      integer :: ierr2,i1count
      integer(8),allocatable :: temp(:)
      integer(8) :: i3savedpoi,totolnonzerome,savedme,i3count,memtmp

      ibuf=min(nsd,3000000)
      
c     ------------------------------------------------------------------------
c     ----- bmp() is used in TRANSITIONS to store the final eigenstates. 
c     ------------------------------------------------------------------------
c

      if (associated(ty%nucl)) deallocate(ty%nucl)
      if (associated(ty%protn)) deallocate(ty%protn)
      if (associated(ty%neutn)) deallocate(ty%neutn)
      if (associated(ty%dim_st)) deallocate(ty%dim_st)
      if (associated(ty%dim_st_p)) deallocate(ty%dim_st_p)
      if (associated(ty%dim_st_n)) deallocate(ty%dim_st_n)


      allocate(ty%nucl)
      allocate(ty%protn)
      allocate(ty%neutn)
      allocate(ty%dim_st)
      allocate(ty%dim_st_p)
      allocate(ty%dim_st_n)
      ty%nucl=nucleons
      ty%protn=nprotons
      ty%neutn=nneutrns
      ty%dim_st=nsd
      ty%dim_st_p=nsdp
      ty%dim_st_n=nsdn
      ty%occ=>iloc
      ty%occp=>ilocp
      ty%occn=>ilocn

c     Setup hold variables for do loop optimization
      N12_max_hold = N12_max
c     N12_max is now correct for current nhw space (below)
c     No further mods need to be made in Lanczos_Hash
c      N12_max = nhw-(nhw_max-N12_max) 

      nspsmin=min(Nmin_HO(nprotons-2)+Nmin_HO(nneutrns),
     $     Nmin_HO(nprotons-1)+Nmin_HO(nneutrns-1),
     $     Nmin_HO(nprotons)+Nmin_HO(nneutrns-2))
      if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
      N12_max=nhw-nspsmin

c     test code
      if (iproc==0) write(6,*) ' nhw,N12_max:',nhw,N12_max

      mnop_save=mnop
      if (mnop==-3) then
         mnop=3
      endif
      if (iproc==0) print *,' mnop,mnop_save=',mnop,mnop_save
      
c     MKGK edit - this is called in make_hash_iloc
c      call construct_hash_table

      if (allocated(intbasi)) deallocate(intbasi)
      if (allocated(occi)) deallocate(occi)

      allocate(intbasi(2*nwd))
      allocate(occi(nucleons))
c      do i1=1,nsd
c         prot_st=iloc(1,i1)
c         neut_st=iloc(2,i1)
c         occi(1:nprotons)=ilocp(:,prot_st)
c         occi(nprotons+1:nucleons)=ilocn(:,neut_st)
c         intbasi=0
c         do cr_1=1,nucleons
c            wd=(occi(cr_1)-1)/nbit+1
c            ibit=mod(occi(cr_1)-1,nbit)
cc     print *,' i,wd,ibit=',i,wd,ibit
c            intbasi(wd)=ibset(intbasi(wd),ibit)
c         end do
c         call get_state_index(nucleons,2*nwd,intbasi,ionb)
c         if (i1/=ionb) then
c            print *,'*** error:'
c            if (iproc==0) then
c               print *,' i=',i1,'    index=',ionb
c               print *, iloc(:,i1)
c               print *,' occi=',occi
c            endif
c            stop
c         endif
c      end do

      if (mnop==3) then
         minz=max(0,nneutrns-2)+max(0,nprotons-2)
         nspsmin=max(minz-3,0)
         if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
         nspsmin=min(Nmin_HO(nprotons-3)+Nmin_HO(nneutrns),
     $        Nmin_HO(nprotons-2)+Nmin_HO(nneutrns-1),
     $        Nmin_HO(nprotons-1)+Nmin_HO(nneutrns-2),
     $        Nmin_HO(nprotons)+Nmin_HO(nneutrns-3))
         if (iproc==0) print *,' nwh,nspsmin=',nhw,nspsmin
         N123_max=nhw-nspsmin
         if (iproc==0) print *,' N123_max=',N123_max
         mult_nucl_dim=0
         do ia_1=1,nucleons-2
            do ia_2=ia_1+1,nucleons-1
               do ia_3=ia_2+1,nucleons
                  mult_nucl_dim=mult_nucl_dim+1
               end do
            end do
         end do
         if (iproc==0) print *,' mult_nucl_dim=',mult_nucl_dim
         memoryallocva=memoryallocva+mult_nucl_dim*3*4
         allocate(mult_nucl(3,mult_nucl_dim))
         ij=0
         do ia_1=1,nucleons-2
            do ia_2=ia_1+1,nucleons-1
               do ia_3=ia_2+1,nucleons
                  ij=ij+1
                  mult_nucl(1,ij)=ia_1
                  mult_nucl(2,ij)=ia_2
                  mult_nucl(3,ij)=ia_3
               end do
            end do
         end do
         if (ij/=mult_nucl_dim) then
            print *,'*** error: ii,mult_nucl_dim=',
     $           ij,mult_nucl_dim
            stop
         endif

c         allocate(parity_MT(0:1,0:mnop))
         mshift=1
         do cr_st_1=1,nasps-2
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N123_max) exit
            do cr_st_2=cr_st_1+1,nasps-1
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
c               if (N_1+N_2>N12_max) cycle
               if (N_1+N_2>N12_max) exit
               if (N_1+N_2>N123_max) exit
               do cr_st_3=cr_st_2+1,nasps
                  N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                  if (N_2+N_3>N12_max) exit
                  if (N_1+N_3>N12_max) exit
                  if (N_1+N_2+N_3>nhw-nspsmin) exit
c                  if (N_2+N_3>N12_max) cycle
c                  if (N_1+N_3>N12_max) cycle
c                  if (N_1+N_2+N_3>nhw-nspsmin) cycle
                  mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)+m2_sp(cr_st_3)
                  if (mj>mshift(1)) mshift(1)=mj
               end do
            end do
         end do
         do cr_st_1=1,nasps-1
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N123_max) exit
            do cr_st_2=cr_st_1+1,nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N12_max) exit
               if (N_1+N_2>N123_max) exit
c               if (N_1+N_2>N12_max) cycle
               do cr_st_3=mxsps+1,mxsps+nasps
                  N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                  if (N_2+N_3>N12_max) exit
                  if (N_1+N_3>N12_max) exit
                  if (N_1+N_2+N_3>nhw-nspsmin) exit
c                  if (N_2+N_3>N12_max) cycle
c                  if (N_1+N_3>N12_max) cycle
c                  if (N_1+N_2+N_3>nhw-nspsmin) cycle
                  mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)+m2_sp(cr_st_3)
                  if (mj>mshift(0)) mshift(0)=mj
               end do
            end do
         end do
         if (iproc==0) print *,' mshift=',mshift

         if (iproc==0) print *,' maxval(mshift)=',maxval(mshift)
         allocate(multiple_dim(0:1,0:3,0:maxval(mshift)))
         multiple_dim=0
         allocate(multiple_point(0:1,0:3,0:maxval(mshift)))
         multiple_point=0

         ij=0
         do pi=0,1
            do mt=-3,3,2

               select case(mt)
               case(3)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_min_3=3
                  cr_st_max_1=nasps-2
                  cr_st_max_2=nasps-1
                  cr_st_max_3=nasps
               case(1)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
                  cr_st_min_3=mxsps+1
                  cr_st_max_3=mxsps+nasps
               case(-1)
                  cr_st_min_1=1
                  cr_st_min_2=mxsps+1
                  cr_st_min_3=mxsps+2
                  cr_st_max_1=nasps
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case(-3)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_min_3=mxsps+3
                  cr_st_max_1=mxsps+nasps-2
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
c               allocate(parity_MT(pi,(mt+mnop)/2)%MJ(
c     $              0:mshift(abs(mt)/2)))
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  multiple_point(pi,mt_ind,mj_ind)=ij
                  num_of_mult=0
                  do N123=0,N123_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
                        if (N_1+N_2>N12_max) exit
c                        if (N_1+N_2>N12_max) cycle
                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
     $                       cr_st_max_3
                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                           if (N_1+N_2+N_3>N123) exit
                           if (N_2+N_3>N12_max) exit
                           if (N_1+N_3>N12_max) exit
c                           if (N_1+N_2+N_3>N123_max) exit
                           if (N_1+N_2+N_3/=N123) cycle
c                           if (N_2+N_3>N12_max) cycle
c                           if (N_1+N_3>N12_max) cycle
c                           if (N_1+N_2+N_3>N123_max) cycle
                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                          +mt2_sp(cr_st_3)/=mt) cycle
                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                          +m2_sp(cr_st_3)/=mj) cycle
                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                          +l_sp(cr_st_3),2)/=pi) cycle
                           num_of_mult=num_of_mult+1
                        end do
                     end do
                  end do
                  end do
c                  if (iproc==0) print *,' pi,mt,mj,num_of_mult=',
c     $                 pi,mt,mj,num_of_mult

c                  parity_MT(pi,mt_ind)%MJ(mj_ind)
c     $                 %dim=num_of_mult

                  ij=ij+num_of_mult
c                  if (num_of_mult>ij) ij=num_of_mult
                  multiple_dim(pi,mt_ind,mj_ind)=num_of_mult


                  if (num_of_mult==0) cycle
c                  allocate(parity_MT(pi,mt_ind)
c     $                 %MJ(mj_ind)
c     $                 %multiple(3,num_of_mult))
c                  num_of_mult=0
c                  do N123=0,N123_max
c                  do cr_st_1=cr_st_min_1,cr_st_max_1
c                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
c                     if (N_1>N123) exit
c                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
c                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
c                        if (N_1+N_2>N123) exit
c                        if (N_1+N_2>N12_max) exit
cc                        if (N_1+N_2>N12_max) cycle
c                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
c     $                       cr_st_max_3
c                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
c                           if (N_1+N_2+N_3>N123) exit
c                           if (N_2+N_3>N12_max) exit
c                           if (N_1+N_3>N12_max) exit
c                           if (N_1+N_2+N_3>N123_max) exit
c                           if (N_1+N_2+N_3/=N123) cycle
cc                           if (N_2+N_3>N12_max) cycle
cc                           if (N_1+N_3>N12_max) cycle
cc                           if (N_1+N_2+N_3>N123_max) cycle
c                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
c     $                          +mt2_sp(cr_st_3)/=mt) cycle
c                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
c     $                          +m2_sp(cr_st_3)/=mj) cycle
c                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
c     $                          +l_sp(cr_st_3),2)/=pi) cycle
c                           num_of_mult=num_of_mult+1
cc                           if (iproc==0) print *,
cc     $                      ' cr_st_1,cr_st_2,cr_st_3,num_of_mult=',
cc     $                          cr_st_1,cr_st_2,cr_st_3,num_of_mult

cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(1,num_of_mult)
cc     $                          =cr_st_1
cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(2,num_of_mult)
cc     $                          =cr_st_2
cc                           parity_MT(pi,mt_ind)
cc     $                          %MJ(mj_ind)
cc     $                          %multiple(3,num_of_mult)
cc     $                          =cr_st_3
c                        end do
c                     end do
c                  end do
c                  end do
               end do
            end do
         end do
         
         if (iproc==0) print *,' total num_of_mult=',ij
         memoryallocva=memoryallocva+ij*3*4
         allocate(multiple_ist(3,ij))
         multiple_ist=0
         ij=0
         do pi=0,1
            do mt=-3,3,2

               select case(mt)
               case(3)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_min_3=3
                  cr_st_max_1=nasps-2
                  cr_st_max_2=nasps-1
                  cr_st_max_3=nasps
               case(1)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
                  cr_st_min_3=mxsps+1
                  cr_st_max_3=mxsps+nasps
               case(-1)
                  cr_st_min_1=1
                  cr_st_min_2=mxsps+1
                  cr_st_min_3=mxsps+2
                  cr_st_max_1=nasps
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case(-3)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_min_3=mxsps+3
                  cr_st_max_1=mxsps+nasps-2
                  cr_st_max_2=mxsps+nasps-1
                  cr_st_max_3=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
c                  num_of_mult=0
                  do N123=0,N123_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
                        if (N_1+N_2>N12_max) exit
c                        if (N_1+N_2>N12_max) cycle
                        do cr_st_3=max(cr_st_2+1,cr_st_min_3),
     $                       cr_st_max_3
                           N_3=2*n_sp(cr_st_3)+l_sp(cr_st_3)
                           if (N_1+N_2+N_3>N123) exit
                           if (N_2+N_3>N12_max) exit
                           if (N_1+N_3>N12_max) exit
c                           if (N_1+N_2+N_3>N123_max) exit
                           if (N_1+N_2+N_3/=N123) cycle
c                           if (N_2+N_3>N12_max) cycle
c                           if (N_1+N_3>N12_max) cycle
c                           if (N_1+N_2+N_3>N123_max) cycle
                           if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                          +mt2_sp(cr_st_3)/=mt) cycle
                           if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                          +m2_sp(cr_st_3)/=mj) cycle
                           if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                          +l_sp(cr_st_3),2)/=pi) cycle
c                           num_of_mult=num_of_mult+1
                           ij=ij+1
                           multiple_ist(1,ij)
     $                          =cr_st_1
                           multiple_ist(2,ij)
     $                          =cr_st_2
                           multiple_ist(3,ij)
     $                          =cr_st_3
                        end do
                     end do
                  end do
                  end do
               end do
            end do
         end do

      elseif (mnop==2) then

         mult_nucl_dim=0
         do ia_1=1,nucleons-1
            do ia_2=ia_1+1,nucleons
               mult_nucl_dim=mult_nucl_dim+1
            end do
         end do
         if (iproc==0) print *,' mult_nucl_dim=',mult_nucl_dim
         memoryallocva=memoryallocva+mult_nucl_dim*2*4
         allocate(mult_nucl(2,mult_nucl_dim))
         ij=0
         do ia_1=1,nucleons-1
            do ia_2=ia_1+1,nucleons
               ij=ij+1
               mult_nucl(1,ij)=ia_1
               mult_nucl(2,ij)=ia_2
            end do
         end do
         if (ij/=mult_nucl_dim) then
            print *,'*** error: ii,mult_nucl_dim=',
     $           ij,mult_nucl_dim
            stop
         endif

c         allocate(parity_MT(0:1,0:mnop))
         mshift=0
         do cr_st_1=1,nasps-1
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N12_max) exit
            do cr_st_2=cr_st_1+1,nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)
               if (mj>mshift(1)) mshift(1)=mj
            end do
         end do
         do cr_st_1=1,nasps
            N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
            if (N_1>N12_max) exit
            do cr_st_2=mxsps+1,mxsps+nasps
               N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
               if (N_1+N_2>N12_max) exit
c               if (N_1+N_2>N12_max) cycle
               mj=m2_sp(cr_st_1)+m2_sp(cr_st_2)
               if (mj>mshift(0)) mshift(0)=mj
            end do
         end do

         if (iproc==0) print *,' mshift=',mshift

         if (iproc==0) print *,' maxval(mshift)=',maxval(mshift)
         allocate(multiple_dim(0:1,0:2,0:maxval(mshift)))
         multiple_dim=0
         allocate(multiple_point(0:1,0:2,0:maxval(mshift)))
         multiple_point=0

         ij=0
         do pi=0,1
            do mt=-2,2,2

               select case(mt)
               case(2)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
               case(0)
                  cr_st_min_1=1
                  cr_st_max_1=nasps
                  cr_st_min_2=mxsps+1
                  cr_st_max_2=mxsps+nasps
               case(-2)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_max_1=mxsps+nasps-1
                  cr_st_max_2=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2

c               allocate(parity_MT(pi,(mt+mnop)/2)%MJ(
c     $              0:mshift(abs(mt)/2)))
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  multiple_point(pi,mt_ind,mj_ind)=ij
                  num_of_mult=0
                  do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
c                        if (N_1+N_2>N12_max) exit
                        if (N_1+N_2/=N123) cycle
c                        if (N_1+N_2>N12_max) cycle
                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                       /=mt) cycle
                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                       ,2)/=pi) cycle
                        num_of_mult=num_of_mult+1
                     end do
                  end do
                  end do
c                  parity_MT(pi,mt_ind)%MJ(mj_ind)
c     $                 %dim=num_of_mult

                  ij=ij+num_of_mult
c                  if (num_of_mult>ij) ij=num_of_mult
                  multiple_dim(pi,mt_ind,mj_ind)=num_of_mult

                  if (num_of_mult==0) cycle
c                  allocate(parity_MT(pi,mt_ind)
c     $                 %MJ(mj_ind)
c     $                 %multiple(2,num_of_mult))
c                  num_of_mult=0
c                  do N123=0,N12_max
c                  do cr_st_1=cr_st_min_1,cr_st_max_1
c                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
c                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
c                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
c                        if (N_1+N_2>N12_max) cycle
c                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
c     $                       /=mt) cycle
c                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
c     $                       /=mj) cycle
c                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
c     $                       ,2)/=pi) cycle
c                        if (N_1+N_2/=N123) cycle
c                        num_of_mult=num_of_mult+1
cc                        parity_MT(pi,mt_ind)
cc     $                       %MJ(mj_ind)
cc     $                       %multiple(1,num_of_mult)
cc     $                       =cr_st_1
cc                        parity_MT(pi,mt_ind)
cc     $                       %MJ(mj_ind)
cc     $                       %multiple(2,num_of_mult)
cc     $                       =cr_st_2
c                     end do
c                  end do
c                  end do
               end do
            end do
         end do
         
         if (iproc==0) print *,' total num_of_mult=',ij
         memoryallocva=memoryallocva+ij*2*4
         allocate(multiple_ist(2,ij))
         multiple_ist=0
         ij=0

         do pi=0,1
            do mt=-2,2,2
               select case(mt)
               case(2)
                  cr_st_min_1=1
                  cr_st_min_2=2
                  cr_st_max_1=nasps-1
                  cr_st_max_2=nasps
               case(0)
                  cr_st_min_1=1
                  cr_st_max_1=nasps
                  cr_st_min_2=mxsps+1
                  cr_st_max_2=mxsps+nasps
               case(-2)
                  cr_st_min_1=mxsps+1
                  cr_st_min_2=mxsps+2
                  cr_st_max_1=mxsps+nasps-1
                  cr_st_max_2=mxsps+nasps
               case default
                  cycle
               end select
               mt_ind=(mt+mnop)/2
               do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2
                  do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_sp(cr_st_1)+l_sp(cr_st_1)
                     if (N_1>N123) exit
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_sp(cr_st_2)+l_sp(cr_st_2)
                        if (N_1+N_2>N123) exit
c                        if (N_1+N_2>N12_max) exit
                        if (N_1+N_2/=N123) cycle
c                        if (N_1+N_2>N12_max) cycle
                        if (mt2_sp(cr_st_1)+mt2_sp(cr_st_2)
     $                       /=mt) cycle
                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_sp(cr_st_1)+l_sp(cr_st_2)
     $                       ,2)/=pi) cycle
                        ij=ij+1
                        multiple_ist(1,ij)
     $                       =cr_st_1
                        multiple_ist(2,ij)
     $                       =cr_st_2
                     end do
                  end do
                  end do
               end do
            end do
         end do

      endif

      xxi=0.d0     
c     ----- Initialize amp:

      if (allocated(amp)) deallocate(amp)
      allocate(amp(nsd,2))
      amp = 0.d0           ! New Lanczos vector
c      call pivot(amp(:,1))  ! Initial Lanczos vector

c     make the pivot call
      if (irest/=2.and.irest/=4) then
         if (piv_saved) then
            if (iproc==0) print *,
     $           ' calling read_pivot from Lanczos_hash'
            call read_pivot(nsd,amp(:,1),0)
         else 
      
            exminho=mod(nhw-nminho,2)
            if (nhw-nminho>=4) then
               x00=sqrt(0.65d0/real(nsd_hw(exminho),kind(0.d0)))
               amp(1:nsd_hw(exminho),1) = x00 
               x00=sqrt(0.25d0/real(nsd_hw(exminho+2)-nsd_hw(exminho),
     $              kind(0.d0)))
               do ij=nsd_hw(exminho)+1,nsd_hw(exminho+2)
                  amp(ij,1) = x00 
               end do
               x00=sqrt(0.1d0/real(nsd_hw(exminho+4)-nsd_hw(exminho+2),
     $              kind(0.d0)))
               do ij=nsd_hw(exminho+2)+1,nsd_hw(exminho+4)
                  amp(ij,1) = x00 
               end do
            elseif (nhw-nminho>=2) then
               x00=sqrt(0.65d0/real(nsd_hw(exminho),kind(0.d0)))
               amp(1:nsd_hw(exminho),1) = x00 
               x00=sqrt(0.35d0/real(nsd_hw(exminho+2)-nsd_hw(exminho),
     $              kind(0.d0)))
               do ij=nsd_hw(exminho)+1,nsd_hw(exminho+2)
                  amp(ij,1) = x00 
               end do
            else
               x00=sqrt(1.d0/real(nsd,kind(0.d0)))
               do ij=1,nsd
                  amp(ij,1) = x00*(-1)**ij
               end do
            endif
         end if                 ! piv_saved or not

c     check norm of new pivot
         x00=0.d0
c      print *,' x00=',x00
         do ij=1,nsd
            x00=x00+amp(ij,1)*amp(ij,1)
c         print *,' i1,amp(i1,1)=',i1,amp(i1,1),x00
         end do
         if (iproc==0) print *,' norm=',x00
      endif

      nnmax = min(nite,nsd)     ! Number of Lanczos iterations
      mysave=0

c     MKGK
      if (allocated(alpha)) deallocate(alpha)
      if (allocated(beta)) deallocate(beta)
      allocate(alpha(nnmax),beta(nnmax))
      alpha=0.d0
      beta=0.d0
         
      if(irest/=2.and.irest/=4) then
         if (iproc==0) then
            mysave=mysave+1
            write(12)(amp(i1,1),i1=1,nsd)
         endif
      endif   

      if (irest==2) return

      if (irest==4) then
         open(29,file='lanczme.tmp',status='old')
         iix=0
         do ij=1,nnmax
            read(29,*,end=4574,err=4574) iix,alpha(ij),beta(ij)
         end do
 4574    continue
c**         iterest=max(1,iix-nproc)
c**         iterest=max(1,iix-2)
cc         iterest=iix
         iterest=iix+1
         if (iproc==0) then
            print *,' iterest=',iterest
            write(9,*) ' Restarting Lanczos iterations from'
     +           ,' #',iterest
            write(8,*) ' Restarting Lanczos iterations from'
     +           ,' #',iterest
         endif   
         rewind(29)
         if (iterest>1) then
            if (iproc==0) then
               read(12) (amp(i1,1),i1=1,nsd)
               mysave=mysave+1
            endif
         endif   
      endif   
c-------------------------------------------------------------------------
c     ----- THE MAIN LOOP: LANCZOS ITERATIONS
c     ----- The matrix elements of H are EITHER read from the files 
c     ----- mfdp.tmp1 and mfdp.tmp2 OR calculated in the first iteration.


      if (allocated(intbasf)) deallocate(intbasf)
cc      if (allocated(intbas)) deallocate(intbas)
      if (allocated(intbasim3)) deallocate(intbasim3)
      if (allocated(occim2)) deallocate(occim2)
      if (allocated(occim3)) deallocate(occim3)
      if (allocated(occf)) deallocate(occf)
      if (allocated(energy)) deallocate(energy)

      allocate(intbasf(2*nwd))
cc      allocate(intbas(2*nwd))
      allocate(intbasim3(2*nwd))
      if (.not.allocated(locz)) allocate(locz(2*mxsps))
      allocate(occim2(nucleons-2))
      allocate(occim3(nucleons-3))
      allocate(occf(nucleons))

      allocate(energy(nnmax))
      energy=0.d0

c---------save in buffer allocations
!!!! iloc, amp(2),bmp      
      memoryallocva=memoryallocva+2*4*nsd+3*8*nsd+bufferdimi1saved+1
     $     +nprotons*2*ty%dim_st_p+nneutrns*2*ty%dim_st_n
      if (iproc==0) then
         print *,' memoryalloc  =',memoryalloc/1024/1024,' MB'
         print *,' memoryallocva=',memoryallocva/1024/1024,' MB'
         print *,' memavail     =',memavail,' MB'
      endif
      memtmp=memavail*1024
      memtmp=memtmp*1024
c      if (iproc==0)  print *,' memavail     =',memavail
!      bufferhsaved=(memtmp-memoryalloc-memoryallocva)/8/2
      bufferhsaved=(memtmp-memoryalloc-memoryallocva)/4/2
      if (iproc==0) then
         print *,' calculated bufferhsaved=',bufferhsaved
      endif
!      bufferhsaved=(memtmp*95/100-memoryalloc/95*100
!     $     -memoryallocva/95*100)/8/2
      bufferhsaved=(memtmp*90/100-memoryalloc/90*100
     $     -memoryallocva/90*100)/4/2
      if (bufferhsaved<0) bufferhsaved=0
c      bufferhsaved=bufferhsaved*7/10
      if (iproc==0) then
         print *,' using bufferhsaved=',bufferhsaved
      endif
      
      if (iproc==0) then
         print *,' saving H in buffer: bufferdimi1saved=',
     $        bufferdimi1saved
      endif
      if (allocated(dimi1saved)) deallocate(dimi1saved)
      allocate(dimi1saved(0:bufferdimi1saved),stat=ierr)
      if (ierr/=0) then
         bufferhsaved=0
         if (iproc==0) print *,' dimi1saved alloc failed'
      else
         dimi1saved(0)=0
      endif
      do
         if (bufferhsaved==0) exit
         if (allocated(i3saved)) deallocate(i3saved)
         if (allocated(hsaved)) deallocate(hsaved)
         allocate(i3saved(bufferhsaved),stat=ierr)
         allocate(hsaved(bufferhsaved),stat=ierr2)
         if (ierr/=0.or.ierr2/=0) then
            bufferhsaved=bufferhsaved/2
            if (iproc==0) then
               print *,' i3saved & hsaved alloc failed'
               print *,' reducing bufferhsaved to ',
     $              bufferhsaved
            endif
         else
            exit
         endif
      end do
      if (bufferhsaved>0.and.iproc==0) then
         print *,' hsaved & i3saved allocated: bufferhsaved=',
     $        bufferhsaved
      endif
      i3savedpoi=0
      
      do 6000 iter = 1, nnmax   ! nnmax = number of Lanczos iterations
         call clock(cputime,walltime,0,0)
         timeused = cputime
         na = mod(iter,2)
         if (na==0) na = 2
         nb = 1
         if (na==1) nb = 2

c         if (iproc==0) then
c            print *,' iter=',iter,'  na=',na,'  nb=',nb
c            print *,' amp(1,na)=',amp(1,na),'  amp(1,nb)=',amp(1,nb)
c            print *,' alpha=',alpha(iter),'  beta=',beta(iter)
c         endif

         if (irest==4.and.iter<iterest) then
            read(29,*) iix,alpha(iter),beta(iter)

c         print *,'iix=',iix,' alpha=',alpha(iter),'  beta=',beta(iter)

            if (iix/=iter) then
               if (iproc==0) print *,'*** error: iix,iter=',iix,iter
               call MPI_Abort(icomm,612,ierr)
               stop
            endif  
            if (iter<iterest-2) then
               if (mod(iter,nproc)==iproc) then
                  read(12) (amp(i1,1),i1=1,nsd)
                  mysave=mysave+1
               endif   
            elseif (iter==iterest-2) then
               if (mod(iter,nproc)==iproc) then
                  read(12) (amp(i1,nb),i1=1,nsd)
                  mysave=mysave+1
               endif      
               if (nproc>1) then
                  call MPI_Barrier(icomm,ierr)
                  if (nsd<=ibuf) then
                     call MPI_Bcast(amp(1,nb),nsd,
     +                    MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                  else
                     iallredbuf=nsd/ibuf
                     iiy=1
                     do iallred=1,iallredbuf
                        call MPI_Bcast(amp(iiy,nb),ibuf,
     +                       MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                        iiy=iiy+ibuf
                     end do
                     iallred=mod(nsd,ibuf)
                     if (iallred/=0) then
                        call MPI_Bcast(amp(iiy,nb),iallred,
     +                       MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                     endif
                  endif
               endif   
            elseif (iter==iterest-1) then
               if (mod(iter,nproc)==iproc) then
                  read(12) (amp(i1,nb),i1=1,nsd)
                  mysave=mysave+1
               endif

               betapre = beta(iter)
               if (nproc>1) then
                  call MPI_Barrier(icomm,ierr)
                  if (nsd<=ibuf) then
                     call MPI_Bcast(amp(1,nb),nsd,
     +                    MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                  else
                     iallredbuf=nsd/ibuf
                     iiy=1
                     do iallred=1,iallredbuf
                        call MPI_Bcast(amp(iiy,nb),ibuf,
     +                       MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                        iiy=iiy+ibuf
                     end do
                     iallred=mod(nsd,ibuf)
                     if (iallred/=0) then
                        call MPI_Bcast(amp(iiy,nb),iallred,
     +                       MPI_REAL8,mod(iter,nproc),icomm,ierr) 
                     endif
                  endif
               endif   
            endif
            cycle
         endif
     
         if (irest==6) iterest=1
         
c------------------------------------------------------------------------
         if (iter>1) then
            if(iproc==0)then
               amp(:,nb) = -betapre*amp(:,nb)
            else
               do i1 = 1, nsd
                  amp(i1,nb) = 0.d0 ! To all other nodes.
               enddo
            endif
         endif

c*** Standard calculation ***********************************

         do 5000 i1 = iproc+1, nsd, nproc
            prot_st=iloc(1,i1)
            neut_st=iloc(2,i1)
            occi(1:nprotons)=ilocp(:,prot_st)
            occi(nprotons+1:nucleons)=ilocn(:,neut_st)

            if (ispe==1) then
               espe=0.d0
               do ia_1=1,nprotons
                  an_st_1=occi(ia_1)
                  espe=espe+Esp(nobt_sp(an_st_1),1)
               end do
               do ia_1=nprotons+1,nucleons
                  an_st_1=occi(ia_1)
                  espe=espe+Esp(nobt_sp(an_st_1),-1)
               end do
               espe=espe+E_core
            elseif (no2bv3n) then
               espe=V_3N_no0b
            else
               espe=0.d0
            endif

            N_HO_i=0
            do ia_1=1,nucleons
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_sp(an_st_1)+l_sp(an_st_1)
            end do

            xamp = amp(i1,na)

            amp(i1,nb) = amp(i1,nb) + espe * xamp

            if (((irest==0.and.iter>1).or.(irest>0.and.iter>iterest))
     $           .and.i1<=i1maxsaved) then
               i1count=(i1-iproc-1)/nproc+1
               do i3count=dimi1saved(i1count-1)+1,dimi1saved(i1count)
                  i3=i3saved(i3count)
                  xcoef=hsaved(i3count)
                  prot_st=iloc(1,i3)
                  neut_st=iloc(2,i3)
                  occf(1:nprotons)=ilocp(:,prot_st)
                  occf(nprotons+1:nucleons)=ilocn(:,neut_st)
                  N_HO_f=0
                  do ia_1=1,nucleons
                     an_st_1=occf(ia_1)
                     N_HO_f=N_HO_f+2*n_sp(an_st_1)+l_sp(an_st_1)
                  end do
                  amp(i3,nb) = amp(i3,nb) + xcoef * xamp
                  if (N_HO_f<N_HO_i) amp(i1,nb) = amp(i1,nb)
     $                 +xcoef*amp(i3,na)
               end do
               cycle
            endif

            nesd = N_HO_i - nespmin
            
            intbasi=0
            do cr_1=1,nucleons
               wd=(occi(cr_1)-1)/nbit+1
               ibit=mod(occi(cr_1)-1,nbit)
c     print *,' i,wd,ibit=',i,wd,ibit
               intbasi(wd)=ibset(intbasi(wd),ibit)
            end do
            locz(1:occi(1))=0
            do cr_1=2,nucleons
               locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
            end do
            locz(occi(nucleons)+1:2*mxsps)=nucleons

            
            if (i1.eq.10000*(i1/10000).and.nproc==1) then 
               print *, ' doing i1=',i1 
            endif
c*pn* 
            if (mnop==3) then

c               print *,' i1',i1
c               print *,' occi=',occi

c               occf=occi ! test
c               locz=0 !test

*C$OMP  PARALLEL DO DEFAULT(SHARED)
*C$OMP& PRIVATE(ij,ia_1,ia_2,ia_3,an_st_1,an_st_2,an_st_3,intbasim3,wd,
*C$OMP& ibit,mt,mj,pi,interm_energy,mt_ind,mj_ind,cr_1,num_of_mult,
*C$OMP& cr_st_1,cr_st_2,cr_st_3,N_HO_f,intbasf,i3,iphase,xcoef)
*C$OMP& SCHEDULE(DYNAMIC)
               do ij=1,mult_nucl_dim

                  ia_1=mult_nucl(1,ij)
                  ia_2=mult_nucl(2,ij)
                  ia_3=mult_nucl(3,ij)

                  an_st_1=occi(ia_1)
                  an_st_2=occi(ia_2)
                  an_st_3=occi(ia_3)

                  intbasim3=intbasi
                  wd=(an_st_1-1)/nbit+1
                  ibit=mod(an_st_1-1,nbit)
                  intbasim3(wd)=ibclr(intbasim3(wd),ibit)
                  wd=(an_st_2-1)/nbit+1
                  ibit=mod(an_st_2-1,nbit)
                  intbasim3(wd)=ibclr(intbasim3(wd),ibit)
                  wd=(an_st_3-1)/nbit+1
                  ibit=mod(an_st_3-1,nbit)
                  intbasim3(wd)=ibclr(intbasim3(wd),ibit)
                  
                  mt=mt2_sp(an_st_1)
     $                 +mt2_sp(an_st_2)+mt2_sp(an_st_3)
                  mj=m2_sp(an_st_1)
     $                 +m2_sp(an_st_2)+m2_sp(an_st_3)
                  pi=mod(l_sp(an_st_1)+l_sp(an_st_2)
     $                 +l_sp(an_st_3),2)
                  
                  interm_energy=-2*n_sp(an_st_1)-l_sp(an_st_1)
     $                 -2*n_sp(an_st_2)-l_sp(an_st_2)
     $                 -2*n_sp(an_st_3)-l_sp(an_st_3)+N_HO_i


                  mt_ind=(mt+mnop)/2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2


                  cr_1=multiple_point(pi,mt_ind,mj_ind)

                  do num_of_mult=1,multiple_dim(pi,mt_ind,mj_ind)

                     cr_1=cr_1+1
                     cr_st_1=multiple_ist(1,cr_1)
                     cr_st_2=multiple_ist(2,cr_1)
                     cr_st_3=multiple_ist(3,cr_1)

                     N_HO_f=2*n_sp(cr_st_1)+l_sp(cr_st_1)
     $                    +2*n_sp(cr_st_2)+l_sp(cr_st_2)
     $                    +2*n_sp(cr_st_3)+l_sp(cr_st_3)
     +                    +interm_energy

                     if (N_HO_f>N_HO_i) exit
                           
                     intbasf=intbasim3
                     wd=(cr_st_1-1)/nbit+1
                     ibit=mod(cr_st_1-1,nbit)
                     if (btest(intbasf(wd),ibit)) cycle
                     intbasf(wd)=ibset(intbasf(wd),ibit)
                     wd=(cr_st_2-1)/nbit+1
                     ibit=mod(cr_st_2-1,nbit)
                     if (btest(intbasf(wd),ibit)) cycle
                     intbasf(wd)=ibset(intbasf(wd),ibit)
                     wd=(cr_st_3-1)/nbit+1
                     ibit=mod(cr_st_3-1,nbit)
                     if (btest(intbasf(wd),ibit)) cycle
                     intbasf(wd)=ibset(intbasf(wd),ibit)

c                     call get_state_index_loc(nucleons,2*nwd,intbasf,i3)
                     call get_state_index(nucleons,2*nwd,intbasf,i3)

                     if (i3==-1) cycle
                     iphase=ia_1+ia_2+ia_3
     $                    +locz(cr_st_1)+locz(cr_st_2)
     $                    +locz(cr_st_3)
                     if (cr_st_1>an_st_1) iphase=iphase+1
                     if (cr_st_1>an_st_2) iphase=iphase+1
                     if (cr_st_1>an_st_3) iphase=iphase+1
                     if (cr_st_2>an_st_1) iphase=iphase+1
                     if (cr_st_2>an_st_2) iphase=iphase+1
                     if (cr_st_2>an_st_3) iphase=iphase+1
                     if (cr_st_3>an_st_1) iphase=iphase+1
                     if (cr_st_3>an_st_2) iphase=iphase+1
                     if (cr_st_3>an_st_3) iphase=iphase+1

                     xcoef=real((-1)**iphase,kind(0.d0))
     $                    *ham3b_cJ(cr_st_3,cr_st_2,cr_st_1,
     $                    an_st_3,an_st_2,an_st_1)

c                                 print *,' xcoef'
c     MKGK - don't uncomment line below (xamp)
c     xamp = amp(i1,na)

*C$OMP CRITICAL
                     amp(i3,nb) = amp(i3,nb) + xcoef * xamp
                     if (N_HO_f<N_HO_i) amp(i1,nb) = amp(i1,nb)
     $                    +xcoef*amp(i3,na)
*C$OMP END CRITICAL
                     if (iter==1.or.iter==iterest) then
                        if (abs(xcoef)>1.d-14) then
                           i3savedpoi=i3savedpoi+1
                           if (i3savedpoi<=bufferhsaved) then
                              i3saved(i3savedpoi)=i3
                              hsaved(i3savedpoi)=xcoef
                           endif
                        endif
                     endif
                  end do
               end do
*C$OMP END PARALLEL DO

            else

*C$OMP  PARALLEL DO DEFAULT(SHARED)
*C$OMP& PRIVATE(ij,ia_1,ia_2,an_st_1,an_st_2,intbasim3,wd,
*C$OMP& ibit,mt,mj,pi,interm_energy,mt_ind,mj_ind,cr_1,num_of_mult,
*C$OMP& cr_st_1,cr_st_2,N_HO_f,intbasf,i3,iphase,xcoef,kk,ll,ii,jj)
*C$OMP& SCHEDULE(DYNAMIC)
               do ij=1,mult_nucl_dim

                  ia_1=mult_nucl(1,ij)
                  ia_2=mult_nucl(2,ij)

                  an_st_1=occi(ia_1)
                  an_st_2=occi(ia_2)

                  intbasim3=intbasi
                  wd=(an_st_1-1)/nbit+1
                  ibit=mod(an_st_1-1,nbit)
                  intbasim3(wd)=ibclr(intbasim3(wd),ibit)
                  wd=(an_st_2-1)/nbit+1
                  ibit=mod(an_st_2-1,nbit)
                  intbasim3(wd)=ibclr(intbasim3(wd),ibit)

                  mt=mt2_sp(an_st_1)
     $                 +mt2_sp(an_st_2)
                  mj=m2_sp(an_st_1)
     $                 +m2_sp(an_st_2)
                  pi=mod(l_sp(an_st_1)+l_sp(an_st_2),2)

                  interm_energy=-2*n_sp(an_st_1)-l_sp(an_st_1)
     $                 -2*n_sp(an_st_2)-l_sp(an_st_2)+N_HO_i
               
                  mt_ind=(mt+mnop)/2
                  mj_ind=(mshift(abs(mt)/2)+mj)/2

                  cr_1=multiple_point(pi,mt_ind,mj_ind)

                  do num_of_mult=1,multiple_dim(pi,mt_ind,mj_ind)
                     
                     cr_1=cr_1+1
                     
                     cr_st_1=multiple_ist(1,cr_1)
                     cr_st_2=multiple_ist(2,cr_1)

                     N_HO_f=2*n_sp(cr_st_1)+l_sp(cr_st_1)
     $                    +2*n_sp(cr_st_2)+l_sp(cr_st_2)
     +                    +interm_energy

c                     if (N_HO_f>nhw) exit
                     if (N_HO_f>N_HO_i) exit

                     intbasf=intbasim3
                     wd=(cr_st_1-1)/nbit+1
                     ibit=mod(cr_st_1-1,nbit)
                     if (btest(intbasf(wd),ibit)) cycle
                     intbasf(wd)=ibset(intbasf(wd),ibit)
                     wd=(cr_st_2-1)/nbit+1
                     ibit=mod(cr_st_2-1,nbit)
                     if (btest(intbasf(wd),ibit)) cycle
                     intbasf(wd)=ibset(intbasf(wd),ibit)

c                     call get_state_index_loc(nucleons,2*nwd,intbasf,
c     $                    i3)
                     call get_state_index(nucleons,2*nwd,intbasf,i3)
                     

                     if (i3==-1) cycle

                     iphase=ia_1+ia_2
     $                    +locz(cr_st_1)+locz(cr_st_2)+1
                     if (cr_st_1>an_st_1) iphase=iphase+1
                     if (cr_st_1>an_st_2) iphase=iphase+1
                     if (cr_st_2>an_st_1) iphase=iphase+1
                     if (cr_st_2>an_st_2) iphase=iphase+1

c                           iphase=(-1)**(ia_1+ia_2+locz(cr_st_1)
c     $                          +locz(cr_st_2)+1)
                     kk=an_st_1
                     ll=an_st_2
                     ii=cr_st_1
                     jj=cr_st_2
************************************************************
                     if (tbmepn) then 
                        xcoef=real((-1)**iphase,kind(0.d0))
     $                       *ham2b_pn(ii,jj,kk,ll)
                     elseif (tbmeTUD) then 
                        xcoef=real((-1)**iphase,kind(0.d0))
     $                       *ham2b_TUD(ii,jj,kk,ll)
                     else
                        xcoef=real((-1)**iphase,kind(0.d0))
     $                       *ham(ii,jj,kk,ll,nesd)
                     endif

************************************************************
*C$OMP CRITICAL
                     amp(i3,nb) = amp(i3,nb) + xcoef * xamp
                     if (N_HO_f<N_HO_i) amp(i1,nb) = amp(i1,nb)
     $                          +xcoef*amp(i3,na)
*C$OMP END CRITICAL
                     if (iter==1.or.iter==iterest) then
                        if (abs(xcoef)>1.d-14) then
                           i3savedpoi=i3savedpoi+1
                           if (i3savedpoi<=bufferhsaved) then
                              i3saved(i3savedpoi)=i3
                              hsaved(i3savedpoi)=xcoef
                           endif
                        endif
                     endif
                  end do
               end do
*C$OMP END PARALLEL DO
            endif
            if (((irest==0.and.iter==1).or.(irest>0.and.iter==iterest))
     $           .and.i3savedpoi<=bufferhsaved) then
               i1maxsaved=i1
               i1count=(i1-iproc-1)/nproc+1
               if (i1count>bufferdimi1saved) then
                  bufferdimi1saved=bufferdimi1saved+ibuf/3
                  if (allocated(temp)) deallocate(temp)
                  allocate(temp(0:bufferdimi1saved),stat=ierr)
                  if (ierr/=0) then
                     i1maxsaved=i1-1
                     i3savedpoi=bufferhsaved+1
                     print *,'*** iproc=',iproc,
     $          ': Next allocation of bufferdimi1saved failed at:',
     $                    bufferdimi1saved-ibuf/3
                     cycle
                  endif
                  temp(0:bufferdimi1saved-ibuf/3)=dimi1saved
                  call move_alloc(temp,dimi1saved)
                  if (iproc==0)
     $                 print *,' i1=',i1,
     $                 ' bufferdimi1saved increased to ',
     $                 bufferdimi1saved
               endif
               dimi1saved(i1count)=i3savedpoi
c               if (iproc==0) print *,' i1,i1count,dimi1saved=',
c     $              i1,i1count,dimi1saved(i1count)
            endif
 5000    continue

         if ((irest==0.and.iter==1).or.(irest>0.and.iter==iterest)) then
            if (iproc==0) then
               print *,' iproc=0: i1maxsaved=',i1maxsaved
               print *,' iproc=0: i3savedpoi=',i3savedpoi
            endif
            call MPI_Barrier(icomm,ierr)
            if (i3savedpoi>bufferhsaved) then
               savedme=bufferhsaved
            else
               savedme=i3savedpoi
            endif
            call MPI_Reduce(savedme,totolnonzerome,1,MPI_INTEGER8,
     +           MPI_SUM,0,icomm,ierr)
            if (iproc==0) then
               print *,' Total # of non-zero m.e. saved=',totolnonzerome
            endif
            call MPI_Barrier(icomm,ierr)
            call MPI_Reduce(i3savedpoi,totolnonzerome,1,MPI_INTEGER8,
     +           MPI_SUM,0,icomm,ierr)
            if (iproc==0) then
               print *,' Total # of non-zero m.e.      =',totolnonzerome
            endif
         endif
         
c         call MPI_Barrier(icomm,ierr)

         if (nproc>1) then
            call MPI_Barrier(icomm,ierr)
            if (iter==1.or.(irest==4.and.iter==iterest)) then
c**               if (allocated(bmp1)) deallocate(bmp1)
c**               allocate(bmp1(nsd))
               if (allocated(bmp0)) deallocate(bmp0)
               if (memsave==1) then
                  allocate(bmp0(ibuf))
                  bmp0=0.d0
               else
                  allocate(bmp0(nsd),stat=ierr)
                  if (ierr/=0) then
                     memsave=1
                     print *,' bmp0 allocation unsuccessful for iproc=',
     $                    iproc
                     allocate(bmp0(ibuf))
                  endif
                  bmp0=0.d0
               endif
            endif   
            if (nsd<=ibuf) then
               call MPI_Allreduce(amp(1,nb),bmp0(1),nsd,MPI_REAL8,
     +              MPI_SUM,icomm,ierr)
               amp(:,nb)=bmp0(:)
            else
               iallredbuf=nsd/ibuf
c               print *,' iproc,ibuf=',iproc,ibuf
c               print *,' iproc,iallredbuf=',iproc,iallredbuf
               iiy=1
               do iallred=1,iallredbuf
c                  print *,' iproc,iiy=',iproc,iiy
                  call MPI_Allreduce(amp(iiy,nb),bmp0(1),ibuf,MPI_REAL8,
     +              MPI_SUM,icomm,ierr)
                  amp(iiy:iiy+ibuf-1,nb)=bmp0(1:ibuf)
                  iiy=iiy+ibuf
               end do
               iallred=mod(nsd,ibuf)
               if (iallred/=0) then
                  if (iallred/=nsd-iiy+1) then
                     print *,'***error:iallred,nsd,iiy:',iallred,nsd,iiy
                  endif
c                  print *,' iproc,iallred,iiy=',iproc,iallred,iiy
                  call MPI_Allreduce(amp(iiy,nb),bmp0(1),iallred,
     $                 MPI_REAL8,MPI_SUM,icomm,ierr)
                  amp(iiy:nsd,nb)=bmp0(1:iallred)
c                  print *,' iproc,iiy,nsd,nsd-iiy+1,iallred=',iproc,
c     $                 iiy,nsd,nsd-iiy+1,iallred
               endif
            endif
c--------------------------------------------------------------------
c     ----- Compute overlap of phi(na) with phibar(na) [stored in amp(nb)]:
cc            alph=dot_product(amp(:,na),amp(:,nb))

            alph=ddot(nsd,amp(:,na),1,amp(:,nb),1)

c            alph=0.d0
*C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i1)
*C$OMP DO 
*C$OMP& SCHEDULE(DYNAMIC)
*C$OMP& REDUCTION(+:alph)
c            do i1=1,nsd
c               alph=alph+amp(i1,na)*amp(i1,nb)
c            end do
*C$OMP END DO NOWAIT
*C$OMP MASTER
            alpha(iter) = alph
*C$OMP END MASTER
cc            amp(:,nb)=amp(:,nb) - alph * amp(:,na)
*C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i1)
*C$OMP DO
*C$OMP& SCHEDULE(DYNAMIC)
            call daxpy(nsd,-alph,amp(:,na),1,amp(:,nb),1)
c            do i1=1,nsd
c               amp(i1,nb)=amp(i1,nb)-alph*amp(i1,na)
c            end do   
*C$OMP END DO NOWAIT
*C$OMP MASTER
cc            betsq = dot_product(amp(:,nb),amp(:,nb))
 
           betapre=dnrm2(nsd,amp(:,nb),1)

c            betsq=0.d0
*C$OMP END MASTER
*C$OMP DO 
*C$OMP& SCHEDULE(DYNAMIC)
*C$OMP& REDUCTION(+:betsq)
c            do i1=1,nsd
c               betsq=betsq+amp(i1,nb)*amp(i1,nb)
c            end do
*C$OMP END DO NOWAIT 
*C$OMP END MASTER
*C$OMP END PARALLEL
c            betapre = dsqrt(betsq)
            beta(iter) = betapre

            if(iter.eq.nsd) goto 5995 ! To jump to the end
            xamp = 1.d0/betapre

c            call dscal(nsd,xamp,amp(:,nb),1)
            amp(:,nb) = amp(:,nb)*xamp

            if (iproc==0) write(8,*)
     +           '  alph =',alph,'   betapre =',betapre

c            if (iproc==0) print *,
c     +           '  alph =',alph,'   betapre =',betapre

c------------------------------------------------------------------------     
c     ----- Second re-orthogonalization of amp( ,nb)
c     ----- Here the new amp is reorthogonalized to all previous amp's which
c     ----- are stored in different files if nproc > 1. What will be done is
c     ----- that the new amp will be reorthogonalized to the previous amp's.
c     ----- The total effect will be collected by calling global sum.
            
            if (memsave==1) then
               write(line,'(i8)') iproc
               line=adjustl(line)
               line1='amb'//trim(line)//'.tmp'
               line='ovl'//trim(line)//'.tmp'
               save47=.false.
               northg=0
            else
               bmp0=0.d0
            endif
            ovlpmx = 0.d0
            inquire(12,opened=op12)
            if (op12) rewind(12)
            do iro = 1,mysave

               if (iro==mysave.and.mod(iter-1,nproc)==iproc) exit

               if (memsave==1) then
                  open(47,file=trim(line1),access='sequential',
     $                 form='unformatted',status='unknown')
                  if (northg==0) then
                     write(47) (amp(i1,nb), i1=1,nsd)
                     close(47)
                     save47=.true.
                  else
                     read(47) (amp(i1,nb), i1=1,nsd)
                     close(47)
                  endif
                  if (mysave>1) then
                     open(48,file=trim(line),access='sequential',
     $                    form='unformatted',status='unknown')
                  endif
               endif

               read(12) (amp(i1,na), i1=1,nsd)
c               read(12) (bmp1(i1), i1=1,nsd)
cc               ovlp = dot_product(amp(:,nb),bmp1(:))

               ovlp=ddot(nsd,amp(:,nb),1,amp(:,na),1)
c               ovlp=0.d0
*C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i1) 
*C$OMP& SCHEDULE(DYNAMIC)
*C$OMP& REDUCTION(+:ovlp)
c               do i1=1,nsd
c                  ovlp=ovlp+amp(i1,nb)*amp(i1,na)
cc                  ovlp=ovlp+amp(i1,nb)*bmp1(i1)
c               end do
*C$OMP END PARALLEL DO 
               ovlpmx = max(ovlpmx,dabs(ovlp))
cc               bmp0(:) = bmp0(:) - ovlp*bmp1(:)
*C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i1) 
*C$OMP& SCHEDULE(DYNAMIC)

c               print *,' iproc,iro,ovlp,ovlpmx=',iproc,iro,ovlp,ovlpmx


               if (memsave==1) then
                  if (northg==0) then
                     amp(:,nb)=0.d0
                  else
                     read(48) (amp(i1,nb), i1=1,nsd)
                  endif

                  call daxpy(nsd,-ovlp,amp(:,na),1,amp(:,nb),1)
c                  do i1=1,nsd
c                     amp(i1,nb)=amp(i1,nb)-ovlp*amp(i1,na)
cc                  bmp0(i1)=bmp0(i1)-ovlp*amp(i1,na)
cc                  bmp0(i1)=bmp0(i1)-ovlp*bmp1(i1)
c                  end do
*C$OMP END PARALLEL DO 
                  northg=northg+1
                  if (iro<mysave) then
                     rewind(48)
                     write(48) (amp(i1,nb), i1=1,nsd)
                     close(48)
                  endif
               else
                  call daxpy(nsd,-ovlp,amp(:,na),1,bmp0(:),1)
c                  do i1=1,nsd
c                     bmp0(i1)=bmp0(i1)-ovlp*amp(i1,na)
c                  end do
               endif
            end do 
            
            call MPI_Barrier(icomm,ierr)
            if (memsave==1) then
               amp(:,na)=0.d0
               if (nsd<=ibuf) then
                  do ii=1,min(iter-2,nproc)
                     if (iproc==mod(ii,nproc)) then
                        bmp0(:)=amp(:,nb)
                     endif
                     call MPI_Bcast(bmp0(1),nsd,MPI_REAL8,
     +                    mod(ii,nproc),icomm,ierr)
                     do i1=1,nsd
                        amp(i1,na)=amp(i1,na)+bmp0(i1)
                     end do
                  end do 
               else
                  iallredbuf=nsd/ibuf
                  iiy=1
                  do iallred=1,iallredbuf
                     do ii=1,min(iter-2,nproc)
                        if (iproc==mod(ii,nproc)) then
                           bmp0(1:ibuf)=amp(iiy:iiy+ibuf-1,nb)
                        endif
                        call MPI_Bcast(bmp0(1),ibuf,MPI_REAL8,
     +                       mod(ii,nproc),icomm,ierr)
                        do i1=1,ibuf
                           amp(i1+iiy-1,na)=amp(i1+iiy-1,na)+bmp0(i1)
                        end do
                     end do
                     iiy=iiy+ibuf
                  end do
                  iallred=mod(nsd,ibuf)
                  if (iallred/=0) then
                     if (iallred/=nsd-iiy+1) then
                        print *,'***error:iallred,nsd,iiy:',
     $                       iallred,nsd,iiy
                     endif
                     do ii=1,min(iter-2,nproc)
                        if (iproc==mod(ii,nproc)) then
                           bmp0(1:iallred)=amp(iiy:iiy+iallred-1,nb)
                        endif
                        call MPI_Bcast(bmp0(1),iallred,MPI_REAL8,
     +                       mod(ii,nproc),icomm,ierr)
                        do i1=1,iallred
                           amp(i1+iiy-1,na)=amp(i1+iiy-1,na)+bmp0(i1)
                        end do
                     end do
                  endif
               endif
               if (save47) then
                  open(47,file=trim(line1),access='sequential',
     $                 form='unformatted',status='old')
                  read(47) (amp(i1,nb), i1=1,nsd)
                  close(47)
               endif
            else
               if (nsd<=ibuf) then
                  call MPI_Allreduce(bmp0(1),amp(1,na),nsd,MPI_REAL8,
     +                 MPI_SUM,icomm,ierr)
               else
                  iallredbuf=nsd/ibuf
                  iiy=1
                  do iallred=1,iallredbuf
                     call MPI_Allreduce(bmp0(iiy),amp(iiy,na),ibuf,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     iiy=iiy+ibuf
                  end do
                  iallred=mod(nsd,ibuf)
                  if (iallred/=0) then
                     if (iallred/=nsd-iiy+1) then
                        print *,'***error:iallred,nsd,iiy:',
     $                       iallred,nsd,iiy
                     endif
c                  print *,' iproc,iallred,iiy=',iproc,iallred,iiy
                     call MPI_Allreduce(bmp0(iiy),amp(iiy,na),iallred,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                  endif
               endif
            endif
            call daxpy(nsd,1.d0,amp(:,na),1,amp(:,nb),1)
c            do i1=1,nsd
c               amp(i1,nb)=amp(i1,nb)+amp(i1,na)
c            end do
c               amp(:,nb)=amp(:,nb)+amp(:,na)
c            amp(:,nb)=amp(:,nb)+bmp1(:)
            call MPI_Allreduce(ovlpmx,ovlp,1,MPI_REAL8,
     +           MPI_MAX,icomm,ierr)
            ovlpmx=ovlp

c            if (iproc==0) print *,' iter,ovlpmx=',iter,ovlpmx


c     ----- Renormalize amp( ,nb):
cc            alph = dot_product(amp(:,nb),amp(:,nb))
            alph=dnrm2(nsd,amp(:,nb),1)
            alph=alph**2

c            if (iproc==0) print *,' norm of amp(nb): iter,alph=',
c     $           iter,alph

c            alph=0.d0
*C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i1) 
*C$OMP& SCHEDULE(DYNAMIC)
*C$OMP& REDUCTION(+:alph)
c            do i1=1,nsd
c               alph=alph+amp(i1,nb)*amp(i1,nb)
c            end do
c            if (iproc==0) print *,
c     $           ' norm of amp(nb) - no blas: iter,alph=',
c     $           iter,alph
*C$OMP END PARALLEL DO

            if (mod(iter-1,nproc)==iproc) then
               read(12) (amp(i1,na), i1=1,nsd)
               if (iproc/=0) then
                  if (nsd<=ibuf) then
                     call MPI_Send(amp(1,na),nsd,MPI_REAL8,
     +                    0,iter,icomm,ierr)
                  else
                     iallredbuf=nsd/ibuf
                     iiy=1
                     do iallred=1,iallredbuf
                        call MPI_Send(amp(iiy,na),ibuf,MPI_REAL8,
     +                       0,iter,icomm,ierr)
                        iiy=iiy+ibuf
                     end do
                     iallred=mod(nsd,ibuf)
                     if (iallred/=0) then
                        call MPI_Send(amp(iiy,na),iallred,MPI_REAL8,
     +                       0,iter,icomm,ierr)
                     endif
                  endif
               endif
            endif

            if (iproc==0) then
               if (mod(iter-1,nproc)/=0) then
                  if (nsd<=ibuf) then
                     call MPI_Recv(amp(1,na),
     $                    nsd,MPI_REAL8,mod(iter-1,nproc),iter,icomm,
     $                    stat,ierr)
                  else
                     iallredbuf=nsd/ibuf
                     iiy=1
                     do iallred=1,iallredbuf
                     call MPI_Recv(amp(iiy,na),
     $                    ibuf,MPI_REAL8,mod(iter-1,nproc),iter,icomm,
     $                    stat,ierr)
                        iiy=iiy+ibuf
                     end do
                     iallred=mod(nsd,ibuf)
                     if (iallred/=0) then
                        call MPI_Recv(amp(iiy,na),
     $                       iallred,MPI_REAL8,mod(iter-1,nproc),iter,
     $                       icomm,stat,ierr)
                     endif
                  endif
               endif
            endif
 
c               else
c                  iallredbuf=nsd/ibuf
c                  iiy=1
c                  do iallred=1,iallredbuf
c                     if (iproc/=0) then
c                        call MPI_Send(amp(iiy,na),ibuf,MPI_REAL8,
c     +                       0,iter,icomm,ierr)
c                     endif
c                     if (iproc==0) then
c                        if (mod(iter-1,nproc)/=0) call MPI_Recv(
c     $                       amp(iiy,na),ibuf,MPI_REAL8,
c     $                       mod(iter-1,nproc),iter,icomm,stat,ierr)
c                     endif 
c                     iiy=iiy+ibuf
c                  end do
c                  iallred=mod(nsd,ibuf)
c                  if (iallred/=0) then
c                     if (iproc/=0) then
c                        call MPI_Send(amp(iiy,na),iallred,MPI_REAL8,
c     +                       0,iter,icomm,ierr)
c                     endif
c                     if (iproc==0) then
c                        if (mod(iter-1,nproc)/=0) call MPI_Recv(
c     $                       amp(iiy,na),iallred,MPI_REAL8,
c     $                       mod(iter-1,nproc),iter,icomm,stat,ierr)
c                     endif 
c                  endif
c               endif
c            endif

         else

c     ----- Compute overlap of phi(na) with phibar(na) [stored in amp(nb)]:
            alph=dot_product(amp(:,na),amp(:,nb))
            alpha(iter) = alph

c            if (iproc==0) print *,' iter,alph=',iter,alph

            amp(:,nb)=amp(:,nb) - alph * amp(:,na)
c            do i1=1,nsd
c               amp(i1,nb)=amp(i1,nb)-alph*amp(i1,na)
c            end do
   
            betsq = dot_product(amp(:,nb),amp(:,nb))
            betapre = dsqrt(betsq)
            beta(iter) = betapre

c            if (iproc==0) print *,' iter,betapre=',iter,betapre

            if(iter.eq.nsd) goto 5995 ! To jump to the end
            xamp = 1.d0/betapre
            amp(:,nb) = amp(:,nb)*xamp
            if (iproc==0) write(8,*)
     &           '  alph =',alph,'   betapre =',betapre
c     ----- Second re-orthogonalization
            ovlpmx = 0.d0
            inquire(12,opened=op12)
            if (op12) rewind(12)
            do iro = 1, iter
               read(12) (amp(i1,na), i1=1,nsd)
               ovlp = dot_product(amp(:,nb),amp(:,na))
               ovlpmx = max(ovlpmx,dabs(ovlp))
               amp(:,nb) = amp(:,nb) - ovlp*amp(:,na)
c               do i1 = 1, nsd
c                  amp(i1,nb) = amp(i1,nb) - ovlp*amp(i1,na)
c               end do
            end do

c     ----- Renormalize amp( ,nb):
            alph = dot_product(amp(:,nb),amp(:,nb))
         endif

         alph = dsqrt(alph)
         if (iproc==0) write(8,*) '  Norm =', alph

c         if (iproc==0) print *,' Norm=', alph

         alph = 1.d0/alph
         amp(:,nb) = amp(:,nb)*alph
c     ----- Write amp to disk (sequential) or store amp in a node:

         if (mod(iter,nproc)==iproc) then
            mysave=mysave+1
c*            print *,'#',iproc,' mysave=',mysave
            write(12) (amp(i1,nb), i1=1,nsd)

         endif   
c-------------------------------------------------------------------------
         if (iproc==0) write(8,*)

 5995    continue
         if (iproc==0) then
            if (iter==1.and.irest/=4) then
               open(29,file='lanczme.tmp',status='unknown')
            elseif (irest==4.and.iter==iterest) then
               continue
            else
             open(29,file='lanczme.tmp',status='old',position='append')
            endif
            write(29,*) iter,alpha(iter),beta(iter)
            close(29)
         endif
         call Diag(nnmax,iter)

         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
         if(iproc.eq.0) then
            if (iter>1) write(8,*) '  Maximal overlap =',ovlpmx
            write(8,*)
     &        '  CPU time for this Lanczos iteration (sec): ',timeused
         endif
 6000 continue
c      if (allocated(bmp1)) deallocate(bmp1)

c     MKGK
c     reset N12_max incase
      N12_max = N12_max_hold
      mnop=mnop_save
      if (allocated(amp)) deallocate(amp)
      if (allocated(bmp0)) deallocate(bmp0)
      if (allocated(temp)) deallocate(temp)
      if (allocated(i3saved)) deallocate(i3saved)
      if (allocated(hsaved)) deallocate(hsaved)
      if (allocated(dimi1saved)) deallocate(dimi1saved)
      irest=0
      iterest=1
      return
      end subroutine Lanczos_hash


c     from MKGK
      subroutine save_pivot
      use nodeinfo
      use vect
      use ctrl
      use iters
      use nnmaxmod !yes, must be here
      use spb
      use ITvar ! set piv_saved in here
      implicit none
      include 'mpif.h'
      real(kind(0.d0)),allocatable:: bmp0(:),bmp(:,:),amp(:),buffer(:)
      integer :: iro, i1, iter, i2
      integer :: k1max,k00,k1, ierr
      real(kind(0.d0)) :: norm
      integer :: iallred,iallredbuf,iiy
      character(len=80) :: line
      logical :: op12

c     set up some stuff for excited states
      k1max = min(iter_old,nf)
      if(k1max.le.0)return
      k1max = k1max - 2
      k00 = 1
      if(ki.ne.kf)then
         k1max = k1max + 1
         k00 = 0
      endif

      if (iproc==0) then
         write(6,*) 'k1max=',k1max
         write(6,*) 'k00=',k00
         write(6,*) 'iter_old=',iter_old
      end if

c     initializing
      if (allocated(bmp0)) deallocate(bmp0)
      allocate(bmp0(nsd))
      bmp0=0.d0
      if (allocated(amp)) deallocate(amp)
      allocate(amp(nsd))
      amp=0.d0
      if (k1max >= 0) then
         if (allocated(bmp)) deallocate(bmp)
         allocate(bmp(nsd,0:0))
         bmp=0.d0
      end if

c     create the eigenstate
      do k1=-1,k1max
         inquire(12,opened=op12)
         if (op12) rewind(12)
         if (k1==-1) then
            bmp0=0.d0
         else
            bmp=0.d0
         endif
         if (iproc==0) 
     +        write(8,*)
     $        '  Reading Lanczos vectors from unit 12-save_pivot'
         do iro = 1, iter_old
            if (mod(iro-1,nproc)==iproc) then
               read(12) (amp(i1),i1=1,nsd)
               if (k1==-1) then
                  do i1=1,nsd
                     bmp0(i1) = bmp0(i1) + v(iro,ki)*amp(i1)
                  end do
               else
c     do k1=0,k1max
                  do i1=1,nsd
                     bmp(i1,0)= bmp(i1,0)+v(iro,kf+k1+k00)*amp(i1)
                  end do
c     end do
               endif
            endif   
         enddo

         if (nproc>1) then
            call MPI_Barrier(icomm,ierr)
            if (nsd<=ibuf) then
               allocate(buffer(nsd))
               buffer=0.d0
               if (k1==-1) then
                  call MPI_Allreduce(bmp0(1),buffer(1),nsd,MPI_REAL8,
     +                 MPI_SUM,icomm,ierr)
                  bmp0(1:nsd)=buffer(1:nsd)
               else
c            do k1=0,k1max
                  call MPI_Allreduce(bmp(1,0),buffer(1),nsd,MPI_REAL8,
     +                 MPI_SUM,icomm,ierr)
                  bmp(1:nsd,0)=buffer(1:nsd)
c     end do
               endif
            else
               allocate(buffer(ibuf))
               buffer=0.d0
               iallredbuf=nsd/ibuf
               iiy=1
               do iallred=1,iallredbuf
                  if (k1==-1) then
                     call MPI_Allreduce(bmp0(iiy),buffer(1),ibuf,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp0(iiy:iiy+ibuf-1)=buffer(1:ibuf) !amp(iiy:iiy+ibuf-1)
c     do k1=0,k1max
                  else
                     call MPI_Allreduce(bmp(iiy,0),buffer(1),ibuf,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp(iiy:iiy+ibuf-1,0)=buffer(1:ibuf) !amp(iiy:iiy+ibuf-1)
                  endif
c               end do   
                  iiy=iiy+ibuf
               end do
               iallred=mod(nsd,ibuf)
               if (iallred/=0) then
                  if (iallred/=nsd-iiy+1) then
                     print *,'***error:iallred,nsd,iiy:',iallred,nsd,iiy
                  endif
                  if (k1==-1) then
                     call MPI_Allreduce(bmp0(iiy),buffer(1),iallred,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp0(iiy:nsd)=buffer(1:iallred) !amp(iiy:nsd)
c                     do k1=0,k1max
                  else
                     call MPI_Allreduce(bmp(iiy,0),buffer(1),iallred,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp(iiy:nsd,0)=buffer(1:iallred) !amp(iiy:nsd)
                  endif
c               end do   
               endif
            endif
            deallocate(buffer)
         endif
         if (k1==-1) then
c            amp(:)=bmp0(:)
c     save ground state for Itruncate
            if(iproc==0) then
               open(30,file='piv0.tmp',            
     +              access='sequential',
     +              form='unformatted',status='unknown')
               write(30) (nsd)
               write(30) (bmp0(i1),i1=1,nsd)
               close(30)
               norm=0.d0
               do i1=1,nsd
                  norm=norm+bmp0(i1)*bmp0(i1)
               end do
               if (iproc==0) write(6,*) 'bmp0 norm = ',norm,
     $              ' for k1=',k1
            end if

         else

c     ground state has index = 0 - not to be confused with i1=0
c     1st excited state index = 1 = 1+48 inside achar
c     this saves the excited states
            if (iproc==0) then
c               do i1=0,k1max
               write(line,'(i8)') k1+1
               line=adjustl(line)
               open(30,file='piv'//trim(line)//'.tmp',
     +              access='sequential',form='unformatted',
     $              status='unknown')
               write(30) (nsd)
               write(30) (bmp(i2,0),i2=1,nsd)
               close(30)
c     end do
               norm=0.d0
               do i1=1,nsd
                  norm=norm+bmp(i1,0)*bmp(i1,0)
               end do
               if (iproc==0) write(6,*) 'bmp norm = ',norm,' for k1=',k1
            end if
            if (k1max == (nf-2)) then
c     do i1=0,k1max
               if (iproc==0) print *,' amp updated for k1=',k1 
               bmp0(:)=bmp0(:)+bmp(:,0)
c               amp(:)=amp(:)+bmp(:,0)
c            end do
            end if
         endif
      end do
      if (iproc==0) write(8,*)
     +     '  Lanczos vectors have been read from unit 12-save pivot'

c     form the pivot state as a linear combination of
c     all 1+nf states present.
c     the nf-2 is required since k1max is reduced by 2 units above
c     THIS pivot state is used in Lanczos_Hash.
c      amp(:)=bmp0(:)
      if (k1max == (nf-2)) then
c         do i1=0,k1max
c            amp(:)=amp(:)+bmp(:,i1)
c         end do
         bmp0(:)=1.d0/dsqrt(1.d0*nf)*bmp0(:)
      end if

c     check norm
      norm=0.d0
      do i1=1,nsd
         norm=norm+bmp0(i1)*bmp0(i1)
      end do
      if (iproc==0) write(6,*) 'Save pivot amp norm = ',norm,'   nf=',nf
      if (abs(norm-1.d0)>1.d-5) then
         if (iproc==0) print *,
     $        ' something went wrong, renormalizing pivot'
         bmp0=bmp0/dsqrt(norm)
      endif

c save the entire pivot state
c the cmin cut is now done in Itruncate pivot call
      if(iproc==0) then
         open(30,file='piv.tmp',            
     +    access='sequential',
     +    form='unformatted',status='unknown')
         write(30) (nsd)
         write(30) (bmp0(i1),i1=1,nsd)
         close(30)
      end if

      deallocate(amp)
      deallocate(bmp0)

      if (allocated(bmp)) deallocate(bmp)
      if (iproc==0) write(6,*) 'save_pivot returned'
      piv_saved = .true.
      end subroutine save_pivot


      subroutine read_pivot(nsd,vec1,statex)
      use nodeinfo 
      implicit none
      include 'mpif.h'
      integer(4), intent(IN) :: nsd
      integer(4), intent(IN) :: statex
      real(kind=kind(0.d0)), intent(INOUT) :: vec1(nsd)
      integer :: nsd_old, i1
      real(kind(0.d0)) :: norm
      character(len=80) :: line
      integer :: ierr
      integer(4) :: ibuf,iallredbuf,iiy,iallred
      
c     pivot for the state you want (gs, 1st excited...)
c     The "-1" is there because the gs has index 1 in the initial call

c     combo pivot
      if (iproc==0) then
         print *,' read_pivot entered, state=',statex
         vec1(:)=0.d0
         if (statex==0) then
            open(30,file='piv.tmp',
     $           access='sequential',form='unformatted',status='old')
            
            read(30) nsd_old
            if (nsd_old>nsd) nsd_old=nsd
            read(30) (vec1(i1),i1=1,nsd_old)
            close(30)
            print *,' pivot read from piv.tmp' 
         else
            write(line,'(i8)') statex-1
            line=adjustl(line)
            open(30,file='piv'//trim(line)//'.tmp',
     $           access='sequential',form='unformatted',status='old')

            read(30) nsd_old
            if (nsd_old>nsd) nsd_old=nsd
            read(30) (vec1(i1),i1=1,nsd_old)
            close(30)
            print *,' pivot read from ','piv'//trim(line)//'.tmp' 
         end if

c     make sure the pivot is normalized
         norm=0.d0
         do i1=1,nsd
            norm=norm+vec1(i1)*vec1(i1)
         end do
         print *, ' read_pivot norm=',norm
         vec1(:) = vec1(:)/dsqrt(norm)
      endif
      call MPI_Barrier(icomm,ierr)
      ibuf=min(nsd,3000000)
      if (nsd<=ibuf) then
         call MPI_Bcast(vec1(1),nsd,MPI_REAL8,0,icomm,ierr) 
      else
         iallredbuf=nsd/ibuf
         iiy=1
         do iallred=1,iallredbuf
            call MPI_Bcast(vec1(iiy),ibuf,MPI_REAL8,0,icomm,ierr) 
            iiy=iiy+ibuf
         end do
         iallred=mod(nsd,ibuf)
         if (iallred/=0) then
            call MPI_Bcast(vec1(iiy),iallred,MPI_REAL8,0,icomm,ierr)
         endif
      endif
      call MPI_Barrier(icomm,ierr)
      end subroutine read_pivot


      subroutine pivot(vec1)
c*** Initial Lanczos vector ****
      use spb
      use nuconf
      use nodeinfo
      use iters
      implicit none
      real(kind=kind(0.d0)),intent(INOUT) :: vec1(nsd)
      real(kind=kind(0.d0)) :: x00,weight1,weight2,dnormgs,dnorm2hw
      data weight1/0.6d0/,weight2/0.4d0/
      integer :: iene1,iene2,icongs,icon2hw,ii,istartgs,istart2hw
      integer :: nstcongs,nstcon2hw

      if (major>0) then
         iene1=negy(1)
         icongs=1
         if (nconf>1) then
            if (negy(2)>negy(1)) then
               iene2=negy(2)
               icon2hw=2
            else
               iene2=negy(1)
               icon2hw=1
               iene1=negy(2)
               icongs=2
            endif
         else
            iene2=0
         endif
         do ii=2,nconf
            if (negy(ii)<iene1) then
               iene1=negy(ii)
               icongs=ii
            elseif (negy(ii)<iene2.and.negy(ii)>iene1) then
               iene2=negy(ii)
               icon2hw=ii
            endif
         end do
         if (iproc==0) then
            print *,' icongs=',icongs,' iene1=',iene1
            print *,' icon2hw=',icon2hw,' iene2=',iene2
         endif
      endif
      if (major>0.and.iene2/=0.and.iene2>iene1) then
         if (icongs==1) then
            nstcongs=iendconf(icongs)
            istartgs=1
         else
            nstcongs=iendconf(icongs)-iendconf(icongs-1)
            istartgs=iendconf(icongs-1)+1
         endif
         if (icon2hw==1) then
            nstcon2hw=iendconf(icon2hw)
            istart2hw=1
         else
            nstcon2hw=iendconf(icon2hw)-iendconf(icon2hw-1)
            istart2hw=iendconf(icon2hw-1)+1
         endif
         dnormgs=dsqrt(weight1/dble(nstcongs))
         dnorm2hw=dsqrt(weight2/dble(nstcon2hw))
         vec1=0.d0
         do ii=istartgs,istartgs+nstcongs-1
            vec1(ii)=dnormgs
         end do
         do ii=istart2hw,istart2hw+nstcon2hw-1
            vec1(ii)=dnorm2hw
         end do
         x00=0.d0
         do ii=1,nsd
            x00=x00+vec1(ii)**2
         end do
         if (iproc==0) print *,' pivot norm=',x00
      else
         x00=1.d0/dsqrt(dble(nsd))
         do ii=1,nsd
            vec1(ii) = x00*(-1)**ii 
         end do
      endif
      deallocate(negy)
      deallocate(iendconf)
      end


c     MKGK
c     make sure files are ready for use again to avoid 
c     'too much data error'
      subroutine file_open_close
      use nodeinfo
      use iters, only: nite
      implicit none
      character(len=80) :: line

c* PN
      if (iproc<=nite) then
         close(12)

         write(line,'(i8)') iproc
         line=adjustl(line)
         line='lancz'//trim(line)//'.tmp'
         open(12,file=trim(line),access='sequential',form='unformatted',
     $        status='unknown')
      endif
c**
         
      if (iproc==0) write(6,*) 'Files closed and opened'

      end subroutine file_open_close



c     ----- This function returns the two-body matrix element:
c     ----- < ii jj | H | kk ll > where indeces are sequence numbers of the
c     ----- single fermion states passed via the common block /hamc/. They 
c     ----- run from 1 to mxsps2.

      function ham(ii,jj,kk,ll,nesd)
      use parameters
      use nodeinfo
      use multig
      use jsgna
      use spodata
      use consts
      use spb
      use spsdata
      use tbmes
      use the3js
!      use TUD_tbme,only: tbmeTUD
c      use hamc
c      use i1i3
      include 'mpif.h'
      integer,intent(IN) :: ii,jj,kk,ll,nesd
c
!      if (tbmeTUD) then
!         ham=ham2b_TUD(ii,jj,kk,ll)
!         return
!      endif
      ham = 0.0
      m2a   = m2_sp(ii)
      m2b   = m2_sp(jj)
      m2c   = m2_sp(kk)
      m2d   = m2_sp(ll)
      mmjj = m2a + m2b
c
      itz2a = mt2_sp(ii)
      itz2b = mt2_sp(jj)
      itz2c = mt2_sp(kk)
      itz2d = mt2_sp(ll)
      mmtt = itz2a + itz2b
c     
c     ----- To be used in the3j(sps1,sps2,J)
c     ----- sps1 and sps2 run from 1 to nasps.
      ii0 = ii + (itz2a-1)/2*mxsps
      jj0 = jj + (itz2b-1)/2*mxsps
      kk0 = kk + (itz2c-1)/2*mxsps
      ll0 = ll + (itz2d-1)/2*mxsps
c     
c     ----- Selecting pp,nn,pn interaction 
      if (mmtt.eq.2) then
        ifaciso = ippint
        faciso=1.d0
      elseif (mmtt.eq.-2) then
        ifaciso = innint
        faciso=0.d0
      else
        ifaciso = 0
        faciso=0.d0
      endif   
c
      j2a = j2_sp(ii)
      j2b = j2_sp(jj)
      j2c = j2_sp(kk) 
      j2d = j2_sp(ll)
c
      if(nsetm1.eq.0)then
         nesp = 0
      elseif(nsetm1.eq.2) then
         nesp = ifaciso
      else
c     ----- Determining which set of G's to use based on the unperturbed 
c     ----- energy of the spectators:
      nesp = nesd - n2l_sp(kk) - n2l_sp(ll)
c***pn****
      inesp=(nesp+iset1)/2 
**********  
      nesp = min0(inesp, nsetm1)+ifaciso
c***pn****
*      if (inesp.gt.nsetm1) then
*         print *, ' error in m-G: (nesp+iset1)/2=',inesp,
*     + '    nsetm1=',nsetm1
*         write(8,*) ' error in m-G: (nesp+iset1)/2=',inesp,
*     + '    nsetm1=',nsetm1
*      endif
c************
      endif

      
 9    continue
      n1 = nobt_sp(ii)
      n2 = nobt_sp(jj)
      n3 = nobt_sp(kk)
      n4 = nobt_sp(ll)
c----------------------------------------------------
      jjamb = (j2a-j2b)/2
      jjcmd = (j2c-j2d)/2
      jmin = max(iabs(jjamb),iabs(jjcmd),iabs(mmjj)/2)
      phsj = real(jsgn(jjamb+jjcmd))
      jmax = min((j2a+j2b)/2,(j2c+j2d)/2)
      iprty = (1-(-1)**(l_sp(ii)+l_sp(jj)))/2
c
      do 100 jval = jmin, jmax, 1
         jjval = jval + jval
c&&&     if(iabs(mmjj).gt.jjval) goto 100
         t3j = the3j(ii0,jj0,jval)*the3j(kk0,ll0,jval)
         if(t3j.eq.0.0)goto 100
         do 50 it = iabs(mmtt)/2, 1
            if((n1.eq.n2.or.n3.eq.n4)
     &           .and.jsgn(jval+it).eq.1)goto 50
            iit = it + it
c&&&        if(iabs(mmtt).gt.iit) goto 50
            t3t = the3t(it,itz2a,itz2b)*the3t(it,itz2c,itz2d)
            if(t3t.eq.0.0)goto 50
c
            call findindx(n1,n2,n3,n4,jval,it,indx,phse)

 25         ham = ham + phsj*real((iit+1)*(jjval+1))*t3t
*     &           *t3j*gful(indx,nesp)*phse
     &           *t3j*(gful(indx,nesp)+faciso*cful(indx))*phse
 50      continue
 100  continue
      if(n1.eq.n2)ham = ham*sqrt2
      if(n3.eq.n4)ham = ham*sqrt2
      return
      end

      real function ham2b_TUD(ii_in,jj_in,kk_in,ll_in)
      use spsdata
      use paramdef
      use spb, only: nucleons
      use TUD_tbme
      use v3b,only: cgj12
      implicit none
      integer,intent(IN) :: ii_in,jj_in,kk_in,ll_in
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd,J12
      integer :: ia,ib,ic,id,iphase,temp,i_ab,i_cd,ii,jj,kk,ll,i_ad
      integer :: m2a,m2b,m2c,m2d
      integer :: itz2a,itz2b,itz2c,itz2d
      real(kind(0.d0)) :: Habcd,clbab,clbcd,
     $     clbabt(0:1,-1:1),clbcdt(0:1,-1:1),Hno1b

      ii=ii_in
      jj=jj_in
      kk=kk_in
      ll=ll_in
      
      iphase=1
      
      na=n_sp(ii)
      la=l_sp(ii)
      ja=j2_sp(ii)
      nb=n_sp(jj)
      lb=l_sp(jj)
      jb=j2_sp(jj)
      ia=nlj_st_TUD(na)%l(la)%j(ja/2)
      ib=nlj_st_TUD(nb)%l(lb)%j(jb/2)
      if (ib<=ia) then
         i_ab=index_ab(ia,ib)
      else
         temp=ii
         ii=jj
         jj=temp
         iphase=-iphase
         i_ab=index_ab(ib,ia)
      endif
      
      nc=n_sp(kk)
      lc=l_sp(kk)
      jc=j2_sp(kk)
      nd=n_sp(ll)
      ld=l_sp(ll)
      jd=j2_sp(ll)
      ic=nlj_st_TUD(nc)%l(lc)%j(jc/2)
      id=nlj_st_TUD(nd)%l(ld)%j(jd/2)
       if (id<=ic) then
         i_cd=index_ab(ic,id)
      else
         temp=kk
         kk=ll
         ll=temp
         iphase=-iphase
         i_cd=index_ab(id,ic)
      endif

      if (i_cd<=i_ab) then
         ja=j2_sp(ii)
         jb=j2_sp(jj)
         jc=j2_sp(kk)
         jd=j2_sp(ll)
      
         m2a = m2_sp(ii)
         m2b = m2_sp(jj)
         m2c = m2_sp(kk)
         m2d = m2_sp(ll)

         itz2a = mt2_sp(ii)
         itz2b = mt2_sp(jj)
         itz2c = mt2_sp(kk)
         itz2d = mt2_sp(ll)

         ii=index_abcd(i_ab,i_cd)
      else
         ja=j2_sp(kk)
         jb=j2_sp(ll)
         jc=j2_sp(ii)
         jd=j2_sp(jj)
      
         m2a = m2_sp(kk)
         m2b = m2_sp(ll)
         m2c = m2_sp(ii)
         m2d = m2_sp(jj)

         itz2a = mt2_sp(kk)
         itz2b = mt2_sp(ll)
         itz2c = mt2_sp(ii)
         itz2d = mt2_sp(jj)
         ii=index_abcd(i_cd,i_ab)
      endif

      Hno1b=0.d0
      if (no2bv3n) then
         if (ii_in==ll_in) then
            if (ib<=ic) then
               i_ad=ib+ic*(ic-1)/2
            else
               i_ad=ic+ib*(ib-1)/2
            endif
            Hno1b=Hno1b-V_3N_no1b_TUD(i_ad,(mt2_sp(jj_in)+1)/2)
         endif
          if (ii_in==kk_in) then
            if (ib<=id) then
               i_ad=ib+id*(id-1)/2
            else
               i_ad=id+ib*(ib-1)/2
            endif
            Hno1b=Hno1b+V_3N_no1b_TUD(i_ad,(mt2_sp(jj_in)+1)/2)
         endif
         if (jj_in==ll_in) then
            if (ia<=ic) then
               i_ad=ia+ic*(ic-1)/2
            else
               i_ad=ic+ia*(ia-1)/2
            endif
            Hno1b=Hno1b+V_3N_no1b_TUD(i_ad,(mt2_sp(ii_in)+1)/2)
         endif
          if (jj_in==kk_in) then
            if (ia<=id) then
               i_ad=ia+id*(id-1)/2
            else
               i_ad=id+ia*(ia-1)/2
            endif
            Hno1b=Hno1b-V_3N_no1b_TUD(i_ad,(mt2_sp(ii_in)+1)/2)
         endif
         Hno1b=Hno1b/real(nucleons-1,kind(0.d0))
      endif

      clbabt(:,:)=cgt12_tud(:,:,(itz2a+1)/2,(itz2b+1)/2)
      clbcdt(:,:)=cgt12_tud(:,:,(itz2c+1)/2,(itz2d+1)/2)
      
      Habcd = 0.d0      
      do J12=max(abs(ja-jb),abs(jc-jd))/2,min(ja+jb,jc+jd)/2
         clbab=cgj12(J12,ja/2,(ja+m2a)/2,jb/2,(jb+m2b)/2)
         clbcd=cgj12(J12,jc/2,(jc+m2c)/2,jd/2,(jd+m2d)/2)
         
         Habcd=Habcd
     $        +H_NN_TUD(ii)*clbab*clbcd*clbabt(0,0)*clbcdt(0,0)
     $        +H_NN_TUD(ii+1)*clbab*clbcd*clbabt(1,-1)*clbcdt(1,-1)
     $        +H_NN_TUD(ii+2)*clbab*clbcd*clbabt(1,0)*clbcdt(1,0)
     $        +H_NN_TUD(ii+3)*clbab*clbcd*clbabt(1,1)*clbcdt(1,1)         

         ii=ii+4

      end do
      ham2b_TUD=Habcd*real(iphase,kind(0.d0))+Hno1b
      end function ham2b_TUD

      real function ham2b_pn(ii,jj,kk,ll)
      use spodata
      use spsdata
      use pn_tbme
      use v3b,only: cgj12
      implicit none
      integer,intent(IN) :: ii,jj,kk,ll
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd,J12,phase_cd
      integer :: ia,ib,ic,id,phase_ab,i_ab,i_cd,index,ip
      integer :: m2a,m2b,m2c,m2d,m2ab
      integer :: itz2a,itz2b,itz2c,itz2d,tz
      real(kind(0.d0)) :: Habcd,clbab,clbcd !,fac_ab,fac_cd

!      na=n_sp(ii)
      la=l_sp(ii)
      ja=j2_sp(ii)
      m2a=m2_sp(ii)
      itz2a=mt2_sp(ii)
      ia=nobt_sp(ii)

!      nb=n_sp(jj)
      lb=l_sp(jj)
      jb=j2_sp(jj)
      m2b=m2_sp(jj)
      itz2b=mt2_sp(jj)
      ib=nobt_sp(jj)

      ip=mod(la+lb,2)
      m2ab=m2a+m2b
      tz=(itz2a+itz2b)/2
!!! No check on parity,T12z conservation
!      nc=n_sp(kk)
!      lc=l_sp(kk)
      jc=j2_sp(kk)
      m2c=m2_sp(kk)
      itz2c=mt2_sp(kk)
      ic=nobt_sp(kk)

!      nd=n_sp(ll)
!      ld=l_sp(ll)
      jd=j2_sp(ll)
      m2d=m2_sp(ll)
!      itz2d=mt2_sp(ll)
      id=nobt_sp(ll)

!*** factors included in set_pn_H
!      if (abs(tz)==1) then
!         if (ia==ib) then
!            fac_ab=sqrt(2.d0)
!         else
!            fac_ab=1.d0
!         endif
!         if (ic==id) then
!            fac_cd=sqrt(2.d0)
!         else
!            fac_cd=1.d0
!         endif
!      else
!         fac_cd=1.d0
!     endif

!      print *,' ii,jj,kk,ll=',ii,jj,kk,ll
!      print *,' ia,ib,ic,id=',ia,ib,ic,id
!      print *,' ip,tz=',ip,tz
!      print *,' ja,jb,jc,jd=',ja,jb,jc,jd

      Habcd=0.d0
      do J12=max(abs(ja-jb),abs(jc-jd),abs(m2ab))/2,min(ja+jb,jc+jd)/2

!         print *,' J12=',J12
         
         if (abs(tz)==1) then
            if ((ia==ib.or.ic==id).and.mod(J12+1,2)==0) cycle
            if (ib<ia) then
               i_ab=pntbdst(J12,ip)%Tz(tz)%sp1(ib)%sp2(ia)
               phase_ab=(-1)**(J12-(ja-jb)/2)
            else
               i_ab=pntbdst(J12,ip)%Tz(tz)%sp1(ia)%sp2(ib)
               phase_ab=1
            endif
            if (id<ic) then
               i_cd=pntbdst(J12,ip)%Tz(tz)%sp1(id)%sp2(ic)
               phase_cd=(-1)**(J12-(jc-jd)/2)
            else
               i_cd=pntbdst(J12,ip)%Tz(tz)%sp1(ic)%sp2(id)
               phase_cd=1
            endif
         else
            if (itz2a==1) then
               i_ab=pntbdst(J12,ip)%Tz(0)%sp1(ia)%sp2(ib)
               phase_ab=1
            else
               i_ab=pntbdst(J12,ip)%Tz(0)%sp1(ib)%sp2(ia)
               phase_ab=(-1)**(J12-(ja-jb)/2)
            endif
            if (itz2c==1) then
               i_cd=pntbdst(J12,ip)%Tz(0)%sp1(ic)%sp2(id)
               phase_cd=1
            else
               i_cd=pntbdst(J12,ip)%Tz(0)%sp1(id)%sp2(ic)
               phase_cd=(-1)**(J12-(jc-jd)/2)
            endif
         endif

!         print *,' i_ab,icd=',i_ab,i_cd
         
         if (i_ab>=i_cd) then
            index=i2belnpoi_pn(J12,ip,tz)+i_cd+i_ab*(i_ab-1)/2
         else
            index=i2belnpoi_pn(J12,ip,tz)+i_ab+i_cd*(i_cd-1)/2
         endif

!         print *,' index=',index
         
         clbab=cgj12(J12,ja/2,(ja+m2a)/2,jb/2,(jb+m2b)/2)
         clbcd=cgj12(J12,jc/2,(jc+m2c)/2,jd/2,(jd+m2d)/2)

!         print *,' clbab=',clbab
!         print *,' clbcd=',clbcd
         
         Habcd=Habcd+clbab*clbcd*H_NN_pn(tz)%el(index) !*fac_ab*fac_cd ! included in set_pn_H
     $        *real(phase_ab*phase_cd,kind(0.d0))

!         print *,' H_NN_pn(tz)%el(index)=',H_NN_pn(tz)%el(index)
!         print *,' Habcd=',Habcd
         
      end do

      ham2b_pn=real(Habcd,kind(0.0))
      
      end function ham2b_pn

      real function ham3b_cJ(ispcr_1_in,ispcr_2_in,ispcr_3_in,
     $     ispan_1_in,ispan_2_in,ispan_3_in)
      use parameters
      use spb
      use spsdata
      use spbasis, only: isp1ntot
c      use hamc
      use v3b
      use TUD_tbme, only: tbmeTUD
      use pn_tbme, only: tbmepn
      implicit none
      integer,intent(IN) :: ispcr_1_in,ispcr_2_in,
     +           ispcr_3_in,ispan_1_in,ispan_2_in,ispan_3_in
      integer :: ispcr_1,ispcr_2,
     +           ispcr_3,ispan_1,ispan_2,ispan_3
      integer :: iphase,itemp,mjtot2,mttot2,itemp_i !,ipar
      integer:: i_a,i_b,i_c,i_d,i_e,i_f,i_abc,i_def,i_abcdef,iii
c     $     ,order_abc,order_def
      real(kind(0.d0)) :: sum,cg1,cg2,cg3,cg4,cgt1(0:1),cgt2(0:1),
     $     cgt3(0:1,0:1),cgt4(0:1,0:1),clebd,sum3
      integer :: j2_a,j2_b,j2_c,j2_d,j2_e,j2_f,j_ab,j_de,t_ab,t_de,
     $     J_3,T_3,m2_a,m2_b,m2_c,m2_d,m2_e,m2_f,mt2_a,mt2_b,mt2_c,
     $     mt2_d,mt2_e,mt2_f
      integer :: isp1n_i,isp1n_j,isp1n_l,isp1n_m,isp1n_n,isp1n_k
      real(kind=kind(0.0)) :: ham,ham2b_TUD,ham2b_pn

      integer :: ii,jj,kk,ll

      ispcr_1=ispcr_1_in
      ispcr_2=ispcr_2_in
      ispcr_3=ispcr_3_in
      ispan_1=ispan_1_in
      ispan_2=ispan_2_in
      ispan_3=ispan_3_in

      mjtot2=m2_sp(ispan_1)+m2_sp(ispan_2)+m2_sp(ispan_3)
c      if (mjtot2/=m2_sp(ispcr_1)+m2_sp(ispcr_2)+m2_sp(ispcr_3)) then
c         ham3b=0.0
c         return
c      endif
      mttot2=mt2_sp(ispan_1)+mt2_sp(ispan_2)+mt2_sp(ispan_3)
c      if (mttot2/=mt2_sp(ispcr_1)+mt2_sp(ispcr_2)+mt2_sp(ispcr_3)) then
c         ham3b=0.0
c         return
c      endif
c      ipar=mod(l_sp(ispan_1)+l_sp(ispan_2)+l_sp(ispan_3),2)

      isp1n_i=ispan_1
      isp1n_j=ispan_2
      isp1n_l=ispan_3

      isp1n_n=ispcr_1
      isp1n_m=ispcr_2
      isp1n_k=ispcr_3

      sum=0.d0

      if (tbmepn) then
         if (isp1n_k==isp1n_l) then 
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_k==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_k==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_pn(ii,jj,kk,ll)
         endif
      elseif (tbmeTUD) then
         if (isp1n_k==isp1n_l) then 
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_k==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_k==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_n==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
         if (isp1n_m==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham2b_TUD(ii,jj,kk,ll)
         endif
      else
         if (isp1n_k==isp1n_l) then 
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_k==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_k==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_n
            jj=isp1n_m
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_n==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_n==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_n==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_m
            jj=isp1n_k
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_m==isp1n_l) then
            kk=isp1n_i
            ll=isp1n_j
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_m==isp1n_i) then
            kk=isp1n_j
            ll=isp1n_l
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
         if (isp1n_m==isp1n_j) then
            kk=isp1n_l
            ll=isp1n_i
            ii=isp1n_k
            jj=isp1n_n
            sum=sum+ham(ii,jj,kk,ll,0)
         endif
      endif
      sum=sum/real(nucleons-2,kind(0.d0))

c      ham3b_cJ=sum
c      return
c 2034 continue

      i_a=nlj_orb(ispcr_1)
      i_b=nlj_orb(ispcr_2)
      i_c=nlj_orb(ispcr_3)
      i_d=nlj_orb(ispan_1)
      i_e=nlj_orb(ispan_2)
      i_f=nlj_orb(ispan_3)

      iphase=1
      if (i_b>i_a) then
         itemp=ispcr_2
         itemp_i=i_b
         ispcr_2=ispcr_1
         i_b=i_a
         ispcr_1=itemp
         i_a=itemp_i
         iphase=-iphase
      endif
      if (i_c>i_b) then
         itemp=ispcr_3
         itemp_i=i_c
         ispcr_3=ispcr_2
         i_c=i_b
         ispcr_2=itemp
         i_b=itemp_i
         iphase=-iphase
         if (i_b>i_a) then
            itemp=ispcr_2
            itemp_i=i_b
            ispcr_2=ispcr_1
            i_b=i_a
            ispcr_1=itemp
            i_a=itemp_i
            iphase=-iphase
         endif
      endif
      if (i_e>i_d) then
         itemp=ispan_2
         itemp_i=i_e
         ispan_2=ispan_1
         i_e=i_d
         ispan_1=itemp
         i_d=itemp_i
         iphase=-iphase
      endif
      if (i_f>i_e) then
         itemp=ispan_3
         itemp_i=i_f
         ispan_3=ispan_2
         i_f=i_e
         ispan_2=itemp
         i_e=itemp_i
         iphase=-iphase
         if (i_e>i_d) then
            itemp=ispan_2
            itemp_i=i_e
            ispan_2=ispan_1
            i_e=i_d
            ispan_1=itemp
            i_d=itemp_i
            iphase=-iphase
         endif
      endif

      i_abc=index_abc(i_a,i_b,i_c)
      i_def=index_abc(i_d,i_e,i_f)

      if (i_def>i_abc) then
c         itemp_i=i_a
c         i_a=i_d
c         i_d=itemp_i
c         itemp_i=i_b
c         i_b=i_e
c         i_e=itemp_i
c         itemp_i=i_c
c         i_c=i_f
c         i_f=itemp_i
         j2_a=j2_sp(ispan_1)
         j2_b=j2_sp(ispan_2)
         j2_c=j2_sp(ispan_3)
         m2_a=m2_sp(ispan_1)
         m2_b=m2_sp(ispan_2)
         m2_c=m2_sp(ispan_3)
         mt2_a=mt2_sp(ispan_1)
         mt2_b=mt2_sp(ispan_2)
         mt2_c=mt2_sp(ispan_3)
         j2_d=j2_sp(ispcr_1)
         j2_e=j2_sp(ispcr_2)
         j2_f=j2_sp(ispcr_3)
         m2_d=m2_sp(ispcr_1)
         m2_e=m2_sp(ispcr_2)
         m2_f=m2_sp(ispcr_3)
         mt2_d=mt2_sp(ispcr_1)
         mt2_e=mt2_sp(ispcr_2)
         mt2_f=mt2_sp(ispcr_3)
         i_abcdef=index_abcdef(i_def,i_abc)
      else
         j2_d=j2_sp(ispan_1)
         j2_e=j2_sp(ispan_2)
         j2_f=j2_sp(ispan_3)
         m2_d=m2_sp(ispan_1)
         m2_e=m2_sp(ispan_2)
         m2_f=m2_sp(ispan_3)
         mt2_d=mt2_sp(ispan_1)
         mt2_e=mt2_sp(ispan_2)
         mt2_f=mt2_sp(ispan_3)
         j2_a=j2_sp(ispcr_1)
         j2_b=j2_sp(ispcr_2)
         j2_c=j2_sp(ispcr_3)
         m2_a=m2_sp(ispcr_1)
         m2_b=m2_sp(ispcr_2)
         m2_c=m2_sp(ispcr_3)
         mt2_a=mt2_sp(ispcr_1)
         mt2_b=mt2_sp(ispcr_2)
         mt2_c=mt2_sp(ispcr_3)
         i_abcdef=index_abcdef(i_abc,i_def)
      endif
c      print *,' i_abcdef=',i_abcdef

      cgt1(0:1)=cgt12(0:1,(mt2_a+1)/2,(mt2_b+1)/2)
      cgt2(0:1)=cgt12(0:1,(mt2_d+1)/2,(mt2_e+1)/2)

      cgt3(0:1,0:1)=cgt123(0:1,0:1,(mt2_a+mt2_b)/2,(mt2_c+1)/2)
      cgt4(0:1,0:1)=cgt123(0:1,0:1,(mt2_d+mt2_e)/2,(mt2_f+1)/2)

      sum3=0.d0
      iii=start_abcdef(i_abcdef)
      do j_ab=abs(j2_a-j2_b),(j2_a+j2_b),2
c         if (abs(m2_a+m2_b)>j_ab) then
c            cg1=0.d0
c         else
c            cg1=clebd(j2_a,m2_a,j2_b,m2_b,j_ab,m2_a+m2_b)
            cg1=cgj12(j_ab/2,j2_a/2,(j2_a+m2_a)/2,j2_b/2,(j2_b+m2_b)/2)
c         endif
         do j_de=abs(j2_d-j2_e),(j2_d+j2_e),2
c            if (abs(m2_d+m2_e)>j_de.or.cg1==0.d0) then
c               cg2=0.d0
c            else
c               cg2=clebd(j2_d,m2_d,j2_e,m2_e,j_de,m2_d+m2_e)
               cg2=cgj12(j_de/2,j2_d/2,(j2_d+m2_d)/2,
     $              j2_e/2,(j2_e+m2_e)/2)
c            endif
            do J_3=max(abs(j_ab-j2_c),
     $           abs(j_de-j2_f)),
     $           min(j_ab+j2_c,j_de+j2_f),2
c               if (abs(mjtot2)>J_3.or.cg1==0.d0.or.cg2==0.d0) then
c                  cg3=0.d0
c                  cg4=0.d0
c               else
c                  cg3=clebd(j_ab,m2_a+m2_b,j2_c,m2_c,J_3,mjtot2)
c                  cg4=clebd(j_de,m2_d+m2_e,j2_f,m2_f,J_3,mjtot2)
                  cg3=cgj123(J_3/2,j_ab/2,(m2_a+m2_b)/2,
     $                 j2_c/2,(j2_c+m2_c)/2)
                  cg4=cgj123(J_3/2,j_de/2,(m2_d+m2_e)/2,
     $                 j2_f/2,(j2_f+m2_f)/2)
c               endif

               sum3=sum3+cg1*cg2*cg3*cg4*(
     $              v3b_cJ(iii)*cgt1(0)*cgt2(0)*cgt3(0,0)*cgt4(0,0)
     $              +v3b_cJ(iii+1)*cgt1(0)*cgt2(1)*cgt3(0,0)*cgt4(1,0)
     $              +v3b_cJ(iii+2)*cgt1(1)*cgt2(0)*cgt3(1,0)*cgt4(0,0)
     $              +v3b_cJ(iii+3)*cgt1(1)*cgt2(1)*cgt3(1,0)*cgt4(1,0)
     $              +v3b_cJ(iii+4)*cgt1(1)*cgt2(1)*cgt3(1,1)*cgt4(1,1))

               iii=iii+5

            end do
         end do
      end do
      ham3b_cJ=sum+sum3*real(iphase,kind(0.d0))
      end


c     ----- This subroutine diagonalizes the Lanczos matrix.

      subroutine Diag(nnmax,iter)
      use nodeinfo
      use parameters
      use albe
      use vect
c     MKGK
      use ITvar
      use iters
c     iter = the Lanczos iteration number
c      real d(iter),e(iter),e2(iter),wd(iter),rv4(iter),rv5(iter)
c      real rv1(iter),rv2(iter),rv3(iter),rv6(iter)
c      integer indd(iter)

      real,allocatable :: d(:),e(:),e2(:),wd(:),rv4(:),rv5(:),
     $     rv1(:),rv2(:),rv3(:),rv6(:)
      integer,allocatable :: indd(:)
      integer :: ifail

c
c     Prepare matrices to diagonalize

      allocate(d(iter),e(iter),e2(iter),wd(iter),rv4(iter),rv5(iter),
     $     rv1(iter),rv2(iter),rv3(iter),rv6(iter))
      allocate(indd(iter))

      if (allocated(v)) deallocate(v)
      allocate(v(iter,iter))
      do j = 1, iter
         d(j) = alpha(j)
         if (j==iter) exit
         e(j+1) = beta(j)
         e2(j+1) = beta(j)**2
      enddo
c     Prepare identity matrix for the nag routine
      do j = 1, iter, 1
         v(j,j) = 1.0
         do k = j+1, iter, 1
            v(j,k) = 0.0
            v(k,j) = 0.0
         enddo
      enddo
      eps1 = -1.0
      call tridib(iter,eps1,d,e,e2,rbb,ubb,
     &     1,iter,wd,indd,ifail,rv4,rv5)
      if(ifail.ne.0.and.iproc.eq.0) 
     &     write(8,*)'  Error in tridib ', ifail
      call tinvit(iter,iter,d,e,e2,iter,wd,indd,
     &     v,ifail,rv1,rv2,rv3,rv4,rv6)
      if(ifail.ne.0.and.iproc.eq.0) 
     &     write(8,*)'  Error in tinvit ', ifail
      if(iproc.eq.0)write(6,*)'Iteration #',iter,(wd(j),j=1,min(iter,3))
      if(iter.le.10.or.iter.eq.10*(iter/10)
     &     .or.iter.gt.(nnmax-3))then
         if(iproc.eq.0)then
            write(8,599)iter
            write(8,600) (j,wd(j),j=1,iter)
 599        format(/,'Iteration # ',i5)
 600        format(i8,f20.7)
         endif
      endif
      energy(1:iter)=wd(1:iter)
      iter_old = iter

      deallocate(d,e,e2,wd,rv4,rv5,rv1,rv2,rv3,rv6)
      deallocate(indd)

      return
      end


c     clean up some memory which is reallocated when Transitions
c     is called. This is required to run the bootstrap loop.
      subroutine cleaner
      use obmes
      use r2HOme
      use gamalog
      implicit none
      if (allocated(xlp)) then
         deallocate(xlp)
         deallocate(xln)
         deallocate(xsp)
         deallocate(xsn)
         deallocate(xqp)
         deallocate(xqn)
      end if
      if (allocated(r2mea)) deallocate(r2mea)
      if (allocated(dsq)) then
         deallocate(dsq,gamal)
      end if
      end subroutine cleaner


c     ----- This subroutine calculates the expectations and transitions for
c     ------ the initial state ki and the "nf" final states starting from kf.
c     ----- It is called by MAIN after all the Lanczos iterations are completed
c     ----- Some major changes were made here when parallelizing the program.
c     ----- Other major parallelizing changes were made in the subroutine
c     ----- Lanczos.

      subroutine Transitions
      use parameters
      use config
      use ibas
      use nuconf
      use lanczvec
      use nodeinfo
c      use multig
      use ctrl
      use jsgna
      use iters
      use spb
      use spsdata
      use effop
      use bits
      use albe
c      use hamc
c      use hmes
      use nnmaxmod
      use vect
      use obmes
      use occ
      use hash_tables
      use ITvar, only: kappa
      include 'mpif.h'

      integer :: i1,i3
      integer :: ii,jj,kk,ll

      real(8),allocatable:: bmp(:,:),bmp0(:)
c     
      real(8) xfacc,xxi         
c**      integer(4) itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)
      integer(8) itemp1(nwd,2),itemp3(nwd,2),itempz(nwd,2)
      real(8) xj2sm(-1:nf-1),xt2sm(-1:nf-1),xrmsm(-1:nf-1),
     &      xcoefs(-1:nf-1),xr1sm(-1:nf-1),xr2sm(-1:nf-1),
     &     gtpsm(-1:nf-1),gtmsm(-1:nf-1),
     &     x1bqpsm(-1:nf-1),x1bqnsm(-1:nf-1),x1blpsm(-1:nf-1),
     &     x1blnsm(-1:nf-1),x1bspsm(-1:nf-1),x1bsnsm(-1:nf-1),
     &     t1bqpsm(-1:nf-1),t1bqnsm(-1:nf-1),t1blpsm(-1:nf-1),
     &     t1blnsm(-1:nf-1),t1bspsm(-1:nf-1),t1bsnsm(-1:nf-1),
     &     ztmp(-1:nf-1)
      real,allocatable :: prob(:,:)
      real :: exeng(-1:nf-1),xjj(-1:nf-1),xtt(-1:nf-1) !,
c     &     occp(nshll),occn(nshll)
      real :: ham
      real(kind(0.d0)) :: xcoef
      integer :: mysave,iix,iterest,iter,na,nb,nesd,iro,iphase
      integer,allocatable :: occi(:),occf(:),occim1(:),occim2(:)
      integer :: N_HO_i,wd,ibit,interm_energy
      integer(8),allocatable :: intbas(:),intbasf(:)
      integer :: ia_1,an_st_1,cr_1,cr_st_1
      integer :: ia_2,an_st_2,cr_st_2,sp_st_an_2,sp_st_cr_2,
     $     cr_st_min_1,cr_st_max_1,cr_st_min_2,cr_st_max_2
      logical :: ob

      integer :: iallred,iallredbuf,iiy

      real(8) cputime,walltime,timeused
c      integer(8) :: nhmer
      integer(2) :: mcon=0
      integer(2),allocatable :: iloc_tmp(:)
      integer :: prot_st,neut_st
      logical :: kappadepwr,op12
c**      integer popcntt  ! for integer(8) comment out

c     ----- The eigenvectors were given in terms of the Lanczos vectors. The
c     ----- coefficients are stored in the array v(nite,nite), which is passed
c     ----- from Diag to this subroutine via a module. In this
c     ----- subroutine, the initial state and "nf" final state eigenvectors
c     ----- are to be re-expressed in terms of the original basis states 
c     ----- (i.e., Slater determinants). Each of these (1+nf) vectors is
c     ----- then specified by nsd coefficients. When all relevant vectors are
c     ----- transformed into the original basis, the array "bmp" will be 
c     ----- used to store the coefficients with different nodes' bmp's having
c     ----- the same content. (This is to avoid too much data communication
c     ----- in the computation of <Eigenvector_i | Operator | Eigenvector_j>,
c     ----- keeping in mind that the indexes (m,n) for the potentially 
c     ----- non-vanishing matrix elements <SD_m | Opeartor | SD_n> are stored
c     ----- in the arrays "irem" and "xrem" with different nodes having
c     ----- different contents.

      iter = nnmax
      if(irest/=2)then
         if (iproc==0) then
            open(14,file='lzegv.tmp',status='unknown')
c     ----- Write the lowest 200 eigenvectors for possible future use:
            write(6,*)'  writing to unit 14:'
            write(8,*)'  writing to unit 14:'
            do k1 = 1, min(iter,max(nf,min(200,nsd)))
               write(14,71)(v(iro,k1),iro=1,iter)
 71            format(5(1x,e16.9))
            enddo
            write(6,*)'  writing to unit 14 completed.'
            write(8,*)'  writing to unit 14 completed.'
            close(14)
         endif   
      else
c     ----- Read the lowest "nf" eigenvectors from the separate run:
         if (allocated(v)) deallocate(v)
         allocate(v(iter,iter))
         open(14,file='lzegv.tmp',status='unknown')
         if (iproc==0) then
            write(6,*)'  reading from unit 14:'
            write(8,*)'  reading from unit 14:'
         endif   
         do k1 = 1, min(iter,kf+nf-1)
            read(14,71)(v(iro,k1),iro=1,iter)
         enddo
         if (iproc==0) then
            write(6,*)'  reading from unit 14 completed.'
            write(8,*)'  reading from unit 14 completed.'
         endif   
         close (14)
      endif
      k1max = min(iter,nf)
      if(k1max.le.0)return
      k1max = k1max - 2
      k00 = 1
      if(ki.ne.kf)then
         k1max = k1max + 1
         k00 = 0
      endif

      call clock(cputime,walltime,0,0)
      timeused = cputime

      if (allocated(amp)) deallocate(amp)
      if (allocated(bmp0)) deallocate(bmp0)
      if (allocated(bmp)) deallocate(bmp)
      allocate(amp(nsd,1))      
      allocate(bmp0(nsd))
      allocate(bmp(nsd,0:0)) !k1max))    

      bmp0=0.d0
c***** start of the eigenstate loop
      do k1=-1,k1max
c*****
         bmp=0.d0
         inquire(12,opened=op12)
         if (op12) rewind(12)
         if (iproc==0) 
     +        write(8,*)'  Reading Lanczos vectors from unit 12:'
         do iro = 1, iter
            if (mod(iro-1,nproc)==iproc) then
               read(12) (amp(i1,1),i1=1,nsd)
c*               print *,'iro=',iro,' #',iproc
               if (k1==-1) then               
                  do i1=1,nsd
                     bmp0(i1) = bmp0(i1) + v(iro,ki)*amp(i1,1)
                  end do
               else
c            do k1 = 0, k1max
                  do i1=1,nsd
                     bmp(i1,0) = bmp(i1,0)+v(iro,kf+k1+k00)*amp(i1,1)
                  end do
               endif
c            enddo
            endif   
         enddo

         if (nproc>1) then
            call MPI_Barrier(icomm,ierr)
            if (nsd<=ibuf) then
               if (k1==-1) then
                  call MPI_Allreduce(bmp0(1),amp(1,1),nsd,MPI_REAL8,
     +                 MPI_SUM,icomm,ierr)
                  bmp0(:)=amp(:,1)
               else
c            do k1=0,k1max
                  call MPI_Allreduce(bmp(1,0),amp(1,1),nsd,MPI_REAL8,
     +                 MPI_SUM,icomm,ierr)
                  bmp(:,0)=amp(:,1)
c            end do   
               endif
            else
               iallredbuf=nsd/ibuf
               iiy=1
               do iallred=1,iallredbuf
                  if (k1==-1) then            
                     call MPI_Allreduce(bmp0(iiy),amp(iiy,1),ibuf,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp0(iiy:iiy+ibuf-1)=amp(iiy:iiy+ibuf-1,1)
                  else
c                  do k1=0,k1max
                     call MPI_Allreduce(bmp(iiy,0),amp(iiy,1),ibuf,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp(iiy:iiy+ibuf-1,0)=amp(iiy:iiy+ibuf-1,1)
c               end do   
                  endif
                  iiy=iiy+ibuf
               end do
               iallred=mod(nsd,ibuf)
               if (iallred/=0) then
                  if (iallred/=nsd-iiy+1) then
                     print *,'***error:iallred,nsd,iiy:',iallred,nsd,iiy
                  endif
                  if (k1==-1) then            
                     call MPI_Allreduce(bmp0(iiy),amp(iiy,1),iallred,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp0(iiy:nsd)=amp(iiy:nsd,1)
                  else
c               do k1=0,k1max
                     call MPI_Allreduce(bmp(iiy,0),amp(iiy,1),iallred,
     $                    MPI_REAL8,MPI_SUM,icomm,ierr)
                     bmp(iiy:nsd,0)=amp(iiy:nsd,1)
c               end do   
                  endif
               endif
            endif
         endif
         if (iproc==0) write(8,*)
     +        '  Lanczos vectors have been read from unit 12 for k1=',k1
         if (k1==k1max) then
            close(12)
            deallocate(amp)
         endif
c***pn*******************
         if (iproc==0) then
            if (k1==-1) then
               nhmex=0
               write(13) nsd,nhmex,k1max,nucleons,nprotons,nneutrns,
     $              hbomeg,nhw
               write(13) mjtotal,mttotal,nshll,nwd,mxsps,iparity,major
               if (memsave/=1) then
                  write(13) (mcon,i1=1,nsd)
               endif
               write(13) (bmp0(i1),i1=1,nsd)
            else
c      do k1=0,k1max
               write(13) (bmp(i1,0),i1=1,nsd)
c      end do
            endif
            if (k1==k1max) then
               write(13) (n_sp(i),l_sp(i),j2_sp(i),m2_sp(i),mt2_sp(i),
     +              i=1,mxsps2)
               if (memsave/=1) then
                  allocate(iloc_tmp(nsd))
                  do ia=1,nucleons
                     do i1=1,nsd
                        if (ia<=nprotons) then
                           prot_st=iloc(1,i1)
                           iloc_tmp(i1)=ilocp(ia,prot_st)
                        else
                           neut_st=iloc(2,i1)
                           iloc_tmp(i1)=ilocn(ia-nprotons,neut_st)
                        endif
                     end do 
                     write(13) (iloc_tmp(i1),i1=1,nsd)
                  end do
                  deallocate(iloc_tmp)
               else
                  write(13) ((iloc(ia,i1),ia=1,2),i1=1,nsd)
                  write(13) nsdp,nsdn
                  write(13) ((ilocp(ia,i1),ia=1,nprotons),i1=1,nsdp)
                  write(13) ((ilocn(ia,i1),ia=1,nneutrns),i1=1,nsdn)
               endif
            endif

            print *,' eigenvectors saved in unit 13 for k1=',k1   
            write(8,*) ' eigenvectors saved in unit 13 for k1=',k1   
         endif
c*************************

         if (k1==-1) then
c     ----- Call obme to calculate one-body matrix elements of observables:
            call obme
      
c     b^2 = (hbar c)^2/(Mc^2 hbar*Omega):
            bsquare = 197.327**2/(938.9185*hbomeg)
         endif

c      allocate(prob(nconf,-1:nf-1))
c     ---------------------------------------------------------------
c     ----- Use bmp to compute occupation probabilities:
cc      do k1 = -1, k1max,1       ! k1=-1 is for the initial state.
cc         do iconf = 1, nconf
cc            prob(iconf,k1) = 0.0
cc         enddo
cc         do i1 = 1, nsd
cc            iconf = mconf(i1)
cc            if(k1.eq.-1) then
cc               xyz = bmp0(i1)**2
cc            else
cc               xyz = bmp(i1,k1)**2
cc            endif
cc            prob(iconf,k1) = prob(iconf,k1) + xyz
cc         enddo
cc         do iconf = 1, nconf
cc            prob(iconf,k1) = 100.0*prob(iconf,k1)
cc         enddo
cc      enddo

cc      deallocate(mconf)

c     ----- Use bmp0 and bmp( ,k1) to compute expectation values of desired 
c     ----- operators. Treat as all two-body operators and follow the logic 
c     ----- in Lanczos.
c      do k1 = -1, k1max, 1
         xrmsm(k1) = 0.d0
         xr1sm(k1) = 0.d0
         xr2sm(k1) = 0.d0
         xj2sm(k1) = 0.d0
         xt2sm(k1) = 0.d0
         gtpsm(k1) = 0.d0
         gtmsm(k1) = 0.d0
c         xcoefs(k1) = 0.d0
c        ovlp(k1) = 0.d0
         x1bqpsm(k1)=0.d0
         x1bqnsm(k1)=0.d0
         x1blpsm(k1)=0.d0
         x1blnsm(k1)=0.d0
         x1bspsm(k1)=0.d0
         x1bsnsm(k1)=0.d0
c--------------------------------------------------------------------
c     ----- Need the followling lines if subshell occupations are desired.
c        do nzh=0,4
c           do lzh=0,8
c              do j2zh=1,17,2
c                 occupsm(k1,nzh,lzh,j2zh)=0.d0
c                 occunsm(k1,nzh,lzh,j2zh)=0.d0
c              enddo
c           enddo
c        enddo
c--------------------------------------------------------------------
         t1bqpsm(k1)=0.d0
         t1bqnsm(k1)=0.d0
         t1blpsm(k1)=0.d0
         t1blnsm(k1)=0.d0
         t1bspsm(k1)=0.d0
         t1bsnsm(k1)=0.d0
c      enddo


c     MKGK
         if (k1==-1) then
            if (allocated(intbas)) deallocate(intbas)
            if (allocated(intbasf)) deallocate(intbasf)
            if (allocated(occi)) deallocate(occi)
            if (allocated(occim2)) deallocate(occim2)
            if (allocated(occf)) deallocate(occf)

            allocate(intbas(2*nwd))
            allocate(intbasf(2*nwd))
ccc      allocate(locz(2*mxsps))
            allocate(occi(nucleons))
            allocate(occim2(nucleons-2))
            allocate(occf(nucleons))
         endif

         do 5000 i1 = iproc+1, nsd, nproc

            prot_st=iloc(1,i1)
            neut_st=iloc(2,i1)
            occi(1:nprotons)=ilocp(:,prot_st)
            occi(nprotons+1:nucleons)=ilocn(:,neut_st)
c         occi(:)=iloc(:,i1)
            N_HO_i=0
            do ia_1=1,nucleons
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_sp(an_st_1)+l_sp(an_st_1)
            end do
c         nesd = N_HO_i - nespmin

            do ia_1=1,nucleons-1
               an_st_1=occi(ia_1)
               do ia_2=ia_1+1,nucleons
                  an_st_2=occi(ia_2)
               
                  select case(mt2_sp(an_st_1)
     $                 +mt2_sp(an_st_2))
                  case(2)
                     cr_st_min_1=1
                     cr_st_min_2=2
                     cr_st_max_1=nasps-1
                     cr_st_max_2=nasps
                  case(-2)
                     cr_st_min_1=mxsps+1
                     cr_st_min_2=mxsps+2
                     cr_st_max_1=mxsps+nasps-1
                     cr_st_max_2=mxsps+nasps
                  case(0)
                     cr_st_min_1=1
                     cr_st_min_2=mxsps+1
                     cr_st_max_1=nasps
                     cr_st_max_2=mxsps+nasps
                  case default
                     cycle
                  end select

                  if (nucleons>2) then
                     occim2(1:ia_1-1)=occi(1:ia_1-1)
                     occim2(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                     occim2(ia_2-1:nucleons-2)=occi(ia_2+1:nucleons)
                     
                     locz(1:occim2(1))=0
                     do cr_1=2,nucleons-2
                        locz(occim2(cr_1-1)+1:occim2(cr_1))=cr_1-1
                     end do
                     locz(occim2(nucleons-2)+1:2*mxsps)=nucleons-2
                  else
                     locz=0
                  endif
                  
                  intbas=0
                  do cr_1=1,nucleons-2
                     wd=(occim2(cr_1)-1)/nbit+1
                     ibit=mod(occim2(cr_1)-1,nbit)
c     print *,' i,wd,ibit=',i,wd,ibit
                     intbas(wd)=ibset(intbas(wd),ibit)
                  end do
                  interm_energy=-2*n_sp(an_st_1)-l_sp(an_st_1)
     $                 -2*n_sp(an_st_2)-l_sp(an_st_2)+N_HO_i
                  
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     if (2*n_sp(cr_st_1)+l_sp(cr_st_1)
     +                    +interm_energy>nhw) exit
                     wd=(cr_st_1-1)/nbit+1
                     ibit=mod(cr_st_1-1,nbit)
                     if (btest(intbas(wd),ibit)) cycle
                     do cr_st_2=max(cr_st_min_2,cr_st_1+1),
     $                    cr_st_max_2
                        if (2*n_sp(cr_st_1)+l_sp(cr_st_1)
     $                       +2*n_sp(cr_st_2)+l_sp(cr_st_2)
     +                       +interm_energy>nhw) exit
                        if (m2_sp(cr_st_1)+m2_sp(cr_st_2)
     $                    -m2_sp(an_st_1)-m2_sp(an_st_2)
     +                       /=0) cycle
                        if (mod(l_sp(cr_st_1)+l_sp(an_st_1)
     $                       +l_sp(cr_st_2)+l_sp(an_st_2),
     +                       2)/=0) cycle
                        wd=(cr_st_2-1)/nbit+1
                        ibit=mod(cr_st_2-1,nbit)
                        if (btest(intbas(wd),ibit)) cycle
                        if (nucleons>2) then
                           occf(1:locz(cr_st_1))=occim2(1:locz(cr_st_1))
                           occf(locz(cr_st_1)+1)=cr_st_1
                           occf(locz(cr_st_1)+2:locz(cr_st_2)+1)=
     $                          occim2(locz(cr_st_1)+1:locz(cr_st_2))
                           occf(locz(cr_st_2)+2)=cr_st_2
                           occf(locz(cr_st_2)+3:nucleons)=
     +                          occim2(locz(cr_st_2)+1:nucleons-2)
                        else
                           occf(1)=cr_st_1
                           occf(2)=cr_st_2
                        endif
                        
                        intbasf=0
                        do cr_1=1,nucleons
                           wd=(occf(cr_1)-1)/nbit+1
                           ibit=mod(occf(cr_1)-1,nbit)
c     print *,' i,wd,ibit=',i,wd,ibit
                           intbasf(wd)=ibset(intbasf(wd),ibit)
                        end do

                        call get_state_index(nucleons,2*nwd,intbasf,i3)
                        if (i3==-1) cycle
                        iphase=(-1)**(ia_1+ia_2+locz(cr_st_1)
     $                       +locz(cr_st_2)+1)
                        kk=an_st_1
                        ll=an_st_2
                        ii=cr_st_1
                        jj=cr_st_2
************************************************************
                        
c                     xcoef=real(iphase)*ham(nesd)


c     ----- Arriving here means that the Hamiltonian will potentially connect
c     ----- the initial Slater determinant (ibasis(,i1)) and the final one 
c     ----- (ibasis(,i3)) which differ, at most, by two s.p. states.

c     ----- Generate contributions from all possible two-body operators
c     ----- between these two many-body states:
                        xrm = 0.0
                        xr1 = 0.0
                        xr2 = 0.0
                        xj2 = 0.0
                        xt2 = 0.0
                        gtp = 0.0
                        gtm = 0.0
                        x1bqp = 0.0
                        x1bqn = 0.0
                        x1blp = 0.0
                        x1bln = 0.0
                        x1bsp = 0.0
                        x1bsn = 0.0
c--------------------------------------------------------------------

c            kk=ianstate
c            ll=ianstatep
c            ii=icrstate
c            jj=icrstatep
************************************************************

c     ----- The orbitals kk and ll have been changed to ii and jj 
c     ----- respectively.
c            ipht = locz(kk)+locz(ll)-1+locz(ii)+locz(jj)
c            if(jj.gt.kk)ipht = ipht - 1
c            if(jj.gt.ll)ipht = ipht - 1
c            if(ii.gt.kk)ipht = ipht - 1
c            if(ii.gt.ll)ipht = ipht - 1
c            phti = jsgn(ipht)
                        phti=real(iphase,kind(0.d0))
c     
                        call tbmeop(ii,jj,kk,ll,rmop1,xj2op1,t2op1,
     $                       xgtp,xgtm)
                        xrm = phti * rmop1
                        xj2 = phti * xj2op1
                        xt2 = phti * t2op1
                        gtp = phti * xgtp
                        gtm = phti * xgtm
c     
                        xrmo2 = 0.5*xrm
                        if(ii.le.mxsps.and.jj.le.mxsps)xr1 = xrm
                        if(ii.le.mxsps.and.jj.gt.mxsps)xr1 = xrmo2
                        if(ii.gt.mxsps.and.jj.gt.mxsps)xr2 = xrm
                        if(ii.le.mxsps.and.jj.gt.mxsps)xr2 = xrmo2
                        if(ii.gt.mxsps.and.jj.le.mxsps)then
                           print *,'  Node #',iproc,
     &                          '  Error: ii.gt.mxsps.and.jj.le.mxsps.'
                           call MPI_Abort(icomm,112,ierr)
                           stop 'ii.gt.mxsps.and.jj.le.mxsps.'
                        endif
c************************************************************************

c         case(2)
c     ----- Here the initial and final m.b. states differ by one orbital:

c**** anihilated state ***
                        ob=.false.
                        if (an_st_2==cr_st_2) then
                           kk=an_st_1
                           ii=cr_st_1
                           ob=.true.
c     ----- The magnetic dipole moment:
                           x1blp = x1blp + xlp(ii,kk)
                           x1bsp = x1bsp + xsp(ii,kk)
                           x1bln = x1bln + xln(ii,kk)
                           x1bsn = x1bsn + xsn(ii,kk)
c     ----- The electric quadrupole moment:
                           x1bqp = x1bqp+xqp(ii,kk)
                           x1bqn = x1bqn+xqn(ii,kk)
                        endif
                        if (an_st_2==cr_st_1) then
                           kk=an_st_1
                           ii=cr_st_2
                           ob=.true.
c     phti=-phti
c     ----- The magnetic dipole moment:
                           x1blp = x1blp - xlp(ii,kk)
                           x1bsp = x1bsp - xsp(ii,kk)
                           x1bln = x1bln - xln(ii,kk)
                           x1bsn = x1bsn - xsn(ii,kk)
c     ----- The electric quadrupole moment:
                           x1bqp = x1bqp-xqp(ii,kk)
                           x1bqn = x1bqn-xqn(ii,kk)
                        endif
                        if (an_st_1==cr_st_1) then
                           kk=an_st_2
                           ii=cr_st_2
                           ob=.true.
c     ----- The magnetic dipole moment:
                           x1blp = x1blp + xlp(ii,kk)
                           x1bsp = x1bsp + xsp(ii,kk)
                           x1bln = x1bln + xln(ii,kk)
                           x1bsn = x1bsn + xsn(ii,kk)
c     ----- The electric quadrupole moment:
                           x1bqp = x1bqp+xqp(ii,kk)
                           x1bqn = x1bqn+xqn(ii,kk)
                        endif
                        if (an_st_1==cr_st_2) then
                           kk=an_st_2
                           ii=cr_st_1
                           ob=.true.
c     phti=-phti
c     ----- The magnetic dipole moment:
                           x1blp = x1blp - xlp(ii,kk)
                           x1bsp = x1bsp - xsp(ii,kk)
                           x1bln = x1bln - xln(ii,kk)
                           x1bsn = x1bsn - xsn(ii,kk)
c     ----- The electric quadrupole moment:
                           x1bqp = x1bqp-xqp(ii,kk)
                           x1bqn = x1bqn-xqn(ii,kk)
                        endif
c            kk=ianstate
c**** created state ***
c            ii=icrstate
******************************************************************
                        if (ob) then
c     ----- For one-body operators:

c     ----- The magnetic dipole moment:
c                        x1blp = x1blp + xlp(ii,kk)
c                        x1bsp = x1bsp + xsp(ii,kk)
c                        x1bln = x1bln + xln(ii,kk)
c                        x1bsn = x1bsn + xsn(ii,kk)
c     ----- The electric quadrupole moment:
c                        x1bqp = xqp(ii,kk)
c                        x1bqn = xqn(ii,kk)

c     ----- For two-body operators:

c            do ix = 1, nucleons
c               index2=locan(ix)
c               if (index2==ianstate) cycle
c               jj=index2
c               ll = jj
c               call tbmeop(tadd,xj2op1,t2op1,xgtp,xgtm)
c               xrm = xrm + tadd
c               xj2 = xj2 + xj2op1
c               xt2 = xt2 + t2op1
c               gtp = gtp + xgtp
c               gtm = gtm + xgtm
c               if(ii.le.mxsps.and.jj.le.mxsps)xr1 = xr1+tadd
c               if(ii.gt.mxsps.and.jj.gt.mxsps)xr2 = xr2+tadd
c               if(ii.le.mxsps.and.jj.gt.mxsps.or.
c     &              ii.gt.mxsps.and.jj.le.mxsps)then
c                  taddo2=0.5*tadd
c                  xr1 = xr1+taddo2
c                  xr2 = xr2+taddo2
c               endif
c            enddo
c            ipht = locz(ii) - locz(kk)
c            if(ii.gt.kk)ipht = ipht - 1
c            phti = jsgn(ipht)
c            xrm = phti*xrm
c            xr1 = phti*xr1
c            xr2 = phti*xr2
c            xj2 = phti*xj2
c            xt2 = phti*xt2
c            gtp = phti*gtp
c            gtm = phti*gtm
c
                           phti=phti/real(nucleons-1,kind(0.d0))
                           x1blp = phti*x1blp
                           x1bln = phti*x1bln
                           x1bsp = phti*x1bsp
                           x1bsn = phti*x1bsn
                           x1bqp = phti*x1bqp
                           x1bqn = phti*x1bqn
                        endif
c***************************************************************************
c     
c         case(0)
c     ----- Handle cases where no orbitals differ between the two m.b. basis
c     ----- states.

c     ----- For one-body operators:
c            do l = 1, nucleons
c               ii = locan(l) 
c     ----- Calc. expectation values of one-body operators between ii and ii:
c---------------------------------------------------------------------------
c     ----- The magnetic dipole moment:
c               x1blp = x1blp + xlp(ii,ii)
c               x1bsp = x1bsp + xsp(ii,ii)
c               x1bln = x1bln + xln(ii,ii)
c               x1bsn = x1bsn + xsn(ii,ii)
c---------------------------------------------------------------------------
c     ----- The electric quadrupole moment:
c               x1bqp = x1bqp + xqp(ii,ii)
c               x1bqn = x1bqn + xqn(ii,ii)
c--------------------------------------------------------------------
c     ----- Need the followling lines if subshell occupations are desired.
c           nnpp   = n_sp(ii)
c           llpp   = l_sp(ii)
c           j2pp   = j2_sp(ii)
c           itz2pp = mt2_sp(ii)
c           if(itz2pp.eq.1)then
c              occup(nnpp,llpp,j2pp) 
c    &              = occup(nnpp,llpp,j2pp)+1.0
c           else
c              occun(nnpp,llpp,j2pp) 
c    &              = occun(nnpp,llpp,j2pp)+1.0
c           endif
c--------------------------------------------------------------------
c            enddo
c            do l = 2, nucleons
c               jj = locan(l) 
c               ll = jj
c               do k = 1, l-1
c                  ii = locan(k) 
c                  kk = ii
c                  call tbmeop(tadd,xj2op1,t2op1,xgtp,xgtm)
c                  xrm = xrm + tadd
c                  xj2 = xj2 + xj2op1
c                  xt2 = xt2 + t2op1
c                  gtp = gtp + xgtp
c                  gtm = gtm + xgtm
c                  if(ii.le.mxsps.and.jj.le.mxsps)xr1 = xr1+tadd
c                  if(ii.gt.mxsps.and.jj.gt.mxsps)xr2 = xr2+tadd
c                  if(ii.le.mxsps.and.jj.gt.mxsps)then
c                     taddo2 = 0.5*tadd
c                     xr1 = xr1+taddo2
c                     xr2 = xr2+taddo2
c                  endif
c                  if(ii.gt.mxsps.and.jj.le.mxsps)then
c                     write(6,*)'  Node #',iproc,
c     &                    '  Error: ii.gt.mxsps.and.jj.le.mxsps.'
c                     call MPI_Abort(icomm,113,ierr)
c                     stop 'ii.gt.mxsps.and.jj.le.mxsps.'
c                  endif
c               enddo
c            enddo
c***************************************************************************

c---------------------------------------------------------------------------
c     I. Calculate static properties for the initial and "k1max" final states:
c                     do k1 = -1, k1max, 1
                        xi10 = bmp0(i1)
                        xi30 = bmp0(i3)
c     if(k1.eq.-1)then
                        if(k1.eq.-1.or.ki.eq.(kf+k1+k00))then
                           xfacc = xi10 * xi30
                        else
                           xi1 = bmp(i1,0)
                           xi3 = bmp(i3,0)
                           xfacc = xi1 * xi3
                        endif
c           xfacc = xfac * bmp(i1,k1) * bmp(i3,k1)
c                        xcoefs(k1) = xcoefs(k1) + xfacc*xcoef
                        xrmsm(k1) = xrmsm(k1) + xfacc*xrm
                        xr1sm(k1) = xr1sm(k1) + xfacc*xr1
                        xr2sm(k1) = xr2sm(k1) + xfacc*xr2
                        xj2sm(k1) = xj2sm(k1) + xfacc*xj2
                        xt2sm(k1) = xt2sm(k1) + xfacc*xt2
                        gtpsm(k1) = gtpsm(k1) + xfacc*gtp
                        gtmsm(k1) = gtmsm(k1) + xfacc*gtm
                        x1bqpsm(k1) = x1bqpsm(k1) + xfacc*x1bqp
                        x1bqnsm(k1) = x1bqnsm(k1) + xfacc*x1bqn
                        x1blpsm(k1) = x1blpsm(k1) + xfacc*x1blp
                        x1blnsm(k1) = x1blnsm(k1) + xfacc*x1bln
                        x1bspsm(k1) = x1bspsm(k1) + xfacc*x1bsp
                        x1bsnsm(k1) = x1bsnsm(k1) + xfacc*x1bsn
c--------------------------------------------------------------------
c     Need the followling lines if subshell occupations are desired.
c           do nzh=0,4
c              do lzh=0,8
c                 do j2zh=1,17,2
c                    occupsm(k1,nzh,lzh,j2zh)=
c    &                    occupsm(k1,nzh,lzh,j2zh)
c    &                    +xfacc*occup(nzh,lzh,j2zh)
c                    occunsm(k1,nzh,lzh,j2zh)=
c    &                    occunsm(k1,nzh,lzh,j2zh)
c    &                    +xfacc*occun(nzh,lzh,j2zh)
c                 enddo
c              enddo
c           enddo
c--------------------------------------------------------------------
c           if(i1.eq.i3)ovlp(k1) = ovlp(k1)+bmp(i1,k1)*bmp(i3,k1)
c     II. Calc. transition amplitudes from the g.s. to a few excited states:
c     ( i3 >= i1, i1 runs from 1 to nsd. note that the factor "xfac"
c     is not present in the following equation because bmp(0).ne.bmp(k1).
c           if(k1.eq.-1)then
                        if(k1.eq.-1.or.ki.eq.(kf+k1+k00))then
                                ! In this case, the <i| and |f> are same states
                           t1bqpsm(k1) = x1bqpsm(k1)
                           t1bqnsm(k1) = x1bqnsm(k1)
                           t1blpsm(k1) = x1blpsm(k1)
                           t1blnsm(k1) = x1blnsm(k1)
                           t1bspsm(k1) = x1bspsm(k1)
                           t1bsnsm(k1) = x1bsnsm(k1)
                        else    ! In this case, the <i| and |f> are different.
                           if(i1.eq.i3)then
                              xfacc = xi10*xi1
                           else
                              xfacc = xi10*xi3 !+ xi1*xi30
                           endif
                           t1bqpsm(k1) = t1bqpsm(k1) + xfacc*x1bqp
                           t1bqnsm(k1) = t1bqnsm(k1) + xfacc*x1bqn
                           t1blpsm(k1) = t1blpsm(k1) + xfacc*x1blp
                           t1blnsm(k1) = t1blnsm(k1) + xfacc*x1bln
                           t1bspsm(k1) = t1bspsm(k1) + xfacc*x1bsp
                           t1bsnsm(k1) = t1bsnsm(k1) + xfacc*x1bsn
                        endif
c                     enddo
                     end do
                  end do
               end do
            end do
 5000    continue
         if (k1==k1max) then
            if (allocated(bmp)) deallocate(bmp) 
            if (allocated(bmp0)) deallocate(bmp0) 
         endif
c------------------------------------------------------------------------
c         if(nproc.eq.1)goto 5001   
c------------------------------------------------------------------------
c     ----- Do the Global sum here for the various observables: 
      
         if (nproc>1) then
            call MPI_Barrier(icomm,ierr)

c         k1max2 = k1max + 2
c     call MPI_REDUCE(ovlp(-1),ztmp(-1),k1max2,
c    &     MPI_REAL8,MPI_SUM,0,icomm,ierr)
c     do k1=-1,k1max,1
c        ovlp(k1) = ztmp(k1)
c     enddo

c      call MPI_Reduce(xcoefs(-1),ztmp(-1),k1max2,
c     &     MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
c         xcoefs(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(xrmsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c         do k1=-1,k1max,1
            xrmsm(k1) = ztmp(k1)
c         enddo

            call MPI_Reduce(xr1sm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c         do k1=-1,k1max,1
            xr1sm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(xr2sm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            xr2sm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(xj2sm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            xj2sm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(xt2sm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c     do k1=-1,k1max,1
            xt2sm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(gtpsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            gtpsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(gtmsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            gtmsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1bqpsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1bqpsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1bqnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1bqnsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1blpsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1blpsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1blnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1blnsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1bspsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1bspsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(x1bsnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            x1bsnsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1bqpsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1bqpsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1bqnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1bqnsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1blpsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1blpsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1blnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1blnsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1bspsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1bspsm(k1) = ztmp(k1)
c      enddo

            call MPI_Reduce(t1bsnsm(k1),ztmp(k1),1,
     &           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c      do k1=-1,k1max,1
            t1bsnsm(k1) = ztmp(k1)
c      enddo

         endif

c------------------------------------------------------------------------
c 5001 if(iproc.gt.0)return
c      xcoefs(-1:k1max)=energy(1:k1max+2)
         if (k1==-1) then
            xcoefs(-1)=energy(ki)
         else
c      do k1 = 0, k1max
            xcoefs(k1)=energy(kf+k1+k00)
c      end do
         endif

         if (iproc==0) then
c****pn*********
            if (k1==k1max) then
               write(13) (xcoefs(iiy),iiy=-1,k1max)
            endif
c***************

            write(8,*)'  Completed 5000 loop for k1=',k1

c      do k1 = -1, k1max, 1
            if(k1.eq.-1) kfx = ki
            if(k1.ge.0) kfx = kf + k1 + k00
            if(egs.eq.0.0)then
               exeng(k1) = xcoefs(k1) - xcoefs(-1)
            else
               exeng(k1) = xcoefs(k1) - egs
            endif
            write(8,*)'State # ',kfx,'   E =',xcoefs(k1)
            write(8,*)'  xj2sm =',xj2sm(k1),'  xt2sm = ',xt2sm(k1)
            write(8,*)'  xr1sm =',xr1sm(k1),'  xr2sm =',xr2sm(k1),
     &           '  xrmsm =',xrmsm(k1)
         if(xj2sm(k1).lt.-0.25d0.or.xt2sm(k1).lt.-0.25d0)then
            write(8,*)'  Error in J, T computation in TRANSITIONS.'
            write(8,*)'  xj2sm, xt2sm: ',xj2sm(k1), xt2sm(k1)
            call MPI_Abort(icomm,115,ierr)
            stop 'Error in J, T computation in TRANSITIONS.'
         endif
 5010    xj = -0.5 + dsqrt(0.25d0 + xj2sm(k1))
         xt = -0.5 + dsqrt(0.25d0 + xt2sm(k1))
         write(8,*)' State # ',kfx,
     &        ' E =',xcoefs(k1),'  J =',xj,'  T =',xt
         xjj(k1) = xj
         xtt(k1) = xt

c****pn*********
         if (k1==k1max) then
            do iiy=-1,k1max
               write(13) iiy,xjj(iiy),xtt(iiy)
            end do
         endif
c***************

         if(nprotons.ne.0)xr1sm(k1) = 2.d0*real(nucleons)
     &        /real(nprotons)*xr1sm(k1)-xrmsm(k1)
         if(nneutrns.ne.0)xr2sm(k1) = 2.d0*real(nucleons)
     &        /real(nneutrns)*xr2sm(k1)-xrmsm(k1)
         ibad = 1
         if(xr1sm(k1).ge.0.d0.and.xr2sm(k1).ge.0.d0
     &        .and.xrmsm(k1).ge.0.d0)
     &        goto 5015
         write(8,*)
     &        '  Error: bad sums ',xr1sm(k1),xr2sm(k1),xrmsm(k1)
         ibad = -1
         xr1sm(k1) = dabs(xr1sm(k1))
         xr2sm(k1) = dabs(xr2sm(k1))
         xrmsm(k1) = dabs(xrmsm(k1))
 5015    r1s = ibad*dsqrt(xr1sm(k1))
         r2s = ibad*dsqrt(xr2sm(k1))
         rms = ibad*dsqrt(xrmsm(k1))
         write(8,*)' Rp =',r1s,' Rn =',r2s,' Rm =',rms
         jx2 = 2*(xj+0.1)
         mx2 = mjtotal
         itx2 = 2*(xt+0.1)
         itz2 = mttotal
         if(k1.eq.-1)then
            jx21=jx2
            itx21=itx2
         endif

c     ----- Calculate magnetic dipole moment:
         t3jd = threej(jx2,2,jx2,-mx2,0)
         if(abs(t3jd).lt.1.e-6)then
            x1bdsm = -99.0
         else
            rt3jd = jsgn((jx2-mx2)/2)
     &           *threej(jx2,2,jx2,-jx2,0)/t3jd
            x1blpsm(k1) = x1blpsm(k1)*rt3jd
            x1blnsm(k1) = x1blnsm(k1)*rt3jd
            x1bspsm(k1) = x1bspsm(k1)*rt3jd
            x1bsnsm(k1) = x1bsnsm(k1)*rt3jd
            x1bdsm = glp*x1blpsm(k1) + gln*x1blnsm(k1)
     &           + gsp*x1bspsm(k1) + gsn*x1bsnsm(k1)
         endif
         write(8,*)'M1s:',' lp=',x1blpsm(k1),
     &        ' ln=',x1blnsm(k1),' sp=',x1bspsm(k1),
     &        ' sn=',x1bsnsm(k1),' M1=',x1bdsm

c     ----- Calculate electric quadrupole moment:
         t3jq = threej(jx2,4,jx2,-mx2,0)
         if(abs(t3jq).lt.1.e-6)then
            x1bqsm = -99.0
         else
            rt3jq = real(jsgn((jx2-mx2)/2))
     &           *threej(jx2,4,jx2,-jx2,0)/t3jq
            x1bqpsm(k1) = x1bqpsm(k1)*rt3jq*bsquare
            x1bqnsm(k1) = x1bqnsm(k1)*rt3jq*bsquare
            x1bqsm = epi*x1bqpsm(k1) + enu*x1bqnsm(k1)
         endif
         write(8,*)'Qs:',' p=',x1bqpsm(k1),
     &        ' n=',x1bqnsm(k1),' Q=',x1bqsm

c     ----- Calculate B(M1) and B(GT):
         t3jm = threej(jx2,2,jx21,-mx2,0)
         if(abs(t3jm).lt.1.e-6)then
            t1bdsm = -99.0
            t1gtp = -99.0
            t1gtm = -99.0
         else
            t3jm = real(jsgn((jx2-mx2)/2))/t3jm
            t1blpsm(k1) = t1blpsm(k1)*t3jm
     &           *sqrt(3./(4.*3.1415926536*(jx21+1.)))
            t1blnsm(k1) = t1blnsm(k1)*t3jm
     &           *sqrt(3./(4.*3.1415926536*(jx21+1.)))
            t1bspsm(k1) = t1bspsm(k1)*t3jm
     &           *sqrt(3./(4.*3.1415926536*(jx21+1.)))
            t1bsnsm(k1) = t1bsnsm(k1)*t3jm
     &           *sqrt(3./(4.*3.1415926536*(jx21+1.)))
            t1bdsm = glp*t1blpsm(k1) + gln*t1blnsm(k1)
     &           + gsp*t1bspsm(k1) + gsn*t1bsnsm(k1)
c           t1bdsm = 3./(4.*3.1415926536*(jx21+1.))*t1bdsm**2
            t1bdsm = t1bdsm**2
            write(8,*)' lp=',t1blpsm(k1),' ln=',t1blnsm(k1),
     &           ' sp=',t1bspsm(k1),' sn=',t1bsnsm(k1),
     &           ' B(M1)=',t1bdsm
            t3jg = threej(itx2,2,itx21,-itz2,0)
            t1gtp = -99.0
            t1gtm = -99.0
            if(abs(t3jg).gt.1.e-6)then
               t1gt = (t1bspsm(k1) - t1bsnsm(k1))
     &              /sqrt(3./(4.*3.1415926536*(jx21+1.)))
c     -----    The factor sqrt(3./(4.*3.1415926536)) was multiplied previously
c     -----    to get t1bspsm(k1) and t1bsnsm(k1) for M1. So we have to 
c     -----    divide this factor now to make things even.

               rt3jp = -threej(itx2,2,itx21,-(itz2+2),2)/t3jg
               rt3jm = -threej(itx2,2,itx21,-(itz2-2),-2)/t3jg
               t1gtp = t1gt*rt3jp
               t1gtm = t1gt*rt3jm
               t1gtp = 2.*t1gtp**2/real(jx21+1)
               t1gtm = 2.*t1gtm**2/real(jx21+1)
               write(8,*)' proton=',t1gtp,'  neutron=',t1gtm
            endif
         endif

c     ----- Calculate B(E2):
         t3je = threej(jx2,4,jx21,-mx2,0)
         if(abs(t3je).lt.1.e-6)then
            t1bqsm = -99.0
         else
            t3je = real(jsgn((jx2-mx2)/2))/t3je
            t1bqpsm(k1) = t1bqpsm(k1)*t3je*bsquare
     &           *sqrt(5./(16.*3.1415926536*(real(jx21+1))))
            t1bqnsm(k1) = t1bqnsm(k1)*t3je*bsquare
     &           *sqrt(5./(16.*3.1415926536*(real(jx21+1))))
            t1bqsm = epi*t1bqpsm(k1) + enu*t1bqnsm(k1)
c           t1bqsm = 5./(16.*3.1415926536*
c    &           (real(jx21+1)))*t1bqsm**2
            t1bqsm = t1bqsm**2
            write(8,*)' proton=',t1bqpsm(k1),
     &           '  neutron=',t1bqnsm(k1),' B(E2)=',t1bqsm
         endif
         inquire(unit=675,opened=kappadepwr)
         if (kappadepwr) then
            write(675,'(e12.5,1x,f10.4,1x,f8.4,1x,f8.4,i4,i10)') kappa,
     $           xcoefs(k1),xj,xt,kfx,nsd
         endif
         write(9,*)
         write(9,6490)kfx,xcoefs(k1),xj,xt
 6490    format(1x,'State #',i2,'   Energy =',f10.4,
     &        '     J =',f8.4,'      T =',f8.4)
         if(x1bdsm.ne.-99.0)write(9,6501)x1blpsm(k1),
     &        x1blnsm(k1),x1bspsm(k1),x1bsnsm(k1),x1bdsm
 6501    format(2x,'M1:    ',' pL=',f9.4,' nL=',f9.4,
     &        ' pS=',f9.4,' nS=',f9.4,'  M1   =',f9.4)
         if(t1bdsm.ne.-99.0)write(9,6502)t1blpsm(k1),
     &        t1blnsm(k1),t1bspsm(k1),t1bsnsm(k1),t1bdsm
 6502    format(2x,'B(M1): ',' pL=',f9.4,' nL=',f9.4,
     &        ' pS=',f9.4,' nS=',f9.4,'  B(M1)=',f9.4)
         if(x1bqsm.ne.-99.0)write(9,6503)x1bqpsm(k1),
     &        x1bqnsm(k1),x1bqsm
 6503    format(2x,'Q:     ',' proton=',f9.4,
     &        '  neutron=',f9.4,18x,'Q    =',f9.4)
         if(t1bqsm.ne.-99.0)write(9,6504)t1bqpsm(k1),
     &        t1bqnsm(k1),t1bqsm
 6504    format(2x,'B(E2): ',' proton=',f9.3,
     &        '  neutron=',f9.3,18x,'B(E2)=',f9.3)
         write(9,6505)r1s,r2s,rms
 6505    format(2x,'Radius:',' proton=',f9.4,
     &        '  neutron=',f9.4,18x,'mass =',f9.4)
         if(t1gtp.ne.-99.0)write(9,6506)t1gtp,t1gtm 
 6506    format(2x,'B(GT): ',' B(GT+)=',f9.4,'  B(GT-) =',f9.4)
         if(igt.eq.1)write(9,6507)gtpsm(k1),gtmsm(k1)
 6507    format(2x,'Sum GT:',' S(GT+)=',f9.4,'  S(GT-) =',f9.4)
         write(9,*)
c         write(9,*)'    Configurations with probability > 1%:'
c         write(9,6508)(ishl,ishl=1,nshll),(ishl,ishl=1,nshll)
c 6508    format(20x,20(i3))
         
c         do ishl=1,nshll
c            occp(ishl)=0.0
c            occn(ishl)=0.0
c         enddo

c         do iconf = 1, nconf

c            do ishl=1,nshll
c               occp(ishl)=occp(ishl)+
c     &              prob(iconf,k1)*nprtn(iconf,ishl)
c               occn(ishl)=occn(ishl)+
c     &              prob(iconf,k1)*nneut(iconf,ishl)
c            enddo

c            if(prob(iconf,k1).gt.1.0)
c     &           write(9,6509)prob(iconf,k1),
c     &           (nprtn(iconf,ishl),ishl=1,nshll),
c     &           (nneut(iconf,ishl),ishl=1,nshll)
c 6509       format(9x,f7.3,'%:',2x,20(i3))
c         enddo
         
c         write(9,*)
c         write(9,6510)(0.01*occp(ishl),ishl=1,nshll)
c         write(9,6511)(0.01*occn(ishl),ishl=1,nshll)
c 6510    format('   protons:  ',8f8.4)
c 6511    format('   neutrons: ',8f8.4)
         endif

c*** k1 eigenstate loop
      enddo
c***
c***************************
c      write(13) nconf
c      write(13) nprtn,nneut
c      write(13) iendconf
c***************************
c--------------------------------------------------------------------
      if (iproc==0) then
         write(9,*)
         write(9,*)' The energy spectrum is:'
         write(9,*)'  n    J    T        Ex(MeV)'
         do k1 = -1, k1max, 1
            if(k1.eq.-1) kfx = ki
            if(k1.ge.0) kfx = kf + k1 + k00
            write(9,6512)kfx,xjj(k1),xtt(k1),exeng(k1)
 6512       format(1x,i3,2f6.2,f12.4)
         enddo
         call clock(cputime,walltime,0,0)
         timeused = cputime - timeused
c*      timeused = dclock( ) - timeused
         write(8,*)'  CPU time used after the final iteration (sec):',
     &        timeused
      endif
      end


c     ----- This subroutine calculates one-body matrix elements for the
c     ----- one-body operators relevant to M1 and E2.

      subroutine obme
      use parameters
      use nodeinfo
      use jsgna
      use consts
      use spb
      use spsdata
      use obmes
      use r2HOme
      include 'mpif.h'
c
      xxxxx = 1.0/real(nucleons-1) - 1.0
      allocate(xlp(mxsps2,mxsps2),xln(mxsps2,mxsps2),
     +     xsp(mxsps2,mxsps2),xsn(mxsps2,mxsps2),
     +     xqp(mxsps2,mxsps2),xqn(mxsps2,mxsps2))
c     
c     ----- Read in matrix elements of r**2 to be used when calculating
c     ----- B(E2) and Q:
      call r2matel(mxn2l)
      
      do 200 ii1 = 1, mxsps2, 1
         if(ii1.gt.nasps.and.ii1.le.mxsps)goto 200
         if(ii1.gt.(mxsps+nasps))goto 200
         do 100 kk1 = ii1, mxsps2, 1
            if(kk1.gt.nasps.and.kk1.le.mxsps)goto 100
            if(kk1.gt.(mxsps+nasps))goto 100

            xlp0 = 0.0
            xln0 = 0.0
            xsp0 = 0.0
            xsn0 = 0.0
            xqp0 = 0.0
            xqn0 = 0.0

            nnpp   = n_sp(ii1)
            llpp   = l_sp(ii1)
            j2pp   = j2_sp(ii1)
            m2pp   = m2_sp(ii1)
            itz2pp = mt2_sp(ii1)
            nnqq   = n_sp(kk1)
            llqq   = l_sp(kk1)
            j2qq   = j2_sp(kk1)
            m2qq   = m2_sp(kk1)
            itz2qq = mt2_sp(kk1)

c     ----- For the magnetic dipole moment:
c     ----- The M1 matrix element is nonzero when states ii1 and kk1
c     ----- differ only by their j values (j2pp and j2qq):
            if(nnpp.ne.nnqq.or.llpp.ne.llqq.or.m2pp.ne.m2qq
     &           .or.itz2pp.ne.itz2qq)goto 25
            do m2sz = -1,1,2
               m2lz = m2pp-m2sz
               if(iabs(m2lz).le.2*llpp)then
                  cgcg = sqrt((j2pp+1.)*(j2qq+1.))
     &                 *threej(2*llpp,1,j2pp,m2lz,m2sz)
     &                 *threej(2*llpp,1,j2qq,m2lz,m2sz)
                  if(itz2pp.eq.1)then
                     xlp0 = xlp0 + cgcg*m2lz*0.5
                     xsp0 = xsp0 + cgcg*m2sz*0.5
                  else
                     xln0 = xln0 + cgcg*m2lz*0.5
                     xsn0 = xsn0 + cgcg*m2sz*0.5
                  endif
               endif
            enddo
 25         continue
c---------------------------------------------------------------------------
c     ----- For the electric quadrupole moment:
            if(m2pp.ne.m2qq.or.itz2pp.ne.itz2qq
     &           .or.iabs(j2pp-j2qq).gt.4
     &           .or.(j2pp+j2qq).lt.4)goto 50
c     ----- y2=reduced matrix element of y2 between |l1,1/2,j1> and |l2,1/2,j2>
            y2 = jsgn((j2pp+j2qq+m2pp+1)/2)*(1+jsgn(llpp+llqq))
            if(y2.eq.0.)goto 50
            r2me = r2mea(nnpp,llpp,nnqq,llqq)
            if(r2me.eq.0.)goto 50
            y2 = y2*sqrt((j2pp+1.)*(j2qq+1.))
     &           *threej(j2qq,4,j2pp,m2pp,0)
     &           *threej(j2pp,4,j2qq,1,0)
            if(itz2pp.eq.1)then
               xqp0 = y2*r2me
            else
               xqn0 = y2*r2me
            endif
 50         continue
            xlp(ii1,kk1) = xlp0
            xsp(ii1,kk1) = xsp0
            xln(ii1,kk1) = xln0
            xsn(ii1,kk1) = xsn0
            xqp(ii1,kk1) = xqp0
            xqn(ii1,kk1) = xqn0
c
            xlp(kk1,ii1) = xlp0
            xsp(kk1,ii1) = xsp0
            xln(kk1,ii1) = xln0
            xsn(kk1,ii1) = xsn0
            xqp(kk1,ii1) = xqp0
            xqn(kk1,ii1) = xqn0
 100     continue
 200  continue
      deallocate(r2mea)
      return
      end

      subroutine r2matel(mxn2l)
      use gamalog
      use r2HOme
      use nodeinfo
      implicit double precision (a-h,o-z)
      parameter(nstep=6000,rc=1.d-6,rinf=22.d0)
      double precision,allocatable,dimension(:,:,:):: u
      nhom=mxn2l
      lrelm=nhom
      nrmax=nhom/2
      anu=1.d0
c**** call gamasub before first use
      call gamasub

c     MKGK
      if (allocated(u)) deallocate(u)
      allocate(u(0:nstep,0:nrmax,0:lrelm))
      rstep=(rinf-rc)/dble(nstep)

      do lr=0,lrelm
         do nr=0,nrmax
            do i=0,nstep
               rr=rc+i*rstep
               call waver(nr,lr,anu,rr,wave)
               u(i,nr,lr)=wave
            end do
         end do
      end do   
      if (iproc==0) print *,' u calculated'

      if (allocated(r2mea)) deallocate(r2mea)
      allocate(r2mea(0:nrmax,0:lrelm,0:nrmax,0:lrelm))
      r2me=0.d0
      do lra=0,lrelm
         do nra=0,nrmax
            if (2*nra+lra>nhom) cycle
            do lrb=0,lrelm
               if ((-1)**(lra+lrb)/=1) cycle
               do nrb=nra,nrmax
                  if (2*nrb+lrb>nhom) cycle
                  over=0.d0
                  fact=4.d0
                  do i=1,nstep-1
                     rr=rc+i*rstep
                     over=over+rr*rr*u(i,nra,lra)*u(i,nrb,lrb)*fact
                     fact=6.d0-fact
                  end do
                  r2mea(nra,lra,nrb,lrb)=(over
     +  +rc*rc*u(0,nra,lra)*u(0,nrb,lrb)
     +  +rinf*rinf*u(nstep,nra,lra)*u(nstep,nrb,lrb))*rstep/3.d0
                  r2mea(nrb,lrb,nra,lra)=r2mea(nra,lra,nrb,lrb)
               end do
            end do
         end do
      end do   
         
      end


      subroutine gamasub
      use gamalog
      implicit double precision(a-h,o-z)
      maxgam=2*nrmax+2*lrelm+3
      allocate(dsq(nrmax,0:lrelm))
      allocate(gamal(maxgam))
      gamal(2)=0.d0
      gamal(1)=0.5d0*dlog(3.14159265358979312d0)
      do i1=3,maxgam
	 gamal(i1)=dlog(dble(float(i1))/2.d0-1.d0)+gamal(i1-2)
      end do
      do il=0,lrelm
         do in=1,nrmax
            dsq(in,il)=dsqrt(dble(in)*(dble(il+in)+0.5d0))
         end do
      end do
      end


      subroutine waver(n,l,anu,r,wave)
      use gamalog
      implicit double precision (a-h,o-z) 

      dlanu=dlog(anu)
      zz=anu*r*r
      wavel=0.25d0*dlanu-zz/2.d0
     &     +dble(l+1)*(0.5d0*dlanu+dlog(r))
      if (n==0) then
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
      elseif (n==1) then
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5)))
     +                   *(dble(l)+1.5d0-zz)         
      else   
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5)))
     +                   *(dble(l)+1.5d0-zz)         
         a=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
         do nnn=2,n
            b=((dble(l+2*nnn)-0.5d0-zz)*guerp
     +           -a*dsq(nnn-1,l))
     +           /dsq(nnn,l)
            a=guerp
            guerp=b
         end do   
      endif 
      wave=dexp(wavel)*guerp
      end 


c     ----- This subroutine returns two-body matrix elements 
c     ----- rmop, xj2op, t2op etc.
      subroutine tbmeop(ii,jj,kk,ll,rmox,xj2ox,t2ox,xgtp,xgtm)
      use parameters
      use jsgna
      use consts
      use spb
      use spsdata
      use tbmes
      use the3js
      use TUD_tbme, only: tbmeTUD
      use pn_tbme, only: tbmepn
c      use hamc
      integer,intent(IN) :: ii,jj,kk,ll
      real(kind(0.d0)) rmo
c     
      rmox = 0.0
      xj2ox= 0.0
      t2ox = 0.0
      xgtp = 0.0
      xgtm = 0.0
c
      m2a  = m2_sp(ii)
      m2b  = m2_sp(jj)
      m2c  = m2_sp(kk)
      m2d  = m2_sp(ll)
      mmjj = m2a + m2b
c
      itz2a = mt2_sp(ii)
      itz2b = mt2_sp(jj)
      itz2c = mt2_sp(kk)
      itz2d = mt2_sp(ll)
      mmtt = itz2a + itz2b
c
c     ----- To be used in the3j(sps1,sps2,J)
c     ----- sps1 and sps2 run from 1 to nasps.
      ii0 = ii + (itz2a-1)/2*mxsps
      jj0 = jj + (itz2b-1)/2*mxsps
      kk0 = kk + (itz2c-1)/2*mxsps
      ll0 = ll + (itz2d-1)/2*mxsps
c
      n1 = nobt_sp(ii)
      n2 = nobt_sp(jj)
      n3 = nobt_sp(kk)
      n4 = nobt_sp(ll)
c----------------------------------------------------
c     ----- Perform the uncoupling transformation:
      nna   = n_sp(ii)
      nnb   = n_sp(jj)
      nnc   = n_sp(kk)
      nnd   = n_sp(ll)
      j2a   = j2_sp(ii)
      j2b   = j2_sp(jj)
      j2c   = j2_sp(kk)
      j2d   = j2_sp(ll)
      lla   = l_sp(ii)
      llb   = l_sp(jj)
      llc   = l_sp(kk)
      lld   = l_sp(ll)
c
c     ----- For rmox, xj2ox, t2ox:
c
      jjamb = (j2a-j2b)/2
      jjcmd = (j2c-j2d)/2
      jmin = max(iabs(jjamb),iabs(jjcmd))
      phsj = real(jsgn(jjamb+jjcmd))
      jjcd = (j2c + j2d)/2
      jmax = min((j2a+j2b)/2,jjcd)
      iprty = (1-(-1)**(lla+llb))/2
c
      do 100 jval = jmin, jmax, 1
         jjval = jval + jval
         if(iabs(mmjj).gt.jjval)goto 100
         t3j = the3j(ii0,jj0,jval)*the3j(kk0,ll0,jval)
         if(t3j.eq.0.0)goto 100
         do 50 it = 0,1
            iit = it + it
            if(iabs(mmtt).gt.iit)goto 50
            if((n1.eq.n2.or.n3.eq.n4)
     &           .and.jsgn(jval+it).eq.1)goto 50
            t3t = the3t(it,itz2a,itz2b)*the3t(it,itz2c,itz2d)
            if(t3t.eq.0.0)goto 50
c
            call findindx(n1,n2,n3,n4,jval,it,indx,phse)
            fac0 = phsj*real((iit+1)*(jjval+1))*t3t*t3j
            if (.not.tbmeTUD.and..not.tbmepn) then
               rmox = rmox + fac0 * rful(indx)*phse
            endif
c     ----- For J**2 and T**2:
            xjtox = 0.0
            if(n1.eq.n3.and.n2.eq.n4) xjtox = 1.0
            if(n1.eq.n4.and.n2.eq.n3) xjtox = xjtox 
     &           + real(jsgn(jjcd+jval+it))
            if(xjtox.ne.0.0)then
               xj2ox = xj2ox + xjtox*fac0*(real(jval*(jval+1))
     &              + xxxxx*0.25*real(j2a*(j2a+2)+j2b*(j2b+2)))
               t2ox  = t2ox + xjtox*fac0*(real(it*(it+1))
     &              + xxxxx*1.5)
            endif
 50      continue
 100  continue
      if (tbmepn) then
         call twobody_op_pn(ii,jj,kk,ll,rmo)
         rmox=rmo
      elseif (tbmeTUD) then
         call twobody_op_TUD(ii,jj,kk,ll,rmo)
         rmox=rmo
      else
         if(n1.eq.n2)rmox=rmox*sqrt2
         if(n3.eq.n4)rmox=rmox*sqrt2
c     ----- This is not done for J**2 and T**2 because the TBMEs obtained
c     ----- in the above for J**2 and T**2 are not normalized.
      endif

      if(igt.eq.0)goto 500      ! If (igt=0), do not calculate (GT+/-)^2.
c
c     ----- For (gt+ and gt- summed strengths) xgtp and xgtm:
c     ----- Assuming only protons and neutrons are in the system. That is to
c     ----- say, both the spin and isopin of the particles are 1/2.
c     ----- gt+ = sum_{i,j} si.sj ti(-)tj(+)
c     -----     = sum_{i}si^2 ti(-)ti(+) + 2sum_{i<j} si.sj ti(-)tj(+)
c     -----     = sum_{i<j} 3[ti(-)ti(+) + tj(-)tj(+)]/(a-1)   <-- diagonal
c     -----        +2sum_{i<j} si.sj ti(-)tj(+)             <-- off diagonal
c
c
c     ----- (ab| |cd) = 0.5*<ab-ba| |cd-dc>
c                     = 0.5[<ab| |cd> + <ba| |dc> - <ab| |dc> - <ba| |cd>]
c     
      if(ii.eq.ll.and.jj.eq.kk)
     &     write(8,*)'  Node #',iproc,'  Warning: ii=ll and jj=kk'
      if(ii.eq.kk.and.jj.eq.ll)then
c     ----- diagonal contributions:
         xgtp = 1.5*(2.-mmtt)/real(nucleons-1)
         xgtm = 1.5*(2.+mmtt)/real(nucleons-1)
      endif
c     
      if(mmtt.ne.0)goto 500
c     ----- Arriving here means that |ii,jj> and |kk,ll> have the following 
c     ----- structure:  |pn> or |np>.
      l2a=lla+lla
      l2b=llb+llb
      l2c=llc+llc
      l2d=lld+lld
      xgts = 0.0
c
      if(itz2a.eq.itz2c)goto 420
c     ----- <pn|t1(-)t2(+)|pn> = <np|t1(+)t2(-)|np> = 0
c     ----- arriving here means that <ii,jj| |kk,ll> has the following 
c     ----- structure: <np|t-(1)t+(2)|pn> for gt+
c     -----         or <pn|t+(1)t-(2)|np> for gt-.
      if(nna.ne.nnc.or.lla.ne.llc)goto 420
      if(nnb.ne.nnd.or.llb.ne.lld)goto 420

c     ----- For <ab| |cd>
c     ----- evaluate contributions from the "direct" matrix element
c     ----- the following loops are for the evaluation of sigma_1*sigma_2
      do 414 msa = -1,1,2
         mla = m2a - msa
         if(iabs(mla).gt.l2a)goto 414
         cgaa = cgcf(l2a,mla,1,msa,j2a)
         mlc = mla
         if(iabs(mlc).gt.l2c)goto 414
         msc = m2c - mlc
         if(iabs(msc).gt.1)goto 414
         cgcc = cgcf(l2c,mlc,1,msc,j2c)

         do 412 msb = -1,1,2
            mlb = m2b - msb
            if(iabs(mlb).gt.l2b)goto 412
            cgbb = cgcf(l2b,mlb,1,msb,j2b)
            mld = mlb
            if(iabs(mld).gt.l2d)goto 412
            msd = m2d - mld
            if(iabs(msd).gt.1)goto 412
            if((msa+msb).ne.(msc+msd))goto 412
            cgdd = cgcf(l2d,mld,1,msd,j2d)
            if(msa.eq.msc)then
               xgtx = 2.*msa*msb
            else
               xgtx = 4.0
            endif
            xgts = xgts + xgtx*cgaa*cgbb*cgcc*cgdd
 412     continue
 414  continue
c
 420  if(itz2a.eq.itz2d)goto 430
      if(nna.ne.nnd.or.lla.ne.lld)goto 430
      if(nnb.ne.nnc.or.llb.ne.llc)goto 430
c     ----- For <ab| |dc>
c     ----- Evaluate contributions from the "exchange" matrix element
c     ----- the following loops are for the evaluation of sigma_1*sigma_2
      do 424 msa = -1,1,2
         mla = m2a - msa
         if(iabs(mla).gt.l2a)goto 424
         cgaa = cgcf(l2a,mla,1,msa,j2a)
         mld = mla
         if(iabs(mld).gt.l2d)goto 424
         msd = m2d - mld
         if(iabs(msd).gt.1)goto 424
         cgdd = cgcf(l2d,mld,1,msd,j2d)
         do 422 msb = -1,1,2
            mlb = m2b - msb
            if(iabs(mlb).gt.l2b)goto 422
            cgbb = cgcf(l2b,mlb,1,msb,j2b)
            mlc = mlb
            if(iabs(mlc).gt.l2c)goto 422
            msc = m2c - mlc
            if(iabs(msc).gt.1)goto 422
            if((msa+msb).ne.(msd+msc))goto 422
            cgcc = cgcf(l2c,mlc,1,msc,j2c)
            if(msa.eq.msd)then
               xgtx = 2.*msa*msb
            else
               xgtx = 4.0
            endif
            xgts = xgts - xgtx*cgaa*cgbb*cgdd*cgcc
 422     continue
 424  continue
 430  if(itz2a.eq.-1)then
         xgtp = xgtp + 0.5*xgts
      else
         xgtm = xgtm + 0.5*xgts
      endif
c
      xgts = 0.0
      if(itz2b.eq.itz2d)goto 440
      if(nnb.ne.nnd.or.llb.ne.lld)goto 440
      if(nna.ne.nnc.or.lla.ne.llc)goto 440
c
c     ----- For <ba| |dc>:
c     ----- evaluate contributions from the "direct" matrix element
c     ----- the following loops are for the evaluation of sigma_1*sigma_2
      do 434 msb = -1,1,2
         mlb = m2b - msb
         if(iabs(mlb).gt.l2b)goto 434
         cgbb = cgcf(l2b,mlb,1,msb,j2b)
         mld = mlb
         if(iabs(mld).gt.l2d)goto 434
         msd = m2d - mld
         if(iabs(msd).gt.1)goto 434
         cgdd = cgcf(l2d,mld,1,msd,j2d)
         do 432 msa = -1,1,2
            mla = m2a - msa
            if(iabs(mla).gt.l2a)goto 432
            cgaa = cgcf(l2a,mla,1,msa,j2a)
            mlc = mla
            if(iabs(mlc).gt.l2c)goto 432
            msc = m2c - mlc
            if(iabs(msc).gt.1)goto 432
            if((msb+msa).ne.(msd+msc))goto 432
            cgcc = cgcf(l2c,mlc,1,msc,j2c)
            if(msb.eq.msd)then
               xgtx = 2.*msb*msa
            else
               xgtx = 4.0
            endif
            xgts = xgts + xgtx*cgbb*cgaa*cgdd*cgcc
 432     continue
 434  continue
c
 440  if(itz2b.eq.itz2c)goto 450
      if(nnb.ne.nnc.or.llb.ne.llc)goto 450
      if(nna.ne.nnd.or.lla.ne.lld)goto 450

c     ----- For <bc| |ad>:
c     ----- Evaluate contributions from the "exchange" matrix element.
c     ----- The following loops are for the evaluation of sigma_1*sigma_2:
      do 444 msb = -1,1,2
         mlb = m2b - msb
         if(iabs(mlb).gt.l2b)goto 444
         cgbb = cgcf(l2b,mlb,1,msb,j2b)
         mlc = mlb
         if(iabs(mlc).gt.l2c)goto 444
         msc = m2c - mlc
         if(iabs(msc).gt.1)goto 444
         cgcc = cgcf(l2c,mlc,1,msc,j2c)
         do 442 msa = -1,1,2
            mla = m2a - msa
            if(iabs(mla).gt.l2a)goto 442
            cgaa = cgcf(l2a,mla,1,msa,j2a)
            mld = mla
            if(iabs(mld).gt.l2d)goto 442
            msd = m2d - mld
            if(iabs(msd).gt.1)goto 442
            if((msb+msa).ne.(msc+msd))goto 442
            cgdd = cgcf(l2d,mld,1,msd,j2d)
            if(msb.eq.msc)then
               xgtx = 2.*msb*msa
            else
               xgtx = 4.0
            endif
            xgts = xgts - xgtx*cgbb*cgaa*cgcc*cgdd
 442     continue
 444  continue
c
 450  if(itz2b.eq.-1)then
         xgtp = xgtp + 0.5*xgts
      else
         xgtm = xgtm + 0.5*xgts
      endif
 500  continue
      return
      end

      subroutine twobody_op_TUD(ii_in,jj_in,kk_in,ll_in,rmo)
      use spsdata
      use paramdef
      use TUD_tbme
      use v3b,only: cgj12
      use tbmes, only: rful
      implicit none
      integer,intent(IN) :: ii_in,jj_in,kk_in,ll_in
      real(kind(0.d0)),intent(OUT) :: rmo
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd,J12
      integer :: ia,ib,ic,id,iphase,temp,i_ab,i_cd,ii,jj,kk,ll
      integer :: m2a,m2b,m2c,m2d
      integer :: itz2a,itz2b,itz2c,itz2d
      real(kind(0.d0)) :: Rabcd,clbab,clbcd,
     $     clbabt(0:1,-1:1),clbcdt(0:1,-1:1)

      ii=ii_in
      jj=jj_in
      kk=kk_in
      ll=ll_in
      
      iphase=1
      
      na=n_sp(ii)
      la=l_sp(ii)
      ja=j2_sp(ii)
      nb=n_sp(jj)
      lb=l_sp(jj)
      jb=j2_sp(jj)
      ia=nlj_st_TUD(na)%l(la)%j(ja/2)
      ib=nlj_st_TUD(nb)%l(lb)%j(jb/2)
      if (ib<=ia) then
         i_ab=index_ab(ia,ib)
      else
         temp=ii
         ii=jj
         jj=temp
         iphase=-iphase
         i_ab=index_ab(ib,ia)
      endif
      
      nc=n_sp(kk)
      lc=l_sp(kk)
      jc=j2_sp(kk)
      nd=n_sp(ll)
      ld=l_sp(ll)
      jd=j2_sp(ll)
      ic=nlj_st_TUD(nc)%l(lc)%j(jc/2)
      id=nlj_st_TUD(nd)%l(ld)%j(jd/2)
       if (id<=ic) then
         i_cd=index_ab(ic,id)
      else
         temp=kk
         kk=ll
         ll=temp
         iphase=-iphase
         i_cd=index_ab(id,ic)
      endif

      if (i_cd<=i_ab) then
         ja=j2_sp(ii)
         jb=j2_sp(jj)
         jc=j2_sp(kk)
         jd=j2_sp(ll)
      
         m2a = m2_sp(ii)
         m2b = m2_sp(jj)
         m2c = m2_sp(kk)
         m2d = m2_sp(ll)

         itz2a = mt2_sp(ii)
         itz2b = mt2_sp(jj)
         itz2c = mt2_sp(kk)
         itz2d = mt2_sp(ll)

         ii=index_abcd(i_ab,i_cd)
      else
         ja=j2_sp(kk)
         jb=j2_sp(ll)
         jc=j2_sp(ii)
         jd=j2_sp(jj)
      
         m2a = m2_sp(kk)
         m2b = m2_sp(ll)
         m2c = m2_sp(ii)
         m2d = m2_sp(jj)

         itz2a = mt2_sp(kk)
         itz2b = mt2_sp(ll)
         itz2c = mt2_sp(ii)
         itz2d = mt2_sp(jj)
         ii=index_abcd(i_cd,i_ab)
      endif

      clbabt(:,:)=cgt12_tud(:,:,(itz2a+1)/2,(itz2b+1)/2)
      clbcdt(:,:)=cgt12_tud(:,:,(itz2c+1)/2,(itz2d+1)/2)
      
      Rabcd = 0.d0      
      do J12=max(abs(ja-jb),abs(jc-jd))/2,min(ja+jb,jc+jd)/2
         clbab=cgj12(J12,ja/2,(ja+m2a)/2,jb/2,(jb+m2b)/2)
         clbcd=cgj12(J12,jc/2,(jc+m2c)/2,jd/2,(jd+m2d)/2)
         
         Rabcd=Rabcd
     $        +rful(ii)*clbab*clbcd*clbabt(0,0)*clbcdt(0,0)
     $        +rful(ii+1)*clbab*clbcd*clbabt(1,-1)*clbcdt(1,-1)
     $        +rful(ii+2)*clbab*clbcd*clbabt(1,0)*clbcdt(1,0)
     $        +rful(ii+3)*clbab*clbcd*clbabt(1,1)*clbcdt(1,1)         

         ii=ii+4

      end do
      rmo=Rabcd*real(iphase,kind(0.d0))
      end subroutine twobody_op_TUD

      subroutine twobody_op_pn(ii,jj,kk,ll,rmo)
      use spodata
      use spsdata
      use pn_tbme
      use v3b,only: cgj12
      implicit none
      integer,intent(IN) :: ii,jj,kk,ll
      real(kind(0.d0)),intent(OUT) :: rmo
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd,J12,phase_cd
      integer :: ia,ib,ic,id,phase_ab,i_ab,i_cd,index,ip
      integer :: m2a,m2b,m2c,m2d,m2ab
      integer :: itz2a,itz2b,itz2c,itz2d,tz
      real(kind(0.d0)) :: Habcd,clbab,clbcd

      !      na=n_sp(ii)
      la=l_sp(ii)
      ja=j2_sp(ii)
      m2a=m2_sp(ii)
      itz2a=mt2_sp(ii)
      ia=nobt_sp(ii)

!      nb=n_sp(jj)
      lb=l_sp(jj)
      jb=j2_sp(jj)
      m2b=m2_sp(jj)
      itz2b=mt2_sp(jj)
      ib=nobt_sp(jj)

      ip=mod(la+lb,2)
      m2ab=m2a+m2b
      tz=(itz2a+itz2b)/2
!!! No check on parity,T12z conservation
!      nc=n_sp(kk)
!      lc=l_sp(kk)
      jc=j2_sp(kk)
      m2c=m2_sp(kk)
      itz2c=mt2_sp(kk)
      ic=nobt_sp(kk)

!      nd=n_sp(ll)
!      ld=l_sp(ll)
      jd=j2_sp(ll)
      m2d=m2_sp(ll)
!      itz2d=mt2_sp(ll)
      id=nobt_sp(ll)

!      print *,' ii,jj,kk,ll=',ii,jj,kk,ll
!      print *,' ia,ib,ic,id=',ia,ib,ic,id
!      print *,' ip,tz=',ip,tz
!      print *,' ja,jb,jc,jd=',ja,jb,jc,jd

      Habcd=0.d0
      do J12=max(abs(ja-jb),abs(jc-jd),abs(m2ab))/2,min(ja+jb,jc+jd)/2

!         print *,' J12=',J12
         
         if (abs(tz)==1) then
            if ((ia==ib.or.ic==id).and.mod(J12+1,2)==0) cycle
            if (ib<ia) then
               i_ab=pntbdst(J12,ip)%Tz(tz)%sp1(ib)%sp2(ia)
               phase_ab=(-1)**(J12-(ja-jb)/2)
            else
               i_ab=pntbdst(J12,ip)%Tz(tz)%sp1(ia)%sp2(ib)
               phase_ab=1
            endif
            if (id<ic) then
               i_cd=pntbdst(J12,ip)%Tz(tz)%sp1(id)%sp2(ic)
               phase_cd=(-1)**(J12-(jc-jd)/2)
            else
               i_cd=pntbdst(J12,ip)%Tz(tz)%sp1(ic)%sp2(id)
               phase_cd=1
            endif
         else
            if (itz2a==1) then
               i_ab=pntbdst(J12,ip)%Tz(0)%sp1(ia)%sp2(ib)
               phase_ab=1
            else
               i_ab=pntbdst(J12,ip)%Tz(0)%sp1(ib)%sp2(ia)
               phase_ab=(-1)**(J12-(ja-jb)/2)
            endif
            if (itz2c==1) then
               i_cd=pntbdst(J12,ip)%Tz(0)%sp1(ic)%sp2(id)
               phase_cd=1
            else
               i_cd=pntbdst(J12,ip)%Tz(0)%sp1(id)%sp2(ic)
               phase_cd=(-1)**(J12-(jc-jd)/2)
            endif
         endif

!         print *,' i_ab,icd=',i_ab,i_cd
         
         if (i_ab>=i_cd) then
            index=i2belnpoi_pn(J12,ip,tz)+i_cd+i_ab*(i_ab-1)/2
         else
            index=i2belnpoi_pn(J12,ip,tz)+i_ab+i_cd*(i_cd-1)/2
         endif

!         print *,' index=',index
         
         clbab=cgj12(J12,ja/2,(ja+m2a)/2,jb/2,(jb+m2b)/2)
         clbcd=cgj12(J12,jc/2,(jc+m2c)/2,jd/2,(jd+m2d)/2)

!         print *,' clbab=',clbab
!         print *,' clbcd=',clbcd
         
         Habcd=Habcd+clbab*clbcd*rful_pn(tz)%el(index) 
     $        *real(phase_ab*phase_cd,kind(0.d0))

!         print *,' rful_pn(tz)%el(index)=',rful_pn(tz)%el(index)
!         print *,' Habcd=',Habcd
         
      end do
      rmo=Habcd
      end subroutine twobody_op_pn

      double precision FUNCTION CLEBRD(A,B,C,D,E,F)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
      IMPLICIT double precision (A-H,O-Z)
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
      LOGICAL FIRST
*      REAL A,B,C,D,E,F
C      REAL*8 FACLOG
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=idINT(2.d0*A)
      IB=idINT(2.d0*C)
      IC=idINT(2.d0*E)
      ID=idINT(2.d0*B)
      IE=idINT(2.d0*D)
      IF=idINT(2.d0*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBID(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEBD(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL THFACD
      RAC=0.0D0
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1.0D0
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1.0D0)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1.D0)**((IAM-IDM)/2)*DSQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1.D0)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5D0*(DLOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0.D0
      S1=(-1.D0)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*DEXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEBRD=RAC
      RETURN
      END

      SUBROUTINE THFACD
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      double precision FACLOG,FN
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
c*      DATA FIRST/.TRUE./,LFACT/LFACTC/
      FIRST=.FALSE.
      LFACT=200
      FACLOG(1)=0.D0
      FACLOG(2)=0.D0
      FN=1.D0
      DO 10 I=3,LFACT
      FN=FN+1.D0
      FACLOG(I)=FACLOG(I-1)+DLOG(FN)
   10 CONTINUE
      RETURN
      END

      block data
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      double precision FACLOG
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
      DATA FIRST/.TRUE./,LFACT/LFACTC/
      end

      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      real d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      real u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,epslon
      integer ind(m)
c
c     ----- This routine is a translation of the algol procedure bisect,
c     ----- Num. Math. 9, 386-393(1967) by Barth, Martin, and Wilkinson.
c     ----- Handbook for Auto. Comp., Vol.ii - Linear Algebra, 249-256(1971).
c
c     ----- This routine finds those eigenvalues of a tridiagonal
c     ----- symmetric matrix between specified boundary indices,
c     ----- using bisection.
c
c     ----- on input
c
c     -----    n is the order of the matrix.
c
c     -----    eps1 is an absolute error tolerance for the computed
c     -----      eigenvalues.  if the input eps1 is non-positive,
c     -----      it is reset for each submatrix to a default value,
c     -----      namely, minus the product of the relative machine
c     -----      precision and the 1-norm of the submatrix.
c
c     -----    d contains the diagonal elements of the input matrix.
c
c     -----    e contains the subdiagonal elements of the input matrix
c     -----      in its last n-1 positions.  e(1) is arbitrary.
c
c     -----    e2 contains the squares of the corresponding elements of e.
c     -----      e2(1) is arbitrary.
c
c     -----    m11 specifies the lower boundary index for the desired
c     -----      eigenvalues.
c
c     -----    m specifies the number of eigenvalues desired.  the upper
c     -----      boundary index m22 is then obtained as m22=m11+m-1.
c
c     ----- on output
c
c     -----    eps1 is unaltered unless it has been reset to its
c     -----      (last) default value.
c
c     -----    d and e are unaltered.
c
c     -----    elements of e2, corresponding to elements of e regarded
c     -----      as negligible, have been replaced by zero causing the
c     -----      matrix to split into a direct sum of submatrices.
c     -----      e2(1) is also set to zero.
c
c     -----    lb and ub define an interval containing exactly the desired
c     -----      eigenvalues.
c
c     -----    w contains, in its first m positions, the eigenvalues
c     -----      between indices m11 and m22 in ascending order.
c
c     -----    ind contains in its first m positions the submatrix indices
c     -----      associated with the corresponding eigenvalues in w --
c     -----      1 for eigenvalues belonging to the first submatrix from
c     -----      the top, 2 for those belonging to the second submatrix, etc..
c
c     -----    ierr is set to
c     -----      zero       for normal return,
c     -----      3*n+1      if multiple eigenvalues at index m11 make
c     -----                 unique selection impossible,
c     -----      3*n+2      if multiple eigenvalues at index m22 make
c     -----                 unique selection impossible.
c
c     -----    rv4 and rv5 are temporary storage arrays.
c
c     ----- note that routine tql1, imtql1, or tqlrat is generally faster
c     ----- than tridib, if more than n/4 eigenvalues are to be found.
c
c     ----- Questions and comments should be directed to Burton S. Garbow,
c     ----- Mathematics and Computer Science Div, Argonne National Laboratory
c
c     ----- This version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0e0
c     .......... look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0e0
         if (i .ne. n) u = abs(e(i+1))
         xu = amin1(d(i)-(x1+u),xu)
         x0 = amax1(d(i)+(x1+u),x0)
         if (i .eq. 1) goto 20
         tst1 = abs(d(i)) + abs(d(i-1))
         tst2 = tst1 + abs(e(i))
         if (tst2 .gt. tst1) goto 40
   20    e2(i) = 0.0e0
   40 continue
c
      x1 = n
      x1 = x1 * epslon(amax1(abs(xu),abs(x0)))
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     .......... determine an interval containing exactly
c                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) goto 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5e0
      if (x1 .eq. v) goto 980
      goto 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      goto 50
   70 x0 = x1
      goto 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) goto 90
      x0 = t2
      isturm = 2
      goto 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) goto 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0e0
c
      do 120 q = p, n
         x1 = u
         u = 0.0e0
         v = 0.0e0
         if (q .eq. n) goto 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = amin1(d(q)-(x1+u),xu)
         x0 = amax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0e0) goto 140
  120 continue
c
  140 x1 = epslon(amax1(abs(xu),abs(x0)))
      if (eps1 .le. 0.0e0) eps1 = -x1
      if (p .ne. q) goto 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) goto 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      goto 900
  180 x1 = x1 * (q - p + 1)
      lb = amax1(t1,xu-x1)
      ub = amin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      goto 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      goto 320
  220 m2 = s
      if (m1 .gt. m2) goto 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) goto 260
            xu = rv4(i)
            goto 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5e0
         if ((x0 - xu) .le. abs(eps1)) goto 420
         tst1 = 2.0e0 * (abs(xu) + abs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) goto 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0e0
c
         do 340 i = p, q
            if (u .ne. 0.0e0) goto 325
            v = abs(e(i)) / epslon(1.0e0)
            if (e2(i) .eq. 0.0e0) v = 0.0e0
            goto 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0e0) s = s + 1
  340    continue
c
         goto (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) goto 400
         xu = x1
         if (s .ge. m1) goto 380
         rv4(m1) = x1
         goto 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         goto 300
  400    x0 = x1
         goto 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) goto 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) goto 910
         if (k .gt. m2) goto 940
         if (rv5(k) .ge. w(l)) goto 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         goto 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) goto 100
      goto 1001
c     .......... set error -- interval cannot be found containing
c                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end


      Subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      real d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      real u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,
     x       pythag
      integer ind(m)
c
c     ----- This routine is a translation of the inverse iteration technique 
c     ----- in the algol procedure Tristurm by Peters and Wilkinson.
c     ----- Handbook for Auto. Comp., Vol.ii - Linear Algebra, 418-439(1971).
c
c     ----- This routine finds those eigenvectors of a tridiagonal
c     ----- symmetric matrix corresponding to specified eigenvalues,
c     ----- using inverse iteration.
c
c     ----- on input
c
c     -----    nm must be set to the row dimension of two-dimensional
c     -----      array parameters as declared in the calling program
c     -----      dimension statement.
c
c     -----    n is the order of the matrix.
c
c     -----    d contains the diagonal elements of the input matrix.
c
c     -----    e contains the subdiagonal elements of the input matrix
c     -----      in its last n-1 positions.  e(1) is arbitrary.
c
c     -----    e2 contains the squares of the corresponding elements of e,
c     -----      with zeros corresponding to negligible elements of e.
c     -----      e(i) is considered negligible if it is not larger than
c     -----      the product of the relative machine precision and the sum
c     -----      of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c     -----      0.0e0 if the eigenvalues are in ascending order, or 2.0e0
c     -----      if the eigenvalues are in descending order.  if  bisect,
c     -----      tridib, or  imtqlv  has been used to find the eigenvalues,
c     -----      their output e2 array is exactly what is expected here.
c
c     -----    m is the number of specified eigenvalues.
c
c     -----    w contains the m eigenvalues in ascending or descending order.
c
c     -----    ind contains in its first m positions the submatrix indices
c     -----      associated with the corresponding eigenvalues in w --
c     -----      1 for eigenvalues belonging to the first submatrix from
c     -----      the top, 2 for those belonging to the second submatrix, etc.
c
c     ----- on output
c
c     -----    all input arrays are unaltered.
c
c     -----    z contains the associated set of orthonormal eigenvectors.
c     -----      any vector which fails to converge is set to zero.
c
c     -----    ierr is set to
c     -----      zero       for normal return,
c     -----      -r         if the eigenvector corresponding to the r-th
c     -----                 eigenvalue fails to converge in 5 iterations.
c
c     -----    rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     ----- calls pythag for  sqrt(a*a + b*b) .
c
c     ----- Questions and comments should be directed to Burton S. Garbow,
c     ----- Mathematics and Computer Science Div, Argonne National Laboratory
c
c     ----- this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (m .eq. 0) goto 1001
      tag = 0
      order = 1.0e0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
 100  p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) goto 140
         if (e2(q+1) .eq. 0.0e0) goto 140
 120  continue
c     .......... find vectors by inverse iteration ..........
 140  tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) goto 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) goto 510
c     .......... check for isolated root ..........
         xu = 1.0e0
         if (p .ne. q) goto 490
         rv6(p) = 1.0e0
         goto 870
  490    norm = abs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = amax1(norm, abs(d(i))+abs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0e-3 * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / sqrt(uk)
         s = p
 505     group = 0
         goto 520
c     .......... look for close or coincident roots ..........
 510     if (abs(x1-x0) .ge. eps2) goto 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0e0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
 520     v = 0.0e0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) goto 560
            if (abs(e(i)) .lt. abs(u)) goto 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0e0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            goto 580
 540        xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0e0
 560        u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
 580     continue
c
         if (u .eq. 0.0e0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0e0
         rv3(q) = 0.0e0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
 600     do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
 620     continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) goto 700
         j = r
c
         do 680 jj = 1, group
 630        j = j - 1
            if (ind(j) .ne. tag) goto 630
            xu = 0.0e0
c
            do 640 i = p, q
 640        xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
 660        rv6(i) = rv6(i) - xu * z(i,j)
c
 680     continue
c
 700     norm = 0.0e0
c
         do 720 i = p, q
 720     norm = norm + abs(rv6(i))
c
         if (norm .ge. 1.0e0) goto 840
c     .......... forward substitution ..........
         if (its .eq. 5) goto 830
         if (norm .ne. 0.0e0) goto 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         goto 780
 740     xu = eps4 / norm
c
         do 760 i = p, q
 760     rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector iterate ..........
 780     do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) goto 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
 800        rv6(i) = u - rv4(i) * rv6(i-1)
 820     continue
c
         its = its + 1
         goto 600
c     .......... set error -- non-converged eigenvector ..........
 830     ierr = -r
         xu = 0.0e0
         goto 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
 840     u = 0.0e0
c
         do 860 i = p, q
 860     u = pythag(u,rv6(i))
c
         xu = 1.0e0 / u
c
 870     do 880 i = 1, n
 880     z(i,r) = 0.0e0
c
         do 900 i = p, q
 900     z(i,r) = rv6(i) * xu
c
         x0 = x1
 920  continue
c
      if (q .lt. n) goto 100
 1001 return
      end


      real function pythag(a,b)
      real a,b
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
      real p,r,s,t,u
      p = amax1(abs(a),abs(b))
      if (p .eq. 0.0e0) goto 20
      r = (amin1(abs(a),abs(b))/p)**2
 10   continue
      t = 4.0e0 + r
      if (t .eq. 4.0e0) goto 20
      s = r/t
      u = 1.0e0 + 2.0e0*s
      p = u*p
      r = (s/u)**2 * r
      goto 10
 20   pythag = p
      return
      end


      real function epslon (x)
      real x
      real a,b,c,eps
c     ----- estimate unit roundoff in quantities of size x.
c     ----- this program should perform properly on all systems
c     ----- satisfying the following two assumptions,
c     -----    1.  the base used in representing floating point
c     -----        numbers is not a power of three.
c     -----    2.  the quantity  a  in statement 10 is represented to
c     -----        the accuracy used in floating point variables
c     -----        that are stored in memory.
c     ----- the statement number 10 and the go to 10 are intended to
c     ----- force optimizing compilers to generate code satisfying
c     ----- assumption 2.
c     ----- under these assumptions, it should be true that,
c     -----        a  is not exactly equal to four-thirds,
c     -----        b  has a zero for its last bit or digit,
c     -----        c  is not exactly equal to one,
c     -----        eps  measures the separation of 1.0 from
c     -----             the next larger floating point number.
c     ----- the developers of eispack would appreciate being informed
c     ----- about any systems where these assumptions do not hold.
c
c     ----- this version dated 4/6/83.
c
      a = 4.0e0/3.0e0
 10   b = a - 1.0e0
      c = b + b + b
      eps = abs(c-1.0e0)
      if (eps .eq. 0.0e0) goto 10
      epslon = eps*abs(x)
      return
      end


      subroutine clock(cputime,walltime,iunit,iprint)
      real(8),intent(OUT) :: cputime,walltime
      integer,intent(IN) :: iunit,iprint
      real(8) fcputime,fwalltime,cputimex
      data fcputime/0.d0/,fwalltime/0.d0/
      character date *8,timechar *10
c*      integer count,fcount/0/
      integer days,hours,minutes,seconds,milisec

      call cpu_time(cputimex)
      cputime=cputimex
ccc      cputime=0.d0
      call date_and_time(date,timechar)
c*      call system_clock(count)

      days=10*(iachar(date(7:7))-48)+(iachar(date(8:8))-48)
c*      print *,' days=',days
      hours=10*(iachar(timechar(1:1))-48)+(iachar(timechar(2:2))-48)
c*      print *,' hours=',hours
      minutes=10*(iachar(timechar(3:3))-48)+(iachar(timechar(4:4))-48)
c*      print *,' minutes=',minutes
      seconds=10*(iachar(timechar(5:5))-48)+(iachar(timechar(6:6))-48)
c*      print *,' seconds=',seconds
      milisec=100*(iachar(timechar(8:8))-48)
     +     +10*(iachar(timechar(9:9))-48)
     +     +(iachar(timechar(10:10))-48)
c*      print *,' milisec=',milisec

      walltime=dble(3600*(24*days+hours)
     +     +60*minutes+seconds)
     +     +dble(milisec)/1000.d0

      if (iprint==1) then
         write(iunit,*)
         write(iunit,1000) date(5:6),date(7:8),date(1:4),timechar(1:2),
     +        timechar(3:4),timechar(5:6),timechar(8:10)
 1000    format('   *** Date (m-d-y): ',(A2),' ',(A2),' ',(A4),
     +     /,'   *** Time: ',(A2),' h ',(A2),' min ',(A2),'.',(A3),' s')
         if (cputime>0.000d0) then
            write(iunit,2000) cputime
 2000       format('   *** CPU time= ',f14.3,' s')
         endif   
c*      write(iunit,3000) count
c* 3000 format('   *** System count:',i20)
         if (fcputime/=0.d0) then
            if (cputime>fcputime) then
               write(iunit,4000) cputime-fcputime
 4000    format('   *** CPU time for this calculation:',f14.3,' s')
            endif   
         else
            fcputime=cputime
         endif   
         if (fwalltime/=0.d0) then
            if (walltime>fwalltime) then
               write(iunit,5000) walltime-fwalltime
 5000          format('   *** Wall-clock time for this calculation:',
     +              f14.3,' s')
            endif
         else
            fwalltime=walltime
         endif   
c*         if (fcount/=0) then
c*            write(iunit,6000) count-fcount
c* 6000       format(' System count for this calculation:',i20)
c*         else
c*            fcount=count
c*         endif   
         write(iunit,*)
      endif
      
      end

c#############################################################################
c#############################################################################
c#############################################################################
c**      integer function popcntt(x)
c**      implicit integer (a-z)
c     population count of a 32-bit integer
c     ====================================
c**      data m1/z'55555555'/      ! binary 01010101...
c**     &     m2/z'33333333'/      ! binary 00110011...
c**     &     m4/z'0f0f0f0f'/      ! binary 00001111...
c     (if the word size exceeds 32 bits, the 2**32 and higher bits
c     are presumed clear)
c     sum the bits within bytes
c     -------------------------
c     IAND = logical AND function
c     NOT(z) = ones complement of z
c     ISHFT(z,-k) = z shifted right by k bit positions (for k > 0)
c**      y = IAND(x,m1) + ISHFT(IAND(x,NOT(m1)),-1)
c**      y = IAND(y,m2) + ISHFT(IAND(y,NOT(m2)),-2)
c**      y = IAND(y,m4) + ISHFT(IAND(y,NOT(m4)),-4)
c     now sum the counts in the bytes
c     -------------------------------
c**      y = y + ISHFT(y, -8)
c**      y = y + ISHFT(y,-16)
c**      popcntt = IAND(y,255)
c**      return
c**      end


c---------------------------------------------------------------------
c     *******  The following subroutines and functions are     *******
c     *******  to be used in the sequential environment.       *******
c---------------------------------------------------------------------


c      subroutine MPI_INIT(ierr)
c      return
c      end

c      subroutine MPI_COMM_RANK(icomm,iproc,ierr)
c      iproc=0
c      return
c      end

c      subroutine MPI_COMM_SIZE(icomm,nproc,ierr)
c      nproc=1
c      return
c      end

c      subroutine MPI_Barrier(icomm,ierr)
c      return
c      end

c      subroutine MPI_Finalize(ierr)
c      return
c      end

c      subroutine MPI_Bcast()
c      return
c      end

c      subroutine MPI_Reduce()
c      return
c      end

c      subroutine MPI_Allreduce()
c      return
c      end

c      subroutine MPI_Abort()
c      return
c      end

c      subroutine MPI_Send()
c      return
c      end

c      subroutine MPI_Recv()
c      return
c      end

