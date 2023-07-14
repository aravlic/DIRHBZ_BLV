c======================================================================c

      subroutine start_vapor(lpr)
c======================================================================c
c
c     initializes potentials (inin=1) and pairing tensor (inink=1)
c     initialization for the vapor phase
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
        
        logical lpr
        common /initia/ inin,inink
        
        if(inin.eq.1) call startpot_vapor(lpr)
        if(inink.eq.1) call startdel_vapor(lpr)
        
        return
        end
        
c======================================================================c
        subroutine startpot_vapor(lpr)
c======================================================================c
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
c
        logical lpr
c
        character nucnam*2                                 ! common nucnuc
c
        dimension vso(2),r0v(2),av(2),rso(2),aso(2)
        dimension rrv(2),rls(2),vp(2),vls(2)
c
        common /baspar/ hom,hb0,b0
        common /coulmbv/ cou_v(MG),drvp_v(MG)
        common /defbas/ beta0,q,bp,bz
        common /defini/ betai
        common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
        common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
        common /gfvsq / sq(0:IGFV)
        common /initia/ inin,inink
        common /mathco/ zero,one,two,half,third,pi
        common /nucnuc/ amas,npr(3),nucnam
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /physco/ hbc,alphi,r0
        common /potpotv/ vps_v(MG,2),vms_v(MG,2)
        common /deldelv/ de_v(NHHX,NB2X)
        common /pair  / del(2),spk(2),spk0(2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c=======================================================================
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
        data v0/-71.28/,akv/0.4616/
        data r0v/1.2334,1.2496/,av/0.615,0.6124/
        data vso/11.1175,8.9698/
        data rso/1.1443,1.1401/,aso/0.6476,0.6469/

c
        if (lpr) then
        write(l6,*) ' ***** BEGIN START VAPOR*************************'
        endif

c
        betas = betai * half* dsqrt(5.d0/(4.d0*pi))
        fac =  one + betas
        fac = (one + betas*cos(120.d0*pi/180.d0))*fac
        fac = (one + betas*cos(-120.d0*pi/180.d0))*fac
        fac = fac**(-third)
        do ih = 0,NGH
            z  = zb(ih)
            zz = z**2
        do il = 0,NGL
            ihl = 1+ih + il*(NGH+1)
            rr = (rb(il)**2 + zz)
            r  = sqrt(rr)
c
C------- Woods Saxon
            ctet = zz/rr
            p20  = 3*ctet - one
            facb = fac*(one + betas*p20) 
            do it = 1,2
            ita = 3-it
            rrv(it) = r0v(it)*amas**third
            rls(it) = rso(it)*amas**third
            vp(it)  = v0*(one - akv*(npr(it)-npr(ita))/amas)
            vls(it) = vp(it) * vso(it)
c
            argv=(r - rrv(it)*facb) / av(it)
            if (argv.lt.65.d0) then
                u = vp(it) /(one + exp(argv))
            else
                u = 0.d0
            endif
            argo=(r - rls(it)*facb) / aso(it)
            if (argo.lt.65.d0) then
                w  = -vls(it)/(one + exp(argo ))
            else
                w = 0.d0
            endif
c
            vps_v(ihl,it) = 0.d0
            vms_v(ihl,it) = 0.d0
c
            enddo   ! it
c
c------- Coulomb
            cou_v(ihl) = zero
            if (icou.ne.0) then
            rc = rrv(2)
            if (r.lt.rc) then
                c = half*(3/rc-r*r/(rc**3))
            else
                c = one/r
            endif
            cou_v(ihl)   = c*npr(2)/alphi
            vps_v(ihl,2) = vps_v(ihl,2) + cou_v(ihl)*hbc
            vms_v(ihl,2) = vms_v(ihl,2) + cou_v(ihl)*hbc
            endif
        enddo   ! il
        enddo   ! ih
        


  100 format(a,2f10.4)
c
        if (lpr) then
        write(l6,*) ' ***** END START VAPOR***************************'
        endif
c
        return
c-end STARTPOT
        end 
c======================================================================c
        subroutine startdel_vapor(lpr)
c======================================================================c
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
c
        logical lpr
        character tb*6  

        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2) 
        common /mathco/ zero,one,two,half,third,pi
        common /deldelv/ de_v(NHHX,NB2X)
        common /pair  / del(2),spk(2),spk0(2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp    
        
c---- set vapor pairing field to zero
        call mzero(NHHX,NHHX,NB2X,de_v)

c---- CHECK if spk needs to be constrained ?
        
        return
c-end STARTDEL
        end 


C======================================================================c

      subroutine gamma_vapor(it)

c======================================================================c
c
c     calculates the Dirac-Matrix in the Hartee-equation - VAPOR 
c 
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'

        character tb*6                                            ! blokap      
c
        common /baspar/ hom,hb0,b0
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /gamgamv/ hh_v(NHHX,NB2X) ! NOTE
        common /masses/ amu,ames(4)
        common /mathco/ zero,one,two,half,third,pi
        common /physco/ hbc,alphi,r0
        common /potpotv/ vps_v(MG,2),vms_v(MG,2)
        common /single/ sp(NFGX,NBX)
c
        emcc2 = 2*amu*hbc
        f = hbc/b0
c
c      write(*,*) 'Mean-field matrix elements:'
        do ib = 1,nb
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            i0f = ia(ib,1)
            i0g = ia(ib,2)
            m   = ib + (it-1)*NBX
c 
            do n2 = 1,nf
            do n1 = 1,ng
                hh_v(nf+n1+(n2-1)*nh,m) = f*sp(n1+(n2-1)*ng,ib)
            enddo
            enddo
            call pot(i0f,nh,nf,vps_v(1,it),hh_v(1,m))
            call pot(i0g,nh,ng,vms_v(1,it),hh_v(nf+1+nf*nh,m))
            do n = nf+1,nh
            hh_v(n+(n-1)*nh,m) = hh_v(n+(n-1)*nh,m) - emcc2
            enddo
c
c     symmetrize HH
            do n2 =    1,nh
            do n1 = n2+1,nh
            hh_v(n2+(n1-1)*nh,m) = hh_v(n1+(n2-1)*nh,m)
            enddo
            enddo

c--- print matrix - Ravlic
c      do n1 = 1,nh
c        do n2 = 1,nh
c           write(*,*) m,n1, n2, hh(n1+(n2-1)*nh,m)
c        enddo
c      enddo
    
        enddo  !ib


        return
C-end-GAMMA
        end


c======================================================================c

      subroutine dirhb_vapor(it,lpr)

c======================================================================c
c
c     solves the RHB-Equation 
c     IT    = 1 for neutrons
c     IT    = 2 for protons
c
c     Implements method of A. Bjelcic et al. []
c
c     Diagonalization is split into 3 steps:
c
c     1) diagonalize h, get Z and calculate dt
c     2) make H and diagonalize to get u,v
c     3) transform u,v to U, V
c
c     Solves the VAPOR system
c     Note that chemical potential is not constrained
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical lpr,lprl
c
        character*1 bbb
        character*8 tbb(NHBX)
        character tp*1,tis*1,tit*8,tl*1                           ! textex
        character tb*6                                            ! blokap
        character tt*11                                            ! quaosc
        character nucnam*2                                        ! nucnuc

c---- Ravlic, LAPACK diagonalization
        integer info, lwork
        external dsyev
        dimension W_HH_1(NHX)
        dimension W_HH_2(2*NFX)
        integer lwmax
        parameter ( lwmax = 10000000 )
        dimension HH_work(lwmax)
        dimension work(lwmax)
        dimension W(2*NFX)

c----------------
        dimension hb(NHBQX),e(NHBX),ez(NHBX)
c---- Ravlic, BCS
        dimension elsp(KX), dksp(KX)
c---- Ravlic      
        dimension hh_hh(NHHX), zz_hh(NHHX)
        dimension HH_1(NHX,NHX,NBX)
        dimension HH_2(2*NFX,2*NFX,NBX)
        dimension e_work(NHX), e_hh(NHX), ee_hh(NFFX,NB2X)
        dimension ez_hh(NHX)
        dimension de_hh(NFFX,NB2X)
        dimension hb_hh(4*NFFX)
        dimension fguv_temp(2*NFX,KX,4)
        dimension de_new(NFFX,NB2X)
        dimension Delta(NFX,NFX,NB2X)
        dimension tempmat(NFX,NFX)
        dimension Zh(NFX,NFX, NB2X)

c
        common /blodir/ ka(NBX,4),kd(NBX,4) !CHECK
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2) ! CHECK
        common /deldelv/ de_v(NHHX,NB2X)
        common /fermi / ala(2),tz(2)
        common /gamgamv/ hh_v(NHHX,NB2X)
        common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
        common /mathco/ zero,one,two,half,third,pi
        common /nucnuc/ amas,npr(3),nucnam
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /quaosc/ nt,nnn(NTX,5),tt(NTX)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)
        common /waveuvv/ fguv_v(NHBX,KX,4),equ_v(KX,4)
        common /temp/ temp
c
        data maxl/200/,epsl/1.d-8/,bbb/'-'/,lprl/.false./
        data fm10/1.0d-10/
    
c
        if (.true.) then
        write(l6,*) ' ****** BEGIN DIRHB VAPOR***********************'
        endif

            al    = ala(it)  
c
            sn  = zero
            klp = 0
            kla = 0
c======================================================================c
        do ib = 1,nb            ! loop over differnt blocks
c======================================================================c
            mul  = mb(ib)
            nf   = id(ib,1)
            ng   = id(ib,2)
            nh   = nf + ng
            nhb  = nh + nh
            m    = ib + (it-1)*NBX

c======================================================================c
c        construct hh_hh matrix, (nf+ng) x (nf+ng) 
c======================================================================c
            do n2 = 1,nh
            do n1 = n2,nh
                HH_1(n1,n2,ib) = hh_v(n1+(n2-1)*nh,m)
            enddo !n2
            enddo !n1



c======================================================================c
c        diagonalization: e_work has nh eigenvalues
c                         hh_hh has nh x nh eigenvectors
c======================================================================c
            
c      call sdiag(nh,nh,hh_hh,e_work,hh_hh,ez_hh,+1)
c--- try LAPACK function
c    Query the optimal workspace.
*
        lwork = -1
        CALL DSYEV( 'V', 'L', nh, HH_1(1,1,ib), NHX
     &     , W_HH_1, work, lwork, info)
        lwork = MIN( lwmax, INT( work( 1 ) ) )
c      write(*,*) 'Dimension of problem: ', lwork, lwmax, info
c      read*

c     Solve eigenproblem.

        CALL DSYEV( 'V', 'L', nh, HH_1(1,1,ib), NHX
     &   , W_HH_1, work, lwork, info )

c     Check for convergence.

        IF( info.GT.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF

c======================================================================c
c        store eigenvalues and eigenvectors
c        
c        Notice: anti-particles decouple at this point
c        e_hh - matrix containing nf eigenvalues corresponding
c        to particles
c        zz_hh are the eigenvectors
c======================================================================c

            do k = 1,nf
            e_hh(k) = W_HH_1(ng+k) !e_work(ng+k)
            do n = 1,nh
                zz_hh(n + (k-1)*nf) = HH_1(n,ng+k,ib)!hh_hh(n+(ng+k-1)*nh)
            enddo
            enddo !k
            
            do k = 1,nf
            do n = 1,nf
                Zh(n,k,m) = HH_1(n,ng+k,ib) !hh_hh(n+(ng+k-1)*nh)
            enddo
            enddo !k

c======================================================================c
c        store anti-particle contribution
c======================================================================c

            ka(ib,it+2) = kla
            do k = 1,ng
            kla = kla + 1
            equ_v(kla,it+2) = W_HH_1(k) !e_work(k) 
            do n = 1,nh
                fguv_v(n,kla,it+2) = HH_1(n,k,ib) !hh_hh(n+(k-1)*nh)
            enddo
            enddo
            kd(ib,it+2) = kla - ka(ib,it+2)


c======================================================================c
c        calculate \tilde{de} = Z^T de Z
c        \tilde{de} has dimenson nf x nf
c======================================================================c


            do j = 1 , nf
                do i = j , nf
                    Delta(i,j,m) = de_v( i + (j-1)*nh , m )
                enddo
            enddo

c--- using LAPACK functions
            call dsymm('L','L',  nf,nf,
     &              +1.0d+0,       Delta(1,1,m),NFX,
     &                             Zh(1,1,m),NFX,
     &              +0.0d+0,        tempmat(1,1),NFX  );
            call dgemm('T','N',  nf,nf,nf,
     &              +1.0d+0,       Zh(1,1,m),NFX,
     &                              tempmat(1,1),NFX,
     &              +0.0d+0,       Delta(1,1,m),NFX  );

            do n2 = 1, nf
            do n1 = 1 ,nf
                de_hh(n1 + (n2-1)*nf,m) = Delta(n1,n2,m)
            enddo
            enddo
c======================================================================c
c        construct diagonal matrix ee_hh(nf x nf)
c======================================================================c
            do n1 = 1,nf
            do n2 = 1,nf
                if (n1.eq.n2) then
                ee_hh(n1 + (n2-1)*nf,m) = e_hh(n1)
                endif
            enddo !n2
            enddo !n1
c======================================================================c
c        construct H_\lambda 2nf x 2nf
c======================================================================c


c------- construct the new matrix (note that it is symmetric)
c------- also this saves the LOWER triangle
            do n2 = 1,nf
            do n1 = n2,nf
                HH_2(n1,n2,ib) = ee_hh(n1+(n2-1)*nf,m)
                HH_2(n1+nf,n2+nf,ib) = -ee_hh(n1+(n2-1)*nf,m)
                HH_2(n1+nf,n2,ib) = de_hh(n1+(n2-1)*nf,m)
                HH_2(n2+nf,n1,ib) = de_hh(n2+(n1-1)*nf,m)
            enddo !n2
            HH_2(n2,n2,ib) = HH_2(n2,n2,ib) - al
            HH_2(n2+nf,n2+nf,ib) = -HH_2(n2,n2,ib)
            enddo !n1

c======================================================================c
c        diagonalize H_\lambda 2nf x 2nf
c======================================================================c
c      write(*,*) 2*nf, 2*NFX, NHBX
c      call sdiag(2*nf,2*nf,hb_hh,e,hb_hh,ez,+1)
        lwork = -1
        CALL DSYEV( 'V', 'L', 2*nf, HH_2(1,1,ib), 2*NFX
     &     , W, work, lwork, info)
        lwork = MIN( lwmax, INT( work( 1 ) ) )
c      write(*,*) 'Dimension of problem 2: ', lwork, lwmax, info
c      read*

c     Solve eigenproblem.

        CALL DSYEV( 'V', 'L', 2*nf, HH_2(1,1,ib), 2*NFX
     &   , W, work, lwork, info )

c     Check for convergence.

        IF( info.GT.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF


c======================================================================c
c        store eigenvalues and wave functions
c        particles, u and v
c======================================================================c

            ka(ib,it) = klp
            do k = 1,nf
            klp = klp + 1
            equ_v(klp,it) = W(nf+k) !e(nf+k)
            do n = 1,2*nf
                fguv_temp(n,klp,it) = HH_2(n,nf+k,ib) !hb_hh(n+(nf+k-1)*2*nf)
            enddo
            enddo
            kd(ib,it) = klp - ka(ib,it)


c======================================================================c
c        transformation from u,v --> U,V
c        U = Z u, V = Z v,  where Z is hh_hh
c======================================================================c

c--- U matrix
            do i = 1,nf
                iik = i + ka(ib,it)
                do k = 1,nh
                fguv_v(k,iik,it) = zero
                do n = 1,nf
                    fguv_v(k,iik,it) = fguv_v(k,iik,it)  
     &          +  HH_1(k,ng+n,ib)*fguv_temp(n,iik,it)

                enddo !n
            enddo !k
        enddo !i
c--- V matrix
        do i = 1,nf
            iik = i + ka(ib,it)
            do k = 1,nh
                fguv_v(nh+k,iik,it) = zero
                do n = 1,nf
                    fguv_v(nh+k,iik,it) = fguv_v(nh+k,iik,it)  
     &          +  HH_1(k,ng+n,ib)*fguv_temp(nf+n,iik,it)

                enddo !n
            enddo !k 
        enddo !i

c


c======================================================================c
        enddo   ! ib
c======================================================================c 

c======================================================================c
c        check transformation
c======================================================================c 
            sn = zero
            sn2 = zero
            un = zero
            un2 = zero
            do ib = 1,nb  
            nf   = id(ib,1)
            ng   = id(ib,2)
            nh = nf + ng
            do k = 1,nf
                do n = 1,nf
                sn = sn + fguv_temp(nf+n,ka(ib,it)+k,it)**2
                un = un + fguv_temp(n,ka(ib,it)+k,it)**2
            enddo
            enddo
            do k = 1,nf
                do n = 1,nh
                sn2 = sn2 + fguv_v(nh+n,ka(ib,it)+k,it)**2
                un2 = un2 + fguv_v(n,ka(ib,it)+k,it)**2
            enddo
            enddo
    
            enddo !ib
         write(l6,*) 'Particle number check: ', sn, sn2
         write(l6,*) 'U norm check:', un, un2
c======================================================================c
c        check normalization, UU + VV = 1
c======================================================================c 
            do ib = 1,nb  
            nf   = id(ib,1)
            ng   = id(ib,2)
            nh = nf + ng
            do k = 1,nf
                kk = ka(ib,it) + k
                snorm = zero
                do n1 = 1,nf
                snorm = snorm + fguv_temp(n1,kk,it)*fguv_temp(n1,kk,it) 
     & + fguv_temp(nf+n1,kk,it)*fguv_temp(nf+n1,kk,it)
            enddo !n1

c            write(*,*) kk, snorm
            if (abs(snorm-1.d0).gt.1e-8) then
                    write(*,*) 'Wrong NORM 1 !', ib,it,kk, snorm
                    stop
            endif

            enddo !k
            enddo !ib
            
        do ib = 1,nb  
            nf   = id(ib,1)
            ng   = id(ib,2)
            nh = nf + ng
            do k = 1,nf
                kk = ka(ib,it) + k
                snorm = zero
                do n1 = 1,nh
                snorm = snorm + fguv_v(n1,kk,it)*fguv_v(n1,kk,it) 
     & + fguv_v(nh+n1,kk,it)*fguv_v(nh+n1,kk,it)
            enddo !n1

c            write(*,*) kk, snorm
            if (abs(snorm-1.d0).gt.1e-8) then
                    write(*,*) 'Wrong NORM 2 !',ib,it, kk, snorm
                    stop
            endif

            enddo !k
            enddo !ib


        if (.true.) then
        write(l6,*) ' ****** END DIRHB VAPOR*************************'
        endif
c
  101 format(i4,a,i4,3f13.8)
  113 format(i4,a,a1,3f13.8)
c
c      read*
        return
C-end-DIRHB
        end


c=====================================================================c

      subroutine denssh_vapor(it,lpr)

c=====================================================================c
C
c     calculates the densities in oscillator basis
c     for VAPOR solution 
C
c---------1---------2---------3---------4---------5---------6---------7-
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
c
        logical lpr
c
        character tb*6                                         ! blokap
        character tt*11                                        ! quaosc
        character tp*1,tis*1,tit*8,tl*1                        ! textex
c
        common /blodir/ ka(NBX,4),kd(NBX,4)      
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /pair  / del(2),spk(2),spk0(2)
        common /pairv / spk_v(2)
      common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
        common /rokaosv/ rosh_v(NHHX,NB2X),aka_v(MVX,2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)
        common /waveuvv/ fguv_v(NHBX,KX,4),equ_v(KX,4)     
        common /temp/ temp
c 
        if (lpr) then
        write(l6,*) ' ****** BEGIN DENSSH VAPOR**********************'
        endif
c
c
        sp  = zero
        il = 0
c       write(*,*) 'Single-particle density:'
c=====================================================================c
        do ib = 1,nb            ! loop over the blocks
c=====================================================================c
c         write(*,*) 'it, ib', it, ib
            nf  = id(ib,1)
            ng  = id(ib,2)
            nh  = nf + ng
            k1  = ka(ib,it) + 1
            ke  = ka(ib,it) + kd(ib,it)
            k1a = ka(ib,it+2) + 1
            kea = ka(ib,it+2) + kd(ib,it+2)
            mul = mb(ib)
            m   = ib + (it-1)*NBX
            if (lpr.and.ib.eq.1) write(l6,'(/,a,1x,a)') tb(ib),tis(it)
c 
c------- calculation of ro
            do n2 =   1,nh
            do n1 =  n2,nh
            sr = zero
            do k = k1,ke
c------- Ravlic: finite-temperature
                if (temp.lt.1e-6) then
                        ftemp = 0
                else
                        ftemp = 1.d0/(1.d0+dexp(equ_v(k,it)/temp))
                endif
        sr = sr + fguv_v(nh+n1,k,it)*fguv_v(nh+n2,k,it)*(1.d0-ftemp)
     &                 + fguv_v(n1,k,it)*fguv_v(n2,k,it)*ftemp
            enddo    ! k
c----------------------------------
c            do k = k1a,kea                 ! no-sea approximation
c               sr = sr + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
c            enddo                          ! no-sea approximation
            sr = mul*sr
            rosh_v(n1+(n2-1)*nh,m) = sr
            rosh_v(n2+(n1-1)*nh,m) = sr
            enddo   ! n1
            enddo   ! n2

c--- Ravlic CHECK:
c         do n1 = 1,nh
c           do n2 = 1,nh
c            write(*,*) m,n1,n2,rosh(n1+(n2-1)*nh,m)
c           enddo
c         enddo

c
c------- contributions of large components f*f to kappa
            i0  = ia(ib,1)
            do n2 =   1,nf
            do n1 =  n2,nf
            ml1=ml(i0+n1)
            ml2=ml(i0+n2)
            i12 = 2 - n2/n1
            il  = il + 1
            sk = zero
            do k = k1,ke
c------- Ravlic: finite-temperature
                if (temp.lt.1e-6) then
                        ftemp = 0
                else
                        ftemp = 1.d0/(1.d0+dexp(equ_v(k,it)/temp))
                endif
        sk = sk - fguv_v(nh+n1,k,it)*fguv_v(n2,k,it)*ftemp
     &        + fguv_v(n1,k,it)*fguv_v(nh+n2,k,it)*(1.d0-ftemp)
            enddo    ! k
c----------------------------------
c            do k = k1a,kea                   ! no-sea approximation
c               sk = sk + fguv(nh+n1,k,it+2)*fguv(n2,k,it+2)
c            enddo                            ! no-sea approximation
            sk = mul*sk 
            if (n1.eq.n2) sp = sp + i12*sk
            if (ml1.ne.ml2) sk = zero      ! remove question ???
            aka_v(il,it) = i12*sk
c            write(*,*) it, il, aka(il,it)
            if (aka_v(il,it).ne.aka_v(il,it)) then
                    stop 'wrong in aka'
            endif
            enddo   ! n1
            enddo   ! n2
            spk_v(it)=half*sp

c
        if (lpr.and.ib.eq.1) then
            i0 = ia(ib,1) + 1
            write(l6,100) 'block ',ib,': ',tb(ib),' ',tit(it)
            nx = nh
            nx = 5
            call aprint(1,3,1,nh,nx,nx,rosh_v(1,m),tt(i0),tt(i0),'RO')
            nx = nf
            nx = 5
            call aprint(3,3,1,nh,nx,nx,aka_v(1,it),tt(i0),tt(i0),'AK')
        endif   
c
        enddo   ! ib

c
        if (lpr) then
        write(l6,*) ' ****** END DENSSH VAPOR************************'
        endif
c
  100 format(//,a,i2,4a)
c
c      read*
        return
C-end-DENSSH
        end


c=====================================================================c

      subroutine densit_vapor(lpr)

c=====================================================================c
C
c     calculates the densities in r-space at Gauss-meshpoints
c     this is for VAPOR
C
c---------1---------2---------3---------4---------5---------6---------7-
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
c
        logical lpr
c
        character tp*1,tis*1,tit*8,tl*1                         ! textex
        character nucnam*2                                      ! nucnuc
        character tb*6                                          ! blokap
        character tt*11                                         ! quaosc
c
        dimension rsh_v(2)
        dimension drs_v(MG,2),drv_v(MG,2)
c
        common /baspar/ hom,hb0,b0
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /coulmbv/ cou_v(MG),drvp_v(MG)
        common /defbas/ beta0,q,bp,bz
        common /densv / ro_v(MG,4),dro_v(MG,4)
        common /gaucor/ ww(MG)
        common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
        common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
        common /gfviv / iv(-IGFV:IGFV)
        common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
        common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
        common /mathco/ zero,one,two,half,third,pi
        common /nucnuc/ amas,npr(3),nucnam
        common /optopt/ itx,icm,icou,ipc,inl,idd
      common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
        common /rhorhov/ rs_v(MG,2),rv_v(MG,2)
        common /rokaosv/ rosh_v(NHHX,NB2X),aka_v(MVX,2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)      
c
        ihl(ih,il) = 1+ih +il*(NGH+1)
c
        if (lpr) then
        write(l6,*) ' ****** BEGIN DENSIT VAPOR***********************'
        endif
c
        ap2  = one/(b0*bp)**2
        az2  = one/(b0*bz)**2
c
        do it = 1,2
        do i = 1,MG
            rs_v(i,it)  = zero
            rv_v(i,it)  = zero
            drs_v(i,it) = zero
            drv_v(i,it) = zero
        enddo   ! i
        enddo   ! it
c
c---- loop over K-parity-blocks
        il = 0
        do ib = 1,nb
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
c 
c------- loop over contributions from large and small components
            do ifg = 1,2
            n0  = (ifg-1)*nf
            nd  = id(ib,ifg)
            i0  = ia(ib,ifg)
            ivv = iv(ifg)
c
c------- loop of oscillator basis states n2 and n1
            do n2 =  1,nd
            nz2 = nz(n2+i0)
            nr2 = nr(n2+i0)
            ml2 = ml(n2+i0)
            do 20 n1 =  n2,nd
            nz1 = nz(n1+i0)
            nr1 = nr(n1+i0)
            ml1 = ml(n1+i0)
            if (ml1.ne.ml2) goto 20
            ll = ml1**2
            i12 = 2 - n2/n1
            rsh_v(1) = rosh_v(n0+n1+(n0+n2-1)*nh,ib)*i12
            rsh_v(2) = rosh_v(n0+n1+(n0+n2-1)*nh,ib+NBX)*i12

c-------    loop over the mesh-points
            ben  = az2*(nz1+nz2+1) + ap2*2*(nr1+nr2+ml1+1)
            do il = 0,NGL
                qlab  = ql(nr1,ml1,il)*ql(nr2,ml1,il)
                qltab = ap2*(ql1(nr1,ml1,il)*ql1(nr2,ml1,il) + 
     &                      ll*qlab/xl(il))
                bfl   = ap2*xl(il) - ben
c
            do ih = 0,NGH
                qhab  = qh(nz1,ih)*qh(nz2,ih)
                qh1ab = az2*qh1(nz1,ih)*qh1(nz2,ih)
                sro   = qlab*qhab
                stau  = qh1ab*qlab + qhab*qltab
                sdro  = sro*(az2*xh(ih)**2 + bfl)
c
c-------       scalar and vector and pairing density
                do it = 1,2
                fgr = rsh_v(it)*sro
                rs_v(ihl(ih,il),it) = rs_v(ihl(ih,il),it) - ivv*fgr
                rv_v(ihl(ih,il),it) = rv_v(ihl(ih,il),it) + fgr
c
c-------          delta rho
                sdt = 2*rsh_v(it)*(sdro+stau)
                drs_v(ihl(ih,il),it) = drs_v(ihl(ih,il),it) - ivv*sdt
                drv_v(ihl(ih,il),it) = drv_v(ihl(ih,il),it) + sdt
                enddo   ! it
            enddo   ! ih
            enddo   ! il

   20    continue   ! n1
            enddo   ! n2
            enddo   ! ifg
        enddo   ! ib
c
c---- check, whether integral over dro vanishes
        s1 = zero
        s2 = zero
        do i = 1,MG
            s1 = s1 + drv_v(i,1)
            s2 = s2 + drv_v(i,2)
        enddo
        if (lpr) write(l6,*) ' integral over dro (V)',s1,s2
c
c
c---- normalization and renormalization to particle number
        do it = 1,2
            s  = zero
            do i = 1,MG
            s  =  s + rv_v(i,it)
            enddo
        if (lpr) write(l6,'(a,i3,2f15.8)') '  Integral over rv:',it,s
c            s  = npr(it)/s
c            do i = 1,MG
c            rv(i,it)  = s*rv(i,it)
c            rs(i,it)  = s*rs(i,it)
c            drs(i,it) = s*drs(i,it)
c            drv(i,it) = s*drv(i,it)
c            enddo
c
c------- printout of the density in configuration space
c         if (lpr) then
c            call prigh(2,rs(1,it),one,'RS  '//tis(it))
c            call prigh(2,rv(1,it),one,'RV  '//tis(it))
c            call prigh(2,drs(1,it),one,'DROS'//tis(it))
c            call prigh(2,drv(1,it),one,'DROV'//tis(it))
c         endif
        enddo  ! it
c
        do i = 1,MG
            f        = one/ww(i)
            ro_v(i,1)  = f*( + rs_v(i,1) + rs_v(i,2) )
            ro_v(i,2)  = f*( + rv_v(i,1) + rv_v(i,2) )
            ro_v(i,3)  = f*( - rs_v(i,1) + rs_v(i,2) )
            ro_v(i,4)  = f*( - rv_v(i,1) + rv_v(i,2) )
            dro_v(i,1) = f*(   drs_v(i,1) + drs_v(i,2) )
            dro_v(i,2) = f*(   drv_v(i,1) + drv_v(i,2) )
            dro_v(i,3) = f*( - drs_v(i,1) + drs_v(i,2) )
            dro_v(i,4) = f*( - drv_v(i,1) + drv_v(i,2) )
            drvp_v(i)  = drv_v(i,2)
        enddo   ! i
c      if (lpr) then
c         call prigh(1,ro(1,1),one,'RO-sig')
c         call prigh(1,ro(1,2),one,'RO-ome')
c         call prigh(1,ro(1,3),one,'RO-del')
c         call prigh(1,ro(1,4),one,'RO-rho')
c         ix = 5
c         do i = 1,ix
c            write(6,100) i,(ro(i,m),m=1,4)
c         enddo   ! i
c      endif   
c
  100 format(i5,4f15.10)
c
        if (lpr) then
        write(l6,*) ' ****** END DENSIT VAPOR********************'
        endif
c
        return
C-end-DENSIT
        end


c======================================================================c
c
      subroutine gdd_vapor(lpr)
c
c----------------------------------------------------------------------c
c
c     calculates 
c        fmes(x,1)       density-dependent coupling constants
c     and
c        fmes(x,2)       their derivatives
c
c     this is for VAPOR
c======================================================================c
        implicit real*8(a-h,o-z)
c
        include 'dirhb.par'
c     
        logical lpr

        common /couplfv/ ff_v(MG,4,2)
        common /densv / ro_v(MG,4),dro_v(MG,4)
        common /couplg/ ggmes(4),lmes(4)
        common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c     
        fun1(x,a,b,c,d) = a*(1+b*(x+d)**2)/(1+c*(x+d)**2)   ! Typel-Wolter
        dun1(x,a,b,c,d) = 2*a*(b-c)*(x+d)/(1+c*(x+d)**2)**2 ! Typel-Wolter
    
        fun2(x,a,b,c,d) = exp(-a*(x-one))                   ! Typel-Wolter
        dun2(x,a,b,c,d) = -a*exp(-a*(x-one))                ! Typel-Wolter
        
        fun3(x,a,b,c,d) = a+(b+c*x)*exp(-d*x)               ! DD-PC1
        dun3(x,a,b,c,d) = c*exp(-d*x)-d*(b+c*x)*exp(-d*x)   ! DD-PC1

        if (lpr) then
        write(l6,*) '****** BEGIN GDD VAPOR**************************'
        endif

        if (ipc.eq.0) then
            do i = 1,MG
            x=ro_v(i,2)/dsat
            do m = 1,2
                ff_v(i,m,1) = fun1(x,a_m(m),b_m(m),c_m(m),d_m(m))
                ff_v(i,m,2) = dun1(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
            do m = 3,4
                ff_v(i,m,1) = fun2(x,a_m(m),b_m(m),c_m(m),d_m(m))
                ff_v(i,m,2) = dun2(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
            enddo   ! i
        elseif (ipc.eq.1) then
            do i = 1,MG
            x=ro_v(i,2)/dsat
            do m = 1,4
                ff_v(i,m,1) = fun3(x,a_m(m),b_m(m),c_m(m),d_m(m))
                ff_v(i,m,2) = dun3(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
            enddo   ! i
        endif   ! ipc
c
c        if (lpr) then
c            ix = 5
c            do i = 1,ix
c            write(6,100) i,(ff_v(i,m,1),m=1,4)
c            enddo   ! i
c        endif   ! lpr
c 
        if (lpr) then
        write(l6,*) '****** END GDD VAPOR*****************************'
        endif
c
  100 format(i4,4f12.6)
c
        return
c-end-GDD
        end

c======================================================================c

      subroutine field_vapor(lpr)

c======================================================================c
c
c     calculation of the meson-fields in the oscillator basis
c     the fields are given in (fm^-1)
c
c     meson fields:  sig(i) = phi(i,1)*ggsig/gsig
c                    ome(i) = phi(i,2)*ggome/gome
c                    del(i) = phi(i,3)*ggdel/gdel
c                    rho(i) = phi(i,4)*ggrho/grho
c
c       VAPOR fields
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical lpr,lprs  
c
        dimension gm(4),so(MG),ph(MG)
c
        common /couplfv/ ff_v(MG,4,2)
        common /couplg/ ggmes(4),lmes(4)
        common /couplm/ gmes(4)
        common /densv / ro_v(MG,4),dro_v(MG,4)
        common /fieldsv/ phi_v(MG,4)
        common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /physco/ hbc,alphi,r0
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
        if (ipc.eq.1) return

        if (lpr) then
        write(l6,*) '****** BEGIN FIELD VAPOR***********************'
        endif

c---- loop over mesons
        do imes = 1,4
            do i = 1,MG
            so(i) = ff_v(i,imes,1)*ro_v(i,imes)
        enddo
            call gordon(imes,so,phi_v(1,imes))
        enddo   ! imes

c        if (lpr) then
c            write(l6,102)
c102    format(/,6x,'sigma',10x,'omega',10x,'rho  ')
c            do i = 1,MG
c            write(*,103) i,phi(i,1)*gm(1),phi(i,2)*gm(2),
c     &                     phi(i,3)*gm(3)
c            enddo   ! i
c        endif
C   
        if (lpr) then
        write(l6,*) '****** END FIELD VAPOR*******************'
        endif
cc
c     CALL zeit(2,'FIELD___',6,.false.)
c
  103 format(i3,4f15.10)
c
        return
C-end-FIELD
        end


c======================================================================c
c
       subroutine efield_vapor(emes,er,ecou)
c
c======================================================================c
c
c     calculates  field energies
c
c     VAPOR phase
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        dimension emes(4)
c
        common /coulmbv/ cou_v(MG),drvp_v(MG)
        common /coupld/ ddmes(4)
        common /couplfv/ ff_v(MG,4,2)
        common /couplg/ ggmes(4),lmes(4)
        common /densv / ro_v(MG,4),dro_v(MG,4)
        common /fieldsv/ phi_v(MG,4)
        common /gaucor/ ww(MG)
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /rhorhov/ rs_v(MG,2),rv_v(MG,2)
        common /rhorho/ rs(MG,2),rv(MG,2)
        common /physco/ hbc,alphi,r0
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
c======================================================================c
c---- field energies
c======================================================================c
c---- meson-fields
        if (ipc.eq.0) then
            er = zero
        do m = 1,4
            s = zero
            do i = 1,MG
            s  = s  + ggmes(m)*ff_v(i,m,1)*phi_v(i,m)*ro_v(i,m)*ww(i)
        er = er + ggmes(m)*ff_v(i,m,2)*phi_v(i,m)*ro_v(i,m)*ro_v(i,2)
     &                   *ww(i)
            enddo   ! i
            emes(m) = half*hbc*s
        enddo   ! m
        er = hbc*er        ! rearrangement term
c
c     point coupling
        elseif (ipc.eq.1) then
            er = zero
            do m = 1,4
            s = zero
            do i = 1,MG
            s  = s  + ggmes(m)*ff_v(i,m,1)*ro_v(i,m)**2*ww(i)
            er = er + ggmes(m)*ff_v(i,m,2)*ro_v(i,m)**2*ro_v(i,2)*ww(i)
            enddo   ! i
c
c          derivative terms
            do i = 1,MG
                s = s + ddmes(m)*ro_v(i,m)*dro_v(i,m)*ww(i)
            enddo   ! i
            emes(m) = half*hbc*s
        enddo   ! m
        er = half*hbc*er    ! rearrangment term
c
        else
        stop 'in EFIELD: ipc not properly defined'
        endif   ! ipc
c
c======================================================================c
c---- Coulomb energy
c======================================================================c
        ecou  = zero
        if (icou.ne.0) then
            do i = 1,MG
            ecou  = ecou + cou_v(i)*(rv(i,2)-rv_v(i,2))
            enddo   ! i
        endif   ! icou
        ecou  = half*hbc*ecou
c
        return
c-end-EFIELD
        end


c======================================================================c

        subroutine poten_vapor(lpr)

c======================================================================c
c
c     CALCULATION OF THE POTENTIALS AT GAUSS-MESHPOINTS
c     for VAPOR phase
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical lpr
c
        dimension glt(4)
c
        common /baspar/ hom,hb0,b0
        common /constr/ vc(MG,2)
        common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
        common /coulmbv/ cou_v(MG),drvp_v(MG)
        common /coupld/ ddmes(4)
        common /couplfv/ ff_v(MG,4,2)
        common /couplg/ ggmes(4),lmes(4)
        common /densv / ro_v(MG,4),dro_v(MG,4)
        common /fieldsv/ phi_v(MG,4)
        common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /physco/ hbc,alphi,r0
        common /potpotv/ vps_v(MG,2),vms_v(MG,2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
        common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)

        if (lpr) then
        write(l6,*) ' ****** BEGIN POTEN VAPOR**********************'
        endif
c
        do i = 1,MG
c------- meson-fields
            if (ipc.eq.0) then
            do m = 1,4
                glt(m) = ggmes(m)*ff_v(i,m,1)*phi_v(i,m)
            enddo   ! m
c
c           rearangement field
        re = zero
        do m = 1,4
            re = re + ggmes(m)*ff_v(i,m,2)*phi_v(i,m)*ro_v(i,m)
        enddo   ! m
            glt(2) = glt(2) + re
c
c------- point-coupling models
            elseif (ipc.eq.1) then
            do m = 1,4
                glt(m) = ggmes(m)*ff_v(i,m,1)*ro_v(i,m)
            enddo   ! m
c
c           derivative terms
            do m = 1,4
                glt(m) = glt(m) + ddmes(m)*dro_v(i,m)
            enddo   ! m
c
c           rearangement field
            re = zero
            do m = 1,4
                re = re + ggmes(m)*ff_v(i,m,2)*ro_v(i,m)**2
            enddo   ! m
            glt(2) = glt(2) + half*re
c
            else
            stop 'in POTEN: ipc not properly defined'
            endif   ! ipc
c
            s1 = hbc*(glt(1) - glt(3))                 ! neutron scalar
            s2 = hbc*(glt(1) + glt(3))                 ! proton  scalar
            v1 = hbc*(glt(2) - glt(4))                 ! neutron vector
            v2 = hbc*(glt(2) + glt(4) + cou_v(i))        ! proton  vector
c
c------- constraining potential
           if (icstr.gt.0) then
            v1 = v1 + vc(i,1)
            v2 = v2 + vc(i,2)
           endif   ! icstr
c
            vps_v(i,1) = v1 + s1
            vps_v(i,2) = v2 + s2
            vms_v(i,1) = v1 - s1
            vms_v(i,2) = v2 - s2

        enddo   ! i

c---- write potentials to a file, assume sph. symmetry
!       open(88, file = 'vps_vapor.out', status = 'unknown')
!       open(89, file = 'vms_vapor.out', status = 'unknown')

!       do il = 0, NGL
!         do ih = 0,NGH

!             ihl = 1+ ih +il*(NGH+1)
!             r = rb(il)
!             z = zb(ih)

!             !if (z.eq.r) then
!         write(88,'(4E15.7)') r, z, vps_v(ihl,1), vps_v(ihl,2)
!         write(89,'(4E15.7)') r, z, vms_v(ihl,1), vms_v(ihl,2)
!             !endif

!         enddo
!       enddo

!       close(88)
!       close(89)

        if (lpr) then
        write(l6,*) ' ****** END POTEN VAPOR*********************'
        endif
c
        return
c-end-POTEN
        end



C======================================================================c

      subroutine delta_vapor(it,lpr)

c======================================================================c
c
c     calculates the pairing field 
c     for separable pairing: Tian,Ma,Ring, PRB 676, 44 (2009)
c 
c     this is for VAPOR
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical  lpr
c
        character tb*6                                            ! blokap
        character tt*11                                            ! quaosc
        character tp*1,tis*1,tit*8,tl*1                           ! textex
c
        dimension pnn(NNNX)
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /quaosc/ nt,nnn(NTX,5),tt(NTX)
        common /deldelv/ de_v(NHHX,NB2X)
        common /mathco/ zero,one,two,half,third,pi
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /pair  / del(2),spk(2),spk0(2)
        common /rokaosv/ rosh_v(NHHX,NB2X),aka_v(MVX,2)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
        common /tmrwnn/ wnn(MVX,NNNX),nnmax
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)
        common /tmrpar/ gl(2),ga

        if (lpr) then
        write(l6,*) ' ****** BEGIN DELTA VAPOR***********************'
        endif
c
c----------------------------------------------------------------------c
c     calculation of PNN
c----------------------------------------------------------------------c
        do nn = 1,nnmax
c        nx = 5
c        call aprint(3,1,1,nh,nx,nx,aka(1,it),' ',' ','KA++')
c        call aprint(3,1,1,nh,nx,nx,wnn(1,nn),' ',' ','WNN++')
            s = zero
            do i = 1,mv
            if (wnn(i,nn).ne.wnn(i,nn)) then
                    write(*,*) 'NaN in wnn !'
            endif
            s = s + wnn(i,nn)*aka_v(i,it)
            enddo   ! i
            pnn(nn) = s
            if (pnn(nn).ne.pnn(nn)) then
                    stop 'NaN in pnn !'
            endif
        enddo  ! nn
c      write(*,*) 'First part: ', mv, nnmax, NNNX      
c     if (lpr) then
c        write(6,*) 'PNN ',mv,nnmax
c        do nn = 1,nnmax
c           write(6,101) ' nn =',nn,' PNN =',pnn(nn)
c        enddo   ! nn
c        write(l6,100) 'pnn',nnmax,(pnn(nn),nn=1,8)
c     endif
c
c      write(*,*) 'Pairing-field matrix elements:'
        g = half*gl(it)
        i12 = 0
        do ib = 1,nb
c         write(*,*) ib
            i0 = ia(ib,1)
            nf = id(ib,1)
            ng = id(ib,2)
            nh = nf + ng
            m  = ib + (it-1)*NBX
            do n2 =  1,nf
            do n1 = n2,nf
            i12 = i12 + 1
            s = zero
            do nn = 1,nnmax
                s = s + wnn(i12,nn)*pnn(nn)
c               write(*,*) ib, n1, n2, nn, wnn(i12,nn),pnn(nn)  
            enddo   ! nn
            de_v(n1+(n2-1)*nh,m) = -g*s
            de_v(n2+(n1-1)*nh,m) = -g*s
            enddo  ! n1
            enddo  ! n2 

c         do n1 = 1,nf
c           do n2 = 1,nf
c              write(*,*) m, n1, n2, de(n1+(n2-1)*nh,m)
c           enddo
c         enddo
c
c        if (lpr.and.ib.eq.1) then
c            k0 = ia(ib,1)+1
c            nx = nf
c            nx = 5
c            call aprint(1,3,1,nh,nx,nx,de(1,m),tt(k0),tt(k0),'DE++')
c        endif
        enddo   ! ib
c
        if (lpr) then
        write(l6,*) ' ****** END DELTA VAPOR********************'
        endif
c
  100 format(a,i6,8f10.6)
  101 format(a,i5,a,f15.6)
c
c      read*
        return
C-end-DELTA
        end


c======================================================================c

        subroutine dirhb_vapor_original(it,lpr)

c======================================================================c
c
c     solves the RHB-Equation 
c     IT    = 1 for neutrons
c     IT    = 2 for protons
c
c     solves for chem. pot. as in the orig. code 
c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical lpr,lprl
c
        character*1 bbb
        character*8 tbb(NHBX)
        character tp*1,tis*1,tit*8,tl*1                           ! textex
        character tb*6                                            ! blokap
        character tt*8                                            ! quaosc
        character nucnam*2                                        ! nucnuc
c
        dimension hb(NHBQX),e(NHBX),ez(NHBX)
c
        common /blodir/ ka(NBX,4),kd(NBX,4)
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /deldelv/ de_v(NHHX,NB2X)
        common /fermi / ala(2),tz(2)
        common /gamgamv/ hh_v(NHHX,NB2X)
        common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
        common /mathco/ zero,one,two,half,third,pi
        common /nucnuc/ amas,npr(3),nucnam
        common /optopt/ itx,icm,icou,ipc,inl,idd
        common /quaosc/ nt,nnn(NTX,5),tt(NTX)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)
        common /temp/ temp
c---- vapor phase solutions:
        common /waveuvv/ fguv_v(NHBX,KX,4),equ_v(KX,4)

c
        data maxl/200/,epsl/1.d-8/,bbb/'-'/,lprl/.false./
c
        if (.true.) then
        write(l6,*) ' ****** BEGIN DIRHB VAPOR ORIGINAL******'
        write(*,*) ' ****** BEGIN DIRHB VAPOR ORIGINAL*******'
        endif

        al = ala(it)
        sn = zero

        klp = 0
        kla = 0
c======================================================================c
        do ib = 1,nb            ! loop over differnt blocks
c======================================================================c
                mul  = mb(ib)
                nf   = id(ib,1)
                ng   = id(ib,2)
                nh   = nf + ng
                nhb  = nh + nh
                m    = ib + (it-1)*NBX
c
c------- calculation of the RHB-Matrix:
                do n2 = 1,nh
                do n1 = n2,nh
                hb(   n1+(   n2-1)*nhb) =  hh_v(n1+(n2-1)*nh,m) 
                hb(nh+n1+(nh+n2-1)*nhb) = -hh_v(n1+(n2-1)*nh,m) 
                hb(nh+n1+(   n2-1)*nhb) =  de_v(n1+(n2-1)*nh,m)
                hb(nh+n2+(   n1-1)*nhb) =  de_v(n2+(n1-1)*nh,m) 
                enddo
                hb(   n2+(   n2-1)*nhb) =  hb(n2+(n2-1)*nhb) - al
                hb(nh+n2+(nh+n2-1)*nhb) = -hb(n2+(n2-1)*nhb)
                enddo
c
c------- Diagonalization:
                if (lpr) then
                i0f = ia(ib,1)  
                do n = 1,nh
                        tbb(n)    = tt(i0f+n)
                        tbb(nh+n) = tbb(n)
                enddo
                write(l6,'(/i3,a,1x,a)') ib,'. Block ',tb(ib)
                nx = nhb
                nx = 5
                call aprint(2,3,6,nhb,nx,nx,hb,tbb,tbb,'HB')
                endif
                call sdiag(nhb,nhb,hb,e,hb,ez,+1)
c
c------- store eigenvalues and wave functions
c------- particles
                ka(ib,it) = klp
                do k = 1,nf
                klp = klp + 1
                equ_v(klp,it) = e(nh+k)
                do n = 1,nhb
                fguv_v(n,klp,it) = hb(n+(nh+k-1)*nhb)
                enddo
c------- Ravlic: finite-temperature
                if (temp.lt.1e-6) then
                ftemp = 0
                else
                ftemp = 1.d0/(1.d0+dexp(equ_v(klp,it)/temp))
                endif
                v2 = zero
                do n = 1,nh 
                v2 = v2 + fguv_v(nh+n,klp,it)**2*(1.d0-ftemp)
     &            + fguv_v(n,klp,it)**2*ftemp 
                enddo
c----------------------------------
                if (v2.lt.zero) v2 = zero
                if (v2.gt.one)  v2 = one

                sn = sn +v2*mul
                enddo
                kd(ib,it) = klp - ka(ib,it)

c
c------- anti-particles - CHECK!
                ka(ib,it+2) = kla
                do k = 1,ng
                kla = kla + 1
                equ_v(kla,it+2) = e(ng-k+1) 
                do n = 1,nhb
                fguv_v(n,kla,it+2) = hb(n+(ng-k)*nhb)
                enddo
                v2 = zero
c            do n = 1,nh                   ! no-sea approximation
c               v2 = v2 + fguv(nh+n,kla,it+2)**2
c            enddo                         ! no-sea approximation
c            sn = sn + v2*mul
                enddo
                kd(ib,it+2) = kla - ka(ib,it+2)
c         write(*,*) 'N (anti-particles) = ', sn
c
c======================================================================c
        enddo   ! ib
c======================================================================c
         write(*,*) 'N (particles) in VAPOR = ', it,  sn
  101 format(i4,a,i4,3f13.8)
  113 format(i4,a,a1,3f13.8)
c

        if (.true.) then
        write(l6,*) ' ****** END DIRHB VAPOR ORIGINAL******'
        write(*,*) ' ****** END DIRHB VAPOR ORIGINAL*******'
        endif
        return
C-end-DIRHB
        end

c======================================================================c

        subroutine emit_neutron()
c
c
c
c       calculates lifetime for emitting neutron
c       NOTE: canonical transf. of vapor states is required
c
c======================================================================c
        implicit real*8 (a-h,o-z)
c
        include 'dirhb.par'
c
        logical lpr

        
        common /erwar / ea,rms,betg,gamg
        
        common /temp/ temp 
        common /fermi / ala(2),tz(2)
        common /mathco/ zero,one,two,half,third,pi
        common /rhorhov/ rs_v(MG,2),rv_v(MG,2)
        common /gaucor/ ww(MG)
        common /basnnn/ n0f,n0b
        common /baspar/ hom,hb0,b0
        common /masses/ amu,amsig,amome,amdel,amrho
        common /physco/ hbc,alphi,r0
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /liflif/ time_life, part_dens
        
c---- canonical transformation of vapor states
        common /eeecanv/ eecan_v(KX,4),decan_v(KX,4),vvcan_v(KX,4),
     &                fgcan_v(NHX,KX,4),ibkcan_v(KX,4)
        common /blocanv/ kacan_v(nbx,4),kdcan_v(nbx,4),nkcan_v(4)
c-----------------------------------------------

        hbar = 6.582119569e-22 ! MeV s
c---- calculate cross-section
        cross = pi*rms**2

c---- calculate gas density (only neutrons)
        do it = 1,1
        s  = zero

        do ih = 1,MG
                s  = s  + rv_v(ih,it)
        enddo

        enddo !it

        part_vap =  s

        

        rrmax = sqrt(2*dble(n0f))*b0 ! Note

        write(*,*) '--------------------------------'
        write(*,*) 'n0f = ', n0f
        write(*,*) 'b0 = ', b0
        write(*,*) 'rmax = ', rrmax
        write(*,*) '--------------------------------'

        volume = 4.d0/3.d0*pi*rrmax**3

        part_dens = part_vap/volume

        write(l6,*) '--------------------------------'
        write(l6,*) 'Particles in vapor: ', part_vap
        write(l6,*) 'Particle density: ', part_dens
        write(l6,*) '--------------------------------'



c---- calculate <v>

        part1 = 0.d0
        part2 = 0.d0

        do it = 1,1

        do n = 1, nkcan_v(it)

                if (eecan_v(n,it).gt.0.d0) then
                ek = (eecan_v(n,it)-ala(it))/temp

                if (ek.lt.-15.0) ftemp = 1.d0
                if (ek.gt.+15.0) ftemp = 0.d0
                if (ek.gt.-15.0.and.ek.lt.+15.0) then
                        ftemp = 1.d0/(exp(ek)+1)
                endif

                vel = sqrt(2.d0*eecan_v(n,it)/(amu*hbc))

                part1 = part1 + ftemp*vel*sqrt(eecan_v(n,it))
                part2 = part2 + ftemp*sqrt(eecan_v(n,it))


c            write(*,*) n, eecan_v(n,it),ftemp, vel, part1, part2
                endif



        enddo !n



        enddo !it

        av_vel = part1/part2

        width = part_dens*hbc*cross*av_vel

        write(l6,*) 'cross (fm^-2) = ', cross
        write(l6,*) 'average velocity = ', av_vel
        write(l6,*) '--------------------------------'
        write(l6,*) 'Calculated width (MeV) = ', width
        write(l6,*) 'Calculated lifetime (s) = ', hbar/width
        write(l6,*) '--------------------------------'
        time_life = hbar/width


        return
        end

c======================================================================c

        subroutine canon_vapor(lpr)

c======================================================================c
c
c     transforms to the canonical basis
c     version for RHB
c     NOTE that this is for vapor states

c----------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
        include 'dirhb.par'
c
        logical lpr,lpr1
c
        character tb*6                                           ! blokap
        character tt*11                                          ! quaosc
        character tp*1,tis*1,tit*8,tl*1                          ! textex
c
        dimension aa(NHHX),dd(NHHX),v2(NHX),z(NHX),eb(NHX),h(NHX),d(NHX)
c      
        common /blodir/ ka(NBX,4),kd(NBX,4)
        common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
        common /bloosc/ ia(NBX,2),id(NBX,2)
        common /gamgamv/ hh_v(NHHX,NB2X)
        common /deldelv/ de_v(NHHX,NB2X)
        common /eeecanv/ eecan_v(KX,4),decan_v(KX,4),vvcan_v(KX,4),
     &                fgcan_v(NHX,KX,4),ibkcan_v(KX,4)
        common /blocanv/ kacan_v(nbx,4),kdcan_v(nbx,4),nkcan_v(4)
        common /mathco/ zero,one,two,half,third,pi
        common /fermi / ala(2),tz(2)
        common /optopt/ itx,icm,icou,ipc,inl,idd
      common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
        common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
        common /textex/ tp(2),tis(2),tit(2),tl(0:30)
        common /waveuvv/ fguv_v(NHBX,KX,4),equ_v(KX,4)
        common /temp/ temp     
        
        data ash/100.d0/
c
        if (lpr) then
        write(l6,*) '****** BEGIN CANON VAPOR************************'
        endif 
c
c======================================================================c
        do it = 1,2   ! loop over neutrons and protons
c======================================================================c

                if (lpr) then
                write(l6,100) tit(it)
                write(l6,101) 'K pi','[nz,nr,ml]','smax',
     &                    'eecan','vvcan','decan'

                endif
                klp = 0
                kla = 0
c======================================================================c
        do ib = 1,nb   ! loop over the blocks
c======================================================================c
                nf  = id(ib,1)
                ng  = id(ib,2)
                nh  = nf + ng
                nhb = 2*nh 
                i0f = ia(ib,1)
                i0g = ia(ib,2)
                m   = ib + (it-1)*NBX
                kf  = kd(ib,it)
                kg  = kd(ib,it+2)
                k0f = ka(ib,it)
                k0g = ka(ib,it+2)
c
c------- transformation to the canonical basis
c------- calculation of the generalized density V*VT
                do n2 = 1,nh 
                do n1 = 1,nh 
                s = zero  
                do k = k0f+1,k0f+kf
                if (temp.lt.1e-6) then
                ftemp = 0
                else
                ftemp = 1.d0/(1.d0+dexp(equ_v(k,it)/temp))
                endif
        s = s + fguv_v(nh+n1,k,it)*fguv_v(nh+n2,k,it)*(1.d0-ftemp)
     &             + fguv_v(n1,k,it)*fguv_v(n2,k,it)*ftemp
                enddo ! k  
                s1 = zero
c------------ no-sea approximation ----------------------------
c                do k = k0g+1,k0g+kg
c                s1 = s1 + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
c                enddo
c--------------------------------------------------------------
                aa(n1+(n2-1)*nh) = s + ash*s1
                enddo   ! n1         
                enddo   ! n2
c       
c------- diagonalizaton
                call sdiag(nh,nh,aa,v2,dd,z,1)
                eps=1.0e-6
                call degen(nh,nh,v2,dd,hh_v(1,m),eb,eps,aa,z)
c
c        major component of the wave function should be > 0
                do k = 1,nh     
                cmax = zero
                do n = 1,nh
                s = abs(dd(n+(k-1)*nh))
                if (s.gt.abs(cmax)) cmax = dd(n+(k-1)*nh)
                enddo   ! n
                if (cmax.lt.zero) then
                do n = 1,nh
                        dd(n+(k-1)*nh) = - dd(n+(k-1)*nh)
                enddo   ! n
                endif
                enddo   ! k
c------- diagonal matrix elements of HH and DE in the canonical basis
                do k = 1,nh     
                hk = zero
                dk = zero
                do n2 = 1,nh
                h2 = zero
                d2 = zero
                do n1 = 1,nh
                h2 = h2 + dd(n1+(k-1)*nh)*hh_v(n1+(n2-1)*nh,m)
                d2 = d2 + dd(n1+(k-1)*nh)*de_v(n1+(n2-1)*nh,m)
                enddo   
                hk = hk + h2*dd(n2+(k-1)*nh)
                dk = dk + d2*dd(n2+(k-1)*nh)
                enddo
                h(k) = hk
                d(k) = dk
                enddo   ! k
c
c------- reordering according to the energy h(k)
                call ordx(nh,h,d,v2,dd)
c
                do k = 1,nh     
        if (v2(k).lt.zero .or. v2(k).gt.2.d0) v2(k) = zero
                if (v2(k).gt.one) v2(k) = one
                enddo ! k
c                  
                kacan_v(ib,it)   = klp
                do k=1,nf     
                klp=klp+1
                eecan_v(klp,it)=h(ng+k)
                decan_v(klp,it)=d(ng+k)
                vvcan_v(klp,it)=v2(ng+k)
                ibkcan_v(klp,it)=ib
                do n = 1,nh
                fgcan_v(n,klp,it)=dd(n+(ng+k-1)*nh)
                enddo
                enddo
                kdcan_v(ib,it)   = klp - kacan_v(ib,it)
c
                kacan_v(ib,it+2) = kla
                do k=1,ng
                kla=kla+1
                eecan_v(kla,it+2)=h(k)
                decan_v(kla,it+2)=d(k)
                vvcan_v(kla,it+2)=v2(k)
                ibkcan_v(kla,it+2)=ib
                do n=1,nh
                fgcan_v(n,kla,it+2)=dd(n+(k-1)*nh)
                enddo
                enddo
                kdcan_v(ib,it+2) = kla - kacan_v(ib,it+2)          
                
c------- printout for particles
                if (lpr) then
                if (ib.eq.1) e0 = h(ng+1)
                k1 = kacan_v(ib,it)+1
                k2 = kacan_v(ib,it)+kdcan_v(ib,it)
                if (kb(ib) .gt. 0) ip=1
                if (kb(ib) .lt. 0) ip=2
                do k = k1,k2
                e1 = eecan_v(k,it)
                v1 = vvcan_v(k,it)
                d1 = decan_v(k,it)
c           search for the main oscillator component
                smax = zero
                do i = 1,nf
                s = abs(fgcan_v(i,k,it))
                if (s.gt.smax) then
                        smax = s
                        imax = i
                endif
                enddo
                fx = fgcan_v(imax,k,it)**2
c
                write(l6,102) k,tb(ib),tt(i0f+imax),fx,e1,v1,d1

                enddo  ! k
                endif   ! lpr
c======================================================================c
   10   enddo      ! ib      end loop over the blocks
c======================================================================c
                nkcan_v(it)   = klp
                nkcan_v(it+2) = kla         
c======================================================================c
        enddo      ! it      end loop over neutrons and protons --------c
c======================================================================c
c
        if (lpr) then
        write(l6,*) '****** END CANON VAPOR**************************'
        endif
c
  100 format(' single-particle energies and gaps ',1x,
     &       'in the canonical basis: ',a,/1x,66(1h-))  

  101 format(7x,a,a,2x,a,5x,a,5x,a,5x,a)
  102 format(i4,3x,a6,a10,f7.2,5f10.4)
c
        return
C-end-CANON
        end
