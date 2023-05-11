       program ptm
c analyze multiplicity data for SIDIS
c do do:

c finish giving wide narrow files different names!

c put smooth interpolation into getrho

c diff. in HMS yptar < 0 and > 0 by 6%.

c multiplicities bigger for phi = pi/2 and 3pi/3 by
c about 10%! (averaged over low pt) when using wide
c acceptance limits.
c Try narrow limits. Could be incorrect HMS offset?

      implicit none

      integer i,j,jp,k,npt,ifit,kk,ipm,kkk,jj,itt,ix,iy,nplt
      integer ikin,ipt,iphi,iz,it,iacc,ikinsv,izsv,ikinp
      integer ngood,ikinx,ipdf,iff,i1,i2,ikinpp
      integer ikinw(8)/29, 13, 9, 5, 11, 1, 7, 3/
      real*8 xpltmin,xpltmax,ypltmin,ypltmax
      real*8 mcsv(3,8,20,10),mcsver(3,8,20,10) 
      real*8 fu1sv(3,8,20),fd1sv(3,8,20),avdelu(3,8),avdeld(3,8)
      real*8 avdeluer(3,8),avdelder(3,8),doverp,doverpp,doverpper
      real*8 mcsvp(3,8,20,12),mcsvper(3,8,20,12),dopsv(3,8,20)
      real*8 mcsv5p(3,8,20,5),mcsv5per(3,8,20,5),uplusd(3,8)
      real*8 delu, deluer, deld, delder,rnp,rsv(20)
      real*8 Mm, Mp, Mmer, Mper, Bm, Bp,rnew,dprop,wp2,fact,du
      real*8 mprat, rd1, rd2, drd, rseas(16,20), rseans(16,20)
      real*8 csvprop(3,20,3,2), csvproper(3,20,3,2),f,Ctq5Pdf
      real*8 BigS, chk1er1, chk1er2,chk1er3,chk1er4,BigT
      real*8 rcsv, rcsver,fu1j,fd1j
      real*8 sumdifffit(56,20,2,2,2)
      real*4 qgev,x4
      logical lastframeused,usewide,okplt
      character*2 ttit(4)/'p+','d+','p-','d-'/
      real*8 x,q2,zmax,w,pt,phicm,z,mmpi2,phi,ratt,wpcut,chk1,chk2
      real*8 sighadp,sighadm,rdm,scsv(56,4,20,4,3),xp,zp,b
      real*8 sighad1, sighad2, sighad3, sighad4,pt2,sig,siger
      real*8 scsver(56,4,20,4,3),xmax,xmin,ymin,ymax,sighad5
      real*8 avcsv(16,3,2),avcsver(16,3,2),pi/3.1415928/
      real*8 ffrat(20,4),ffrater(20,4)
      real*8 ffrat4(20,4),ffrat4er(20,4)
      real*8 yymin,yymax, zforwp225,zforwp23,u1,d1
      real*8 q2n,q2er,xn,xer,zn,zer,rd,rder,rdrho
      real*8 mltnorc, mltnorcer,mlt,mlter,mrho,xx
      real*8 sighad,u,ub,d,db,s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu
      real*8 w2,f1p,f2p,f1d,f2d,a,zz,rplus,rminus,zcut
      real*8 rplusd,rminusd,rat(6),ratrho(6),rater(6)
      real*8 rv(5,20),rver(5,20)
      real*8 avm(56,20,6),avmer(56,20,6),avk(56,20,5)
      real*8 avmod(56,20,6),mmod,mratio
      real*8 avker(56,20),avmrho(56,20,6)
      real*8 avmk(6),avmker(6),lund(16,20,4)
      real*8 avmkr(6),avmkrer(6)
      real*8 avms(56,20,6),rats(6),diff,diffmax
      real*8 srat(56,20,6),srater(56,20,6)
      real*8 sratrho(56,20,6)
      real*8 sq2(10,2),sx(10,2),sz(10,2),smx(10,2)
      real*8 pq2(10,2),px(10,2),pz(10,2),pmx(10,2)
      real*8 pkin(56,2),pw(10,2),rhofactp
      integer rhorat(25,2),usedup(56,20,6)
      integer bnpt,bnpt0,bitv(90000)
      integer srath(40,6),sratrhoh(40,6),rdh(40,2)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000),avmit(5)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      real*8 bchi2x(15),bdfx(15),avit(5),avitrho(5),aviter(5) 
      real*8 ffsv(3,3,20,4),ffsver(3,3,20,4)
      real*8 avffsv(20,4,0:1),avffsver(20,4,0:1)
      real*8 avffsv4(20,4,0:1),avffsv4er(20,4,0:1)
      real*8 avffsvcsv(20,4),avffsvcsver(20,4)
      real*8 avffsvex(20,4),avffsvexer(20,4)
      real*8 avffsv4lund(20,4),avffsv4lunder(20,4)

      integer ncallsfit
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      real*8 rhofact,fbest,chibest,fnorm(36),fnormsv(36)
      integer bkin(90000),ncall/0/
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt,ncallsfit
        integer ikinfit,itfit,npts
        real*8 zfit
        common/sstuff/ zfit,ikinfit,itfit,npts
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag
        integer sn
        real*8 sphi(1000), spt(1000),sy(10000),syer(1000)
        common/sforplot/ sphi,spt,sy,syer,sn

        real*8 ep0sv(28)/
     >    5.24,   3.31,  5.24,  4.49,  5.98,
     >    4.94,   5.98,  6.36,  5.24,  6.36,
     >    4.72,   5.24,  5.71,  4.33,  5.24,
     >    4.87,   3.98,  4.78, 6.590, 5.262,
     >    4.686, 4.183, 3.253, 
     >    2.178, 0.888, 2.320,
     >    1.816, 0.962/
        real*8 the0sv(28)/
     >    13.50, 19.69, 16.31, 16.63, 14.24,
     >    17.26, 15.74, 15.30, 17.21, 13.96,
     >    19.06, 15.38, 17.65, 20.24, 18.51,
     >    19.10, 19.68, 19.69, 11.910, 11.159,
     >    17.08, 14.93,23.005, 
     >    27.25, 36.15, 27.77,
     >    25.895, 49.308 /
      real*8 e0,ep,th,sin2,cos2,nu,q2kin(56),xkin(56),wkin(56)
      real*8 e0sv(56),am/0.938/
      real*8 qu,qd,qs,sum_sqp,sum_sqn,w2chk
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)

      real*8 coef(100)  ! final values of coeficients
      integer nparam    ! number of parameters
      integer arglis(10) ! to pass info to Minuit
      real*8 std(100),chi2,grad(100),p(100)
      integer  nerror_migrad,ierflg,npar
      real*8 zero/0./  ! to pass to Minuit a 0.
      character*10 pname(100)
      external bfit_fcn,sfit_fcn,mfit_fcn,mfit5_fcn,fffit_fcn 
      external fffit4_fcn,fffitcsv_fcn,fffitex_fcn 
      real*8 futil ! auxially function
      external futil
      real*8 mmpi2cut,q2cut,wcut,ptcut,ptcutlo
c inclusive d/p from bcm.out using scalers
      real*8 rrsv(8)/0.86, 0.86, 0.81, 0.86,
     >  0.80, 0.81, 0.78, 0.77/
c shuo data
      character*200 string
       real*8 avq2,avq2er,avx,avxer,zset,avz,avzer,yrd,q2set,xset,zbin,
     >  yrder,yrdrho,ry,ryrho,ryrho2,rynoex,rynodel,ryerr,xbin
       real*8 sjryav(32,20),sjryaver(32,20),v1,v2,v3,v4,v5,v6
       real*8 sjryrhoav(32,20)
       integer group,ncallff,iswap,ioff

       real*8 sum1,sum2,phimin,phimax
       real*8 mvlund(4),chi2ff
       real*8 mv(4),mver(4),muf,mubf,mdf,mdbf,msf,msbf,mvf(4)
       common/ffstuff/ mv,mver,muf,mdf,mubf,mdbf,msf,
     >  msbf,fs1,fsb,mvf,ncallff,iswap,chi2ff

c this is for DSS
      integer IHdss,ICdss,IOdss, fini 
      real*8 fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1
      COMMON / FRAGINI / FINI
      fini=0

c note: low cut values!
      mmpi2cut = 2.0
      wcut=1.8
      q2cut=1.8
      zcut=0.96
c note: very low cut values!
      mmpi2cut = 1.1
c more reasonable cut
c      mmpi2cut = 2.56
      wcut=2.0
      q2cut=1.0
      zcut=0.95
c  can chang from default 
c      ptcut = 3./16.
c range of pt to use for averages
c large value!
      ptcut = 4./16.
      ptcutlo = 2./16.
c wide or narrow acceptance?
      usewide=.true.
c      usewide = .false.

! ASSIGN UNITS (NORMALLY 5,6,7 FOR INTERACTIVE)
      CALL MNINIT(61,62,63)

      if(usewide) then
       open(unit=6,file='ptmw.out')
       open(unit=16,file='ptmsfitw.out')
       open(unit=17,file='ptmnfitw.out')
      else
       open(unit=6,file='ptm.out')
       open(unit=16,file='ptmsfit.out')
       open(unit=17,file='ptmnfit.out')
      endif

c read in Lund multiplicities
      open(unit=9,file='ptm.lund')
      do ikin=1,8
       read(9,'(a)') string
       read(9,'(a)') string
       read(9,'(a)') string
       read(9,'(a)') string
       do iz=1,19
        read(9,*) i1,i2,lund(ikin,iz,1), ! p pi+
     >    lund(ikin,iz,2),! p pi-
     >    v1,lund(ikin,iz,3), ! n pi+
     >    lund(ikin,iz,4) ! n pi-
       enddo
      enddo
      close(unit=9)

c check of CSV formula
      do i=1,3
       rnew = 0.25 * i
       do j=1,6
        du = 0.1 * (j-1)
        mprat = (4.*(1-du)*rnew +  (1 + du))/
     >          (4.*(1-du) + (1 +du)*rnew)
        Dprop = (1. - rnew) / (1. + rnew)
        RD = (4.*mprat - 1.) / (1. - mprat)
c this gives delu / (u + d + ub + db), if s=sbar=0
        chk1 = (4*rnew + 1 - mprat * (4 + rnew)) /
     >         (4*rnew - 1 - mprat * (4 - rnew))
c this is -4/3 (delta u - delta d) / (u + d)
c or -8/3 delta u / (u+d)
        chk2 = 2.5  - Dprop * (2.5 + RD)
        write(6,'(''chkcsv'',6f7.2)') 
     >   rnew, du, chk1, 
     >    -3./8.*chk2
        enddo
       enddo
      do i=1,9
       z = 0.1*i
       q2 = 2.
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
       write(6,'(''fdss'',8f7.4)') z,q2, fu1,fd1/fu1,
     >  fub/fu1, fdb/fu1, fs1/fu1, fsb/fu1
c ub is same as d1
c db is a bit bigger than u1
c s1 and sb are about same as d1
      enddo

      do ikin=1,56
        e0sv(ikin)=10.600
        if(i kin.ge.25 .and. ikin.le.32) e0sv(ikin)=10.214
        if(ikin.ge.47 .and. ikin.le.50) e0sv(ikin)=6.190
        if(ikin.ge.51 .and. ikin.le.56) e0sv(ikin)=8.209
      enddo
c get kinematics by setting
      do ikin=1,56
       e0 = e0sv(ikin)
       ep = ep0sv( (ikin+1) / 2)
       th = 3.14159 / 180. * the0sv( (ikin+1) / 2)
       if( (ikin/2)*2 .eq. ikin) then
        th = th + 0.010
       else
        th = th - 0.010
       endif
       sin2 = sin(th/2)**2
       cos2 = 1. - sin2
       q2kin(ikin) = 4. * e0 * ep * sin2
       nu = e0 - ep
       xkin(ikin) = q2kin(ikin) / 2. / am / nu
       x = xkin(ikin)
       q2 = q2kin(ikin)
       w = sqrt(am**2 + 2.*am*nu -q2)
       wkin(ikin) = w
      enddo

c read in Shuo's results for yield ratios from d (pi-/pi+) in grid of
c 20 x by 20 z bins for each kin. setting
c      open(unit=9, file='sjcsv_simadd.txt')
      open(unit=9, file='csv_datasub.csv')
      read(9,*) string
      do i=1,669
       read(9,*) q2set,avq2,avq2er,xset,xbin,avx,avxer,zset,zbin,avz,
     >  avzer,group,yrd,yrder,yrdrho,v1,v2,
     >  ry,ryrho,ryrho2,rynoex,rynodel,
     >  ryerr
       ikin=0.
       do j=1,32,2
        if(abs(q2set - (q2kin(j)+q2kin(j+1))/2.).lt.0.3 .and.
     >     abs(xset - (xkin(j)+xkin(j+1))/2.).lt.0.03) ikin=j
       enddo
       if(ikin.eq.0) write(6,'(''error,ikin=0'')')
       iz = int(20.*avz)+1
       if(ryerr.lt.0.2) then
        sjryav(ikin,iz) = sjryav(ikin,iz) + ry/ryerr**2 
        sjryrhoav(ikin,iz) = sjryrhoav(ikin,iz) + 
     >    ryrho/ryerr**2 
        sjryaver(ikin,iz) = sjryaver(ikin,iz) + 1./ryerr**2 
       endif
       write(6,'(2i4,i3,10f7.2)') group,ikin,iz,xset,q2set,xbin,zbin,
     >  avz,ry,ryrho,ryerr
      enddo
      close(unit=9)
      do ikin=1,32
       do iz=1,20
        if(sjryaver(ikin,iz).ne.0.) then
         sjryav(ikin,iz) = sjryav(ikin,iz)/sjryaver(ikin,iz)
         sjryrhoav(ikin,iz) = sjryrhoav(ikin,iz) /
     >     sjryaver(ikin,iz)
         sjryaver(ikin,iz) = sqrt(1./sjryaver(ikin,iz))
        endif
       enddo
      enddo

c wide limits
      if(usewide) then
       open(unit=277,file='ptm.inpw')
c narrow limits
      else
       open(unit=277,file='ptm.inpn')
      endif
      bnpt = 0
       do i=1,1000000
       read(277,'(i2,3i3,2i2,6f6.3,6f10.4,e12.4)',end=11) 
     >  ikin,ipt,iphi,iz,it,iacc,x,q2,
     >  z,pt,phicm,mmpi2,
     >  mlt,mlter,mltnorc,mltnorcer,mrho,
     >  mratio,mmod
       w = sqrt(0.938**2 + q2 * (1./x -1))
       if(mmpi2.gt.mmpi2cut .and. it.ne.3 .and. it.ne.6
     >  .and. ikin.le.56 .and. z.lt.zcut
c for test
c     >  .and. ikin.gt.2
     >  .and. mlter.gt.0.
     >  .and. q2.gt.q2cut .and. w.gt.wcut) then
c special corr. for f1f2in21 being 5% begger than
c f1f2in09 for deuteron.
c took out because now using f1f2in09 in SIMC
c        if(it.eq.2 .or. it.eq.5) then
c         mlt = mlt * 1.05
c         mlter = mlter * 1.05
c        endif
c for test
c        mlt = mltnorc
c        mlter = mltnorcer

c for test
c        if(it.eq.2) then
c         mlt = mlt * 1.065
c         mlter = mlter * 1.065
c        endif
             bnpt = bnpt + 1
             bkin(bnpt) = ikin
             bitv(bnpt) = it
             bq2v(bnpt) = q2
             bxv(bnpt) = x
             bmmpi2(bnpt) = mmpi2
c epsilon not used
c             berv(bnpt) = erkin(ikin)
             bzv(bnpt) = z
             bptv(bnpt) = pt
             phi = 2. * 3.1415928 / 15. * (iphi-0.5)
             bphiv(bnpt) = phi
             byv(bnpt) = mlt
             byerv(bnpt) = mlter
             brho(bnpt) = mrho
        call getrho(xkin(ikin),q2kin(ikin),z,pt,phi,rplus,rminus,
     >    rplusd, rminusd)
        if(it.eq.1) brhop(bnpt) = mlt * rplus
        if(it.eq.4) brhop(bnpt) = mlt * rminus
        if(it.eq.2) brhop(bnpt) = mlt * rplusd
        if(it.eq.5) brhop(bnpt) = mlt * rminusd

        if(it.eq.1 .and. rplus .gt.0.02) then
         kkk = min(25,max(1,int(mrho/mlt/rplus*4 .)+1))
         rhorat(kkk,1) = rhorat(kkk,1)+1
        endif
        if(it.eq.4 .and. rminus.gt.0.02) then
         kkk = min(25,max(1,int(mrho/mlt/rminus*4.)+1))
         rhorat(kkk,2) = rhorat(kkk,2)+1
        endif

        ipdf = 1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c        write(6,'(''tst'',i6,4f10.3)') i,mlt,mlter,
c     >   sighad,mlt/sighad
        buv(bnpt) = u
        bubv(bnpt) = ub
        bdv(bnpt) = d
        bdbv(bnpt) = db
        bsv(bnpt) = s
        bsbv(bnpt) = sb
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        bffu(bnpt) = fu1 / z
        bffd(bnpt) = fd1 / z
        bffs(bnpt) = fs1 / z
        bffub(bnpt) = fub / z
        bffdb(bnpt) = fdb / z
        bffsb(bnpt) = fsb / z
        if(mrho/mlt.gt.0.2) write(6,'(''rho'',
     >   i2,9f6.2,f7.3)') 
     >   it,x,q2,z,pt,phicm,mmpi2,mlt, mlter,mlt/sighad,
     >   mrho/mlt
c ratio of pi-/pi+ deuteron, etc. 
c these are averages of FF (dsigdz, based on M0), not sighad
c they should match fDSS
        if(pt.lt.ptcut .and. pt.gt.ptcutlo) then
c used in SIMC. see v2.pdf
c         b = 1./ (0.120 * z**2 + 0.200)
c rough fit to pt-SIDIS results. see v1.pdf
         zmax = min(0.6, z)
         b = 1./ sqrt(0.16**2 + zmax**2 * 0.30**2)
         if(it.eq.1.or.it.eq.2) then
c          b = 1./ sqrt(0.15**2 + zmax**2 * 0.40**2)
         endif
        fact = 1./ (b * exp(-b * pt**2)/ 2. / pi)
c for test of systematics see v3.pdf
c         fact = 1.
         write(6,'(''tst'',4i3,5f8.3)') ikin,iz,ipt,iphi,
     >    mlt,mlter,
     >    sighad,mlt/sighad,mmod/sighad
        avm(ikin,iz,it) = avm(ikin,iz,it) +
     >    fact * mlt / mlter**2 / fact**2
        avms(ikin,iz,it) = avms(ikin,iz,it) +
     >    fact * sighad / mlter**2 / fact**2
        avmrho(ikin,iz,it) = avmrho(ikin,iz,it) +
     >    fact * mrho/ mlter**2 / fact**2
        avmer(ikin,iz,it) = avmer(ikin,iz,it) +
     >    1. / mlter**2 / fact**2
        avmod(ikin,iz,it) = avmod(ikin,iz,it) +
     >    fact * mmod / mlter**2 / fact**2
        avk(ikin,iz,1) = avk(ikin,iz,1) + x/ mlter**2 / fact**2 
        avk(ikin,iz,2) = avk(ikin,iz,2) +q2/ mlter**2 / fact**2 
        avk(ikin,iz,3) = avk(ikin,iz,3) + w/ mlter**2 / fact**2 
        avk(ikin,iz,4) = avk(ikin,iz,4) + z/ mlter**2 / fact**2 
        avk(ikin,iz,5) = avk(ikin,iz,5) + mmpi2/mlter**2 / fact**2 
        avker(ikin,iz) = avker(ikin,iz) + 1./mlter**2  / fact**2
        endif ! pt cut
       endif
      enddo
 11   continue

c results for pi-/pi+ deuteron vs z
      if(usewide) then
       open(unit=14,file='ptmratzw.txt')
       open(unit=15,file='ptmratzrhow.txt')
      else
       open(unit=14,file='ptmratz.txt')
       open(unit=15,file='ptmratzrho.txt')
      endif
      do ikin=1,56
       do iz=1,20
c super ratios of data to model 
        do it=1,6
         srat(ikin,iz,it)=0.
         srater(ikin,iz,it)=0.
         usedup(ikin,iz,it)=0
        enddo
        do it=1,6
          if(avmer(ikin,iz,it).ne.0.) then
           avm(ikin,iz,it) = avm(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmod(ikin,iz,it) = avmod(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avms(ikin,iz,it) = avms(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmrho(ikin,iz,it) = avmrho(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmer(ikin,iz,it) = 1./
     >       sqrt(avmer(ikin,iz,it))
          endif
        enddo
        if(avmer(ikin,iz,2).ne.0. .and.
     >      avker(ikin,iz).ne.0.) then
         do it=1,6
          rat(it)=0.
          rats(it)=0.
          ratrho(it)=0.
          rater(it)=0.
          if(avmer(ikin,iz,it).ne.0 .and.
     >       avmer(ikin,iz,2).ne.0.) then
           rat(it) = avm(ikin,iz,it) / 
     >               avm(ikin,iz,2)
           rats(it) = avms(ikin,iz,it) / 
     >                avms(ikin,iz,2)
           ratrho(it) = avmrho(ikin,iz,it) / 
     >               avmrho(ikin,iz,2)
           rater(it) = rat(it) * sqrt(
     >      (avmer(ikin,iz,it)/avm(ikin,iz,it))**2 +
     >      (avmer(ikin,iz,2)/avm(ikin,iz,2))**2)
           sratrho(ikin,iz,it) = ratrho(it)/rats(it)
           srat(ikin,iz,it) = rat(it)/rats(it)
           srater(ikin,iz,it) = rater(it)/rats(it)
           diff = (srat(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
           diff = max(-7.,min(6.99,diff))
           k = int((diff + 7.)/14.*40.)+1
           if(k.gt.40) write(6,'(i5,2f8.3)') k,
     >      srat(ikin,iz,it),srater(ikin,iz,it)
           srath(k,it) = srath(k,it)+1
           diff = (sratrho(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
           diff = max(-7.,min(6.99,diff))
           k = int((diff + 7.)/14.*40.)+1
           sratrhoh(k,it) = sratrhoh(k,it)+1

           if(it.eq.5) then
            diff = (sratrho(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
c            write(6,
c     >      '(''err diff'',4e10.2)') sratrho(ikin,iz,it),
c     >       srat(ikin,iz,it),srater(ikin,iz,it),diff
            x = avk(ikin,iz,1)/avker(ikin,iz)
            q2 = avk(ikin,iz,2)/avker(ikin,iz)
            w = avk(ikin,iz,3)/avker(ikin,iz)
            z = avk(ikin,iz,4)/avker(ikin,iz)
            mmpi2 = avk(ikin,iz,5)/avker(ikin,iz)
       k=max(1,min(10,int((q2-3.)/3.5*10)+1))
       pq2(k,1) = pq2(k,1) + diff
       pq2(k,2) = pq2(k,2) + 1.
       k=max(1,min(10,int((x-0.25)/0.35*10)+1))
       px(k,1) = px(k,1) + diff
       px(k,2) = px(k,2) + 1.
       k=max(1,min(10,int((z-0.3)/0.5*10)+1))
       pz(k,1) = pz(k,1) + diff
       pz(k,2) = pz(k,2) + 1.   
       k=max(1,min(10,int((w-1.8)/1.5*10)+1))
       pw(k,1) = pw(k,1) + diff
       pw(k,2) = pw(k,2) + 1.   
       k=max(1,min(10,int((mmpi2-2.)/5*10)+1))
       pmx(k,1) = pmx(k,1) + diff
       pmx(k,2) = pmx(k,2) + 1.   
           endif
          endif
         enddo
         if(rater(5).lt.0.2) then
          write(14,'(2i3,5f6.3,9f7.3)') ikin,iz,
     >     (avk(ikin,iz,jj)/avker(ikin,iz),jj=1,5),
     >     rat(5),rater(5),rat(1),rater(1),
     >     rat(4),rater(4)
          write(15,'(2i3,5f7.3,6f7.3)') ikin,iz,
     >     (avk(ikin,iz,jj)/avker(ikin,iz),jj=1,5),
     >     ratrho(5),rater(5),ratrho(1),rater(1),
     >     ratrho(4),rater(4)
         endif
        endif
       enddo
      enddo

c write super ratios in order of worst first
      if(usewide) then
       open(unit=7,file='ptmsratw.txt')
      else
       open(unit=7,file='ptmsrat.txt')
      endif
      do it=1,6
       if(it.eq.1.or.it.eq.4.or.it.eq.5) then
        do i=1,100000
         diffmax=0.
         ikinsv=0
         izsv=0
         do ikin=1,32
          do iz=1,20
           if(usedup(ikin,iz,it).eq.0.and.
     >      srater(ikin,iz,it).ne.0.) then
            diff = abs(srat(ikin,iz,it) - 1.) /
     >       srater(ikin,iz,it)
            if(diff.gt.diffmax) then
             diffmax = diff
             ikinsv = ikin
             izsv = iz
            endif
           endif
          enddo
         enddo ! ikin
         if(diffmax.eq.0.) goto 12
         usedup(ikinsv,izsv,it)=1
         write(7,'(i1,2i3,5f6.3,2f7.3)')
     >    it,ikinsv,izsv,
     >    (avk(ikinsv,izsv,jj)/avker(ikinsv,izsv),jj=1,5),
     >    srat(ikinsv,izsv,it),
     >    srater(ikinsv,izsv,it)
        enddo
 12     continue
       endif
      enddo

c table of average multiplicites (<ptcut) versus z
      close(unit=22)
      close(unit=23)
      if(usewide) then
       open(unit=22,file='ptmavmw.txt')
       open(unit=23,file='ptmavmkltw.txt')
      else
       open(unit=22,file='ptmavmn.txt')
       open(unit=23,file='ptmavmkltn.txt')
      endif
      do ikin=1,56,2
       x = (xkin(ikin) + xkin(ikin+1))/2.
       q2 = (q2kin(ikin) + q2kin(ikin+1))/2.
       w2 = am**2 + q2*(1/x-1)
       do iz=4,19
        do it=1,5
         avit(it)=0.
         avitrho(it)=0.
         aviter(it)=0.
         If(avmer(ikin ,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).gt.0.) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v3 = avmod(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avmod(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avit(it) = v1 / v2
          avmit(it) = v3 / v2
          aviter(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >         (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >          avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avitrho(it) = v1 / v2
         endif
        enddo ! it
        z = .05 * (iz-0.5)
        write(22,'(i2,3f6.3,12f8.4)') (ikin+1)/2,x,q2,z,
     >  avit(1),avitrho(1),aviter(1),
     >  avit(2),avitrho(2),aviter(2),
     >  avit(4),avitrho(4),aviter(4),
     >  avit(5),avitrho(5),aviter(5) 
c        write(22,'(20x,4(f8.4,16x))')
c     >   avmit(1),avmit(2),avmit(4),avmit(5)
        if(ikin.gt.0) then
         write(23,'(i2,3f6.3,8f8.4)') 
     >    (ikin+1)/2,x,q2,z,
     >    avit(1)/avmit(1),aviter(1)/avmit(1),
     >    avit(2)/avmit(2),aviter(2)/avmit(2),
     >    avit(4)/avmit(4),aviter(4)/avmit(4),
     >    avit(5)/avmit(5),aviter(5)/avmit(5)
        endif
       enddo ! iz
      enddo ! ikin

c plot multiplicity ratios over deuteron pi+ versus z  and compare
c to SHuo if available.
c also, get two ratio of FF and plot those two
      close(unit=21)
      close(unit=22)
      close(unit=23)
      open(unit=21,file='ptmavm.top')
      write(21,'(''set device postscript'')')
      open(unit=22,file='ptmavmrat.top')
      write(22,'(''set device postscript'')')
      open(unit=23,file='ptmavmratff.top')
      write(23,'(''set device postscript'')')
      open(unit=25,file='ptmavmratff4.top')
      write(25,'(''set device postscript'')')
      open(unit=31,file='ptmavmratff4no.top')
      write(31,'(''set device postscript'')')
      open(unit=24,file='ptmavmratffrm.top')
      write(24,'(''set device postscript'')')
      open(unit=26,file='ptmavmratffr.top')
      write(26,'(''set device postscript'')')
      open(unit=27,file='ptmavmratcsv.top')
      write(27,'(''set device postscript'')')
      open(unit=28,file='ptmavmratex.top')
      write(28,'(''set device postscript'')')
      close(unit=29)
      open(unit=29,file='ptmavmratcsv2.top')
      write(29,'(''set device postscript'')')
      close(unit=30)
      open(unit=30,file='ptmdiffrat.top')
      write(30,'(''set device postscript'')')
      ix=0
      iy=3
c only do cases with proton
      do ikinpp=1,8
       ikin = 2*ikinpp - 1
       if(ikinpp.eq.8) ikin=29
c change order to increasing W
       ikin = ikinw(ikinpp)
c      do ikin=1,32,2
       ix=ix+1
       if(ix.gt.4) then
        ix=1
        iy = iy-1
       endif
       xpltmin = 1.1 + 2.95*(ix-1)
       xpltmax = xpltmin + 2.95
       ypltmax = 9.5 - 4.0*(3-iy)
       ypltmin = ypltmax - 4.0
       x = (xkin(ikin) + xkin(ikin+1))/2.
       q2 = (q2kin(ikin) + q2kin(ikin+1))/2.
       w2 = am**2 + q2*(1/x-1)
c       write(21,2188) ix,iy,x,q2,sqrt(w2)
       do kk=1,12
        if(ix.eq.1) write(20+kk,'(''set labels left on'')')
        if(ix.ne.1) write(20+kk,'(''set labels left off'')')
        if(iy.eq.2) write(20+kk,'(''set labels bottom on'')')
        if(iy.ne.2) write(20+kk,'(''set labels bottom off'')')
        write(20+kk,2088) xpltmin,xpltmax,ypltmin,ypltmax,
     >   xpltmin+0.25, ypltmax-0.35,x,q2,sqrt(w2)
 2088   format(
     >   'set intensity 4'/
     >   'set color white'/
     >   'set window x ',2f8.3,' y ', 2f8.3/
     >   'title ',2f8.3,' size 1.6 ',1h',
     >   'x=',f4.2,'  Q223=',f3.1,' W=',f3.1,1h'/
     >   'case ',1h','         X X',1H'/
     >   'set symbol 9O size 1.0 ; set bar size 0.'/
     >   'set ticks size 0.04 ; set order x y dy sym')
       enddo ! kk

       if(ix.eq. 1.and.iy.eq.2) then
        write(21,8801)
 8801   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','zM(z)',1h')
        write(23,8803)
 8803   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z)',1h')
        write(24,8804)
 8804   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
c     >          1h','D(z) / D0DSS3(z)',1h'/
     >          1h','D(z) / D0JAM3(z)',1h'/
     >  'CASE ',1H','        X   X   ',1H')
        write(25,8805)
        write(31,8805)
 8805   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','zD(z)',1h')
        write(26,8806)
 8806   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >           1h','D(z) / D0u12+3(z)',1h'/
     >  'case  ',1h','        X XX X   ',1h')
        write(27,8807)
 8807   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z)',1h')
        write(28,8808)
 8808   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z) or excl factor',1h')
        write(29,8809)
 8809   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','CSV(x)',1h')
        write(30,8810)
 8810   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','Ratio',1h')
        endif

       write(21,2188)
 2188  format(
     > 'set limits x 0.25 0.75 y 0.0 0.49'/
     > 'plot axes')
       write(22,2288)
 2288  format(
     > 'set limits x 0.25 0.75 y 0.35 1.4'/
     > 'plot axes')
       write(23,2372)
 2372  format(
     > 'set limits x 0.25 0.75 y 0.0 0.57'/
     > 'plot axes')
       write(24,2472)
 2472  format(
     > 'set limits x 0.25 0.75 y 0.0 1.9'/
     > ' 0. 1. ; 1. 1. ; set pattern .03 .03 .03 .03 ; join pattern'/
     > 'plot axes')
       write(25,2572)
 2572  format(
     > 'set limits x 0.25 0.75 y  0.  0.57'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(31,3172)
 3172  format(
     > 'set limits x 0.25 0.75 y  -0.2  0.57'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(26,2562)
 2562  format(
     > 'set limits x 0.25 0.75 y -0.5 1.5'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(27,2762)
 2762  format(
     > 'set limits x 0.25 0.75 y 0. 0.59'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(28,2862)
 2862  format(
     > 'set limits x 0.25 0.75 y 0 1.5'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(29,2962)
 2962  format(
     > 'set limits x 0.25 0.75 y -0.6 0.59'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(30,3062)
 3062  format(
     > 'set limits x 0.25 0.75 y 0.3 1.25'/
c     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
c     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')

c my multiplicity results
       do it=1,5
       if(it.eq.1.or.it.eq.4.or.it.eq.5.or.it.eq.2) then
       if(it.eq.1) write(21,'(''set color blue'')')
       if(it.eq.2) write(21,'(''set color red'')')
       if(it.eq.4) write(21,'(''set color green'')')
       if(it.eq.5) write(21,'(''set color white'')')
       nplt = 0
       do iz=1,20
        z = .05 * (iz-0.5)
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.3 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.3) then
         v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >        avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         write(21,'(3f10.4,'' 9O'')')  
     >    0.05*(iz-0.5), z*v1,z*v2
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) then
            write(21,'(''set symbol size 1.0 ; plot'')')
            write(21,'(''set symbol size 0.8 ; plot'')')
            write(21,'(''set symbol size 0.6 ; plot'')')
            write(21,'(''set symbol size 0.4 ; plot'')')
            write(21,'(''set symbol size 0.2 ; plot'')')
       endif
! with rho subtraction
       nplt = 0
       do iz=1,20
        z = .05 * (iz-0.5)
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.3 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.3) then
         v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         write(21,'(3f10.4,'' 9O'')')  
     >    0.05*(iz-0.5)-0.01, z*v1,z*v2
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(21,'(''set symbol size 1.0 ; plot'')')
       write(21,'(''set intensity 2'')')

c plot various simc models
c now only plot dss (my fit too big) and one pdf
c (they give almost identical results)
c now plot dss with/without swap, jam (2, 3, 4)
       do iff=2,4
        do ipdf=2,2
         okplt = .true.
         if(iff.eq.3 .and. (it.eq.1 .or. it.eq.4)) okplt=.false.
         do iz=6,18
          z = 0.05*iz
          pt=0.
          mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >      sighadm,u,ub,d,db,u1,d1,
     >      s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          if(okplt) write(21,'(f8.3,f10.4,2f8.3,2i2)') z,z*dsigdz
         enddo
         if(okplt) then
          if(iff.eq.2.and.ipdf.eq.2) 
     >     write(21,'(''set pattern .05 .05 .05 .05'')')
          if(iff.eq.3.and.ipdf.eq.2) 
     >     write(21,'(''set pattern .015 .015 .015 .015'')')
          if(iff.le.3) write(21,'(''join pattern'')')
          if(iff.eq.4) write(21,'(''join'')')
         endif
        enddo ! pdf
       enddo ! iff
       endif ! it
       enddo ! it

! my ratio to D+ results
       do it=1,5
       do iz=1,20
        rver(it,iz)=0.
       enddo
       if(it.eq.1.or.it.eq.4.or.it.eq.5) then
       if(it.eq.1) write(22,'(''set color blue'')')
       if(it.eq.4) write(22,'(''set color green'')')
       if(it.eq.5) write(22,'(''set color white'')')
       nplt = 0
       do iz=1,20
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz, 2).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin,iz, 2).lt.0.2.and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz, 2).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz, 2).lt.0.2) then 
         v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >        avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v3 = avm(ikin,  iz, 2)/avmer(ikin,  iz, 2)**2 +
     >        avm(ikin+1,iz, 2)/avmer(ikin+1,iz, 2)**2 
         v4 =              1./avmer(ikin,  iz, 2)**2 +
     >                     1./avmer(ikin+1,iz, 2)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         v3 = v3 / v4
         v4 = sqrt(1./v4)
         v5 = v1 / v3
         v6 = sqrt( (v2/v3)**2 + (v5/v3*v4)**2)
         rv(it,iz) = v1
         rver(it,iz) = v2
         rv(2,iz) = v3
         rver(2,iz) = v4
         write(22,'(3f10.4,'' 9O'')')  0.05*(iz-0.5), v5, v6
c do again with rho subtraction
         v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v3 = (avm(ikin,  iz, 2)-avmrho(ikin,  iz,2))/
     >          avmer(ikin,  iz, 2)**2 +
     >        (avm(ikin+1,iz, 2)-avmrho(ikin+1,iz,2))/
     >          avmer(ikin+1,iz, 2)**2 
         v4 =              1./avmer(ikin,  iz, 2)**2 +
     >                     1./avmer(ikin+1,iz, 2)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         v3 = v3 / v4
         v4 = sqrt(1./v4)
         v5 = v1 / v3
         v6 = sqrt( (v2/v3)**2 + (v5/v3*v4)**2)
         write(22,'(3f10.4,'' 8O'')')  0.05*(iz-0.5)-0.01, v5, v6
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(22,'(''plot'')')
c plot model
       nplt = 0
       do iz=4,18
        if(avms(ikin,iz,it).gt.0. .and.
     >     avms(ikin,iz,2).gt.0.) then
         write(22,'(3f10.4)')  0.05*(iz-0.5),
     >    avms(ikin,iz,it)/avms(ikin,iz,2)
         nplt = nplt + 1
        endif
       enddo
       if(nplt.gt.0) write(22,'(''join'')')
       write(22,'(''set order x y dum dum'')')
c plot various simc models
       do iff=1,2
        do ipdf=1,2
         do iz=4,18
          z = 0.05*iz
          pt=0.
          mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >     sighadp,u,ub,d,db,u1,d1,
     >     s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >      sighadm,u,ub,d,db,u1,d1,
     >      s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          write(22,'(f8.3,3f10.4,2f8.3,2i2)') z,sighadm/sighadp,
     >      sighadm,sighadp,u,d,iff,ipdf
         enddo
         if(iff.eq.1.and.ipdf.eq.1) 
     >    write(22,'(''set pattern .01 .01 .01 .01'')')
         if(iff.eq.1.and.ipdf.eq.2) 
     >    write(22,'(''set pattern .01 .05 .01 .05'')')
         if(iff.eq.2.and.ipdf.eq.1) 
     >    write(22,'(''set pattern .02 .02 .02 .02'')')
         if(iff.eq.2.and.ipdf.eq.2) 
     >    write(22,'(''set pattern .05 .05 .05 .05'')')
         write(22,'(''join pattern'')')
        enddo
       enddo
c plot lund predictions
c values are from low x, but plot all of them
       if(xkin(ikin).lt.0.60) then
        write(22,'(''set order x y '')')
        call FNP_NMC(X,Q2,rnp)
        ikinp = max(1,min(8,int((x-0.275)/0.05)+1))
c changed to use lowest x bin
        ikinp = 1
        do iz=3,18
         v2 = (lund(ikinp,iz,1)  + rnp * 
     >         lund(ikinp,iz,3)) / (1 + rnp)
         if(it.eq.1) v1 = lund(ikinp,iz,1)
         if(it.eq.4) v1 = lund(ikinp,iz,2)
         if(it.eq.5) v1  =(lund(ikinp,iz,2)  + rnp * 
     >         lund(ikinp,iz,4)) / (1 + rnp)
         mvlund(1) = lund(ikinp,iz,1)
         mvlund(2) = v2
         mvlund(3) = lund(ikinp,iz,2)
         mvlund(4) = (lund(ikinp,iz,2)  + rnp * 
     >         lund(ikinp,iz,4)) / (1 + rnp)
         write(22,'(2f8.4)') 0.05*(iz-0.5),v1/v2
        enddo
        write(22,'(''set pattern .09 .09 .09 .09'')')
        write(22,'(''join pattern'')')
c        write(22,'(''join ; set color white'')')
        write(22,'(''set order x y dy sym'')')
       endif
       write(6,'(''xx3, it, ikin'',3i3)') it,ikin
      endif ! it
      write(6,'(''xx3, it, ikin'',3i3)') it,ikin
      enddo ! it

      write(6,'(''xx3, it, ikin'',3i3)') it,ikin
c fit 3 D(z) and then 4 D(z) for this kin
c either swapping the two favoreds or not
      do iz=3,19
       do j=1,4
        do iswap=0,1
         avffsver(iz,j,iswap)= 0.
         avffsv4er(iz,j,iswap)= 0.
        enddo
        ffrater(iz,j)=0.
        ffrat4er(iz,j)=0.
        avffsvcsver(iz,j)= 0.
       enddo
       if(rver(1,iz).gt.0. .and.
     >    rver(2,iz).gt.0. .and.
     >    rver(4,iz).gt.0. .and.
     >    rver(5,iz).gt.0.) then
        z = 0.05 * (iz-0.5)
! fit to data with NO rho subratcion
        mv(1) = rv(1,iz) * z
        mv(2) = rv(2,iz) * z
        mv(3) = rv(4,iz) * z
        mv(4) = rv(5,iz) * z
        mver(1) = rver(1,iz) * z
        mver(2) = rver(2,iz) * z
        mver(3) = rver(4,iz) * z
        mver(4) = rver(5,iz) * z
        it = 1
        pt = 0.
        phicm = 0.
        ipdf = 1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,muf,mubf,mdf,mdbf,u1,d1,
     >   msf,msbf,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        write(6,'(''rfffit'',3i3,8f7.3)') ikin,iz,ifit,
     >   (mv(it),mver(it),it=1,4)
c fit with two fav. one unfavored
        do iswap=0,1
         nparam=3
         npar = 3
         p(1) = fu1 ! fav u->pi+
         p(2) = fub ! unfav u->pi+
         p(3) = fdb
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
         arglis(1)=0
         ncallff=0
         call mnexcm(fffit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          avffsv  (iz,j,iswap)= coef(j)
          avffsver(iz,j,iswap)= std(j)
c          write(6,'(''rfffit'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
c     >     j,p(j),coef(j),std(j)
         enddo ! j
         write(6,'(''3param'',2i3,i2,5f7.3)') ikin,iz,
     >    iswap,chi2ff,
     >    (mv(1) - mv(3))/(mv(2) - mv(4)),
     >    (mvf(1) - mvf(3))/(mvf(2) - mvf(4)),
     >    (mv(1) + mv(3))/(mv(2) + mv(4)),
     >    (mvf(1) + mvf(3))/(mvf(2) + mvf(4))
         do it=1,4
          ffrat(iz,it) = mv(it)/mvf(it)
          ffrater(iz,it) = mver(it)/mvf(it)
         enddo
        enddo ! iswap
c fit with two fav. two unfavor21,ed
c using swap or not 
        do iswap=0,1
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = fub ! unfav u-> pi-
        p(4) = fdb ! fav d-> pi-+
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffit4_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsv4  (iz,j,iswap)= coef(j)
         avffsv4er(iz,j,iswap)= std(j)
         write(6,'(''rfffit4'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        do it=1,4
         ffrat4(iz,it) = mv(it)/mvf(it)
         ffrat4er(iz,it) = mver(it)/mvf(it)
         write(6,'(''ffrat4'',3i3,2f8.3)') iswap,iz,it,ffrat4(iz,it),
     >    ffrat4er(iz,it)
        enddo
        enddo ! iswap
c fit with one favored, one unfavored, and two cv params
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = 0. ! delta d / d
        p(4) = 0. ! delta u / u
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        iswap = 1.
        call mnexcm(fffitcsv_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsvcsv  (iz,j)= coef(j)
         avffsvcsver(iz,j)= std(j)
         write(6,'(''rfffitcsv'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        write(6,'(''csv chk'',f8.3)') 
     >   (4.*mdf*coef(3) + muf*coef(4)) /
     >   (4.*mdf + muf)
c fit with one favored, one unfavored, and two excl. 
c enhancement factors
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = 1. ! scales up pi+ from p
        p(4) = 1. ! scales up pi- from n
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffitex_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsvex  (iz,j)= coef(j)
         avffsvexer(iz,j)= std(j)
         write(6,'(''rfffitex'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
       endif ! check on rver>0
      enddo ! iz
c plot the three D(z) results for this kin
        do j=1,3
         if(j.eq.2) write(23,'(''set color blue'')')
         if(j.eq.3) write(23,'(''set color red'')')
         do iswap=0,1
          nplt = 0
          do iz=1,20
           if(avffsver(iz,j,iswap).ne.0. .and.
     >        avffsver(iz,j,iswap).lt.0.1) then
            nplt = nplt + 1
            fact=1.
            write(23,136) 
     >       0.05*(iz-0.5) + 0.005*(j)-0.02*iswap,
     >       fact*avffsv  (iz,j,iswap),
     >       fact*avffsver(iz,j,iswap)
           endif
          enddo
          if(nplt.gt.0) then
           write(23,'(''plot'')')
           if(iswap.eq.1) then
            write(23,'(''set symbol size 0.8 ; plot'')')
            write(23,'(''set symbol size 0.6 ; plot'')')
            write(23,'(''set symbol size 0.4 ; plot'')')
            write(23,'(''set symbol size 0.2 ; plot'')')
           endif
          endif
         enddo ! iswap
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(23,'(2f8.4)') z, fu1
          if(j.eq.2) write(23,'(2f8.4)') z, fd1
          if(j.eq.3) write(23,'(2f8.4)') z, fdb
         enddo
         write(23,'(''join'')')
        enddo ! j
c plot the four D(z) results for this kin
c also plot ratios to u+
c on unit-27, 3 D(z) and one CSV params on unit 29
        do j=1,4
         if(j.eq.2) write(25,'(''set color blue'')')
         if(j.eq.3) write(25,'(''set color red'')')
         if(j.eq.4) write(25,'(''set color green'')')
         if(j.eq.2) write(31,'(''set color blue'')')
         if(j.eq.3) write(31,'(''set color red'')')
         if(j.eq.4) write(31,'(''set color green'')')
         if(j.eq.2) write(26,'(''set color blue'')')
         if(j.eq.3) write(26,'(''set color red'')')
         if(j.eq.4) write(26,'(''set color white'')')
         if(j.eq.2) write(27,'(''set color blue'')')
         if(j.eq.3) write(27,'(''set color red'')')
         if(j.eq.4) write(29,'(''set color magenta'')')
         if(j.eq.2) write(28,'(''set color blue'')')
         if(j.eq.3) write(28,'(''set color red'')')
         if(j.eq.4) write(28,'(''set color green'')')
! ff4 noswap on unit 31
         nplt = 0
         do iz=1,20
          if(avffsv4er(iz,j,0).ne.0. .and.
     >       avffsv4er(iz,j,0).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,0),
     >        fact*avffsv4er(iz,j,0)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''set symbol size 1.0 ; plot'')')
          write(31,'(''set symbol size 0.8 ; plot'')')
          write(31,'(''set symbol size 0.6 ; plot'')')
          write(31,'(''set symbol size 0.4 ; plot'')')
          write(31,'(''set symbol size 0.2 ; plot'')')
         endif
         do iz=1,20
          if(avffsv4er(iz,j,1).ne.0. .and.
     >       avffsv4er(iz,j,1).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           write(25,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,1),
     >        fact*avffsv4er(iz,j,1)
           if(j.gt.1) write(26,136) 
     >       0.05*(iz-0.5) + 0.005*(j-2),
     >       avffsv4  (iz,j,1)/avffsv4(iz,1,1),
     >       avffsv4er(iz,j,1)/avffsv4(iz,1,1)
          endif
         enddo
         if(nplt.gt.0) then
          write(25,'(''plot'')')
          write(25,'(''set symbol size 1.0 ; plot'')')
          write(25,'(''set symbol size 0.8 ; plot'')')
          write(25,'(''set symbol size 0.6 ; plot'')')
          write(25,'(''set symbol size 0.4 ; plot'')')
          write(25,'(''set symbol size 0.2 ; plot'')')
          write(26,'(''plot'')')
          write(26,'(''set symbol size 0.8 ; plot'')')
          write(26,'(''set symbol size 0.6 ; plot'')')
          write(26,'(''set symbol size 0.4 ; plot'')')
          write(26,'(''set symbol size 0.2 ; plot'')')
         endif
         if(j.le.3) ioff=0
         if(j.eq.4) ioff=2
         nplt = 0
         do iz=1,20
          if(avffsvcsver(iz,j).ne.0. .and.
     >       avffsvcsver(iz,j).lt.0.4) then
           nplt = nplt + 1
           fact=1.
           write(27+ioff,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >       fact*avffsvcsv  (iz,j),
     >       fact*avffsvcsver(iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(27+ioff,'(''plot'')')
          write(27+ioff,'(''set symbol size 0.8 ; plot'')')
          write(27+ioff,'(''set symbol size 0.6 ; plot'')')
          write(27+ioff,'(''set symbol size 0.4 ; plot'')')
          write(27+ioff,'(''set symbol size 0.2 ; plot'')')
         endif
c limits on csv(x) from MRST
         if(j.eq.4) then
          v1 = 0.8 * (1. - x)**4 / sqrt(x) *
     >     (x - 0.0909)
          write(29,'(''0. '',f8.4,'' ; 1. '',f8.4,
     >     '' ; join'')') v1,v1
          v1 = -0.65 * (1. - x)**4 / sqrt(x) *
     >     (x - 0.0909)
          write(29,'(''0. '',f8.4,'' ; 1. '',f8.4,
     >     '' ; join'')') v1,v1
         endif ! j.eq.4
         nplt = 0
         do iz=1,20
          if(avffsvex er(iz,j).ne.0. .and.
     >       avffsvex er(iz,j).lt.0.2) then
           nplt = nplt + 1
           fact=1.
           write(28,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsvex   (iz,j),
     >        fact*avffsvex er(iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(28,'(''plot'')')
          write(28,'(''set symbol size 0.8 ; plot'')')
          write(28,'(''set symbol size 0.6 ; plot'')')
          write(28,'(''set symbol size 0.4 ; plot'')')
          write(28,'(''set symbol size 0.2 ; plot'')')
         endif
c plot my fit to D(z)
c xxx need to update with new fit
         do iz=6,19
          z = 0.05 * iz
          iff=1
          ipdf=1
          it=1
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
c          if(j.eq.1) write(25,'(2f8.4)') z, u1
c          if(j.eq.2) write(25,'(2f8.4)') z, d1
         enddo
c         if(j.le.2) write(25,'(''join dash'')')
c plot DSS
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(25,'(2f8.4)') z, fu1
          if(j.eq.2) write(25,'(2f8.4)') z, fd1
          if(j.eq.3) write(25,'(2f8.4)') z, fub
          if(j.eq.4) write(25,'(2f8.4)') z, fdb
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fub
          if(j.eq.4) write(31,'(2f8.4)') z, fdb
          if(j.eq.2) write(26,'(2f8.4)') z, fd1/fu1
          if(j.eq.3) write(26,'(2f8.4)') z, fub/fu1
          if(j.eq.4) write(26,'(2f8.4)') z, fdb/fu1
          if(j.eq.1) write(27,'(2f8.4)') z, fu1
          if(j.eq.2) write(27,'(2f8.4)') z, fd1
          if(j.eq.1) write(28,'(2f8.4)') z, fu1
          if(j.eq.2) write(28,'(2f8.4)') z, fd1
         enddo
         write(25,'(''join dash'')')
         write(31,'(''join dash'')')
         if(j.gt.1) write(26,'(''join dash'')')
         if(j.le.2) write(27,'(''join dash'')')
         if(j.le.2) write(28,'(''join dash'')')
c plot JAM
         do iz=6,19
          z = 0.05 * iz
          call jamff(z,q2kin(ikin),fd1,fu1)
          if(j.eq.1) write(25,'(2f8.4)') z, fu1
          if(j.eq.2) write(25,'(2f8.4)') z, fd1
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.1) write(27,'(2f8.4)') z, fu1
          if(j.eq.2) write(27,'(2f8.4)') z, fd1
          if(j.eq.2) write(26,'(2f8.4)') z, fd1/fu1
         enddo
         if(j.le.2) write(25,'(''join'')')
         if(j.le.2) write(27,'(''join'')')
         if(j.eq.2) write(26,'(''join'')')
c plot LUND ratios
c xxx need to get correct values from fitlund.f
         if(j.gt.1) then
          do iz=7,14
           z = 0.05 * iz
c           write(26,'(2f8.4)') 
c     >      z, avffsv4lund(iz,j)/
c     >      avffsv4lund(iz,1)
           write(6,'(''chklund''2i3,3f8.4)') ikin,iz, 
     >      z, avffsv4lund(iz,j),
     >      avffsv4lund(iz,1)
          enddo
c          write(26,'(''join dotdash'')')
         endif
        enddo ! j
c plot the D(z) / D(z)_fit results for this kin
c changed to be D(z) / FDSS (using averge of the
c two favored ones from DSS)
        do j=1,4
         if(j.eq.2) write(24,'(''set color blue'')')
         if(j.eq.3) write(24,'(''set color red'')')
         if(j.eq.4) write(24,'(''set color green'')')
         nplt=0
         do iz=1,20
          if(avffsv4er(iz,j,1).ne.0. .and.
     >       avffsv4er(iz,j,1).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           z = 0.05 * (iz-0.5)
           call fDSS (1,1,0, Z, Q2, 
     >       fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
           call jamff(z,q2kin(ikin),fd1j,fu1j)
           if(iz.eq.10 .and. j.eq.1) write(6,
     >      '(''fdsschk4'',6f8.3)') fu1,fdb,fu1j,fub,fd1,fd1j
           if(j.eq.1.or.j.eq.4) fact=2./(fu1 + fdb)
           if(j.eq.2.or.j.eq.3) fact = 2./(fub + fd1)
c ratios to JAM
           if(j.eq.1.or.j.eq.4) fact=1./fu1j
           if(j.eq.2.or.j.eq.3) fact = 1./fd1j
           write(24,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,1),
     >        fact*avffsv4er(iz,j,1)
          endif
         enddo
         if(nplt.gt.0) then
          write(24,'(''plot'')')
          write(24,'(''set symbol size 0.8 ; plot'')')
          write(24,'(''set symbol size 0.6 ; plot'')')
          write(24,'(''set symbol size 0.4 ; plot'')')
          write(24,'(''set symbol size 0.2 ; plot'')')
         endif
        enddo ! j
c sum,diff ratios on unit 30 for this ikin
       nplt = 0
       do iz=1,20
        do it=1,6
         avmk(it)=0.
         avmkr(it)=0.
         If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avmk(it) = v1 / v2
          avmker(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 

          avmkr(it) = v1 / v2
          avmkrer(it) = sqrt(1./v2)
          write(6,'(''diffdbg'',i3,i3,i2,6f7.3)') ikin,iz,it,
     >     avm(ikin,  iz,it),avmer(ikin,  iz,it),
     >     avm(ikin+1,iz,it),avmer(ikin+1,iz,it),avmk(it),avmker(it)
         endif
        enddo ! it
        if(avmk(1).ne.0. .and. avmk(2).ne.0. .and.
     >     avmk(4).ne.0. .and. avmk(5).ne.0.) then
c diff ratio without rho subtraction (d+ - d-) / (p+ - p-)
         v1 = 1./ (avmk(1) - avmk(4)) *
     >        (avmk(2) - avmk(5))
         v2 = 1./(avmk(1) + avmker(1) - avmk(4)) *
     >        (avmk(2) - avmk(5))
         v3 = 1./(avmk(1) - avmk(4) - avmker(4)) *
     >        (avmk(2) - avmk(5))
         v4 = 1./(avmk(1) - avmk(4)) *
     >        (avmk(2) + avmker(2) - avmk(5))
         v5 = 1./(avmk(1) - avmk(4)) *
     >        (avmk(2) - avmk(5) - avmker(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
         write(30,'(3f10.4,'' 3O'')')  
c xxx    >    0.05*(iz-0.5)+0.006, v1,v6
     >    0.05*(iz-0.5)-0.006, v1,v6
! sum ratio without rho subtraction. (d+ + d-) / (p+ + p-)
         v1 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmk(5))
         v2 = 1./(avmk(1) + avmker(1) + avmk(4)) *
     >        (avmk(2) + avmk(5))
         v3 = 1./(avmk(1) + avmk(4) + avmker(4)) *
     >        (avmk(2) + avmk(5))
         v4 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmker(2) + avmk(5))
         v5 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmk(5) + avmker(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
         write(30,'(3f10.4,'' 9O'')')  
cxxx     >    0.05*(iz-0.5)+0.006, v1,v6
     >    0.05*(iz-0.5)-0.006, v1,v6
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) then
          write(30,'(''set symbol size 1.0 ; plot'')')
          write(30,'(''set symbol size 0.8 ; plot'')')
          write(30,'(''set symbol size 0.6 ; plot'')')
          write(30,'(''set symbol size 0.4 ; plot'')')
          write(30,'(''set symbol size 0.2 ; plot'')')
       endif
c diff ratio with rho subtraction
       nplt = 0
       do iz=1,20
        do it=1,6
         avmk(it)=0.
         avmkr(it)=0.
         If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avmk(it) = v1 / v2
          avmker(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 

          avmkr(it) = v1 / v2
          avmkrer(it) = sqrt(1./v2)
          write(6,'(''diffdbg'',i3,i3,i2,6f7.3)') ikin,iz,it,
     >     avm(ikin,  iz,it),avmer(ikin,  iz,it),
     >     avm(ikin+1,iz,it),avmer(ikin+1,iz,it),avmk(it),avmker(it)
         endif
        enddo ! it
        if(avmk(1).ne.0. .and. avmk(2).ne.0. .and.
     >     avmk(4).ne.0. .and. avmk(5).ne.0.) then
         v1 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5))
         v2 = 1./(avmkr(1) + avmkrer(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5))
         v3 = 1./(avmkr(1) - avmkr(4) - avmkrer(4)) *
     >        (avmkr(2) - avmkr(5))
         v4 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) + avmkrer(2) - avmkr(5))
         v5 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5) - avmkrer(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
c         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
c         write(30,'(3f10.4,'' 3O'')')  0.05*(iz-0.5)-0.006, v1,v6
         v1 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5))
         v2 = 1./(avmkr(1) + avmkrer(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5))
         v3 = 1./(avmkr(1) + avmkr(4) + avmkrer(4)) *
     >        (avmkr(2) + avmkr(5))
         v4 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkrer(2) + avmkr(5))
         v5 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5) + avmkrer(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
c         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
c         write(30,'(3f10.4,'' 9O'')')  0.05*(iz-0.5)-0.006, v1,v6
c         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(30,'(''set symbol size 1.0 ; plot'')')
       pt=0.
       mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
       call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >     sighadp,u,ub,d,db,u1,d1,
     >     s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,2)
       v1 = (d - db) / (u - ub)
       v2 = (4.*(u-ub) - (d-db)) / 3. / (u - ub + d - db)
       v2 = v2 * 5. * (u + ub + d + db) /
     >           (4.*(u + ub) + (d + db))
c v2 = (4u - d) / 3 / (u+d) * 5 * (u+d) / (4u + d)
c v2 = (5/3) (4u-d) / (4u + d)
       v2 = 1./v2
c (3/5) (4u + d) / (4u -d)      
        write(30,'(''0.2 '',f10.4,'' ; 1.0 '',f10.4)') 
     >   v2,v2
       write(30,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
c diff ratio with d - db = 0.
c       write(22,'(''set pattern .02 .02 .02 .02'')')
c       write(22,'(''0.2 1.66 ; 1.0 1.66 ; join 1'')') 
c sum rattio
       write(30,'(''0.2 1.00 ; 1.0 1.00'')') 
       write(30,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
c simc model using 2 Frag Fun.
       do iff=2,4,2
       if(iff.eq.2) write(30,'(''set pattern .05 .05 .05 .05'')')
       if(iff.eq.4) write(30,'(''set pattern .02 .02 .02 .02'')')
       if(iff.eq.3) write(30,'(''set pattern .01 .01 .01  .01 '')')
       do iz=4,18
        z =0.05*iz
        ipdf=1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,1,
     >   sighad1,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >   sighad2,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,4,
     >   sighad4,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,5,
     >   sighad5,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
         write(30,'(2f8.4)') z,1./(sighad1-sighad4)*
     >    (sighad2 - sighad5)
         rsv(iz) = 1./(sighad1 + sighad4)*
     >    (sighad2 + sighad5)
        enddo
        if(iff.eq.4) then
         write(30,'(''join '')')
        else
         write(30,'(''join pattern'')')
        endif
        do iz=3,18
         z =0.05*iz
         write(30,'(2f8.4)') z,rsv(iz)
        enddo
        if(iff.eq.4) then
         write(30,'(''join '')')
        else
         write(30,'(''join pattern'')')
        endif
       enddo ! iff
c       write(22,'(8f7.2)') x,q2,v1,v2
      enddo ! ikin

c debug simcmodel
        do iff=1,2
         do ipdf=1,2
          do j=1,5
           do iz=1,3
             z= 0.25 * iz
             x = 0.2 + 0.1 * j
             q2=3.
             pt=0.
             phicm=0.
             it=1
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v1 = (4 * u * u1 + d * d1)/
     >        (4.*u + d) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v1,v1/dsigdz,
c     >        sighad1/dsigdz,(ub + db + s + sb)/(u+d)
             it=2
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v2 = (4 * (u+d)* u1 + (d+u) * d1)/
     >        (4.*(u+d) + (d+u)) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v2,v2/dsigdz,
c     >        sighad2/dsigdz,(ub + db + s + sb)/(u+d)
             it=4
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v3 = (4 * u * u1 + d * d1)/
     >        (4.*u + d) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v3,v3/dsigdz,
c     >        sighad3/dsigdz,(ub + db + s + sb)/(u+d)
             it=5
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v4 = (4 * (u+d)* u1 + (d+u) * d1)/
     >        (4.*(u+d) + (d+u)) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v4,v4/dsigdz,
c     >        sighad4/dsigdz,(ub + db + s + sb)/(u+d)
             v5 = (d - db) / (u - ub)
             write(6,'(''dbgs'',2i2,f6.2,10f7.3)') iff,ipdf,x,
     >        sighad1, sighad2, sighad3, sighad4,
     >        (sighad2 - sighad4) / (sighad1 - sighad3),
     >        (sighad2 + sighad4) / (sighad1 + sighad3),
     >        0.6 * (1 + v5/4.)/(1 - v5/4.),v5,
     >        (v2 - v4) / (v1 - v3)
            enddo ! iz
           enddo ! j (x)
          enddo ! iff
         enddo ! ipdf

       do j=3,8
        z=0.1*j
        q2=3.
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
         write(6,'(''FDSS'',f6.2,6F7.4)') z,
     >      fU1, fUB, fD1, fDB, fS1, fSB
       enddo

c       goto 999
ccccccccccccccccc
c look at shuo's results. OLD ones
      open(unit=8,file='shuoratlarge.txt')
      do i=1,386
       read(8,*) q2n,q2,q2er,xn,x,xer,zn,z,zer,j,
     >   rd,rder,rdrho
       pt=0.1
       phicm=0.
       it=2
       mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighadp,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
       it=5
       write(6,'(''sh mx'',5f7.3)') q2,x,z,mmpi2
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighadm,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c for test
c       sighadp = sighadp * 0.95
       rdm = (4*sighadm - sighadp) /
     >       (sighadp - sighadm)
       diff = (rd/rdm - 1.0) /
     >       (rder/rdm)
       diff = max(-5.,min(4.999,diff))
       k = int((diff + 5.)/0.25) +1
       rdh(k,1) = rdh(k,1)+1
       diff = (rdrho/rdm - 1.0) /
     >       (rder/rdm)
       diff = max(-5.,min(4.999,diff))
       k = int((diff + 5.)/0.25) +1
       rdh(k,2) = rdh(k,2)+1
        write(6,'(''shuo'',3f6.3,6f7.3)') q2,x,z,
     >  rd/rdm,rder/rdm,rd,rder,rdrho,rdm
c average diff versus various variables
       diff = (rdrho/rdm - 1.0) /
     >       (rder/rdm)
       k=max(1,min(10,int((q2-3.)/3.5*10)+1))
       sq2(k,1) = sq2(k,1) + diff
       sq2(k,2) = sq2(k,2) + 1.
       k=max(1,min(10,int((x-0.25)/0.35*10)+1))
       sx(k,1) = sx(k,1) + diff
       sx(k,2) = sx(k,2) + 1.
       k=max(1,min(10,int((z-0.3)/0.4*10)+1))
       sz(k,1) = sz(k,1) + diff
       sz(k,2) = sz(k,2) + 1.   
       k=max(1,min(10,int((mmpi2-2.)/5*10)+1))
       smx(k,1) = smx(k,1) + diff
       smx(k,2) = smx(k,2) + 1.   
      enddo

      do k=1,10
       write(6,'(''sh kin'',4(f6.1,f5.0))') 
     >  sq2(k,1)/max(1.,sq2(k,2)),sq2(k,2),
     >  sx(k,1)/max(1.,sx(k,2)),sx(k,2),
     >  sz(k,1)/max(1.,sz(k,2)),sz(k,2),
     >  smx(k,1)/max(1.,smx(k,2)),smx(k,2)
      enddo
      do k=1,10
       write(6,'(''pb'',5(f6.1,f5.0))') 
     >  pq2(k,1)/max(1.,pq2(k,2)),pq2(k,2),
     >  px(k,1)/max(1.,px(k,2)),px(k,2),
     >  pz(k,1)/max(1.,pz(k,2)),pz(k,2),
     >  pw(k,1)/max(1.,pw(k,2)),pw(k,2),
     >  pmx(k,1)/max(1.,pmx(k,2)),pmx(k,2) 
      enddo
      do k=1,10
       write(6,'(''pb'',e12.4,f7.0)') 
     >  pq2(k,1),pq2(k,2)
      enddo

      open(unit=22,file='ptmsrat.top')
      write(22,'(''set device postscript'')')
      do i=1,2
       do j=1,3
        if(j.eq.1) it=1
        if(j.eq.2) it=4
        if(j.eq.3) it=5
        write(22,222) j,i,it,i
 222    format('set window x ',i1,' of 3 y ',i1,' of 2'/
     >   'title top ',1h','it,i=',2i3,1h')
        do k=1,40
         if(i.eq.1) write(22,'(f7.2,i6)') 
     >    -7. + 14./40.*(k-0.5),
     >    srath(k,it)
         if(i.eq.2) write(22,'(f7.2,i6)') 
     >    -7. + 14./40.*(k-0.5),
     >    sratrhoh(k,it)
        enddo
        write(22,'(''hist'')')
       enddo
      enddo
      write(22,'(''new frame'')')
      do i=1,2
       write(22,223) i
 223   format('set window x ',i1,' of 2 y 1 of 1')
       do k=1,40
        write(22,'(f7.2,i6)') 
     >    -5. + 0.25*(k-0.5),rdh(k,i)
       enddo
       write(22,'(''hist'')')
      enddo

      close(unit=22)

      write(6,'(''bnpt='',i6)') bnpt
c add in inclusiev
      bnpt0 = bnpt
      open(unit=44,file='ptm.incltxt')
      do i=1,100
       read(44,'(4x,6f10.4)') q2,w2,f1p,f2p,f1d,f2d
       xx = (w2 - am**2) / q2 
       x = 1./(xx + 1.)
       w2chk = am**2 + q2 * (1/x -1)
       z=0.4
       pt=0.
       phicm=0.
       mmpi2 = 4.
       it=1
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
       if(x.lt.0.7 .and. q2.gt.2.5) then
        bnpt = bnpt + 1
        bq2v(bnpt) = q2
        bxv(bnpt) = x
        buv(bnpt) = u
        bubv(bnpt) = ub
        bdv(bnpt) = d
        bdbv(bnpt) = db
        bsv(bnpt) = s
        bsbv(bnpt) = sb
c took out f1f2in21
c       byv(bnpt) = f2d / f2p
c use this instad
        call FNP_NMC(X,Q2,ratt)
        byv(bnpt) = 1. + ratt
        byerv(bnpt) = byv(bnpt) * 0.03
        bkin(bnpt) = 0 ! special code
        delu=0.
        deld=0.
        sum_sqp = qu**2*(u+ub) + qd**2*(d+db) + 
     >     qs**2*(s+sb)
        sum_sqn = qu**2*(d + deld + db) + 
     >         qd**2*(u + delu + ub) + qs**2*(s+sb)
        write(6,'(''f1f2'',3f5.2,5f7.3)') x,q2,w2,
     >   f2p, f2d, f2d/f2p, 1.+sum_sqn/sum_sqp, 1.+ratt
       endif
      enddo

c take out incl
      bnpt = bnpt0

c       goto 99

      if(usewide) then
       open(unit=22,file='ptmfitsw.top')
       open(unit=23,file='ptmfitpw.top')
      else
       open(unit=22,file='ptmfits.top')
       open(unit=23,file='ptmfitp.top')
      endif
      write(22,'(''set device postscript'')')
      write(23,'(''set device postscript'')')

c start of huge loop over ifit
      do ifit=1,3
      rhofact = 0.0
      rhofactp = 0.0
      if(ifit.eq.2) rhofact = 1.
c      if(ifit.eq.3) rhofact = 10.
      if(ifit.eq.2) rhofactp = 1.
c      if(ifit.eq.3) rhofactp = 0.
      ipdf = 1
      if(ifit.eq.3) ipdf = 2
      write(6,'(''using ipdf, rhofact = '',2i2,f5.1)') 
     >  ifit,ipdf,rhofact

      do j=1,19
        write(pname(j),'(''P'',i1)') j
      enddo

c fit each kin and target individually
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo
      lastframeused = .false.
      do ikinfit=1,32,2
       do iz=6,18
        if(lastframeused.and.ifit.eq.1) write(22,'(''new frame'')')
        if(lastframeused.and.ifit.eq.1) write(23,'(''new frame'')')
        lastframeused = .false.
        do itt=1,4
         itfit=itt
         if(itt.gt.2) itfit = itt + 1
         zfit = 0.05 * (iz - 0.5)
c         nparam=4
c         npar = 4
         nparam=3
         npar = 3
         p(1) = 0.1
         p(2) = 0.2
         p(3) = 0.0
         p(4) = 0.
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         ncallsfit=0
         call sfit_fcn(npar,grad,chi2,p,ierflg,futil)
         if(npts.gt.40) then
          write(6,'(''initial chi2'',2i4,f6.2,i6,f8.2)') 
     >     ikinfit,itfit,zfit,npts,
     >     chi2/float(npts)
          arglis(1)=0
          call mnexcm(sfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
          do j=1,nparam
           call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
           write(6,'(''sf '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
           scsv((ikinfit+1)/2,itt,iz,j,ifit)= coef(j)
           scsver((ikinfit+1)/2,itt,iz,j,ifit)= std(j)
          enddo
          call sfit_fcn(npar,grad,chi2,coef,ierflg,futil)
          write(6,'(''final chi2'',2i4,f6.2,i6,f8.2)') 
     >     ikinfit,itfit,zfit,npts,
     >     chi2/float(npts)
          write(16,'(i2,i2,f5.2,i4,9f7.3)') ikinfit,itt,
     >     zfit,npts,chi2/float(npts),
     >     (coef(j),std(j),j=1,4)
          if((ikinfit.eq.1 .or. ikinfit.eq.5) .and.
     >     ifit.eq.1.and.iz.gt.4 .and. iz.lt.18) then
           ix=1
           iy=1
           if(itt.eq.2.or.itt.eq.4) ix=2
           if(itt.eq.1.or.itt.eq.2) iy=2
           pt2 = 0.
           ymax = 1.5 * (coef(1)/zfit  ) / coef(2) * 
     >        exp(-pt2/coef(2)) 
           pt2 = 0.7**2
           ymin = 0.8 * (coef(1)/zfit  ) / coef(2) * 
     >        exp(-pt2/coef(2)) 
           x = (xkin(ikinfit) + xkin(ikinfit+1))/2.
           q2 = (q2kin(ikinfit) + q2kin(ikinfit+1))/2.
           write(22,833) ix,iy,ttit(itt),x,q2,zfit,
     >      npts,chi2/float(npts),ymin,ymax
 833       format('set window x ',i1,' of 2 y ',
     >                            i1,' of 2'/
     >      'set color white'/ 
     >      'title top ',1h',a2,' x=',f4.2,' Q2=',f3.1,
     >      ' z=',f4.2,' npt=',i3,' chi2/df=',f5.1,1h'/
     >      'title bottom ',1h','phi*',1h'/
     >      'title left ',1h','Multiplicity',1h'/
     >      'set limits x 0. 6.3 y ',2f10.3/
     >      'set scale y log ; set bar size 0.'/
     >      'set symbol 9O size 1.0 ; set order x y dy')
           write(23,834) ix,iy,ttit(itt),x,q2,zfit,
     >      npts,chi2/float(npts),ymax/3.,ymax
 834       format('set window x ',i1,' of 2 y ',
     >                            i1,' of 2'/
     >      'title top ',1h',a2,' x=',f4.2,' Q2=',f3.1,
     >      ' z=',f4.2,' npt=',i4,' chi2/df=',f5.1,1h'/
     >      'title bottom ',1h','phi*',1h'/
     >      'title left ',1h','Multiplicity',1h'/
     >      'set limits x 0. 6.3 y ',2f10.3/
     >      'set bar size 0.'/
     >      'set symbol 9O size 1.0 ; set order x y dy')
           do ipt=1,10
            if(ipt.eq.1) write(22,'(''set color white'')')
            if(ipt.eq.2) write(22,'(''set color magenta'')')
            if(ipt.eq.3) write(22,'(''set color red'')')
            if(ipt.eq.4) write(22,'(''set color green'')')
            if(ipt.eq.5) write(22,'(''set color blue'')')
            if(ipt.eq.6) write(22,'(''set color cyan'')')
            if(ipt.eq.7) write(22,'(''set color white'')')
            if(ipt.eq.8) write(22,'(''set color magenta'')')
            if(ipt.eq.9) write(22,'(''set color red'')')
            if(ipt.eq.10) write(22,'(''set color green'')')
            if(ipt.eq.11) write(22,'(''set color blue'')')
            nplt=0
            do iphi=1,14,2
             phimin = 2.*3.1415928/15.*(iphi-1)
             phimax = 2.*3.1415928/15.*(iphi+1)
             sum1=0.
             sum2=0.
             do jj=1,npts
c              if(spt(jj).gt.0.06*(ipt-1) .and.
c     >           spt(jj).le.0.06*(ipt+1).and.
c use all pt bins
              if(spt(jj).gt.0.06*(ipt-0.5) .and.
     >           spt(jj).le.0.06*(ipt+0.5).and.
     >           sphi(jj).gt.phimin .and.
     >           sphi(jj).lt.phimax) then
               sum1 = sum1 + sy(jj)/syer(jj)**2
               sum2 = sum2 + 1./syer(jj)**2
              endif
             enddo
             if(sum2.gt.0.) then
              sum1 = sum1 / sum2
              sum2 = 1./sqrt(sum2)
              write(22,'(3f7.3)') (phimin + phimax)/2.,sum1,sum2
              nplt = nplt + 1
             endif
            enddo
            if(nplt.gt.0) then
             lastframeused = .true.
             write(22,'(''plot'')')
             pt = 0.06*(ipt)
             pt2 = pt**2
             do jj=1,21
              phi = 2.*3.1415928*(jj-1)/20.
              sig = (coef(1)/zfit  ) / coef(2) * 
     >         exp(-pt2/coef(2)) *
     >         (1. + coef(3) * pt * cos(phi))
              write(22,'(3f8.3)') phi,sig
             enddo
             write(22,'(''join'')')
            endif
           enddo

c plot mult versus phi averaged over pt<0.25
           do iphi=1,15
            sig = 0.
            siger = 0.
            do jj=1,npts
             if(sphi(jj).gt. 2.*3.1415928/15.*(iphi-1) .and.
     >          sphi(jj).le. 2.*3.1415928/15*iphi .and. 
     >          spt(jj).le.0.25) then
              sig = sig + sy(jj)/syer(jj)**2
              siger = siger + 1./syer(jj)**2
             endif
            enddo
            sig = sig / siger
            siger = 1./sqrt(siger)
            if(siger.ne.0.) write(23,'(3f10.4)') 
     >       2.*3.1415928/15.*(iphi-0.5),sig,siger
           enddo
           write(23,'(''plot'')')
           pt = 0.12
           pt2 = pt**2
           do jj=1,21
            phi = 2.*3.1415928*(jj-1)/20.
            sig = (coef(1)/zfit  ) / coef(2) * 
     >       exp(-pt2/coef(2)) /2. / pi **
     >       (1. + coef(3) * pt * cos(phi))
            write(23,'(3f8.3)') phi,sig
           enddo
           write(23,'(''join'')')
          endif
         endif
        enddo
       enddo
      enddo
c xxx for test
      goto 77
! get csv from d only using modified fdds plus Geiger and both total integral
! and values at pt=0 for ratio
      do ikin=1,3
       q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
       x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
       do iz=5,18
        z= 0.05*(iz-0.5)
        wp2 = am**2 + q2*(1/x - 1)*(1-z)
        Mm = scsv(ikin,4,iz,1,ifit)
        Mp = scsv(ikin,2,iz,1,ifit)
        Mmer = scsver(ikin,4,iz,1,ifit)
        Mper = scsver(ikin,2,iz,1,ifit)
        Bm = scsv(ikin,4,iz,2,ifit)
        Bp = scsv(ikin,2,iz,2,ifit)
        if(Mm.ne.0. .and. Mp.ne.0.) then
         pt=0.
         phicm=0.
         it=4
         call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >    sighad,u,ub,d,db,u1,d1,
     >    s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
         Rseans(ikin,ikin) = 5. * (ub + db) / 
     >   (u - ub + d - db)
         Rseas(ikin,iz) = (s + sb) / (u - ub + d - db)
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        rnew = fd1 / fu1
c for test
          rnew = rnew * 1.35
c swith to geiger: works better
         call rgeiger(z,rnew)
         Dprop = (1. - rnew) / (1. + rnew)
         mprat = Mm / Mp 
         RD = (4.*mprat - 1.) / (1. - mprat)
         mprat = (Mm + Mmer) / Mp 
         RD1 = (4.*mprat - 1.) / (1. - mprat)
         mprat = Mm / (Mp + Mper) 
         RD2 = (4.*mprat - 1.) / (1. - mprat)
         dRd = sqrt((rd1-rd)**2 + (rd2-rd)**2)
         CSVprop(ikin,iz,ifit,1) = 2.5 + rseans(ikin,iz) + 
     >    rseas(ikin,iz) - Dprop * (2.5 + RD)
c this is delta u / (u + d) if delta u = -1 *delta d
c ignoring sea
         chk1 = (4*rnew + 1 - mprat * (4 + rnew)) /
     >         (4*rnew - 1 - mprat * (4 - rnew))
c this is -4/3 (delta u - delta d) / (u + d)
c or -8/3 delta u / (u+d)
         chk2 = 2.5  - Dprop * (2.5 + RD)
         CSVproper(ikin,iz,ifit,1) = Dprop * dRD
         mprat = (Mm / Bm) / (Mp / Bp) 
         RD = (4.*mprat - 1.) / (1. - mprat)
         CSVprop(ikin,iz,ifit,2) = 2.5 + rseans(ikin,iz) + 
     >    rseas(ikin,iz) - Dprop * (2.5 + RD)
         CSVproper(ikin,iz,ifit,2) = 
     >    CSVproper(ikin,iz,ifit,1) 
         write(6,'(''csv'',2i3,6f7.2)') ikin,iz,
     >    (csvprop(ikin,iz,ifit,k),
     >    CSVproper(ikin,iz,ifit,k),k=1,2),chk1,
     >    -3./8.*chk2
         do k=1,2
          if(wp2 .gt. 3.0) then
           avcsv(ikin,ifit,k) = avcsv(ikin,ifit,k) +
     >      csvprop(ikin,iz,ifit,k) /
     >      csvproper(ikin,iz,ifit,k)**2
           avcsver(ikin,ifit,k) = avcsver(ikin,ifit,k) +
     >      1./csvproper(ikin,iz,ifit,k)**2
          endif
         enddo
        endif
       enddo ! iz
       do k=1,2
        avcsv(ikin,ifit,k) = avcsv(ikin,ifit,k) /
     >    avcsver(ikin,ifit,k)
        avcsver(ikin,ifit,k) = 1. /
     >   sqrt(avcsver(ikin,ifit,k))
       enddo
       write(6,'(''avcsv'',2i3,4f8.3)') ikin,ifit,
     >  (avcsv(ikin,ifit,k),avcsver(ikin,ifit,k),k=1,2)
      enddo ! ikin

! get csv, fav, fav/unfav from all four it cases
c not fitting d/u because results are not sensitive to this:
c would have to constrain inclusive to remain unchanged
      do ikinx=1,8
       ikin = ikinx
       if(ikinx.eq.8) ikin=15
       q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
       x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
c get pdfs 
       do iz=5,18
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,mu,mub,md,mdb,u1,d1,
     >   ms,msb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        z = 0.05*(iz-0.5)
        wp2 = am**2 + q2*(1/x - 1)*(1-z)
        ngood=0
        BigS = 2.*(ms + msb) / (mu + md)
        BigT = (mub + mdb) / (mu + md)
        uplusd(ifit,ikinx) = mu + md
        do it=1,4
         Mfit(it) = scsv(ikin,it,iz,1,ifit) * 2. * 3.1415928
c try using M0 * b instead of M0. This makes CSV even bigger
c and also diff ratio worse too
c     >     * scsv(ikin,it,iz,2,ifit) / 0.2
         Mfiter(it) = scsver(ikin,it,iz,1,ifit) * 2. * 3.1415928
c     >     * scsv(ikin,it,iz,2,ifit) / 0.2
         if(Mfiter(it).ne.0 .and. 
     >    Mfiter(it).lt.0.1) ngood = ngood + 1
c values using average mult.
         mfitp(it)=0.
         mfitper(it)=0.
          itt=it
          if(it.gt.2) itt=it+1
          if(avmer(2*ikin-1,iz,itt).ne.0. .and.
     >      avmer(2*ikin  ,iz,itt).ne.0.) then
           do ikinp = 2*ikin-1, 2*ikin
            mfitp(it) = mfitp(it) + 
     >      (avm(ikinp,iz,itt) - 
     >       rhofactp * avmrho(ikinp,iz,itt)) /
     >       avmer(ikinp,iz,itt)**2
            write(6,'(''dbg'',4i3,f6.2,4f8.3)') ikinp,iz,itt,
     >       ifit,rhofactp,avm(ikinp,iz,itt), 
     >       avmer(ikinp,iz,itt),
     >       rhofactp * avmrho(ikinp,iz,itt),
     >       avmrho(ikinp,iz,itt)
            mfitper(it) = mfitper(it) + 1./
     >       avmer(ikinp,iz,itt)**2
           enddo
           mfitp(it) = mfitp(it) / mfitper(it)
           mfitper(it) = 1./sqrt(mfitper(it))
         endif
         write(6,'(''mfitp'',3i3,8f7.3)') ikin,iz,it,
     >    mfit(it),mfiter(it),mfitp(it),mfitper(it),
     >    avm  (2*ikin-1,iz,itt),
     >    avmer(2*ikin-1,iz,itt),
     >    avm  (2*ikin  ,iz,itt),
     >    avmer(2*ikin  ,iz,itt) 
 
c  for test, replace with averaged mult.
c         mfit(it) = mfitp(it)
c         mfiter(it) = mfitper(it)
        enddo !it

c for test
c this lowers funfav / ffav, but doesn't change delta u
c        mfit(3) = mfit(3) * 0.9
c        mfit(4) = mfit(4) * 0.9
c for test can change inclusive d/p ratio
        fact = 1.00
c        fact = 1.05
        mfit(2) = mfit(2) * fact
        mfiter(2) = mfiter(2) * fact
        mfit(4) = mfit(4) * fact
        mfiter(4) = mfiter(4) * fact
c  for test, increase d pi+ only
        fact = 1.00
        mfit(2) = mfit(2) * fact
        mfiter(2) = mfiter(2) * fact
c  for test, increase d pi+ only
        fact = 1.00
        mfitp(2) = mfitp(2) * fact
        mfitper(2) = mfitper(2) * fact
        if(ngood.eq.4) then
c get starting z*FF/z (matches mfit which is also scaled by z)
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
c         nparam=4
c         npar = 4
         nparam=3
         npar = 3
c for test, try varying csv by hand
         do j=1,10
          p(1) = fu1
          p(2) = fd1 / fu1
          p(3) = 0.05 * float(j-5)
          p(4) = -p(3)
          call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
          write(6,'(''mt'',10f7.3)') 
     >     p(3),chi2,(mfit(it),mfitr(it),it=1,4)
         enddo
! fit using mfit
         mfitflag = 0
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
c         p(4) = 0.0 ! deltd
c         p(3) = md / mu ! d/u
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit I'',2i2,i3,f8.2/8f8.3/3f8.3/6f7.3)') 
     >     itfit,ikin,iz,chi2,(mfit(it),mfitr(it),it=1,4),
     >     (p(it),it=1,3),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         arglis(1)=0
         call mnexcm(mfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
          mcsv(ifit,ikinx,iz,j)= coef(j)
          mcsver(ifit,ikinx,iz,j)= std(j)
         enddo
         fu1sv(ifit,ikinx,iz)= fu1
         fd1sv(ifit,ikinx,iz)= fd1
         fact=1.
c delu, deld from sum, diff
         call getcsv(delu,deluer,deld,delder,doverp,
     >     doverpp, doverpper,fact)
         mcsv(ifit,ikinx,iz, 9)= delu / (mu + md)
         mcsver(ifit,ikinx,iz, 9)= deluer / (mu + md)
         mcsv(ifit,ikinx,iz,10)= deld / (mu + md)
         mcsver(ifit,ikinx,iz,10)= delder / (mu + md)

! fit using mfitp now
         mfitflag = 1
         nparam=3
         npar = 3
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
c         p(4) = 0.0 ! deltd
c        p(4) = md / mu         ! d/u
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit Ip'',2i2,i3,f8.2/8f8.3/3f8.3/6f7.3)') 
     >     itfit,ikin,iz,chi2,(mfit(it),mfitr(it),it=1,4),
     >     (p(it),it=1,3),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         arglis(1)=0
         call mnexcm(mfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
          mcsvp(ifit,ikinx,iz,j)= coef(j)
          mcsvper(ifit,ikinx,iz,j)= std(j)
         enddo
c influeance of d/u ratio (usng mfitp)
         do j=1,5
          fact = 0.7 + 0.1*j
          call getcsv(delu,deluer,deld,delder,doverp,
     >    doverpp, doverpper,fact )
          write(6,'(''duchk'',3i3,7f7.3)') ifit,ikinx,iz,
     >     fact,delu,deluer,deld,delder,doverpp/doverp,doverp
         enddo
         fact=1.
         call getcsv(delu,deluer,deld,delder,doverp,
     >    doverpp, doverpper,fact)
         mcsvp(ifit,ikinx,iz, 9)= delu / (mu + md)
         mcsvper(ifit,ikinx,iz, 9)= deluer / (mu + md)
         mcsvp(ifit,ikinx,iz,10)= deld / (mu + md)
         mcsvper(ifit,ikinx,iz,10)= delder / (mu + md)
         dopsv(ifit,ikinx,iz)= doverp          
         mcsvp(ifit,ikinx,iz,11)= doverpp          
         mcsvper(ifit,ikinx,iz,11)= doverpper
         write(6,'(''cmpcsv'',8f7.3)') delu,deluer,deld,delder,
     >    coef(3)/delu,std(3)/deluer,coef(4)/deld,std(4)/delder
         write(6,'(''doverp'',3i3,3f8.3)') ifit,ikinx,iz,
     >    doverp, doverpp, doverpper
         avdelu(ifit,ikinx) = avdelu(ifit,ikinx)  +
     >    delu / (mu + md) / (deluer / (mu + md))**2
         avdeluer(ifit,ikinx) = avdelu(ifit,ikinx)  +
     >    1. / (deluer / (mu + md))**2
         avdeld(ifit,ikinx) = avdeld(ifit,ikinx)  +
     >    deld / (mu + md) / (delder / (mu + md))**2
         avdelder(ifit,ikinx) = avdeld(ifit,ikinx)  +
     >    1. / (delder / (mu + md))**2

c also get delu from deuteron only for delu + deld = 0
c this gives delu / (u + d), using mfit
         call rgeiger(z,rnew)
         mprat = mfit(4) / mfit(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mprat = (mfit(4) + mfiter(4)) / mfit(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >             (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mprat = mfit(4) / (mfit(2) + mfiter(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >             (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mcsv(ifit,ikinx,iz,4) = chk1
         mcsver(ifit,ikinx,iz,4) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c also get delu from deuteron only for delu + deld = 0
c this gives delu / (u + d), using mfitp
         a = 1/4
         call rgeiger(z,rnew)
         mprat = mfitp(4) / mfitp(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mcsvp(ifit,ikinx,iz,4) = chk1
         mcsvper(ifit,ikinx,iz,4) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c same thing but with DSS
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
         rnew = fd1 / fu1
         mprat = mfitp(4) / mfitp(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mcsvp(ifit,ikinx,iz,12) = chk1
         mcsvper(ifit,ikinx,iz,12) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c get fu / fd from proton data only using mfit
         mprat = mfit(3) / mfit(1)
         chk1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.* mub + ms + msb))
         chk2 = (mprat * (4.*mu + mdb) - md*1.1 - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md*1.1 + 4.* mub + ms + msb))
         write(6,'(''fup 1'',3f8.3)') mprat,chk1,chk2
         mprat = (mfit(3) + mfiter(3)) / mfit(1)
         chk1er1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfit(3) / (mfit(1) + mfiter(1))
         chk1er2 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsv(ifit,ikinx,iz,5) = chk1
         mcsver(ifit,ikinx,iz,5) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,5),mcsver(ifit,ikinx,iz,5)
c get fu / fd from proton data only using mfitp
         mprat = mfitp(3) / mfitp(1)
         chk1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.* mub + ms + msb))
         chk2 = (mprat * (4.*mu + mdb) - md*1.1 - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md*1.1 + 4.* mub + ms + msb))
         write(6,'(''fup 1'',3f8.3)') mprat,chk1,chk2
         mprat = (mfitp(3) + mfitper(3)) / mfitp(1)
         chk1er1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
c         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfitp(3) / (mfitp(1) + mfitper(1))
         chk1er2 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
c         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsvp(ifit,ikinx,iz,5) = chk1
         mcsvper(ifit,ikinx,iz,5) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fupp '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,5),mcsvper(ifit,ikinx,iz,5)
c get fu / fd from deuteron data only assuming csv=0 using mfit
         mprat = mfit(4 ) / mfit(2)
         chk1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 1'',2f8.3)') mprat,chk1
         mprat = (mfit(4) + mfiter(4)) / mfit(2)
         chk1er1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfit(4) / (mfit(2) + mfiter(2))
         chk1er2 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsv(ifit,ikinx,iz,6) = chk1
         mcsver(ifit,ikinx,iz,6) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,6),mcsver(ifit,ikinx,iz,6)
c get fu / fd from deuteron data only w/csv=0 using mfitp
         mprat = mfitp(4 ) / mfitp(2)
         chk1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
c         write(6,'(''fup 1'',2f8.3)') mprat,chk1
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
c         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsvp(ifit,ikinx,iz,6) = chk1
         mcsvper(ifit,ikinx,iz,6) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,6),mcsvper(ifit,ikinx,iz,6)
c 5-param fit using mfitp. Includes d/u with constraint on 
c inclusive p/d remaining unchanged
         nparam=5
         npar = 5
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
         p(4) = 0.0 ! deltd
         p(5) = 1.0 ! d/p
         call mfit5_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit5 Ip'',2i2,i3,f8.2/8f8.3/5f8.3/6f7.3)') 
     >     ifit,ikin,iz,chi2,(mfitp(it),mfitr(it),it=1,4),
     >     (p(it),it=1,5),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
c get better guess for p1
         p(1) = p(1) * 
     >   (mfitp(1)+mfitp(2)+mfitp(3)+mfitp(4)) / 
     >   (mfitr(1)+mfitr(2)+mfitr(3)+mfitr(4))
         call mfit5_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit5 IIp'',2i2,i3,f8.2/8f8.3/5f8.3/6f7.3)') 
     >     ifit,ikin,iz,chi2,(mfitp(it),mfitr(it),it=1,4),
     >     (p(it),it=1,5),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         call mnparm( 1,"P1 ",p(1), 0.001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.001D0,zero,zero,ierflg)
         call mnparm( 4,"P4 ",p(4), 0.001D0,zero,zero,ierflg)
         call mnparm( 5,"P5 ",p(5), 0.001D0,zero,zero,ierflg)
         arglis(1)=0
         call mnexcm(mfit5_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit5 '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
           mcsv5p(ifit,ikinx,iz,j)= coef(j)
          mcsv5per(ifit,ikinx,iz,j)= std(j)
         enddo
c get sum and diff ratios using mfit
         chk1 = (mfit(2) - mfit(4))/
     >          (mfit(1) - mfit(3))
         chk1er1 = (mfit(2)+mfiter(2) - mfit(4))/
     >             (mfit(1)           - mfit(3))
         chk1er2 = (mfit(2)           - mfit(4))/
     >             (mfit(1)+mfiter(1) - mfit(3))
         chk1er3 = (mfit(2) - mfit(4) - mfiter(4))/
     >             (mfit(1) - mfit(3))
         chk1er4 = (mfit(2) - mfit(4))/
     >             (mfit(1) - mfit(3) - mfiter(3))
         mcsv(ifit,ikinx,iz,7) = chk1
         mcsver(ifit,ikinx,iz,7) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fdif'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,7),mcsver(ifit,ikinx,iz,7)
         chk1 = (mfit(2) + mfit(4))/
     >          (mfit(1) + mfit(3))
         chk1er1 = (mfit(2)+mfiter(2) + mfit(4))/
     >             (mfit(1)           + mfit(3))
         chk1er2 = (mfit(2)           + mfit(4))/
     >             (mfit(1)+mfiter(1) + mfit(3))
         chk1er3 = (mfit(2) + mfit(4) + mfiter(4))/
     >             (mfit(1) + mfit(3))
         chk1er4 = (mfit(2) + mfit(4))/
     >             (mfit(1) + mfit(3) + mfiter(3))
         mcsv(ifit,ikinx,iz,8) = chk1
         mcsver(ifit,ikinx,iz,8) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fsum'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,8),mcsver(ifit,ikinx,iz,8)
c get sum and diff ratios using mfitp
         chk1 = (mfitp(2) - mfitp(4))/
     >          (mfitp(1) - mfitp(3))
         chk1er1 = (mfitp(2)+mfitper(2) - mfitp(4))/
     >             (mfitp(1)           - mfitp(3))
         chk1er2 = (mfitp(2)           - mfitp(4))/
     >             (mfitp(1)+mfitper(1) - mfitp(3))
         chk1er3 = (mfitp(2) - mfitp(4) - mfitper(4))/
     >             (mfitp(1) - mfitp(3))
         chk1er4 = (mfitp(2) - mfitp(4))/
     >             (mfitp(1) - mfitp(3) - mfitper(3))
         mcsvp(ifit,ikinx,iz,7) = chk1
         mcsvper(ifit,ikinx,iz,7) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fdif p'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,7),mcsvper(ifit,ikinx,iz,7)
         chk1 = (mfitp(2) + mfitp(4))/
     >          (mfitp(1) + mfitp(3))
         chk1er1 = (mfitp(2)+mfitper(2) + mfitp(4))/
     >             (mfitp(1)           + mfitp(3))
         chk1er2 = (mfitp(2)           + mfitp(4))/
     >             (mfitp(1)+mfitper(1) + mfitp(3))
         chk1er3 = (mfitp(2) + mfitp(4) + mfitper(4))/
     >             (mfitp(1) + mfitp(3))
         chk1er4 = (mfitp(2) + mfitp(4))/
     >             (mfitp(1) + mfitp(3) + mfitper(3))
         mcsvp(ifit,ikinx,iz,8) = chk1
         mcsvper(ifit,ikinx,iz,8) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fsum p'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,8),mcsvper(ifit,ikinx,iz,8)

c back to multi-param fit again using mfit
         call mfit_fcn(npar,grad,chi2,coef,ierflg,futil)
         write(6,'(''mfit F'',2i2,i3,f8.2/8f8.3)') 
     >     itfit,ikinx,iz,chi2,(mfit(it),mfitr(it),it=1,4)
         write(17,'(i2,i2,i3,11f7.3)') ifit,ikinx,iz,
     >     chi2,coef(1)/fu1,std(1)/fu1,
     >    coef(2) / (fd1 / fu1),
     >    std(2) / (fd1 / fu1),
     >    coef(3),std(3),
     >    mcsv(ifit,ikinx,iz,4),
     >    mcsver(ifit,ikinx,iz,4)
        endif
       enddo ! iz
       avdelu(ifit,ikinx) = avdelu(ifit,ikinx) / avdeluer(ifit,ikinx)
       avdeld(ifit,ikinx) = avdeld(ifit,ikinx) / avdelder(ifit,ikinx)
      enddo ! ikin
c if skipping all the csv stuff for test
 77   continue
c end of big loop over ifit
      enddo ! ifit

      goto 99

c big global fit to all of the data
      do ifit= 1,2
      rhofact = 1.0
      write(6,'(''using rho fact = '',i2,f5.1)') ifit,rhofact
      if(ifit.eq.2) rhofact = 0.
c xxx fitx this in future?
      if(ifit.eq.3) then
       do i=1,bnpt
        brho(i) = brhop(i)
       enddo
      endif
      nparam=19
      npar = 19
      p(1) = 0.2
      p(2) = 0.2
      p(3) = 0.2
      p(4) = 0.2
      p(5) = 0.
      p(6) = 0.
      p(7) = 0.
      p(8) = 0.
      p(9) = 0.0
      p(10) = 0.
      do i=11,19
       p(i) = 0.
      enddo
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo
      call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
      write(6,'(''initial chi2'',f8.2)') chi2/float(bnpt)
c find best norm. factor for each kin.
      do ikin=1,32
       chibest=10000000.
       fbest=0.
       do kk=1,32
        fnorm(kk)=1.0
       enddo
       do k=1,40
        fnorm(ikin) = 0.8 + 0.01 * k
        call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
        if(chi2.lt.chibest) then
         fbest=fnorm(ikin)
         chibest = chi2
        endif
        fnormsv(ikin)=fbest
       enddo
      enddo
      do ikin=1,32
       fnorm(ikin)=fnormsv(ikin)
      enddo
      do ikin=1,32,2
       write(6,'(''ikin fnorm='',i3,3f6.2,3f7.3)') ikin,
     >  (xkin(ikin) + xkin(ikin+1))/2.,
     >  (q2kin(ikin) + q2kin(ikin+1))/2.,
     >  (wkin(ikin) + wkin(ikin+1))/2.,
     >  fnorm(ikin),fnorm(ikin+1),
     >  fnorm(ikin)/fnorm(ikin+1)
      enddo

c set norm factors back to 1.0
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo

      call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
      write(6,'(''initial chi2'',f8.2)') chi2/float(bnpt)
      do ikin=0,36
       if(bdfk(ikin).gt.0) then
        if(ikin.gt.0) then
         write(6,'(''kin,df,chi2/df'',i3,3f6.2,f7.0,f7.2)') ikin,
     >   xkin(ikin),q2kin(ikin),wkin(ikin),
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
        else
         write(6,'(''kin,df,chi2/df'',i3,18x,f7.0,f7.2)') ikin,
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
        endif
       endif
      enddo
      write(6,'(''inital x,q2,w,z,pt,mx'')')
      do k=1,15
       write(6,'(6(i5,f4.1))') 
     >  int(bdfx(k)),bchi2x(k)/max(1.,bdfx(k)),
     >  int(bdfq2(k)),bchi2q2(k)/max(1.,bdfq2(k)),
     >  int(bdfw(k)),bchi2w(k)/max(1.,bdfw(k)),
     >  int(bdfz(k)),bchi2z(k)/max(1.,bdfz(k)),
     >  int(bdfpt(k)),bchi2pt(k)/max(1.,bdfpt(k)),
     >  int(bdfmx(k)),bchi2mx(k)/max(1.,bdfmx(k))
      enddo
      do k=1,6
       write(6,'(i2,i7,f6.2)') k,int(bdft(k)),
     >   bchi2t(k)/max(1.,bdft(k))
      enddo

! initial parameters
      call mnparm( 1,"P1 ",p(1), 0.001D0,.19D0,.21D0,ierflg)
      call mnparm( 2,"P2 ",p(2), 0.001D0,.19D0,.21D0,ierflg)
      call mnparm( 3,"P3 ",p(3), 0.001D0,zero,zero,ierflg)
      call mnparm( 4,"P4 ",p(4), 0.001D0,zero,zero,ierflg)
      call mnparm( 5,"P5 ",p(5), 0.001D0,zero,zero,ierflg)
      call mnparm( 6,"P6 ",p(6), 0.001D0,zero,zero,ierflg)
      call mnparm( 7,"P7 ",p(7), 0.001D0,zero,zero,ierflg)
      call mnparm( 8,"P8 ",p(8), 0.001D0,zero,zero,ierflg)
      call mnparm( 9,"P9 ",p(9), 0.001D0,zero,zero,ierflg)
      call mnparm( 10,"P10 ",p(10), 0.001D0,zero,zero,ierflg)
      call mnparm( 11,"P11 ",p(11), 0.001D0,zero,zero,ierflg)
      call mnparm( 12,"P12 ",p(12), 0.001D0,zero,zero,ierflg)
      call mnparm( 13,"P13 ",p(13), 0.001D0,zero,zero,ierflg)
      call mnparm( 14,"P14 ",p(14), 0.001D0,zero,zero,ierflg)
      call mnparm( 15,"P15 ",p(15), 0.001D0,zero,zero,ierflg)
      call mnparm( 16,"P16 ",p(16), 0.001D0,zero,zero,ierflg)
      call mnparm( 17,"P17 ",p(17), 0.001D0,zero,zero,ierflg)
      call mnparm( 18,"P18 ",p(18), 0.001D0,zero,zero,ierflg)
      call mnparm( 19,"P19 ",p(19), 0.001D0,zero,zero,ierflg)
      ARGLIS(1) = 5.
      CALL MNEXCM(bfit_fcn,'FIX', ARGLIS ,1,IERFLG)
      arglis(1)=0
      call mnexcm(bfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
      
      do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         write(6,'(''bf '',i2,8f7.3)') j,coef(j),std(j)
      enddo
      call bfit_fcn(npar,grad,chi2,coef,ierflg,futil)
      write(6,'(''final chi2/d.f.'',f8.2)') chi2/float(bnpt)
      do ikin=1,36
       if(bdfk(ikin).gt.0) then
        write(6,'(''kin,df,chi2/df'',i3,3f6.2,f7.0,f7.2)') ikin,
     >   xkin(ikin),q2kin(ikin),wkin(ikin),
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
       endif
      enddo
      write(6,'(''final x,q2,w,z,pt,mx'')')
      do k=1,15
       write(6,'(6(i5,f4.1))') 
     >  int(bdfx(k)),bchi2x(k)/max(1.,bdfx(k)),
     >  int(bdfq2(k)),bchi2q2(k)/max(1.,bdfq2(k)),
     >  int(bdfw(k)),bchi2w(k)/max(1.,bdfw(k)),
     >  int(bdfz(k)),bchi2z(k)/max(1.,bdfz(k)),
     >  int(bdfpt(k)),bchi2pt(k)/max(1.,bdfpt(k)),
     >  int(bdfmx(k)),bchi2mx(k)/max(1.,bdfmx(k))
      enddo
      do k=1,6
       write(6,'(i2,i7,f6.2)') k,int(bdft(k)),
     >   bchi2t(k)/max(1.,bdft(k))
      enddo

      do i=1,bnpt
       if(bkin(i).eq.0) then
        write(6,'(''incl'',5f7.3)') bxv(i),bq2v(i),byv(i),
     >   bsigv(i),byv(i)/bsigv(i)
       endif
      enddo

      enddo ! j fits
 99   continue

      close(unit=22)
      open(unit=22,file='ptmrho.top')
      write(22,'(''set device postscript'')')
      do ipm=1,2
       write(22,122) ipm
 122   format('set window x ',i1,' of 2'/
     >  'title bottom',1h','hepgen/SIMC',1h')
       do kk=1,25
        write(22,'(f7.2,i8)') 0.25*(kk-0.5),rhorat(kk,ipm)
       enddo
       write(22,'(''hist'')')
      enddo

c plot 3-parameter fit values
      close(unit=31)
      open(unit=31,file='ptmmult.top')
      write(31,'(''set device postscript'')')
      close(unit=32)
      open(unit=32,file='ptmpt2.top')
      write(32,'(''set device postscript'')')
      close(unit=33)
      open(unit=33,file='ptmcos.top')
      write(33,'(''set device postscript'')')
      close(unit=34)
      open(unit=34,file='ptmcos2.top')
      write(34,'(''set device postscript'')')
      ikin=0
      do iy=1,4
       do ix=1,4
        ikin = ikin + 1
        xmin=1.0 + 2.5*(ix-1)
        xmax=xmin+2.5
        ymax=9.0-2.0*(iy-1)
        ymin=ymax-2.0
        do j=1,4
         if(j.eq.1) yymin = 0.
         if(j.eq.1) yymax = 0.099
         if(j.eq.2) yymin = 0.
         if(j.eq.2) yymax = 0.499
         if(j.eq.3) yymin = -1.
         if(j.eq.3) yymax = 0.999
         if(j.eq.4) yymin = -1.
         if(j.eq.4) yymax = 0.999
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         if(iy.eq.4) write(30+j,'(1x,''set labels bottom on'')')
         if(iy.lt.4) write(30+j,'(1x,''set labels bottom off'')')
         write(30+j,131) xmin,xmax,ymin,ymax,yymin,
     >    yymax,yymin + 0.9*(yymax-yymin),
     >    (xkin(2*ikin-1) + xkin(2*ikin))/2.,
     >    (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
 131     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.79 y',2f8.3/
     >     1x,'set ticks size 0.03 ; set labels size 1.0'/
     >     1x,' set sym 9O size 0.2'/
     >     1x,'title 5. 0.6 size 1.5',1h','z',1h'/
     >     1x,'title 0.3 ',f8.3,' data size 1.0 ',1h',
     >     '(',f4.2,',',f3.1,')',1h')
         if(j.eq.1) write(31,132)
 132     format('title 0.1 4.6 angle=90 size=1.5',1h','zM',1h')
         if(j.eq.1) write(32,133)
 133     format('title 0.1 4.6 angle=90 size=1.5',1h','<P0t1>223',1h'/
     >          'case                          ',1h','  X X X X',1h')
         if(j.eq.3) write(33,134)
 134     format('title 0.5 4. angle=90 size=1.5',1h','A',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.4) write(34,135)
 135     format('title 0.5 4. angle=90 size=1.5',1h','B',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         write(30+j,'(''plot axes'')')
         do itt=1,4
          if(itt.eq.1) write(30+j,'(''set color red'')')
          if(itt.eq.2) write(30+j,'(''set color green'')')
          if(itt.eq.3) write(30+j,'(''set color blue'')')
          if(itt.eq.4) write(30+j,'(''set color white'')')
          nplt=0
          ifit=1
          do iz=1,20
           if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >        scsver(ikin,itt,iz,j,ifit).lt.0.2) then
            nplt = nplt + 1
            write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >       scsv(ikin,itt,iz,j,ifit),
     >       scsver(ikin,itt,iz,j,ifit)
 136        format(5f10.4)
           endif
          enddo
          if(nplt.gt.0) write(30+j,'(''plot'')')
         enddo
        enddo
       enddo
      enddo
c plot results from ptsidis runs only
c plot 3-parameter fit values
c "n" is for narrow delta cuts
c results look quite similar to wide cuts
      close(unit=31)
      open(unit=31,file='ptmmults.top')
c      open(unit=31,file='ptmmultsn.top')
      write(31,'(''set device postscript'')')
      close(unit=32)
      open(unit=32,file='ptmpt2s.top')
c      open(unit=32,file='ptmpt2sn.top')
      write(32,'(''set device postscript'')')
      close(unit=33)
      open(unit=33,file='ptmcoss.top')
c      open(unit=33,file='ptmcossn.top')
      write(33,'(''set device postscript'')')
      close(unit=34)
      open(unit=34,file='ptmcos2s.top')
c      open(unit=34,file='ptmcos2sn.top')
      write(34,'(''set device postscript'')')
      open(unit=35,file='ptmffs.top')
c      open(unit=35,file='ptmffn.top')
      write(35,'(''set device postscript'')')
      open(unit=36,file='ptmfus.top')
c      open(unit=36,file='ptmfun.top')
      write(36,'(''set device postscript'')')
      open(unit=37,file='ptmcsvs.top')
c      open(unit=37,file='ptmcsvn.top')
      write(37,'(''set device postscript'')')
      open(unit=38,file='ptmcsvsd.top')
c      open(unit=38,file='ptmcsvnd.top')
      write(38,'(''set device postscript'')')
      open(unit=39,file='ptmfups.top')
c      open(unit=39,file='ptmfupsn.top')
      write(39,'(''set device postscript'')')
      open(unit=40,file='ptmfuds.top')
c      open(unit=40,file='ptmfudsn.top')
      write(40,'(''set device postscript'')')
      open(unit=41,file='ptmdiff.top')
c      open(unit=41,file='ptmdiffn.top')
      write(41,'(''set device postscript'')')
      open(unit=42,file='ptmsums.top')
c      open(unit=42,file='ptmsumsn.top')
      write(42,'(''set device postscript'')')
      ikin=0
      do iy=1,1
       do ix=1,3
        ikin = ikin + 1
        xmin=1.0 + 3.5*(ix-1)
        xmax=xmin + 3.5
        ymax=9.0-4.5*(iy-1)
        ymin=ymax-4.5
        do j=1,12
         if(j.eq.1) yymin = 0.
         if(j.eq.1) yymax = 0.15
         if(j.eq.2) yymin = 0.
         if(j.eq.2) yymax = 0.5
         if(j.eq.3) yymin = -1.
         if(j.eq.3) yymax = 1.
         if(j.eq.4) yymin = -1.
         if(j.eq.4) yymax = 1.0
         if(j.eq.5) yymin = 0.5
         if(j.eq.5) yymax = 1.5
         if(j.eq.6) yymin = 0. 
         if(j.eq.6) yymax = 1.0
         if(j.eq.7) yymin = -1.0
         if(j.eq.7) yymax = 1.0
         if(j.eq.8) yymin = -1.0
         if(j.eq.8) yymax = 1.0
         if(j.eq.9) yymin = 0. 
         if(j.eq.9) yymax = 1.0
         if(j.eq.10) yymin = 0. 
         if(j.eq.10) yymax = 1.0
         if(j.eq.11) yymin = 0. 
         if(j.eq.11) yymax = 2.0
         if(j.eq.12) yymin = 0. 
         if(j.eq.12) yymax = 2.0
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
         q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
         w2 = am**2 + q2*(1/x-1)
         zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
         zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
         WRITE(30+J,231) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,x,q2
 231     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.89 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 0.6'/
     >     1x,'title 6.2 4.1 size 1.5',1h','z',1h'/
     >     1x,'title top size 1.5 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X',1h'/
     >     'title ',1h','VERY PRELIMINARY',1H')
         if(j.eq.1) write(31,232)
 232     format('title 0.35 6.8 angle=90 size=1.8',
     >     1h','z223 M001',1h'/'case ',1h',' X X  X X',1h')
         if(j.eq.2) write(32,233)
 233     format('title 0.5 6.8 angle=90 size=1.5',1h','<P0t1>223',1h'/
     >          'case                          ',1h','  X X X X',1h')
         if(j.eq.3) write(33,234)
 234     format('title 0.5 6.8 angle=90 size=1.5',1h','A',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.4) write(34,235)
 235     format('title 0.5 6.8 angle=90 size=1.5',1h','B',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.5) write(35,239)
 239     format('title 0.5 6.8 angle=90 size=1.5',1h','F/F0DSS1',1h'/
     >    ' case                                ',1h','   X   X',1h'/
     >    '0. 1. ; 1. 1. ; join dash')
         if(j.eq.6 .or. j.eq.9 .or. j.eq.10) write(36,237)
 237     format('title 0.5 6.8 angle=90 size=1.5',1h','F0u1/F0f1',1h'/
     >     'case                                ',1h',' X X  X X',1h')
         if(j.eq.7) write(37,238)
         if(j.eq.8) write(38,238)
 238     format('title 0.5 6.8 angle=90 size=1.5',1h','Du/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.11) write(41,288)
 288     format('title 0.5 6.8 angle=90 size=1.5',
     >            1h','(d2+3-d2-3)/(p2+3-p2-3)',1h'/
     >    'case ',1h','  X X  X X    X X  X X ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.12) write(42,289)
 289     format('title 0.5 6.8 angle=90 size=1.5',
     >            1h','(d2+3+d2-3)/(p2+3+p2-3)',1h'/
     >    'case ',1h','  X X  X X    X X  X X ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         write(30+j,'(''plot axes'')')
         wpcut=2.5
         write(30+j,236) zforwp225,yymin+(yymax-yymin)*0.06,
     >    zforwp225,yymin,
     >    zforwp225-0.01,yymin+(yymax-yymin)*0.09,wpcut
         wpcut=3.0
         write(30+j,236) zforwp23,yymin+(yymax-yymin)*0.06,
     >    zforwp23,yymin,
     >    zforwp23-0.01,yymin+(yymax-yymin)*0.09,wpcut
 236     format('arrow from ',2f8.3,' data to ',
     >    2f8.3,' data size 1.5'/
     >    'title ',2f8.3,' data size 0.8 ',1h',f3.1,1h')
         if(j.le.4) then
          do itt=1,4
           if(itt.eq.1) write(30+j,'(''set color red'')')
           if(itt.eq.2) write(30+j,'(''set color green'')')
           if(itt.eq.3) write(30+j,'(''set color blue'')')
           if(itt.eq.4) write(30+j,'(''set color white'')')
           ifit=1
           nplt=0
           do iz=1,20
            if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >       (scsver(ikin,itt,iz,j,ifit).lt.0.007 .or. j.gt.1).and.
     >         scsver(ikin,itt,iz,j,ifit).lt.0.2) then
             nplt = nplt + 1
             fact=1.
             if(j.eq.1) fact = 0.05*(iz-0.5) * 2. * 3.1415928
             write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >        fact*scsv(ikin,itt,iz,j,ifit),
     >        fact*scsver(ikin,itt,iz,j,ifit)
            endif
           enddo
           if(nplt.gt.0) then
             write(30+j,'(''plot'')')
             write(30+j,'(''set symbol size 0.5 ; plot'')')
             write(30+j,'(''set symbol size 0.4 ; plot'')')
             write(30+j,'(''set symbol size 0.3 ; plot'')')
            write(30+j,'(''set symbol size 0.2 ; plot'')')
            write(30+j,'(''set symbol size 0.1 ; plot'')')
           endif
c          write(30+j,'(''set symbol 1O size 0.5'')')
           ifit=2
           nplt=0
           do iz=1,20
            if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >       (scsver(ikin,itt,iz,j,ifit).lt.0.007 .or. j.gt.1).and.
     >         scsver(ikin,itt,iz,j,ifit).lt.0.2) then
            nplt = nplt + 1
            write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >       scsv(ikin,itt,iz,j,ifit),
     >       scsver(ikin,itt,iz,j,ifit)/100000.
            endif
           enddo
           if(nplt.gt.0) then
             write(30+j,'(''plot'')')
           endif
           if(j.eq.1) then
            it=itt
            if(itt.gt.2) it = itt+1
            do iz=4,18
             z = 0.05 * iz
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c correct b ratio of actual pt slope to one assumed
c in the fit
             dsigdz = dsigdz * sqrt(0.16**2 + z**2 * 0.25**2) /
     >                        sqrt(0.20**2 + z**2 * 0.12**2)
c            write(30+j,'(2f8.4)') z, z**2 * dsigdz/2./3.1415928
            enddo
c           write(30+j,'(''join'')')
c using DSS, still for ptmmult plot
            do iz=4,18
             z = 0.05 * iz
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,2,ipdf)
c             write(30+j,'(2f8.4)') z, z**2 * dsigdz/2./3.1415928
             write(30+j,'(2f8.4)') z, z**2 * dsigdz
            enddo
c           write(31,'(''set pattern .01 .03 .01 .03 ; join pattern'')')
            write(30+j,'(''join'')')
           endif

           if(j.eq.2 .and. itt.eq.4) then
            do iz=1,20
             z = 0.05 * iz
c            write(30+j,'(2f8.4)') z, 0.2 * (1. + z**2)
             write(30+j,'(2f8.4)') z,sqrt(0.16**2 + z**2 * 0.40**2)
            enddo
            write(30+j,'(''join'')')
           endif
           if(itt.eq.2) write(30+j,242) yymin +
     >      (yymax - yymin) * 0.17
 242       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' d P2+3',1h'/
     >      'case ',1h','   GX X',1h')
          if(itt.eq.3) write(30+j,243) yymin +
     >      (yymax - yymin) * 0.12
 243       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' p P2-3',1h'/
     >      'case ',1h','   GX X',1h')
           if(itt.eq.4) write(30+j,244) yymin +
     >     (yymax - yymin) * 0.07
 244       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' d P2-3',1h'/
     >      'case ',1h','   GX X',1h')
           if(itt.eq.1) write(30+j,241) yymin +
     >      (yymax - yymin) * 0.22
 241       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' p P2+3',1h'/
     >      'case ',1h','   GX X',1h')
          enddo ! itt
         endif ! j.le.4

         if(j.gt.4) then
          jp = j - 4
          do ifit=1,3
           nplt=0
           do iz=1,20
            if(mcsver(ifit,ikin,iz,jp).ne.0) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsv(ifit,ikin,iz,jp)*fact,
     >        mcsver(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
          enddo ! ifit
c plot both results from p only and d only on j=10 plot
          if(j.eq.10) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsver(ifit,ikin,iz,5).ne.0) then
              nplt = nplt + 1
              fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsv(ifit,ikin,iz,5)*fact,
     >         mcsver(ifit,ikin,iz,5)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.10
          if(jp.eq.1) then
           do iz=4,19
            z = 0.05*iz
            fact = 1.3 - 0.5 * (1-z)**4 / 0.7**4 -
     >                   0.2 * z**4 
            write(30+j,'(2f10.3)') z,fact
           enddo
           write(30+j,'(''join'')')
          endif
          if(jp.eq.2 .or. jp.eq.5 .or. jp.eq.6) then
c ratios from dss
           ifit=1
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             fact = fd1sv(ifit,ikin,iz)/fu1sv(ifit,ikin,iz)
             write(30+j,'(2f10.3)') z,fact
            endif
           enddo
           write(30+j,'(''join'')')
c ratios from geiger
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             call rgeiger(z,rnew)
             write(30+j,'(2f10.3)') z,rnew
            endif
           enddo
           write(30+j,'(''join dash'')')
          endif
          if(jp.eq.3.or.jp.eq.4) then
           x4 = x
           qgev = sqrt(q2)
           u = ctq5pdf(1,x4,qgev)
           d = ctq5pdf(2,x4,qgev)
! MRST fit. Note they have kappa = -0.2, but definition
! of delu is oppositie to mine
           f = 0.8 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
           f = -0.65 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
 251       format('0. ',f10.4,' ; 1. ',f10.4,' ; join dash')
          endif
          if(jp.eq.7.or.jp.eq.8) then
            do iz=3,18
             z = 0.05 * iz
             pt=0.
             phicm=0.
             call simcmodel(x,q2,z,pt,phicm,mmpi2,1 ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,2 ,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,4 ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,5 ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             if(jp.eq.7) write(30+j,'(2f10.4)') z,
     >        (sighad2 - sighad4) / (sighad1 - sighad3)
             if(jp.eq.8) write(30+j,'(2f10.4)') z,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
             v1 = (d - db) / (u - ub)
             write(6,'(''dbgd'',i3,2f6.2,8f7.2)') ikin,x,q2,
     >        sighad1,sighad2,
     >        sighad3,sighad4,
     >        (sighad2 - sighad4) / (sighad1 - sighad3),
     >        0.6 * (1 - v1/4.)/(1+v1/4.),v1,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
            enddo
            write(30+j,'(''join'')')
          endif
         endif ! j.gt.4
        enddo ! j
       enddo ! ix
      enddo ! iy

c for ptSIDIS runs ony, get four FF and plot
      CALL MNINIT(861,862,863)
      do ikin=1,3
       do iz=6,16
       do ifit=1,3
        z = 0.05 * (iz-0.5)
        q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
        x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
        it = 1
        pt = 0.
        phicm = 0.
        ipdf = 1
        if(ifit.eq.3) ipdf=2
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,muf,mubf,mdf,mdbf,u1,d1,
     >   msf,msbf,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        do it=1,4
         mv(it) = scsv(ikin,it,iz,1,ifit) * 2. * pi
         mver(it) = scsver(ikin,it,iz,1,ifit) * 2. * pi
        enddo
        write(6,'(''fffit'',3i3,8f7.3)') ikin,iz,ifit,
     >   (mv(it),mver(it),it=1,4)
c with four FF, errors too big
c        nparam=4
c        npar = 4
c chang to fit with two fav. one unfavored
        nparam=3
        npar = 3
        p(1) = fu1 ! fav u->pi+
        p(2) = fub ! unfav u->pi+
c        p(3) = fd1
c        p(4) = fdb
        p(3) = fdb
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         ffsv  (ifit,ikin,iz,j)= coef(j)
         ffsver(ifit,ikin,iz,j)= std(j)
         write(6,'(''fffit'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        enddo ! ifit
       enddo ! iz
      enddo ! ikin

c plot 3-paam FF results from pt-SIDIS runs
      close(unit=31)
      open(unit=31,file='ptmff3.top')
      write(31,'(''set device postscript'')')
      ikin=0
      do iy=1,1
       do ix=1,3
        ikin = ikin + 1
        xmin=1.0 + 3.5*(ix-1)
        xmax=xmin + 3.5
        ymax=9.0-4.5*(iy-1)
        ymin=ymax-4.5
        yymin=0.
        yymax=0.6
        x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
        q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
        w2 = am**2 + q2*(1/x-1)
        zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
        zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
        if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
        if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
        WRITE(31,2231) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,x,q2
 2231     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.79 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 1.0'/
     >     1x,'title 6.2 4.1 size 1.5',1h','z',1h'/
     >     1x,'title top size 1.5 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X',1h'/
     >     'title ',1h','VERY PRELIMINARY',1H')
        write(31,2232)
 2232   format('title 0.35 6.8 angle=90 size=1.8',
     >     1h','D(z)',1h')
        do j=1,3
         if(j.eq.2) write(31,'(''set color blue'')')
         if(j.eq.3) write(31,'(''set color red'')')
         ifit=1
         do iz=1,20
          if(ffsver(ifit,ikin,iz,j).ne.0. .and.
     >       ffsver(ifit,ikin,iz,j).lt.0.1) then
           nplt = nplt + 1
           fact=1.
c           fact = 0.05*(iz-0.5) * 2. * 3.1415928
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*ffsv  (ifit,ikin,iz,j),
     >        fact*ffsver(ifit,ikin,iz,j),
     >        fact*ffsv  (3,ikin,iz,j),
     >        fact*ffsver(3,ikin,iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''plot'')')
          write(31,'(''set symbol size 0.8 ; plot'')')
          write(31,'(''set symbol size 0.6 ; plot'')')
          write(31,'(''set symbol size 0.4 ; plot'')')
          write(31,'(''set symbol size 0.2 ; plot'')')
         endif
         ifit=2
         do iz=1,20
          if(ffsver(ifit,ikin,iz,j).ne.0. .and.
     >       ffsver(ifit,ikin,iz,j).lt.0.1) then
           nplt = nplt + 1
           fact=1.
c           fact = 0.05*(iz-0.5) * 2. * 3.1415928
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*ffsv  (ifit,ikin,iz,j),
     >        0.0 * fact*ffsver(ifit,ikin,iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''set symbol size 1.0 ; plot'')')
         endif
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fdb
         enddo
         write(31,'(''join'')')
c plot fdss again using zprime and NLO
         do iz=6,19
          z = 0.05 * iz
          xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
          zp = (z / 2.) * (xp / x) *(1. + sqrt(1 - 4 * x**2 * am**2 *  
     >       (0.02 + 0.1) / z**2 / q2**2))
          call fDSS (1,1,1, Zp, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fdb
         enddo
         write(31,'(''join dotdash'')')
        enddo ! j
       enddo ! ix
      enddo ! iy

c plot results from runs with both d and p (8 of them)
c using both intercept and averaged values
      close(unit=35)
      close(unit=36)
      close(unit=37)
      close(unit=38)
      close(unit=39)
      close(unit=40)
      close(unit=41)
      close(unit=42)
      close(unit=43)
      close(unit=44)
      close(unit=45)
      open(unit=35,file='ptmffp.top')
      write(35,'(''set device postscript'')')
      open(unit=36,file='ptmfup.top')
      write(36,'(''set device postscript'')')
      open(unit=37,file='ptmcsvp.top')
      write(37,'(''set device postscript'')')
      open(unit=38,file='ptmcsvpd.top')
      write(38,'(''set device postscript'')')
      open(unit=39,file='ptmfupp.top')
      write(39,'(''set device postscript'')')
      open(unit=40,file='ptmfudp.top')
      write(40,'(''set device postscript'')')
      open(unit=41,file='ptmdiffp.top')
      write(41,'(''set device postscript'')')
      open(unit=42,file='ptmsump.top')
      write(42,'(''set device postscript'')')
      open(unit=43,file='ptmdelup.top')
      write(43,'(''set device postscript'')')
      open(unit=44,file='ptmdeldp.top')
      write(44,'(''set device postscript'')')
      open(unit=45,file='ptmdoverp.top')
      write(45,'(''set device postscript'')')
      ikin=0
      do iy=1,2
       do ix=1,4
        ikin = ikin + 1
        ikinx = ikin
        if(ikin.eq.8) ikinx=15
        xmin=0.95 + 3.0*(ix-1)
        xmax=xmin + 3.0
        ymax=9.0-3.3*(iy-1)
        ymin=ymax-3.3
        do j=5,15
         if(j.eq.6) yymin = 0. 
         if(j.eq.6) yymax = 1.0
         if(j.eq.7) yymin = -0.9
         if(j.eq.7) yymax = 0.9
         if(j.eq.8) yymin = -0.9
         if(j.eq.8) yymax = 0.9
         if(j.eq.9) yymin = 0. 
         if(j.eq.9) yymax = 1.0
         if(j.eq.10) yymin = 0. 
         if(j.eq.10) yymax = 1.0
         if(j.eq.11) yymin = 0.0
         if(j.eq.11) yymax = 1.0
         if(j.eq.12) yymin = 0.5 
         if(j.eq.12) yymax = 1.5
         if(j.eq.13) yymin = -0.9
         if(j.eq.13) yymax =  0.9
         if(j.eq.14) yymin = -0.9
         if(j.eq.14) yymax =  0.9
         if(j.eq.15) yymin =  1.0
         if(j.eq.15) yymax =  1.99
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         if(iy.eq.1) write(30+j,'(1x,''set labels bottom off'')')
         if(iy.eq.2) write(30+j,'(1x,''set labels bottom on'')')
         x = (xkin(2*ikinx-1) + xkin(2*ikinx))/2.
         q2 = (q2kin(2*ikinx-1) + q2kin(2*ikinx))/2.
         w2 = am**2 + q2*(1/x-1)
         zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
         zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
         WRITE(30+J,431) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,yymin + 0.87*(yymax-yymin),x,q2,sqrt(w2)
 431     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.89 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 0.6'/
     >     1x,'title 0.3 ',f8.3,' data size 1.1 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223 W=',f3.1,' GeV',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X ',1h')
         if(j.eq.5.or.j.ge.8) write(30+j,433)
 433     format('0. 1. ; 1. 1. ; join dash')

         if(ix.eq.1 .and. iy.eq.1) then
          write(30+j,432)
 432      format('title 6.6 2.0 size 2.0',1h','z',1h')
          if(j.eq.5) write(35,439)
 439      format('title 0.2 4.8 angle=90 size=2.0',1h','F/F0DSS1',1h'/
     >     ' case                                ',1h','   X   X',1h')
         endif
         if(j.eq.6 .or. j.eq.9 .or. j.eq.10) write(36,437)
 437     format('title 0.2 4.8 angle=90 size=1.5',1h','F0u1/F0f1',1h'/
     >     'case                                ',1h',' X X  X X',1h')
         if(j.eq.7) write(37,438)
         if(j.eq.8) write(38,438)
 438     format('title 0.2 4.8 angle=90 size=2.0',1h','Du/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h')
         if(j.eq.11) write(41,488)
 488     format('title 0.2 4.8 angle=90 size=2.0',
     >            1h','(d2+3 - d2-3) / (p2+3 - p2-3)',1h'/
     >    'case ',1h','  X X    X X      X X    X X ',1h')
         if(j.eq.12) write(42,489)
 489     format('title 0.2 4.8 angle=90 size=2.0',
     >            1h','(d2+3 + d2-3) / (p2+3 + p2-3)',1h'/
     >    'case ',1h','  X X    X X      X X    X X ',1h')
         if(j.eq.13) write(30+j,438)
         if(j.eq.14) write(30+j,458)
 458     format('title 0.2 4.8 angle=90 size=2.0',1h','Dd/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h')
         if(j.eq.15) write(30+j,459)
 459     format('title 0.2 4.8 angle=90 size=2.0',1h',
     >    'inclusive  d/p',1h')
         write(30+j,'(''plot axes'')')
         wpcut=2.5
         write(30+j,236) zforwp225,yymin+(yymax-yymin)*0.06,
     >    zforwp225,yymin,
     >    zforwp225-0.01,yymin+(yymax-yymin)*0.09,wpcut
         wpcut=3.0
         write(30+j,236) zforwp23,yymin+(yymax-yymin)*0.06,
     >    zforwp23,yymin,
     >    zforwp23-0.01,yymin+(yymax-yymin)*0.09,wpcut

         if(j.gt.4) then
          jp = j - 4
          do ifit=1,3
           nplt=0
           do iz=1,20
            if(mcsvper(ifit,ikin,iz,jp).ne.0 .and.
     >         mcsvper(ifit,ikin,iz,jp).lt.0.3) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             if(jp.eq.3) fact = 1. / uplusd(ifit,ikin)
             write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsvp(ifit,ikin,iz,jp)*fact,
     >        mcsvper(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
          enddo ! ifit
c add pt, cosphi fit results for pt-sidsis runs
          if(ikin.le.3 .and. j.ne.15) then
           write(30+j,'(''set color blue'')')
           ifit=1
           do iz=1,20
            if(mcsver(ifit,ikin,iz,jp).ne.0) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             if(jp.eq.3) fact = 1. / uplusd(ifit,ikin)
             write(30+j,136) 0.01+0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsv(ifit,ikin,iz,jp)*fact,
     >        mcsver(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
           write(30+j,'(''set color white'')')
          endif
c plot both results from p only and d only on j=10 plot
          if(j.eq.10) then
           write(30+j,'(''set color blue'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsvper(ifit,ikin,iz,5).ne.0) then
              nplt = nplt + 1
              fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsvp(ifit,ikin,iz,5)*fact,
     >         mcsvper(ifit,ikin,iz,5)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.10
c plot both results from geiger and dss if j=8
          if(j.eq.8) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsvper(ifit,ikin,iz,12).ne.0) then
              nplt = nplt + 1
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsvp(ifit,ikin,iz,12),
     >         mcsvper(ifit,ikin,iz,12)
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.8
c plot delu, deld from 5-param fits
          if(j.eq.13 .or. j.eq.14) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsv5per(ifit,ikin,iz,j-10).ne.0) then
              nplt = nplt + 1
              fact = 1. / uplusd(ifit,ikin)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsv5p(ifit,ikin,iz,j-10)*fact,
     >         mcsv5per(ifit,ikin,iz,j-10)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.13 or 14
c plot d/p with delu = deld = 0. if ctegq (1), jam (3)
c also ratio from scalers
          if(j.eq.15) then
           write(30+j,'(''set color blue'')')
           do ifit=1,3,2
            nplt=0
            do iz=1,20
             if(dopsv(ifit,ikin,iz).ne.0) then
              nplt = nplt + 1
              write(30+j,136) 0.05*(iz-0.5),dopsv(ifit,ikin,iz)
             endif
            enddo
            if(nplt.gt.0.and.ifit.eq.1) write(30+j,'(''join'')')
            if(nplt.gt.0.and.ifit.eq.3) write(30+j,'(''join dash'')')
           enddo ! ifit
           write(30+j,'(''set color white'')')
           write(30+j,623) 2.*rrsv(ikin),2.*rrsv(ikin)
 623       format('0. ',f8.3, ' ; 1. ',f8.3, ' ; join')
           write(30+j,'(''set color green'')')
c checked in bcm.f that this is almost identical to calling
c f2_NMC_NEW with p and d targets 
          call FNP_NMC(X,Q2,ratt)
           write(30+j,623) 1.+ratt, 1.+ratt
           write(30+j,'(''set color white'')')
          endif ! j.eq.15
          if(jp.eq.1) then
           do iz=4,19
            z = 0.05*iz
            fact = 1.3 - 0.5 * (1-z)**4 / 0.7**4 -
     >                   0.2 * z**4 
            write(30+j,'(2f10.3)') z,fact
           enddo
           write(30+j,'(''join'')')
          endif
          if(jp.eq.2 .or. jp.eq.5 .or. jp.eq.6) then
c ratios from dss
           ifit=1
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             fact = fd1sv(ifit,ikin,iz)/fu1sv(ifit,ikin,iz)
             write(30+j,'(2f10.3)') z,fact
            endif
           enddo
           write(30+j,'(''join'')')
c ratios from geiger
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             call rgeiger(z,rnew)
             write(30+j,'(2f10.3)') z,rnew
            endif
           enddo
           write(30+j,'(''join dash'')')
          endif
          if(jp.eq.3.or.jp.eq.4.or.jp.eq.9.or.jp.eq.10) then
           x4 = x
           qgev = sqrt(q2)
           u = ctq5pdf(1,x4,qgev)
           d = ctq5pdf(2,x4,qgev)
! MRST fit. Note they have kappa = -0.2, but definition
! of delu is oppositie to mine
           f = 0.8 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
           f = -0.65 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
          endif
          if(jp.eq.7.or.jp.eq.8) then
            do iz=3,18
             z = 0.05 * iz
             pt=0.
             phicm=0.
             iff = 2
             call simcmodel(x,q2,z,pt,phicm,mmpi2,1 ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,2 ,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,4 ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,5 ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             if(jp.eq.7) write(30+j,'(2f10.4)') z,
     >        (sighad2 - sighad4) / (sighad1 - sighad3)
             if(jp.eq.8) write(30+j,'(2f10.4)') z,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
            enddo
            write(30+j,'(''join'')')
          endif
         endif ! j.gt.4
        enddo ! j
       enddo ! ix
      enddo ! iy
c end of 2 by 4 plots

      do ikin=1,3
       do ifit=1,3
        write(6,'(''avcsv'',2i3,4f8.3)') ikin,ifit,
     >  (avcsv(ikin,ifit,k),avcsver(ikin,ifit,k),k=1,2)
       enddo
      enddo ! ikin

      do ikin=1,8
       do ifit=1,3
        do iz=1,20
         if(mcsv5per(ifit,ikin,iz,1).ne.0.) then
           write(6,'(''m5'',3i3,10f7.2)') ikin,ifit,iz,
     >      (mcsv5p(ifit,ikin,iz,j),
     >       mcsv5per(ifit,ikin,iz,j),j=1,5)
         endif
        enddo
       enddo
      enddo ! ikin

 999  stop
      end

c modified Geiger ratio to give best agreement with our data
      subroutine rgeiger(z,rnew)
      implicit none
      real*8 z,rnew

      rnew = (1.0 -z)**0.083583 / (1.0 +z)**1.9838
c     > * (1.33 - 0.6 * z)
      return
      end


      subroutine getrho(x,q2,z,pt,phi,rplus,rminus,
     >  rplusd,rminusd)
      implicit none
      integer i,ix,iq,iz,ipt,ipm,ncall/0/
      real*8 abst(2),rat(2),asv(10,10,10,10,2)
      real*8 dp(2),dpsv(10,10,10,10,2),rplusd,rminusd
      real*8 ratsv(10,10,10,10,2),x,q2,z,pt,phi,rplus,rminus
      real*8 e0,ep,theta,nu,sin2,am/0.938/
      logical first/.true./

      if(first) then
       open(unit=11,file='rhoratio.txt')
       do i=1,9200
        read(11,'(4i3,6f8.3)') ix,iq,iz,ipt,
     >      abst(1),abst(2),rat(1),rat(2),dp(1),dp(2)
        do ipm=1,2
         asv(ix,iq,iz,ipt,ipm)=abst(ipm)
         ratsv(ix,iq,iz,ipt,ipm)=rat(ipm)
c protect against zero valuess
         if(dp(ipm).lt.0.001) dp(ipm)=1.0
         dpsv(ix,iq,iz,ipt,ipm)=dp(ipm)
        enddo
       enddo
       first=.false.
      endif

      e0 = 10.6
      nu = q2 / 2. / am / x
      ep = e0 - nu
      sin2 = q2 / 4. / e0 / ep
      theta = 57.3 * 2. * asin(sqrt(sin2))
      ix = min(10,max(1,int((ep - 3) / 4. * 10) + 1))
      iq = min(10,max(1,int((theta - 12.) / 10.0  * 10) + 1))
      iz = min(10,max(1,int((z-0.25)/0.6 * 10)+1))
      ipt = min(10,max(1,int(pt / 0.7 * 10.) + 1))

      rplus = ratsv(ix,iq,iz,ipt,1) * 
     >  (1. + asv(ix,iz,iz,ipt,1) * cos(phi)) 
      rminus = ratsv(ix,iq,iz,ipt,2) * 
     >  (1. + asv(ix,iz,iz,ipt,2) * cos(phi)) 
      rplusd = 0.
      rminusd = 0.
      if(dpsv(ix,iz,iz,ipt,1).gt.0 .and.
     >   dpsv(ix,iz,iz,ipt,2).gt.0.) then
       rplusd = rplus / dpsv(ix,iz,iz,ipt,1)
       rminusd = rminus / dpsv(ix,iz,iz,ipt,2)
      endif
c for checking
      ncall = ncall + 1
      if((ncall/100)*100.eq.ncall) then
       write(6,'(''getrho'',4i3,10f6.2)') ix,iq,iz,ipt,
     >  x,q2,ep,theta,z,pt,rplus,rminus,rplusd,rminusd
      endif

      return
      end

      subroutine mfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deld
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag

      fav = p(1)
      unf = p(2) * fav
      u = mu
      d = md
      ub = mub
      db = mdb
      s = ms
      sb = msb
c      delu = p(3) * (u + d)
c      deld = -1.* delu 
      delu = -p(3)
c      deld = -p(4)
      deld = -1./4. *delu

      sqp = 4 * (u + ub) + (d + db) + (s + sb)
c  for test: this makes delu about 50% bigger
c      sqn = 4 * (d + deld + db) + 
c     >               (u + delu + ub) + (s + sb)
c took delu out of denominator
      sqn = 4 * (d + db) + (u + ub) + (s + sb)

c change u and d to preseve incl. p and n cross sections
c took out: results in HUGE values of delu and
c more than 10% change to u (relative)
c      u = mu - delu/5.
c      d = md + 4.*delu/5.
      sqpchk = 4 * (u + ub) + (d + db) + (s + sb)
      sqnchk = 4 * (d + deld + db) + (u + delu +ub) + (s + sb)
c      If(abs(sqp - sqpchk).gt.0.001*mu .or.
c     >   abs(sqn - sqnchk).gt.0.001*mu) write(6,
c     >   '(''error sqp,sqn'',4f10.4)') sqp,sqpchk,sqn,sqnchk

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp = (4 * d + ub + 4.*deld) * fav + 
     >       (4 * db + u  + delu + s + sb) * unf 

      mnm = (4 * db + u + delu ) * fav + 
     >      (4 * d  + 4.*deld + ub  + s + sb) * unf 

      mfitr(1) = mpp/sqp
      mfitr(2) = (mpp + mnp)/(sqp + sqn)
      mfitr(3) = mpm/sqp
      mfitr(4) = (mpm + mnm)/(sqp + sqn)

      chi2 = 0.
      do i=1,4
       if(mfitflag.eq.0) then
        chi2 = chi2 + (mfit(i) -  mfitr(i))**2/mfiter(1)**2
       else
        chi2 = chi2 + (mfitp(i) - mfitr(i))**2/mfitper(1)**2
       endif
      enddo

      fval = chi2

      return
      end

      subroutine mfit5_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deld,rr,sqp0,sqn0
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag

      fav = p(1)
      unf = p(2) * fav
      u = mu
      d = md * p(5)
      ub = mub
      db = mdb
      s = ms
      sb = msb
c      delu = p(3) * (u + d)
c      deld = -1.* delu 
      delu = -p(3)
      deld = -p(4)

c values with modified d/u
      sqp = 4 * (u + ub) + (d + db) + (s + sb)
      sqn = 4 * (d + deld + db) + 
     >               (u + delu + ub) + (s + sb)

c values with original values of d
      sqp0 = 4 * (u + ub) + (md + db) + (s + sb)
      sqn0 = 4 * (md + db) + (u + ub) + (s + sb)

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp = (4 * d + ub + 4.*deld) * fav + 
     >       (4 * db + u  + delu + s + sb) * unf 

      mnm = (4 * db + u + delu ) * fav + 
     >      (4 * d  + 4.*deld + ub  + s + sb) * unf 

c fit doesn't converge with this form
      mfitr(1) = mpp/sqp
      mfitr(2) = (mpp + mnp)/(sqp + sqn)
      mfitr(3) = mpm/sqp
      mfitr(4) = (mpm + mnm)/(sqp + sqn)

c try this
      mfitr(1) = mpp/sqp0
      mfitr(2) = (mpp + mnp)/(sqp0 + sqn0)
      mfitr(3) = mpm/sqp0
      mfitr(4) = (mpm + mnm)/(sqp0 + sqn0)

      chi2 = 0.
      do i=1,4
       if(mfitflag.eq.0) then
        chi2 = chi2 + (mfit(i) -  mfitr(i))**2/mfiter(1)**2
       else
        chi2 = chi2 + (mfitp(i) - mfitr(i))**2/mfitper(1)**2
       endif
      enddo

c added contraint: p/d inclusive must remain unchanged
      rr = (sqp + sqn) / (sqp0 + sqn0)
      chi2 = chi2 + (rr-1)**2/0.05**2

      fval = chi2

      write(6,'(''mm5p'',5f7.2,f7.3,f7.1)') (p(i),i=1,5),rr,chi2

      return
      end

      subroutine getcsv(delu, deluer, deld, delder,doverp,
     >  doverpp, doverpper,fact)
! get delu and deld from sum and diff ratios assuming
c d/u, sea quarks known
      implicit none

      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deluer, deld, delder
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4),doverp,doverpp, doverpper,sqdp
      real*8 dop(5),sqnp
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag
      real*8 v(4),dv(4),sqd,mnp0,mnm0,srat,drat,drat0,C,E
      real*8 dd(5),du(5),sratp,dratp,srat0,fact

      fav = 1.
      unf = 0.5
      u = mu
c fact modifies d quark 
      d = md * fact
      ub = mub
      db = mdb
      s = ms
      sb = msb

      delu = 0.
      deld = 0. 

      sqp = 4 * (u + ub) + (d + db) + (s + sb)
      sqn = 4 * (d + db) + (u + ub) + (s + sb)
      sqd = sqp + sqn
! nominal ratio of d / p inclusive
      doverp = (sqp + sqn ) / sqp

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp0 = (4 * d + ub - 4.*deld) * fav + 
     >       (4 * db + u  - delu + s + sb) * unf 

      mnm0 = (4 * db  + u - delu ) * fav + 
     >      (4 * d  - 4.*deld + ub  + s + sb) * unf 

      do it=1,4
       v(it) = mfit(it)
       dv(it) = mfiter(it)
       if(mfitflag.ne.0) then
        v(it) = mfitp(it)
        dv(it) = mfitper(it)
       endif
      enddo

      do i=1,5
       if(i.eq.1) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1)          + v(3)          )
        drat = (v(2)          - v(4)          ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.2) then
        srat = (v(2) + dv(2)  + v(4)          ) /
     >         (v(1)          + v(3)          )
        drat = (v(2) + dv(2)  - v(4)          ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.3) then
        srat = (v(2)          + v(4) + dv(4)  ) /
     >         (v(1)          + v(3)          )
        drat = (v(2)          - v(4) - dv(4)  ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.4) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1) + dv(1)  + v(3)          )
        drat = (v(2)          - v(4)          ) /
     >         (v(1) + dv(1)  - v(3)          )
       endif
       if(i.eq.5) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1)          + v(3) + dv(3)  )
        drat = (v(2)          - v(4)          ) /
     >         (v(1)          - v(3) - dv(3)  )
       endif
       srat0 = ((sqp / sqd) * (mpp + mnp0 + mpm + mnm0) / 
     >                    (mpp        + mpm))
       C = (srat0 - srat) * (sqd / sqp) * (mpp + mpm) / (fav + unf) 
       drat0 = ((sqp / sqd) * (mpp + mnp0 - mpm - mnm0) / 
     >                    (mpp        - mpm))
       E = (drat0 - drat) * (sqd / sqp) * (mpp - mpm) / (fav - unf)
       dd(i) = (C + E)/8.
       du(i) = (C - E)/2.

! ratio with delu and deld
       sqnp = 4 * (d  - dd(i) + db) + (u - du(i) + ub) + (s + sb)
       dop(i) = (sqp + sqnp) / sqp

c check
       delu = du(i)
       deld = dd(i)
       mnp = (4 * d + ub - 4.*deld) * fav + 
     >       (4 * db + u  - delu + s + sb) * unf 
       mnm = (4 * db + u - delu ) * fav + 
     >      (4 * d  - 4.*deld + ub  + s + sb) * unf 
       dratp = ((sqp / sqd) * (mpp + mnp - mpm - mnm) / 
     >                    (mpp        - mpm))
       sratp = ((sqp / sqd) * (mpp + mnp + mpm + mnm) / 
     >                    (mpp        + mpm))
       if(abs(srat/sratp - 1.).gt.0.001 .or.
     >    abs(drat/dratp - 1.).gt.0.001) write(6,
     >    '(''csverr'',8f8.3)') delu,deld,
     <   drat,drat0,dratp,srat,srat0,sratp

      enddo

      delu = du(1)
      deld = dd(1)
      deluer = sqrt((du(1) - du(2))**2 + 
     >              (du(1) - du(3))**2 +
     >              (du(1) - du(4))**2 +
     >              (du(1) - du(5))**2)
      delder = sqrt((dd(1) - dd(2))**2 + 
     >              (dd(1) - dd(3))**2 +
     >              (dd(1) - dd(4))**2 +
     >              (dd(1) - dd(5))**2)

      doverpp = dop(1)
      doverpper = sqrt((dop(1) - dop(2))**2 + 
     >                (dop(1) - dop(3))**2 +
     >                (dop(1) - dop(4))**2 +
     >                (dop(1) - dop(5))**2)


      return
      end

      subroutine sfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 u,d,rpp,rnp,rpm,rnm,b1,b2,ds, w,mmpi2,w2
      real*8 qu,qd,qs,nu,c0,c1,c2,c3,c4,dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 sum_sqp,sum_sqn,fu,fub,fd,fdb,fs,fsb
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/,pi/3.141593/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      integer bkin(90000),ncall
      real*8 rhofact,fnorm(36)
      real*8 bchi2x(15),bdfx(15)
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt,ncall
        integer ikinfit,itfit,npts
        real*8 zfit
        common/sstuff/ zfit,ikinfit,itfit,npts
        integer sn
        real*8 sphi(1000), spt(1000),sy(10000),syer(1000)
        common/sforplot/ sphi,spt,sy,syer,sn

        ncall = ncall + 1
        chi2 = 0.
        npts = 0
        sn = 0

      do i=1,bnpt
       if(bitv(i).eq.itfit .and. bkin(i).ge.ikinfit .and.
     >   bkin(i).le.ikinfit + 1.and.
     >   abs(bzv(i) - zfit).lt.0.025) then
        npts = npts + 1
        sn = sn + 1
        phi = bphiv(i)
        sphi(sn) = phi
c        if(bkin(i).eq.ikinfit) sphi(sn) = sphi(sn)+0.05
c        if(bkin(i).eq.ikinfit+1) sphi(sn) = sphi(sn)-0.05
        pt = bptv(i)
        spt(sn) = pt
        pt2 = pt**2
        z = bzv(i)
c using p1 / z
        b = 1./p(2)
ccc        sig = (p(1)/z  ) * 2.*sqrt(b/pi) * exp(-b *pt2) *
        sig = (p(1)/z  ) * b * exp(-b *pt2) / 2. / pi **
     >  (1. + p(3) * pt * cos(phi))
c     >  (1. + p(3) * pt * cos(phi)  + 
c     >        p(4) * pt2 * cos(2.*phi))
        bsigv(i) = sig
! how much rho subtraction
        y = byv(i)  - rhofact * brho(i)
        sy(sn) = y
        syer(sn) = byerv(i)
! special code to use brhop. Not currently in use, I hope
        if(rhofact.eq.10.) y = byv(i)  - brhop(i)
        if(ncall.eq.1 .and. ikinfit.eq.1 .and. zfit.gt.0.5
     >   .and. zfit.lt.0.55) then
         write(6,'(''sdbg'',f4.1,i3,3f8.3)') rhofact,sn,y,byerv(i),
     >    rhofact * brho(i)
        endif
        chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       endif
      enddo

      if((ncall/100)*100.eq.ncall .and. ncall.gt.0 .and.
     >  npts.gt.0)
     >  write(6,'(''ncall small, chi2'',i8,i6,2i3,f6.2,f8.3)')
     >  ncall,npts,ikinfit,itfit,zfit,chi2/float(npts)

      fval = chi2

      return
      end

      subroutine bfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 u,d,rpp,rnp,rpm,rnm,b1,b2,ds, w,mmpi2,w2
      real*8 qu,qd,qs,nu,c0,c1,c2,c3,c4,dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 sum_sqp,sum_sqn,fu,fub,fd,fdb,fs,fsb
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      integer bkin(90000),ncall/0/
      real*8 rhofact,fnorm(36)
      real*8 bchi2x(15),bdfx(15)
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt

      ncall = ncall + 1
      chi2 = 0.
      bchi2k(0)=0.
      bdfk(0)=0.
      do i=1,56
       bchi2k(i)=0.
       bdfk(i)=0.
       if(i.le.15) bchi2x(i)=0.
       if(i.le.15) bdfx(i)=0.
       if(i.le.15) bchi2q2(i)=0.
       if(i.le.15) bdfq2(i)=0.
       if(i.le.15) bchi2w(i)=0.
       if(i.le.15) bdfw(i)=0.
       if(i.le.15) bchi2z(i)=0.
       if(i.le.15) bdfz(i)=0.
       if(i.le.15) bchi2f(i)=0.
       if(i.le.15) bdff(i)=0.
       if(i.le.15) bchi2pt(i)=0.
       if(i.le.15) bdfpt(i)=0.
       if(i.le.15) bchi2mx(i)=0.
       if(i.le.15) bdfmx(i)=0.
       if(i.le.6) bchi2t(i)=0.
       if(i.le.6) bdft(i)=0.
      enddo

      do i=1,bnpt
       it = bitv(i)
       u = buv(i) 
       d = bdv(i)
       s = bsv(i)
       ubar = bubv(i)
       dbar = bdbv(i)
       sbar = bsbv(i)
       delu = p(19) * (1 - x)**4 / sqrt(x) * (x-0.0909)
       deld = -1. * delu
       sum_sqp = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)
       sum_sqn = qu**2*(d + deld + dbar) + 
     >         qd**2*(u + delu + ubar) + qs**2*(s+sbar)
       sum_sq = sum_sqp
       if (it.eq.2 .or. it.eq.5) then
        sum_sq = sum_sqp + sum_sqn
       endif
       if(bkin(i).gt.0.) then
       phi = bphiv(i)
       pt = bptv(i)
       pt2 = pt**2
       q2 = bq2v(i)
       x = bxv(i)
       mmpi2 = bmmpi2(i)
       w = sqrt(0.938**2 + q2 * (1./x -1))
       w2 = w**2
       nu = q2 / 2. / 0.938 / x
       y = nu / 10.6
       xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
       z0 = bzv(i)
       zp = (z0 / 2.) * (xp / x) *(1. + 
     >     sqrt(1 - 4 * x**2 * am**2 *  
     >     (ampi**2 + pt2) / z0**2 / q2**2))
c       z = z0
c use zp instead of z
       z = zp
       z2 = z**2
c       if(i .lt. 100 .and. ncall.eq.1)
c     >   write(6,'(10f7.3)') xp/x,zp/z0,x,xp,z0,z,zp,q2,am,ampi

c start off with Binneweis frag. func. 
c        sv = log( log(Q2/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
c        N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
c        a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
c        a2 = 1.430 + 0.543*sv - 0.023*sv**2
c        N = n * p(5)
c        a1 = a1 * p(6)
c        a2 = a2 * p(7)
c	Ns = 4.250 - 3.147*sv + 0.755*sv**2
c	a1s = -0.770 -0.573*sv + 0.117*sv**2
c	a2s = 4.48 + 0.890*sv - 0.138*sv**2
cxx	zhadm = min(z, 0.75)
c        zhadm = min(z, 0.99)
c	D_sum = N * zhadm**a1 * (1.0-zhadm)**a2
c	D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
c        b1 = 0.083583 * p(8)
c        b2 = 1.983 * p(9)
c	R_D = (1.0-zhadm)**b1 / (1.0+zhadm)**b2
c	D_fav = D_sum/(1.0+R_D)
c	D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
c	D_s = D_sum_s/2.0
c        Dp = d_fav
c        Dm = D_unfav
c new way using DSS
       fu = bffu(i) * (1. + p(5)/w2 +
     >  p(6)/mmpi2 + p(7)/q2)
       fdb= bffdb(i) * (1. + p(5)/w2 +
     >  p(6)/mmpi2 + p(7)/q2)
       fd = bffd(i) * (1. + p(8)/w2 +
     >  p(9)/mmpi2 + p(10)/q2)
       fub= bffub(i) * (1. + p(8)/w2 +
     >  p(9)/mmpi2 + p(10)/q2)
       fs = bffs(i)
       fsb = bffsb(i)
       c1 = 1. / (z2*p(1) + p(3)) / 9.
       c2 = 1. / (z2*p(2) + p(4)) / 9.
       c3 = 1. / (z2*p(1) + p(4)) / 9.
       c4 = 1. / (z2*p(2) + p(3)) / 9.

c need to add strange sea!
       rpp = 4.*u*fu * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >          d*fd * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >       dbar*fdb* c1*exp(-pt2/(z2*p(1) + p(3)))+
     >    4.*ubar*fub* c2*exp(-pt2/(z2*p(2) + p(4)))
     >   + (s + sbar)*fs * c2*exp(-pt2/(z2*p(2) + p(4)))

       rpm = 4.*u*fub* c3*exp(-pt2/(z2*p(1) + p(4)))+
     >          d*fdb* c4*exp(-pt2/(z2*p(2) + p(3)))+
     >       dbar*fd * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >    4.*ubar*fu * c4*exp(-pt2/(z2*p(2) + p(3)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))
       rnp = 4.*(d+deld)*fu * c4*exp(-pt2/(z2*p(2) + p(3)))+
     >          (u+delu)*fd * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >       ubar*fdb* c4*exp(-pt2/(z2*p(2) + p(3)))+
     >    4.*dbar*fub* c3*exp(-pt2/(z2*p(1) + p(4)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))
       rnm = 4.*(d+deld)*fub* c2*exp(-pt2/(z2*p(2) + p(4)))+
     >          (u+delu)*fdb* c1*exp(-pt2/(z2*p(1) + p(3)))+
     >       ubar*fd * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >    4.*dbar*fu * c1*exp(-pt2/(z2*p(1) + p(3)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))

       if(it.eq.1) sig = rpp * 
     >   (1. + p(11)*sqrt(pt2/q2)*cos(phi) +
     >       p(12)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.2) sig = (rpp + rnp) *
     >   (1. + p(13)*sqrt(pt2/q2)*cos(phi) +
     >       p(14)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.4) sig = rpm *
     >   (1. + p(15)*sqrt(pt2/q2)*cos(phi) +
     >       p(16)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.5) sig = (rpm + rnm) * 
     >   (1. + p(17)*sqrt(pt2/q2)*cos(phi) +
     >       p(18)*(pt2/q2)*cos(2.*phi)) 

       sig = sig / sum_sq / 2. / 3.1415928

       bsigv(i) = sig

! how much rho subtraction
       y = byv(i)  - rhofact * brho(i)
c normalization factor by ikin
       y = y * fnorm(bkin(i))
c inclusive ratio here
       else
        sig = 1. + sum_sqn / sum_sqp
       endif
       bsigv(i) = sig
       chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       j=bkin(i)
       bchi2k(j) = bchi2k(j) + (sig - y)**2 / byerv(i)**2
       bdfk(j) = bdfk(j) + 1.
       if(j.ne.0) then
        j = min(15,max(1,int((x-0.2)/0.6*15)+1))
        bchi2x(j) = bchi2x(j) + (sig - y)**2 / byerv(i)**2
        bdfx(j) = bdfx(j) + 1.
        j = min(15,max(1,int((w-1.6)/1.8*15)+1))
        bchi2w(j) = bchi2w(j) + (sig - y)**2 / byerv(i)**2
        bdfw(j) = bdfw(j) + 1.
        j = min(15,max(1,int((q2-1.)/6.*15)+1))
        bchi2q2(j) = bchi2q2(j) + (sig - y)**2 / byerv(i)**2
        bdfq2(j) = bdfq2(j) + 1.
        j = min(15,max(1,int(phi/2./3.1415928*15)+1))
        bchi2f(j) = bchi2f(j) + (sig - y)**2 / byerv(i)**2
        bdff(j) = bdff(j) + 1.
        j = min(15,max(1,int((z-0.2)/0.75*15)+1))
        bchi2z(j) = bchi2z(j) + (sig - y)**2 / byerv(i)**2
        bdfz(j) = bdfz(j) + 1.
        j = min(15,max(1,int(pt/0.95*15)+1))
        bchi2pt(j) = bchi2pt(j) + (sig - y)**2 / byerv(i)**2
        bdfpt(j) = bdfpt(j) + 1.
        j = min(15,max(1,int((sqrt(mmpi2)-1.4)/1.4*15)+1))
        bchi2mx(j) = bchi2mx(j) + (sig - y)**2 / byerv(i)**2
        bdfmx(j) = bdfmx(j) + 1.
        j = it
        bchi2t(j) = bchi2t(j) + (sig - y)**2 / byerv(i)**2
        bdft(j) = bdft(j) + 1.
       endif
      enddo

      if((ncall/100)*100.eq.ncall .and. ncall.gt.0) 
     >  write(6,'(''ncall, chi2'',i8,f8.3)')
     >  ncall,chi2/float(bnpt)

      fval = chi2

      return
      end

      subroutine bfit_fcnOLD(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 fm,p1p,p2p,p3p,p4p,u,d,rpp,rnp,rpm,rnm,rppp,rpmp,b1,b2
      real*8 qu,qd,qs,nu,er,c0,c1,c2,c3,c4,c1p,c2p,c3p,c4p, dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000)
      real*8 byv(90000),byerv(90000)
      real*8 buv(90000),bubv(90000),brho(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(56),bdfk(56)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      integer bkin(90000),ncall/0/
      real*8 rhofact
      common/bstuffold/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,rhofact,bnpt

      ncall = ncall + 1
      chi2 = 0.
      do i=1,56
       bchi2k(i)=0.
       bdfk(i)=0.
      enddo

      fm = 0.005
      p1p = p(1) + fm
      p2p = p(2) + fm
      p3p = p(3) + fm
      p4p = p(4) + fm

      do i=1,bnpt
       it = bitv(i)
       phi = bphiv(i)
       pt = bptv(i)
       pt2 = pt**2
       z0 = bzv(i)
       q2 = bq2v(i)
       x = bxv(i)
       er = berv(i) ! epsilon * R=sigl/sigt
       nu = q2 / 2. / 0.938 / x
       y = nu / 10.6
       u = buv(i) 
       xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
       zp = (z0 / 2.) * (xp / x) *(1. + 
     >     sqrt(1 - 4 * x**2 * am**2 *  
     >     (ampi**2 + pt2) / z0**2 / q2**2))
       z = z0
c use zp instead of z
       z = zp
       z2 = z**2
       if(i .lt. 100 .and. ncall.eq.1)
     >   write(6,'(10f7.3)') xp/x,zp/z0,x,xp,z0,z,zp,q2,am,ampi
c old use of p10
c       d = bdv(i) * p(10)
       d = bdv(i)
       s = bsv(i)
       ubar = bubv(i)
       dbar = bdbv(i)
       sbar = bsbv(i)
       sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)
       if (it.eq.2 .or. it.eq.5) then
          sum_sq = sum_sq + qu**2*(d+dbar) + 
     >         qd**2*(u+ubar) + qs**2*(s+sbar)
       endif

       c0 = -4. * (2. - y) * sqrt(1. - y) * z / sqrt(q2) /  
     >      (1. + (1. - y)**2) * sqrt(pt2) * cos(phi)
       c0 = c0 * p(10)
cxx turned Cahn teerm off!
       c0 = 0.
       c1 = 1. + c0 * p(1) / (z2*p(1) + p(3))
       c2 = 1. + c0 * p(2) / (z2*p(2) + p(4))
       c3 = 1. + c0 * p(1) / (z2*p(1) + p(4))
       c4 = 1. + c0 * p(2) / (z2*p(2) + p(3))
       c1 = c1 / (z2*p(1) + p(3)) / 9.
       c2 = c2 / (z2*p(2) + p(4)) / 9.
       c3 = c3 / (z2*p(1) + p(4)) / 9.
       c4 = c4 / (z2*p(2) + p(3)) / 9.
c these are for deuteron (have Fermi motion added)
       c1p = 1. + c0 * p1p / (z2*p1p + p(3))
       c2p = 1. + c0 * p2p / (z2*p2p + p(4))
       c3p = 1. + c0 * p1p / (z2*p1p + p(4))
       c4p = 1. + c0 * p2p / (z2*p2p + p(3))
       c1p = c1p / (z2*p1p + p(3)) / 9.
       c2p = c2p / (z2*p2p + p(4)) / 9.
       c3p = c3p / (z2*p1p + p(4)) / 9.
       c4p = c4p / (z2*p2p + p(3)) / 9.

c start off with Binneweis frag. func. 
        sv = log( log(Q2/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
        N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
        a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
        a2 = 1.430 + 0.543*sv - 0.023*sv**2
        N = n * p(5)
        a1 = a1 * p(6)
        a2 = a2 * p(7)
	Ns = 4.250 - 3.147*sv + 0.755*sv**2
	a1s = -0.770 -0.573*sv + 0.117*sv**2
	a2s = 4.48 + 0.890*sv - 0.138*sv**2
cxx	zhadm = min(z, 0.75)
        zhadm = min(z, 0.99)
	D_sum = N * zhadm**a1 * (1.0-zhadm)**a2
	D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
        b1 = 0.083583 * p(8)
        b2 = 1.983 * p(9)
	R_D = (1.0-zhadm)**b1 / (1.0+zhadm)**b2
	D_fav = D_sum/(1.0+R_D)
	D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
	D_s = D_sum_s/2.0
        Dp = d_fav
        Dm = D_unfav

c need to add sea!
       rpp = 4.*u*Dp * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >          d*Dm * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >       dbar*Dp * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >    4.*ubar*Dm * c2*exp(-pt2/(z2*p(2) + p(4)))
       rpm = 4.*u*Dm * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >          d*Dp * c4*exp(-pt2/(z2*p(2) + p(3)))+
     >       dbar*Dm * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >    4.*ubar*Dp * c4*exp(-pt2/(z2*p(2) + p(3)))
       rppp= 4.*u*Dp * c1p*exp(-pt2/(z2*p1p + p(3)))+
     >          d*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))+
     >       dbar*Dp * c1p*exp(-pt2/(z2*p2p + p(3)))+
     >    4.*ubar*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))
       rpmp= 4.*u*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >          d*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))+
     >       dbar*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >    4.*ubar*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))
       delu = p(19) * (1 - x)**4 / sqrt(x) * (x-0.0909)
       deld = -1. * delu
       rnp = 4.*(d+deld)*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))+
     >          (u+delu)*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >       ubar*Dp * c1p*exp(-pt2/(z2*p2p + p(3)))+
     >    4.*dbar*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))
       rnm = 4.*(d+deld)*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))+
     >          (u+delu)*Dp * c1p*exp(-pt2/(z2*p1p + p(3)))+
     >       ubar*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >    4.*dbar*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))

       if(it.eq.1) sig = rpp * 
     >   (1. + p(11)*sqrt(pt2/q2)*cos(phi) +
     >       p(12)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.2) sig = (rppp + rnp) *
     >   (1. + p(13)*sqrt(pt2/q2)*cos(phi) +
     >       p(14)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.4) sig = rpm *
     >   (1. + p(15)*sqrt(pt2/q2)*cos(phi) +
     >       p(16)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.5) sig = (rpmp + rnm) * 
     >   (1. + p(17)*sqrt(pt2/q2)*cos(phi) +
     >       p(18)*sqrt(pt2/q2)*cos(2.*phi)) 

       sig = sig / sum_sq / 2. / 3.1415928

! eps * R term
       sig = sig * (1. + er)
       bsigv(i) = sig

! new use of p(10): how much rho subtraction
       y = byv(i) *(1. - p(10)*brho(i))
       chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       j=bkin(i)
       bchi2k(j) = bchi2k(j) + (sig - y)**2 / byerv(i)**2
       bdfk(j) = bdfk(j) + 1.
      enddo

      if((ncall/100)*100.eq.ncall) write(6,'(''ncall, chi2'',i8,f8.3)')
     >  ncall,chi2/float(bnpt)

      fval = chi2

      return
      end

      subroutine simcmodel(x,q2,z,pt,phi,mmpi2,
     >  it,sighad,u,ubar,d,dbar,u1,d1,
     >  s,sbar,ff,fu,fs,dsigdz,yf,yu,zp,rfu,iFF,ipdf)
      implicit none

      real*8 x,z,pt,phi,wsq,w
      integer it,ncall,iFF,ipdf

	integer iset  !which set (1=cteq5m)
	integer ipart !particle u=1, ubar=-1, d=2, dbar=-2, s=3, sbar=-3
	real*8 u,d,ubar,dbar,s,sbar,fav,unfav,rgeiger
	real*8 qu,qd,qs ! u, d, and s quark charges
	real*8 D_fav, D_unfav, D_sum, R_D !favored,unfavored,sum,ratio of FFs
	real*8 D_sum_s, D_s  ! strange frag. functions
	real*8 lambda, Q2zero ! scales for FF param
	real*8 sv !scaling variable for FF param
	real*8 N,a1,a2 !parameters for FF param
	real*8 Ns, a1s, a2s !parameters for strange FF param
        real*8 ff,fu,fs,chk,rfu,q28
C Some local kinematic variables
        real*4 xbj,q2gev, Qgev
	real*8 sx, pt2gev !unitless or GeV
	real*8 b  ! pt2 parameter for FFs

	real*8 nu,qx,qy,qz,mtar,Q2,Eb,Eprime,Epx,Epy,Epz  !MeV
	real*8 pt2,zhad,Ehad,phad,mhad,zhadm ! all in MeV
	real*8 cthpq,pi/3.1415928/

	real*8 kcent,klo,khi,mmpi2,mtargev,nugev
	integer i,nwritten/0/,nprint/0/

	real*8 sum_sq, dsigdz, sigsemi, jacobian, fac, sigma_eepiX
	real*8 sighad, sige, dsigdzold,dsigdzp,dsigdzn
	real*8 dsigdzpold, dsigdznold

	real*8 F1,F2,W1,W2,sin2th2,cos2th2,W2coeff
	real*8 Ctq5Pdf	

c parameters for PB fit of 9/11/2020
c versus zp
c       real*8 pf(8)/   1.2803,   0.0851,   0.8379,   0.1586,
c     >                 0.0140,   0.2133,  -4.4985,   4.1285/
c       real*8 pu(8)/   0.8290,  -0.1416,   0.9869,   0.2559,
c     >                 0.0090,  -1.2306,  -1.5292,   2.4169/

c from fit of 12/5/2020
c      real*8 pf(8)/ 2.1676,   0.4117,   1.4852,   0.1480,
c     >              0.1406,   1.6525,  -8.9585,   8.1429/
c      real*8 pu(8)/ 2.1322,   0.3988,   1.9398,   0.1897,
c     >              0.1337,  -0.7782,  -3.2305,   4.3361/

c from fit versus zp of 12/15/2020, which uses
c z-dependant values of b as below
c      real*8 pf(8)/2.4367,   0.6459,   1.6099,   0.1529,
c     >             0.1439,   1.8692,  -9.9865,   9.3691/
c      real*8 pu(8)/1.7180,   0.3416,   1.8522,   0.2327,
c     >             0.1332,  -1.1574,  -2.0046,   3.4474/
c fit of June 20, 2021
      real*8 pf(12)/   1.0424,  -0.1714,   1.8960,  -0.0307,
     >                 0.1636,  -0.1272,  -4.2093,   5.0103,
     >                 2.7406,  -0.5778,   3.5292,   7.3910/
      real*8 pu(12)/   0.7840,   0.2369,   1.4238,   0.1484,
     >                 0.1518,  -1.2923,  -1.5710,   3.0305,
     >                 1.1995,   1.3553,   2.5868,   8.0666/

	logical first, doinghkns

	parameter (iset=1)
	parameter (qu=2./3.)
	parameter (qd=-1./3.)
	parameter (qs=-1./3.)

	parameter (lambda=0.227) !0.227 GeV for NLO
	parameter (Q2zero=2.0)   !Gev^2 for u,d,s,g

	data first /.TRUE./

c this is for DSS
        integer IHdss,ICdss,IOdss, fini 
        real*8  U1, UB, D1, DB, S1, SB, C1, B1, GL1
        COMMON / FRAGINI / FINI

! this is for HKNS frag. fun.
	integer ixx,iqq,iswap
        real*8 hkusv(100,100),hkdsv(100,100)
        real*8 xp,zp,yu,yf
c parameter from March 2022 fit
        real*8 p(10),ffav,funf,am/0.938/,ampi/0.139/,mpi/0.139/
c
cf  5  0.713  0.061
cf  6  0.455  0.024
cf  7 -0.192  0.018
cf  8 -1.582  0.122
cf  9  1.208  0.056
cf 10  0.791  0.044
       p(5)=  0.715
       p(6)=  0.455
       p(7)= -0.192
       p(8)= -1.582
       p(9)=  1.208
       p(10)=  0.791

        if(first) fini=0

        if(first) then
         open(unit=23,file='FFHKNS07.dat')
         do iqq=1,100
          do ixx=1,100
           read(23,'(f5.1,f5.2,2f9.5)') q2gev,zhad,
     >       hkusv(iqq,ixx),hkdsv(iqq,ixx)
          enddo
         enddo
         close(unit=23)
        endif
	if(first) then
	   call SetCtq5(iset)	! initialize Cteq5 (we're using cteq5m)
	   first=.FALSE.
	endif

	b = 3.8
	mhad = Mpi
	doinghkns = .true.
        xbj = x
        q2gev = q2
	Qgev = sqrt(Q2gev)
        zhad = z
        pt2gev = pt**2

! using xp has almost no effet compared to x
        xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
        xbj = xp

        if(ipdf.eq.1) then
 	ipart=1
	u = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-1
	ubar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=2
	d = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-2
	dbar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=3
	s = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-3
	sbar = Ctq5pdf (ipart , xbj, Qgev)

        else
c         call getjam(xp, q2, u,ubar,d,dbar,s,sbar)
         call getjam(x, q2, u,ubar,d,dbar,s,sbar)
c for test, make d/u smaller at higher x
c this makes p+ / d+ better, but p- / d+ worse!
c         if(x.gt.0.35) d = 0.
        endif

	sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)

	if (it.eq.2 .or. it.eq.5) then
	   sum_sq = sum_sq + qu**2*(d+dbar) + 
     >       qd**2*(u+ubar) + qs**2*(s+sbar)
	endif

	sv = log( log(Q2gev/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
	   N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
	   a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
	   a2 = 1.430 + 0.543*sv - 0.023*sv**2
	   Ns = 4.250 - 3.147*sv + 0.755*sv**2
	   a1s = -0.770 -0.573*sv + 0.117*sv**2
	   a2s = 4.48 + 0.890*sv - 0.138*sv**2
   	   zhadm = min(zhad, 0.75)
cxx	   zhadm = zhad
	   D_sum = N*zhadm**a1*(1.0-zhadm)**a2
c correction to get better agreement with data
	   D_sum = D_sum * (1.5 - 0.5 * zhad)
	   D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
	   R_D = (1.0-zhad)**0.083583/(1.0+zhad)**1.9838
	   R_D = (1.0-zhadm)**0.083583/(1.0+zhadm)**1.983
           rgeiger = r_d
c New from Peter's fit of August 27,2020
           R_D = (1. - 0.8 * zhad) * (1. + 
     >      zhad    * 1.988 + 
     >      zhad**2 * (-9.014) +
     >      zhad**3 * 8.219)
	   D_fav = D_sum/(1.0+R_D)
	   D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
	   D_s = D_sum_s/2.0

c new PB it using zp for pions. This is z * D
        wsq = am**2 + q2 * (1/x -1)
        w = sqrt(wsq) 
         xp = 2.*xbj / (1. + sqrt(1. + 4. * xbj**2 * am**2 / q2gev))
         zp = (zhad / 2.) * (xp / xbj) *(1. + 
     >     sqrt(1 - 4 * xbj**2 * am**2 *  
     >     (ampi**2 + pt2gev) / zhad**2 / q2gev**2))
         sv = log(q2gev/2.)
c         yf = pf(1) * zp**(pf(2) + pf(4)*sv) * 
c     >        (1.-zp)**(pf(3) + pf(5)*sv) 
c         yf = yf * (1. + pf(6)*zp + pf(7)*zp**2 + pf(8)*zp**3)
c         yu = pu(1) * zp**(pu(2) + pu(4)*sv) * 
c     >        (1.-zp)**(pu(3) + pu(5)*sv) 
c         yu = yu * (1. + pu(6)*zp + pu(7)*zp**2 + pu(8)*zp**3)
       yf = pf(1) * zp**(pf(2) + pf(4)*sv + pf(9)/w) * 
     >       (1.-zp)**(pf(3) + pf(5)*sv + pf(10)/w) 
       yf = yf * 
     >      (1. + pf(6)*zp + pf(7)*zp**2 + pf(8)*zp**3) *
     >      (1. + pf(11)/w + pf(12)/w**2)
       yu = pu(1) * zp**(pu(2) + pu(4)*sv + pu(9)/w) * 
     >       (1.-zp)**(pu(3) + pu(5)*sv + pu(10)/w) 
       yu = yu * 
     >      (1. + pu(6)*zp + pu(7)*zp**2 + pu(8)*zp**3) *
     >      (1. + pu(11)/w + pu(12)/w**2)

         if(it.le.3) then
c forget about hkns for now
c 	  u1 = D_fav * zhad
c	  d1 = D_unfav * zhad
c new fit
	  u1 = yf
	  d1 = yu
	 else
c	  d1 = zhad * D_fav
c	  u1 = zhad * D_unfav
c new fit
	  u1 = yu
	  d1 = yf
	 endif
         rfu = yu / yf

	 ub = d1
	 db = u1
c	 s1 = zhad * D_s
c	 sb = s1
c new fit
	 s1 = yu
	 sb = s1

c override with JAM
        if(iFF.eq.4) then
         q28 = q2gev
         IHdss=1
         ICdss=1
         if(it.gt.3) ICdss = -1
         IOdss=0
! use fDSS for sea quarks
         call fDSS (IHdss, ICdss, IOdss, Zp, Q28, 
     >       U1, UB, D1, DB, S1, SB, C1, B1, GL1)
! JAM for u1, d1
         if(it.le.3) then
          call jamff(zp,q28,d1,u1) ! pi+
         else
          call jamff(zp,q28,u1,d1) ! pi-
         endif
        endif

c use DSS
        if(iFF.eq.2.or.iff.eq.3) then
         q28 = q2gev
         IHdss=1
         ICdss=1
         if(it.gt.3) ICdss = -1
         IOdss=0
         call fDSS (IHdss, ICdss, IOdss, Zp, Q28, 
     >       U1, UB, D1, DB, S1, SB, C1, B1, GL1)
        endif
        iswap=0
        if(iff.eq.3) iswap=1

c put in first-pass corr. from fit
        if(it.le.3) then
         ffav = 1. + p(5)/wsq + p(6)/mmpi2 + p(7)/q2
         funf = 1. + p(8)/wsq + p(9)/mmpi2 + p(10)/q2
        else
         funf = 1. + p(5)/wsq + p(6)/mmpi2 + p(7)/q2
         ffav = 1. + p(8)/wsq + p(9)/mmpi2 + p(10)/q2
        endif

        if(ncall.lt.10) then
         write(6,'(''ffav'',i2,6f7.3)') it,wsq,mmpi2,
     >    q2,ffav,funf
         ncall = ncall + 1
        endif

cturn off for test
        ffav = 1.
        funf = 1.
           
c for test. Doesnt fix pi+/d
c        if(it.le.3) d1 = d1 * 1.15
c        if(it.gt.3) u1 = u1 * 1.15
        dsigdz = (qu**2 * u    * u1 * ffav + 
     >            qu**2 * ubar * ub * funf +
     >  	  qd**2 * d    * d1 * funf + 
     >            qd**2 * dbar * db * ffav + 
     >  	  qs**2 * s    * s1 + 
     >            qs**2 * sbar * sb)/sum_sq/zhad
        if(it.le.3) then
 	 ff = (qu**2 * u + qd**2 * dbar)/sum_sq/zhad
	 fu = (qd**2 * d + qu**2 * ubar)/sum_sq/zhad
        else
 	 fu = (qu**2 * u + qd**2 * dbar)/sum_sq/zhad
	 ff = (qd**2 * d + qu**2 * ubar)/sum_sq/zhad
        endif
	fs = (qs**2 * s + qs**2 * sbar)/sum_sq/zhad
	dsigdzp = dsigdz
	if(it.eq.2.or.it.eq.5) then ! Deut
         if(iswap.ne.1) then
         dsigdzn =(qu**2 * d    * u1 * ffav + 
     >             qu**2 * dbar * ub * funf +
     >  	   qd**2 * u    * d1 * funf + 
     >             qd**2 * ubar * db * ffav + 
     >  	   qs**2 * s    * s1 + 
     >             qs**2 * sbar * sb)/sum_sq/zhad
         else
         dsigdzn =(qu**2 * d    * d1 * ffav + 
     >             qu**2 * dbar * ub * funf +
     >  	   qd**2 * u    * u1 * funf + 
     >             qd**2 * ubar * db * ffav + 
     >  	   qs**2 * s    * s1 + 
     >             qs**2 * sbar * sb)/sum_sq/zhad
         endif
	 dsigdz = dsigdzp + dsigdzn
         if(it.le.3) then
    	 ff = ff + (qu**2 * d + qd**2 * ubar)/sum_sq/zhad
	 fu = fu + (qd**2 * u + qu**2 * dbar)/sum_sq/zhad
         else
    	 fu = fu + (qu**2 * d + qd**2 * ubar)/sum_sq/zhad
	 ff = ff + (qd**2 * u + qu**2 * dbar)/sum_sq/zhad
         endif
	 fs = fs +(qs**2 * s + qs**2 * sbar)/sum_sq/zhad
        endif
 	chk = ff*u1 + fu*d1 + fs*s1
        nprint = nprint+1
        if(nprint.lt.10) then
         write(6,'(/''chk'',9f7.3)') 
     >    chk/dsigdz,chk,dsigdz,
     >    ff,fu,fs,u1,d1,s1    
         write(6,'(i2,9f7.3)') 
     >    it,xbj,q2gev,zhad,u,ubar,d,dbar,s,sbar    
        endif

c b as function of z
        b = 1./ (0.200 * zhad**2 + 0.200)
c used in SIMC
	sighad = dsigdz * b * exp(-b * pt2gev)/ 2. / pi

      return
      end
C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.

C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      real*8 Function Ctq5Pdf (Iparton, X, Q)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      real*8 partonx
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
        Print *, 'Setting to Alambda'
        Q=Alambda
c	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

      real*8 FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
c      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      implicit integer (I-N)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
c         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
c     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq5 (Iset)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      Parameter (Isetmax=9)

      Character Flnm(Isetmax)*17, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5/cteq5m.tbl', 'cteq5/cteq5d.tbl', 'cteq5/cteq5l.tbl', 
     >   'cteq5/cteq5hj.tbl'
     > , 'cteq5/cteq5hq.tbl', 'cteq5/cteq5f3.tbl', 'cteq5/cteq5f4.tbl'
     > , 'cteq5/cteq5m1.tbl', 'cteq5/ctq5hq1.tbl'  /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
C                             ********************
      End

      Subroutine ReadTbl (Nu)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
      integer N, NextUn
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
c      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) stop
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      real*8 function futil()
      return
      end

********************************************************************
*                                                                  *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
*  D.de Florian, R.Sassot, M.Stratmann   Phys.Rev.D75 114010 2007  *
*                                 *and*  Phys.Rev.D76 074033 2007  *
*                                                                  *
*     CALL fDSS (IH,IC,IO, X, Q2, U, UB, D, DB, S, SB, C, B, GL)   *
*                                                                  *	
*  INPUT:                                                          *
*  IH = hadron type    1: PION                                     *
*                      2: KAON                                     *
*                      3: PROTON                                   *
*                      4: CHARGED HADRONS                          *
*                                                                  *
*  IC = Hadron Charge  0: 0 (as average of + and -)                *
*                      1: +                                        *
*                     -1: -                                        *
*                                                                  *
*  IO= Order           0: LO                                       *
*                      1: NLO                                      *
*                                                                  *
*            X                    (between  0.05   and  1.0)       *
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
*             (for values outside the allowed range the program    *
*              writes a warning and extrapolates to the x and      *
*              Q2 values requested)                                *
*                                                                  *
*   OUTPUT: U, UB, D, DB, S, SB,   C,           B,       GL        *
*           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon      *
*           Always X times the distribution is returned            *
*                                                                  *
*                                                                  *
*   COMMON:  The main program or the calling routine has to have   *
*            a common block  COMMON / FRAGINI / FINI , and  FINI   *
*            has always to be zero when DSS is called for the      *
*            first time or when the SET has been changed.          *
*                                                                  *
********************************************************************

      SUBROUTINE fDSS (IH,IC,IO, X, Q2, U, UB, D, DB, 
     >  S, SB, C, B, GL)
      implicit real*8 (A-H,O-Z)
      PARAMETER (NPART=9, NX=35, NQ=24, NARG=2)
      DIMENSION XUTOTF(NX,NQ), XDTOTF(NX,NQ), XSTOTF(NX,NQ)
      DIMENSION XUVALF(NX,NQ), XDVALF(NX,NQ), XSVALF(NX,NQ)
      DIMENSION XCTOTF(NX,NQ), XBTOTF(NX,NQ)
      DIMENSION XGF(NX,NQ), PARTON (NPART,NQ,NX-1)
      DIMENSION QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      integer fini
      COMMON / FRAGINI / FINI
      SAVE XUTOTF, XDTOTF, XSTOTF, XCTOTF, XBTOTF, XGF, NA, ARRF
      SAVE XUVALF, XDVALF, XSVALF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.d0, 1.25D0, 1.5D0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
       DATA XB /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6        0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9 , 0.93, 1.0/
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.05D0) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C          STOP
       ENDIF
       IF ( (Q2.LT.1.D0) .OR. (Q2.GT.1.D5) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
C          STOP
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
      IF (FINI.NE.0) GOTO 16
      IF ((IH.EQ.1).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='PINLO.GRID')
      ELSEIF ((IH.EQ.1).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='PILO.GRID')
      ELSEIF ((IH.EQ.2).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='KANLO.GRID')
       ELSEIF ((IH.EQ.2).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='KALO.GRID')
      ELSEIF ((IH.EQ.3).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='PRONLO.GRID')       
      ELSEIF ((IH.EQ.3).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='PROLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='HNLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='HLO.GRID')       
	ELSE
         WRITE(6,93)
 93      FORMAT (2X,' WRONG SET')
         STOP
      END IF
C
       DO 15 M = 1, NX-1 
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M), PARTON(9,N,M)
  90   FORMAT (9(1PE10.3))
c       write(6,'(2i4,8f8.2)') m,n,PARTON(1,N,M), PARTON(2,N,M), 
c     >   PARTON(3,N,M),PARTON(4,N,M) 
  15   CONTINUE
       CLOSE(IIREAD)
C
      FINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUTOTF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**4 * XB0**0.5)
        XDTOTF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.5)
        XSTOTF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**4 * XB0**0.5) 
        XCTOTF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0**0.3) 
        XBTOTF(IX,IQ) = PARTON(5,IQ,IX) / (XB1**7 * XB0**0.3)
        XGF(IX,IQ)    = PARTON(6,IQ,IX) / (XB1**4 * XB0**0.3)
        XUVALF(IX,IQ) = PARTON(7,IQ,IX) / (XB1**4 * XB0**0.5)
        XDVALF(IX,IQ) = PARTON(8,IQ,IX) / (XB1**4 * XB0**0.5)
        XSVALF(IX,IQ) = PARTON(9,IQ,IX) / (XB1**4 * XB0**0.5)
  20  CONTINUE
        XUTOTF(NX,IQ) = 0.D0
        XDTOTF(NX,IQ) = 0.D0
        XSTOTF(NX,IQ) = 0.D0
        XCTOTF(NX,IQ) = 0.D0
        XBTOTF(NX,IQ) = 0.D0
        XGF(NX,IQ)    = 0.D0
        XUVALF(NX,IQ) = 0.D0
        XDVALF(NX,IQ) = 0.D0
        XSVALF(NX,IQ) = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = LOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = LOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = LOG(X)
      XT(2) = LOG(Q2)
      UTOT = FINT(NARG,XT,NA,ARRF,XUTOTF) * (1.D0-X)**4 * X**0.5
      DTOT = FINT(NARG,XT,NA,ARRF,XDTOTF) * (1.D0-X)**4 * X**0.5 
      STOT = FINT(NARG,XT,NA,ARRF,XSTOTF) * (1.D0-X)**4 * X**0.5
      CTOT = FINT(NARG,XT,NA,ARRF,XCTOTF) * (1.D0-X)**7 * X**0.3
      BTOT = FINT(NARG,XT,NA,ARRF,XBTOTF) * (1.D0-X)**7 * X**0.3
      GL   = FINT(NARG,XT,NA,ARRF,XGF)    * (1.D0-X)**4 * X**0.3
      UVAL = FINT(NARG,XT,NA,ARRF,XUVALF) * (1.D0-X)**4 * X**0.5
      DVAL = FINT(NARG,XT,NA,ARRF,XDVALF) * (1.D0-X)**4 * X**0.5 
      SVAL = FINT(NARG,XT,NA,ARRF,XSVALF) * (1.D0-X)**4 * X**0.5
       
       Up  = (UTOT+UVAL)/2.
       UBp = (UTOT-UVAL)/2.
       Dp  = (DTOT+DVAL)/2.
       DBp = (DTOT-DVAL)/2.
       Sp  = (STOT+SVAL)/2.
       SBp = (STOT-SVAL)/2.
       Cp  =  CTOT/2.
       Bp  =  BTOT/2.
              
       IF (IC.EQ.1) THEN
       U  = Up
       UB = UBp
       D  = Dp 
       DB = DBp
       S  = Sp
       SB = SBp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.-1) THEN
       U  = UBp
       UB = Up
       D  = DBp 
       DB = Dp
       S  = SBp
       SB = Sp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.0) THEN
       U  = (UBp+Up)/2.
       UB =  U
       D  = (DBp+Dp)/2. 
       DB =  D
       S  = (SBp+Sp)/2.
       SB =  S
       C  =  Cp
       B  =  Bp 
       ELSE
         WRITE(6,94)
 94      FORMAT (2X,' WRONG CHARGE')
         STOP
       END IF 
c       write(6,'(i2,8f7.3)') ic,up,u,(utot+uval)/2.,
c     >  d, dp, (DTOT+DVAL)/2.
 60   RETURN
       END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(59),TABLE(840)
      DIMENSION D(2),NCOMB(2),IENT(2)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
!---------------------------------------------------------------------
      subroutine FNP_NMC(X,QSQ,rat)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,QSQ,A,B,X2,X3,rat
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      rat = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END

      subroutine getjam(x, q2, u,ubar,d,dbar,s,sbar)
      implicit none
      logical first/.true./
      real*8 x, q2, Q, u,ubar,d,dbar,s,sbar
      real*8 usv(51,13),ubsv(51,13)
      real*8 dsv(51,13),dbsv(51,13)
      real*8 ssv(51,13),sbsv(51,13)
      real*8 xsv(51),qsv(13)
      integer ix,iq,j

      if(first) then
       first=.false.
       open(unit=77, file='jampdf.dat')
       do ix=1,51
        do iq=1,13
         read(77,'(2f6.3,11f8.4)') xsv(ix),qsv(iq),
     >     sbsv(ix,iq),ubsv(ix,iq),dbsv(ix,iq),
     >     dsv(ix,iq), usv(ix,iq),ssv(ix,iq)
        enddo
       enddo
       close(unit=77)
      endif

      ix=0
      iq=0
      do j=1,50
       if(x.ge.xsv(j) .and. x.lt.xsv(j+1)) ix=j
      enddo
      q = sqrt(q2)
      do j=1,12
       if(q.ge.qsv(j) .and. q.lt.qsv(j+1)) iq=j
      enddo

      if(ix.eq.0 .or. iq.eq.0) then
       write(6,'(''error in jampdf'',2f8.3)') x,q
       stop
      endif

      u = usv(ix,iq)
      d = dsv(ix,iq)
      s = ssv(ix,iq)
      ubar = ubsv(ix,iq)
      dbar = dbsv(ix,iq)
      sbar = sbsv(ix,iq)

      return
      end

      SUBROUTINE FFFIT_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! fit for two favored, one unfavored FF (3 params)
! CALCULATE CHISQ FOR MINUIT
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call,iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 fs1, fsb
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
c      mpp = 4 * u * p(1) + 
c     >          d * p(3) + 
c     >      4.*ub * p(2) +
c     >         db * p(4) + 
c     >          s * ffs + 
c     >         sb * ffsb 

c      mpm = 4 * u * p(2) + 
c     >          d * p(4) + 
c     >      4.*ub * p(1) +
c     >         db * p(3) + 
c     >          s * ffsb + 
c     >         sb * ffs 
c      mnp = 4 * d * p(1) + 
c     >          u * p(3) + 
c     >      4.*db * p(2) +
c     >         ub * p(4) + 
c     >          s * ffs + 
c     >         sb * ffsb 

c      mnm = 4 * d * p(2) + 
c     >          u * p(4) + 
c     >      4.*db * p(1) +
c     >         ub * p(3) + 
c     >          s * ffsb + 
c     >         sb * ffs 

      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 

      mpm = 4 * u * p(2) + 
     >          d * p(3) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 

      if(iswap.ne.1) then
      mnp = 4 * d * p(1) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * d * p(2) + 
     >          u * p(3) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      else
      mnp = 4 * d * p(3) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >      4.*db * p(3) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db) + u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mpm,
     >   msv(3),msver(3),mdp,msv(4),msver(4),mdm
       write(6,'(''fchk'',i4,12f7.3)') ncall,p(1),p(2),p(2),p(3),
     >   ffs,ffsb
       write(6,'(''fchk'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFIT4_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fits (two favored, two unfavored)
! depending on "swap", interchange for p, n or not
! CALCULATE CHISQ FOR MINUIT
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 
      mpm = 4 * u * p(3) + 
     >          d * p(4) + 
     >      4.*ub * p(4) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 
c regular (FF don't depend on target)
      if(iswap.eq.0) then
      mnp = 4 * d * p(1) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(3) + 
     >          u * p(4) + 
     >      4.*db * p(4) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
c swap the two favored FF and the two unfavored FF
c as it looks like happens in LUND
      else
      mnp = 4 * d * p(4) + 
     >          u * p(3) + 
     >      4.*db * p(3) +
     >         ub * p(4) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >         db * p(1) +
     >      4.*ub * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db)+ u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFITcsv_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fit assumes twofavored, one
! unfavored FF, and delta u = -delta d
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2,delta_u,delta_d
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 

      mpm = 4 * u * p(2) + 
     >          d * p(3) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 

      delta_d = p(4) * (u + d) * 3. / 8.
      delta_u = -1.0 * delta_d
      if(iswap.ne.1) then       program ptm
c analyze multiplicity data for SIDIS
c do do:

c finish giving wide narrow files different names!

c put smooth interpolation into getrho

c diff. in HMS yptar < 0 and > 0 by 6%.

c multiplicities bigger for phi = pi/2 and 3pi/3 by
c about 10%! (averaged over low pt) when using wide
c acceptance limits.
c Try narrow limits. Could be incorrect HMS offset?

      implicit none

      integer i,j,jp,k,npt,ifit,kk,ipm,kkk,jj,itt,ix,iy,nplt
      integer ikin,ipt,iphi,iz,it,iacc,ikinsv,izsv,ikinp
      integer ngood,ikinx,ipdf,iff,i1,i2,ikinpp
      integer ikinw(8)/29, 13, 9, 5, 11, 1, 7, 3/
      real*8 xpltmin,xpltmax,ypltmin,ypltmax
      real*8 mcsv(3,8,20,10),mcsver(3,8,20,10) 
      real*8 fu1sv(3,8,20),fd1sv(3,8,20),avdelu(3,8),avdeld(3,8)
      real*8 avdeluer(3,8),avdelder(3,8),doverp,doverpp,doverpper
      real*8 mcsvp(3,8,20,12),mcsvper(3,8,20,12),dopsv(3,8,20)
      real*8 mcsv5p(3,8,20,5),mcsv5per(3,8,20,5),uplusd(3,8)
      real*8 delu, deluer, deld, delder,rnp,rsv(20)
      real*8 Mm, Mp, Mmer, Mper, Bm, Bp,rnew,dprop,wp2,fact,du
      real*8 mprat, rd1, rd2, drd, rseas(16,20), rseans(16,20)
      real*8 csvprop(3,20,3,2), csvproper(3,20,3,2),f,Ctq5Pdf
      real*8 BigS, chk1er1, chk1er2,chk1er3,chk1er4,BigT
      real*8 rcsv, rcsver,fu1j,fd1j
      real*8 sumdifffit(56,20,2,2,2)
      real*4 qgev,x4
      logical lastframeused,usewide,okplt
      character*2 ttit(4)/'p+','d+','p-','d-'/
      real*8 x,q2,zmax,w,pt,phicm,z,mmpi2,phi,ratt,wpcut,chk1,chk2
      real*8 sighadp,sighadm,rdm,scsv(56,4,20,4,3),xp,zp,b
      real*8 sighad1, sighad2, sighad3, sighad4,pt2,sig,siger
      real*8 scsver(56,4,20,4,3),xmax,xmin,ymin,ymax,sighad5
      real*8 avcsv(16,3,2),avcsver(16,3,2),pi/3.1415928/
      real*8 ffrat(20,4),ffrater(20,4)
      real*8 ffrat4(20,4),ffrat4er(20,4)
      real*8 yymin,yymax, zforwp225,zforwp23,u1,d1
      real*8 q2n,q2er,xn,xer,zn,zer,rd,rder,rdrho
      real*8 mltnorc, mltnorcer,mlt,mlter,mrho,xx
      real*8 sighad,u,ub,d,db,s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu
      real*8 w2,f1p,f2p,f1d,f2d,a,zz,rplus,rminus,zcut
      real*8 rplusd,rminusd,rat(6),ratrho(6),rater(6)
      real*8 rv(5,20),rver(5,20)
      real*8 avm(56,20,6),avmer(56,20,6),avk(56,20,5)
      real*8 avmod(56,20,6),mmod,mratio
      real*8 avker(56,20),avmrho(56,20,6)
      real*8 avmk(6),avmker(6),lund(16,20,4)
      real*8 avmkr(6),avmkrer(6)
      real*8 avms(56,20,6),rats(6),diff,diffmax
      real*8 srat(56,20,6),srater(56,20,6)
      real*8 sratrho(56,20,6)
      real*8 sq2(10,2),sx(10,2),sz(10,2),smx(10,2)
      real*8 pq2(10,2),px(10,2),pz(10,2),pmx(10,2)
      real*8 pkin(56,2),pw(10,2),rhofactp
      integer rhorat(25,2),usedup(56,20,6)
      integer bnpt,bnpt0,bitv(90000)
      integer srath(40,6),sratrhoh(40,6),rdh(40,2)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000),avmit(5)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      real*8 bchi2x(15),bdfx(15),avit(5),avitrho(5),aviter(5) 
      real*8 ffsv(3,3,20,4),ffsver(3,3,20,4)
      real*8 avffsv(20,4,0:1),avffsver(20,4,0:1)
      real*8 avffsv4(20,4,0:1),avffsv4er(20,4,0:1)
      real*8 avffsvcsv(20,4),avffsvcsver(20,4)
      real*8 avffsvex(20,4),avffsvexer(20,4)
      real*8 avffsv4lund(20,4),avffsv4lunder(20,4)

      integer ncallsfit
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      real*8 rhofact,fbest,chibest,fnorm(36),fnormsv(36)
      integer bkin(90000),ncall/0/
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt,ncallsfit
        integer ikinfit,itfit,npts
        real*8 zfit
        common/sstuff/ zfit,ikinfit,itfit,npts
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag
        integer sn
        real*8 sphi(1000), spt(1000),sy(10000),syer(1000)
        common/sforplot/ sphi,spt,sy,syer,sn

        real*8 ep0sv(28)/
     >    5.24,   3.31,  5.24,  4.49,  5.98,
     >    4.94,   5.98,  6.36,  5.24,  6.36,
     >    4.72,   5.24,  5.71,  4.33,  5.24,
     >    4.87,   3.98,  4.78, 6.590, 5.262,
     >    4.686, 4.183, 3.253, 
     >    2.178, 0.888, 2.320,
     >    1.816, 0.962/
        real*8 the0sv(28)/
     >    13.50, 19.69, 16.31, 16.63, 14.24,
     >    17.26, 15.74, 15.30, 17.21, 13.96,
     >    19.06, 15.38, 17.65, 20.24, 18.51,
     >    19.10, 19.68, 19.69, 11.910, 11.159,
     >    17.08, 14.93,23.005, 
     >    27.25, 36.15, 27.77,
     >    25.895, 49.308 /
      real*8 e0,ep,th,sin2,cos2,nu,q2kin(56),xkin(56),wkin(56)
      real*8 e0sv(56),am/0.938/
      real*8 qu,qd,qs,sum_sqp,sum_sqn,w2chk
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)

      real*8 coef(100)  ! final values of coeficients
      integer nparam    ! number of parameters
      integer arglis(10) ! to pass info to Minuit
      real*8 std(100),chi2,grad(100),p(100)
      integer  nerror_migrad,ierflg,npar
      real*8 zero/0./  ! to pass to Minuit a 0.
      character*10 pname(100)
      external bfit_fcn,sfit_fcn,mfit_fcn,mfit5_fcn,fffit_fcn 
      external fffit4_fcn,fffitcsv_fcn,fffitex_fcn 
      real*8 futil ! auxially function
      external futil
      real*8 mmpi2cut,q2cut,wcut,ptcut,ptcutlo
c inclusive d/p from bcm.out using scalers
      real*8 rrsv(8)/0.86, 0.86, 0.81, 0.86,
     >  0.80, 0.81, 0.78, 0.77/
c shuo data
      character*200 string
       real*8 avq2,avq2er,avx,avxer,zset,avz,avzer,yrd,q2set,xset,zbin,
     >  yrder,yrdrho,ry,ryrho,ryrho2,rynoex,rynodel,ryerr,xbin
       real*8 sjryav(32,20),sjryaver(32,20),v1,v2,v3,v4,v5,v6
       real*8 sjryrhoav(32,20)
       integer group,ncallff,iswap,ioff

       real*8 sum1,sum2,phimin,phimax
       real*8 mvlund(4),chi2ff
       real*8 mv(4),mver(4),muf,mubf,mdf,mdbf,msf,msbf,mvf(4)
       common/ffstuff/ mv,mver,muf,mdf,mubf,mdbf,msf,
     >  msbf,fs1,fsb,mvf,ncallff,iswap,chi2ff

c this is for DSS
      integer IHdss,ICdss,IOdss, fini 
      real*8 fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1
      COMMON / FRAGINI / FINI
      fini=0

c note: low cut values!
      mmpi2cut = 2.0
      wcut=1.8
      q2cut=1.8
      zcut=0.96
c note: very low cut values!
      mmpi2cut = 1.1
c more reasonable cut
c      mmpi2cut = 2.56
      wcut=2.0
      q2cut=1.0
      zcut=0.95
c  can chang from default 
c      ptcut = 3./16.
c range of pt to use for averages
c large value!
      ptcut = 4./16.
      ptcutlo = 2./16.
c wide or narrow acceptance?
      usewide=.true.
c      usewide = .false.

! ASSIGN UNITS (NORMALLY 5,6,7 FOR INTERACTIVE)
      CALL MNINIT(61,62,63)

      if(usewide) then
       open(unit=6,file='ptmw.out')
       open(unit=16,file='ptmsfitw.out')
       open(unit=17,file='ptmnfitw.out')
      else
       open(unit=6,file='ptm.out')
       open(unit=16,file='ptmsfit.out')
       open(unit=17,file='ptmnfit.out')
      endif

c read in Lund multiplicities
      open(unit=9,file='ptm.lund')
      do ikin=1,8
       read(9,'(a)') string
       read(9,'(a)') string
       read(9,'(a)') string
       read(9,'(a)') string
       do iz=1,19
        read(9,*) i1,i2,lund(ikin,iz,1), ! p pi+
     >    lund(ikin,iz,2),! p pi-
     >    v1,lund(ikin,iz,3), ! n pi+
     >    lund(ikin,iz,4) ! n pi-
       enddo
      enddo
      close(unit=9)

c check of CSV formula
      do i=1,3
       rnew = 0.25 * i
       do j=1,6
        du = 0.1 * (j-1)
        mprat = (4.*(1-du)*rnew +  (1 + du))/
     >          (4.*(1-du) + (1 +du)*rnew)
        Dprop = (1. - rnew) / (1. + rnew)
        RD = (4.*mprat - 1.) / (1. - mprat)
c this gives delu / (u + d + ub + db), if s=sbar=0
        chk1 = (4*rnew + 1 - mprat * (4 + rnew)) /
     >         (4*rnew - 1 - mprat * (4 - rnew))
c this is -4/3 (delta u - delta d) / (u + d)
c or -8/3 delta u / (u+d)
        chk2 = 2.5  - Dprop * (2.5 + RD)
        write(6,'(''chkcsv'',6f7.2)') 
     >   rnew, du, chk1, 
     >    -3./8.*chk2
        enddo
       enddo
      do i=1,9
       z = 0.1*i
       q2 = 2.
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
       write(6,'(''fdss'',8f7.4)') z,q2, fu1,fd1/fu1,
     >  fub/fu1, fdb/fu1, fs1/fu1, fsb/fu1
c ub is same as d1
c db is a bit bigger than u1
c s1 and sb are about same as d1
      enddo

      do ikin=1,56
        e0sv(ikin)=10.600
        if(i kin.ge.25 .and. ikin.le.32) e0sv(ikin)=10.214
        if(ikin.ge.47 .and. ikin.le.50) e0sv(ikin)=6.190
        if(ikin.ge.51 .and. ikin.le.56) e0sv(ikin)=8.209
      enddo
c get kinematics by setting
      do ikin=1,56
       e0 = e0sv(ikin)
       ep = ep0sv( (ikin+1) / 2)
       th = 3.14159 / 180. * the0sv( (ikin+1) / 2)
       if( (ikin/2)*2 .eq. ikin) then
        th = th + 0.010
       else
        th = th - 0.010
       endif
       sin2 = sin(th/2)**2
       cos2 = 1. - sin2
       q2kin(ikin) = 4. * e0 * ep * sin2
       nu = e0 - ep
       xkin(ikin) = q2kin(ikin) / 2. / am / nu
       x = xkin(ikin)
       q2 = q2kin(ikin)
       w = sqrt(am**2 + 2.*am*nu -q2)
       wkin(ikin) = w
      enddo

c read in Shuo's results for yield ratios from d (pi-/pi+) in grid of
c 20 x by 20 z bins for each kin. setting
c      open(unit=9, file='sjcsv_simadd.txt')
      open(unit=9, file='csv_datasub.csv')
      read(9,*) string
      do i=1,669
       read(9,*) q2set,avq2,avq2er,xset,xbin,avx,avxer,zset,zbin,avz,
     >  avzer,group,yrd,yrder,yrdrho,v1,v2,
     >  ry,ryrho,ryrho2,rynoex,rynodel,
     >  ryerr
       ikin=0.
       do j=1,32,2
        if(abs(q2set - (q2kin(j)+q2kin(j+1))/2.).lt.0.3 .and.
     >     abs(xset - (xkin(j)+xkin(j+1))/2.).lt.0.03) ikin=j
       enddo
       if(ikin.eq.0) write(6,'(''error,ikin=0'')')
       iz = int(20.*avz)+1
       if(ryerr.lt.0.2) then
        sjryav(ikin,iz) = sjryav(ikin,iz) + ry/ryerr**2 
        sjryrhoav(ikin,iz) = sjryrhoav(ikin,iz) + 
     >    ryrho/ryerr**2 
        sjryaver(ikin,iz) = sjryaver(ikin,iz) + 1./ryerr**2 
       endif
       write(6,'(2i4,i3,10f7.2)') group,ikin,iz,xset,q2set,xbin,zbin,
     >  avz,ry,ryrho,ryerr
      enddo
      close(unit=9)
      do ikin=1,32
       do iz=1,20
        if(sjryaver(ikin,iz).ne.0.) then
         sjryav(ikin,iz) = sjryav(ikin,iz)/sjryaver(ikin,iz)
         sjryrhoav(ikin,iz) = sjryrhoav(ikin,iz) /
     >     sjryaver(ikin,iz)
         sjryaver(ikin,iz) = sqrt(1./sjryaver(ikin,iz))
        endif
       enddo
      enddo

c wide limits
      if(usewide) then
       open(unit=277,file='ptm.inpw')
c narrow limits
      else
       open(unit=277,file='ptm.inpn')
      endif
      bnpt = 0
       do i=1,1000000
       read(277,'(i2,3i3,2i2,6f6.3,6f10.4,e12.4)',end=11) 
     >  ikin,ipt,iphi,iz,it,iacc,x,q2,
     >  z,pt,phicm,mmpi2,
     >  mlt,mlter,mltnorc,mltnorcer,mrho,
     >  mratio,mmod
       w = sqrt(0.938**2 + q2 * (1./x -1))
       if(mmpi2.gt.mmpi2cut .and. it.ne.3 .and. it.ne.6
     >  .and. ikin.le.56 .and. z.lt.zcut
c for test
c     >  .and. ikin.gt.2
     >  .and. mlter.gt.0.
     >  .and. q2.gt.q2cut .and. w.gt.wcut) then
c special corr. for f1f2in21 being 5% begger than
c f1f2in09 for deuteron.
c took out because now using f1f2in09 in SIMC
c        if(it.eq.2 .or. it.eq.5) then
c         mlt = mlt * 1.05
c         mlter = mlter * 1.05
c        endif
c for test
c        mlt = mltnorc
c        mlter = mltnorcer

c for test
c        if(it.eq.2) then
c         mlt = mlt * 1.065
c         mlter = mlter * 1.065
c        endif
             bnpt = bnpt + 1
             bkin(bnpt) = ikin
             bitv(bnpt) = it
             bq2v(bnpt) = q2
             bxv(bnpt) = x
             bmmpi2(bnpt) = mmpi2
c epsilon not used
c             berv(bnpt) = erkin(ikin)
             bzv(bnpt) = z
             bptv(bnpt) = pt
             phi = 2. * 3.1415928 / 15. * (iphi-0.5)
             bphiv(bnpt) = phi
             byv(bnpt) = mlt
             byerv(bnpt) = mlter
             brho(bnpt) = mrho
        call getrho(xkin(ikin),q2kin(ikin),z,pt,phi,rplus,rminus,
     >    rplusd, rminusd)
        if(it.eq.1) brhop(bnpt) = mlt * rplus
        if(it.eq.4) brhop(bnpt) = mlt * rminus
        if(it.eq.2) brhop(bnpt) = mlt * rplusd
        if(it.eq.5) brhop(bnpt) = mlt * rminusd

        if(it.eq.1 .and. rplus .gt.0.02) then
         kkk = min(25,max(1,int(mrho/mlt/rplus*4 .)+1))
         rhorat(kkk,1) = rhorat(kkk,1)+1
        endif
        if(it.eq.4 .and. rminus.gt.0.02) then
         kkk = min(25,max(1,int(mrho/mlt/rminus*4.)+1))
         rhorat(kkk,2) = rhorat(kkk,2)+1
        endif

        ipdf = 1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c        write(6,'(''tst'',i6,4f10.3)') i,mlt,mlter,
c     >   sighad,mlt/sighad
        buv(bnpt) = u
        bubv(bnpt) = ub
        bdv(bnpt) = d
        bdbv(bnpt) = db
        bsv(bnpt) = s
        bsbv(bnpt) = sb
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        bffu(bnpt) = fu1 / z
        bffd(bnpt) = fd1 / z
        bffs(bnpt) = fs1 / z
        bffub(bnpt) = fub / z
        bffdb(bnpt) = fdb / z
        bffsb(bnpt) = fsb / z
        if(mrho/mlt.gt.0.2) write(6,'(''rho'',
     >   i2,9f6.2,f7.3)') 
     >   it,x,q2,z,pt,phicm,mmpi2,mlt, mlter,mlt/sighad,
     >   mrho/mlt
c ratio of pi-/pi+ deuteron, etc. 
c these are averages of FF (dsigdz, based on M0), not sighad
c they should match fDSS
        if(pt.lt.ptcut .and. pt.gt.ptcutlo) then
c used in SIMC. see v2.pdf
c         b = 1./ (0.120 * z**2 + 0.200)
c rough fit to pt-SIDIS results. see v1.pdf
         zmax = min(0.6, z)
         b = 1./ sqrt(0.16**2 + zmax**2 * 0.30**2)
         if(it.eq.1.or.it.eq.2) then
c          b = 1./ sqrt(0.15**2 + zmax**2 * 0.40**2)
         endif
        fact = 1./ (b * exp(-b * pt**2)/ 2. / pi)
c for test of systematics see v3.pdf
c         fact = 1.
         write(6,'(''tst'',4i3,5f8.3)') ikin,iz,ipt,iphi,
     >    mlt,mlter,
     >    sighad,mlt/sighad,mmod/sighad
        avm(ikin,iz,it) = avm(ikin,iz,it) +
     >    fact * mlt / mlter**2 / fact**2
        avms(ikin,iz,it) = avms(ikin,iz,it) +
     >    fact * sighad / mlter**2 / fact**2
        avmrho(ikin,iz,it) = avmrho(ikin,iz,it) +
     >    fact * mrho/ mlter**2 / fact**2
        avmer(ikin,iz,it) = avmer(ikin,iz,it) +
     >    1. / mlter**2 / fact**2
        avmod(ikin,iz,it) = avmod(ikin,iz,it) +
     >    fact * mmod / mlter**2 / fact**2
        avk(ikin,iz,1) = avk(ikin,iz,1) + x/ mlter**2 / fact**2 
        avk(ikin,iz,2) = avk(ikin,iz,2) +q2/ mlter**2 / fact**2 
        avk(ikin,iz,3) = avk(ikin,iz,3) + w/ mlter**2 / fact**2 
        avk(ikin,iz,4) = avk(ikin,iz,4) + z/ mlter**2 / fact**2 
        avk(ikin,iz,5) = avk(ikin,iz,5) + mmpi2/mlter**2 / fact**2 
        avker(ikin,iz) = avker(ikin,iz) + 1./mlter**2  / fact**2
        endif ! pt cut
       endif
      enddo
 11   continue

c results for pi-/pi+ deuteron vs z
      if(usewide) then
       open(unit=14,file='ptmratzw.txt')
       open(unit=15,file='ptmratzrhow.txt')
      else
       open(unit=14,file='ptmratz.txt')
       open(unit=15,file='ptmratzrho.txt')
      endif
      do ikin=1,56
       do iz=1,20
c super ratios of data to model 
        do it=1,6
         srat(ikin,iz,it)=0.
         srater(ikin,iz,it)=0.
         usedup(ikin,iz,it)=0
        enddo
        do it=1,6
          if(avmer(ikin,iz,it).ne.0.) then
           avm(ikin,iz,it) = avm(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmod(ikin,iz,it) = avmod(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avms(ikin,iz,it) = avms(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmrho(ikin,iz,it) = avmrho(ikin,iz,it) / 
     >      avmer(ikin,iz,it)
           avmer(ikin,iz,it) = 1./
     >       sqrt(avmer(ikin,iz,it))
          endif
        enddo
        if(avmer(ikin,iz,2).ne.0. .and.
     >      avker(ikin,iz).ne.0.) then
         do it=1,6
          rat(it)=0.
          rats(it)=0.
          ratrho(it)=0.
          rater(it)=0.
          if(avmer(ikin,iz,it).ne.0 .and.
     >       avmer(ikin,iz,2).ne.0.) then
           rat(it) = avm(ikin,iz,it) / 
     >               avm(ikin,iz,2)
           rats(it) = avms(ikin,iz,it) / 
     >                avms(ikin,iz,2)
           ratrho(it) = avmrho(ikin,iz,it) / 
     >               avmrho(ikin,iz,2)
           rater(it) = rat(it) * sqrt(
     >      (avmer(ikin,iz,it)/avm(ikin,iz,it))**2 +
     >      (avmer(ikin,iz,2)/avm(ikin,iz,2))**2)
           sratrho(ikin,iz,it) = ratrho(it)/rats(it)
           srat(ikin,iz,it) = rat(it)/rats(it)
           srater(ikin,iz,it) = rater(it)/rats(it)
           diff = (srat(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
           diff = max(-7.,min(6.99,diff))
           k = int((diff + 7.)/14.*40.)+1
           if(k.gt.40) write(6,'(i5,2f8.3)') k,
     >      srat(ikin,iz,it),srater(ikin,iz,it)
           srath(k,it) = srath(k,it)+1
           diff = (sratrho(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
           diff = max(-7.,min(6.99,diff))
           k = int((diff + 7.)/14.*40.)+1
           sratrhoh(k,it) = sratrhoh(k,it)+1

           if(it.eq.5) then
            diff = (sratrho(ikin,iz,it) - 1.0) /
     >       srater(ikin,iz,it)
c            write(6,
c     >      '(''err diff'',4e10.2)') sratrho(ikin,iz,it),
c     >       srat(ikin,iz,it),srater(ikin,iz,it),diff
            x = avk(ikin,iz,1)/avker(ikin,iz)
            q2 = avk(ikin,iz,2)/avker(ikin,iz)
            w = avk(ikin,iz,3)/avker(ikin,iz)
            z = avk(ikin,iz,4)/avker(ikin,iz)
            mmpi2 = avk(ikin,iz,5)/avker(ikin,iz)
       k=max(1,min(10,int((q2-3.)/3.5*10)+1))
       pq2(k,1) = pq2(k,1) + diff
       pq2(k,2) = pq2(k,2) + 1.
       k=max(1,min(10,int((x-0.25)/0.35*10)+1))
       px(k,1) = px(k,1) + diff
       px(k,2) = px(k,2) + 1.
       k=max(1,min(10,int((z-0.3)/0.5*10)+1))
       pz(k,1) = pz(k,1) + diff
       pz(k,2) = pz(k,2) + 1.   
       k=max(1,min(10,int((w-1.8)/1.5*10)+1))
       pw(k,1) = pw(k,1) + diff
       pw(k,2) = pw(k,2) + 1.   
       k=max(1,min(10,int((mmpi2-2.)/5*10)+1))
       pmx(k,1) = pmx(k,1) + diff
       pmx(k,2) = pmx(k,2) + 1.   
           endif
          endif
         enddo
         if(rater(5).lt.0.2) then
          write(14,'(2i3,5f6.3,9f7.3)') ikin,iz,
     >     (avk(ikin,iz,jj)/avker(ikin,iz),jj=1,5),
     >     rat(5),rater(5),rat(1),rater(1),
     >     rat(4),rater(4)
          write(15,'(2i3,5f7.3,6f7.3)') ikin,iz,
     >     (avk(ikin,iz,jj)/avker(ikin,iz),jj=1,5),
     >     ratrho(5),rater(5),ratrho(1),rater(1),
     >     ratrho(4),rater(4)
         endif
        endif
       enddo
      enddo

c write super ratios in order of worst first
      if(usewide) then
       open(unit=7,file='ptmsratw.txt')
      else
       open(unit=7,file='ptmsrat.txt')
      endif
      do it=1,6
       if(it.eq.1.or.it.eq.4.or.it.eq.5) then
        do i=1,100000
         diffmax=0.
         ikinsv=0
         izsv=0
         do ikin=1,32
          do iz=1,20
           if(usedup(ikin,iz,it).eq.0.and.
     >      srater(ikin,iz,it).ne.0.) then
            diff = abs(srat(ikin,iz,it) - 1.) /
     >       srater(ikin,iz,it)
            if(diff.gt.diffmax) then
             diffmax = diff
             ikinsv = ikin
             izsv = iz
            endif
           endif
          enddo
         enddo ! ikin
         if(diffmax.eq.0.) goto 12
         usedup(ikinsv,izsv,it)=1
         write(7,'(i1,2i3,5f6.3,2f7.3)')
     >    it,ikinsv,izsv,
     >    (avk(ikinsv,izsv,jj)/avker(ikinsv,izsv),jj=1,5),
     >    srat(ikinsv,izsv,it),
     >    srater(ikinsv,izsv,it)
        enddo
 12     continue
       endif
      enddo

c table of average multiplicites (<ptcut) versus z
      close(unit=22)
      close(unit=23)
      if(usewide) then
       open(unit=22,file='ptmavmw.txt')
       open(unit=23,file='ptmavmkltw.txt')
      else
       open(unit=22,file='ptmavmn.txt')
       open(unit=23,file='ptmavmkltn.txt')
      endif
      do ikin=1,56,2
       x = (xkin(ikin) + xkin(ikin+1))/2.
       q2 = (q2kin(ikin) + q2kin(ikin+1))/2.
       w2 = am**2 + q2*(1/x-1)
       do iz=4,19
        do it=1,5
         avit(it)=0.
         avitrho(it)=0.
         aviter(it)=0.
         If(avmer(ikin ,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).gt.0.) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v3 = avmod(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avmod(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avit(it) = v1 / v2
          avmit(it) = v3 / v2
          aviter(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >         (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >          avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avitrho(it) = v1 / v2
         endif
        enddo ! it
        z = .05 * (iz-0.5)
        write(22,'(i2,3f6.3,12f8.4)') (ikin+1)/2,x,q2,z,
     >  avit(1),avitrho(1),aviter(1),
     >  avit(2),avitrho(2),aviter(2),
     >  avit(4),avitrho(4),aviter(4),
     >  avit(5),avitrho(5),aviter(5) 
c        write(22,'(20x,4(f8.4,16x))')
c     >   avmit(1),avmit(2),avmit(4),avmit(5)
        if(ikin.gt.0) then
         write(23,'(i2,3f6.3,8f8.4)') 
     >    (ikin+1)/2,x,q2,z,
     >    avit(1)/avmit(1),aviter(1)/avmit(1),
     >    avit(2)/avmit(2),aviter(2)/avmit(2),
     >    avit(4)/avmit(4),aviter(4)/avmit(4),
     >    avit(5)/avmit(5),aviter(5)/avmit(5)
        endif
       enddo ! iz
      enddo ! ikin

c plot multiplicity ratios over deuteron pi+ versus z  and compare
c to SHuo if available.
c also, get two ratio of FF and plot those two
      close(unit=21)
      close(unit=22)
      close(unit=23)
      open(unit=21,file='ptmavm.top')
      write(21,'(''set device postscript'')')
      open(unit=22,file='ptmavmrat.top')
      write(22,'(''set device postscript'')')
      open(unit=23,file='ptmavmratff.top')
      write(23,'(''set device postscript'')')
      open(unit=25,file='ptmavmratff4.top')
      write(25,'(''set device postscript'')')
      open(unit=31,file='ptmavmratff4no.top')
      write(31,'(''set device postscript'')')
      open(unit=24,file='ptmavmratffrm.top')
      write(24,'(''set device postscript'')')
      open(unit=26,file='ptmavmratffr.top')
      write(26,'(''set device postscript'')')
      open(unit=27,file='ptmavmratcsv.top')
      write(27,'(''set device postscript'')')
      open(unit=28,file='ptmavmratex.top')
      write(28,'(''set device postscript'')')
      close(unit=29)
      open(unit=29,file='ptmavmratcsv2.top')
      write(29,'(''set device postscript'')')
      close(unit=30)
      open(unit=30,file='ptmdiffrat.top')
      write(30,'(''set device postscript'')')
      ix=0
      iy=3
c only do cases with proton
      do ikinpp=1,8
       ikin = 2*ikinpp - 1
       if(ikinpp.eq.8) ikin=29
c change order to increasing W
       ikin = ikinw(ikinpp)
c      do ikin=1,32,2
       ix=ix+1
       if(ix.gt.4) then
        ix=1
        iy = iy-1
       endif
       xpltmin = 1.1 + 2.95*(ix-1)
       xpltmax = xpltmin + 2.95
       ypltmax = 9.5 - 4.0*(3-iy)
       ypltmin = ypltmax - 4.0
       x = (xkin(ikin) + xkin(ikin+1))/2.
       q2 = (q2kin(ikin) + q2kin(ikin+1))/2.
       w2 = am**2 + q2*(1/x-1)
c       write(21,2188) ix,iy,x,q2,sqrt(w2)
       do kk=1,12
        if(ix.eq.1) write(20+kk,'(''set labels left on'')')
        if(ix.ne.1) write(20+kk,'(''set labels left off'')')
        if(iy.eq.2) write(20+kk,'(''set labels bottom on'')')
        if(iy.ne.2) write(20+kk,'(''set labels bottom off'')')
        write(20+kk,2088) xpltmin,xpltmax,ypltmin,ypltmax,
     >   xpltmin+0.25, ypltmax-0.35,x,q2,sqrt(w2)
 2088   format(
     >   'set intensity 4'/
     >   'set color white'/
     >   'set window x ',2f8.3,' y ', 2f8.3/
     >   'title ',2f8.3,' size 1.6 ',1h',
     >   'x=',f4.2,'  Q223=',f3.1,' W=',f3.1,1h'/
     >   'case ',1h','         X X',1H'/
     >   'set symbol 9O size 1.0 ; set bar size 0.'/
     >   'set ticks size 0.04 ; set order x y dy sym')
       enddo ! kk

       if(ix.eq. 1.and.iy.eq.2) then
        write(21,8801)
 8801   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','zM(z)',1h')
        write(23,8803)
 8803   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z)',1h')
        write(24,8804)
 8804   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
c     >          1h','D(z) / D0DSS3(z)',1h'/
     >          1h','D(z) / D0JAM3(z)',1h'/
     >  'CASE ',1H','        X   X   ',1H')
        write(25,8805)
        write(31,8805)
 8805   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','zD(z)',1h')
        write(26,8806)
 8806   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >           1h','D(z) / D0u12+3(z)',1h'/
     >  'case  ',1h','        X XX X   ',1h')
        write(27,8807)
 8807   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z)',1h')
        write(28,8808)
 8808   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','D(z) or excl factor',1h')
        write(29,8809)
 8809   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','CSV(x)',1h')
        write(30,8810)
 8810   format('title 6.7 0.9 size 2.2 ',1h','z',1h'/
     >   'title 0.12 5.2 size 2.2 angle=90 ',
     >   1h','Ratio',1h')
        endif

       write(21,2188)
 2188  format(
     > 'set limits x 0.25 0.75 y 0.0 0.49'/
     > 'plot axes')
       write(22,2288)
 2288  format(
     > 'set limits x 0.25 0.75 y 0.35 1.4'/
     > 'plot axes')
       write(23,2372)
 2372  format(
     > 'set limits x 0.25 0.75 y 0.0 0.57'/
     > 'plot axes')
       write(24,2472)
 2472  format(
     > 'set limits x 0.25 0.75 y 0.0 1.9'/
     > ' 0. 1. ; 1. 1. ; set pattern .03 .03 .03 .03 ; join pattern'/
     > 'plot axes')
       write(25,2572)
 2572  format(
     > 'set limits x 0.25 0.75 y  0.  0.57'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(31,3172)
 3172  format(
     > 'set limits x 0.25 0.75 y  -0.2  0.57'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(26,2562)
 2562  format(
     > 'set limits x 0.25 0.75 y -0.5 1.5'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(27,2762)
 2762  format(
     > 'set limits x 0.25 0.75 y 0. 0.59'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(28,2862)
 2862  format(
     > 'set limits x 0.25 0.75 y 0 1.5'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(29,2962)
 2962  format(
     > 'set limits x 0.25 0.75 y -0.6 0.59'/
     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')
       write(30,3062)
 3062  format(
     > 'set limits x 0.25 0.75 y 0.3 1.25'/
c     > '0. 0. ; 1. 0. ; set pattern .02 .02 .02 .02 ; join pattern'/
c     > '0. 1. ; 1. 1. ; set pattern .02 .02 .02 .02 ; join pattern'/
     > 'plot axes')

c my multiplicity results
       do it=1,5
       if(it.eq.1.or.it.eq.4.or.it.eq.5.or.it.eq.2) then
       if(it.eq.1) write(21,'(''set color blue'')')
       if(it.eq.2) write(21,'(''set color red'')')
       if(it.eq.4) write(21,'(''set color green'')')
       if(it.eq.5) write(21,'(''set color white'')')
       nplt = 0
       do iz=1,20
        z = .05 * (iz-0.5)
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.3 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.3) then
         v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >        avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         write(21,'(3f10.4,'' 9O'')')  
     >    0.05*(iz-0.5), z*v1,z*v2
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) then
            write(21,'(''set symbol size 1.0 ; plot'')')
            write(21,'(''set symbol size 0.8 ; plot'')')
            write(21,'(''set symbol size 0.6 ; plot'')')
            write(21,'(''set symbol size 0.4 ; plot'')')
            write(21,'(''set symbol size 0.2 ; plot'')')
       endif
! with rho subtraction
       nplt = 0
       do iz=1,20
        z = .05 * (iz-0.5)
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.3 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.3) then
         v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         write(21,'(3f10.4,'' 9O'')')  
     >    0.05*(iz-0.5)-0.01, z*v1,z*v2
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(21,'(''set symbol size 1.0 ; plot'')')
       write(21,'(''set intensity 2'')')

c plot various simc models
c now only plot dss (my fit too big) and one pdf
c (they give almost identical results)
c now plot dss with/without swap, jam (2, 3, 4)
       do iff=2,4
        do ipdf=2,2
         okplt = .true.
         if(iff.eq.3 .and. (it.eq.1 .or. it.eq.4)) okplt=.false.
         do iz=6,18
          z = 0.05*iz
          pt=0.
          mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >      sighadm,u,ub,d,db,u1,d1,
     >      s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          if(okplt) write(21,'(f8.3,f10.4,2f8.3,2i2)') z,z*dsigdz
         enddo
         if(okplt) then
          if(iff.eq.2.and.ipdf.eq.2) 
     >     write(21,'(''set pattern .05 .05 .05 .05'')')
          if(iff.eq.3.and.ipdf.eq.2) 
     >     write(21,'(''set pattern .015 .015 .015 .015'')')
          if(iff.le.3) write(21,'(''join pattern'')')
          if(iff.eq.4) write(21,'(''join'')')
         endif
        enddo ! pdf
       enddo ! iff
       endif ! it
       enddo ! it

! my ratio to D+ results
       do it=1,5
       do iz=1,20
        rver(it,iz)=0.
       enddo
       if(it.eq.1.or.it.eq.4.or.it.eq.5) then
       if(it.eq.1) write(22,'(''set color blue'')')
       if(it.eq.4) write(22,'(''set color green'')')
       if(it.eq.5) write(22,'(''set color white'')')
       nplt = 0
       do iz=1,20
        If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz, 2).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin,iz, 2).lt.0.2.and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz, 2).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz, 2).lt.0.2) then 
         v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >        avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v3 = avm(ikin,  iz, 2)/avmer(ikin,  iz, 2)**2 +
     >        avm(ikin+1,iz, 2)/avmer(ikin+1,iz, 2)**2 
         v4 =              1./avmer(ikin,  iz, 2)**2 +
     >                     1./avmer(ikin+1,iz, 2)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         v3 = v3 / v4
         v4 = sqrt(1./v4)
         v5 = v1 / v3
         v6 = sqrt( (v2/v3)**2 + (v5/v3*v4)**2)
         rv(it,iz) = v1
         rver(it,iz) = v2
         rv(2,iz) = v3
         rver(2,iz) = v4
         write(22,'(3f10.4,'' 9O'')')  0.05*(iz-0.5), v5, v6
c do again with rho subtraction
         v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
         v2 =                1./avmer(ikin,iz,it)**2 +
     >                       1./avmer(ikin+1,iz,it)**2 
         v3 = (avm(ikin,  iz, 2)-avmrho(ikin,  iz,2))/
     >          avmer(ikin,  iz, 2)**2 +
     >        (avm(ikin+1,iz, 2)-avmrho(ikin+1,iz,2))/
     >          avmer(ikin+1,iz, 2)**2 
         v4 =              1./avmer(ikin,  iz, 2)**2 +
     >                     1./avmer(ikin+1,iz, 2)**2 
         v1 = v1 / v2
         v2 = sqrt(1./v2)
         v3 = v3 / v4
         v4 = sqrt(1./v4)
         v5 = v1 / v3
         v6 = sqrt( (v2/v3)**2 + (v5/v3*v4)**2)
         write(22,'(3f10.4,'' 8O'')')  0.05*(iz-0.5)-0.01, v5, v6
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(22,'(''plot'')')
c plot model
       nplt = 0
       do iz=4,18
        if(avms(ikin,iz,it).gt.0. .and.
     >     avms(ikin,iz,2).gt.0.) then
         write(22,'(3f10.4)')  0.05*(iz-0.5),
     >    avms(ikin,iz,it)/avms(ikin,iz,2)
         nplt = nplt + 1
        endif
       enddo
       if(nplt.gt.0) write(22,'(''join'')')
       write(22,'(''set order x y dum dum'')')
c plot various simc models
       do iff=1,2
        do ipdf=1,2
         do iz=4,18
          z = 0.05*iz
          pt=0.
          mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >     sighadp,u,ub,d,db,u1,d1,
     >     s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >      sighadm,u,ub,d,db,u1,d1,
     >      s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
          write(22,'(f8.3,3f10.4,2f8.3,2i2)') z,sighadm/sighadp,
     >      sighadm,sighadp,u,d,iff,ipdf
         enddo
         if(iff.eq.1.and.ipdf.eq.1) 
     >    write(22,'(''set pattern .01 .01 .01 .01'')')
         if(iff.eq.1.and.ipdf.eq.2) 
     >    write(22,'(''set pattern .01 .05 .01 .05'')')
         if(iff.eq.2.and.ipdf.eq.1) 
     >    write(22,'(''set pattern .02 .02 .02 .02'')')
         if(iff.eq.2.and.ipdf.eq.2) 
     >    write(22,'(''set pattern .05 .05 .05 .05'')')
         write(22,'(''join pattern'')')
        enddo
       enddo
c plot lund predictions
c values are from low x, but plot all of them
       if(xkin(ikin).lt.0.60) then
        write(22,'(''set order x y '')')
        call FNP_NMC(X,Q2,rnp)
        ikinp = max(1,min(8,int((x-0.275)/0.05)+1))
c changed to use lowest x bin
        ikinp = 1
        do iz=3,18
         v2 = (lund(ikinp,iz,1)  + rnp * 
     >         lund(ikinp,iz,3)) / (1 + rnp)
         if(it.eq.1) v1 = lund(ikinp,iz,1)
         if(it.eq.4) v1 = lund(ikinp,iz,2)
         if(it.eq.5) v1  =(lund(ikinp,iz,2)  + rnp * 
     >         lund(ikinp,iz,4)) / (1 + rnp)
         mvlund(1) = lund(ikinp,iz,1)
         mvlund(2) = v2
         mvlund(3) = lund(ikinp,iz,2)
         mvlund(4) = (lund(ikinp,iz,2)  + rnp * 
     >         lund(ikinp,iz,4)) / (1 + rnp)
         write(22,'(2f8.4)') 0.05*(iz-0.5),v1/v2
        enddo
        write(22,'(''set pattern .09 .09 .09 .09'')')
        write(22,'(''join pattern'')')
c        write(22,'(''join ; set color white'')')
        write(22,'(''set order x y dy sym'')')
       endif
       write(6,'(''xx3, it, ikin'',3i3)') it,ikin
      endif ! it
      write(6,'(''xx3, it, ikin'',3i3)') it,ikin
      enddo ! it

      write(6,'(''xx3, it, ikin'',3i3)') it,ikin
c fit 3 D(z) and then 4 D(z) for this kin
c either swapping the two favoreds or not
      do iz=3,19
       do j=1,4
        do iswap=0,1
         avffsver(iz,j,iswap)= 0.
         avffsv4er(iz,j,iswap)= 0.
        enddo
        ffrater(iz,j)=0.
        ffrat4er(iz,j)=0.
        avffsvcsver(iz,j)= 0.
       enddo
       if(rver(1,iz).gt.0. .and.
     >    rver(2,iz).gt.0. .and.
     >    rver(4,iz).gt.0. .and.
     >    rver(5,iz).gt.0.) then
        z = 0.05 * (iz-0.5)
! fit to data with NO rho subratcion
        mv(1) = rv(1,iz) * z
        mv(2) = rv(2,iz) * z
        mv(3) = rv(4,iz) * z
        mv(4) = rv(5,iz) * z
        mver(1) = rver(1,iz) * z
        mver(2) = rver(2,iz) * z
        mver(3) = rver(4,iz) * z
        mver(4) = rver(5,iz) * z
        it = 1
        pt = 0.
        phicm = 0.
        ipdf = 1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,muf,mubf,mdf,mdbf,u1,d1,
     >   msf,msbf,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        write(6,'(''rfffit'',3i3,8f7.3)') ikin,iz,ifit,
     >   (mv(it),mver(it),it=1,4)
c fit with two fav. one unfavored
        do iswap=0,1
         nparam=3
         npar = 3
         p(1) = fu1 ! fav u->pi+
         p(2) = fub ! unfav u->pi+
         p(3) = fdb
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
         arglis(1)=0
         ncallff=0
         call mnexcm(fffit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          avffsv  (iz,j,iswap)= coef(j)
          avffsver(iz,j,iswap)= std(j)
c          write(6,'(''rfffit'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
c     >     j,p(j),coef(j),std(j)
         enddo ! j
         write(6,'(''3param'',2i3,i2,5f7.3)') ikin,iz,
     >    iswap,chi2ff,
     >    (mv(1) - mv(3))/(mv(2) - mv(4)),
     >    (mvf(1) - mvf(3))/(mvf(2) - mvf(4)),
     >    (mv(1) + mv(3))/(mv(2) + mv(4)),
     >    (mvf(1) + mvf(3))/(mvf(2) + mvf(4))
         do it=1,4
          ffrat(iz,it) = mv(it)/mvf(it)
          ffrater(iz,it) = mver(it)/mvf(it)
         enddo
        enddo ! iswap
c fit with two fav. two unfavor21,ed
c using swap or not 
        do iswap=0,1
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = fub ! unfav u-> pi-
        p(4) = fdb ! fav d-> pi-+
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffit4_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsv4  (iz,j,iswap)= coef(j)
         avffsv4er(iz,j,iswap)= std(j)
         write(6,'(''rfffit4'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        do it=1,4
         ffrat4(iz,it) = mv(it)/mvf(it)
         ffrat4er(iz,it) = mver(it)/mvf(it)
         write(6,'(''ffrat4'',3i3,2f8.3)') iswap,iz,it,ffrat4(iz,it),
     >    ffrat4er(iz,it)
        enddo
        enddo ! iswap
c fit with one favored, one unfavored, and two cv params
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = 0. ! delta d / d
        p(4) = 0. ! delta u / u
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        iswap = 1.
        call mnexcm(fffitcsv_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsvcsv  (iz,j)= coef(j)
         avffsvcsver(iz,j)= std(j)
         write(6,'(''rfffitcsv'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        write(6,'(''csv chk'',f8.3)') 
     >   (4.*mdf*coef(3) + muf*coef(4)) /
     >   (4.*mdf + muf)
c fit with one favored, one unfavored, and two excl. 
c enhancement factors
        nparam=4
        npar = 4
        p(1) = fu1 ! fav u->pi+
        p(2) = fd1 ! unfav d->pi+
        p(3) = 1. ! scales up pi+ from p
        p(4) = 1. ! scales up pi- from n
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffitex_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         avffsvex  (iz,j)= coef(j)
         avffsvexer(iz,j)= std(j)
         write(6,'(''rfffitex'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
       endif ! check on rver>0
      enddo ! iz
c plot the three D(z) results for this kin
        do j=1,3
         if(j.eq.2) write(23,'(''set color blue'')')
         if(j.eq.3) write(23,'(''set color red'')')
         do iswap=0,1
          nplt = 0
          do iz=1,20
           if(avffsver(iz,j,iswap).ne.0. .and.
     >        avffsver(iz,j,iswap).lt.0.1) then
            nplt = nplt + 1
            fact=1.
            write(23,136) 
     >       0.05*(iz-0.5) + 0.005*(j)-0.02*iswap,
     >       fact*avffsv  (iz,j,iswap),
     >       fact*avffsver(iz,j,iswap)
           endif
          enddo
          if(nplt.gt.0) then
           write(23,'(''plot'')')
           if(iswap.eq.1) then
            write(23,'(''set symbol size 0.8 ; plot'')')
            write(23,'(''set symbol size 0.6 ; plot'')')
            write(23,'(''set symbol size 0.4 ; plot'')')
            write(23,'(''set symbol size 0.2 ; plot'')')
           endif
          endif
         enddo ! iswap
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(23,'(2f8.4)') z, fu1
          if(j.eq.2) write(23,'(2f8.4)') z, fd1
          if(j.eq.3) write(23,'(2f8.4)') z, fdb
         enddo
         write(23,'(''join'')')
        enddo ! j
c plot the four D(z) results for this kin
c also plot ratios to u+
c on unit-27, 3 D(z) and one CSV params on unit 29
        do j=1,4
         if(j.eq.2) write(25,'(''set color blue'')')
         if(j.eq.3) write(25,'(''set color red'')')
         if(j.eq.4) write(25,'(''set color green'')')
         if(j.eq.2) write(31,'(''set color blue'')')
         if(j.eq.3) write(31,'(''set color red'')')
         if(j.eq.4) write(31,'(''set color green'')')
         if(j.eq.2) write(26,'(''set color blue'')')
         if(j.eq.3) write(26,'(''set color red'')')
         if(j.eq.4) write(26,'(''set color white'')')
         if(j.eq.2) write(27,'(''set color blue'')')
         if(j.eq.3) write(27,'(''set color red'')')
         if(j.eq.4) write(29,'(''set color magenta'')')
         if(j.eq.2) write(28,'(''set color blue'')')
         if(j.eq.3) write(28,'(''set color red'')')
         if(j.eq.4) write(28,'(''set color green'')')
! ff4 noswap on unit 31
         nplt = 0
         do iz=1,20
          if(avffsv4er(iz,j,0).ne.0. .and.
     >       avffsv4er(iz,j,0).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,0),
     >        fact*avffsv4er(iz,j,0)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''set symbol size 1.0 ; plot'')')
          write(31,'(''set symbol size 0.8 ; plot'')')
          write(31,'(''set symbol size 0.6 ; plot'')')
          write(31,'(''set symbol size 0.4 ; plot'')')
          write(31,'(''set symbol size 0.2 ; plot'')')
         endif
         do iz=1,20
          if(avffsv4er(iz,j,1).ne.0. .and.
     >       avffsv4er(iz,j,1).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           write(25,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,1),
     >        fact*avffsv4er(iz,j,1)
           if(j.gt.1) write(26,136) 
     >       0.05*(iz-0.5) + 0.005*(j-2),
     >       avffsv4  (iz,j,1)/avffsv4(iz,1,1),
     >       avffsv4er(iz,j,1)/avffsv4(iz,1,1)
          endif
         enddo
         if(nplt.gt.0) then
          write(25,'(''plot'')')
          write(25,'(''set symbol size 1.0 ; plot'')')
          write(25,'(''set symbol size 0.8 ; plot'')')
          write(25,'(''set symbol size 0.6 ; plot'')')
          write(25,'(''set symbol size 0.4 ; plot'')')
          write(25,'(''set symbol size 0.2 ; plot'')')
          write(26,'(''plot'')')
          write(26,'(''set symbol size 0.8 ; plot'')')
          write(26,'(''set symbol size 0.6 ; plot'')')
          write(26,'(''set symbol size 0.4 ; plot'')')
          write(26,'(''set symbol size 0.2 ; plot'')')
         endif
         if(j.le.3) ioff=0
         if(j.eq.4) ioff=2
         nplt = 0
         do iz=1,20
          if(avffsvcsver(iz,j).ne.0. .and.
     >       avffsvcsver(iz,j).lt.0.4) then
           nplt = nplt + 1
           fact=1.
           write(27+ioff,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >       fact*avffsvcsv  (iz,j),
     >       fact*avffsvcsver(iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(27+ioff,'(''plot'')')
          write(27+ioff,'(''set symbol size 0.8 ; plot'')')
          write(27+ioff,'(''set symbol size 0.6 ; plot'')')
          write(27+ioff,'(''set symbol size 0.4 ; plot'')')
          write(27+ioff,'(''set symbol size 0.2 ; plot'')')
         endif
c limits on csv(x) from MRST
         if(j.eq.4) then
          v1 = 0.8 * (1. - x)**4 / sqrt(x) *
     >     (x - 0.0909)
          write(29,'(''0. '',f8.4,'' ; 1. '',f8.4,
     >     '' ; join'')') v1,v1
          v1 = -0.65 * (1. - x)**4 / sqrt(x) *
     >     (x - 0.0909)
          write(29,'(''0. '',f8.4,'' ; 1. '',f8.4,
     >     '' ; join'')') v1,v1
         endif ! j.eq.4
         nplt = 0
         do iz=1,20
          if(avffsvex er(iz,j).ne.0. .and.
     >       avffsvex er(iz,j).lt.0.2) then
           nplt = nplt + 1
           fact=1.
           write(28,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsvex   (iz,j),
     >        fact*avffsvex er(iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(28,'(''plot'')')
          write(28,'(''set symbol size 0.8 ; plot'')')
          write(28,'(''set symbol size 0.6 ; plot'')')
          write(28,'(''set symbol size 0.4 ; plot'')')
          write(28,'(''set symbol size 0.2 ; plot'')')
         endif
c plot my fit to D(z)
c xxx need to update with new fit
         do iz=6,19
          z = 0.05 * iz
          iff=1
          ipdf=1
          it=1
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
c          if(j.eq.1) write(25,'(2f8.4)') z, u1
c          if(j.eq.2) write(25,'(2f8.4)') z, d1
         enddo
c         if(j.le.2) write(25,'(''join dash'')')
c plot DSS
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(25,'(2f8.4)') z, fu1
          if(j.eq.2) write(25,'(2f8.4)') z, fd1
          if(j.eq.3) write(25,'(2f8.4)') z, fub
          if(j.eq.4) write(25,'(2f8.4)') z, fdb
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fub
          if(j.eq.4) write(31,'(2f8.4)') z, fdb
          if(j.eq.2) write(26,'(2f8.4)') z, fd1/fu1
          if(j.eq.3) write(26,'(2f8.4)') z, fub/fu1
          if(j.eq.4) write(26,'(2f8.4)') z, fdb/fu1
          if(j.eq.1) write(27,'(2f8.4)') z, fu1
          if(j.eq.2) write(27,'(2f8.4)') z, fd1
          if(j.eq.1) write(28,'(2f8.4)') z, fu1
          if(j.eq.2) write(28,'(2f8.4)') z, fd1
         enddo
         write(25,'(''join dash'')')
         write(31,'(''join dash'')')
         if(j.gt.1) write(26,'(''join dash'')')
         if(j.le.2) write(27,'(''join dash'')')
         if(j.le.2) write(28,'(''join dash'')')
c plot JAM
         do iz=6,19
          z = 0.05 * iz
          call jamff(z,q2kin(ikin),fd1,fu1)
          if(j.eq.1) write(25,'(2f8.4)') z, fu1
          if(j.eq.2) write(25,'(2f8.4)') z, fd1
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.1) write(27,'(2f8.4)') z, fu1
          if(j.eq.2) write(27,'(2f8.4)') z, fd1
          if(j.eq.2) write(26,'(2f8.4)') z, fd1/fu1
         enddo
         if(j.le.2) write(25,'(''join'')')
         if(j.le.2) write(27,'(''join'')')
         if(j.eq.2) write(26,'(''join'')')
c plot LUND ratios
c xxx need to get correct values from fitlund.f
         if(j.gt.1) then
          do iz=7,14
           z = 0.05 * iz
c           write(26,'(2f8.4)') 
c     >      z, avffsv4lund(iz,j)/
c     >      avffsv4lund(iz,1)
           write(6,'(''chklund''2i3,3f8.4)') ikin,iz, 
     >      z, avffsv4lund(iz,j),
     >      avffsv4lund(iz,1)
          enddo
c          write(26,'(''join dotdash'')')
         endif
        enddo ! j
c plot the D(z) / D(z)_fit results for this kin
c changed to be D(z) / FDSS (using averge of the
c two favored ones from DSS)
        do j=1,4
         if(j.eq.2) write(24,'(''set color blue'')')
         if(j.eq.3) write(24,'(''set color red'')')
         if(j.eq.4) write(24,'(''set color green'')')
         nplt=0
         do iz=1,20
          if(avffsv4er(iz,j,1).ne.0. .and.
     >       avffsv4er(iz,j,1).lt.0.1) then
           nplt = nplt + 1
           fact=1.
           z = 0.05 * (iz-0.5)
           call fDSS (1,1,0, Z, Q2, 
     >       fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
           call jamff(z,q2kin(ikin),fd1j,fu1j)
           if(iz.eq.10 .and. j.eq.1) write(6,
     >      '(''fdsschk4'',6f8.3)') fu1,fdb,fu1j,fub,fd1,fd1j
           if(j.eq.1.or.j.eq.4) fact=2./(fu1 + fdb)
           if(j.eq.2.or.j.eq.3) fact = 2./(fub + fd1)
c ratios to JAM
           if(j.eq.1.or.j.eq.4) fact=1./fu1j
           if(j.eq.2.or.j.eq.3) fact = 1./fd1j
           write(24,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*avffsv4  (iz,j,1),
     >        fact*avffsv4er(iz,j,1)
          endif
         enddo
         if(nplt.gt.0) then
          write(24,'(''plot'')')
          write(24,'(''set symbol size 0.8 ; plot'')')
          write(24,'(''set symbol size 0.6 ; plot'')')
          write(24,'(''set symbol size 0.4 ; plot'')')
          write(24,'(''set symbol size 0.2 ; plot'')')
         endif
        enddo ! j
c sum,diff ratios on unit 30 for this ikin
       nplt = 0
       do iz=1,20
        do it=1,6
         avmk(it)=0.
         avmkr(it)=0.
         If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avmk(it) = v1 / v2
          avmker(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 

          avmkr(it) = v1 / v2
          avmkrer(it) = sqrt(1./v2)
          write(6,'(''diffdbg'',i3,i3,i2,6f7.3)') ikin,iz,it,
     >     avm(ikin,  iz,it),avmer(ikin,  iz,it),
     >     avm(ikin+1,iz,it),avmer(ikin+1,iz,it),avmk(it),avmker(it)
         endif
        enddo ! it
        if(avmk(1).ne.0. .and. avmk(2).ne.0. .and.
     >     avmk(4).ne.0. .and. avmk(5).ne.0.) then
c diff ratio without rho subtraction (d+ - d-) / (p+ - p-)
         v1 = 1./ (avmk(1) - avmk(4)) *
     >        (avmk(2) - avmk(5))
         v2 = 1./(avmk(1) + avmker(1) - avmk(4)) *
     >        (avmk(2) - avmk(5))
         v3 = 1./(avmk(1) - avmk(4) - avmker(4)) *
     >        (avmk(2) - avmk(5))
         v4 = 1./(avmk(1) - avmk(4)) *
     >        (avmk(2) + avmker(2) - avmk(5))
         v5 = 1./(avmk(1) - avmk(4)) *
     >        (avmk(2) - avmk(5) - avmker(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
         write(30,'(3f10.4,'' 3O'')')  
c xxx    >    0.05*(iz-0.5)+0.006, v1,v6
     >    0.05*(iz-0.5)-0.006, v1,v6
! sum ratio without rho subtraction. (d+ + d-) / (p+ + p-)
         v1 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmk(5))
         v2 = 1./(avmk(1) + avmker(1) + avmk(4)) *
     >        (avmk(2) + avmk(5))
         v3 = 1./(avmk(1) + avmk(4) + avmker(4)) *
     >        (avmk(2) + avmk(5))
         v4 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmker(2) + avmk(5))
         v5 = 1./(avmk(1) + avmk(4)) *
     >        (avmk(2) + avmk(5) + avmker(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
         write(30,'(3f10.4,'' 9O'')')  
cxxx     >    0.05*(iz-0.5)+0.006, v1,v6
     >    0.05*(iz-0.5)-0.006, v1,v6
         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) then
          write(30,'(''set symbol size 1.0 ; plot'')')
          write(30,'(''set symbol size 0.8 ; plot'')')
          write(30,'(''set symbol size 0.6 ; plot'')')
          write(30,'(''set symbol size 0.4 ; plot'')')
          write(30,'(''set symbol size 0.2 ; plot'')')
       endif
c diff ratio with rho subtraction
       nplt = 0
       do iz=1,20
        do it=1,6
         avmk(it)=0.
         avmkr(it)=0.
         If(avmer(ikin,iz,it).gt.0. .and. 
     >     avmer(ikin,iz,it).lt.0.2 .and. 
     >     avmer(ikin+1,iz,it).gt.0. .and. 
     >     avmer(ikin+1,iz,it).lt.0.2) then
          v1 = avm(ikin,  iz,it)/avmer(ikin,  iz,it)**2 +
     >         avm(ikin+1,iz,it)/avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 
          avmk(it) = v1 / v2
          avmker(it) = sqrt(1./v2)
          v1 = (avm(ikin,  iz,it)-avmrho(ikin,  iz,it))/
     >          avmer(ikin,  iz,it)**2 +
     >        (avm(ikin+1,iz,it)-avmrho(ikin+1,iz,it))/
     >         avmer(ikin+1,iz,it)**2 
          v2 =                1./avmer(ikin,iz,it)**2 +
     >                        1./avmer(ikin+1,iz,it)**2 

          avmkr(it) = v1 / v2
          avmkrer(it) = sqrt(1./v2)
          write(6,'(''diffdbg'',i3,i3,i2,6f7.3)') ikin,iz,it,
     >     avm(ikin,  iz,it),avmer(ikin,  iz,it),
     >     avm(ikin+1,iz,it),avmer(ikin+1,iz,it),avmk(it),avmker(it)
         endif
        enddo ! it
        if(avmk(1).ne.0. .and. avmk(2).ne.0. .and.
     >     avmk(4).ne.0. .and. avmk(5).ne.0.) then
         v1 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5))
         v2 = 1./(avmkr(1) + avmkrer(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5))
         v3 = 1./(avmkr(1) - avmkr(4) - avmkrer(4)) *
     >        (avmkr(2) - avmkr(5))
         v4 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) + avmkrer(2) - avmkr(5))
         v5 = 1./(avmkr(1) - avmkr(4)) *
     >        (avmkr(2) - avmkr(5) - avmkrer(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
c         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
c         write(30,'(3f10.4,'' 3O'')')  0.05*(iz-0.5)-0.006, v1,v6
         v1 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5))
         v2 = 1./(avmkr(1) + avmkrer(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5))
         v3 = 1./(avmkr(1) + avmkr(4) + avmkrer(4)) *
     >        (avmkr(2) + avmkr(5))
         v4 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkrer(2) + avmkr(5))
         v5 = 1./(avmkr(1) + avmkr(4)) *
     >        (avmkr(2) + avmkr(5) + avmkrer(5))
         v6 = sqrt((v2-v1)**2 + (v3-v1)**2 + (v4-v1)**2 + (v5-v1)**2)
c         write(6,'(6f8.3)') v1,v2,v3,v4,v5,v6
c         write(30,'(3f10.4,'' 9O'')')  0.05*(iz-0.5)-0.006, v1,v6
c         nplt = nplt + 1
        endif
       enddo ! iz
       if(nplt.gt.0) write(30,'(''set symbol size 1.0 ; plot'')')
       pt=0.
       mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
       call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >     sighadp,u,ub,d,db,u1,d1,
     >     s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,2)
       v1 = (d - db) / (u - ub)
       v2 = (4.*(u-ub) - (d-db)) / 3. / (u - ub + d - db)
       v2 = v2 * 5. * (u + ub + d + db) /
     >           (4.*(u + ub) + (d + db))
c v2 = (4u - d) / 3 / (u+d) * 5 * (u+d) / (4u + d)
c v2 = (5/3) (4u-d) / (4u + d)
       v2 = 1./v2
c (3/5) (4u + d) / (4u -d)      
        write(30,'(''0.2 '',f10.4,'' ; 1.0 '',f10.4)') 
     >   v2,v2
       write(30,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
c diff ratio with d - db = 0.
c       write(22,'(''set pattern .02 .02 .02 .02'')')
c       write(22,'(''0.2 1.66 ; 1.0 1.66 ; join 1'')') 
c sum rattio
       write(30,'(''0.2 1.00 ; 1.0 1.00'')') 
       write(30,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
c simc model using 2 Frag Fun.
       do iff=2,4,2
       if(iff.eq.2) write(30,'(''set pattern .05 .05 .05 .05'')')
       if(iff.eq.4) write(30,'(''set pattern .02 .02 .02 .02'')')
       if(iff.eq.3) write(30,'(''set pattern .01 .01 .01  .01 '')')
       do iz=4,18
        z =0.05*iz
        ipdf=1
        call simcmodel(x,q2,z,pt,phicm,mmpi2,1,
     >   sighad1,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,2,
     >   sighad2,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,4,
     >   sighad4,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
        call simcmodel(x,q2,z,pt,phicm,mmpi2,5,
     >   sighad5,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
         write(30,'(2f8.4)') z,1./(sighad1-sighad4)*
     >    (sighad2 - sighad5)
         rsv(iz) = 1./(sighad1 + sighad4)*
     >    (sighad2 + sighad5)
        enddo
        if(iff.eq.4) then
         write(30,'(''join '')')
        else
         write(30,'(''join pattern'')')
        endif
        do iz=3,18
         z =0.05*iz
         write(30,'(2f8.4)') z,rsv(iz)
        enddo
        if(iff.eq.4) then
         write(30,'(''join '')')
        else
         write(30,'(''join pattern'')')
        endif
       enddo ! iff
c       write(22,'(8f7.2)') x,q2,v1,v2
      enddo ! ikin

c debug simcmodel
        do iff=1,2
         do ipdf=1,2
          do j=1,5
           do iz=1,3
             z= 0.25 * iz
             x = 0.2 + 0.1 * j
             q2=3.
             pt=0.
             phicm=0.
             it=1
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v1 = (4 * u * u1 + d * d1)/
     >        (4.*u + d) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v1,v1/dsigdz,
c     >        sighad1/dsigdz,(ub + db + s + sb)/(u+d)
             it=2
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v2 = (4 * (u+d)* u1 + (d+u) * d1)/
     >        (4.*(u+d) + (d+u)) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v2,v2/dsigdz,
c     >        sighad2/dsigdz,(ub + db + s + sb)/(u+d)
             it=4
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v3 = (4 * u * u1 + d * d1)/
     >        (4.*u + d) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v3,v3/dsigdz,
c     >        sighad3/dsigdz,(ub + db + s + sb)/(u+d)
             it=5
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             v4 = (4 * (u+d)* u1 + (d+u) * d1)/
     >        (4.*(u+d) + (d+u)) / z
c             write(6,'(3i2,10f7.3)') iff,ipdf,it,
c     >        u,d,u1,d1,dsigdz,v4,v4/dsigdz,
c     >        sighad4/dsigdz,(ub + db + s + sb)/(u+d)
             v5 = (d - db) / (u - ub)
             write(6,'(''dbgs'',2i2,f6.2,10f7.3)') iff,ipdf,x,
     >        sighad1, sighad2, sighad3, sighad4,
     >        (sighad2 - sighad4) / (sighad1 - sighad3),
     >        (sighad2 + sighad4) / (sighad1 + sighad3),
     >        0.6 * (1 + v5/4.)/(1 - v5/4.),v5,
     >        (v2 - v4) / (v1 - v3)
            enddo ! iz
           enddo ! j (x)
          enddo ! iff
         enddo ! ipdf

       do j=3,8
        z=0.1*j
        q2=3.
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
         write(6,'(''FDSS'',f6.2,6F7.4)') z,
     >      fU1, fUB, fD1, fDB, fS1, fSB
       enddo

c       goto 999
ccccccccccccccccc
c look at shuo's results. OLD ones
      open(unit=8,file='shuoratlarge.txt')
      do i=1,386
       read(8,*) q2n,q2,q2er,xn,x,xer,zn,z,zer,j,
     >   rd,rder,rdrho
       pt=0.1
       phicm=0.
       it=2
       mmpi2 = am**2 + q2*(1./x-1.)*(1.-z)
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighadp,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
       it=5
       write(6,'(''sh mx'',5f7.3)') q2,x,z,mmpi2
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighadm,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c for test
c       sighadp = sighadp * 0.95
       rdm = (4*sighadm - sighadp) /
     >       (sighadp - sighadm)
       diff = (rd/rdm - 1.0) /
     >       (rder/rdm)
       diff = max(-5.,min(4.999,diff))
       k = int((diff + 5.)/0.25) +1
       rdh(k,1) = rdh(k,1)+1
       diff = (rdrho/rdm - 1.0) /
     >       (rder/rdm)
       diff = max(-5.,min(4.999,diff))
       k = int((diff + 5.)/0.25) +1
       rdh(k,2) = rdh(k,2)+1
        write(6,'(''shuo'',3f6.3,6f7.3)') q2,x,z,
     >  rd/rdm,rder/rdm,rd,rder,rdrho,rdm
c average diff versus various variables
       diff = (rdrho/rdm - 1.0) /
     >       (rder/rdm)
       k=max(1,min(10,int((q2-3.)/3.5*10)+1))
       sq2(k,1) = sq2(k,1) + diff
       sq2(k,2) = sq2(k,2) + 1.
       k=max(1,min(10,int((x-0.25)/0.35*10)+1))
       sx(k,1) = sx(k,1) + diff
       sx(k,2) = sx(k,2) + 1.
       k=max(1,min(10,int((z-0.3)/0.4*10)+1))
       sz(k,1) = sz(k,1) + diff
       sz(k,2) = sz(k,2) + 1.   
       k=max(1,min(10,int((mmpi2-2.)/5*10)+1))
       smx(k,1) = smx(k,1) + diff
       smx(k,2) = smx(k,2) + 1.   
      enddo

      do k=1,10
       write(6,'(''sh kin'',4(f6.1,f5.0))') 
     >  sq2(k,1)/max(1.,sq2(k,2)),sq2(k,2),
     >  sx(k,1)/max(1.,sx(k,2)),sx(k,2),
     >  sz(k,1)/max(1.,sz(k,2)),sz(k,2),
     >  smx(k,1)/max(1.,smx(k,2)),smx(k,2)
      enddo
      do k=1,10
       write(6,'(''pb'',5(f6.1,f5.0))') 
     >  pq2(k,1)/max(1.,pq2(k,2)),pq2(k,2),
     >  px(k,1)/max(1.,px(k,2)),px(k,2),
     >  pz(k,1)/max(1.,pz(k,2)),pz(k,2),
     >  pw(k,1)/max(1.,pw(k,2)),pw(k,2),
     >  pmx(k,1)/max(1.,pmx(k,2)),pmx(k,2) 
      enddo
      do k=1,10
       write(6,'(''pb'',e12.4,f7.0)') 
     >  pq2(k,1),pq2(k,2)
      enddo

      open(unit=22,file='ptmsrat.top')
      write(22,'(''set device postscript'')')
      do i=1,2
       do j=1,3
        if(j.eq.1) it=1
        if(j.eq.2) it=4
        if(j.eq.3) it=5
        write(22,222) j,i,it,i
 222    format('set window x ',i1,' of 3 y ',i1,' of 2'/
     >   'title top ',1h','it,i=',2i3,1h')
        do k=1,40
         if(i.eq.1) write(22,'(f7.2,i6)') 
     >    -7. + 14./40.*(k-0.5),
     >    srath(k,it)
         if(i.eq.2) write(22,'(f7.2,i6)') 
     >    -7. + 14./40.*(k-0.5),
     >    sratrhoh(k,it)
        enddo
        write(22,'(''hist'')')
       enddo
      enddo
      write(22,'(''new frame'')')
      do i=1,2
       write(22,223) i
 223   format('set window x ',i1,' of 2 y 1 of 1')
       do k=1,40
        write(22,'(f7.2,i6)') 
     >    -5. + 0.25*(k-0.5),rdh(k,i)
       enddo
       write(22,'(''hist'')')
      enddo

      close(unit=22)

      write(6,'(''bnpt='',i6)') bnpt
c add in inclusiev
      bnpt0 = bnpt
      open(unit=44,file='ptm.incltxt')
      do i=1,100
       read(44,'(4x,6f10.4)') q2,w2,f1p,f2p,f1d,f2d
       xx = (w2 - am**2) / q2 
       x = 1./(xx + 1.)
       w2chk = am**2 + q2 * (1/x -1)
       z=0.4
       pt=0.
       phicm=0.
       mmpi2 = 4.
       it=1
       call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,u,ub,d,db,u1,d1,
     >   s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
       if(x.lt.0.7 .and. q2.gt.2.5) then
        bnpt = bnpt + 1
        bq2v(bnpt) = q2
        bxv(bnpt) = x
        buv(bnpt) = u
        bubv(bnpt) = ub
        bdv(bnpt) = d
        bdbv(bnpt) = db
        bsv(bnpt) = s
        bsbv(bnpt) = sb
c took out f1f2in21
c       byv(bnpt) = f2d / f2p
c use this instad
        call FNP_NMC(X,Q2,ratt)
        byv(bnpt) = 1. + ratt
        byerv(bnpt) = byv(bnpt) * 0.03
        bkin(bnpt) = 0 ! special code
        delu=0.
        deld=0.
        sum_sqp = qu**2*(u+ub) + qd**2*(d+db) + 
     >     qs**2*(s+sb)
        sum_sqn = qu**2*(d + deld + db) + 
     >         qd**2*(u + delu + ub) + qs**2*(s+sb)
        write(6,'(''f1f2'',3f5.2,5f7.3)') x,q2,w2,
     >   f2p, f2d, f2d/f2p, 1.+sum_sqn/sum_sqp, 1.+ratt
       endif
      enddo

c take out incl
      bnpt = bnpt0

c       goto 99

      if(usewide) then
       open(unit=22,file='ptmfitsw.top')
       open(unit=23,file='ptmfitpw.top')
      else
       open(unit=22,file='ptmfits.top')
       open(unit=23,file='ptmfitp.top')
      endif
      write(22,'(''set device postscript'')')
      write(23,'(''set device postscript'')')

c start of huge loop over ifit
      do ifit=1,3
      rhofact = 0.0
      rhofactp = 0.0
      if(ifit.eq.2) rhofact = 1.
c      if(ifit.eq.3) rhofact = 10.
      if(ifit.eq.2) rhofactp = 1.
c      if(ifit.eq.3) rhofactp = 0.
      ipdf = 1
      if(ifit.eq.3) ipdf = 2
      write(6,'(''using ipdf, rhofact = '',2i2,f5.1)') 
     >  ifit,ipdf,rhofact

      do j=1,19
        write(pname(j),'(''P'',i1)') j
      enddo

c fit each kin and target individually
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo
      lastframeused = .false.
      do ikinfit=1,32,2
       do iz=6,18
        if(lastframeused.and.ifit.eq.1) write(22,'(''new frame'')')
        if(lastframeused.and.ifit.eq.1) write(23,'(''new frame'')')
        lastframeused = .false.
        do itt=1,4
         itfit=itt
         if(itt.gt.2) itfit = itt + 1
         zfit = 0.05 * (iz - 0.5)
c         nparam=4
c         npar = 4
         nparam=3
         npar = 3
         p(1) = 0.1
         p(2) = 0.2
         p(3) = 0.0
         p(4) = 0.
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         ncallsfit=0
         call sfit_fcn(npar,grad,chi2,p,ierflg,futil)
         if(npts.gt.40) then
          write(6,'(''initial chi2'',2i4,f6.2,i6,f8.2)') 
     >     ikinfit,itfit,zfit,npts,
     >     chi2/float(npts)
          arglis(1)=0
          call mnexcm(sfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
          do j=1,nparam
           call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
           write(6,'(''sf '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
           scsv((ikinfit+1)/2,itt,iz,j,ifit)= coef(j)
           scsver((ikinfit+1)/2,itt,iz,j,ifit)= std(j)
          enddo
          call sfit_fcn(npar,grad,chi2,coef,ierflg,futil)
          write(6,'(''final chi2'',2i4,f6.2,i6,f8.2)') 
     >     ikinfit,itfit,zfit,npts,
     >     chi2/float(npts)
          write(16,'(i2,i2,f5.2,i4,9f7.3)') ikinfit,itt,
     >     zfit,npts,chi2/float(npts),
     >     (coef(j),std(j),j=1,4)
          if((ikinfit.eq.1 .or. ikinfit.eq.5) .and.
     >     ifit.eq.1.and.iz.gt.4 .and. iz.lt.18) then
           ix=1
           iy=1
           if(itt.eq.2.or.itt.eq.4) ix=2
           if(itt.eq.1.or.itt.eq.2) iy=2
           pt2 = 0.
           ymax = 1.5 * (coef(1)/zfit  ) / coef(2) * 
     >        exp(-pt2/coef(2)) 
           pt2 = 0.7**2
           ymin = 0.8 * (coef(1)/zfit  ) / coef(2) * 
     >        exp(-pt2/coef(2)) 
           x = (xkin(ikinfit) + xkin(ikinfit+1))/2.
           q2 = (q2kin(ikinfit) + q2kin(ikinfit+1))/2.
           write(22,833) ix,iy,ttit(itt),x,q2,zfit,
     >      npts,chi2/float(npts),ymin,ymax
 833       format('set window x ',i1,' of 2 y ',
     >                            i1,' of 2'/
     >      'set color white'/ 
     >      'title top ',1h',a2,' x=',f4.2,' Q2=',f3.1,
     >      ' z=',f4.2,' npt=',i3,' chi2/df=',f5.1,1h'/
     >      'title bottom ',1h','phi*',1h'/
     >      'title left ',1h','Multiplicity',1h'/
     >      'set limits x 0. 6.3 y ',2f10.3/
     >      'set scale y log ; set bar size 0.'/
     >      'set symbol 9O size 1.0 ; set order x y dy')
           write(23,834) ix,iy,ttit(itt),x,q2,zfit,
     >      npts,chi2/float(npts),ymax/3.,ymax
 834       format('set window x ',i1,' of 2 y ',
     >                            i1,' of 2'/
     >      'title top ',1h',a2,' x=',f4.2,' Q2=',f3.1,
     >      ' z=',f4.2,' npt=',i4,' chi2/df=',f5.1,1h'/
     >      'title bottom ',1h','phi*',1h'/
     >      'title left ',1h','Multiplicity',1h'/
     >      'set limits x 0. 6.3 y ',2f10.3/
     >      'set bar size 0.'/
     >      'set symbol 9O size 1.0 ; set order x y dy')
           do ipt=1,10
            if(ipt.eq.1) write(22,'(''set color white'')')
            if(ipt.eq.2) write(22,'(''set color magenta'')')
            if(ipt.eq.3) write(22,'(''set color red'')')
            if(ipt.eq.4) write(22,'(''set color green'')')
            if(ipt.eq.5) write(22,'(''set color blue'')')
            if(ipt.eq.6) write(22,'(''set color cyan'')')
            if(ipt.eq.7) write(22,'(''set color white'')')
            if(ipt.eq.8) write(22,'(''set color magenta'')')
            if(ipt.eq.9) write(22,'(''set color red'')')
            if(ipt.eq.10) write(22,'(''set color green'')')
            if(ipt.eq.11) write(22,'(''set color blue'')')
            nplt=0
            do iphi=1,14,2
             phimin = 2.*3.1415928/15.*(iphi-1)
             phimax = 2.*3.1415928/15.*(iphi+1)
             sum1=0.
             sum2=0.
             do jj=1,npts
c              if(spt(jj).gt.0.06*(ipt-1) .and.
c     >           spt(jj).le.0.06*(ipt+1).and.
c use all pt bins
              if(spt(jj).gt.0.06*(ipt-0.5) .and.
     >           spt(jj).le.0.06*(ipt+0.5).and.
     >           sphi(jj).gt.phimin .and.
     >           sphi(jj).lt.phimax) then
               sum1 = sum1 + sy(jj)/syer(jj)**2
               sum2 = sum2 + 1./syer(jj)**2
              endif
             enddo
             if(sum2.gt.0.) then
              sum1 = sum1 / sum2
              sum2 = 1./sqrt(sum2)
              write(22,'(3f7.3)') (phimin + phimax)/2.,sum1,sum2
              nplt = nplt + 1
             endif
            enddo
            if(nplt.gt.0) then
             lastframeused = .true.
             write(22,'(''plot'')')
             pt = 0.06*(ipt)
             pt2 = pt**2
             do jj=1,21
              phi = 2.*3.1415928*(jj-1)/20.
              sig = (coef(1)/zfit  ) / coef(2) * 
     >         exp(-pt2/coef(2)) *
     >         (1. + coef(3) * pt * cos(phi))
              write(22,'(3f8.3)') phi,sig
             enddo
             write(22,'(''join'')')
            endif
           enddo

c plot mult versus phi averaged over pt<0.25
           do iphi=1,15
            sig = 0.
            siger = 0.
            do jj=1,npts
             if(sphi(jj).gt. 2.*3.1415928/15.*(iphi-1) .and.
     >          sphi(jj).le. 2.*3.1415928/15*iphi .and. 
     >          spt(jj).le.0.25) then
              sig = sig + sy(jj)/syer(jj)**2
              siger = siger + 1./syer(jj)**2
             endif
            enddo
            sig = sig / siger
            siger = 1./sqrt(siger)
            if(siger.ne.0.) write(23,'(3f10.4)') 
     >       2.*3.1415928/15.*(iphi-0.5),sig,siger
           enddo
           write(23,'(''plot'')')
           pt = 0.12
           pt2 = pt**2
           do jj=1,21
            phi = 2.*3.1415928*(jj-1)/20.
            sig = (coef(1)/zfit  ) / coef(2) * 
     >       exp(-pt2/coef(2)) /2. / pi **
     >       (1. + coef(3) * pt * cos(phi))
            write(23,'(3f8.3)') phi,sig
           enddo
           write(23,'(''join'')')
          endif
         endif
        enddo
       enddo
      enddo
c xxx for test
      goto 77
! get csv from d only using modified fdds plus Geiger and both total integral
! and values at pt=0 for ratio
      do ikin=1,3
       q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
       x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
       do iz=5,18
        z= 0.05*(iz-0.5)
        wp2 = am**2 + q2*(1/x - 1)*(1-z)
        Mm = scsv(ikin,4,iz,1,ifit)
        Mp = scsv(ikin,2,iz,1,ifit)
        Mmer = scsver(ikin,4,iz,1,ifit)
        Mper = scsver(ikin,2,iz,1,ifit)
        Bm = scsv(ikin,4,iz,2,ifit)
        Bp = scsv(ikin,2,iz,2,ifit)
        if(Mm.ne.0. .and. Mp.ne.0.) then
         pt=0.
         phicm=0.
         it=4
         call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >    sighad,u,ub,d,db,u1,d1,
     >    s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
         Rseans(ikin,ikin) = 5. * (ub + db) / 
     >   (u - ub + d - db)
         Rseas(ikin,iz) = (s + sb) / (u - ub + d - db)
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        rnew = fd1 / fu1
c for test
          rnew = rnew * 1.35
c swith to geiger: works better
         call rgeiger(z,rnew)
         Dprop = (1. - rnew) / (1. + rnew)
         mprat = Mm / Mp 
         RD = (4.*mprat - 1.) / (1. - mprat)
         mprat = (Mm + Mmer) / Mp 
         RD1 = (4.*mprat - 1.) / (1. - mprat)
         mprat = Mm / (Mp + Mper) 
         RD2 = (4.*mprat - 1.) / (1. - mprat)
         dRd = sqrt((rd1-rd)**2 + (rd2-rd)**2)
         CSVprop(ikin,iz,ifit,1) = 2.5 + rseans(ikin,iz) + 
     >    rseas(ikin,iz) - Dprop * (2.5 + RD)
c this is delta u / (u + d) if delta u = -1 *delta d
c ignoring sea
         chk1 = (4*rnew + 1 - mprat * (4 + rnew)) /
     >         (4*rnew - 1 - mprat * (4 - rnew))
c this is -4/3 (delta u - delta d) / (u + d)
c or -8/3 delta u / (u+d)
         chk2 = 2.5  - Dprop * (2.5 + RD)
         CSVproper(ikin,iz,ifit,1) = Dprop * dRD
         mprat = (Mm / Bm) / (Mp / Bp) 
         RD = (4.*mprat - 1.) / (1. - mprat)
         CSVprop(ikin,iz,ifit,2) = 2.5 + rseans(ikin,iz) + 
     >    rseas(ikin,iz) - Dprop * (2.5 + RD)
         CSVproper(ikin,iz,ifit,2) = 
     >    CSVproper(ikin,iz,ifit,1) 
         write(6,'(''csv'',2i3,6f7.2)') ikin,iz,
     >    (csvprop(ikin,iz,ifit,k),
     >    CSVproper(ikin,iz,ifit,k),k=1,2),chk1,
     >    -3./8.*chk2
         do k=1,2
          if(wp2 .gt. 3.0) then
           avcsv(ikin,ifit,k) = avcsv(ikin,ifit,k) +
     >      csvprop(ikin,iz,ifit,k) /
     >      csvproper(ikin,iz,ifit,k)**2
           avcsver(ikin,ifit,k) = avcsver(ikin,ifit,k) +
     >      1./csvproper(ikin,iz,ifit,k)**2
          endif
         enddo
        endif
       enddo ! iz
       do k=1,2
        avcsv(ikin,ifit,k) = avcsv(ikin,ifit,k) /
     >    avcsver(ikin,ifit,k)
        avcsver(ikin,ifit,k) = 1. /
     >   sqrt(avcsver(ikin,ifit,k))
       enddo
       write(6,'(''avcsv'',2i3,4f8.3)') ikin,ifit,
     >  (avcsv(ikin,ifit,k),avcsver(ikin,ifit,k),k=1,2)
      enddo ! ikin

! get csv, fav, fav/unfav from all four it cases
c not fitting d/u because results are not sensitive to this:
c would have to constrain inclusive to remain unchanged
      do ikinx=1,8
       ikin = ikinx
       if(ikinx.eq.8) ikin=15
       q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
       x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
c get pdfs 
       do iz=5,18
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,mu,mub,md,mdb,u1,d1,
     >   ms,msb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        z = 0.05*(iz-0.5)
        wp2 = am**2 + q2*(1/x - 1)*(1-z)
        ngood=0
        BigS = 2.*(ms + msb) / (mu + md)
        BigT = (mub + mdb) / (mu + md)
        uplusd(ifit,ikinx) = mu + md
        do it=1,4
         Mfit(it) = scsv(ikin,it,iz,1,ifit) * 2. * 3.1415928
c try using M0 * b instead of M0. This makes CSV even bigger
c and also diff ratio worse too
c     >     * scsv(ikin,it,iz,2,ifit) / 0.2
         Mfiter(it) = scsver(ikin,it,iz,1,ifit) * 2. * 3.1415928
c     >     * scsv(ikin,it,iz,2,ifit) / 0.2
         if(Mfiter(it).ne.0 .and. 
     >    Mfiter(it).lt.0.1) ngood = ngood + 1
c values using average mult.
         mfitp(it)=0.
         mfitper(it)=0.
          itt=it
          if(it.gt.2) itt=it+1
          if(avmer(2*ikin-1,iz,itt).ne.0. .and.
     >      avmer(2*ikin  ,iz,itt).ne.0.) then
           do ikinp = 2*ikin-1, 2*ikin
            mfitp(it) = mfitp(it) + 
     >      (avm(ikinp,iz,itt) - 
     >       rhofactp * avmrho(ikinp,iz,itt)) /
     >       avmer(ikinp,iz,itt)**2
            write(6,'(''dbg'',4i3,f6.2,4f8.3)') ikinp,iz,itt,
     >       ifit,rhofactp,avm(ikinp,iz,itt), 
     >       avmer(ikinp,iz,itt),
     >       rhofactp * avmrho(ikinp,iz,itt),
     >       avmrho(ikinp,iz,itt)
            mfitper(it) = mfitper(it) + 1./
     >       avmer(ikinp,iz,itt)**2
           enddo
           mfitp(it) = mfitp(it) / mfitper(it)
           mfitper(it) = 1./sqrt(mfitper(it))
         endif
         write(6,'(''mfitp'',3i3,8f7.3)') ikin,iz,it,
     >    mfit(it),mfiter(it),mfitp(it),mfitper(it),
     >    avm  (2*ikin-1,iz,itt),
     >    avmer(2*ikin-1,iz,itt),
     >    avm  (2*ikin  ,iz,itt),
     >    avmer(2*ikin  ,iz,itt) 
 
c  for test, replace with averaged mult.
c         mfit(it) = mfitp(it)
c         mfiter(it) = mfitper(it)
        enddo !it

c for test
c this lowers funfav / ffav, but doesn't change delta u
c        mfit(3) = mfit(3) * 0.9
c        mfit(4) = mfit(4) * 0.9
c for test can change inclusive d/p ratio
        fact = 1.00
c        fact = 1.05
        mfit(2) = mfit(2) * fact
        mfiter(2) = mfiter(2) * fact
        mfit(4) = mfit(4) * fact
        mfiter(4) = mfiter(4) * fact
c  for test, increase d pi+ only
        fact = 1.00
        mfit(2) = mfit(2) * fact
        mfiter(2) = mfiter(2) * fact
c  for test, increase d pi+ only
        fact = 1.00
        mfitp(2) = mfitp(2) * fact
        mfitper(2) = mfitper(2) * fact
        if(ngood.eq.4) then
c get starting z*FF/z (matches mfit which is also scaled by z)
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
c         nparam=4
c         npar = 4
         nparam=3
         npar = 3
c for test, try varying csv by hand
         do j=1,10
          p(1) = fu1
          p(2) = fd1 / fu1
          p(3) = 0.05 * float(j-5)
          p(4) = -p(3)
          call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
          write(6,'(''mt'',10f7.3)') 
     >     p(3),chi2,(mfit(it),mfitr(it),it=1,4)
         enddo
! fit using mfit
         mfitflag = 0
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
c         p(4) = 0.0 ! deltd
c         p(3) = md / mu ! d/u
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit I'',2i2,i3,f8.2/8f8.3/3f8.3/6f7.3)') 
     >     itfit,ikin,iz,chi2,(mfit(it),mfitr(it),it=1,4),
     >     (p(it),it=1,3),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         arglis(1)=0
         call mnexcm(mfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
          mcsv(ifit,ikinx,iz,j)= coef(j)
          mcsver(ifit,ikinx,iz,j)= std(j)
         enddo
         fu1sv(ifit,ikinx,iz)= fu1
         fd1sv(ifit,ikinx,iz)= fd1
         fact=1.
c delu, deld from sum, diff
         call getcsv(delu,deluer,deld,delder,doverp,
     >     doverpp, doverpper,fact)
         mcsv(ifit,ikinx,iz, 9)= delu / (mu + md)
         mcsver(ifit,ikinx,iz, 9)= deluer / (mu + md)
         mcsv(ifit,ikinx,iz,10)= deld / (mu + md)
         mcsver(ifit,ikinx,iz,10)= delder / (mu + md)

! fit using mfitp now
         mfitflag = 1
         nparam=3
         npar = 3
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
c         p(4) = 0.0 ! deltd
c        p(4) = md / mu         ! d/u
         call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c         call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
         call mfit_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit Ip'',2i2,i3,f8.2/8f8.3/3f8.3/6f7.3)') 
     >     itfit,ikin,iz,chi2,(mfit(it),mfitr(it),it=1,4),
     >     (p(it),it=1,3),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         arglis(1)=0
         call mnexcm(mfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
          mcsvp(ifit,ikinx,iz,j)= coef(j)
          mcsvper(ifit,ikinx,iz,j)= std(j)
         enddo
c influeance of d/u ratio (usng mfitp)
         do j=1,5
          fact = 0.7 + 0.1*j
          call getcsv(delu,deluer,deld,delder,doverp,
     >    doverpp, doverpper,fact )
          write(6,'(''duchk'',3i3,7f7.3)') ifit,ikinx,iz,
     >     fact,delu,deluer,deld,delder,doverpp/doverp,doverp
         enddo
         fact=1.
         call getcsv(delu,deluer,deld,delder,doverp,
     >    doverpp, doverpper,fact)
         mcsvp(ifit,ikinx,iz, 9)= delu / (mu + md)
         mcsvper(ifit,ikinx,iz, 9)= deluer / (mu + md)
         mcsvp(ifit,ikinx,iz,10)= deld / (mu + md)
         mcsvper(ifit,ikinx,iz,10)= delder / (mu + md)
         dopsv(ifit,ikinx,iz)= doverp          
         mcsvp(ifit,ikinx,iz,11)= doverpp          
         mcsvper(ifit,ikinx,iz,11)= doverpper
         write(6,'(''cmpcsv'',8f7.3)') delu,deluer,deld,delder,
     >    coef(3)/delu,std(3)/deluer,coef(4)/deld,std(4)/delder
         write(6,'(''doverp'',3i3,3f8.3)') ifit,ikinx,iz,
     >    doverp, doverpp, doverpper
         avdelu(ifit,ikinx) = avdelu(ifit,ikinx)  +
     >    delu / (mu + md) / (deluer / (mu + md))**2
         avdeluer(ifit,ikinx) = avdelu(ifit,ikinx)  +
     >    1. / (deluer / (mu + md))**2
         avdeld(ifit,ikinx) = avdeld(ifit,ikinx)  +
     >    deld / (mu + md) / (delder / (mu + md))**2
         avdelder(ifit,ikinx) = avdeld(ifit,ikinx)  +
     >    1. / (delder / (mu + md))**2

c also get delu from deuteron only for delu + deld = 0
c this gives delu / (u + d), using mfit
         call rgeiger(z,rnew)
         mprat = mfit(4) / mfit(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mprat = (mfit(4) + mfiter(4)) / mfit(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >             (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mprat = mfit(4) / (mfit(2) + mfiter(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >             (4         *rnew - 1 - 
     >    mprat * (4 - rnew))
         mcsv(ifit,ikinx,iz,4) = chk1
         mcsver(ifit,ikinx,iz,4) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c also get delu from deuteron only for delu + deld = 0
c this gives delu / (u + d), using mfitp
         a = 1/4
         call rgeiger(z,rnew)
         mprat = mfitp(4) / mfitp(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mcsvp(ifit,ikinx,iz,4) = chk1
         mcsvper(ifit,ikinx,iz,4) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c same thing but with DSS
         call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
         rnew = fd1 / fu1
         mprat = mfitp(4) / mfitp(2)
         chk1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = ((4 + BigS)*rnew + 1  + BigT*(rnew + 4) - 
     >    mprat * (4 + rnew + BigS*rnew + BigT*(4 + rnew))) /
     >          (4.*a*rnew - 1 - mprat * (4.*a - rnew))
         mcsvp(ifit,ikinx,iz,12) = chk1
         mcsvper(ifit,ikinx,iz,12) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
c get fu / fd from proton data only using mfit
         mprat = mfit(3) / mfit(1)
         chk1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.* mub + ms + msb))
         chk2 = (mprat * (4.*mu + mdb) - md*1.1 - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md*1.1 + 4.* mub + ms + msb))
         write(6,'(''fup 1'',3f8.3)') mprat,chk1,chk2
         mprat = (mfit(3) + mfiter(3)) / mfit(1)
         chk1er1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfit(3) / (mfit(1) + mfiter(1))
         chk1er2 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsv(ifit,ikinx,iz,5) = chk1
         mcsver(ifit,ikinx,iz,5) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,5),mcsver(ifit,ikinx,iz,5)
c get fu / fd from proton data only using mfitp
         mprat = mfitp(3) / mfitp(1)
         chk1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.* mub + ms + msb))
         chk2 = (mprat * (4.*mu + mdb) - md*1.1 - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md*1.1 + 4.* mub + ms + msb))
         write(6,'(''fup 1'',3f8.3)') mprat,chk1,chk2
         mprat = (mfitp(3) + mfitper(3)) / mfitp(1)
         chk1er1 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
c         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfitp(3) / (mfitp(1) + mfitper(1))
         chk1er2 = (mprat * (4.*mu + mdb) - md - 4.*mub) /
     >     (4.*mu + mdb + ms + msb - 
     >      mprat * (md + 4.*mub + ms + msb))
c         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsvp(ifit,ikinx,iz,5) = chk1
         mcsvper(ifit,ikinx,iz,5) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fupp '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,5),mcsvper(ifit,ikinx,iz,5)
c get fu / fd from deuteron data only assuming csv=0 using mfit
         mprat = mfit(4 ) / mfit(2)
         chk1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 1'',2f8.3)') mprat,chk1
         mprat = (mfit(4) + mfiter(4)) / mfit(2)
         chk1er1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfit(4) / (mfit(2) + mfiter(2))
         chk1er2 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsv(ifit,ikinx,iz,6) = chk1
         mcsver(ifit,ikinx,iz,6) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,6),mcsver(ifit,ikinx,iz,6)
c get fu / fd from deuteron data only w/csv=0 using mfitp
         mprat = mfitp(4 ) / mfitp(2)
         chk1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
c         write(6,'(''fup 1'',2f8.3)') mprat,chk1
         mprat = (mfitp(4) + mfitper(4)) / mfitp(2)
         chk1er1 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
c         write(6,'(''fup 2'',2f8.3)') mprat,chk1er1
         mprat = mfitp(4) / (mfitp(2) + mfitper(2))
         chk1er2 = (mprat * (4.*(mu+md) + (mdb+mub)) - 
     >     (md+mu) - 4.*(mub+mdb)) /
     >     (4.*(mu+md) + (mdb+mub) + 2.*ms + 2.*msb - 
     >      mprat * ((md+mu) + 4.*(mub+mdb) + 2.*ms + 2.*msb))
         write(6,'(''fup 3'',2f8.3)') mprat,chk1er2
         mcsvp(ifit,ikinx,iz,6) = chk1
         mcsvper(ifit,ikinx,iz,6) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2)
         write(6,'(''fup '',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,6),mcsvper(ifit,ikinx,iz,6)
c 5-param fit using mfitp. Includes d/u with constraint on 
c inclusive p/d remaining unchanged
         nparam=5
         npar = 5
         p(1) = fu1
         p(2) = fd1 / fu1
         p(3) = 0.0 ! deltu
         p(4) = 0.0 ! deltd
         p(5) = 1.0 ! d/p
         call mfit5_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit5 Ip'',2i2,i3,f8.2/8f8.3/5f8.3/6f7.3)') 
     >     ifit,ikin,iz,chi2,(mfitp(it),mfitr(it),it=1,4),
     >     (p(it),it=1,5),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
c get better guess for p1
         p(1) = p(1) * 
     >   (mfitp(1)+mfitp(2)+mfitp(3)+mfitp(4)) / 
     >   (mfitr(1)+mfitr(2)+mfitr(3)+mfitr(4))
         call mfit5_fcn(npar,grad,chi2,p,ierflg,futil)
         write(6,'(''mfit5 IIp'',2i2,i3,f8.2/8f8.3/5f8.3/6f7.3)') 
     >     ifit,ikin,iz,chi2,(mfitp(it),mfitr(it),it=1,4),
     >     (p(it),it=1,5),mu,mub/mu,md/mu,mdb/mu,ms/mu,msb/mu
         call mnparm( 1,"P1 ",p(1), 0.001D0,zero,zero,ierflg)
         call mnparm( 2,"P2 ",p(2), 0.001D0,zero,zero,ierflg)
         call mnparm( 3,"P3 ",p(3), 0.001D0,zero,zero,ierflg)
         call mnparm( 4,"P4 ",p(4), 0.001D0,zero,zero,ierflg)
         call mnparm( 5,"P5 ",p(5), 0.001D0,zero,zero,ierflg)
         arglis(1)=0
         call mnexcm(mfit5_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
         do j=1,nparam
          call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
          write(6,'(''mfit5 '',2i2,8f7.3)') ifit,j,
     >     coef(j),std(j)
           mcsv5p(ifit,ikinx,iz,j)= coef(j)
          mcsv5per(ifit,ikinx,iz,j)= std(j)
         enddo
c get sum and diff ratios using mfit
         chk1 = (mfit(2) - mfit(4))/
     >          (mfit(1) - mfit(3))
         chk1er1 = (mfit(2)+mfiter(2) - mfit(4))/
     >             (mfit(1)           - mfit(3))
         chk1er2 = (mfit(2)           - mfit(4))/
     >             (mfit(1)+mfiter(1) - mfit(3))
         chk1er3 = (mfit(2) - mfit(4) - mfiter(4))/
     >             (mfit(1) - mfit(3))
         chk1er4 = (mfit(2) - mfit(4))/
     >             (mfit(1) - mfit(3) - mfiter(3))
         mcsv(ifit,ikinx,iz,7) = chk1
         mcsver(ifit,ikinx,iz,7) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fdif'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,7),mcsver(ifit,ikinx,iz,7)
         chk1 = (mfit(2) + mfit(4))/
     >          (mfit(1) + mfit(3))
         chk1er1 = (mfit(2)+mfiter(2) + mfit(4))/
     >             (mfit(1)           + mfit(3))
         chk1er2 = (mfit(2)           + mfit(4))/
     >             (mfit(1)+mfiter(1) + mfit(3))
         chk1er3 = (mfit(2) + mfit(4) + mfiter(4))/
     >             (mfit(1) + mfit(3))
         chk1er4 = (mfit(2) + mfit(4))/
     >             (mfit(1) + mfit(3) + mfiter(3))
         mcsv(ifit,ikinx,iz,8) = chk1
         mcsver(ifit,ikinx,iz,8) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fsum'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsv(ifit,ikinx,iz,8),mcsver(ifit,ikinx,iz,8)
c get sum and diff ratios using mfitp
         chk1 = (mfitp(2) - mfitp(4))/
     >          (mfitp(1) - mfitp(3))
         chk1er1 = (mfitp(2)+mfitper(2) - mfitp(4))/
     >             (mfitp(1)           - mfitp(3))
         chk1er2 = (mfitp(2)           - mfitp(4))/
     >             (mfitp(1)+mfitper(1) - mfitp(3))
         chk1er3 = (mfitp(2) - mfitp(4) - mfitper(4))/
     >             (mfitp(1) - mfitp(3))
         chk1er4 = (mfitp(2) - mfitp(4))/
     >             (mfitp(1) - mfitp(3) - mfitper(3))
         mcsvp(ifit,ikinx,iz,7) = chk1
         mcsvper(ifit,ikinx,iz,7) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fdif p'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,7),mcsvper(ifit,ikinx,iz,7)
         chk1 = (mfitp(2) + mfitp(4))/
     >          (mfitp(1) + mfitp(3))
         chk1er1 = (mfitp(2)+mfitper(2) + mfitp(4))/
     >             (mfitp(1)           + mfitp(3))
         chk1er2 = (mfitp(2)           + mfitp(4))/
     >             (mfitp(1)+mfitper(1) + mfitp(3))
         chk1er3 = (mfitp(2) + mfitp(4) + mfitper(4))/
     >             (mfitp(1) + mfitp(3))
         chk1er4 = (mfitp(2) + mfitp(4))/
     >             (mfitp(1) + mfitp(3) + mfitper(3))
         mcsvp(ifit,ikinx,iz,8) = chk1
         mcsvper(ifit,ikinx,iz,8) = sqrt(
     >    (chk1 - chk1er1)**2 +
     >    (chk1 - chk1er2)**2 +
     >    (chk1 - chk1er3)**2 +
     >    (chk1 - chk1er4)**2)
         write(6,'(''fsum p'',3i3,2f8.3)') ifit,ikinx,iz,
     >    mcsvp(ifit,ikinx,iz,8),mcsvper(ifit,ikinx,iz,8)

c back to multi-param fit again using mfit
         call mfit_fcn(npar,grad,chi2,coef,ierflg,futil)
         write(6,'(''mfit F'',2i2,i3,f8.2/8f8.3)') 
     >     itfit,ikinx,iz,chi2,(mfit(it),mfitr(it),it=1,4)
         write(17,'(i2,i2,i3,11f7.3)') ifit,ikinx,iz,
     >     chi2,coef(1)/fu1,std(1)/fu1,
     >    coef(2) / (fd1 / fu1),
     >    std(2) / (fd1 / fu1),
     >    coef(3),std(3),
     >    mcsv(ifit,ikinx,iz,4),
     >    mcsver(ifit,ikinx,iz,4)
        endif
       enddo ! iz
       avdelu(ifit,ikinx) = avdelu(ifit,ikinx) / avdeluer(ifit,ikinx)
       avdeld(ifit,ikinx) = avdeld(ifit,ikinx) / avdelder(ifit,ikinx)
      enddo ! ikin
c if skipping all the csv stuff for test
 77   continue
c end of big loop over ifit
      enddo ! ifit

      goto 99

c big global fit to all of the data
      do ifit= 1,2
      rhofact = 1.0
      write(6,'(''using rho fact = '',i2,f5.1)') ifit,rhofact
      if(ifit.eq.2) rhofact = 0.
c xxx fitx this in future?
      if(ifit.eq.3) then
       do i=1,bnpt
        brho(i) = brhop(i)
       enddo
      endif
      nparam=19
      npar = 19
      p(1) = 0.2
      p(2) = 0.2
      p(3) = 0.2
      p(4) = 0.2
      p(5) = 0.
      p(6) = 0.
      p(7) = 0.
      p(8) = 0.
      p(9) = 0.0
      p(10) = 0.
      do i=11,19
       p(i) = 0.
      enddo
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo
      call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
      write(6,'(''initial chi2'',f8.2)') chi2/float(bnpt)
c find best norm. factor for each kin.
      do ikin=1,32
       chibest=10000000.
       fbest=0.
       do kk=1,32
        fnorm(kk)=1.0
       enddo
       do k=1,40
        fnorm(ikin) = 0.8 + 0.01 * k
        call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
        if(chi2.lt.chibest) then
         fbest=fnorm(ikin)
         chibest = chi2
        endif
        fnormsv(ikin)=fbest
       enddo
      enddo
      do ikin=1,32
       fnorm(ikin)=fnormsv(ikin)
      enddo
      do ikin=1,32,2
       write(6,'(''ikin fnorm='',i3,3f6.2,3f7.3)') ikin,
     >  (xkin(ikin) + xkin(ikin+1))/2.,
     >  (q2kin(ikin) + q2kin(ikin+1))/2.,
     >  (wkin(ikin) + wkin(ikin+1))/2.,
     >  fnorm(ikin),fnorm(ikin+1),
     >  fnorm(ikin)/fnorm(ikin+1)
      enddo

c set norm factors back to 1.0
      do ikin=1,36
       fnorm(ikin)=1.0
      enddo

      call bfit_fcn(npar,grad,chi2,p,ierflg,futil)
      write(6,'(''initial chi2'',f8.2)') chi2/float(bnpt)
      do ikin=0,36
       if(bdfk(ikin).gt.0) then
        if(ikin.gt.0) then
         write(6,'(''kin,df,chi2/df'',i3,3f6.2,f7.0,f7.2)') ikin,
     >   xkin(ikin),q2kin(ikin),wkin(ikin),
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
        else
         write(6,'(''kin,df,chi2/df'',i3,18x,f7.0,f7.2)') ikin,
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
        endif
       endif
      enddo
      write(6,'(''inital x,q2,w,z,pt,mx'')')
      do k=1,15
       write(6,'(6(i5,f4.1))') 
     >  int(bdfx(k)),bchi2x(k)/max(1.,bdfx(k)),
     >  int(bdfq2(k)),bchi2q2(k)/max(1.,bdfq2(k)),
     >  int(bdfw(k)),bchi2w(k)/max(1.,bdfw(k)),
     >  int(bdfz(k)),bchi2z(k)/max(1.,bdfz(k)),
     >  int(bdfpt(k)),bchi2pt(k)/max(1.,bdfpt(k)),
     >  int(bdfmx(k)),bchi2mx(k)/max(1.,bdfmx(k))
      enddo
      do k=1,6
       write(6,'(i2,i7,f6.2)') k,int(bdft(k)),
     >   bchi2t(k)/max(1.,bdft(k))
      enddo

! initial parameters
      call mnparm( 1,"P1 ",p(1), 0.001D0,.19D0,.21D0,ierflg)
      call mnparm( 2,"P2 ",p(2), 0.001D0,.19D0,.21D0,ierflg)
      call mnparm( 3,"P3 ",p(3), 0.001D0,zero,zero,ierflg)
      call mnparm( 4,"P4 ",p(4), 0.001D0,zero,zero,ierflg)
      call mnparm( 5,"P5 ",p(5), 0.001D0,zero,zero,ierflg)
      call mnparm( 6,"P6 ",p(6), 0.001D0,zero,zero,ierflg)
      call mnparm( 7,"P7 ",p(7), 0.001D0,zero,zero,ierflg)
      call mnparm( 8,"P8 ",p(8), 0.001D0,zero,zero,ierflg)
      call mnparm( 9,"P9 ",p(9), 0.001D0,zero,zero,ierflg)
      call mnparm( 10,"P10 ",p(10), 0.001D0,zero,zero,ierflg)
      call mnparm( 11,"P11 ",p(11), 0.001D0,zero,zero,ierflg)
      call mnparm( 12,"P12 ",p(12), 0.001D0,zero,zero,ierflg)
      call mnparm( 13,"P13 ",p(13), 0.001D0,zero,zero,ierflg)
      call mnparm( 14,"P14 ",p(14), 0.001D0,zero,zero,ierflg)
      call mnparm( 15,"P15 ",p(15), 0.001D0,zero,zero,ierflg)
      call mnparm( 16,"P16 ",p(16), 0.001D0,zero,zero,ierflg)
      call mnparm( 17,"P17 ",p(17), 0.001D0,zero,zero,ierflg)
      call mnparm( 18,"P18 ",p(18), 0.001D0,zero,zero,ierflg)
      call mnparm( 19,"P19 ",p(19), 0.001D0,zero,zero,ierflg)
      ARGLIS(1) = 5.
      CALL MNEXCM(bfit_fcn,'FIX', ARGLIS ,1,IERFLG)
      arglis(1)=0
      call mnexcm(bfit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
      
      do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         write(6,'(''bf '',i2,8f7.3)') j,coef(j),std(j)
      enddo
      call bfit_fcn(npar,grad,chi2,coef,ierflg,futil)
      write(6,'(''final chi2/d.f.'',f8.2)') chi2/float(bnpt)
      do ikin=1,36
       if(bdfk(ikin).gt.0) then
        write(6,'(''kin,df,chi2/df'',i3,3f6.2,f7.0,f7.2)') ikin,
     >   xkin(ikin),q2kin(ikin),wkin(ikin),
     >   bdfk(ikin),bchi2k(ikin)/bdfk(ikin)
       endif
      enddo
      write(6,'(''final x,q2,w,z,pt,mx'')')
      do k=1,15
       write(6,'(6(i5,f4.1))') 
     >  int(bdfx(k)),bchi2x(k)/max(1.,bdfx(k)),
     >  int(bdfq2(k)),bchi2q2(k)/max(1.,bdfq2(k)),
     >  int(bdfw(k)),bchi2w(k)/max(1.,bdfw(k)),
     >  int(bdfz(k)),bchi2z(k)/max(1.,bdfz(k)),
     >  int(bdfpt(k)),bchi2pt(k)/max(1.,bdfpt(k)),
     >  int(bdfmx(k)),bchi2mx(k)/max(1.,bdfmx(k))
      enddo
      do k=1,6
       write(6,'(i2,i7,f6.2)') k,int(bdft(k)),
     >   bchi2t(k)/max(1.,bdft(k))
      enddo

      do i=1,bnpt
       if(bkin(i).eq.0) then
        write(6,'(''incl'',5f7.3)') bxv(i),bq2v(i),byv(i),
     >   bsigv(i),byv(i)/bsigv(i)
       endif
      enddo

      enddo ! j fits
 99   continue

      close(unit=22)
      open(unit=22,file='ptmrho.top')
      write(22,'(''set device postscript'')')
      do ipm=1,2
       write(22,122) ipm
 122   format('set window x ',i1,' of 2'/
     >  'title bottom',1h','hepgen/SIMC',1h')
       do kk=1,25
        write(22,'(f7.2,i8)') 0.25*(kk-0.5),rhorat(kk,ipm)
       enddo
       write(22,'(''hist'')')
      enddo

c plot 3-parameter fit values
      close(unit=31)
      open(unit=31,file='ptmmult.top')
      write(31,'(''set device postscript'')')
      close(unit=32)
      open(unit=32,file='ptmpt2.top')
      write(32,'(''set device postscript'')')
      close(unit=33)
      open(unit=33,file='ptmcos.top')
      write(33,'(''set device postscript'')')
      close(unit=34)
      open(unit=34,file='ptmcos2.top')
      write(34,'(''set device postscript'')')
      ikin=0
      do iy=1,4
       do ix=1,4
        ikin = ikin + 1
        xmin=1.0 + 2.5*(ix-1)
        xmax=xmin+2.5
        ymax=9.0-2.0*(iy-1)
        ymin=ymax-2.0
        do j=1,4
         if(j.eq.1) yymin = 0.
         if(j.eq.1) yymax = 0.099
         if(j.eq.2) yymin = 0.
         if(j.eq.2) yymax = 0.499
         if(j.eq.3) yymin = -1.
         if(j.eq.3) yymax = 0.999
         if(j.eq.4) yymin = -1.
         if(j.eq.4) yymax = 0.999
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         if(iy.eq.4) write(30+j,'(1x,''set labels bottom on'')')
         if(iy.lt.4) write(30+j,'(1x,''set labels bottom off'')')
         write(30+j,131) xmin,xmax,ymin,ymax,yymin,
     >    yymax,yymin + 0.9*(yymax-yymin),
     >    (xkin(2*ikin-1) + xkin(2*ikin))/2.,
     >    (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
 131     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.79 y',2f8.3/
     >     1x,'set ticks size 0.03 ; set labels size 1.0'/
     >     1x,' set sym 9O size 0.2'/
     >     1x,'title 5. 0.6 size 1.5',1h','z',1h'/
     >     1x,'title 0.3 ',f8.3,' data size 1.0 ',1h',
     >     '(',f4.2,',',f3.1,')',1h')
         if(j.eq.1) write(31,132)
 132     format('title 0.1 4.6 angle=90 size=1.5',1h','zM',1h')
         if(j.eq.1) write(32,133)
 133     format('title 0.1 4.6 angle=90 size=1.5',1h','<P0t1>223',1h'/
     >          'case                          ',1h','  X X X X',1h')
         if(j.eq.3) write(33,134)
 134     format('title 0.5 4. angle=90 size=1.5',1h','A',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.4) write(34,135)
 135     format('title 0.5 4. angle=90 size=1.5',1h','B',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         write(30+j,'(''plot axes'')')
         do itt=1,4
          if(itt.eq.1) write(30+j,'(''set color red'')')
          if(itt.eq.2) write(30+j,'(''set color green'')')
          if(itt.eq.3) write(30+j,'(''set color blue'')')
          if(itt.eq.4) write(30+j,'(''set color white'')')
          nplt=0
          ifit=1
          do iz=1,20
           if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >        scsver(ikin,itt,iz,j,ifit).lt.0.2) then
            nplt = nplt + 1
            write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >       scsv(ikin,itt,iz,j,ifit),
     >       scsver(ikin,itt,iz,j,ifit)
 136        format(5f10.4)
           endif
          enddo
          if(nplt.gt.0) write(30+j,'(''plot'')')
         enddo
        enddo
       enddo
      enddo
c plot results from ptsidis runs only
c plot 3-parameter fit values
c "n" is for narrow delta cuts
c results look quite similar to wide cuts
      close(unit=31)
      open(unit=31,file='ptmmults.top')
c      open(unit=31,file='ptmmultsn.top')
      write(31,'(''set device postscript'')')
      close(unit=32)
      open(unit=32,file='ptmpt2s.top')
c      open(unit=32,file='ptmpt2sn.top')
      write(32,'(''set device postscript'')')
      close(unit=33)
      open(unit=33,file='ptmcoss.top')
c      open(unit=33,file='ptmcossn.top')
      write(33,'(''set device postscript'')')
      close(unit=34)
      open(unit=34,file='ptmcos2s.top')
c      open(unit=34,file='ptmcos2sn.top')
      write(34,'(''set device postscript'')')
      open(unit=35,file='ptmffs.top')
c      open(unit=35,file='ptmffn.top')
      write(35,'(''set device postscript'')')
      open(unit=36,file='ptmfus.top')
c      open(unit=36,file='ptmfun.top')
      write(36,'(''set device postscript'')')
      open(unit=37,file='ptmcsvs.top')
c      open(unit=37,file='ptmcsvn.top')
      write(37,'(''set device postscript'')')
      open(unit=38,file='ptmcsvsd.top')
c      open(unit=38,file='ptmcsvnd.top')
      write(38,'(''set device postscript'')')
      open(unit=39,file='ptmfups.top')
c      open(unit=39,file='ptmfupsn.top')
      write(39,'(''set device postscript'')')
      open(unit=40,file='ptmfuds.top')
c      open(unit=40,file='ptmfudsn.top')
      write(40,'(''set device postscript'')')
      open(unit=41,file='ptmdiff.top')
c      open(unit=41,file='ptmdiffn.top')
      write(41,'(''set device postscript'')')
      open(unit=42,file='ptmsums.top')
c      open(unit=42,file='ptmsumsn.top')
      write(42,'(''set device postscript'')')
      ikin=0
      do iy=1,1
       do ix=1,3
        ikin = ikin + 1
        xmin=1.0 + 3.5*(ix-1)
        xmax=xmin + 3.5
        ymax=9.0-4.5*(iy-1)
        ymin=ymax-4.5
        do j=1,12
         if(j.eq.1) yymin = 0.
         if(j.eq.1) yymax = 0.15
         if(j.eq.2) yymin = 0.
         if(j.eq.2) yymax = 0.5
         if(j.eq.3) yymin = -1.
         if(j.eq.3) yymax = 1.
         if(j.eq.4) yymin = -1.
         if(j.eq.4) yymax = 1.0
         if(j.eq.5) yymin = 0.5
         if(j.eq.5) yymax = 1.5
         if(j.eq.6) yymin = 0. 
         if(j.eq.6) yymax = 1.0
         if(j.eq.7) yymin = -1.0
         if(j.eq.7) yymax = 1.0
         if(j.eq.8) yymin = -1.0
         if(j.eq.8) yymax = 1.0
         if(j.eq.9) yymin = 0. 
         if(j.eq.9) yymax = 1.0
         if(j.eq.10) yymin = 0. 
         if(j.eq.10) yymax = 1.0
         if(j.eq.11) yymin = 0. 
         if(j.eq.11) yymax = 2.0
         if(j.eq.12) yymin = 0. 
         if(j.eq.12) yymax = 2.0
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
         q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
         w2 = am**2 + q2*(1/x-1)
         zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
         zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
         WRITE(30+J,231) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,x,q2
 231     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.89 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 0.6'/
     >     1x,'title 6.2 4.1 size 1.5',1h','z',1h'/
     >     1x,'title top size 1.5 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X',1h'/
     >     'title ',1h','VERY PRELIMINARY',1H')
         if(j.eq.1) write(31,232)
 232     format('title 0.35 6.8 angle=90 size=1.8',
     >     1h','z223 M001',1h'/'case ',1h',' X X  X X',1h')
         if(j.eq.2) write(32,233)
 233     format('title 0.5 6.8 angle=90 size=1.5',1h','<P0t1>223',1h'/
     >          'case                          ',1h','  X X X X',1h')
         if(j.eq.3) write(33,234)
 234     format('title 0.5 6.8 angle=90 size=1.5',1h','A',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.4) write(34,235)
 235     format('title 0.5 6.8 angle=90 size=1.5',1h','B',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.5) write(35,239)
 239     format('title 0.5 6.8 angle=90 size=1.5',1h','F/F0DSS1',1h'/
     >    ' case                                ',1h','   X   X',1h'/
     >    '0. 1. ; 1. 1. ; join dash')
         if(j.eq.6 .or. j.eq.9 .or. j.eq.10) write(36,237)
 237     format('title 0.5 6.8 angle=90 size=1.5',1h','F0u1/F0f1',1h'/
     >     'case                                ',1h',' X X  X X',1h')
         if(j.eq.7) write(37,238)
         if(j.eq.8) write(38,238)
 238     format('title 0.5 6.8 angle=90 size=1.5',1h','Du/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.11) write(41,288)
 288     format('title 0.5 6.8 angle=90 size=1.5',
     >            1h','(d2+3-d2-3)/(p2+3-p2-3)',1h'/
     >    'case ',1h','  X X  X X    X X  X X ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         if(j.eq.12) write(42,289)
 289     format('title 0.5 6.8 angle=90 size=1.5',
     >            1h','(d2+3+d2-3)/(p2+3+p2-3)',1h'/
     >    'case ',1h','  X X  X X    X X  X X ',1h'/
     >    '0. 0. ; 1. 0. ; join dash')
         write(30+j,'(''plot axes'')')
         wpcut=2.5
         write(30+j,236) zforwp225,yymin+(yymax-yymin)*0.06,
     >    zforwp225,yymin,
     >    zforwp225-0.01,yymin+(yymax-yymin)*0.09,wpcut
         wpcut=3.0
         write(30+j,236) zforwp23,yymin+(yymax-yymin)*0.06,
     >    zforwp23,yymin,
     >    zforwp23-0.01,yymin+(yymax-yymin)*0.09,wpcut
 236     format('arrow from ',2f8.3,' data to ',
     >    2f8.3,' data size 1.5'/
     >    'title ',2f8.3,' data size 0.8 ',1h',f3.1,1h')
         if(j.le.4) then
          do itt=1,4
           if(itt.eq.1) write(30+j,'(''set color red'')')
           if(itt.eq.2) write(30+j,'(''set color green'')')
           if(itt.eq.3) write(30+j,'(''set color blue'')')
           if(itt.eq.4) write(30+j,'(''set color white'')')
           ifit=1
           nplt=0
           do iz=1,20
            if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >       (scsver(ikin,itt,iz,j,ifit).lt.0.007 .or. j.gt.1).and.
     >         scsver(ikin,itt,iz,j,ifit).lt.0.2) then
             nplt = nplt + 1
             fact=1.
             if(j.eq.1) fact = 0.05*(iz-0.5) * 2. * 3.1415928
             write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >        fact*scsv(ikin,itt,iz,j,ifit),
     >        fact*scsver(ikin,itt,iz,j,ifit)
            endif
           enddo
           if(nplt.gt.0) then
             write(30+j,'(''plot'')')
             write(30+j,'(''set symbol size 0.5 ; plot'')')
             write(30+j,'(''set symbol size 0.4 ; plot'')')
             write(30+j,'(''set symbol size 0.3 ; plot'')')
            write(30+j,'(''set symbol size 0.2 ; plot'')')
            write(30+j,'(''set symbol size 0.1 ; plot'')')
           endif
c          write(30+j,'(''set symbol 1O size 0.5'')')
           ifit=2
           nplt=0
           do iz=1,20
            if(scsver(ikin,itt,iz,j,ifit).ne.0. .and.
     >       (scsver(ikin,itt,iz,j,ifit).lt.0.007 .or. j.gt.1).and.
     >         scsver(ikin,itt,iz,j,ifit).lt.0.2) then
            nplt = nplt + 1
            write(30+j,136) 0.05*(iz-0.5)+0.005*(itt-2.5),
     >       scsv(ikin,itt,iz,j,ifit),
     >       scsver(ikin,itt,iz,j,ifit)/100000.
            endif
           enddo
           if(nplt.gt.0) then
             write(30+j,'(''plot'')')
           endif
           if(j.eq.1) then
            it=itt
            if(itt.gt.2) it = itt+1
            do iz=4,18
             z = 0.05 * iz
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
c correct b ratio of actual pt slope to one assumed
c in the fit
             dsigdz = dsigdz * sqrt(0.16**2 + z**2 * 0.25**2) /
     >                        sqrt(0.20**2 + z**2 * 0.12**2)
c            write(30+j,'(2f8.4)') z, z**2 * dsigdz/2./3.1415928
            enddo
c           write(30+j,'(''join'')')
c using DSS, still for ptmmult plot
            do iz=4,18
             z = 0.05 * iz
             call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >       sighad,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,2,ipdf)
c             write(30+j,'(2f8.4)') z, z**2 * dsigdz/2./3.1415928
             write(30+j,'(2f8.4)') z, z**2 * dsigdz
            enddo
c           write(31,'(''set pattern .01 .03 .01 .03 ; join pattern'')')
            write(30+j,'(''join'')')
           endif

           if(j.eq.2 .and. itt.eq.4) then
            do iz=1,20
             z = 0.05 * iz
c            write(30+j,'(2f8.4)') z, 0.2 * (1. + z**2)
             write(30+j,'(2f8.4)') z,sqrt(0.16**2 + z**2 * 0.40**2)
            enddo
            write(30+j,'(''join'')')
           endif
           if(itt.eq.2) write(30+j,242) yymin +
     >      (yymax - yymin) * 0.17
 242       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' d P2+3',1h'/
     >      'case ',1h','   GX X',1h')
          if(itt.eq.3) write(30+j,243) yymin +
     >      (yymax - yymin) * 0.12
 243       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' p P2-3',1h'/
     >      'case ',1h','   GX X',1h')
           if(itt.eq.4) write(30+j,244) yymin +
     >     (yymax - yymin) * 0.07
 244       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' d P2-3',1h'/
     >      'case ',1h','   GX X',1h')
           if(itt.eq.1) write(30+j,241) yymin +
     >      (yymax - yymin) * 0.22
 241       format('title 0.32 ',f7.4,' data size 1.0',
     >      1h',' p P2+3',1h'/
     >      'case ',1h','   GX X',1h')
          enddo ! itt
         endif ! j.le.4

         if(j.gt.4) then
          jp = j - 4
          do ifit=1,3
           nplt=0
           do iz=1,20
            if(mcsver(ifit,ikin,iz,jp).ne.0) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsv(ifit,ikin,iz,jp)*fact,
     >        mcsver(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
          enddo ! ifit
c plot both results from p only and d only on j=10 plot
          if(j.eq.10) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsver(ifit,ikin,iz,5).ne.0) then
              nplt = nplt + 1
              fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsv(ifit,ikin,iz,5)*fact,
     >         mcsver(ifit,ikin,iz,5)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.10
          if(jp.eq.1) then
           do iz=4,19
            z = 0.05*iz
            fact = 1.3 - 0.5 * (1-z)**4 / 0.7**4 -
     >                   0.2 * z**4 
            write(30+j,'(2f10.3)') z,fact
           enddo
           write(30+j,'(''join'')')
          endif
          if(jp.eq.2 .or. jp.eq.5 .or. jp.eq.6) then
c ratios from dss
           ifit=1
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             fact = fd1sv(ifit,ikin,iz)/fu1sv(ifit,ikin,iz)
             write(30+j,'(2f10.3)') z,fact
            endif
           enddo
           write(30+j,'(''join'')')
c ratios from geiger
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             call rgeiger(z,rnew)
             write(30+j,'(2f10.3)') z,rnew
            endif
           enddo
           write(30+j,'(''join dash'')')
          endif
          if(jp.eq.3.or.jp.eq.4) then
           x4 = x
           qgev = sqrt(q2)
           u = ctq5pdf(1,x4,qgev)
           d = ctq5pdf(2,x4,qgev)
! MRST fit. Note they have kappa = -0.2, but definition
! of delu is oppositie to mine
           f = 0.8 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
           f = -0.65 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
 251       format('0. ',f10.4,' ; 1. ',f10.4,' ; join dash')
          endif
          if(jp.eq.7.or.jp.eq.8) then
            do iz=3,18
             z = 0.05 * iz
             pt=0.
             phicm=0.
             call simcmodel(x,q2,z,pt,phicm,mmpi2,1 ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,2 ,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,4 ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,5 ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
             if(jp.eq.7) write(30+j,'(2f10.4)') z,
     >        (sighad2 - sighad4) / (sighad1 - sighad3)
             if(jp.eq.8) write(30+j,'(2f10.4)') z,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
             v1 = (d - db) / (u - ub)
             write(6,'(''dbgd'',i3,2f6.2,8f7.2)') ikin,x,q2,
     >        sighad1,sighad2,
     >        sighad3,sighad4,
     >        (sighad2 - sighad4) / (sighad1 - sighad3),
     >        0.6 * (1 - v1/4.)/(1+v1/4.),v1,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
            enddo
            write(30+j,'(''join'')')
          endif
         endif ! j.gt.4
        enddo ! j
       enddo ! ix
      enddo ! iy

c for ptSIDIS runs ony, get four FF and plot
      CALL MNINIT(861,862,863)
      do ikin=1,3
       do iz=6,16
       do ifit=1,3
        z = 0.05 * (iz-0.5)
        q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
        x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
        it = 1
        pt = 0.
        phicm = 0.
        ipdf = 1
        if(ifit.eq.3) ipdf=2
        call simcmodel(x,q2,z,pt,phicm,mmpi2,it,
     >   sighad,muf,mubf,mdf,mdbf,u1,d1,
     >   msf,msbf,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,1,ipdf)
        call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
        do it=1,4
         mv(it) = scsv(ikin,it,iz,1,ifit) * 2. * pi
         mver(it) = scsver(ikin,it,iz,1,ifit) * 2. * pi
        enddo
        write(6,'(''fffit'',3i3,8f7.3)') ikin,iz,ifit,
     >   (mv(it),mver(it),it=1,4)
c with four FF, errors too big
c        nparam=4
c        npar = 4
c chang to fit with two fav. one unfavored
        nparam=3
        npar = 3
        p(1) = fu1 ! fav u->pi+
        p(2) = fub ! unfav u->pi+
c        p(3) = fd1
c        p(4) = fdb
        p(3) = fdb
        call mnparm( 1,"P1 ",p(1), 0.0001D0,zero,zero,ierflg)
        call mnparm( 2,"P2 ",p(2), 0.0001D0,zero,zero,ierflg)
        call mnparm( 3,"P3 ",p(3), 0.0001D0,zero,zero,ierflg)
c        call mnparm( 4,"P4 ",p(4), 0.0001D0,zero,zero,ierflg)
        arglis(1)=0
        ncallff=0
        call mnexcm(fffit_fcn,'MIGRAD',arglis,0,nerror_migrad,0)
        do j=1,nparam
         call mnpout(j,pname(j),coef(j),std(j),zero,zero,ierflg)
         ffsv  (ifit,ikin,iz,j)= coef(j)
         ffsver(ifit,ikin,iz,j)= std(j)
         write(6,'(''fffit'',i2,i2,i3,i2,3f8.3)') ikin,iz,ifit,
     >    j,p(j),coef(j),std(j)
        enddo ! j
        enddo ! ifit
       enddo ! iz
      enddo ! ikin

c plot 3-paam FF results from pt-SIDIS runs
      close(unit=31)
      open(unit=31,file='ptmff3.top')
      write(31,'(''set device postscript'')')
      ikin=0
      do iy=1,1
       do ix=1,3
        ikin = ikin + 1
        xmin=1.0 + 3.5*(ix-1)
        xmax=xmin + 3.5
        ymax=9.0-4.5*(iy-1)
        ymin=ymax-4.5
        yymin=0.
        yymax=0.6
        x = (xkin(2*ikin-1) + xkin(2*ikin))/2.
        q2 = (q2kin(2*ikin-1) + q2kin(2*ikin))/2.
        w2 = am**2 + q2*(1/x-1)
        zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
        zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
        if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
        if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
        WRITE(31,2231) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,x,q2
 2231     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
c     >     1x,'set scale x log'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.79 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 1.0'/
     >     1x,'title 6.2 4.1 size 1.5',1h','z',1h'/
     >     1x,'title top size 1.5 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X',1h'/
     >     'title ',1h','VERY PRELIMINARY',1H')
        write(31,2232)
 2232   format('title 0.35 6.8 angle=90 size=1.8',
     >     1h','D(z)',1h')
        do j=1,3
         if(j.eq.2) write(31,'(''set color blue'')')
         if(j.eq.3) write(31,'(''set color red'')')
         ifit=1
         do iz=1,20
          if(ffsver(ifit,ikin,iz,j).ne.0. .and.
     >       ffsver(ifit,ikin,iz,j).lt.0.1) then
           nplt = nplt + 1
           fact=1.
c           fact = 0.05*(iz-0.5) * 2. * 3.1415928
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*ffsv  (ifit,ikin,iz,j),
     >        fact*ffsver(ifit,ikin,iz,j),
     >        fact*ffsv  (3,ikin,iz,j),
     >        fact*ffsver(3,ikin,iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''plot'')')
          write(31,'(''set symbol size 0.8 ; plot'')')
          write(31,'(''set symbol size 0.6 ; plot'')')
          write(31,'(''set symbol size 0.4 ; plot'')')
          write(31,'(''set symbol size 0.2 ; plot'')')
         endif
         ifit=2
         do iz=1,20
          if(ffsver(ifit,ikin,iz,j).ne.0. .and.
     >       ffsver(ifit,ikin,iz,j).lt.0.1) then
           nplt = nplt + 1
           fact=1.
c           fact = 0.05*(iz-0.5) * 2. * 3.1415928
           write(31,136) 0.05*(iz-0.5) + 0.005*(j-2),
     >        fact*ffsv  (ifit,ikin,iz,j),
     >        0.0 * fact*ffsver(ifit,ikin,iz,j)
          endif
         enddo
         if(nplt.gt.0) then
          write(31,'(''set symbol size 1.0 ; plot'')')
         endif
         do iz=6,19
          z = 0.05 * iz
          call fDSS (1,1,0, Z, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fdb
         enddo
         write(31,'(''join'')')
c plot fdss again using zprime and NLO
         do iz=6,19
          z = 0.05 * iz
          xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
          zp = (z / 2.) * (xp / x) *(1. + sqrt(1 - 4 * x**2 * am**2 *  
     >       (0.02 + 0.1) / z**2 / q2**2))
          call fDSS (1,1,1, Zp, Q2, 
     >      fU1, fUB, fD1, fDB, fS1, fSB, fC1, fB1, fGL1)
          if(j.eq.1) write(31,'(2f8.4)') z, fu1
          if(j.eq.2) write(31,'(2f8.4)') z, fd1
          if(j.eq.3) write(31,'(2f8.4)') z, fdb
         enddo
         write(31,'(''join dotdash'')')
        enddo ! j
       enddo ! ix
      enddo ! iy

c plot results from runs with both d and p (8 of them)
c using both intercept and averaged values
      close(unit=35)
      close(unit=36)
      close(unit=37)
      close(unit=38)
      close(unit=39)
      close(unit=40)
      close(unit=41)
      close(unit=42)
      close(unit=43)
      close(unit=44)
      close(unit=45)
      open(unit=35,file='ptmffp.top')
      write(35,'(''set device postscript'')')
      open(unit=36,file='ptmfup.top')
      write(36,'(''set device postscript'')')
      open(unit=37,file='ptmcsvp.top')
      write(37,'(''set device postscript'')')
      open(unit=38,file='ptmcsvpd.top')
      write(38,'(''set device postscript'')')
      open(unit=39,file='ptmfupp.top')
      write(39,'(''set device postscript'')')
      open(unit=40,file='ptmfudp.top')
      write(40,'(''set device postscript'')')
      open(unit=41,file='ptmdiffp.top')
      write(41,'(''set device postscript'')')
      open(unit=42,file='ptmsump.top')
      write(42,'(''set device postscript'')')
      open(unit=43,file='ptmdelup.top')
      write(43,'(''set device postscript'')')
      open(unit=44,file='ptmdeldp.top')
      write(44,'(''set device postscript'')')
      open(unit=45,file='ptmdoverp.top')
      write(45,'(''set device postscript'')')
      ikin=0
      do iy=1,2
       do ix=1,4
        ikin = ikin + 1
        ikinx = ikin
        if(ikin.eq.8) ikinx=15
        xmin=0.95 + 3.0*(ix-1)
        xmax=xmin + 3.0
        ymax=9.0-3.3*(iy-1)
        ymin=ymax-3.3
        do j=5,15
         if(j.eq.6) yymin = 0. 
         if(j.eq.6) yymax = 1.0
         if(j.eq.7) yymin = -0.9
         if(j.eq.7) yymax = 0.9
         if(j.eq.8) yymin = -0.9
         if(j.eq.8) yymax = 0.9
         if(j.eq.9) yymin = 0. 
         if(j.eq.9) yymax = 1.0
         if(j.eq.10) yymin = 0. 
         if(j.eq.10) yymax = 1.0
         if(j.eq.11) yymin = 0.0
         if(j.eq.11) yymax = 1.0
         if(j.eq.12) yymin = 0.5 
         if(j.eq.12) yymax = 1.5
         if(j.eq.13) yymin = -0.9
         if(j.eq.13) yymax =  0.9
         if(j.eq.14) yymin = -0.9
         if(j.eq.14) yymax =  0.9
         if(j.eq.15) yymin =  1.0
         if(j.eq.15) yymax =  1.99
         if(ix.eq.1) write(30+j,'(1x,''set labels left on'')')
         if(ix.ge.2) write(30+j,'(1x,''set labels left off'')')
         if(iy.eq.1) write(30+j,'(1x,''set labels bottom off'')')
         if(iy.eq.2) write(30+j,'(1x,''set labels bottom on'')')
         x = (xkin(2*ikinx-1) + xkin(2*ikinx))/2.
         q2 = (q2kin(2*ikinx-1) + q2kin(2*ikinx))/2.
         w2 = am**2 + q2*(1/x-1)
         zforwp225 = 1. - (2.5 - am**2) / (w2 - am**2)
         zforwp23 = 1. - (3.0 - am**2) / (w2 - am**2)
         WRITE(30+J,431) XMIN,XMAX,YMIN,YMAX,yymin,
     >     yymax,yymin + 0.87*(yymax-yymin),x,q2,sqrt(w2)
 431     format(1x,'set window x ',2f6.2,' y ',2f6.2/
     >     1x,'set bar size 0.0 ; set order x y dy'/
     >     1x,'set intensity 4'/
     >     1x,'set color white'/
     >     1x,'set limits x 0.25 0.89 y',2f8.3/
     >     1x,'set ticks size 0.04 ; set labels size 1.20'/
     >     1x,' set sym 9O size 0.6'/
     >     1x,'title 0.3 ',f8.3,' data size 1.1 ',1h',
     >     'x=',f4.2,' Q223=',f3.1,' GeV223 W=',f3.1,' GeV',1h'/
     >     1x,'case               ',1h',
     >     '        X X        X X ',1h')
         if(j.eq.5.or.j.ge.8) write(30+j,433)
 433     format('0. 1. ; 1. 1. ; join dash')

         if(ix.eq.1 .and. iy.eq.1) then
          write(30+j,432)
 432      format('title 6.6 2.0 size 2.0',1h','z',1h')
          if(j.eq.5) write(35,439)
 439      format('title 0.2 4.8 angle=90 size=2.0',1h','F/F0DSS1',1h'/
     >     ' case                                ',1h','   X   X',1h')
         endif
         if(j.eq.6 .or. j.eq.9 .or. j.eq.10) write(36,437)
 437     format('title 0.2 4.8 angle=90 size=1.5',1h','F0u1/F0f1',1h'/
     >     'case                                ',1h',' X X  X X',1h')
         if(j.eq.7) write(37,438)
         if(j.eq.8) write(38,438)
 438     format('title 0.2 4.8 angle=90 size=2.0',1h','Du/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h')
         if(j.eq.11) write(41,488)
 488     format('title 0.2 4.8 angle=90 size=2.0',
     >            1h','(d2+3 - d2-3) / (p2+3 - p2-3)',1h'/
     >    'case ',1h','  X X    X X      X X    X X ',1h')
         if(j.eq.12) write(42,489)
 489     format('title 0.2 4.8 angle=90 size=2.0',
     >            1h','(d2+3 + d2-3) / (p2+3 + p2-3)',1h'/
     >    'case ',1h','  X X    X X      X X    X X ',1h')
         if(j.eq.13) write(30+j,438)
         if(j.eq.14) write(30+j,458)
 458     format('title 0.2 4.8 angle=90 size=2.0',1h','Dd/(u+d)',1h'/
     >    'case                                 ',1h','G       ',1h')
         if(j.eq.15) write(30+j,459)
 459     format('title 0.2 4.8 angle=90 size=2.0',1h',
     >    'inclusive  d/p',1h')
         write(30+j,'(''plot axes'')')
         wpcut=2.5
         write(30+j,236) zforwp225,yymin+(yymax-yymin)*0.06,
     >    zforwp225,yymin,
     >    zforwp225-0.01,yymin+(yymax-yymin)*0.09,wpcut
         wpcut=3.0
         write(30+j,236) zforwp23,yymin+(yymax-yymin)*0.06,
     >    zforwp23,yymin,
     >    zforwp23-0.01,yymin+(yymax-yymin)*0.09,wpcut

         if(j.gt.4) then
          jp = j - 4
          do ifit=1,3
           nplt=0
           do iz=1,20
            if(mcsvper(ifit,ikin,iz,jp).ne.0 .and.
     >         mcsvper(ifit,ikin,iz,jp).lt.0.3) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             if(jp.eq.3) fact = 1. / uplusd(ifit,ikin)
             write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsvp(ifit,ikin,iz,jp)*fact,
     >        mcsvper(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
          enddo ! ifit
c add pt, cosphi fit results for pt-sidsis runs
          if(ikin.le.3 .and. j.ne.15) then
           write(30+j,'(''set color blue'')')
           ifit=1
           do iz=1,20
            if(mcsver(ifit,ikin,iz,jp).ne.0) then
             nplt = nplt + 1
             fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
             if(jp.eq.3) fact = 1. / uplusd(ifit,ikin)
             write(30+j,136) 0.01+0.05*(iz-0.5)+0.005*(ifit-1.5),
     >        mcsv(ifit,ikin,iz,jp)*fact,
     >        mcsver(ifit,ikin,iz,jp)*fact
            endif
           enddo
           if(nplt.gt.0) then
            write(30+j,'(''plot'')')
           endif
           write(30+j,'(''set color white'')')
          endif
c plot both results from p only and d only on j=10 plot
          if(j.eq.10) then
           write(30+j,'(''set color blue'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsvper(ifit,ikin,iz,5).ne.0) then
              nplt = nplt + 1
              fact=1.
             if(jp.eq.1) fact = 1. / fu1sv(ifit,ikin,iz)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsvp(ifit,ikin,iz,5)*fact,
     >         mcsvper(ifit,ikin,iz,5)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.10
c plot both results from geiger and dss if j=8
          if(j.eq.8) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsvper(ifit,ikin,iz,12).ne.0) then
              nplt = nplt + 1
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsvp(ifit,ikin,iz,12),
     >         mcsvper(ifit,ikin,iz,12)
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.8
c plot delu, deld from 5-param fits
          if(j.eq.13 .or. j.eq.14) then
           write(30+j,'(''set color green'')')
           do ifit=1,3
            nplt=0
            do iz=1,20
             if(mcsv5per(ifit,ikin,iz,j-10).ne.0) then
              nplt = nplt + 1
              fact = 1. / uplusd(ifit,ikin)
              write(30+j,136) 0.05*(iz-0.5)+0.005*(ifit-1.5),
     >         mcsv5p(ifit,ikin,iz,j-10)*fact,
     >         mcsv5per(ifit,ikin,iz,j-10)*fact
             endif
            enddo
            if(nplt.gt.0) then
             write(30+j,'(''plot'')')
            endif
           enddo ! ifit
           write(30+j,'(''set color white'')')
          endif ! j.eq.13 or 14
c plot d/p with delu = deld = 0. if ctegq (1), jam (3)
c also ratio from scalers
          if(j.eq.15) then
           write(30+j,'(''set color blue'')')
           do ifit=1,3,2
            nplt=0
            do iz=1,20
             if(dopsv(ifit,ikin,iz).ne.0) then
              nplt = nplt + 1
              write(30+j,136) 0.05*(iz-0.5),dopsv(ifit,ikin,iz)
             endif
            enddo
            if(nplt.gt.0.and.ifit.eq.1) write(30+j,'(''join'')')
            if(nplt.gt.0.and.ifit.eq.3) write(30+j,'(''join dash'')')
           enddo ! ifit
           write(30+j,'(''set color white'')')
           write(30+j,623) 2.*rrsv(ikin),2.*rrsv(ikin)
 623       format('0. ',f8.3, ' ; 1. ',f8.3, ' ; join')
           write(30+j,'(''set color green'')')
c checked in bcm.f that this is almost identical to calling
c f2_NMC_NEW with p and d targets 
          call FNP_NMC(X,Q2,ratt)
           write(30+j,623) 1.+ratt, 1.+ratt
           write(30+j,'(''set color white'')')
          endif ! j.eq.15
          if(jp.eq.1) then
           do iz=4,19
            z = 0.05*iz
            fact = 1.3 - 0.5 * (1-z)**4 / 0.7**4 -
     >                   0.2 * z**4 
            write(30+j,'(2f10.3)') z,fact
           enddo
           write(30+j,'(''join'')')
          endif
          if(jp.eq.2 .or. jp.eq.5 .or. jp.eq.6) then
c ratios from dss
           ifit=1
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             fact = fd1sv(ifit,ikin,iz)/fu1sv(ifit,ikin,iz)
             write(30+j,'(2f10.3)') z,fact
            endif
           enddo
           write(30+j,'(''join'')')
c ratios from geiger
           do iz=4,19
            z = 0.05*(iz-0.5)
            if(fu1sv(ifit,ikin,iz).gt.0.) then
             call rgeiger(z,rnew)
             write(30+j,'(2f10.3)') z,rnew
            endif
           enddo
           write(30+j,'(''join dash'')')
          endif
          if(jp.eq.3.or.jp.eq.4.or.jp.eq.9.or.jp.eq.10) then
           x4 = x
           qgev = sqrt(q2)
           u = ctq5pdf(1,x4,qgev)
           d = ctq5pdf(2,x4,qgev)
! MRST fit. Note they have kappa = -0.2, but definition
! of delu is oppositie to mine
           f = 0.8 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
           f = -0.65 * (1 - x)**4 / sqrt(x) * (x-0.0909) / 
     >      uplusd(1,ikin)
           write(30+j,251) f,f
          endif
          if(jp.eq.7.or.jp.eq.8) then
            do iz=3,18
             z = 0.05 * iz
             pt=0.
             phicm=0.
             iff = 2
             call simcmodel(x,q2,z,pt,phicm,mmpi2,1 ,
     >       sighad1,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,2 ,
     >       sighad2,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,4 ,
     >       sighad3,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             call simcmodel(x,q2,z,pt,phicm,mmpi2,5 ,
     >       sighad4,u,ub,d,db,u1,d1,
     >       s,sb,ff,fu,fs,dsigdz,fff,ffu,zpm,rfu,iff,ipdf)
             if(jp.eq.7) write(30+j,'(2f10.4)') z,
     >        (sighad2 - sighad4) / (sighad1 - sighad3)
             if(jp.eq.8) write(30+j,'(2f10.4)') z,
     >        (sighad2 + sighad4) / (sighad1 + sighad3)
            enddo
            write(30+j,'(''join'')')
          endif
         endif ! j.gt.4
        enddo ! j
       enddo ! ix
      enddo ! iy
c end of 2 by 4 plots

      do ikin=1,3
       do ifit=1,3
        write(6,'(''avcsv'',2i3,4f8.3)') ikin,ifit,
     >  (avcsv(ikin,ifit,k),avcsver(ikin,ifit,k),k=1,2)
       enddo
      enddo ! ikin

      do ikin=1,8
       do ifit=1,3
        do iz=1,20
         if(mcsv5per(ifit,ikin,iz,1).ne.0.) then
           write(6,'(''m5'',3i3,10f7.2)') ikin,ifit,iz,
     >      (mcsv5p(ifit,ikin,iz,j),
     >       mcsv5per(ifit,ikin,iz,j),j=1,5)
         endif
        enddo
       enddo
      enddo ! ikin

 999  stop
      end

c modified Geiger ratio to give best agreement with our data
      subroutine rgeiger(z,rnew)
      implicit none
      real*8 z,rnew

      rnew = (1.0 -z)**0.083583 / (1.0 +z)**1.9838
c     > * (1.33 - 0.6 * z)
      return
      end


      subroutine getrho(x,q2,z,pt,phi,rplus,rminus,
     >  rplusd,rminusd)
      implicit none
      integer i,ix,iq,iz,ipt,ipm,ncall/0/
      real*8 abst(2),rat(2),asv(10,10,10,10,2)
      real*8 dp(2),dpsv(10,10,10,10,2),rplusd,rminusd
      real*8 ratsv(10,10,10,10,2),x,q2,z,pt,phi,rplus,rminus
      real*8 e0,ep,theta,nu,sin2,am/0.938/
      logical first/.true./

      if(first) then
       open(unit=11,file='rhoratio.txt')
       do i=1,9200
        read(11,'(4i3,6f8.3)') ix,iq,iz,ipt,
     >      abst(1),abst(2),rat(1),rat(2),dp(1),dp(2)
        do ipm=1,2
         asv(ix,iq,iz,ipt,ipm)=abst(ipm)
         ratsv(ix,iq,iz,ipt,ipm)=rat(ipm)
c protect against zero valuess
         if(dp(ipm).lt.0.001) dp(ipm)=1.0
         dpsv(ix,iq,iz,ipt,ipm)=dp(ipm)
        enddo
       enddo
       first=.false.
      endif

      e0 = 10.6
      nu = q2 / 2. / am / x
      ep = e0 - nu
      sin2 = q2 / 4. / e0 / ep
      theta = 57.3 * 2. * asin(sqrt(sin2))
      ix = min(10,max(1,int((ep - 3) / 4. * 10) + 1))
      iq = min(10,max(1,int((theta - 12.) / 10.0  * 10) + 1))
      iz = min(10,max(1,int((z-0.25)/0.6 * 10)+1))
      ipt = min(10,max(1,int(pt / 0.7 * 10.) + 1))

      rplus = ratsv(ix,iq,iz,ipt,1) * 
     >  (1. + asv(ix,iz,iz,ipt,1) * cos(phi)) 
      rminus = ratsv(ix,iq,iz,ipt,2) * 
     >  (1. + asv(ix,iz,iz,ipt,2) * cos(phi)) 
      rplusd = 0.
      rminusd = 0.
      if(dpsv(ix,iz,iz,ipt,1).gt.0 .and.
     >   dpsv(ix,iz,iz,ipt,2).gt.0.) then
       rplusd = rplus / dpsv(ix,iz,iz,ipt,1)
       rminusd = rminus / dpsv(ix,iz,iz,ipt,2)
      endif
c for checking
      ncall = ncall + 1
      if((ncall/100)*100.eq.ncall) then
       write(6,'(''getrho'',4i3,10f6.2)') ix,iq,iz,ipt,
     >  x,q2,ep,theta,z,pt,rplus,rminus,rplusd,rminusd
      endif

      return
      end

      subroutine mfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deld
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag

      fav = p(1)
      unf = p(2) * fav
      u = mu
      d = md
      ub = mub
      db = mdb
      s = ms
      sb = msb
c      delu = p(3) * (u + d)
c      deld = -1.* delu 
      delu = -p(3)
c      deld = -p(4)
      deld = -1./4. *delu

      sqp = 4 * (u + ub) + (d + db) + (s + sb)
c  for test: this makes delu about 50% bigger
c      sqn = 4 * (d + deld + db) + 
c     >               (u + delu + ub) + (s + sb)
c took delu out of denominator
      sqn = 4 * (d + db) + (u + ub) + (s + sb)

c change u and d to preseve incl. p and n cross sections
c took out: results in HUGE values of delu and
c more than 10% change to u (relative)
c      u = mu - delu/5.
c      d = md + 4.*delu/5.
      sqpchk = 4 * (u + ub) + (d + db) + (s + sb)
      sqnchk = 4 * (d + deld + db) + (u + delu +ub) + (s + sb)
c      If(abs(sqp - sqpchk).gt.0.001*mu .or.
c     >   abs(sqn - sqnchk).gt.0.001*mu) write(6,
c     >   '(''error sqp,sqn'',4f10.4)') sqp,sqpchk,sqn,sqnchk

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp = (4 * d + ub + 4.*deld) * fav + 
     >       (4 * db + u  + delu + s + sb) * unf 

      mnm = (4 * db + u + delu ) * fav + 
     >      (4 * d  + 4.*deld + ub  + s + sb) * unf 

      mfitr(1) = mpp/sqp
      mfitr(2) = (mpp + mnp)/(sqp + sqn)
      mfitr(3) = mpm/sqp
      mfitr(4) = (mpm + mnm)/(sqp + sqn)

      chi2 = 0.
      do i=1,4
       if(mfitflag.eq.0) then
        chi2 = chi2 + (mfit(i) -  mfitr(i))**2/mfiter(1)**2
       else
        chi2 = chi2 + (mfitp(i) - mfitr(i))**2/mfitper(1)**2
       endif
      enddo

      fval = chi2

      return
      end

      subroutine mfit5_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deld,rr,sqp0,sqn0
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4)
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag

      fav = p(1)
      unf = p(2) * fav
      u = mu
      d = md * p(5)
      ub = mub
      db = mdb
      s = ms
      sb = msb
c      delu = p(3) * (u + d)
c      deld = -1.* delu 
      delu = -p(3)
      deld = -p(4)

c values with modified d/u
      sqp = 4 * (u + ub) + (d + db) + (s + sb)
      sqn = 4 * (d + deld + db) + 
     >               (u + delu + ub) + (s + sb)

c values with original values of d
      sqp0 = 4 * (u + ub) + (md + db) + (s + sb)
      sqn0 = 4 * (md + db) + (u + ub) + (s + sb)

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp = (4 * d + ub + 4.*deld) * fav + 
     >       (4 * db + u  + delu + s + sb) * unf 

      mnm = (4 * db + u + delu ) * fav + 
     >      (4 * d  + 4.*deld + ub  + s + sb) * unf 

c fit doesn't converge with this form
      mfitr(1) = mpp/sqp
      mfitr(2) = (mpp + mnp)/(sqp + sqn)
      mfitr(3) = mpm/sqp
      mfitr(4) = (mpm + mnm)/(sqp + sqn)

c try this
      mfitr(1) = mpp/sqp0
      mfitr(2) = (mpp + mnp)/(sqp0 + sqn0)
      mfitr(3) = mpm/sqp0
      mfitr(4) = (mpm + mnm)/(sqp0 + sqn0)

      chi2 = 0.
      do i=1,4
       if(mfitflag.eq.0) then
        chi2 = chi2 + (mfit(i) -  mfitr(i))**2/mfiter(1)**2
       else
        chi2 = chi2 + (mfitp(i) - mfitr(i))**2/mfitper(1)**2
       endif
      enddo

c added contraint: p/d inclusive must remain unchanged
      rr = (sqp + sqn) / (sqp0 + sqn0)
      chi2 = chi2 + (rr-1)**2/0.05**2

      fval = chi2

      write(6,'(''mm5p'',5f7.2,f7.3,f7.1)') (p(i),i=1,5),rr,chi2

      return
      end

      subroutine getcsv(delu, deluer, deld, delder,doverp,
     >  doverpp, doverpper,fact)
! get delu and deld from sum and diff ratios assuming
c d/u, sea quarks known
      implicit none

      integer i,j,it
      real*8 chi2, u, d, ub, db, s, sb, fav, unf, sqpchk,sqnchk
      real*8 mpp, mpm, mnp, mnm, sqp, sqn, delu, deluer, deld, delder
      real*8 mfit(4),mfiter(4),mu,md,mub,mdb,ms,msb,mfitr(4)
      real*8 mfitp(4),mfitper(4),doverp,doverpp, doverpper,sqdp
      real*8 dop(5),sqnp
      integer mfitflag
      common/mstuff/ mfit,mfiter,mfitr,mfitp,mfitper,
     > mu,md,mub,mdb,ms,msb,mfitflag
      real*8 v(4),dv(4),sqd,mnp0,mnm0,srat,drat,drat0,C,E
      real*8 dd(5),du(5),sratp,dratp,srat0,fact

      fav = 1.
      unf = 0.5
      u = mu
c fact modifies d quark 
      d = md * fact
      ub = mub
      db = mdb
      s = ms
      sb = msb

      delu = 0.
      deld = 0. 

      sqp = 4 * (u + ub) + (d + db) + (s + sb)
      sqn = 4 * (d + db) + (u + ub) + (s + sb)
      sqd = sqp + sqn
! nominal ratio of d / p inclusive
      doverp = (sqp + sqn ) / sqp

c get sidis multiplicities
      mpp = (4 * u + db) * fav + 
     >       (4 * ub + d  + s + sb) * unf 

      mpm = (4 * ub + d ) * fav + 
     >       (4 * u  + db  + s + sb) * unf 

      mnp0 = (4 * d + ub - 4.*deld) * fav + 
     >       (4 * db + u  - delu + s + sb) * unf 

      mnm0 = (4 * db  + u - delu ) * fav + 
     >      (4 * d  - 4.*deld + ub  + s + sb) * unf 

      do it=1,4
       v(it) = mfit(it)
       dv(it) = mfiter(it)
       if(mfitflag.ne.0) then
        v(it) = mfitp(it)
        dv(it) = mfitper(it)
       endif
      enddo

      do i=1,5
       if(i.eq.1) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1)          + v(3)          )
        drat = (v(2)          - v(4)          ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.2) then
        srat = (v(2) + dv(2)  + v(4)          ) /
     >         (v(1)          + v(3)          )
        drat = (v(2) + dv(2)  - v(4)          ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.3) then
        srat = (v(2)          + v(4) + dv(4)  ) /
     >         (v(1)          + v(3)          )
        drat = (v(2)          - v(4) - dv(4)  ) /
     >         (v(1)          - v(3)          )
       endif
       if(i.eq.4) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1) + dv(1)  + v(3)          )
        drat = (v(2)          - v(4)          ) /
     >         (v(1) + dv(1)  - v(3)          )
       endif
       if(i.eq.5) then
        srat = (v(2)          + v(4)          ) /
     >         (v(1)          + v(3) + dv(3)  )
        drat = (v(2)          - v(4)          ) /
     >         (v(1)          - v(3) - dv(3)  )
       endif
       srat0 = ((sqp / sqd) * (mpp + mnp0 + mpm + mnm0) / 
     >                    (mpp        + mpm))
       C = (srat0 - srat) * (sqd / sqp) * (mpp + mpm) / (fav + unf) 
       drat0 = ((sqp / sqd) * (mpp + mnp0 - mpm - mnm0) / 
     >                    (mpp        - mpm))
       E = (drat0 - drat) * (sqd / sqp) * (mpp - mpm) / (fav - unf)
       dd(i) = (C + E)/8.
       du(i) = (C - E)/2.

! ratio with delu and deld
       sqnp = 4 * (d  - dd(i) + db) + (u - du(i) + ub) + (s + sb)
       dop(i) = (sqp + sqnp) / sqp

c check
       delu = du(i)
       deld = dd(i)
       mnp = (4 * d + ub - 4.*deld) * fav + 
     >       (4 * db + u  - delu + s + sb) * unf 
       mnm = (4 * db + u - delu ) * fav + 
     >      (4 * d  - 4.*deld + ub  + s + sb) * unf 
       dratp = ((sqp / sqd) * (mpp + mnp - mpm - mnm) / 
     >                    (mpp        - mpm))
       sratp = ((sqp / sqd) * (mpp + mnp + mpm + mnm) / 
     >                    (mpp        + mpm))
       if(abs(srat/sratp - 1.).gt.0.001 .or.
     >    abs(drat/dratp - 1.).gt.0.001) write(6,
     >    '(''csverr'',8f8.3)') delu,deld,
     <   drat,drat0,dratp,srat,srat0,sratp

      enddo

      delu = du(1)
      deld = dd(1)
      deluer = sqrt((du(1) - du(2))**2 + 
     >              (du(1) - du(3))**2 +
     >              (du(1) - du(4))**2 +
     >              (du(1) - du(5))**2)
      delder = sqrt((dd(1) - dd(2))**2 + 
     >              (dd(1) - dd(3))**2 +
     >              (dd(1) - dd(4))**2 +
     >              (dd(1) - dd(5))**2)

      doverpp = dop(1)
      doverpper = sqrt((dop(1) - dop(2))**2 + 
     >                (dop(1) - dop(3))**2 +
     >                (dop(1) - dop(4))**2 +
     >                (dop(1) - dop(5))**2)


      return
      end

      subroutine sfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 u,d,rpp,rnp,rpm,rnm,b1,b2,ds, w,mmpi2,w2
      real*8 qu,qd,qs,nu,c0,c1,c2,c3,c4,dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 sum_sqp,sum_sqn,fu,fub,fd,fdb,fs,fsb
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/,pi/3.141593/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      integer bkin(90000),ncall
      real*8 rhofact,fnorm(36)
      real*8 bchi2x(15),bdfx(15)
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt,ncall
        integer ikinfit,itfit,npts
        real*8 zfit
        common/sstuff/ zfit,ikinfit,itfit,npts
        integer sn
        real*8 sphi(1000), spt(1000),sy(10000),syer(1000)
        common/sforplot/ sphi,spt,sy,syer,sn

        ncall = ncall + 1
        chi2 = 0.
        npts = 0
        sn = 0

      do i=1,bnpt
       if(bitv(i).eq.itfit .and. bkin(i).ge.ikinfit .and.
     >   bkin(i).le.ikinfit + 1.and.
     >   abs(bzv(i) - zfit).lt.0.025) then
        npts = npts + 1
        sn = sn + 1
        phi = bphiv(i)
        sphi(sn) = phi
c        if(bkin(i).eq.ikinfit) sphi(sn) = sphi(sn)+0.05
c        if(bkin(i).eq.ikinfit+1) sphi(sn) = sphi(sn)-0.05
        pt = bptv(i)
        spt(sn) = pt
        pt2 = pt**2
        z = bzv(i)
c using p1 / z
        b = 1./p(2)
ccc        sig = (p(1)/z  ) * 2.*sqrt(b/pi) * exp(-b *pt2) *
        sig = (p(1)/z  ) * b * exp(-b *pt2) / 2. / pi **
     >  (1. + p(3) * pt * cos(phi))
c     >  (1. + p(3) * pt * cos(phi)  + 
c     >        p(4) * pt2 * cos(2.*phi))
        bsigv(i) = sig
! how much rho subtraction
        y = byv(i)  - rhofact * brho(i)
        sy(sn) = y
        syer(sn) = byerv(i)
! special code to use brhop. Not currently in use, I hope
        if(rhofact.eq.10.) y = byv(i)  - brhop(i)
        if(ncall.eq.1 .and. ikinfit.eq.1 .and. zfit.gt.0.5
     >   .and. zfit.lt.0.55) then
         write(6,'(''sdbg'',f4.1,i3,3f8.3)') rhofact,sn,y,byerv(i),
     >    rhofact * brho(i)
        endif
        chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       endif
      enddo

      if((ncall/100)*100.eq.ncall .and. ncall.gt.0 .and.
     >  npts.gt.0)
     >  write(6,'(''ncall small, chi2'',i8,i6,2i3,f6.2,f8.3)')
     >  ncall,npts,ikinfit,itfit,zfit,chi2/float(npts)

      fval = chi2

      return
      end

      subroutine bfit_fcn(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 u,d,rpp,rnp,rpm,rnm,b1,b2,ds, w,mmpi2,w2
      real*8 qu,qd,qs,nu,c0,c1,c2,c3,c4,dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 sum_sqp,sum_sqn,fu,fub,fd,fdb,fs,fsb
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000),bmmpi2(90000)
      real*8 byv(90000),byerv(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bffub(90000),bffdb(90000),bffsb(90000)
      real*8 buv(90000),bubv(90000),brho(90000),brhop(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(0:56),bdfk(0:56)
      integer bkin(90000),ncall/0/
      real*8 rhofact,fnorm(36)
      real*8 bchi2x(15),bdfx(15)
      real*8 bchi2q2(15),bdfq2(15)
      real*8 bchi2w(15),bdfw(15)
      real*8 bchi2z(15),bdfz(15)
      real*8 bchi2f(15),bdff(15)
      real*8 bchi2pt(15),bchi2mx(15),bdfpt(15),bdfmx(15)
      real*8 bchi2t(6),bdft(6)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      common/bstuff/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,
     >  bchi2x,bdfx,bffub,bffdb,bffsb,
     >  bchi2q2,bdfq2,bmmpi2,brhop,
     >  bchi2w,bdfw,
     >  bchi2z,bdfz,
     >  bchi2f,bdff,
     >  bchi2pt,bdfpt,
     >  bchi2mx,bdfmx,
     >  bchi2t,bdft,fnorm,
     >  rhofact,bnpt

      ncall = ncall + 1
      chi2 = 0.
      bchi2k(0)=0.
      bdfk(0)=0.
      do i=1,56
       bchi2k(i)=0.
       bdfk(i)=0.
       if(i.le.15) bchi2x(i)=0.
       if(i.le.15) bdfx(i)=0.
       if(i.le.15) bchi2q2(i)=0.
       if(i.le.15) bdfq2(i)=0.
       if(i.le.15) bchi2w(i)=0.
       if(i.le.15) bdfw(i)=0.
       if(i.le.15) bchi2z(i)=0.
       if(i.le.15) bdfz(i)=0.
       if(i.le.15) bchi2f(i)=0.
       if(i.le.15) bdff(i)=0.
       if(i.le.15) bchi2pt(i)=0.
       if(i.le.15) bdfpt(i)=0.
       if(i.le.15) bchi2mx(i)=0.
       if(i.le.15) bdfmx(i)=0.
       if(i.le.6) bchi2t(i)=0.
       if(i.le.6) bdft(i)=0.
      enddo

      do i=1,bnpt
       it = bitv(i)
       u = buv(i) 
       d = bdv(i)
       s = bsv(i)
       ubar = bubv(i)
       dbar = bdbv(i)
       sbar = bsbv(i)
       delu = p(19) * (1 - x)**4 / sqrt(x) * (x-0.0909)
       deld = -1. * delu
       sum_sqp = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)
       sum_sqn = qu**2*(d + deld + dbar) + 
     >         qd**2*(u + delu + ubar) + qs**2*(s+sbar)
       sum_sq = sum_sqp
       if (it.eq.2 .or. it.eq.5) then
        sum_sq = sum_sqp + sum_sqn
       endif
       if(bkin(i).gt.0.) then
       phi = bphiv(i)
       pt = bptv(i)
       pt2 = pt**2
       q2 = bq2v(i)
       x = bxv(i)
       mmpi2 = bmmpi2(i)
       w = sqrt(0.938**2 + q2 * (1./x -1))
       w2 = w**2
       nu = q2 / 2. / 0.938 / x
       y = nu / 10.6
       xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
       z0 = bzv(i)
       zp = (z0 / 2.) * (xp / x) *(1. + 
     >     sqrt(1 - 4 * x**2 * am**2 *  
     >     (ampi**2 + pt2) / z0**2 / q2**2))
c       z = z0
c use zp instead of z
       z = zp
       z2 = z**2
c       if(i .lt. 100 .and. ncall.eq.1)
c     >   write(6,'(10f7.3)') xp/x,zp/z0,x,xp,z0,z,zp,q2,am,ampi

c start off with Binneweis frag. func. 
c        sv = log( log(Q2/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
c        N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
c        a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
c        a2 = 1.430 + 0.543*sv - 0.023*sv**2
c        N = n * p(5)
c        a1 = a1 * p(6)
c        a2 = a2 * p(7)
c	Ns = 4.250 - 3.147*sv + 0.755*sv**2
c	a1s = -0.770 -0.573*sv + 0.117*sv**2
c	a2s = 4.48 + 0.890*sv - 0.138*sv**2
cxx	zhadm = min(z, 0.75)
c        zhadm = min(z, 0.99)
c	D_sum = N * zhadm**a1 * (1.0-zhadm)**a2
c	D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
c        b1 = 0.083583 * p(8)
c        b2 = 1.983 * p(9)
c	R_D = (1.0-zhadm)**b1 / (1.0+zhadm)**b2
c	D_fav = D_sum/(1.0+R_D)
c	D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
c	D_s = D_sum_s/2.0
c        Dp = d_fav
c        Dm = D_unfav
c new way using DSS
       fu = bffu(i) * (1. + p(5)/w2 +
     >  p(6)/mmpi2 + p(7)/q2)
       fdb= bffdb(i) * (1. + p(5)/w2 +
     >  p(6)/mmpi2 + p(7)/q2)
       fd = bffd(i) * (1. + p(8)/w2 +
     >  p(9)/mmpi2 + p(10)/q2)
       fub= bffub(i) * (1. + p(8)/w2 +
     >  p(9)/mmpi2 + p(10)/q2)
       fs = bffs(i)
       fsb = bffsb(i)
       c1 = 1. / (z2*p(1) + p(3)) / 9.
       c2 = 1. / (z2*p(2) + p(4)) / 9.
       c3 = 1. / (z2*p(1) + p(4)) / 9.
       c4 = 1. / (z2*p(2) + p(3)) / 9.

c need to add strange sea!
       rpp = 4.*u*fu * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >          d*fd * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >       dbar*fdb* c1*exp(-pt2/(z2*p(1) + p(3)))+
     >    4.*ubar*fub* c2*exp(-pt2/(z2*p(2) + p(4)))
     >   + (s + sbar)*fs * c2*exp(-pt2/(z2*p(2) + p(4)))

       rpm = 4.*u*fub* c3*exp(-pt2/(z2*p(1) + p(4)))+
     >          d*fdb* c4*exp(-pt2/(z2*p(2) + p(3)))+
     >       dbar*fd * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >    4.*ubar*fu * c4*exp(-pt2/(z2*p(2) + p(3)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))
       rnp = 4.*(d+deld)*fu * c4*exp(-pt2/(z2*p(2) + p(3)))+
     >          (u+delu)*fd * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >       ubar*fdb* c4*exp(-pt2/(z2*p(2) + p(3)))+
     >    4.*dbar*fub* c3*exp(-pt2/(z2*p(1) + p(4)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))
       rnm = 4.*(d+deld)*fub* c2*exp(-pt2/(z2*p(2) + p(4)))+
     >          (u+delu)*fdb* c1*exp(-pt2/(z2*p(1) + p(3)))+
     >       ubar*fd * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >    4.*dbar*fu * c1*exp(-pt2/(z2*p(1) + p(3)))
     >   + (s + sbar)*Ds * c2*exp(-pt2/(z2*p(2) + p(4)))

       if(it.eq.1) sig = rpp * 
     >   (1. + p(11)*sqrt(pt2/q2)*cos(phi) +
     >       p(12)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.2) sig = (rpp + rnp) *
     >   (1. + p(13)*sqrt(pt2/q2)*cos(phi) +
     >       p(14)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.4) sig = rpm *
     >   (1. + p(15)*sqrt(pt2/q2)*cos(phi) +
     >       p(16)*(pt2/q2)*cos(2.*phi)) 
       if(it.eq.5) sig = (rpm + rnm) * 
     >   (1. + p(17)*sqrt(pt2/q2)*cos(phi) +
     >       p(18)*(pt2/q2)*cos(2.*phi)) 

       sig = sig / sum_sq / 2. / 3.1415928

       bsigv(i) = sig

! how much rho subtraction
       y = byv(i)  - rhofact * brho(i)
c normalization factor by ikin
       y = y * fnorm(bkin(i))
c inclusive ratio here
       else
        sig = 1. + sum_sqn / sum_sqp
       endif
       bsigv(i) = sig
       chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       j=bkin(i)
       bchi2k(j) = bchi2k(j) + (sig - y)**2 / byerv(i)**2
       bdfk(j) = bdfk(j) + 1.
       if(j.ne.0) then
        j = min(15,max(1,int((x-0.2)/0.6*15)+1))
        bchi2x(j) = bchi2x(j) + (sig - y)**2 / byerv(i)**2
        bdfx(j) = bdfx(j) + 1.
        j = min(15,max(1,int((w-1.6)/1.8*15)+1))
        bchi2w(j) = bchi2w(j) + (sig - y)**2 / byerv(i)**2
        bdfw(j) = bdfw(j) + 1.
        j = min(15,max(1,int((q2-1.)/6.*15)+1))
        bchi2q2(j) = bchi2q2(j) + (sig - y)**2 / byerv(i)**2
        bdfq2(j) = bdfq2(j) + 1.
        j = min(15,max(1,int(phi/2./3.1415928*15)+1))
        bchi2f(j) = bchi2f(j) + (sig - y)**2 / byerv(i)**2
        bdff(j) = bdff(j) + 1.
        j = min(15,max(1,int((z-0.2)/0.75*15)+1))
        bchi2z(j) = bchi2z(j) + (sig - y)**2 / byerv(i)**2
        bdfz(j) = bdfz(j) + 1.
        j = min(15,max(1,int(pt/0.95*15)+1))
        bchi2pt(j) = bchi2pt(j) + (sig - y)**2 / byerv(i)**2
        bdfpt(j) = bdfpt(j) + 1.
        j = min(15,max(1,int((sqrt(mmpi2)-1.4)/1.4*15)+1))
        bchi2mx(j) = bchi2mx(j) + (sig - y)**2 / byerv(i)**2
        bdfmx(j) = bdfmx(j) + 1.
        j = it
        bchi2t(j) = bchi2t(j) + (sig - y)**2 / byerv(i)**2
        bdft(j) = bdft(j) + 1.
       endif
      enddo

      if((ncall/100)*100.eq.ncall .and. ncall.gt.0) 
     >  write(6,'(''ncall, chi2'',i8,f8.3)')
     >  ncall,chi2/float(bnpt)

      fval = chi2

      return
      end

      subroutine bfit_fcnOLD(npar,grad,fval,p,iflag,futil)
! Calculate Chisq for Minuit
      implicit none
      integer npar,iflag
      real*8 grad(*)
      real*8 p(*) ! vector of parameters
      real*8 fval  ! chisq
      real*8 futil ! auxially function
      external futil
      integer i,j,it
      real*8 chi2, y, phi, pt, pt2, b0, b, x, q2, z, z2, sig
      real*8 fm,p1p,p2p,p3p,p4p,u,d,rpp,rnp,rpm,rnm,rppp,rpmp,b1,b2
      real*8 qu,qd,qs,nu,er,c0,c1,c2,c3,c4,c1p,c2p,c3p,c4p, dp, dm
      real*8 sv,n,a1,a2,n2,a1s,a2s,d_fav,d_unfav,r_d,d_sum,d_sum_s
      real*8 lambda, Q2zero,sum_sq,ubar,dbar,s,sbar,zhadm,d_s,Ns
      real*8 deld,delu,xp,z0,zp,ampi/0.14/,am/0.938/
      parameter (qu=2./3.)
      parameter (qd=-1./3)
      parameter (qs=-1./3.)
      parameter (lambda=0.227)  !0.227 GeV for NLO
      parameter (Q2zero=2.0)    !Gev^2 for u,d,s,g
      integer bnpt,bitv(90000)
      real*8 bptv(90000),bzv(90000),bphiv(90000)
      real*8 bq2v(90000),bxv(90000)
      real*8 byv(90000),byerv(90000)
      real*8 buv(90000),bubv(90000),brho(90000)
      real*8 bffu(90000),bffd(90000),bffs(90000)
      real*8 bdv(90000),bdbv(90000),bchi2k(56),bdfk(56)
      real*8 bsv(90000),bsbv(90000),berv(90000),bsigv(90000)
      integer bkin(90000),ncall/0/
      real*8 rhofact
      common/bstuffold/ bptv,bzv,bphiv,bq2v,bxv,byv,byerv,
     >  buv,bubv,bdv,bdbv,bsv,bsbv,berv,bsigv,bitv,brho,
     >  bkin,bffu,bffd,bffs,bchi2k,bdfk,rhofact,bnpt

      ncall = ncall + 1
      chi2 = 0.
      do i=1,56
       bchi2k(i)=0.
       bdfk(i)=0.
      enddo

      fm = 0.005
      p1p = p(1) + fm
      p2p = p(2) + fm
      p3p = p(3) + fm
      p4p = p(4) + fm

      do i=1,bnpt
       it = bitv(i)
       phi = bphiv(i)
       pt = bptv(i)
       pt2 = pt**2
       z0 = bzv(i)
       q2 = bq2v(i)
       x = bxv(i)
       er = berv(i) ! epsilon * R=sigl/sigt
       nu = q2 / 2. / 0.938 / x
       y = nu / 10.6
       u = buv(i) 
       xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
       zp = (z0 / 2.) * (xp / x) *(1. + 
     >     sqrt(1 - 4 * x**2 * am**2 *  
     >     (ampi**2 + pt2) / z0**2 / q2**2))
       z = z0
c use zp instead of z
       z = zp
       z2 = z**2
       if(i .lt. 100 .and. ncall.eq.1)
     >   write(6,'(10f7.3)') xp/x,zp/z0,x,xp,z0,z,zp,q2,am,ampi
c old use of p10
c       d = bdv(i) * p(10)
       d = bdv(i)
       s = bsv(i)
       ubar = bubv(i)
       dbar = bdbv(i)
       sbar = bsbv(i)
       sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)
       if (it.eq.2 .or. it.eq.5) then
          sum_sq = sum_sq + qu**2*(d+dbar) + 
     >         qd**2*(u+ubar) + qs**2*(s+sbar)
       endif

       c0 = -4. * (2. - y) * sqrt(1. - y) * z / sqrt(q2) /  
     >      (1. + (1. - y)**2) * sqrt(pt2) * cos(phi)
       c0 = c0 * p(10)
cxx turned Cahn teerm off!
       c0 = 0.
       c1 = 1. + c0 * p(1) / (z2*p(1) + p(3))
       c2 = 1. + c0 * p(2) / (z2*p(2) + p(4))
       c3 = 1. + c0 * p(1) / (z2*p(1) + p(4))
       c4 = 1. + c0 * p(2) / (z2*p(2) + p(3))
       c1 = c1 / (z2*p(1) + p(3)) / 9.
       c2 = c2 / (z2*p(2) + p(4)) / 9.
       c3 = c3 / (z2*p(1) + p(4)) / 9.
       c4 = c4 / (z2*p(2) + p(3)) / 9.
c these are for deuteron (have Fermi motion added)
       c1p = 1. + c0 * p1p / (z2*p1p + p(3))
       c2p = 1. + c0 * p2p / (z2*p2p + p(4))
       c3p = 1. + c0 * p1p / (z2*p1p + p(4))
       c4p = 1. + c0 * p2p / (z2*p2p + p(3))
       c1p = c1p / (z2*p1p + p(3)) / 9.
       c2p = c2p / (z2*p2p + p(4)) / 9.
       c3p = c3p / (z2*p1p + p(4)) / 9.
       c4p = c4p / (z2*p2p + p(3)) / 9.

c start off with Binneweis frag. func. 
        sv = log( log(Q2/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
        N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
        a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
        a2 = 1.430 + 0.543*sv - 0.023*sv**2
        N = n * p(5)
        a1 = a1 * p(6)
        a2 = a2 * p(7)
	Ns = 4.250 - 3.147*sv + 0.755*sv**2
	a1s = -0.770 -0.573*sv + 0.117*sv**2
	a2s = 4.48 + 0.890*sv - 0.138*sv**2
cxx	zhadm = min(z, 0.75)
        zhadm = min(z, 0.99)
	D_sum = N * zhadm**a1 * (1.0-zhadm)**a2
	D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
        b1 = 0.083583 * p(8)
        b2 = 1.983 * p(9)
	R_D = (1.0-zhadm)**b1 / (1.0+zhadm)**b2
	D_fav = D_sum/(1.0+R_D)
	D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
	D_s = D_sum_s/2.0
        Dp = d_fav
        Dm = D_unfav

c need to add sea!
       rpp = 4.*u*Dp * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >          d*Dm * c2*exp(-pt2/(z2*p(2) + p(4)))+
     >       dbar*Dp * c1*exp(-pt2/(z2*p(1) + p(3)))+
     >    4.*ubar*Dm * c2*exp(-pt2/(z2*p(2) + p(4)))
       rpm = 4.*u*Dm * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >          d*Dp * c4*exp(-pt2/(z2*p(2) + p(3)))+
     >       dbar*Dm * c3*exp(-pt2/(z2*p(1) + p(4)))+
     >    4.*ubar*Dp * c4*exp(-pt2/(z2*p(2) + p(3)))
       rppp= 4.*u*Dp * c1p*exp(-pt2/(z2*p1p + p(3)))+
     >          d*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))+
     >       dbar*Dp * c1p*exp(-pt2/(z2*p2p + p(3)))+
     >    4.*ubar*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))
       rpmp= 4.*u*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >          d*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))+
     >       dbar*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >    4.*ubar*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))
       delu = p(19) * (1 - x)**4 / sqrt(x) * (x-0.0909)
       deld = -1. * delu
       rnp = 4.*(d+deld)*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))+
     >          (u+delu)*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >       ubar*Dp * c1p*exp(-pt2/(z2*p2p + p(3)))+
     >    4.*dbar*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))
       rnm = 4.*(d+deld)*Dm * c2p*exp(-pt2/(z2*p2p + p(4)))+
     >          (u+delu)*Dp * c1p*exp(-pt2/(z2*p1p + p(3)))+
     >       ubar*Dm * c3p*exp(-pt2/(z2*p1p + p(4)))+
     >    4.*dbar*Dp * c4p*exp(-pt2/(z2*p2p + p(3)))

       if(it.eq.1) sig = rpp * 
     >   (1. + p(11)*sqrt(pt2/q2)*cos(phi) +
     >       p(12)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.2) sig = (rppp + rnp) *
     >   (1. + p(13)*sqrt(pt2/q2)*cos(phi) +
     >       p(14)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.4) sig = rpm *
     >   (1. + p(15)*sqrt(pt2/q2)*cos(phi) +
     >       p(16)*sqrt(pt2/q2)*cos(2.*phi)) 
       if(it.eq.5) sig = (rpmp + rnm) * 
     >   (1. + p(17)*sqrt(pt2/q2)*cos(phi) +
     >       p(18)*sqrt(pt2/q2)*cos(2.*phi)) 

       sig = sig / sum_sq / 2. / 3.1415928

! eps * R term
       sig = sig * (1. + er)
       bsigv(i) = sig

! new use of p(10): how much rho subtraction
       y = byv(i) *(1. - p(10)*brho(i))
       chi2 = chi2 + (sig - y)**2 / byerv(i)**2
       j=bkin(i)
       bchi2k(j) = bchi2k(j) + (sig - y)**2 / byerv(i)**2
       bdfk(j) = bdfk(j) + 1.
      enddo

      if((ncall/100)*100.eq.ncall) write(6,'(''ncall, chi2'',i8,f8.3)')
     >  ncall,chi2/float(bnpt)

      fval = chi2

      return
      end

      subroutine simcmodel(x,q2,z,pt,phi,mmpi2,
     >  it,sighad,u,ubar,d,dbar,u1,d1,
     >  s,sbar,ff,fu,fs,dsigdz,yf,yu,zp,rfu,iFF,ipdf)
      implicit none

      real*8 x,z,pt,phi,wsq,w
      integer it,ncall,iFF,ipdf

	integer iset  !which set (1=cteq5m)
	integer ipart !particle u=1, ubar=-1, d=2, dbar=-2, s=3, sbar=-3
	real*8 u,d,ubar,dbar,s,sbar,fav,unfav,rgeiger
	real*8 qu,qd,qs ! u, d, and s quark charges
	real*8 D_fav, D_unfav, D_sum, R_D !favored,unfavored,sum,ratio of FFs
	real*8 D_sum_s, D_s  ! strange frag. functions
	real*8 lambda, Q2zero ! scales for FF param
	real*8 sv !scaling variable for FF param
	real*8 N,a1,a2 !parameters for FF param
	real*8 Ns, a1s, a2s !parameters for strange FF param
        real*8 ff,fu,fs,chk,rfu,q28
C Some local kinematic variables
        real*4 xbj,q2gev, Qgev
	real*8 sx, pt2gev !unitless or GeV
	real*8 b  ! pt2 parameter for FFs

	real*8 nu,qx,qy,qz,mtar,Q2,Eb,Eprime,Epx,Epy,Epz  !MeV
	real*8 pt2,zhad,Ehad,phad,mhad,zhadm ! all in MeV
	real*8 cthpq,pi/3.1415928/

	real*8 kcent,klo,khi,mmpi2,mtargev,nugev
	integer i,nwritten/0/,nprint/0/

	real*8 sum_sq, dsigdz, sigsemi, jacobian, fac, sigma_eepiX
	real*8 sighad, sige, dsigdzold,dsigdzp,dsigdzn
	real*8 dsigdzpold, dsigdznold

	real*8 F1,F2,W1,W2,sin2th2,cos2th2,W2coeff
	real*8 Ctq5Pdf	

c parameters for PB fit of 9/11/2020
c versus zp
c       real*8 pf(8)/   1.2803,   0.0851,   0.8379,   0.1586,
c     >                 0.0140,   0.2133,  -4.4985,   4.1285/
c       real*8 pu(8)/   0.8290,  -0.1416,   0.9869,   0.2559,
c     >                 0.0090,  -1.2306,  -1.5292,   2.4169/

c from fit of 12/5/2020
c      real*8 pf(8)/ 2.1676,   0.4117,   1.4852,   0.1480,
c     >              0.1406,   1.6525,  -8.9585,   8.1429/
c      real*8 pu(8)/ 2.1322,   0.3988,   1.9398,   0.1897,
c     >              0.1337,  -0.7782,  -3.2305,   4.3361/

c from fit versus zp of 12/15/2020, which uses
c z-dependant values of b as below
c      real*8 pf(8)/2.4367,   0.6459,   1.6099,   0.1529,
c     >             0.1439,   1.8692,  -9.9865,   9.3691/
c      real*8 pu(8)/1.7180,   0.3416,   1.8522,   0.2327,
c     >             0.1332,  -1.1574,  -2.0046,   3.4474/
c fit of June 20, 2021
      real*8 pf(12)/   1.0424,  -0.1714,   1.8960,  -0.0307,
     >                 0.1636,  -0.1272,  -4.2093,   5.0103,
     >                 2.7406,  -0.5778,   3.5292,   7.3910/
      real*8 pu(12)/   0.7840,   0.2369,   1.4238,   0.1484,
     >                 0.1518,  -1.2923,  -1.5710,   3.0305,
     >                 1.1995,   1.3553,   2.5868,   8.0666/

	logical first, doinghkns

	parameter (iset=1)
	parameter (qu=2./3.)
	parameter (qd=-1./3.)
	parameter (qs=-1./3.)

	parameter (lambda=0.227) !0.227 GeV for NLO
	parameter (Q2zero=2.0)   !Gev^2 for u,d,s,g

	data first /.TRUE./

c this is for DSS
        integer IHdss,ICdss,IOdss, fini 
        real*8  U1, UB, D1, DB, S1, SB, C1, B1, GL1
        COMMON / FRAGINI / FINI

! this is for HKNS frag. fun.
	integer ixx,iqq,iswap
        real*8 hkusv(100,100),hkdsv(100,100)
        real*8 xp,zp,yu,yf
c parameter from March 2022 fit
        real*8 p(10),ffav,funf,am/0.938/,ampi/0.139/,mpi/0.139/
c
cf  5  0.713  0.061
cf  6  0.455  0.024
cf  7 -0.192  0.018
cf  8 -1.582  0.122
cf  9  1.208  0.056
cf 10  0.791  0.044
       p(5)=  0.715
       p(6)=  0.455
       p(7)= -0.192
       p(8)= -1.582
       p(9)=  1.208
       p(10)=  0.791

        if(first) fini=0

        if(first) then
         open(unit=23,file='FFHKNS07.dat')
         do iqq=1,100
          do ixx=1,100
           read(23,'(f5.1,f5.2,2f9.5)') q2gev,zhad,
     >       hkusv(iqq,ixx),hkdsv(iqq,ixx)
          enddo
         enddo
         close(unit=23)
        endif
	if(first) then
	   call SetCtq5(iset)	! initialize Cteq5 (we're using cteq5m)
	   first=.FALSE.
	endif

	b = 3.8
	mhad = Mpi
	doinghkns = .true.
        xbj = x
        q2gev = q2
	Qgev = sqrt(Q2gev)
        zhad = z
        pt2gev = pt**2

! using xp has almost no effet compared to x
        xp = 2.*x / (1. + sqrt(1. + 4. * x**2 * am**2 / q2))
        xbj = xp

        if(ipdf.eq.1) then
 	ipart=1
	u = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-1
	ubar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=2
	d = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-2
	dbar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=3
	s = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-3
	sbar = Ctq5pdf (ipart , xbj, Qgev)

        else
c         call getjam(xp, q2, u,ubar,d,dbar,s,sbar)
         call getjam(x, q2, u,ubar,d,dbar,s,sbar)
c for test, make d/u smaller at higher x
c this makes p+ / d+ better, but p- / d+ worse!
c         if(x.gt.0.35) d = 0.
        endif

	sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar) + 
     >     qs**2*(s+sbar)

	if (it.eq.2 .or. it.eq.5) then
	   sum_sq = sum_sq + qu**2*(d+dbar) + 
     >       qd**2*(u+ubar) + qs**2*(s+sbar)
	endif

	sv = log( log(Q2gev/lambda**2)/log(Q2zero/lambda**2) )
C Form of parameterization is D = N z^a1 (1-z)^a2
	   N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
	   a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
	   a2 = 1.430 + 0.543*sv - 0.023*sv**2
	   Ns = 4.250 - 3.147*sv + 0.755*sv**2
	   a1s = -0.770 -0.573*sv + 0.117*sv**2
	   a2s = 4.48 + 0.890*sv - 0.138*sv**2
   	   zhadm = min(zhad, 0.75)
cxx	   zhadm = zhad
	   D_sum = N*zhadm**a1*(1.0-zhadm)**a2
c correction to get better agreement with data
	   D_sum = D_sum * (1.5 - 0.5 * zhad)
	   D_sum_s = Ns*zhadm**a1s*(1.0-zhadm)**a2s
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
	   R_D = (1.0-zhad)**0.083583/(1.0+zhad)**1.9838
	   R_D = (1.0-zhadm)**0.083583/(1.0+zhadm)**1.983
           rgeiger = r_d
c New from Peter's fit of August 27,2020
           R_D = (1. - 0.8 * zhad) * (1. + 
     >      zhad    * 1.988 + 
     >      zhad**2 * (-9.014) +
     >      zhad**3 * 8.219)
	   D_fav = D_sum/(1.0+R_D)
	   D_unfav = D_sum/(1.0+1.0/R_D)
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
	   D_s = D_sum_s/2.0

c new PB it using zp for pions. This is z * D
        wsq = am**2 + q2 * (1/x -1)
        w = sqrt(wsq) 
         xp = 2.*xbj / (1. + sqrt(1. + 4. * xbj**2 * am**2 / q2gev))
         zp = (zhad / 2.) * (xp / xbj) *(1. + 
     >     sqrt(1 - 4 * xbj**2 * am**2 *  
     >     (ampi**2 + pt2gev) / zhad**2 / q2gev**2))
         sv = log(q2gev/2.)
c         yf = pf(1) * zp**(pf(2) + pf(4)*sv) * 
c     >        (1.-zp)**(pf(3) + pf(5)*sv) 
c         yf = yf * (1. + pf(6)*zp + pf(7)*zp**2 + pf(8)*zp**3)
c         yu = pu(1) * zp**(pu(2) + pu(4)*sv) * 
c     >        (1.-zp)**(pu(3) + pu(5)*sv) 
c         yu = yu * (1. + pu(6)*zp + pu(7)*zp**2 + pu(8)*zp**3)
       yf = pf(1) * zp**(pf(2) + pf(4)*sv + pf(9)/w) * 
     >       (1.-zp)**(pf(3) + pf(5)*sv + pf(10)/w) 
       yf = yf * 
     >      (1. + pf(6)*zp + pf(7)*zp**2 + pf(8)*zp**3) *
     >      (1. + pf(11)/w + pf(12)/w**2)
       yu = pu(1) * zp**(pu(2) + pu(4)*sv + pu(9)/w) * 
     >       (1.-zp)**(pu(3) + pu(5)*sv + pu(10)/w) 
       yu = yu * 
     >      (1. + pu(6)*zp + pu(7)*zp**2 + pu(8)*zp**3) *
     >      (1. + pu(11)/w + pu(12)/w**2)

         if(it.le.3) then
c forget about hkns for now
c 	  u1 = D_fav * zhad
c	  d1 = D_unfav * zhad
c new fit
	  u1 = yf
	  d1 = yu
	 else
c	  d1 = zhad * D_fav
c	  u1 = zhad * D_unfav
c new fit
	  u1 = yu
	  d1 = yf
	 endif
         rfu = yu / yf

	 ub = d1
	 db = u1
c	 s1 = zhad * D_s
c	 sb = s1
c new fit
	 s1 = yu
	 sb = s1

c override with JAM
        if(iFF.eq.4) then
         q28 = q2gev
         IHdss=1
         ICdss=1
         if(it.gt.3) ICdss = -1
         IOdss=0
! use fDSS for sea quarks
         call fDSS (IHdss, ICdss, IOdss, Zp, Q28, 
     >       U1, UB, D1, DB, S1, SB, C1, B1, GL1)
! JAM for u1, d1
         if(it.le.3) then
          call jamff(zp,q28,d1,u1) ! pi+
         else
          call jamff(zp,q28,u1,d1) ! pi-
         endif
        endif

c use DSS
        if(iFF.eq.2.or.iff.eq.3) then
         q28 = q2gev
         IHdss=1
         ICdss=1
         if(it.gt.3) ICdss = -1
         IOdss=0
         call fDSS (IHdss, ICdss, IOdss, Zp, Q28, 
     >       U1, UB, D1, DB, S1, SB, C1, B1, GL1)
        endif
        iswap=0
        if(iff.eq.3) iswap=1

c put in first-pass corr. from fit
        if(it.le.3) then
         ffav = 1. + p(5)/wsq + p(6)/mmpi2 + p(7)/q2
         funf = 1. + p(8)/wsq + p(9)/mmpi2 + p(10)/q2
        else
         funf = 1. + p(5)/wsq + p(6)/mmpi2 + p(7)/q2
         ffav = 1. + p(8)/wsq + p(9)/mmpi2 + p(10)/q2
        endif

        if(ncall.lt.10) then
         write(6,'(''ffav'',i2,6f7.3)') it,wsq,mmpi2,
     >    q2,ffav,funf
         ncall = ncall + 1
        endif

cturn off for test
        ffav = 1.
        funf = 1.
           
c for test. Doesnt fix pi+/d
c        if(it.le.3) d1 = d1 * 1.15
c        if(it.gt.3) u1 = u1 * 1.15
        dsigdz = (qu**2 * u    * u1 * ffav + 
     >            qu**2 * ubar * ub * funf +
     >  	  qd**2 * d    * d1 * funf + 
     >            qd**2 * dbar * db * ffav + 
     >  	  qs**2 * s    * s1 + 
     >            qs**2 * sbar * sb)/sum_sq/zhad
        if(it.le.3) then
 	 ff = (qu**2 * u + qd**2 * dbar)/sum_sq/zhad
	 fu = (qd**2 * d + qu**2 * ubar)/sum_sq/zhad
        else
 	 fu = (qu**2 * u + qd**2 * dbar)/sum_sq/zhad
	 ff = (qd**2 * d + qu**2 * ubar)/sum_sq/zhad
        endif
	fs = (qs**2 * s + qs**2 * sbar)/sum_sq/zhad
	dsigdzp = dsigdz
	if(it.eq.2.or.it.eq.5) then ! Deut
         if(iswap.ne.1) then
         dsigdzn =(qu**2 * d    * u1 * ffav + 
     >             qu**2 * dbar * ub * funf +
     >  	   qd**2 * u    * d1 * funf + 
     >             qd**2 * ubar * db * ffav + 
     >  	   qs**2 * s    * s1 + 
     >             qs**2 * sbar * sb)/sum_sq/zhad
         else
         dsigdzn =(qu**2 * d    * d1 * ffav + 
     >             qu**2 * dbar * ub * funf +
     >  	   qd**2 * u    * u1 * funf + 
     >             qd**2 * ubar * db * ffav + 
     >  	   qs**2 * s    * s1 + 
     >             qs**2 * sbar * sb)/sum_sq/zhad
         endif
	 dsigdz = dsigdzp + dsigdzn
         if(it.le.3) then
    	 ff = ff + (qu**2 * d + qd**2 * ubar)/sum_sq/zhad
	 fu = fu + (qd**2 * u + qu**2 * dbar)/sum_sq/zhad
         else
    	 fu = fu + (qu**2 * d + qd**2 * ubar)/sum_sq/zhad
	 ff = ff + (qd**2 * u + qu**2 * dbar)/sum_sq/zhad
         endif
	 fs = fs +(qs**2 * s + qs**2 * sbar)/sum_sq/zhad
        endif
 	chk = ff*u1 + fu*d1 + fs*s1
        nprint = nprint+1
        if(nprint.lt.10) then
         write(6,'(/''chk'',9f7.3)') 
     >    chk/dsigdz,chk,dsigdz,
     >    ff,fu,fs,u1,d1,s1    
         write(6,'(i2,9f7.3)') 
     >    it,xbj,q2gev,zhad,u,ubar,d,dbar,s,sbar    
        endif

c b as function of z
        b = 1./ (0.200 * zhad**2 + 0.200)
c used in SIMC
	sighad = dsigdz * b * exp(-b * pt2gev)/ 2. / pi

      return
      end
C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.

C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      real*8 Function Ctq5Pdf (Iparton, X, Q)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      real*8 partonx
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
        Print *, 'Setting to Alambda'
        Q=Alambda
c	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

      real*8 FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
c      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      implicit integer (I-N)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
c         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
c     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq5 (Iset)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      Parameter (Isetmax=9)

      Character Flnm(Isetmax)*17, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5/cteq5m.tbl', 'cteq5/cteq5d.tbl', 'cteq5/cteq5l.tbl', 
     >   'cteq5/cteq5hj.tbl'
     > , 'cteq5/cteq5hq.tbl', 'cteq5/cteq5f3.tbl', 'cteq5/cteq5f4.tbl'
     > , 'cteq5/cteq5m1.tbl', 'cteq5/ctq5hq1.tbl'  /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
C                             ********************
      End

      Subroutine ReadTbl (Nu)
c      Implicit Double Precision (A-H,O-Z)
      implicit integer (I-N)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
      integer N, NextUn
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
 
c      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) stop
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      real*8 function futil()
      return
      end

********************************************************************
*                                                                  *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
*  D.de Florian, R.Sassot, M.Stratmann   Phys.Rev.D75 114010 2007  *
*                                 *and*  Phys.Rev.D76 074033 2007  *
*                                                                  *
*     CALL fDSS (IH,IC,IO, X, Q2, U, UB, D, DB, S, SB, C, B, GL)   *
*                                                                  *	
*  INPUT:                                                          *
*  IH = hadron type    1: PION                                     *
*                      2: KAON                                     *
*                      3: PROTON                                   *
*                      4: CHARGED HADRONS                          *
*                                                                  *
*  IC = Hadron Charge  0: 0 (as average of + and -)                *
*                      1: +                                        *
*                     -1: -                                        *
*                                                                  *
*  IO= Order           0: LO                                       *
*                      1: NLO                                      *
*                                                                  *
*            X                    (between  0.05   and  1.0)       *
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
*             (for values outside the allowed range the program    *
*              writes a warning and extrapolates to the x and      *
*              Q2 values requested)                                *
*                                                                  *
*   OUTPUT: U, UB, D, DB, S, SB,   C,           B,       GL        *
*           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon      *
*           Always X times the distribution is returned            *
*                                                                  *
*                                                                  *
*   COMMON:  The main program or the calling routine has to have   *
*            a common block  COMMON / FRAGINI / FINI , and  FINI   *
*            has always to be zero when DSS is called for the      *
*            first time or when the SET has been changed.          *
*                                                                  *
********************************************************************

      SUBROUTINE fDSS (IH,IC,IO, X, Q2, U, UB, D, DB, 
     >  S, SB, C, B, GL)
      implicit real*8 (A-H,O-Z)
      PARAMETER (NPART=9, NX=35, NQ=24, NARG=2)
      DIMENSION XUTOTF(NX,NQ), XDTOTF(NX,NQ), XSTOTF(NX,NQ)
      DIMENSION XUVALF(NX,NQ), XDVALF(NX,NQ), XSVALF(NX,NQ)
      DIMENSION XCTOTF(NX,NQ), XBTOTF(NX,NQ)
      DIMENSION XGF(NX,NQ), PARTON (NPART,NQ,NX-1)
      DIMENSION QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      integer fini
      COMMON / FRAGINI / FINI
      SAVE XUTOTF, XDTOTF, XSTOTF, XCTOTF, XBTOTF, XGF, NA, ARRF
      SAVE XUVALF, XDVALF, XSVALF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.d0, 1.25D0, 1.5D0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
       DATA XB /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6        0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9 , 0.93, 1.0/
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.05D0) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C          STOP
       ENDIF
       IF ( (Q2.LT.1.D0) .OR. (Q2.GT.1.D5) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
C          STOP
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
      IF (FINI.NE.0) GOTO 16
      IF ((IH.EQ.1).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='PINLO.GRID')
      ELSEIF ((IH.EQ.1).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='PILO.GRID')
      ELSEIF ((IH.EQ.2).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='KANLO.GRID')
       ELSEIF ((IH.EQ.2).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='KALO.GRID')
      ELSEIF ((IH.EQ.3).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='PRONLO.GRID')       
      ELSEIF ((IH.EQ.3).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='PROLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='HNLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='HLO.GRID')       
	ELSE
         WRITE(6,93)
 93      FORMAT (2X,' WRONG SET')
         STOP
      END IF
C
       DO 15 M = 1, NX-1 
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M), PARTON(9,N,M)
  90   FORMAT (9(1PE10.3))
c       write(6,'(2i4,8f8.2)') m,n,PARTON(1,N,M), PARTON(2,N,M), 
c     >   PARTON(3,N,M),PARTON(4,N,M) 
  15   CONTINUE
       CLOSE(IIREAD)
C
      FINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUTOTF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**4 * XB0**0.5)
        XDTOTF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.5)
        XSTOTF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**4 * XB0**0.5) 
        XCTOTF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0**0.3) 
        XBTOTF(IX,IQ) = PARTON(5,IQ,IX) / (XB1**7 * XB0**0.3)
        XGF(IX,IQ)    = PARTON(6,IQ,IX) / (XB1**4 * XB0**0.3)
        XUVALF(IX,IQ) = PARTON(7,IQ,IX) / (XB1**4 * XB0**0.5)
        XDVALF(IX,IQ) = PARTON(8,IQ,IX) / (XB1**4 * XB0**0.5)
        XSVALF(IX,IQ) = PARTON(9,IQ,IX) / (XB1**4 * XB0**0.5)
  20  CONTINUE
        XUTOTF(NX,IQ) = 0.D0
        XDTOTF(NX,IQ) = 0.D0
        XSTOTF(NX,IQ) = 0.D0
        XCTOTF(NX,IQ) = 0.D0
        XBTOTF(NX,IQ) = 0.D0
        XGF(NX,IQ)    = 0.D0
        XUVALF(NX,IQ) = 0.D0
        XDVALF(NX,IQ) = 0.D0
        XSVALF(NX,IQ) = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = LOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = LOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = LOG(X)
      XT(2) = LOG(Q2)
      UTOT = FINT(NARG,XT,NA,ARRF,XUTOTF) * (1.D0-X)**4 * X**0.5
      DTOT = FINT(NARG,XT,NA,ARRF,XDTOTF) * (1.D0-X)**4 * X**0.5 
      STOT = FINT(NARG,XT,NA,ARRF,XSTOTF) * (1.D0-X)**4 * X**0.5
      CTOT = FINT(NARG,XT,NA,ARRF,XCTOTF) * (1.D0-X)**7 * X**0.3
      BTOT = FINT(NARG,XT,NA,ARRF,XBTOTF) * (1.D0-X)**7 * X**0.3
      GL   = FINT(NARG,XT,NA,ARRF,XGF)    * (1.D0-X)**4 * X**0.3
      UVAL = FINT(NARG,XT,NA,ARRF,XUVALF) * (1.D0-X)**4 * X**0.5
      DVAL = FINT(NARG,XT,NA,ARRF,XDVALF) * (1.D0-X)**4 * X**0.5 
      SVAL = FINT(NARG,XT,NA,ARRF,XSVALF) * (1.D0-X)**4 * X**0.5
       
       Up  = (UTOT+UVAL)/2.
       UBp = (UTOT-UVAL)/2.
       Dp  = (DTOT+DVAL)/2.
       DBp = (DTOT-DVAL)/2.
       Sp  = (STOT+SVAL)/2.
       SBp = (STOT-SVAL)/2.
       Cp  =  CTOT/2.
       Bp  =  BTOT/2.
              
       IF (IC.EQ.1) THEN
       U  = Up
       UB = UBp
       D  = Dp 
       DB = DBp
       S  = Sp
       SB = SBp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.-1) THEN
       U  = UBp
       UB = Up
       D  = DBp 
       DB = Dp
       S  = SBp
       SB = Sp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.0) THEN
       U  = (UBp+Up)/2.
       UB =  U
       D  = (DBp+Dp)/2. 
       DB =  D
       S  = (SBp+Sp)/2.
       SB =  S
       C  =  Cp
       B  =  Bp 
       ELSE
         WRITE(6,94)
 94      FORMAT (2X,' WRONG CHARGE')
         STOP
       END IF 
c       write(6,'(i2,8f7.3)') ic,up,u,(utot+uval)/2.,
c     >  d, dp, (DTOT+DVAL)/2.
 60   RETURN
       END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(59),TABLE(840)
      DIMENSION D(2),NCOMB(2),IENT(2)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
!---------------------------------------------------------------------
      subroutine FNP_NMC(X,QSQ,rat)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,QSQ,A,B,X2,X3,rat
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      rat = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END

      subroutine getjam(x, q2, u,ubar,d,dbar,s,sbar)
      implicit none
      logical first/.true./
      real*8 x, q2, Q, u,ubar,d,dbar,s,sbar
      real*8 usv(51,13),ubsv(51,13)
      real*8 dsv(51,13),dbsv(51,13)
      real*8 ssv(51,13),sbsv(51,13)
      real*8 xsv(51),qsv(13)
      integer ix,iq,j

      if(first) then
       first=.false.
       open(unit=77, file='jampdf.dat')
       do ix=1,51
        do iq=1,13
         read(77,'(2f6.3,11f8.4)') xsv(ix),qsv(iq),
     >     sbsv(ix,iq),ubsv(ix,iq),dbsv(ix,iq),
     >     dsv(ix,iq), usv(ix,iq),ssv(ix,iq)
        enddo
       enddo
       close(unit=77)
      endif

      ix=0
      iq=0
      do j=1,50
       if(x.ge.xsv(j) .and. x.lt.xsv(j+1)) ix=j
      enddo
      q = sqrt(q2)
      do j=1,12
       if(q.ge.qsv(j) .and. q.lt.qsv(j+1)) iq=j
      enddo

      if(ix.eq.0 .or. iq.eq.0) then
       write(6,'(''error in jampdf'',2f8.3)') x,q
       stop
      endif

      u = usv(ix,iq)
      d = dsv(ix,iq)
      s = ssv(ix,iq)
      ubar = ubsv(ix,iq)
      dbar = dbsv(ix,iq)
      sbar = sbsv(ix,iq)

      return
      end

      SUBROUTINE FFFIT_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! fit for two favored, one unfavored FF (3 params)
! CALCULATE CHISQ FOR MINUIT
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call,iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 fs1, fsb
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
c      mpp = 4 * u * p(1) + 
c     >          d * p(3) + 
c     >      4.*ub * p(2) +
c     >         db * p(4) + 
c     >          s * ffs + 
c     >         sb * ffsb 

c      mpm = 4 * u * p(2) + 
c     >          d * p(4) + 
c     >      4.*ub * p(1) +
c     >         db * p(3) + 
c     >          s * ffsb + 
c     >         sb * ffs 
c      mnp = 4 * d * p(1) + 
c     >          u * p(3) + 
c     >      4.*db * p(2) +
c     >         ub * p(4) + 
c     >          s * ffs + 
c     >         sb * ffsb 

c      mnm = 4 * d * p(2) + 
c     >          u * p(4) + 
c     >      4.*db * p(1) +
c     >         ub * p(3) + 
c     >          s * ffsb + 
c     >         sb * ffs 

      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 

      mpm = 4 * u * p(2) + 
     >          d * p(3) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 

      if(iswap.ne.1) then
      mnp = 4 * d * p(1) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * d * p(2) + 
     >          u * p(3) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      else
      mnp = 4 * d * p(3) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >      4.*db * p(3) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db) + u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mpm,
     >   msv(3),msver(3),mdp,msv(4),msver(4),mdm
       write(6,'(''fchk'',i4,12f7.3)') ncall,p(1),p(2),p(2),p(3),
     >   ffs,ffsb
       write(6,'(''fchk'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFIT4_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fits (two favored, two unfavored)
! depending on "swap", interchange for p, n or not
! CALCULATE CHISQ FOR MINUIT
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 
      mpm = 4 * u * p(3) + 
     >          d * p(4) + 
     >      4.*ub * p(4) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 
c regular (FF don't depend on target)
      if(iswap.eq.0) then
      mnp = 4 * d * p(1) + 
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(3) + 
     >          u * p(4) + 
     >      4.*db * p(4) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
c swap the two favored FF and the two unfavored FF
c as it looks like happens in LUND
      else
      mnp = 4 * d * p(4) + 
     >          u * p(3) + 
     >      4.*db * p(3) +
     >         ub * p(4) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >         db * p(1) +
     >      4.*ub * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db)+ u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFITcsv_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fit assumes twofavored, one
! unfavored FF, and delta u = -delta d
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2,delta_u,delta_d
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 

      mpm = 4 * u * p(2) + 
     >          d * p(3) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 

      delta_d = p(4) * (u + d) * 3. / 8.
      delta_u = -1.0 * delta_d
      if(iswap.ne.1) then
      mnp = 4 * (d + delta_d) * p(1) + 
     >          (u + delta_u) * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * (d + delta_d) * p(2) + 
     >          (u + delta_u) * p(3) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      else
      mnp = 4 * (d + delta_d) * p(3) + 
     >          (u + delta_u) * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * (d + delta_d) * p(2) + 
     >          (u + delta_u) * p(1) + 
     >      4.*db * p(3) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
c      sumn = 4.*(d * (1. + p(3)) + db)+ 
c     >           u * (1. + p(4)) + ub + s + sb
      sumn = 4.*(d + db)+ 
     >           u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFITex_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fit assume only 1 favored, one
! unfavored FF,and one enhancemeent factor for ech
! of pi+ from p and pi- from n
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 
      mpm = 4 * u * p(2) + 
     >          d * p(1) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 
c regular (FF don't depend on target)
c p3 is delta d / d, p4 is delta u/u
      mnp = 4 * d * p(1) +
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      mpp = mpp  * p(3)
      mnm = mnm  * p(4)

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db)+ u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      subroutine jamff(z,q2,d1,u1)
c frag. func. from JAM 2020. Table generaated from their files
c using jamtest.f
c Note: in their model ubar=u=u1, dbar=d=d1, s=sbar=d1
c d1 is "unfavored", u1 is "favored"
      implicit none
      integer j,k,jj,kk,i
      real*8 zsv(51),qsv(13),ff(51,13,2),ffz(2,2),q
      real*8 z,q2,d1,u1
      logical first/.true./

      if(first) then
       open(unit=777,file='jamff.dat')
       do j=1,44
        do k=1,13
         read(777,'(4f9.5)') zsv(j),qsv(k),ff(j,k,1),ff(j,k,2)
        enddo
       enddo
       first = .false.
      endif

      jj=0
      do j=1,43
       if(z.gt.zsv(j).and.z.le.zsv(j+1)) jj=j
      enddo
      if(z.le.zsv(1)) jj=1
      if(z.gt.zsv(44)) jj=44

      q = sqrt(q2)
      kk=0
      do k=1,12
       if(q.gt.qsv(k).and.q.le.qsv(k+1)) kk=k
      enddo
      if(q.le.qsv(1)) kk=1
      if(q.gt.qsv(12)) kk=12

c interpolate in z for 2 q2 bins
      do i=1,2
       ffz(1,i) = ff(jj, kk, i) + 
     >  (ff(jj+1, kk, i) - ff(jj, kk, i)) *
     >  (z - zsv(jj)) / 
     >  (zsv(jj+1) - zsv(jj))
       ffz(2,i) = ff(jj, kk+1, i) + 
     >  (ff(jj+1, kk+1, i) - ff(jj, kk+1, i)) *
     >  (z - zsv(jj)) / 
     >  (zsv(jj+1) - zsv(jj))
      enddo
c interpolate in q
       d1 = ffz(1,1) + 
     >  (ffz(2,1) - ffz(1,1)) *
     >  (q - qsv(kk)) / 
     >  (qsv(kk+1) - qsv(kk))
       u1 = ffz(1,2) + 
     >  (ffz(2,2) - ffz(1,2)) *
     >  (q - qsv(kk)) / 
     >  (qsv(kk+1) - qsv(kk))

       return
       end

      mnp = 4 * (d + delta_d) * p(1) + 
     >          (u + delta_u) * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(3) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * (d + delta_d) * p(2) + 
     >          (u + delta_u) * p(3) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      else
      mnp = 4 * (d + delta_d) * p(3) + 
     >          (u + delta_u) * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 

      mnm = 4 * (d + delta_d) * p(2) + 
     >          (u + delta_u) * p(1) + 
     >      4.*db * p(3) +
     >         ub * p(2) + 
     >          s * ffsb + 
     >         sb * ffs 
      endif

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
c      sumn = 4.*(d * (1. + p(3)) + db)+ 
c     >           u * (1. + p(4)) + ub + s + sb
      sumn = 4.*(d + db)+ 
     >           u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      SUBROUTINE FFFITex_FCN(NPAR,GRAD,FVAL,P,IFLAG,FUTIL)
! four-parameter fit assume only 1 favored, one
! unfavored FF,and one enhancemeent factor for ech
! of pi+ from p and pi- from n
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD(*)
      REAL*8 P(*) ! VECTOR OF PARAMETERS
      REAL*8 FVAL  ! CHISQ
      REAL*8 FUTIL ! AUXIALLY FUNCTION
      EXTERNAL FUTIL
      INTEGER I,J,IT,n call, iswap
      REAL*8 CHI2, U, D, UB, DB, S, SB, FFS, FFSB, sump, sumd
      REAL*8 mpp,mpm,mnp,mnm,mdp,mdm,msv(4),msver(4),sumn
      real*8 msvfit(4),sf1,sf2
      common/ffstuff/ msv,msver,u,d,ub,db,s,sb,ffs,ffsb,
     >   msvfit,ncall,iswap,chi2

c get sidis multiplicities
      mpp = 4 * u * p(1) + 
     >          d * p(2) + 
     >      4.*ub * p(2) +
     >         db * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf1=     (s * ffs + 
     >         sb * ffsb)/mpp 
      mpm = 4 * u * p(2) + 
     >          d * p(1) + 
     >      4.*ub * p(1) +
     >         db * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      sf2=     (s * ffs + 
     >         sb * ffsb)/mpm 
c regular (FF don't depend on target)
c p3 is delta d / d, p4 is delta u/u
      mnp = 4 * d * p(1) +
     >          u * p(2) + 
     >      4.*db * p(2) +
     >         ub * p(1) + 
     >          s * ffs + 
     >         sb * ffsb 
      mnm = 4 * d * p(2) + 
     >          u * p(1) + 
     >      4.*db * p(1) +
     >         ub * p(2) + 
     >          s * ffs + 
     >         sb * ffsb 
      mpp = mpp  * p(3)
      mnm = mnm  * p(4)

      mdp = (mpp + mnp) 
      mdm = (mpm + mnm) 

      sump = 4.*(u + ub) + d + db + s + sb
      sumn = 4.*(d + db)+ u + ub + s + sb
      sumd = sump + sumn

! these should be same as z * Mult (integrated over phi*)
      mpp = mpp / sump
      mpm = mpm / sump
      mdp = mdp / sumd
      mdm = mdm / sumd

      chi2 = (msv(1) - mpp)**2 / msver(1)**2 +
     >       (msv(2) - mdp)**2 / msver(2)**2 +
     >       (msv(3) - mpm)**2 / msver(3)**2 +
     >       (msv(4) - mdm)**2 / msver(4)**2 

      msvfit(1) = mpp
      msvfit(2) = mdp
      msvfit(3) = mpm
      msvfit(4) = mdm

      if(ncall.eq.0 .and. chi2.lt.0.1) then
       ncall = ncall + 1

       write(6,'(/''fchk4'',i4,12f7.1)') ncall,chi2
       write(6,'(''fchk4'',i4,12f7.3)') ncall,
     >   msv(1),msver(1),mpp,msv(2),msver(2),mdp,
     >   msv(3),msver(3),mpm,msv(4),msver(4),mdm
       write(6,'(''fchk4'',i4,12f7.3)') ncall,p(1),p(2),p(3),p(4),
     >   ffs,ffsb
       write(6,'(''fchk4'',i4,12f7.3)') ncall,u,ub,d,db,s,sb
       write(6,'(''sf'',2f9.5)') sf1, sf2
      endif


       fval = chi2
       return
       end

      subroutine jamff(z,q2,d1,u1)
c frag. func. from JAM 2020. Table generaated from their files
c using jamtest.f
c Note: in their model ubar=u=u1, dbar=d=d1, s=sbar=d1
c d1 is "unfavored", u1 is "favored"
      implicit none
      integer j,k,jj,kk,i
      real*8 zsv(51),qsv(13),ff(51,13,2),ffz(2,2),q
      real*8 z,q2,d1,u1
      logical first/.true./

      if(first) then
       open(unit=777,file='jamff.dat')
       do j=1,44
        do k=1,13
         read(777,'(4f9.5)') zsv(j),qsv(k),ff(j,k,1),ff(j,k,2)
        enddo
       enddo
       first = .false.
      endif

      jj=0
      do j=1,43
       if(z.gt.zsv(j).and.z.le.zsv(j+1)) jj=j
      enddo
      if(z.le.zsv(1)) jj=1
      if(z.gt.zsv(44)) jj=44

      q = sqrt(q2)
      kk=0
      do k=1,12
       if(q.gt.qsv(k).and.q.le.qsv(k+1)) kk=k
      enddo
      if(q.le.qsv(1)) kk=1
      if(q.gt.qsv(12)) kk=12

c interpolate in z for 2 q2 bins
      do i=1,2
       ffz(1,i) = ff(jj, kk, i) + 
     >  (ff(jj+1, kk, i) - ff(jj, kk, i)) *
     >  (z - zsv(jj)) / 
     >  (zsv(jj+1) - zsv(jj))
       ffz(2,i) = ff(jj, kk+1, i) + 
     >  (ff(jj+1, kk+1, i) - ff(jj, kk+1, i)) *
     >  (z - zsv(jj)) / 
     >  (zsv(jj+1) - zsv(jj))
      enddo
c interpolate in q
       d1 = ffz(1,1) + 
     >  (ffz(2,1) - ffz(1,1)) *
     >  (q - qsv(kk)) / 
     >  (qsv(kk+1) - qsv(kk))
       u1 = ffz(1,2) + 
     >  (ffz(2,2) - ffz(1,2)) *
     >  (q - qsv(kk)) / 
     >  (qsv(kk+1) - qsv(kk))

       return
       end
