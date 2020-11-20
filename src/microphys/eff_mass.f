      subroutine eff_mass (nn,np,effmn,effmp)
c     suffix n,p stands for neutron or proton
c     nn, np input are number densities (/cm**3)
c     effm (output) are effective mass/bare mass.
c     Coding is primarily Table Look-up and interpolation

      real*8 nn,np,effmn,effmp
      real*4 knt(42),efnt(42),kpt(42),efpt(42),kn,kp,d,omd
      integer*4 in,ip,inn,ipn

      save knt,efnt,kpt,efpt

      data kpt /
     $ 0.0000, 0.0053, 0.0374, 0.0695, 0.1016, 0.1342, 0.1675, 0.2017, 
     $ 0.2373, 0.2744, 0.3133, 0.3544, 0.3979, 0.4441, 0.4933, 0.5457, 
     $ 0.6017, 0.6615, 0.7255, 0.7938, 0.8668, 0.9448, 1.0187, 1.0867, 
     $ 1.1489, 1.2060, 1.2583, 1.3061, 1.3500, 1.3903, 1.4275, 1.4618, 
     $ 1.4938, 1.5239, 1.5524, 1.5798, 1.6064, 1.6327, 1.6592, 1.6861, 
     $ 1.7139, 1.7430/ 

      data efpt /
     $ 1.0000, 0.9947, 0.9518, 0.9112, 0.8728, 0.8366, 0.8026, 0.7706, 
     $ 0.7408, 0.7129, 0.6869, 0.6629, 0.6407, 0.6204, 0.6018, 0.5849, 
     $ 0.5697, 0.5561, 0.5440, 0.5335, 0.5244, 0.5168, 0.5106, 0.5053, 
     $ 0.5008, 0.4971, 0.4941, 0.4918, 0.4900, 0.4888, 0.4881, 0.4878, 
     $ 0.4879, 0.4884, 0.4891, 0.4901, 0.4912, 0.4924, 0.4937, 0.4950, 
     $ 0.4962, 0.4973/

      data knt /
     $ 0.0000, 0.0159, 0.0683, 0.1172, 0.1628, 0.2056, 0.2458, 0.2839, 
     $ 0.3201, 0.3548, 0.3883, 0.4210, 0.4532, 0.4852, 0.5175, 0.5502, 
     $ 0.5839, 0.6187, 0.6551, 0.6934, 0.7339, 0.7770, 0.8227, 0.8703, 
     $ 0.9197, 0.9705, 1.0224, 1.0752, 1.1284, 1.1819, 1.2353, 1.2883, 
     $ 1.3406, 1.3920, 1.4420, 1.4904, 1.5370, 1.5813, 1.6232, 1.6622, 
     $ 1.6982, 1.7307/ 

      data efnt /
     $ 1.0000, 0.9947, 0.9923, 0.9898, 0.9873, 0.9846, 0.9817, 0.9787, 
     $ 0.9754, 0.9719, 0.9681, 0.9639, 0.9594, 0.9545, 0.9492, 0.9435, 
     $ 0.9372, 0.9304, 0.9231, 0.9151, 0.9066, 0.8973, 0.8877, 0.8779, 
     $ 0.8680, 0.8582, 0.8483, 0.8387, 0.8292, 0.8200, 0.8111, 0.8027, 
     $ 0.7947, 0.7872, 0.7804, 0.7743, 0.7689, 0.7644, 0.7607, 0.7580, 
     $ 0.7563, 0.7558/

c           3 pi**2/1.e+39
      kn = (2.96088d-38*nn)**0.333333333333333
      kp = (2.96088d-38*np)**0.333333333333333
      
      if ( kn.le.0.) then
          effmn = efnt(1)
      elseif ( kn.ge.knt(42) ) then
          effmn = efnt(42)
      else
          do 1 in = 1,41
              if ( knt(in).le.kn.and.kn.le.knt(in+1) ) then
                  inn = in
                  go to 2
              endif
    1     continue
          write(6,'(a)') 'Got to end of 1 loop in eff_mass'
          inn = 1
    2     continue
          d = (kn - knt(inn))/(knt(inn+1) - knt(inn))
          omd = 1. - d
          effmn = d*efnt(inn+1) + omd*efnt(inn)
      endif

      if ( kp.le.0.) then
          effmp = efpt(1)
      elseif ( kp.ge.kpt(42) ) then
          effmp = efpt(42)
      else
          do 3 ip = 1,41
              if ( kpt(ip).le.kp.and.kp.le.kpt(ip+1) ) then
                  ipn = ip
                  go to 4
              endif
    3     continue
          write(6,'(a)') 'Got to end of 3 loop in eff_mass'
          ipn = 1
    4     continue
          d = (kp - kpt(ipn))/(kpt(ipn+1) - kpt(ipn))
          omd = 1. - d
          effmp = d*efpt(ipn+1) + omd*efpt(ipn)
      endif

      return
      end
