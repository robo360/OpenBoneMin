#PROB]

# Authored by Roberto Berwa

# 24-Oct-2020

#[GLOBAL]
#macro max(a,b) ((a) > (b) ? (a) : (b))
#macro F11 T85
#macro PicOCkin Pic0
#macro SETINIT if(NEWIND <=1)
#macro PKCP (PKCENT/PKVC)
#macro DENCP (DENCENT/DENVC)
#macro CaConc0 (P_0/V1)
#macro PTHconc0 (PTH_0/V1)
#macro OB (OBfast + OBslow)
#macro PTHconc (PTH/V1)
#macro CaConc (P/V1)
#macro C1 (P/V1)
#macro C2 (ECCPhos/V1)
#macro C8 (B/V1)
#macro D  ROB1
#macro Osteoclast OC
#macro OB0 (OBfast_0 + OBslow_0)
#macro Calcitriol0 (B_0/V1)
#denosumab concentration (mol)
#macro DENMOL (DENCENT/DENVC/150000)*1000*1000
# [MAIN]
#
# TGFB_0 = Pic0*1000.0
# TGFBact_0 = Pic0
# OBfast_0 = OBtot0*FracOBfast
# OBslow_0 = OBtot0*(1-FracOBfast)
# OB0 = (OBfast_0 + OBslow_0)
# M_0 = k3*RNK_0*L_0/k4
# N_0 = k1*O_0*L_0/k2
# AOH_0 = B_0/10.0
# F_DENSC = DENF*1E6
# F_TERISC = TERIF
#
#
# [PARAM]
using DifferentialEquations

function modelODE!(u, p, t)
  OBtot0 = 0.00501324
  k1 = 0.00000624
  k2 = 0.112013
  k3 = 0.00000624
  k4 = 0.112013
  V1= 14.0
  CaDay = 88.0
  FracJ14 = 0.107763
  J14OCmax = 0.543488
  J14OCgam = 1.6971
  FracJ15 = 0.114376
  kinRNKgam = 0.151825
  koutRNK = 0.00323667
  MOCratioGam = 0.603754
  Da = 0.7/24.0
  OBtgfGAM = 0.0111319
  koutTGF0 = 0.0000298449
  koutTGFGam = 0.919131
  OCtgfGAM = 0.593891
  EmaxPicROB = 3.9745
  PicROBgam = 1.80968
  FracPicROB = 0.883824
  PicOBgam = 0.122313
  FracPicOB = 0.000244818
  EmaxPicOB = 0.251636
  E0Meff = 0.388267
  EmaxMeffOC = 3.15667
  kinOCgam = 8.53065
  EmaxPicOC = 1.9746
  FracPicOC = 0.878215
  PicOCgam = 1.0168
  E0RANKL = 3.80338
  EmaxL = 0.469779
  T16 = 1.06147
  T64 = 0.05
  T65 = 6.3
  T67 = 1.54865
  AlphOHgam =  0.111241
  k14a =  0.0000244437
  HApMRT = 3.60609
  koutL = 0.00293273
  OsteoEffectGam = 0.173833
  opgPTH50 = 3.85
  RX2Kout0 = 0.693
  E0rx2Kout = 0.125
  EmaxPTHRX2x = 5.0
  E0crebKin = 0.5
  EmaxPTHcreb = 3.39745
  crebKout = 0.00279513
  bcl2Kout = 0.693
  ScaEffGam = 0.9
  PhosEff0 = 1.52493
  PhosEff50 = 1.3021
  PhosEffGam = 8.25229
  PO4inhPTHgam = 0
  T69 = 0.10
  Reabs50 = 1.57322
  T7 = 2.0
  T9 = 90.0
  T70 = 0.01
  T71 = 0.03
  T33 = 0.003
  T34 = 0.037
  T35 = 90.0
  CaPOgam = 1.0
  T46 = 1.142
  T52 =  0.365
  OralPhos = 10.5/24
  F12 = 0.7
  T49 = 51.8
  T55 = 0.019268
  PicOBgamkb = 2.92375
  MultPicOBkb = 3.11842
  FracPic0kb = 0.764028
  E0RUNX2kbEffFACT = 1.01
  RUNkbGAM = 3.81644
  T43 = 1.03856
  T45 = 3.85
  T84 = 0.25
  RUNkbMaxFact = 0.638114
  RUNX20 = 10.0
  Frackb = 0.313186
  T81 = 0.75
  T87 = 0.0495
  T28 = 0.9
  OralCa = 24.055/24
  T310 = 0.105929
  T77 = 0.909359
  T80 = 4.0
  CtriolPTgam = 12.5033
  CtriolMax = 4.1029
  CtriolMin = 0.9
  PTout = 0.0001604
  kout = 100/14
  T58 = 6249.09
  T59 = 11.7387
  T61 = 96.25
  IPTHint = 0
  IPTHinf = 0
  Pic0 = 0.228142
  LsurvOCCgam = 3.0923
  EmaxLpth = 1.30721
  kO = 15.8885
  kb = 0.000605516
  LsurvOCgam =3.09023
  FracOBfast = 0.797629
  TERIKA = 10.4
  TERIVC = 94.4
  TERIVD  = 7.0
  TERICL = 62.2
  TERIF = 1.0
  PKKA = 0
  PKVC = 10.0
  PKQ1 = 0
  PKQ2 = 0
  PKVP1 = 1.0
  PKVP2 = 1.0
  PKCL = 0
  PKVMAX=0
  PKKM=1.0
  koutBMDls = 0.000397
  koutBMDlsDEN = 0.000145
  koutBMDfnBAS = 0.000005651
  gamOB = 0.0793
  gamOCls = 0.14
  gamOClsDEN = 0.0679
  gamOCfnBAS = 0.3101
  ETHN=0
  BMI=28.0
  kdenosl = 2.0e-06
  E2scalePicB1 = 0.0000116832
  # Newly created not sure what they mean
  RNK_0 =  10.0
  L_0 = 0.4
  O_0 = 4.0


  #= Denosumab PK parameters from...
  //  ##  M. Peterson B. Stouch D. Chen S. Baughman D. Holloway P. Bekker and S. Martin.
  //  ##  A population pk/pd model describes the rapid profound and sustained suppression of urinary
  //  ##  n-telopeptide following administration of amg 162 a fully human monoclonal antibody against
  //  ##  rankl to healthy postmenopausal women.
  //  ## The AAPS Journal 24(6 Abstract W4340) 2004.
  //
  //  ## DENK(1112) = 0.0141
  //  ## DENK(1211) = 0.00798
  =#
  DENVMAX = 3110.0
  DENKM = 188.0
  DENVC = 2340.0
  DENVP = 1324.0  # = Q/K(1211)
  DENCL = 2.75
  DENQ = 18.67 # = K(1211) * VC
  DENKA = 0.00592
  DENF = 0.729

  # This code block feeds in to control estrogen decline during menopause transition
  ESTON = 0.0
  koutEST=0.05776227
  menoDUR = 8736*1.66 #as.hour(as.year(1.66))
  ageGAM = -2.3
  age50 = 0.64
  ageENTER = 8736*41 #as.hour(as.year(41))
  ageDONE = 8736*51 #as.hour(as.year(51))
  tgfbGAM = 0.0374
  tgfbactGAM = 0.045273
  robGAM = 0.16
  obGAM = 0.000012
  maxTmESTkid = 0.923737

  # To use to get GFR to decline with time
  GFRtau=10.0 ## years over which GFR declines
  GFRdelta=0 ## ml/min

#$CMTN DENSC, TERISC

#$INIT
  DENSC= u[1]
  PTH = u[2]
  S = u[3]
  PTmax = u[4]
  B = u[5]
  P = u[6]
  ECCPhos = u[7]
  T = u[8]
  R = u[9]
  HAp = u[9]
  PhosGut = u[10]
  IntraPO = u[11]
  OC = u[12]
  ROB1 = u[13]
  L = u[14]
  RNK = u[15]
  O = u[16]
  Q = u[17]
  Qbone = u[18]
  RX2 = u[19]
  CREB = u[20]
  BCL2 = u[21]
  TERISC = u[22]
  TERICENT= u[23]
  PKGUT= u[24]
  PKCENT= u[25]
  PKPER1 = u[26]
  PKPER2 = u[27]
  DENCENT= u[28]
  DENPER = u[29]
  UCA= u[30]
  TGFB= u[31]
  TGFBact= u[32]
  OBfast= u[33]
  OBslow= u[34]
  M= u[35]
  N= u[36]
  AOH= u[37]
  EST  = u[38]
  BMDls = u[39]
  BMDlsDEN = u[40]
  BMDfn = u[41]
  GFR = u[42]
  #was in main

 # Initial conditions to build the rest
  DENSC_0= 0.0
  PTH_0 = 53.90
  S_0 = 0.5
  PTmax_0 = 1.0
  B_0 = 1260.0
  P_0 = 32.90
  ECCPhos_0 = 16.8
  T_0 = 1.58471
  R_0 = 0.50
  HAp_0 = 1.00
  PhosGut = 0.839
  IntraPO_0 = 3226.0
  OC_0 = 0.00115398
  ROB1_0 = 0.00104122
  L_0 = 0.4
  RNK_0 = 10.0
  O_0 = 4.0
  Q_0 = 100.0
  Qbone_0 = 24900.0
  RX2_0 = 10.0
  CREB_0 = 10.0
  BCL2_0 = 100.0
  TERISC_0 = 0.0
  TERICENT_0 = 0.0
  PKGUT_0 =0.0
  PKCENT_0 =0.0
  PKPER1_0 = 0.0
  PKPER2_0 = 0.0
  DENCENT_0 =0.0
  DENPER_0 = 0.0
  UCA_0=0.0
  TGFB_0=0.0
  TGFBact_0=0.0
  OBfast_0=0.0
  OBslow_0=0.0
  M_0=0.0
  N_0=0.0
  AOH_0=126.0
  EST_0  = 1.0
  BMDls_0 = 1.0
  BMDlsDEN_0 = 1.0
  BMDfn_0 = 1.0
  GFR_0 = 100.0/16.667

  TGFB_0 = Pic0*1000.0
  TGFBact_0 = Pic0
  OBfast_0 = OBtot0*FracOBfast
  OBslow_0 = OBtot0*(1-FracOBfast)
  OB0 = (OBfast_0 + OBslow_0)
  M_0 = k3*RNK_0*L_0/k4
  N_0 = k1*O_0*L_0/k2
  AOH_0 = B_0/10.0
  F_DENSC = DENF*1E6
  F_TERISC = TERIF

#[ODE]
  #=**************************************************************
  //   Calcium / bone model algebraic relationships and
  //   differential equations
  **************************************************************=#

  PhosEffect = 0
  J48 = 0
  J27 = 0
  RUNX2 = 0
  kinEST = 0

  # parameters derived from SS initial conditions #
  T13 = (CaDay/24)/Q_0

  CaConc0 = (P_0/V1)

  T15 = CaDay/(CaConc0*V1*24)

  J14OC50= exp(log((J14OCmax*pow(OC_0,J14OCgam)/T13) - pow(OC_0,J14OCgam))/J14OCgam)

  OCeqn = (J14OCmax*pow(Osteoclast,J14OCgam))/(pow(Osteoclast,J14OCgam) + pow(J14OC50,J14OCgam))

  kinRNK = (koutRNK*RNK_0 + k3*RNK_0*L_0 - k4*M_0) / pow(TGFBact_0,kinRNKgam)

  MOCratio = M/Osteoclast

  MOCratio0 = M_0/OC_0

  MOCratioEff = pow((MOCratio/MOCratio0), MOCratioGam)

  J14OCdepend = OCeqn*Q_0*FracJ14*MOCratioEff

  J14 = T13*Q_0*(1-FracJ14) + J14OCdepend
  # 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. #
  J41 = 0.464*J14

  bigDb = kb*OB0*Pic0/ROB1_0

  kinTGF = koutTGF0*TGFB_0

  koutTGF = koutTGF0*(pow((TGFB/TGFB_0),koutTGFGam))

  koutTGFact = koutTGF0*1000

  koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC_0), OCtgfGAM))

  E0PicROB = FracPicROB*Pic0

  EC50PicROBparen= (EmaxPicROB*pow(TGFBact_0,PicROBgam) / (Pic0 - E0PicROB)) - pow(TGFBact_0,PicROBgam)

  EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam)

  Dr = kb*OB0/Pic0

  PicROB = E0PicROB + EmaxPicROB*pow(TGFBact,PicROBgam)/(pow(TGFBact,PicROBgam) + pow(EC50PicROB,PicROBgam))

  ROBin = Dr*PicROB

  E0PicOB = FracPicOB*Pic0

  EC50PicOBparen = (EmaxPicOB*pow(TGFBact_0,PicOBgam)/(Pic0 - E0PicOB)) - pow(TGFBact_0,PicOBgam)

  EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam)

  PicOB = E0PicOB + EmaxPicOB*pow(TGFBact,PicOBgam) / (pow(TGFBact,PicOBgam) + pow(EC50PicOB,PicOBgam))

  KPT =1*(bigDb/PicOB)

  EC50MeffOC = exp(log(pow(M_0, kinOCgam)*EmaxMeffOC/(1-E0Meff) - pow(M_0, kinOCgam))/kinOCgam)

  MeffOC = E0Meff + (EmaxMeffOC * pow(M, kinOCgam)/(pow(M, kinOCgam) + pow(EC50MeffOC,kinOCgam)))

  PicOCkin = Pic0

  kinOC2 = Da*PicOCkin*MeffOC*OC_0

  E0PicOC = FracPicOC*Pic0

  EC50PicOCparen = (EmaxPicOC*pow(TGFBact_0, PicOCgam)/(Pic0 - E0PicOC)) - pow(TGFBact_0, PicOCgam)

  EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam)

  PicOC = E0PicOC + ((EmaxPicOC*pow(TGFBact, PicOCgam))/(pow(TGFBact, PicOCgam) + pow(EC50PicOC, PicOCgam)))

  PiL0 = (k3/k4)*L_0

  PiL = M/10

  EC50survInPar = (E0RANKL - EmaxL)*(pow(PiL0, LsurvOCgam)/(E0RANKL - 1)) - pow(PiL0, LsurvOCgam)

  EC50surv = exp(log(EC50survInPar)/LsurvOCgam)

  LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(pow(PiL, LsurvOCgam)/(pow(PiL, LsurvOCgam) + pow(EC50surv, LsurvOCgam)))

  KLSoc = Da*PicOC*LsurvOC

  PTHconc0 = (PTH_0/V1)

  T66 = (pow(T67, AlphOHgam) + pow(PTHconc0, AlphOHgam) )/pow(PTHconc0, AlphOHgam)

  k15a = k14a*Qbone_0/Q_0

  J14a = k14a*Qbone

  J15a = k15a*Q

  # Hydroxy-apatite #
  kLShap = 1/HApMRT

  kHApIn = kLShap/OB0

  # Calcium flux from plasma into bone #
  J15 = (T15*P*(1-FracJ15) + T15*P*FracJ15*HAp)

  # 0.464, reported as the molar ratio of P / Ca in hydroxyapatite. #
  J42 = 0.464*J15

  Osteoblast = OBfast + OBslow

  kinLbase = koutL*L_0

  OsteoEffect = pow((Osteoblast/OB0), OsteoEffectGam)

  PTH50 = EmaxLpth*PTHconc0 - PTHconc0

  LpthEff = EmaxLpth*(PTHconc) / ((PTH50*OsteoEffect) + (PTHconc))

  kinL = kinLbase*(OsteoEffect)*LpthEff

  pObase = kO*O_0

  pO = pObase*(D/ROB1_0)*((PTHconc+(opgPTH50*(D/ROB1_0)))/(2*PTHconc))

  RX2Kin = RX2Kout0*RX2_0

  EC50PTHRX2x = ((EmaxPTHRX2x*PTHconc0)/(RX2Kout0 - E0rx2Kout)) - PTHconc0

  RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x)


  #*****************************************************#
  # START CREB-RELATED EQUATIONS                        #
  #*****************************************************#
  EC50PTHcreb = ((EmaxPTHcreb*PTHconc0)/(1-E0crebKin)) -  PTHconc0

  crebKin0= crebKout*CREB_0

  crebKin = crebKin0* (E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb))

  bcl2Kin = RX2*CREB*bcl2Kout
  #*****************************************************#

  #*****************************************************#
  # START PHOS-RELATED EQUATIONS                        #
  #*****************************************************#

  # C2 is extracellular phosphate concentration #

  C2 = (ECCPhos/V1)

  PO4inhPTH = pow((C2/1.2),PO4inhPTHgam)

  PhosEffTop = (PhosEff0 - 1)*(pow(1.2, PhosEffGam) + pow(PhosEff50, PhosEffGam) )

  PhosEffBot =PhosEff0 * pow(1.2, PhosEffGam)

  PhosEffMax =  PhosEffTop / PhosEffBot

  PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * pow(C2, PhosEffGam) /(pow(C2, PhosEffGam)  + pow(PhosEff50, PhosEffGam)))

  if (C2 > 1.2)
    PhosEffect = PhosEff
  else
    PhosEffect = 1
  end

  T68 = T66*pow(PTHconc, AlphOHgam)/(pow(T67, AlphOHgam)*PO4inhPTH+pow(PTHconc, AlphOHgam))

  SE = T65*T68*PhosEffect

  # Equations relating to calcitriol-dependent calcium absorption #
  C8 = (B/V1)
  T36 = T33 + (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)))
  T37 = T34 - (T34-T33)*(pow(C8,CaPOgam)/(pow(T35,CaPOgam)+ pow(C8,CaPOgam)))

  #======================================
  //   RENAL CALCIUM HANDLING
  //   ======================================
  //   Calcium filtration rate in kidney
  //   We assume that 50% of filtered calcium is reabsorbed in a PTH-independent manner
  //   ... and 50% is reabsorbed in a PTH-dependent manner
  //   Fraction unbound assumed to be 0.6
  //======================================#

  CaConc = (P/V1)
  CaFilt = 0.6*0.5*GFR*CaConc

  # Maximum calcium reabsorption in the kidney - PTH sensitive#
  mtmEST = (1-maxTmESTkid)/(1-0.1)
  tmEST = 1 - mtmEST + mtmEST*EST

  ReabsMax = tmEST * (0.3*GFR*CaConc0 - 0.149997)*(Reabs50 + CaConc0) / CaConc0

  # Effect of PTH on calcium reabsorption #
  T17 = PTHconc0*T16 - PTHconc0

  ReabsPTHeff = (T16*PTHconc)/(PTHconc + T17)

  # PTH-sensitive calcium reabsorption in kidney #
  # Reabs50 = 1.573 = H(4-u)-delta #
  C1 = (P/V1)
  CaReabsActive =  (ReabsMax*C1/(Reabs50 + C1))*ReabsPTHeff

  T20 = CaFilt - CaReabsActive

  T10 = T7*C8/(C8+T9)

  # Temporary calcium excretion rate #
  J27a = (2-T10)*T20

  # J27 will be the flux of calcium out of the plasma via the kidney #
  if (J27a<0)
    J27 = 0
  else
    J27 = J27a
  end

  ScaEff = pow( (CaConc0/CaConc), ScaEffGam)

  T72 = 90 * ScaEff

  T73 = T71 * (C8 - T72)

  T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))

  T75 = T70 * (0.85 * (1 + T74) + 0.15)

  T76 = T70 * (0.85 * (1 - T74) + 0.15)

  # phosphate renal excretion #
  T47 = T46*0.88*GFR

  J48a = 0.88*GFR*C2 - T47

  if (J48a < 0)
    J48 = 0
  else
    J48 = J48a
  end

  # phosphate oral absorption #
  J53 = T52*PhosGut

  J54 = T49*C2

  J56 = T55*IntraPO

  # Parameters describing TGF-beta effects on Osteoblast and clast differentiation and apoptosis #
  E0PicOBkb = MultPicOBkb*Pic0

  EmaxPicOBkb = FracPic0kb*Pic0

  EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*pow(TGFBact_0,PicOBgamkb)) / (E0PicOBkb - Pic0)  - pow(TGFBact_0,PicOBgamkb)

  EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb)

  PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*pow(TGFBact,PicOBgamkb) / (pow(TGFBact,PicOBgamkb) + pow(EC50PicOBkb,PicOBgamkb))

  # Estrogen effect that propogates through to OB apoptosis#
  PicOBkbEff = (PicOBkb/Pic0)*(1/(pow(EST,E2scalePicB1)))

  # Parameters describing osteoblast apoptosis as affected by PTH (continuous vs intermitent) #
  E0RUNX2kbEff= E0RUNX2kbEffFACT*kb

  if (BCL2 > 105)
    RUNX2 = BCL2 - 90
  else
    RUNX2 = 10
  end

  RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact

  INparen = (RUNkbMax * pow(RUNX20,RUNkbGAM)) / (E0RUNX2kbEff - kb) - pow(RUNX20,RUNkbGAM)

  RUNkb50 = exp(log(INparen)/RUNkbGAM)

  RUNX2kbPrimeEff = RUNkbMax*pow(RUNX2,RUNkbGAM) / (pow(RUNX2,RUNkbGAM) + pow(RUNkb50,RUNkbGAM))

  kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff

  kbslow = kbprime*Frackb

  kbfast = (kb*OB0 + kbslow*OBfast_0 - kbslow*OB0) / OBfast_0

  Frackb2 = kbfast/kbprime

  #*********************************************************#
  # Equations relating to calcium movement to/from the gut #
  #*******************************************************#
  T29 = (T28*T_0 - T310*T_0)/T310

  T31 = T28*T/(T+T29)

  # R is calcitriol-dependent gut Ca2+ absorption #
  T83 = R/0.5

  # J40 = calcium flux from gut to plasma #
  J40 = T31*T*T83/(T + T81) + T87*T

  # T85 relates to extent of absorption of orally-administered dose #

  T85Rpart = pow(R, T80)/(pow(R,T80) + pow(T81,T80))
  T85 = T77*T85Rpart

  #= ****************************#
  # Calcitriol equations     #
  #*************************# =#

  Calcitriol0 = (B_0/V1)

  INparenCtriol =((CtriolMax - CtriolMin) * pow(Calcitriol0, CtriolPTgam)) / (CtriolMax - 1)- pow(Calcitriol0,CtriolPTgam)

  Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam)

  CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * pow(C8, CtriolPTgam) / (pow(C8, CtriolPTgam) + pow(Ctriol50, CtriolPTgam))

  PTin = PTout * CtriolPTeff

  # S is the PTH gland pool cmt
  FCTD = (S / 0.5) * PTmax

  INparenCa =(T58 - T61) * pow(CaConc0, T59) / (T58 - 385) - pow(CaConc0, T59)
  T60 = exp(log(INparenCa) / T59)
  T63 =  T58 - (T58 - T61) * pow((CaConc), T59) / (pow((CaConc), T59) + pow(T60, T59))

  # Total PTH input (production) rate #
  SPTH = T63*FCTD

  # Teriparatide pk #
  TERIPKIN = TERISC*TERICL/TERIVC

  # Plasma PTH (pmol)
  #   SPTH = PTH input rate
  #   kout = PTH first order elimination rate
  #   TERIPKIN = first order rate from tpar subq dosing into plasma

  du[2] = dxdt_PTH = SPTH - kout*PTH + TERIPKIN

  # PTH gland pool
  du[3] = dxdt_S = (1 - S) * T76 - (S* T75)

  # PT gland max capacity #
  du[4] = dxdt_PTmax = PTin - PTout * PTmax

  du[] = dxdt_B = AOH - T69 * B

  # 1-alpha-hydroxylase (AOH) cmt
  du[] = dxdt_AOH = SE - T64*AOH

  #= J14 = nu(12-4) calcium flux from bone into plasma
  //   J15 = nu(4-12) calcium flux from plasma into bone
  //   J27 = nu(4-u)  calcium flux from plasma to urine
  //   J40 = nu(1-4)  calcium flux from gut to plasma
  // =#
  du[] = dxdt_P = J14 - J15- J27 + J40

  # Extracelluar phosphate (mmol) #
  # d/dt(extracellular phosphate) = J41 -  J42 - J48 + J53 - J54 + J56
  #=   d/dt(intracellular phosphate) = J54  -  J56

  //      The exchange fluxes of PO4 between ECF and bone (J41 and J42)
  //      set same as respective Ca fluxes but multiplied by the stoichiometric
  //      factor of 0.464, reported as the molar ratio of P / Ca in hydroxyapatite.
  //
  //      d/dt(dietary phosphate) = OralPhos*F12 -J53
  =#

  du[] = dxdt_PhosGut = OralPhos *F12 - J53
  du[] = dxdt_IntraPO = J54 - J56
  du[] = dxdt_ECCPhos = J41  - J42 - J48 + J53 - J54 + J56

  # Oral calcium #
  # CMT: T  UNITS: mmol #
  # J40 --> flux from gut to plasma #
  # F11 == T85 by definition #
  du[] = dxdt_T = OralCa*T85 - J40

  # Calcitriol-dependent Ca2+ absorption #
  # CMT: R #
  du[] = dxdt_R = T36*(1 - R) - T37*R

  # Hydroxyapatite #
  du[] = dxdt_HAp = kHApIn*Osteoblast - kLShap*HAp

  #*****************#
  # Estrogen piece#
  #*****************#

  AGE = ageENTER + SOLVERTIME

  ageONSET = ageDONE-menoDUR

  if (AGE < ageONSET)
    kinEST = koutEST * pow((AGE/ageENTER),ageGAM)
  end

  if (AGE >= ageONSET)
    kinEST = koutEST * pow((AGE/ageENTER),ageGAM) * (1 - age50 * (pow((AGE-ageONSET),2)/(pow((menoDUR/2),2) + pow((AGE-ageONSET),2))))
  end

  du[] = dxdt_EST = (kinEST - koutEST * EST)*ESTON

  # Osteoblasts were considered to exist as two populations: fast and slow.
  #   Fast and slow refer to removal rates (kbfast, kbslow) input assumed to
  #   be the same for each.  Total osteoblasts = OBfast + OBslow #

  # Estrogen effect added: E2dosePicB1 --> PicOBkbEff --> kbprime --> kbslow, kbfast and Frackb2#
  du[] = dxdt_OBfast = (bigDb/PicOB)*D*FracOBfast*Frackb2  - kbfast*OBfast

  du[] = dxdt_OBslow = (bigDb/PicOB)*D*(1-FracOBfast)*Frackb - kbslow*OBslow

  # OC: Active Osteoclasts #
  du[] = dxdt_OC = kinOC2 - KLSoc*OC

  # D = ROB1 Responding Osteoblasts #
  du[] = dxdt_ROB1 = ROBin * pow(1/EST,robGAM) - KPT*ROB1

  # Latent TGF-beta pool, production dependent on osteoblast function #
  du[] = dxdt_TGFB = kinTGF*(pow((Osteoblast/OB0),OBtgfGAM)) * pow(1/EST,tgfbGAM) - koutTGFeqn * pow(EST,tgfbactGAM)

  # active TGF-beta pool, production dependent on osteoclast function
  #   koutTGFeqn = koutTGF*TGFB*(pow((Osteoclast/OC_0), OCtgfGAM))

  du[] = dxdt_TGFBact = koutTGFeqn * pow(EST,tgfbactGAM) - koutTGFact*TGFBact

  #********************************#
  # BMD BONE MINERAL DENSITY #
  #********************************#
  # Solve for kinBMD based on OC=OC_0, OB=OB0, and BMD=BMD0 #

  #Lumbar spine#
  kinBMDls =  koutBMDls*BMDls_0

  du[] = dxdt_BMDls = kinBMDls * pow(OB/OB0,gamOB) - koutBMDls * pow(OC/OC_0,gamOCls) * BMDls

  #Lumbar spine with DENOSUMAB#
  kinBMDlsDEN =  koutBMDlsDEN*BMDlsDEN_0

  du[] = dxdt_BMDlsDEN = kinBMDlsDEN * pow(OB/OB0,gamOB) - koutBMDlsDEN * pow(OC/OC_0,gamOClsDEN) * BMDlsDEN

  #Femoral neck#
  gamOCfn = gamOCfnBAS
  koutBMDfn = koutBMDfnBAS

  kinBMDfn =  koutBMDfn*BMDfn_0

  du[] = dxdt_BMDfn = kinBMDfn * pow(OB/OB0,gamOB) - koutBMDfn * pow(OC/OC_0,gamOCfn) * BMDfn

  #**************************************************#
  #Treatment effects#
  #**************************************************#

  DENMOL =  (DENCENT/DENVC/150000)*1000*1000

  # L: RANK-L #
  du[] = dxdt_L = kinL- koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M -  kdenosl*DENMOL*L

  # RNK: RANK #
  du[] = dxdt_RNK = kinRNK*pow(TGFBact,kinRNKgam) - koutRNK*RNK - k3*RNK*L  + k4*M

  # M: RANK - RANK-L complex #
  du[] = dxdt_M = k3*RNK*L - k4*M

  # N: OPG - RANK-L complex #
  du[] = dxdt_N = k1*O*L - k2*N

  # O: OPG #
  du[] = dxdt_O = pO - k1*O*L + k2*N - kO*O

  #= d/dt(Q) = J15-J14+ J14a-J15a  Q=exchangeable bone Ca
  //   d/dt(Qbone) = -J14a + J15a    Qbone=non(immediately)-exchangeable bone Ca
  //   ~ 99% of the total Ca stored in bone approximately 100 mmol of 25000 - 30000 mmol
  //   of total skeletal Ca considered immediately exchangeable with plasma Ca
  =#
  du[] = dxdt_Q = J15 - J14 + J14a - J15a

  du[] = dxdt_Qbone = J15a - J14a

  #=
  //      Initial conditions set empirically:
  //      both RX2 and CREB started at 10, BCL2 was the product (10*10=100).
  //
  //      BCL2 assumed half-life of 1 hour: rate const set to 0.693 (bcl2Kout)
  //
  //      BCL2 affected osteoblast survival:
  //      decreases elimination rate constant for OB (kbprime) =#

  du[] = dxdt_RX2 = RX2Kin - RX2Kout*RX2

  du[] = dxdt_CREB = crebKin - crebKout*CREB

  du[] = dxdt_BCL2 = bcl2Kout*CREB*RX2 - bcl2Kout*BCL2

  # Teriparatide PK info #
  # TERIPKIN = TERISC * TERICL/TERIVC #
  du[] = dxdt_TERISC = -TERIPKIN
  du[] = dxdt_TERICENT  = TERIPKIN - TERICENT*TERIKA

  # GENERAL PK COMPARTMENT #

  # NONLINEAR PIECE #
  PKCP = (PKCENT/PKVC)
  PKCLNL = PKVMAX/(PKKM+PKCP)

  # PKGUT #
  du[] = dxdt_PKGUT = -PKKA*PKGUT

  # PKCENT #
  du[] = dxdt_PKCENT = PKKA*PKGUT + PKQ1*PKPER1/PKVP1 + PKQ2*PKPER2/PKVP2 -
    (PKQ1+PKQ2+PKCL+PKCLNL)*PKCENT/PKVC

  # PKPER1 #
  du[] = dxdt_PKPER1 = PKQ1*PKCENT/PKVC - PKQ1*PKPER1/PKVP1
  # PKPER2 #
  du[] = dxdt_PKPER2 = PKQ2*PKCENT/PKVC - PKQ2*PKPER2/PKVP2
  # DENOSUMAB PK #
  DENCP = (DENCENT/DENVC)
  DENCLNL =  (DENVMAX/(DENKM+DENCP))
  du[1] = dxdt_DENSC =  -DENKA*DENSC
  du[2] = dxdt_DENCENT = DENKA*DENSC + DENQ*DENPER/DENVP -
    (DENQ+DENCL+DENCLNL)*DENCENT/DENVC
  du[] = dxdt_DENPER = DENQ*DENCENT/DENVC - DENQ*DENPER/DENVP
  # DENOSUMAB PK #
  GFRend = GFR_0 - GFRdelta/16.667
  GFRtau_ = GFRtau*8766
  kGFR = -log(GFRend/GFR_0)/GFRtau_
  du[] = dxdt_GFR = -kGFR*GFR
  # Collects calcium in urine - cumulative rate #
  # This is a differential equation #
  du[] = dxdt_UCA = J27
end

uo = (0.0, 53.90, 0.5, 1.00, 1260.0, 32.90, 16.8, 1.58471, 0.50, 1.00, 0.839, 3226.0, 0.00115398, 0.00104122, 0.4, 10.0, 4.0,
100.0, 24900.0, 10.0, 10.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 126.0, 1.0, 1.0, 1.0, 1.0, 100.0/16.667)

prob = ODEProblem(modelODE!, u0, (0.0,10.0))
sol = solve(prob)

# [TABLE]
# capture PTHpm = (PTH/V1)   # PTH conc in pM #
# capture PTHc = (PTH/V1*9.4)      # PTH conc in pg/mL #
# capture CAchange = (100*P/P_0)    /# Calcium (total, serum), %change from baseline #
# capture OBchange = 100*OB/OB0  # bone-specific alkaline phosphatase, %change from baseline #
# capture OCchange = 100*OC/OC_0  # serum CTx, %change from baseline #
# capture CaC =  P/V1    # Calcium (total, serum) in mM #
# capture BMDlsDENchange = (BMDlsDEN-1)*100
# capture OBtot = OBfast + OBslow
# capture OCfrac = OC/OC_0
# capture OBfrac = (OBfast+OBslow)/(OBfast_0 + OBslow_0)
#
# $CAPTURE DENMOL T43 DENCP
