import Conversion as Conv

SimNum = 1

numin = 0.4
numax = 0.9
MaxSim = 5.

nu = numin + SimNum*(numax-numin)/MaxSim

P = Conv.AnalyticToSimul(Gamma = 0.5, nu = nu, l = 5., epsilon= 0.01,writting= False,ParticleType='Hexagon')
#  ____                  _
# / ___|   _   _   ___  | |_    ___   _ __ ___
# \___ \  | | | | / __| | __|  / _ \ | '_ ` _ \
#  ___) | | |_| | \__ \ | |_  |  __/ | | | | | |
# |____/   \__, | |___/  \__|  \___| |_| |_| |_|
#          |___/
ParticleType = P.ParticleType#'Hexagon'
Kmain=P.k#1.
Kcoupling=P.kc#0.04629
Eps=P.epsilon#0.01
KVOL=P.kA#0.288675
#----------------
J=P.J#2.774*10**(-5)
#----------------
SizeX=35
SizeY=35
<<<<<<< HEAD
NumberOfParticle=100
=======
NumberOfParticle=300
>>>>>>> daf6d0a2156fb9ed1817149e5fea07fbf208c177
Expansion = True
Output=True
#  ____                                              _
# |  _ \    __ _   _ __    __ _   _ __ ___     ___  | |_    ___   _ __   ___
# | |_) |  / _` | | '__|  / _` | | '_ ` _ \   / _ \ | __|  / _ \ | '__| / __|
# |  __/  | (_| | | |    | (_| | | | | | | | |  __/ | |_  |  __/ | |    \__ \
# |_|      \__,_| |_|     \__,_| |_| |_| |_|  \___|  \__|  \___| |_|    |___/

<<<<<<< HEAD
TimeStepTot=1000
StatTime=100#TimeStepTot//100
=======
TimeStepTot=500000
StatTime=TimeStepTot//100
>>>>>>> daf6d0a2156fb9ed1817149e5fea07fbf208c177
BetaInitial=0
BetaFinal=1.6*10**2
Seed=None
DEG=0.0125
