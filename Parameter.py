import Conversion as Conv

SimNum = 1
SimMin = 0
SimMax = 9

numin = 0.4
numax = 0.9
GammaMin = 0.03
GammaMax = 0.8

#nu = numin + SimNum*(numax-numin)/MaxSim
Gamma = GammaMin + (SimNum-SimMax)*(GammaMax-GammaMin)/(SimMax-SimMin)

P = Conv.AnalyticToSimul(Gamma = 0.5, nu = 0.95, l = 5., epsilon= 0.01,writting= False,ParticleType='Hexagon')
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
SizeX=30
SizeY=30
NumberOfParticle=200
Expansion = True
Output=True
#  ____                                              _
# |  _ \    __ _   _ __    __ _   _ __ ___     ___  | |_    ___   _ __   ___
# | |_) |  / _` | | '__|  / _` | | '_ ` _ \   / _ \ | __|  / _ \ | '__| / __|
# |  __/  | (_| | | |    | (_| | | | | | | | |  __/ | |_  |  __/ | |    \__ \
# |_|      \__,_| |_|     \__,_| |_| |_| |_|  \___|  \__|  \___| |_|    |___/

TimeStepTot=1000000
StatTime=TimeStepTot//100
BetaInitial=0
BetaFinal=1.6*10**2
Seed=None
DEG=0.0125
