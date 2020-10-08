#  ____                  _
# / ___|   _   _   ___  | |_    ___   _ __ ___
# \___ \  | | | | / __| | __|  / _ \ | '_ ` _ \
#  ___) | | |_| | \__ \ | |_  |  __/ | | | | | |
# |____/   \__, | |___/  \__|  \___| |_| |_| |_|
#          |___/

Kmain=1.
Kcoupling=21
Eps=0.01
KVOL=0.
#----------------
J=0.0003
#----------------
SizeX=35
SizeY=35
NumberOfParticle=150

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
#  _____                           _                   _
# |_   _|   ___    _ __     ___   | |   ___     __ _  (_)   ___
#   | |    / _ \  | '_ \   / _ \  | |  / _ \   / _` | | |  / _ \
#   | |   | (_) | | |_) | | (_) | | | | (_) | | (_| | | | |  __/
#   |_|    \___/  | .__/   \___/  |_|  \___/   \__, | |_|  \___|
#                 |_|                          |___/
TopologieDown = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieUp = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
#TopologieDown = [(1,0),(-1,0),(0,1)]
#TopologieUp = [(1,0),(-1,0),(0,-1)]
