import Rasterizer
import GameLogic

cont = GameLogic.getCurrentController()
cont.owner.suspendDynamics()
cont.owner.restoreDynamics()

