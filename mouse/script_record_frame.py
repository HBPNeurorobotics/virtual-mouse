
import Rasterizer
# ffmpeg -b 64000 -r 20 -sameq -i game_%04d.png -vcodec libx264 ../video.mp4

rec_i = -1

def record():
    global rec_i
    scr = Rasterizer.makeScreenshot("/tmp/game_"+str('%04d' % rec_i)+".jpg")
    rec_i+=1
    

