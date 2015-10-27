import sys
while len(sys.meta_path)>3:
    sys.meta_path.pop(0)


import os

filename_tmp = "script_importpaths.py"
exec(compile(open(filename_tmp).read(), filename_tmp, 'exec'))

filename_tmp = "script_fcts.py"
exec(compile(open(filename_tmp).read(), filename_tmp, 'exec'))

import numpy as np
import math
#~ import scipy.misc
import ctypes

import bpy
import bge
import GameLogic as g
import VideoTexture
scn = g.getCurrentScene()



firstTimeCsX = True
try:
    bpy.ExistentialVariable
    firstTimeCsX = False
except:
    bpy.ExistentialVariable = 42




# =================================================================================================================================================
#                                       High-level functions to access from controller.py



def getVisual(camera_name, max_dimensions):
    global scn
    global camera_array
    src1 = VideoTexture.ImageRender(scn, scn.objects[camera_name])
    #src1 = VideoTexture.ImageViewport(scn, scn.objects[camera_name])
    src1.capsize = max_dimensions[0],max_dimensions[1]
    imX,imY = src1.size
    image_np = np.copy(VideoTexture.imageToArray(src1) ).reshape((imY, imX, 4))
    image_np = image_np[::-1,:,:]
    return image_np


def getOlfactory(olfactory_object_name, receptor_names, diffusivity = 10.0):
    olfactory_outputs = np.zeros(len(receptor_names), dtype=np.float32)
    for objkeys in bpy.data.objects.keys():
        for irn,rn in enumerate(receptor_names):
            smellampl = None
            try:
                smellampl = bpy.data.objects[objkeys][rn]
            except:
                smellampl=None
            if smellampl!=None:
                olfactory_outputs[irn] += smellampl * np.exp(-(np.linalg.norm(scn.objects[objkeys].worldPosition-scn.objects[olfactory_object_name].worldPosition)**2.0)/diffusivity)
    return olfactory_outputs


# TODO: make it touch-dependent rather than distance-dependent
def getTaste(taste_object_name, receptor_names, distance_to_object):
    taste_outputs = np.zeros(len(receptor_names), dtype=np.float32)
    for objkeys in bpy.data.objects.keys():
        for irn,rn in enumerate(receptor_names):
            tasteampl = None
            try:
                tasteampl = bpy.data.objects[objkeys][rn]
            except:
                tasteampl=None
            if tasteampl!=None:
                if np.linalg.norm(scn.objects[objkeys].worldPosition-scn.objects[taste_object_name].worldPosition) <= distance_to_object:
                    taste_outputs[irn] = tasteampl
    return taste_outputs


velocities_old = {}
def getVestibular(vestibular_object_name = "obj_head"):
    global velocities_old
    av_ = scn.objects[vestibular_object_name].localAngularVelocity
    lv_ = scn.objects[vestibular_object_name].localLinearVelocity
    velocities_tmp = np.array([av_[0], av_[1], av_[2], lv_[0], lv_[1], lv_[2]])
    gravitational_acceleration = newref([0.0,0.0,-9.80], scn.objects[vestibular_object_name].worldOrientation.transposed())
    # if first step
    if vestibular_object_name not in velocities_old.keys():
        velocities_old[vestibular_object_name] = np.copy(velocities_tmp)
    accelerations_tmp = (velocities_tmp - velocities_old[vestibular_object_name]) * float(bpy.context.scene.game_settings.fps) 
    accelerations_tmp[3:] -= gravitational_acceleration
    return accelerations_tmp










controllable_muscles = []
object_rotation_sum     = {}


def setPositionServo(reference_object_name, attached_object_name, P, maxV = None, DoF_vector = None):
    muscle_type = "position_servo"
    global controllable_muscles
    if DoF_vector==None:
        bpy.context.scene.objects.active = bpy.data.objects[reference_object_name]
        vecc = np.array([1.0,0.0,0.0])
        vecY = rotate_around_vector(vec = np.array([0.0,1.0,0.0]), vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = vecY,                    angle = bpy.context.object.constraints["Rigid Body Joint"].axis_y)
        DoF_vector = vecc
    controllable_muscles.append( {"type": muscle_type, "ref_name": reference_object_name, "att_name": attached_object_name, "P": P, "maxV": maxV, "orientation":  np.array(DoF_vector), "limit_angles": [bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x, bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x], "vel": 0.0, "rot": -10000.0, "rotSpindle": -10000.0  } )
    return len(controllable_muscles)-1

def setVelocityServo(reference_object_name, attached_object_name, maxV, DoF_vector = None):
    muscle_type = "velocity_servo"
    global controllable_muscles
    if DoF_vector==None:
        bpy.context.scene.objects.active = bpy.data.objects[reference_object_name]
        vecc = np.array([1.0,0.0,0.0])
        vecY = rotate_around_vector(vec = np.array([0.0,1.0,0.0]), vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = vecY,                    angle = bpy.context.object.constraints["Rigid Body Joint"].axis_y)
        DoF_vector = vecc
    controllable_muscles.append( {"type": muscle_type, "ref_name": reference_object_name, "att_name": attached_object_name, "maxV": maxV, "orientation":  np.array(DoF_vector), "limit_angles": [bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x, bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x], "vel": 0.0, "rot": -10000.0  } )
    return len(controllable_muscles)-1

def setTorqueMusclePair(reference_object_name, attached_object_name, maxF, DoF_vector = None):
    muscle_type = "torque-based"
    global controllable_muscles
    if DoF_vector==None:
        bpy.context.scene.objects.active = bpy.data.objects[reference_object_name]
        vecc = np.array([1.0,0.0,0.0])
        vecY = rotate_around_vector(vec = np.array([0.0,1.0,0.0]), vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = np.array([0.0,0.0,1.0]), angle = bpy.context.object.constraints["Rigid Body Joint"].axis_z)
        vecc = rotate_around_vector(vec = vecc,                    vec_ref = vecY,                    angle = bpy.context.object.constraints["Rigid Body Joint"].axis_y)
        DoF_vector = vecc
    controllable_muscles.append( {"type": muscle_type, "ref_name": reference_object_name, "att_name": attached_object_name, "isflex": 1, "maxF": maxF, "orientation":  np.array(DoF_vector), "limit_angles": [bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x, bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x], "vel": 0.0, "rot": -20000.0  } )
    controllable_muscles.append( {"type": muscle_type, "ref_name": reference_object_name, "att_name": attached_object_name, "isflex": 0, "maxF": maxF, "orientation": -np.array(DoF_vector), "limit_angles": [bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x, bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x], "vel": 0.0, "rot": -20000.0  } )
    return [len(controllable_muscles)-2, len(controllable_muscles)-1]

def setSingleLinearMuscle(reference_object_name, reference_object_application_point, attached_object_name, attached_object_application_point, maxF):
    muscle_type = "single-linear"
    global controllable_muscles
    controllable_muscles.append( {"type": muscle_type, "ref_name": reference_object_name, "ref_point": reference_object_application_point, "att_name": attached_object_name, "att_point": attached_object_application_point, "maxF": maxF, "vel": 0.0, "muscledist": -20000.0  } )
    return len(controllable_muscles)-1





def getMuscleSpindle(control_id):
    global scn
    if controllable_muscles[control_id]["type"]=="position_servo":
        a_ = 1.0/(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0])
        b_ = -a_ * controllable_muscles[control_id]["limit_angles"][0]
        rot_n = (a_*( scn.objects[controllable_muscles[control_id]["ref_name"]].getRigidBodyJointAngles(0) )  + b_) 
        vel_n = 0.0
        if controllable_muscles[control_id]["rotSpindle"]>-10000.0:
            vel_n = np.fabs(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0]) * (rot_n - controllable_muscles[control_id]["rotSpindle"])/dt_bl
        controllable_muscles[control_id]["rotSpindle"] = rot_n
        return [rot_n, vel_n, 0.0]
    if controllable_muscles[control_id]["type"]=="torque-based":
        a_ = 1.0/(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0])
        b_ = -a_ * controllable_muscles[control_id]["limit_angles"][0]
        rot_n = (a_*( scn.objects[controllable_muscles[control_id]["ref_name"]].getRigidBodyJointAngles(0) )  + b_) 
        if controllable_muscles[control_id]["isflex"]==1:
            rot_n = 1.0 - rot_n
        vel_n = 0.0
        if controllable_muscles[control_id]["rot"]>-10000.0:
            vel_n = np.fabs(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0]) * (rot_n - controllable_muscles[control_id]["rot"])/dt_bl
        controllable_muscles[control_id]["rot"] = rot_n
        return [rot_n, vel_n, 0.0]
    if controllable_muscles[control_id]["type"]=="single-linear":
        muscle_pos1_W = np.array(newref(controllable_muscles[control_id]["ref_point"], scn.objects[controllable_muscles[control_id]["ref_name"]].worldOrientation))
        muscle_pos2_W = np.array(newref(controllable_muscles[control_id]["att_point"], scn.objects[controllable_muscles[control_id]["att_name"]].worldOrientation))
        dir_W = (np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) + muscle_pos1_W) - (np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) + muscle_pos2_W)
        muscledist = np.linalg.norm(dir_W)
        musclevel = 0.0
        if controllable_muscles[control_id]["muscledist"]>-10000.0:
            musclevel = (muscledist - controllable_muscles[control_id]["muscledist"])/dt_bl
        controllable_muscles[control_id]["muscledist"] = muscledist
        return [muscledist, musclevel, 0.0]
    else:
        print("No spindle return implemented...")
        return None


def controlActivity(control_id, control_activity):
    global scn
    global object_rotation_sum
    if control_activity>1.0:
        control_activity = 1.0
        print("Activity should be between 0 and 1 ("+str(control_activity)+" in this case)")
    if control_activity<0.0:
        control_activity = 0.0
        print("Activity should be between 0 and 1 ("+str(control_activity)+" in this case)")
    # ------------------------------------- Positional servo ----------------------------------------------------------------------------------------------------
    if controllable_muscles[control_id]["type"]=="position_servo":
        # define rotation
        a_ = 1.0/(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0])
        b_ = -a_ * controllable_muscles[control_id]["limit_angles"][0]
        rot_n_old = controllable_muscles[control_id]["rot"]
        rot_n = (a_*( scn.objects[controllable_muscles[control_id]["ref_name"]].getRigidBodyJointAngles(0) )  + b_) 
        controllable_muscles[control_id]["rot"] = rot_n
        #~ print(controllable_muscles[control_id]["limit_angles"][0], controllable_muscles[control_id]["limit_angles"][1], rot_n)
        ROT = (control_activity-rot_n)*controllable_muscles[control_id]["P"]
        rotV_joint = np.fabs(controllable_muscles[control_id]["limit_angles"][1] - controllable_muscles[control_id]["limit_angles"][0]) * np.fabs(rot_n - rot_n_old)/dt_bl
        if controllable_muscles[control_id]["maxV"]!=None:
            if np.fabs(ROT)>controllable_muscles[control_id]["maxV"]:
                ROT = np.sign(ROT) * controllable_muscles[control_id]["maxV"]
        #~ printCsX(str(control_id)+('   ----------------------------- '), 'cyan')
        #~ printCsX(np.fabs(rot_n - rot_n_old)/dt_bl, 'yellow')
        #~ if np.fabs(rot_n - rot_n_old)/dt_bl>maxV and rot_n_old>-10000.0:
        #~ if rotV_joint>maxV:
            #~ printCsX(ROT, 'red')
            #~ ROT = np.sign(ROT) * maxV
            #~ printCsX(ROT, 'green')
            #~ ROT = 0.00001
            #~ scn.objects[controllable_muscles[control_id]["ref_name"]].worldAngularVelocity = [0.0,0.0,0.0]
            #~ scn.objects[controllable_muscles[control_id]["ref_name"]].worldLinearVelocity = [0.0,0.0,0.0]
            #~ scn.objects[controllable_muscles[control_id]["ref_name"]].worldAngularVelocity.magnitude = .000001
            #~ scn.objects[controllable_muscles[control_id]["ref_name"]].worldAngularVelocity.magnitude = 0.000001
        #~ else:
            #~ printCsX("Noes!!!", 'purple')
        ROT_W = np.array(newref(ROT*controllable_muscles[control_id]["orientation"], scn.objects[controllable_muscles[control_id]["ref_name"]].worldOrientation))
        #~ printCsX(scn.objects[controllable_muscles[control_id]["ref_name"]].worldAngularVelocity.magnitude, 'blue')
        #~ printCsX(ROT, 'cyan')
        #~ printCsX( dt_bl, 'red')
        #~ printCsX( controllable_muscles[control_id]["limit_angles"][1], 'red')
        #~ printCsX( rotV_joint, 'yellow')
        #~ printCsX(scn.objects[controllable_muscles[control_id]["ref_name"]].worldAngularVelocity, 'green')
        # apply rotation
        try:     object_rotation_sum[controllable_muscles[control_id]["ref_name"]]
        except:  object_rotation_sum[controllable_muscles[control_id]["ref_name"]] = np.array([0.0,0.0,0.0])
        try:     object_rotation_sum[controllable_muscles[control_id]["att_name"]]
        except:  object_rotation_sum[controllable_muscles[control_id]["att_name"]] = np.array([0.0,0.0,0.0])
        scn.objects[controllable_muscles[control_id]["ref_name"]].setAngularVelocity( object_rotation_sum[controllable_muscles[control_id]["ref_name"]]+ROT_W, False)
        scn.objects[controllable_muscles[control_id]["att_name"]].setAngularVelocity( object_rotation_sum[controllable_muscles[control_id]["att_name"]]-ROT_W, False)
        object_rotation_sum[controllable_muscles[control_id]["ref_name"]] +=  ROT_W
        object_rotation_sum[controllable_muscles[control_id]["att_name"]] += -ROT_W
        
        # TEST !!!
        #~ torque_ = 
        #~ scn.objects[controllable_muscles[control_id]["ref_name"]].setTorque()
        
    # ------------------------------------- Velocity servo ----------------------------------------------------------------------------------------------------
    if controllable_muscles[control_id]["type"]=="velocity_servo":
        # define rotation
        ROT = (control_activity-0.5)*controllable_muscles[control_id]["maxV"]*controllable_muscles[control_id]["orientation"]
        ROT_W = np.array(newref(ROT, scn.objects[controllable_muscles[control_id]["ref_name"]].worldOrientation))
        # apply rotation
        try:     object_rotation_sum[controllable_muscles[control_id]["ref_name"]]
        except:  object_rotation_sum[controllable_muscles[control_id]["ref_name"]] = np.array([0.0,0.0,0.0])
        try:     object_rotation_sum[controllable_muscles[control_id]["att_name"]]
        except:  object_rotation_sum[controllable_muscles[control_id]["att_name"]] = np.array([0.0,0.0,0.0])
        scn.objects[controllable_muscles[control_id]["ref_name"]].setAngularVelocity( object_rotation_sum[controllable_muscles[control_id]["ref_name"]]+ROT_W, False)
        scn.objects[controllable_muscles[control_id]["att_name"]].setAngularVelocity( object_rotation_sum[controllable_muscles[control_id]["att_name"]]-ROT_W, False)
        object_rotation_sum[controllable_muscles[control_id]["ref_name"]] +=  ROT_W
        object_rotation_sum[controllable_muscles[control_id]["att_name"]] += -ROT_W
    # ------------------------------------- Torque-Based ----------------------------------------------------------------------------------------------------
    if controllable_muscles[control_id]["type"]=="torque-based":
        # printCsX(control_activity, 'green')
        torqueN = controllable_muscles[control_id]["orientation"] * controllable_muscles[control_id]["maxF"] * control_activity
        torqueN_W = np.array(newref(torqueN, scn.objects[controllable_muscles[control_id]["ref_name"]].worldOrientation))
        scn.objects[controllable_muscles[control_id]["ref_name"]].applyTorque( torqueN_W, False )
        scn.objects[controllable_muscles[control_id]["att_name"]].applyTorque(-torqueN_W, False )
        bge.render.drawLine(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition, np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) + ( torqueN_W/20.0), [0.5,0.0,1.0], 1.0)
        bge.render.drawLine(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition, np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) + (-torqueN_W/20.0), [0.5,0.5,1.0], 1.0)
    # ------------------------------------- Single-Linear ----------------------------------------------------------------------------------------------------
    if controllable_muscles[control_id]["type"]=="single-linear":
        muscle_pos1_W = np.array(newref(controllable_muscles[control_id]["ref_point"], scn.objects[controllable_muscles[control_id]["ref_name"]].worldOrientation))
        muscle_pos2_W = np.array(newref(controllable_muscles[control_id]["att_point"], scn.objects[controllable_muscles[control_id]["att_name"]].worldOrientation))
        #~ impulse_W = normal_part_of((np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) ) - (np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) )) * controllable_muscles[control_id]["maxF"] * control_activity / float(bpy.context.scene.game_settings.fps)
        #~ impulse_W = np.array([0.0,0.0,1.0])
        impulse_W = normal_part_of((np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) + muscle_pos1_W) - (np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) + muscle_pos2_W)) * controllable_muscles[control_id]["maxF"] * control_activity / float(bpy.context.scene.game_settings.fps)
        #~ impulse_W = normal_part_of((np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) + np.array([0.0,0.0,1.0])) - (np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) + np.array([0.0,0.0,1.0]) )) * controllable_muscles[control_id]["maxF"] * control_activity / float(bpy.context.scene.game_settings.fps)
        #~ impulse_W = normal_part_of((muscle_pos1_W) - (muscle_pos2_W)) * controllable_muscles[control_id]["maxF"] * control_activity / float(bpy.context.scene.game_settings.fps)
        #~ scn.objects[controllable_muscles[control_id]["ref_name"]].applyForce( -impulse_W*60.0, False )
        #~ scn.objects[controllable_muscles[control_id]["att_name"]].applyForce(  impulse_W*60.0, False )
        #~ scn.objects[controllable_muscles[control_id]["ref_name"]].applyImpulse( controllable_muscles[control_id]["ref_point"], -impulse_W )
        #~ scn.objects[controllable_muscles[control_id]["att_name"]].applyImpulse( controllable_muscles[control_id]["att_point"],  impulse_W )
        #~ scn.objects[controllable_muscles[control_id]["ref_name"]].applyImpulse( [0.0,0.0,0.0], -impulse_W )
        #~ scn.objects[controllable_muscles[control_id]["att_name"]].applyImpulse( [0.0,0.0,0.0],  impulse_W )
        scn.objects[controllable_muscles[control_id]["ref_name"]].applyImpulse( muscle_pos1_W, -impulse_W )
        scn.objects[controllable_muscles[control_id]["att_name"]].applyImpulse( muscle_pos2_W,  impulse_W )
        bge.render.drawLine( np.array(scn.objects[controllable_muscles[control_id]["ref_name"]].worldPosition) + muscle_pos1_W , np.array(scn.objects[controllable_muscles[control_id]["att_name"]].worldPosition) + muscle_pos2_W, [1.0,0.0,0.0], 1.0 + 3.0*control_activity)        
    # ------------------------------------- OTHER ------------------------------------------------------------------------------------------------------------
    # torque-based muscle 
    #~ torqueN = controllable_muscles[control_id][4] * controllable_muscles[control_id][3] * control_activity * MuscleLengthGain_NONE(rot_n) * MuscleVelocityGain_NONE(vel_n)
    #~ torque_ = [
    #~ fsf*Ix*rot[0] - bV*(vel[0]) * Ix + torqueN_W[0], 
    #~ fsf*Iy*rot[1] - bV*(vel[1]) * Iy + torqueN_W[1], 
    #~ fsf*Iz*rot[2] - bV*(vel[2]) * Iz + torqueN_W[2]
    #~ ]








# =================================================================================================================================================
#                                       Useful functions


def progressbar(progress, pbnum=40):
    loadbar = ""
    for i in range(pbnum):
        if i<=int(progress*float(pbnum)):
            loadbar = loadbar + "="
        else:
            loadbar = loadbar + " "
    sys.stdout.write("\r"+"["+loadbar+"] "+str(int(progress*100.0))+"%")




def bezierPowCsX(t):
	if t<0.0: return 0.0;
	if t>1.0: return 1.0;
	return 3.0*(t)*(t) - 2.0*(t)*(t)*(t);
	

def bezierPowNpCsX(t):
    ret_ = np.zeros(t.shape)
    ret_[np.nonzero(t>=1.0)] = 1.0
    in_ids = np.nonzero((t>=0.0) * (t<=1.0))
    t_in = t[in_ids]
    ret_[in_ids] = 3.0*(t_in)*(t_in) - 2.0*(t_in)*(t_in)*(t_in)
    return ret_


def newref(vec, mat):
    return [
    vec[0]*mat[0][0] + vec[1]*mat[0][1] + vec[2]*mat[0][2],
    vec[0]*mat[1][0] + vec[1]*mat[1][1] + vec[2]*mat[1][2],
    vec[0]*mat[2][0] + vec[1]*mat[2][1] + vec[2]*mat[2][2]
    ]




def deg_to_rad(deg):
    return math.pi*deg/180.0

def rad_to_deg(rad):
    return rad*180.0/math.pi



def sign_triangle(v1, v2, v3):
    return (v1[0]-v3[0])*(v2[1]-v3[1]) - (v2[0]-v3[0])*(v1[1]-v3[1])



def point_in_triangle(pt, v1, v2, v3):
    st1 = sign_triangle(pt, v1, v2)
    st2 = sign_triangle(pt, v2, v3)
    st3 = sign_triangle(pt, v3, v1)
    if   st1<=0.0 and st2<=0.0 and st3<=0.0:
        return True;
    elif st1>=0.0 and st2>=0.0 and st3>=0.0:
        return True;
    else:
        return False;



def inverse2x2(Tmat):
    det_ = Tmat[0,0]*Tmat[1,1] - Tmat[0,1]*Tmat[1,0]
    if det_==0.0:return Tmat*0.0;
    return np.array([[Tmat[1,1], -Tmat[0,1]], [-Tmat[1,0], Tmat[0,0]]])/det_
    

def matrix_vec_mult_2x2(mat, vec):
    return np.array([ mat[0,0]*vec[0]+mat[0,1]*vec[1] ,  mat[1,0]*vec[0]+mat[1,1]*vec[1] ])


def simple_gauss3D(vec1, vec2, sigma):
    dist = norm(sub(vec1, vec2))
    return math.exp(-dist*dist/sigma)
    


def sigmoid_CsX(x, Tk, dec):
    return 1.0/(1.0+np.exp(-(x-dec)/Tk))





def rotate_around_vector(vec, vec_ref, angle):
    vec_final = np.copy(vec)
    vec_ref_N = normal_part_of(vec_ref)
    circleY = np.cross(vec_ref_N, vec)
    rotation_norm = np.linalg.norm(circleY)
    circleY = normal_part_of(circleY)
    #~ print(circleX)
    if np.linalg.norm(circleY)==0.0: 
        return vec_final
    circleX = normal_part_of(np.cross(circleY, vec_ref))
    vec_around = vec_ref_N * np.dot(vec_ref_N, vec)
    #~ print(vec_around)
    #~ print(circleX)
    #~ print(circleY)
    #~ print(np.dot(circleX, circleY))
    #~ print(np.linalg.norm(vec_around))
    vec_final  = vec_around + rotation_norm*circleX*(np.cos(angle)) + rotation_norm*circleY*(np.sin(angle))
    return vec_final



def normal_part_of(arr):
    arr_norm = np.linalg.norm(arr)
    if arr_norm>0.0: 
        return arr/arr_norm
    else:
        return arr

def set_joint_(child_object_name, DoF_vector, limit_angles, maxforce):
    bpy.context.scene.objects.active = bpy.data.objects[child_object_name]
    #~ bpy.ops.object.constraint_add(type='RIGID_BODY_JOINT')
    bpy.context.object.constraints["Rigid Body Joint"].show_pivot = True
    rotv_xy = [ DoF_vector[0], DoF_vector[1] ]
    rotv_xy = normal_part_of(rotv_xy)
    rotv_  = DoF_vector
    rotv_ = normal_part_of(rotv_)
    theta_ =  math.acos(rotv_xy[0])
    if rotv_xy[1]<0.0:theta_ *= -1.0
    phi_   = -math.asin(rotv_[2])
    #~ print(rotv_, rotv_xy, theta_, phi_)
    bpy.context.object.constraints["Rigid Body Joint"].axis_z = theta_
    bpy.context.object.constraints["Rigid Body Joint"].axis_y = phi_
    bpy.context.object.constraints["Rigid Body Joint"].axis_x = 0.0
    #joint_registry[bone_[0]] = [pivot_x_, pivot_y_, pivot_z_]
    #joint_bvectors_registry[bone_[0]] = []
    #joint_bvectors_registry[bone_[0]].append( [0.0, 0.0, 0.15] )
    bpy.context.object.constraints["Rigid Body Joint"].use_linked_collision = True
    bpy.context.object.constraints["Rigid Body Joint"].pivot_type = 'HINGE'
    #bpy.context.object.constraints["Rigid Body Joint"].pivot_type = 'GENERIC_6_DOF'
    bpy.context.object.constraints["Rigid Body Joint"].use_angular_limit_x = True
    bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x = deg_to_rad(limit_angles[0])
    bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x = deg_to_rad(limit_angles[1])     










# =================================================================================================================================================
#                                       Calling controller.py


dt_bl = 1.0/float(bpy.context.scene.game_settings.fps)
#~ dt_bl = 1.0/float(bpy.context.scene.game_settings.fps + bpy.context.scene.game_settings.physics_step_sub)
t_bl  = 0.0
i_bl  = 0
reinitialized = 0
gamestate = None
def evolve_():
    global gamestate
    global reinitialized, t_bl, dt_bl, i_bl, scn
    global object_rotation_sum
    #~ if reinitialized==0:
        #~ print("Restarting")
        #~ reinitialized += 1
        #~ scn.restart()
    #~ elif reinitialized==1:
        #~ scn = g.getCurrentScene()
    for orsk in object_rotation_sum.keys():
        object_rotation_sum[orsk] = np.array([0.0,0.0,0.0])
    if gamestate==None:
        gamestate = evolve()
    elif gamestate=="restart":
        bge.logic.restartGame()
    elif gamestate=="end":
        bge.logic.endGame()
    else: 
        gamestate=None
    t_bl += dt_bl
    i_bl += 1




class MyImporter(object):
    def find_module(self, module_name, package_path):
        # Return a loader
        #~ global modulenames_CsX
        global bpy
        global sys
        print("finding module", module_name)
        #~ if module_name not in sys.modules.keys():
            #~ modulenames_CsX.append(module_name)
        if module_name not in sys.modules.keys():
            bpy.modules_CsX[module_name] = None
        return None
    def load_module(self, module_name):
        # Return a module
        print("loading module")
        return self





if firstTimeCsX:
    bpy.modules_CsX = {}
else:
    for mnC in bpy.modules_CsX.keys():
        sys.modules[mnC] = bpy.modules_CsX[mnC]    


sys.meta_path.insert(0,MyImporter())

# Calling controller
exec(compile(open(controller_name).read(), controller_name, 'exec'))

sys.meta_path.pop(0)


for mnC in bpy.modules_CsX.keys():
    #~ print(mnC)
    #~ if sys.modules[mnC] not in modules_CsX:
    #~ modules_CsX.append( sys.modules[mnC] )
    if mnC in sys.modules.keys():
        bpy.modules_CsX[mnC] = sys.modules[mnC] 














'''

filename = "../../SCRIPTS/script_importpaths.py"
exec(compile(open(filename).read(), filename, 'exec'))

filename = "../../SCRIPTS/script_fcts.py"
exec(compile(open(filename).read(), filename, 'exec'))


import numpy as np
#~ from pylab import *
import scipy

# needs python3-imaging
import scipy.misc

from mathutils import Euler
import bpy
import bge
import math
import VideoTexture


from rg import accelerometer_registry

from rg import olfactory_object

from rg import camera_list

from rg import joint_registry

from rg import joint_bvectors_registry

from rg import sensory_map
from rg import sensory_map_pos
from rg import sensory_map_objname
from rg import sensory_map_isused
from rg import sensory_ids_used
from rg import objstr_ids

from rg import vertex_registry
from rg import uv_registry
from rg import face_local_registry
from rg import familial_registry



controllername = "controller.py"
exec(compile(open(filename).read(), controllername, 'exec'))
#~ init_WPs()




position_save_tab = {}
#~ position_save_tab["forearm.L.008"] = []



Pos_tmp = None
Work_tmp = 0.0

game_ended = False

first_unique_run = True

img1_global = None
muscle_pos_global = None
muscle_inputs = None
spindle_outputs = None
vestibular_outputs = None
tmp_velocities = None
olfactory_outputs = None

total_spikes = {}
total_spikes["senders"] = []
total_spikes["times"  ] = []

scn = None

muscle_indices = {}
joint_forces = {}
muscle_tension = {}
dynamic_functions = {}

i_evo = 0




filename = "../../SCRIPTS/script_collision.py"
exec(compile(open(filename).read(), filename, 'exec'))



filename = "../../SCRIPTS/script_muscle.py"
exec(compile(open(filename).read(), filename, 'exec'))


# ============================================================================================================================================================================

rot1s = []
rot2s = []
f1s = []
vel1s = []
f2s = []
vel2s = []

def evolve_():
    global rot1s, rot2s, f1s, vel1s, f2s, vel2s 
    
    global i_evo, i_evo_sock, s, joint_forces, muscle_tension, dynamic_functions, scn, sensory_map, src1, accelerometer_registry, MFM1, img1_global, muscle_pos_global, muscle_inputs, spindle_outputs, vestibular_outputs, sensory_ids_used
    global tmp_velocities, olfactory_outputs, camera_list, muscle_indices, Pos_tmp, Work_tmp, scene_ini, first_unique_run
    global game_ended
    global total_spikes
    if game_ended==False:
        import bge
        import GameLogic as g
        scn = g.getCurrentScene()
        cont = g.getCurrentController()
        own = cont.owner
        
        # OLFACTORY INPUT
        if olfactory_object:
            num_smells = 3
            olfactory_outputs_tmp = np.zeros(num_smells, dtype=np.float32)
            for objkeys in bpy.data.objects.keys():
                for ismell in range(num_smells):
                    smellampl = None
                    try:
                        smellampl = bpy.data.objects[objkeys]['smell'+str(ismell+1)]
                    except:
                        smellampl=None
                    if smellampl!=None:
                        pos1 = scn.objects[objkeys].worldPosition
                        pos2 = scn.objects["obj_"+olfactory_object].worldPosition
                        #olfactory_outputs_tmp[ismell] += smellampl * np.exp(-norm(pos1-pos2)/5.0)
                        olfactory_outputs_tmp[ismell] += smellampl * np.exp(-(norm(pos1-pos2)**2.0)/10.0)
            olfactory_outputs = np.copy(olfactory_outputs_tmp)
            MFM1.set_input_1D(olfactory_outputs, 0, 0.0, 1.0, 0.0, 1.0, -60.0, -43.0, -40.0, -23.0, 0)
        
        # VESTIBULAR INPUT
        if len(accelerometer_registry)>0:
            vel_ = []
            steady_acc_ = []
            for acc_name in accelerometer_registry:
                av_ = scn.objects["obj_"+acc_name].localAngularVelocity
                lv_ = scn.objects["obj_"+acc_name].localLinearVelocity
                vel_.append([av_[0], av_[1], av_[2], lv_[0], lv_[1], lv_[2]])
                steady_acc_tmp = newref([0.0,0.0,-9.80], scn.objects["obj_"+acc_name].worldOrientation.transposed())
                steady_acc_.append([steady_acc_tmp[0], steady_acc_tmp[1], steady_acc_tmp[2]])
            tmp_velocities_new = np.array(vel_, dtype=np.float32)
            if tmp_velocities==None:tmp_velocities = np.copy(tmp_velocities_new)
            vestibular_outputs_tmp = ((tmp_velocities_new - tmp_velocities) * float(bpy.context.scene.game_settings.fps) )
            vestibular_outputs_tmp[:,3:] -= np.array(steady_acc_, dtype=np.float32)
            vestibular_outputs_tmp[:,:] *= 0.1
            vestibular_outputs_tmp_tmp = np.zeros((vestibular_outputs_tmp.shape[0], vestibular_outputs_tmp.shape[1]*2), dtype=np.float32)
            vestibular_outputs_tmp_tmp[:,0::2] = -vestibular_outputs_tmp[:,:] * np.float32(vestibular_outputs_tmp[:,:]<0.0)
            vestibular_outputs_tmp_tmp[:,1::2] = +vestibular_outputs_tmp[:,:] * np.float32(vestibular_outputs_tmp[:,:]>0.0)
            vestibular_outputs = np.copy(vestibular_outputs_tmp_tmp).T
            tmp_velocities = np.copy(tmp_velocities_new)
            #print(vestibular_outputs)
            MFM1.set_input_2D(vestibular_outputs, 1, 0.4, 0.0, 1.0, 1.0,-60.0, -65.0, -40.0, -45.0, 0)
        
        # MOTOR OUTPUT & SPINDLE INPUT
        if first_unique_run==True:
            muscle_indices = {}; mi = 0
            positions_ = []
            for kk in joint_registry.keys():
                for ijs in range(len(joint_registry[kk])):
                    if(joint_registry[kk][ijs][1]>0.0):
                        for iside in [-1,1]:
                            vecc = scn.objects["obj_"+kk].worldPosition
                            positions_.append( [vecc[0], vecc[1] + 0.05*float(iside), vecc[2]] )
                            muscle_indices[kk+"_"+str(ijs+1)+"_"+str(iside)] = mi
                            mi += 1
            muscle_pos_global = np.array(positions_, dtype=np.float32)
            muscle_pos_global /= (np.max(muscle_pos_global[:,0]) - np.min(muscle_pos_global[:,0]))
            muscle_pos_global = muscle_pos_global - np.mean(muscle_pos_global, axis=0)
            muscle_pos_global *= 8.0
            muscle_inputs = np.zeros(muscle_pos_global.shape[0], dtype=np.float32)
        MFM1.set_input_1D( muscle_inputs, 5, 0.5, 0.5, 0.2, 1.0, -19.0, -21.0, 15.0, -1.0+20.0, 1 )
        i_MI = 0
        spindle_outputs_list = []
        for kk in joint_registry.keys():
            for ijs in range(len(joint_registry[kk])):
                if(joint_registry[kk][ijs][1]>0.0):
                    for iss,iside in enumerate([-1,1]):
                        joint_registry[kk][ijs][5+iss] += (1.0 / float(bpy.context.scene.game_settings.fps)) * (1.0/0.01) * (muscle_inputs[i_MI] - joint_registry[kk][ijs][5+iss])
                        #~ joint_registry[kk][ijs][5+iss] = muscle_inputs[i_MI]
                        #~ joint_registry[kk][ijs][5+iss] = 0.0
                        #~ if iss==0 and len(vel1s)>0:
                            #~ joint_registry[kk][ijs][5+iss] = (0.2)*vel1s[len(vel1s)-1]
                        #~ if iss==1 and len(vel2s)>0:
                            #~ joint_registry[kk][ijs][5+iss] = (0.2)*vel2s[len(vel2s)-1]
                        [vel_n, tension_, rot_n] = muscle_transfer_function(
                        name_1 = "obj_"+kk, name_2 = "obj_"+familial_registry[kk], 
                        muscle_direction = joint_registry[kk][ijs][0], 
                        muscle_force = float(iside) * joint_registry[kk][ijs][1], neural_activity = joint_registry[kk][ijs][5+iss], 
                        objects_ = scn.objects,
                        jr = joint_registry[kk][ijs],
                        oldrot = joint_registry[kk][ijs][7+iss],
                        dt = 1.0/float(bpy.context.scene.game_settings.fps))
                        joint_registry[kk][ijs][7+iss] = rot_n
                        Work_tmp += np.abs(tension_);
                        i_MI+=1
                        # computing spindles
                        spindle_outputs_list.append([ vel_n, tension_, rot_n ])
                        #~ if iss==0:
                            #~ rot1s.append((rot_n - 0.5) * 1.0)
                            #~ vel1s.append(vel_n * 0.05)
                            #~ f1s.append(  joint_registry[kk][ijs][5+iss])
                        #~ if iss==1:
                            #~ rot2s.append((rot_n - 0.5) * 1.0)
                            #~ vel2s.append(vel_n * 0.05)
                            #~ f2s.append(  joint_registry[kk][ijs][5+iss])
        spindle_outputs = np.array(spindle_outputs_list, dtype=np.float32)
        MFM1.set_input_2D( spindle_outputs, 4, 0.5, 0.5, 0.2, 1.0, -60.0, -21.0, -40.0, -1.0+20.0, 0 )
        
        # SENSORY INPUT
        #~ if first_unique_run==True:
        if i_evo==0:
            for kk in familial_registry.keys():
                scn.objects["obj_"+kk].collisionCallbacks.clear()
                dynamic_functions["obj_"+kk] = collision_function_builder("obj_"+kk)
                try:
                    scn.objects["obj_"+kk].collisionCallbacks.append(dynamic_functions["obj_"+kk])
                except:
                    print("Bone ", kk, " has no associated polygons...")
        sensory_map[np.nonzero(sensory_map>1.0)] = 1.0
        #scipy.misc.imsave("../../../HighD/S1/SENSORY_INPUT/test_"+('%05d' % (i_evo+1))+".png", sensory_map)
        
        #~ MFM1.set_input_3D(sensory_map, 2, 1.0, 1.0, 1.0, 1.0, -60.0, 23.0, -40.0, 43.0, 0 )
        ssmap_tmp = np.copy(sensory_map[:,:,:3])
        MFM1.set_input_3D(ssmap_tmp, 2, 1.0, 1.0, 1.0, 1.0, -60.0, 23.0, -40.0, 43.0, 0 )
        
        
        #~ # VISUAL INPUT
        #~ tvcam = scn.objects[camera_list[0]]
        #~ src1 = VideoTexture.ImageRender(scn,tvcam)
        #~ #src1.capsize = 512,512
        #~ src1.capsize = 128,128
        #~ imX,imY = src1.size
        #~ #img1 = np.array(src1.image).reshape((imY, imX, 4))
        #~ #img1_global = np.float32(img1[:,:,:])/256.0
        #~ img1_global = np.array(src1.image, dtype=np.float32)/256.0
        #~ #scipy.misc.imsave("../../HighD/V1/VISUAL_INPUT/test_"+('%05d' % (i_evo+1))+".png", img1[::-1,:,:])
        #~ if going_MFM:
            #~ MFM1.set_input_splat(img1_global, 3, 0.0, 0.6, 1.0, 1.0, -60.0, 1.0, -40.0, 21.0, 0, imX, imY, 4)
        
        
        # EVOLUTIONARY ALGORITHM
        numsecs = 10.0
        if i_evo>(int(float(bpy.context.scene.game_settings.fps)*numsecs)) and evo_algorithm:
            Pos_tmp_2 = [scn.objects["obj_head"].worldPosition[0], scn.objects["obj_head"].worldPosition[1], scn.objects["obj_head"].worldPosition[2]] 
            print(Pos_tmp, Pos_tmp_2)
            worksc = 1.0 - sigmoid_CsX(x = Work_tmp/numsecs, dec = 700.0, Tk = 200.0)  # decreasing with energy consumptiongame_ended
            print("Work:", Work_tmp/numsecs, "  Sig:", worksc)
            Rew_tmp = (-1.0)*(Pos_tmp_2[1] - Pos_tmp[1]) * np.exp(-np.abs(Pos_tmp_2[0] - Pos_tmp[0])/1.0) * worksc
            Work_tmp = 0.0
            next_generation( Rew_tmp )
            #~ MFM1.sys_clear_CsX();
            i_evo=0
            if current_generation>200:
                bge.logic.endGame()                
        if i_evo==0:
            MFM1.deleteAllCells(); 
            MFM1.deleteRectangles(); 
            MFM1.deleteMessages(); 
            genesis(MFM1)
        t_to_sim = 1.0/float(bpy.context.scene.game_settings.fps)
        if usePtNeu:
            nest.SetStatus(cells, "I_e", 0.0)
            for i_in_,in_ in enumerate(input_neurons_):
                if i_in_<sensory_map.shape[2]:
                    nest.SetStatus(list(in_+1), "I_e", 0.0 + (500.0*np.float64(sensory_map[vtxss[i_in_], vtyss[i_in_], i_in_])))        
            I_e_old = np.array(nest.GetStatus(cells, "I_e"))
            nest.SetStatus(cells, "I_e", I_e_old + extraInput)
            #~ nest.Simulate(10.0)
            nest.Simulate(float(int(t_to_sim*1000.0)))
            spikes = nest.GetStatus(rec, 'events')[0]
            nest.SetStatus(rec, [{"n_events": 0}])
            #~ print(spikes)
        if usePtNeu:
            MFM1.sendSpikesPtNeu(np.int64(spikes["senders"]-1))
            total_spikes["senders"] += list(spikes["senders"])
            total_spikes["times"  ] += list(spikes["times"  ])
            del spikes
        # simulate MFM
        #~ print(int(t_to_sim*1000.0))
        #~ for iMFMevo in range(15):
        for iMFMevo in range(int(t_to_sim*1000.0)):
            MFM1.evolve(0.001)
        #~ MFM1.draw()
        if makevideo:
            zoom_ = 300.0;
            #~ Rspeed = 0.002
            #~ theta = Rspeed*float(i_evo)
            #~ phi   = np.pi*(0.35)
            #~ MFM1.draw(export_image_name = ("../../../HighD/VIDEO/img"+str('%05d' % i_evo)+".png").encode("utf-8"), camera_coord = [zoom_*np.cos(theta)*np.sin(phi), zoom_*np.sin(theta)*np.sin(phi), zoom_*np.cos(phi)+100.0, 0.0,0.0,0.0+100.0, 0.0,0.0,1.0], draw_only_PtNeus = 1)
            #~ theta = -np.pi * 0.4
            #~ phi   = np.pi*(0.25)
            #~ camera_coord_initial = np.array([0.0,-300.0,300.0, 0.0,-1.3,93.5, 0.0, 0.0, 1.0])
            #~ camera_coord_final   = np.array([147.8293,-23.7114,226.512, 0.0,-1.3,73.5, 0.0,0.0,1.0])
            camera_coord_initial = np.array([105.0,-180.0,300.0, 0.0,-1.3,83.5, 0.0,0.0,1.0])
            camera_coord_final   = np.array([105.0,-180.0,300.0, 0.0,-1.3,83.5, 0.0,0.0,1.0])
            MFM1.draw(export_image_name = ("../../../HighD/VIDEO/img"+str('%05d' % i_evo)).encode("utf-8"), export_array_name = ("../../../HighD/VIDEO/input"+str('%05d' % i_evo)).encode("utf-8"), 
            camera_coord = list( camera_coord_initial + bezierPowCsX(float(i_evo)/240.0) * (camera_coord_final - camera_coord_initial) ), 
            draw_only_PtNeus = 1)
        else:
            MFM1.draw()
            #~ MFM1.draw(draw_only_PtNeus = 1)
        use_controller = MFM1.events(viewtype = 1)
        #~ sensory_map[sensory_ids_used[0], sensory_ids_used[1], sensory_ids_used[2]] = 0.0
        sensory_map[:,:,:] *= 0.0
        for kk in position_save_tab.keys():
            position_save_tab[kk].append( [scn.objects["obj_"+kk].worldPosition[0], scn.objects["obj_"+kk].worldPosition[1], scn.objects["obj_"+kk].worldPosition[2]] )
        i_evo += 1
        #~ print(i_evo)
        if first_unique_run==True:first_unique_run=False;
        if looplimit!=None:
            if i_evo>=looplimit:
                
                #~ import pylab as plb
                #~ plb.figure();
                #~ plb.plot(range(len(rot1s)), rot1s, '-', color=[0.0,0.0,1.0]); 
                #~ plb.plot(range(len(vel1s)), vel1s, '-', color=[0.0,1.0,0.0]); 
                #~ plb.plot(range(len(f1s))  , f1s  , '-', color=[1.0,0.0,0.0]); 
                #~ plb.legend(['angle', 'velocity', 'force'])
                #~ plb.figure();
                #~ plb.plot(range(len(rot2s)), rot2s, '-', color=[0.0,0.0,1.0]); 
                #~ plb.plot(range(len(vel2s)), vel2s, '-', color=[0.0,1.0,0.0]); 
                #~ plb.plot(range(len(f2s))  , f2s  , '-', color=[1.0,0.0,0.0]); 
                #~ plb.legend(['angle', 'velocity', 'force'])
                #~ plb.show()
                
                if usePtNeu and game_ended==False and len(total_spikes["times"])>0:
                    #~ eventsM = nest.GetStatus(multiM)[iddd-1]['events']
                    spikes_to_Fr_fMRI("../../../HighD/MouseVfMRI", timestep_limits = [0, looplimit+10], dt = 1000.0/float(bpy.context.scene.game_settings.fps), boundaries_min=[-6600.0, -4000.0, -5700.0], boundaries_max=[6600.0, 4000.0, 5700.0], resolution = [132,80,114], positions = np.float32(np.vstack((circuit["x"], circuit["y"], circuit["z"])).T), spikes = total_spikes, numtype = np.float)
                    spikes_to_rasterplot(filename = "raster", spikes = total_spikes, regions = np.array(circuit["Larea"])[:len(cells)], colors = (np.float32(np.vstack((circuit["colorx"], circuit["colory"], circuit["colorz"])).T) / 255.0)[:,:len(cells)] )
                    #~ spikes_to_rasterplot(filename = "raster", spikes = total_spikes, regions = whiskmapnumber, colors = (np.float32(np.vstack((circuit["colorx"], circuit["colory"], circuit["colorz"])).T) / 255.0)[:,:len(cells)] )
                    
                game_ended = True
                print("Game Ended !!!")
                bge.logic.endGame()













'''


