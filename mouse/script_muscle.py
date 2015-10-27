
# ----------------------------------------------------- Stimulus Gain ------------------------------------------------------------------------

def MuscleStimulusGain_NONE(inp_):
    return inp_


def MuscleStimulusGain_sigmoid(inp_, A, B, C):
    return ((m_fltB/(1.0+np.exp(C*(A-inp_)))) + D);



# ----------------------------------------------------- Velocity Gain -----------------------------------------------------------------------


def MuscleVelocityGain_NONE(inp_):
    return 1.0



def MuscleVelocityGain_Hill1(vel, a = 0.25, b = 0.25):
    if vel<0.0: 
        return 1.0;
    else: 
        return b*((1.0 + a)/(vel + b)) - a


# ------------------------------------------------------ Length Gain ------------------------------------------------------------------------

def MuscleLengthGain_NONE(phi):
    return 1.0

def MuscleLengthGain_square(phi, phi_rest, phi_width):
    Lnew = L - Lrest;
    scFac = (-((Lnew**2.0)/Lwidth) + 1.0);
    if scFac<0.0: scFac = 0.0;
    return scFac


# -----------------------------------------------------------------------------------------------------------------------------------------------

def muscle_transfer_function(name_1, name_2, muscle_direction, muscle_force, neural_activity, objects_, jr, oldrot, dt):
    phi1 = objects_[name_1].worldOrientation.transposed()
    phi2 = objects_[name_2].worldOrientation
    rot = (phi1 * phi2).to_euler()
    #~ phiV1 = objects_[name_1].worldAngularVelocity
    #~ phiV2 = objects_[name_2].worldAngularVelocity
    #rot = (objects_["obj_"+kk].worldOrientation.transposed() * objects_["obj_"+familial_registry[kk]].worldOrientation).to_euler()
    #~ vel = dot((phiV1 - phiV2), jr[0])
    #vel *= jr[1]/15.0
    a_ = 1.0/(jr[2][1] - jr[2][0])
    b_ = -a_ * jr[2][0]
    rot_n = (a_*rad_to_deg( dot(rot, jr[0]) ) * (-1.0) + b_) 
    if muscle_force>0.0:
        rot_n = 1.0 - rot_n
        #~ vel_n *= -1.0    
    vel_n = (rot_n - oldrot)/dt
    if oldrot<-10000.0:vel_n = 0.0
    vel = mult( jr[0], vel_n )
    #~ print(objects_[name_1].localInertia, objects_[name_2].localInertia)
    I1x = objects_[name_1].localInertia[0]; I1y = objects_[name_1].localInertia[1]; I1z = objects_[name_1].localInertia[2]
    I2x = objects_[name_2].localInertia[0]; I2y = objects_[name_2].localInertia[1]; I2z = objects_[name_2].localInertia[2]
    Ix = I1x; Iy = I1y; Iz = I1z
    if I2x<I1x:Ix=I2x
    if I2y<I1y:Iy=I2y
    if I2z<I1z:Iz=I2z
    #~ fsf = 30000.0
    #~ bV  = 800.0
    fsf = jr[3]
    bV  = jr[4]
    if muscle_force>0.0:
        vel = mult(vel, -1.0)
    torqueN = mult(muscle_direction, muscle_force * MuscleStimulusGain_NONE( neural_activity ) * MuscleLengthGain_NONE(rot_n) * MuscleVelocityGain_NONE(vel_n) )
    #~ torqueN = mult(muscle_direction, muscle_force * MuscleStimulusGain_NONE( neural_activity ) * MuscleLengthGain_NONE(rot_n) * MuscleVelocityGain_Hill1(-vel_n/20.0, 0.5, 0.5) )
    torqueN_W = newref(torqueN, objects_[name_1].worldOrientation)
    torque_ = [
    fsf*Ix*rot[0] - bV*(vel[0]) * Ix + torqueN_W[0], 
    fsf*Iy*rot[1] - bV*(vel[1]) * Iy + torqueN_W[1], 
    fsf*Iz*rot[2] - bV*(vel[2]) * Iz + torqueN_W[2]
    ]
    torque_m = [-torque_[0],-torque_[1],-torque_[2]]
    objects_[name_1].applyTorque( torque_ , False )
    objects_[name_2].applyTorque( torque_m, False )
    bge.render.drawLine(objects_[name_1].worldPosition, add(objects_[name_1].worldPosition, div(torque_ , 40.0)), [0.5,0.0,1.0])
    bge.render.drawLine(objects_[name_2].worldPosition, add(objects_[name_2].worldPosition, div(torque_m, 40.0)), [0.5,0.5,1.0])
    #
    return [vel_n, norm(torque_), rot_n]











