# =================================================================================================================================================
#                                       Import modules



# =================================================================================================================================================
#                                       Creating muscles

#~ muscle_ids = {}
#~ [ muscle_ids["forearm.L_FLEX"], muscle_ids["forearm.L_EXT"] ] = setTorqueMusclePair(reference_object_name = "obj_forearm.L",  attached_object_name = "obj_upper_arm.L",  maxF = 5.0)

#~ servo_ids = {}
#~ servo_ids["forearm.L"] = setVelocityServo(reference_object_name = "obj_forearm.L",  attached_object_name = "obj_upper_arm.L",  maxV = 10.0)

#~ servo_ids = {}
#~ servo_ids["forearm.L"]    = setPositionServo(reference_object_name = "obj_forearm.L",    attached_object_name = "obj_upper_arm.L", P = 50.0)
#~ servo_ids["forearm.R"]    = setPositionServo(reference_object_name = "obj_forearm.R",    attached_object_name = "obj_upper_arm.R", P = 50.0)
#~ servo_ids["upper_arm.L"]  = setPositionServo(reference_object_name = "obj_upper_arm.L",  attached_object_name = "obj_shoulder.L", P = 50.0)
#~ servo_ids["upper_arm.R"]  = setPositionServo(reference_object_name = "obj_upper_arm.R",  attached_object_name = "obj_shoulder.R", P = 50.0)
#~ servo_ids["shin_lower.L"] = setPositionServo(reference_object_name = "obj_shin_lower.L", attached_object_name = "obj_shin.L", P = 50.0)
#~ servo_ids["shin_lower.R"] = setPositionServo(reference_object_name = "obj_shin_lower.R", attached_object_name = "obj_shin.R", P = 50.0)
#~ servo_ids["shin.L"]       = setPositionServo(reference_object_name = "obj_shin.L",       attached_object_name = "obj_thigh.L", P = 50.0)
#~ servo_ids["shin.R"]       = setPositionServo(reference_object_name = "obj_shin.R",       attached_object_name = "obj_thigh.R", P = 50.0)


# =================================================================================================================================================
#                                       Network creation





# =================================================================================================================================================
#                                       Evolve function
def evolve():
    #~ print("Step:", i_bl, "  Time:", t_bl)
    # ------------------------------------- Visual ------------------------------------------------------------------------------------------------
    #visual_array     = getVisual(camera_name = "Meye", max_dimensions = [256,256])
    #scipy.misc.imsave("test_"+('%05d' % (i_bl+1))+".png", visual_array)
    # ------------------------------------- Olfactory ---------------------------------------------------------------------------------------------
    olfactory_array  = getOlfactory(olfactory_object_name = "obj_nose", receptor_names = ["smell1", "plastic1"])
    # ------------------------------------- Taste -------------------------------------------------------------------------------------------------
    taste_array      = getTaste(    taste_object_name =     "obj_mouth", receptor_names = ["smell1", "plastic1"], distance_to_object = 1.0)
    # ------------------------------------- Vestibular --------------------------------------------------------------------------------------------
    vestibular_array = getVestibular(vestibular_object_name = "obj_head")
    # ------------------------------------- Sensory -----------------------------------------------------------------------------------------------
    # ------------------------------------- Proprioception ----------------------------------------------------------------------------------------
    #~ spindle_FLEX = getMuscleSpindle(control_id = muscle_ids["forearm.L_FLEX"])
    #~ spindle_EXT  = getMuscleSpindle(control_id = muscle_ids["forearm.L_EXT"])
    # ------------------------------------- Neural Simulation -------------------------------------------------------------------------------------
    # nest.Simulate()
    # ------------------------------------- Muscle Activation -------------------------------------------------------------------------------------
    act_tmp      = 0.5 + 0.5*np.sin(2.0*t_bl)
    anti_act_tmp = 1.0 - act_tmp
    # Torque-based muscles
    #~ controlActivity(control_id = muscle_ids["forearm.L_FLEX"], control_activity = act_tmp)
    #~ controlActivity(control_id = muscle_ids["forearm.L_EXT"] , control_activity = anti_act_tmp)
    # Servos
    #~ controlActivity(control_id = servo_ids["forearm.L"], control_activity = act_tmp)
    #~ controlActivity(control_id = servo_ids["forearm.R"], control_activity = act_tmp)
    #~ controlActivity(control_id = servo_ids["obj_shin_lower.L"], control_activity = act_tmp)
    #~ controlActivity(control_id = servo_ids["obj_shin_lower.R"], control_activity = act_tmp)
    #~ controlActivity(control_id = servo_ids["forearm.L"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["forearm.R"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["upper_arm.L"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["upper_arm.R"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["shin_lower.L"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["shin_lower.R"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["shin.L"], control_activity = 0.5)
    #~ controlActivity(control_id = servo_ids["shin.R"], control_activity = 0.5)


















#bge.logic.endGame()
#bge.logic.restartGame()
#scn.reset()


