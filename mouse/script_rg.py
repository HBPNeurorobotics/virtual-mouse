# Required simulation-specific variables
# armaturename, meshname


import bpy
from mathutils import Euler
import math

# allows blender to use external python modules
filename = "../../SCRIPTS/script_importpaths.py"
exec(compile(open(filename).read(), filename, 'exec'))


# use numpy 
import numpy as np

# some basic functions
filename = "../../SCRIPTS/script_fcts.py"
exec(compile(open(filename).read(), filename, 'exec'))



# olfactory sensor
olfactory_object = None

# acceleration sensors
accelerometer_registry = []

# list of eye cameras
camera_list = []

# joints
joint_bvectors_registry = {}
joint_registry = {}





# for each bone tells his parent (useful after bones have been disconnected)
familial_registry = {}
# meshes for rigid bodies
vertex_registry = {}
uv_registry     = {}
face_id_registry_name = []
face_id_registry_num  = []
face_id_registry_name_2nd  = []
face_id_registry_num_2nd  = []
face_local_registry    = {}
face_material_registry = {}
def init():
    import gc
    gc.disable()
    global familial_registry, vertex_registry, uv_registry, face_id_registry_name, face_id_registry_num, face_id_registry_name_2nd, face_id_registry_num_2nd, face_local_registry, face_material_registry
    familial_registry = {}
    face_id_registry_name = []
    face_id_registry_num  = []
    face_id_registry_name_2nd = []
    face_id_registry_num_2nd  = []
    for bone_ in bpy.data.objects[armaturename].pose.bones.items():
        #~ print(bone_[0])
        vertex_registry[bone_[0]]     = []
        uv_registry[    bone_[0]]     = []
        face_local_registry[   bone_[0]] = []
        face_material_registry[bone_[0]] = []
        if bone_[1].parent!=None:
            parentname = str(bone_[1].parent)
            #~ print (parentname[parentname.find("(\"")+2:parentname.find("\")")])
            familial_registry[bone_[0]] = parentname[parentname.find("(\"")+2:parentname.find("\")")]
            joint_registry[   bone_[0]] = []
        else:
            #~ print ("!!!")
            familial_registry[bone_[0]] = "None"
    bpy.ops.object.mode_set(mode='EDIT')
    #bpy.ops.object.select_all(action='TOGGLE')
    bpy.ops.armature.parent_clear(type='CLEAR')
    bpy.ops.armature.parent_clear(type='DISCONNECT')
    for poly in bpy.data.objects[meshname].data.polygons:
        max_group_id = -1
        max_weight   = -0.0
        for vert in poly.vertices:
            for group in bpy.data.objects[meshname].data.vertices[vert].groups:
                if group.weight>max_weight:
                    max_weight   = group.weight
                    max_group_id = group.group
            #print(vert)
        #~ print(max_group_id)
        bone_name = bpy.data.objects[meshname].vertex_groups[max_group_id].name
        #~ bone_name = bpy.data.objects[armaturename].pose.bones[max_group_id].name
        #~ print(bone_name)
        flr = []
        for li in poly.loop_indices:
            uv_registry[    bone_name].append(bpy.data.objects[meshname].data.uv_layers.active.data[li].uv)
        #~ print(poly.loop_indices)
        for vert in poly.vertices:
            vertex_registry[bone_name].append(bpy.data.objects[meshname].data.vertices[vert].co)
            flr.append(len(vertex_registry[bone_name])-1)
        #print(bone_name, flr)
        face_local_registry[   bone_name].append(flr)
        face_material_registry[bone_name].append(poly.material_index)
        #print(bone_name)




def rg():
    global joint_registry, joint_bvectors_registry, vertex_registry
    rot = bpy.data.objects[armaturename].rotation_euler
    sc_mesh  = bpy.data.objects[meshname].scale
    sc  = bpy.data.objects[armaturename].scale
    loc = bpy.data.objects[armaturename].location
    # set mesh origin to armature origin
    #~ bpy.context.scene.cursor_location = (loc[0], loc[1], loc[2]) 
    #~ bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    #~ bpy.data.objects[meshname].origin_set(type='ORIGIN_CURSOR')
    #~ print(sc)
    #~ print(loc)
    for eb in bpy.context.scene.objects["Armature"].data.edit_bones:
        eb.roll = 0.0
    additional_objects_to_merge = {}
    for bone_ in bpy.data.objects[armaturename].pose.bones.items():
        # bone location
        vecloc_head = (loc[0] + sc[0]*bone_[1].head[0], 
        loc[1] + math.cos(rot[0]) * sc[1]*bone_[1].head[1] - math.sin(rot[0]) * sc[2]*bone_[1].head[2], 
        loc[2] + math.sin(rot[0]) * sc[1]*bone_[1].head[1] + math.cos(rot[0]) * sc[2]*bone_[1].head[2])
        vecloc_tail = (loc[0] + sc[0]*bone_[1].tail[0], 
        loc[1] + math.cos(rot[0]) * sc[1]*bone_[1].tail[1] - math.sin(rot[0]) * sc[2]*bone_[1].tail[2], 
        loc[2] + math.sin(rot[0]) * sc[1]*bone_[1].tail[1] + math.cos(rot[0]) * sc[2]*bone_[1].tail[2])
        vecdir = ((vecloc_head[0]-vecloc_tail[0]), (vecloc_head[1]-vecloc_tail[1]), (vecloc_head[2]-vecloc_tail[2]))
        vecloc = ((vecloc_head[0]+vecloc_tail[0])*0.5, (vecloc_head[1]+vecloc_tail[1])*0.5, (vecloc_head[2]+vecloc_tail[2])*0.5)
        vecdir_length = math.sqrt(vecdir[0]*vecdir[0] + vecdir[1]*vecdir[1] + vecdir[2]*vecdir[2])
        v_origin = bpy.data.objects[meshname].location
        vertices_ = []
        mass_ratio = 0.0
        #~ vecloc_head[0] - bpy.data.objects["obj_"+bone_[0]].location[0]
        for vert in vertex_registry[bone_[0]]:
            vertices_.append([
            (vert[0]) * sc_mesh[0] + v_origin[0], 
            (vert[1]) * sc_mesh[1] + v_origin[1], 
            (vert[2]) * sc_mesh[2] + v_origin[2]
            ])
            #mass_ratio += norm(bpy.data.objects[meshname].data.vertices[ivid].co)
            #mass_ratio += norm(sub(vertices_[len(vertices_)-1], vecloc))
            mass_ratio += norm(sub(vertices_[len(vertices_)-1], vecloc))**2.0           
        #mass_ratio = math.sqrt(mass_ratio)
        try:
            mass_ratio /= float(len(vertex_registry[bone_[0]])) 
        except:
            mass_ratio = 0.1
        me = bpy.data.meshes.new("mesh_"+bone_[0])
        me.from_pydata(vertices_, [], face_local_registry[bone_[0]])
        me.update()
        # merge double vertices
        #~ me.remDoubles(0.0001)
        #~ print(dir(me))
        #~ bpy.context.active_object.name = 
        ob = bpy.data.objects.new("obj_"+bone_[0], me)
        scn = bpy.context.scene
        scn.objects.link(ob)
        # decimate        
        bpy.context.area.type = "VIEW_3D"
        bpy.context.scene.objects.active = bpy.data.objects["obj_"+bone_[0]]
        bpy.ops.object.editmode_toggle()
        #~ print(bpy.context.area.type)
        try:
            bpy.ops.mesh.remove_doubles()
        except BaseException as e:
            a_tmp = 0
            #~ print(str(e))
        bpy.ops.object.editmode_toggle()
        bpy.context.area.type = "CONSOLE"
        mod = bpy.data.objects["obj_"+bone_[0]].modifiers.new(name='decimate', type='DECIMATE')
        mod.ratio=1.0
        #~ mod.ratio=0.2
        #~ mod.ratio=0.03
        # adding material for friction properties
        mat = bpy.data.materials.new("mat_"+bone_[0])
        mat.physics.friction = 3.0
        ob.data.materials.append(mat)
        # setting physical shape
        ob.game.physics_type = "RIGID_BODY"
        ob.game.use_collision_bounds = True
        #~ ob.game.collision_bounds_type = "TRIANGLE_MESH"
        ob.game.collision_bounds_type = "CONVEX_HULL"
        #~ ob.game.collision_bounds_type = "BOX"
        ob.game.collision_margin = 0.040
        #~ ob.game.collision_margin = 0.010
        #~ ob.game.collision_margin = 0.001
        #~ ob.game.collision_margin = 0.0
        ob.hide_render = True
        if friction_coefficients_!=None:
            ob.game.use_anisotropic_friction = True
            ob.game.friction_coefficients = friction_coefficients_
        ob.game.damping = 0.5
        ob.game.rotation_damping = 0.5
        #obj_mass = vecdir_length * mass_ratio * mass_ratio
        #obj_mass =  mass_ratio * mass_ratio
        obj_mass =  mass_ratio
        #obj_mass =  np.sqrt(mass_ratio)
        #Njoints = len(joint_registry[bone_[0]])
        #Njoints = 1.0
        #if Njoints>=1:obj_mass /= Njoints
        ob.game.mass = obj_mass * 1.0 * mass_scaling
        #ob.game.mass = 1.0
        ob.game.form_factor = 1.0
        ob.game.use_sleep = True
        #ob.game.use_ghost = True
        ob.game.collision_group[1] = True
        ob.game.collision_group[0] = False
        ob.game.collision_mask[0]  = True
        ob.game.collision_mask[1]  = False
        bpy.ops.logic.sensor_add(    type='ALWAYS', object="obj_"+bone_[0])
        bpy.ops.logic.controller_add(type='PYTHON', object="obj_"+bone_[0])
        ob.game.controllers[-1].type   = 'PYTHON'
        ob.game.controllers[-1].states = 1
        ob.game.controllers[-1].text = bpy.data.texts['ghostit.py']
        ob.game.sensors[-1].link(ob.game.controllers[-1])
        bpy.ops.object.mode_set(mode='OBJECT')
        #bpy.context.active_object.name = "obj_"+bone_[0]
        v1 = [-(vecloc_head[0] - vecloc_tail[0]), -(vecloc_head[1] - vecloc_tail[1]), -(vecloc_head[2] - vecloc_tail[2])]
        v1 = div(v1, norm(v1))
        crossed = [-v1[1], v1[0], 0.0]
        quat = [math.sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) + v1[2], crossed[0], crossed[1], crossed[2]]
        quat = div(quat, norm(quat))
        bpy.ops.object.armature_add(location = vecloc_head)
        bpy.context.active_object.name = "CuBone_"+bone_[0]
        bpy.context.object.rotation_mode = "QUATERNION"
        #bpy.context.object.rotation_mode = "XYZ"
        bpy.ops.object.mode_set(mode='EDIT')
        #vecdir2 = [bone_[1].head[0] - bone_[1].tail[0], bone_[1].head[1] - bone_[1].tail[1], bone_[1].head[2] - bone_[1].tail[2]]
        bpy.ops.transform.translate(value=(-vecdir[0]*0.1, -vecdir[1]*0.1, -vecdir[2]*0.1 - 1.0), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        #print(dir(bpy.ops))
        #bpy.ops.transform.translate(value=(0.0,0.0,-1.0), constraint_axis=(False, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
        bpy.ops.object.mode_set(mode='OBJECT')
        #bpy.context.active_object.scale = [0.1,0.1,0.1]
        bpy.data.objects["CuBone_"+bone_[0]].select = True
        bpy.data.objects["obj_"+bone_[0]].select = True
        bpy.context.scene.objects.active = bpy.data.objects["obj_"+bone_[0]]
        bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
        #bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
        for ivert in range(len(vertex_registry[bone_[0]])):
            vertex_registry[bone_[0]][ivert] = sub( add(v_origin, vertex_registry[bone_[0]][ivert]), bpy.data.objects["obj_"+bone_[0]].location)
        bpy.ops.object.parent_set()
        bpy.data.objects["CuBone_"+bone_[0]].select = False
        bpy.data.objects["obj_"+bone_[0]].select = False
        # Constraints
        bone_[1].constraints.new(type="COPY_LOCATION")
        bone_[1].constraints["Copy Location"].target = bpy.data.objects["CuBone_"+bone_[0]]
        bone_[1].constraints["Copy Location"].subtarget = "Bone"
        bone_[1].constraints.new(type="COPY_ROTATION")
        bone_[1].constraints["Copy Rotation"].target = bpy.data.objects["CuBone_"+bone_[0]]
        bone_[1].constraints["Copy Rotation"].subtarget = "Bone"
        if familial_registry[bone_[0]]!="None":
            Njoints = len(joint_registry[bone_[0]])
            if Njoints==0:
                if bone_[1].lock_rotation[0]==True:                
                    joint_registry[bone_[0]].append([ [1.0,0.0,0.0], [0.0, 0.0] ])
                else:
                    joint_registry[bone_[0]].append([ [1.0,0.0,0.0], [0.0, 0.0] ])
                Njoints = 1
            for isj in range(1, Njoints):
                mejoint = bpy.data.meshes.new("joint_mesh_"+bone_[0])
                #mejoint.from_pydata([[0.0,0.0,0.0]], [], [])
                mejoint.from_pydata([vecloc_head], [], [])
                mejoint.update()
                objoint = bpy.data.objects.new("joint_"+str(isj)+"_obj_"+bone_[0], mejoint)
                scn.objects.link(objoint)
                objoint.game.physics_type = "RIGID_BODY"
                objoint.hide_render = True
                objoint.game.damping = 0.5
                objoint.game.rotation_damping = 0.5
                #objoint.game.mass = 0.5 * vecdir_length
                #objoint.game.mass = 0.5 * vecdir_length / float(Njoints)
                objoint.game.mass = obj_mass * 1.0 * mass_scaling
                # WHY ???????
                #objoint.game.mass = 0.01
                objoint.game.form_factor = 0.4
                #objoint.game.use_sleep = True
                #objoint.game.use_ghost = True
                #~ ob.game.collision_margin = 0.04
                ob.game.collision_margin = 0.0
                objoint.game.collision_group[1] = True
                objoint.game.collision_group[0] = False
                #objoint.game.collision_mask[0]  = True
                for i_mask in range(0,8):
                    objoint.game.collision_mask[i_mask]  = False
                bpy.ops.logic.sensor_add(    type='ALWAYS', object="joint_"+str(isj)+"_obj_"+bone_[0])
                bpy.ops.logic.controller_add(type='PYTHON', object="joint_"+str(isj)+"_obj_"+bone_[0])
                objoint.game.controllers[-1].type   = 'PYTHON'
                objoint.game.controllers[-1].states = 1
                objoint.game.controllers[-1].text = bpy.data.texts['ghostit.py']
                #objoint.game.sensors[-1].use_pulse_true_level = True
                objoint.game.sensors[-1].link(objoint.game.controllers[-1])
                bpy.ops.object.mode_set(mode='OBJECT')
                bpy.context.scene.objects.active = bpy.data.objects["joint_"+str(isj)+"_obj_"+bone_[0]]
                bpy.data.objects["joint_"+str(isj)+"_obj_"+bone_[0]].select = True
                #bpy.ops.transform.translate(value=(vecloc_head[0], vecloc_head[1], vecloc_head[2]))
                bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')  
                bpy.data.objects["joint_"+str(isj)+"_obj_"+bone_[0]].select = False              
            for isj in range(Njoints):
                if isj==0:bpy.context.scene.objects.active = bpy.data.objects["obj_"+bone_[0]]
                else:     bpy.context.scene.objects.active = bpy.data.objects["joint_"+str(isj)+"_obj_"+bone_[0]]
                bpy.ops.object.constraint_add(type='RIGID_BODY_JOINT')
                if isj==Njoints-1:bpy.context.object.constraints["Rigid Body Joint"].target = bpy.data.objects["obj_"+familial_registry[bone_[0]]]
                else:             bpy.context.object.constraints["Rigid Body Joint"].target = bpy.data.objects["joint_"+str(isj+1)+"_obj_"+bone_[0]]
                bpy.context.object.constraints["Rigid Body Joint"].show_pivot = True
                if isj==0:
                    pivot_x_ = (vecloc_head[0] - bpy.data.objects["obj_"+bone_[0]].location[0])
                    pivot_y_ = (vecloc_head[1] - bpy.data.objects["obj_"+bone_[0]].location[1])
                    pivot_z_ = (vecloc_head[2] - bpy.data.objects["obj_"+bone_[0]].location[2])
                else:
                    pivot_x_ = 0.0; pivot_y_ = 0.0; pivot_z_ = 0.0; 
                bpy.context.object.constraints["Rigid Body Joint"].pivot_x = pivot_x_
                bpy.context.object.constraints["Rigid Body Joint"].pivot_y = pivot_y_
                bpy.context.object.constraints["Rigid Body Joint"].pivot_z = pivot_z_
                rotv_xy = [joint_registry[bone_[0]][isj][0][0], joint_registry[bone_[0]][isj][0][1]]
                rotv_xy = div(rotv_xy, norm(rotv_xy))
                rotv_  = joint_registry[bone_[0]][isj][0]
                rotv_ = div(rotv_, norm(rotv_))
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
                bpy.context.object.constraints["Rigid Body Joint"].limit_angle_min_x = deg_to_rad(joint_registry[bone_[0]][isj][1][0])
                bpy.context.object.constraints["Rigid Body Joint"].limit_angle_max_x = deg_to_rad(joint_registry[bone_[0]][isj][1][1])






sensory_map         = None
sensory_map_pos     = None
sensory_map_objname = None
sensory_map_isused  = None
sensory_ids_used = None
objstr_ids = {}
def prepare_sensory_map():
    global sensory_map, sensory_map_pos, sensory_map_objname, sensory_map_isused, objstr_ids
    dimZ = 0
    for iflrk, flrk in enumerate(face_material_registry.keys()):
        for iflr,flr in enumerate(face_material_registry[flrk]):
            #~ print(flr+1, dimZ)
            if(flr+1>dimZ):
                dimZ = flr+1    
    sensory_map = np.zeros((dimX, dimY, dimZ), dtype = np.float32)
    sensory_map_pos = np.ones((dimX, dimY, dimZ, 3), dtype = np.float32)
    sensory_map_objname = np.zeros((dimX, dimY, dimZ), dtype = '|S128')
    sensory_map_isused = np.zeros((dimX, dimY, dimZ), dtype = np.int16)
    vertex_array = {}
    tot_polys = 0
    for iflrk, flrk in enumerate(face_local_registry.keys()):
        vertices_tmp = []
        quads_tmp    = []
        uv_tmp       = []
        polypoints_tmp = []
        polyuv_tmp = []
        for iflr,flr in enumerate(face_local_registry[flrk]):
            polypoints_tmp_tmp = []
            polyuv_tmp_tmp     = []
            #for i_vertex_id in range(3):
            for i_vertex_id in range(len(flr)):
                polypoints_tmp_tmp.append(vertex_registry[flrk][flr[i_vertex_id]])
                polyuv_tmp_tmp.append(    uv_registry[    flrk][flr[i_vertex_id]])
            polypoints_tmp.append(polypoints_tmp_tmp)
            polyuv_tmp.append(polyuv_tmp_tmp)
            for vertex_id in flr:
                vertices_tmp.append(vertex_registry[flrk][vertex_id])
                quads_tmp.append(   iflr )
                uv_tmp.append(      uv_registry[flrk][vertex_id] )
        vertex_array[flrk]            = np.array(vertices_tmp)
        vertex_array[flrk+"_polynum"] = np.array(quads_tmp)
        vertex_array[flrk+"_uv"]      = np.array(uv_tmp)
        vertex_array[flrk+"_polypoints"]= polypoints_tmp
        vertex_array[flrk+"_polyuv"]    = polyuv_tmp
        #~ print(vertex_array[flrk].shape)
        #~ print(vertex_array[flrk+"_polynum"].shape)
        #~ print(vertex_array[flrk+"_uv"].shape)
        #~ print("========")
    for flrk in face_local_registry.keys():
        #~ print(len(face_material_registry[flrk]))
        #~ print(len(vertex_array[flrk+"_polyuv"]))
        for ipoly in range(len(vertex_array[flrk+"_polyuv"])):
            #get material (for dimZ)
            idz = face_material_registry[flrk][ipoly]
            ist = np.zeros((len(vertex_array[flrk+"_polyuv"][ipoly]),2), dtype = np.int32)
            for ipolypoint in range(ist.shape[0]):
                id_sens = np.int32(np.array(vertex_array[flrk+"_polyuv"][ipoly][ipolypoint]) * np.float64(sensory_map.shape[:2]))
                ist[ipolypoint,:] = ( id_sens[:] )
            # for all pixels, check if it is in
            min_ix = np.min(ist[:,0])
            min_iy = np.min(ist[:,1])
            max_ix = np.max(ist[:,0])
            max_iy = np.max(ist[:,1])
            if min_ix<0:min_ix = 0
            if min_iy<0:min_iy = 0
            if max_ix>=dimX:max_ix = dimX-1
            if max_iy>=dimY:max_iy = dimY-1
            Tmat = np.array([[ist[0,0]-ist[2,0],ist[1,0]-ist[2,0]],[ist[0,1]-ist[2,1],ist[1,1]-ist[2,1]]])
            Tmat_inv = inverse2x2(Tmat)
            for ix in range(min_ix, max_ix):
                for iy in range(min_iy, max_iy):
                    itisin = False
                    for i_ist in range(2, ist.shape[0]):
                        if point_in_triangle([ix,iy], ist[0], ist[i_ist-1], ist[i_ist]):
                            itisin = True;
                    if itisin==True:
                        lambda_ = matrix_vec_mult_2x2(Tmat_inv, sub([ix, iy], ist[2]))
                        lambda3 = 1.0 -lambda_[0] - lambda_[1]
                        #sensory_map_pos[   ix, iy, :] = np.array(vertex_array[flrk+"_polypoints"][ipoly][0])
                        v1 = np.array(vertex_array[flrk+"_polypoints"][ipoly][0])
                        v2 = np.array(vertex_array[flrk+"_polypoints"][ipoly][1])
                        v3 = np.array(vertex_array[flrk+"_polypoints"][ipoly][2])
                        #~ print(sensory_map_pos.shape)
                        #~ print(dimZ)
                        #~ print(dimY-iy-1, ix, idz)
                        sensory_map_pos[    dimY-iy-1, ix, idz, :] = add( add(mult(v1, lambda_[0]), mult(v2, lambda_[1])), mult(v3, lambda3) )
                        sensory_map_objname[dimY-iy-1, ix, idz]    = flrk
                        sensory_map_isused[ dimY-iy-1, ix, idz]    = 1
    for flrk in face_local_registry.keys():
        objstr_ids[flrk] = np.nonzero(sensory_map_objname==flrk.encode())
    sensory_ids_used = np.nonzero(sensory_map_isused==1)
    sensory_map[:,:,:] = 0.0
    #~ sensory_map[sensory_ids_used[0], sensory_ids_used[1], sensory_ids_used[2]] = 0.0




children_registry = {}

def write_children(f, highest_object, numgen):
    global children_registry, familial_registry
    global joint_registry, joint_bvectors_registry, vertex_registry
    if numgen==0:
        f.write(("  "*numgen)+"Robot{ \n")
        f.write(("  "*numgen)+"rotation 1.0 0.0 0.0 -1.54 \n")
        f.write(("  "*numgen)+"translation 0.0 2.0 0.0 \n")
        f.write(("  "*numgen)+"children[ \n")
    else:
        f.write(("  "*numgen)+"DEF "+highest_object.replace(".", "_")+" HingeJoint{ \n")
        f.write(("  "*numgen)+"device RotationalMotor { \n")
        f.write(("  "*numgen)+"name \""+highest_object+"\" \n")
        f.write(("  "*numgen)+"maxVelocity 4.0 \n")
        f.write(("  "*numgen)+"minPosition -1.0 \n")
        f.write(("  "*numgen)+"maxPosition  1.0 \n")
        f.write(("  "*numgen)+"} \n")
        f.write(("  "*numgen)+"jointParameters HingeJointParameters { \n")
        f.write(("  "*numgen)+"axis 1.0 0.0 0.0  \n")
        pp_ = np.array(bpy.data.objects["obj_"+highest_object].location) - np.array(bpy.data.objects["obj_"+familial_registry[highest_object]].location)
        rbj = bpy.data.objects["obj_"+highest_object].constraints["Rigid Body Joint"]
        pivotpos = pp_ + np.array([rbj.pivot_x, rbj.pivot_y, rbj.pivot_z]) 
        #~ f.write(("  "*numgen)+"anchor 0 -0.0 -0.0  \n")
        f.write(("  "*numgen)+"anchor "+str(pivotpos[0])+" "+str(pivotpos[1])+" "+str(pivotpos[2])+" \n")
        f.write(("  "*numgen)+"} \n")
        f.write(("  "*numgen)+"endPoint Solid { \n")
        f.write(("  "*numgen)+"translation "+str(pp_[0])+" "+str(pp_[1])+" "+str(pp_[2])+" \n")
        #~ f.write(("  "*numgen)+"translation 0 0 0 \n")
        f.write(("  "*numgen)+" children [ \n")
    f.write(("  "*numgen)+"Shape{ \n")
    f.write(("  "*numgen)+"appearance Appearance{ \n")
    f.write(("  "*numgen)+"material Material{ \n")
    f.write(("  "*numgen)+"ambientIntensity 0.4 \n")
    f.write(("  "*numgen)+"diffuseColor 0.5 0.5 0.5 \n")
    f.write(("  "*numgen)+"} \n")
    f.write(("  "*numgen)+"} \n")
    
    # ========================================================================================================================================================================================================
    #~ f.write(("  "*numgen)+"geometry Cylinder { \n")
    #~ f.write(("  "*numgen)+"height 0.4 \n")
    #~ f.write(("  "*numgen)+"radius 0.02 \n")
    #~ f.write(("  "*numgen)+"subdivision 10 \n")
    #~ f.write(("  "*numgen)+"} \n")
    f.write(("  "*numgen)+"geometry IndexedFaceSet { \n")
    f.write(("  "*numgen)+"  coord Coordinate { \n")
    f.write(("  "*numgen)+"    point [ \n")
    for vert in bpy.data.objects["obj_"+highest_object].data.vertices:
        f.write( ("  "*numgen)+str(vert.co[0])+" "+str(vert.co[1])+" "+str(vert.co[2])+"  \n")
    f.write(("  "*numgen)+"    ] \n")
    f.write(("  "*numgen)+"  } \n")
    f.write(("  "*numgen)+"  coordIndex [ \n")
    for poly in bpy.data.objects["obj_"+highest_object].data.polygons:
        verts_ = ""
        for vert in poly.vertices:
            verts_+=str(vert)
            verts_+=" "
        f.write( ("  "*numgen)+verts_+" -1  \n")
    f.write(("  "*numgen)+"  ] \n")
    f.write(("  "*numgen)+"  creaseAngle 2.0  \n")
    f.write(("  "*numgen)+"} \n")
    # ========================================================================================================================================================================================================
    
    f.write(("  "*numgen)+"} \n")   
    for cr in children_registry[highest_object]:
        write_children(f, cr, numgen+1)
    f.write(("  "*numgen)+" ] \n")
    #~ f.write(("  "*numgen)+"boundingObject Transform { \n")
    f.write(("  "*numgen)+"boundingObject Group { \n")
    #~ f.write(("  "*numgen)+"translation 0.0 0.0 0.0 \n")
    f.write(("  "*numgen)+"children [ \n")
    
    # ========================================================================================================================================================================================================
    max_xy = 0.0
    max_z  = 0.0
    ave_xy = 0.0
    ave_z  = 0.0
    for ivert, vert in enumerate(bpy.data.objects["obj_"+highest_object].data.vertices):
        vertpos = np.array(vert.co)
        dist_xy = np.sqrt(vertpos[0]*vertpos[0] + vertpos[2]*vertpos[2])
        dist_z  = np.sqrt(vertpos[1]*vertpos[1])
        ave_xy += dist_xy
        ave_z  += dist_z
        if dist_xy>max_xy:max_xy = dist_xy;
        if dist_z >max_z :max_z  = dist_z;
    ave_xy /= float(len(bpy.data.objects["obj_"+highest_object].data.vertices))
    ave_z  /= float(len(bpy.data.objects["obj_"+highest_object].data.vertices))
    f.write(("  "*numgen)+"Cylinder{ \n")
    f.write(("  "*numgen)+"height "+str(max_z)+" \n")
    #~ f.write(("  "*numgen)+"radius "+str(max_xy*0.6)+" \n")
    #~ f.write(("  "*numgen)+"height "+str(ave_z)+" \n")
    f.write(("  "*numgen)+"radius "+str(ave_xy)+" \n")
    f.write(("  "*numgen)+"subdivision 8 \n")
    f.write(("  "*numgen)+"} \n")
    #~ f.write(("  "*numgen)+"Shape { \n")
    #~ f.write(("  "*numgen)+"geometry IndexedFaceSet { \n")
    #~ f.write(("  "*numgen)+"  coord Coordinate { \n")
    #~ f.write(("  "*numgen)+"    point [ \n")
    #~ 
    #~ #f.write(("  "*numgen)+"      -0.000088 -0.00000 0.000047    \n")
    #~ #f.write(("  "*numgen)+"      -0.014830 -0.009771 0.013278  \n")
    #~ #f.write(("  "*numgen)+"      -0.011830 -0.019771 0.023278  \n")
    #~ #f.write(("  "*numgen)+"      -0.010830 0.009771 -0.013278  \n")
    #~ print(len(bpy.data.objects["obj_"+highest_object].data.vertices))
    #~ for ivert, vert in enumerate(bpy.data.objects["obj_"+highest_object].data.vertices):
        #~ if ivert<4:
            #~ #f.write( ("  "*numgen)+str(vert.co[0])+" "+str(vert.co[1])+" "+str(vert.co[2])+"  \n")
            #~ f.write( ("  "*numgen)+"      "+('%f' % vert.co[0])+" "+str('%f' % vert.co[1])+" "+str('%f' % vert.co[2])+"  \n")
    #~ 
    #~ f.write(("  "*numgen)+"    ] \n")
    #~ f.write(("  "*numgen)+"  } \n")
    #~ f.write(("  "*numgen)+"  coordIndex [ \n")
    #~ #for poly in bpy.data.objects["obj_"+highest_object].data.polygons:
        #~ #verts_ = ""
        #~ #vertcount = 0
        #~ #for vert in poly.vertices:
            #~ #if vertcount<=2:
                #~ #verts_+=str(vert)
                #~ #verts_+=" "
            #~ #vertcount += 1
        #~ #f.write( ("  "*numgen)+verts_+" -1  \n")
    #~ 
    #~ f.write(("  "*numgen)+"    0 1 2 -1  \n")
    #~ #f.write(("  "*numgen)+"    1 2 3 -1  \n")
    #~ 
    #~ f.write(("  "*numgen)+"  ] \n")
    #~ #f.write(("  "*numgen)+"  creaseAngle 2.0  \n")
    #~ f.write(("  "*numgen)+"} \n")    
    #~ f.write(("  "*numgen)+"} \n")    
    #~ # ========================================================================================================================================================================================================
    f.write(("  "*numgen)+"] \n")
    f.write(("  "*numgen)+"} \n")
    f.write(("  "*numgen)+"physics Physics { \n")
    f.write(("  "*numgen)+"density -1 \n")
    f.write(("  "*numgen)+"mass "+str(bpy.data.objects["obj_"+highest_object].game.mass)+" \n")
    #~ f.write(("  "*numgen)+"mass 0.1 \n")
    f.write(("  "*numgen)+"} \n")
    #~ f.write(("  "*numgen)+"} \n")
    #~ if numgen!=0:
        #~ f.write(("  "*numgen)+"maxForce 110 \n")
        #~ f.write(("  "*numgen)+"controlP 100 \n")
    if numgen==0:
        f.write(("  "*numgen)+"controller \"defaultCsX\" \n")
    else:
        f.write(("  "*numgen)+"} \n")
    f.write(("  "*numgen)+"} \n")


def rg2webots():
    global children_registry
    print("Exporting to webots")
    f = open(export2webots, 'w')
    webots_intro = '''#VRML_SIM V6.0 utf8 \nWorldInfo { \n  info [ \n    "Blender Export Robot" \n    "Authors: Csaba Eroe" \n    "Date: Sept 29, 2014" \n  ] \n  title "Blender2Webots" \n  ERP 0.6 \n  basicTimeStep 2.5 \n  lineScale 1 \n} \nViewpoint { \n  orientation 0.0702507 0.980637 0.182799 3.65007 \n  position -1.33071 1.7029 -3.23726 \n} \nBackground { \n} \nDirectionalLight { \n  ambientIntensity 0.9 \n  direction 0.4 -0.8 -0.2 \n  intensity 0.8 \n  castShadows TRUE \n} \n\n\n '''
    f.write(webots_intro)
    webots_start = ''' \n \n'''
    f.write(webots_start)
    global joint_registry, joint_bvectors_registry, vertex_registry
    # get: children registry, highest order object
    highest_object = ""
    for bone_ in bpy.data.objects[armaturename].pose.bones.items():
        children_registry[bone_[0]] = []
    for bone_ in bpy.data.objects[armaturename].pose.bones.items():
        if familial_registry[bone_[0]]=="None":
            highest_object = bone_[0]
        else:
            if bone_[0].find("Whisker")<0:
                children_registry[familial_registry[bone_[0]]].append(bone_[0])
                #~ print(familial_registry[bone_[0]], bone_[0])
    print("Highest parent: ", highest_object)
    write_children(f, highest_object, 0)
    #~ for bone_ in bpy.data.objects[armaturename].pose.bones.items():
        #~ for poly in bpy.data.objects["obj_"+bone_[0]].data.polygons:
            #~ for vert in poly.vertices:
                #~ f.write( str(bpy.data.objects["obj_"+bone_[0]].data.vertices[vert].co[0]) +", "+str(bpy.data.objects["obj_"+bone_[0]].data.vertices[vert].co[1]) +", "+str(bpy.data.objects["obj_"+bone_[0]].data.vertices[vert].co[2])+"  \n"  )
                #~ flr.append(len(vertex_registry[bone_name])-1)
            #~ face_local_registry[   bone_name].append(flr)
        #~ f.write("--------------------------------------------- \n\n")
    webots_outro = ''' \nDEF FLOOR Solid { \n  children [ \n    Transform { \n      translation -4 0 -4 \n      children [ \n        Shape { \n          appearance Appearance { \n            texture ImageTexture { \n              url [ \n                "wood3.jpg" \n              ] \n            } \n            textureTransform TextureTransform { \n              scale 5 5 \n            } \n          } \n          geometry ElevationGrid { \n            colorPerVertex FALSE \n            xDimension 10 \n            zDimension 10 \n          } \n        } \n      ] \n    } \n  ] \n  boundingObject Plane { \n    size 300 300 \n  } \n  locked TRUE \n}'''
    f.write(webots_outro)
    f.close()





init()

def build():
    rg()
    ew_ = False
    try:
        export2webots
        ew_ = True
    except:
        print("No webots export")
    if ew_:rg2webots()
    prepare_sensory_map()
