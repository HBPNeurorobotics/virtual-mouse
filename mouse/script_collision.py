def collision_function_builder(name_):
    def cb(other, px, py, pz, impulse_):
        global scn, sensory_map
        #print(impulse_)
        #print(name_, "collided with", other)
        #print(other.worldPosition)
        #print(scn.objects[name_].worldPosition)
        #print(px, py, pz)
        clashpos = [px, py, pz]
        sc_cross = impulse_*10.0
        bge.render.drawLine(add(clashpos, [0.0,0.0,-sc_cross]), add(clashpos, [0.0,0.0,sc_cross]), [1.0,1.0,0.0])
        bge.render.drawLine(add(clashpos, [0.0,-sc_cross,0.0]), add(clashpos, [0.0,sc_cross,0.0]), [1.0,1.0,0.0])
        bge.render.drawLine(add(clashpos, [-sc_cross,0.0,0.0]), add(clashpos, [sc_cross,0.0,0.0]), [1.0,1.0,0.0])
        #clashpos_local = np.array( newref( sub(clashpos, scn.objects[name_].worldPosition), scn.objects[name_].worldOrientation.inverted()) )
        
        clashpos_local = np.array( newref( sub(clashpos, scn.objects[name_].worldPosition), scn.objects[name_].worldOrientation.inverted()) )
        distances = np.sqrt(np.sum((sensory_map_pos[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], objstr_ids[name_[4:]][2], :] - clashpos_local)**2.0, axis=1))
        pressures = np.exp( -distances*distances/0.04 )
        #~ pressures = np.exp( -distances*distances/0.6 )
        #~ pressures = 10.0 * impulse_ * np.exp( -distances*distances/0.06 )
        #~ pressures = 60.0 * impulse_ * np.exp( -distances*distances/0.06 )
        sensory_map[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], objstr_ids[name_[4:]][2]] += pressures[:]
        #~ sensory_map[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], objstr_ids[name_[4:]][2]] += 1.0
        #~ print("->", distances)
        
        #~ if pressures.shape[0]>0:
        #    max_pressure_id = np.argmax(pressures[:])
        #    sensory_map[objstr_ids[name_[4:]][0][max_pressure_id], objstr_ids[name_[4:]][1][max_pressure_id] ,0] = 1.0
        
        #for igg in range(len(sensory_map_pos[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], :])):
        #    bge.render.drawLine(add(sensory_map_pos[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], igg], [-1.0,0.0,0.0]), add(sensory_map_pos[objstr_ids[name_[4:]][0], objstr_ids[name_[4:]][1], igg], [1.0,0.0,0.0]), [1.0,1.0,0.0])
        
    return cb
