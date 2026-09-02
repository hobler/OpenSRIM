"""Target related definitions and operations.

Currently, only planar targets are supported.

Available functions:
    get_distance_from_surface: get the distance of a given position from the
        surface, positive if inside, negative if outside
    check_exit: check if the projectile exits the target and take appropriate 
        action
    is_inside_target: check if a given position is inside the target (obsolescent)
    get_layer_index: get the layer index for a given position
    get_element_index: randomly select an element index for a given projectile's 
        layer index
"""
import numpy as np
from numba import jit


@jit(inline = "always")
def is_inside_target(pos, params):
    """Check if a given position is inside the target.

    If the position is within roughness of the top surface, the result is 
    determined by a random number. The probability of being inside is 
    proportional to the distance from the surface.

    Parameters:
        pos (ndarray): Position to check (size 3)
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (bool): True if the position is inside the target, False otherwise
    """
    if pos[0] < params.geometry.x_intf[0]:
        return False
    elif pos[0] < params.geometry.x_intf[0] + params.geometry.roughness:
        dist_surf = pos[0] - params.geometry.x_intf[0]
        return dist_surf > params.geometry.roughness * np.random.rand()
    elif pos[0] > params.geometry.x_intf[params.geometry.nlayers]:
        return False
    else:
        return True
    return (params.geometry.x_intf[0] <= pos[0] 
            <= params.geometry.x_intf[params.geometry.nlayers])


@jit(inline = "always")
def get_distance_from_surface(pos, params):
    """Get the distance of a given position from the surface.

    The distance is positive if the position is inside, and negative if 
    it is outside the target.
    
    Parameters:
        pos (ndarray): Position to check (size 3)
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (float): Distance to the surface
        (bool): True if pos is on the beam side of the target
    """
    if pos[0] < 0.5 * (params.geometry.x_intf[0] 
                       + params.geometry.x_intf[params.geometry.nlayers]):
        return pos[0] - params.geometry.x_intf[0], True
    else:
        return params.geometry.x_intf[params.geometry.nlayers] - pos[0], False


@jit#(inline = "always")
def check_exit_and_move(proj, free_path, params):
    """Check if the projectile exits the target and move it accordingly.

    proj is modified in-place to reflect its motion along the free flight path.
    In addition, refraction or reflection is applied if the edge of the surface 
    layer is reached.

    The edge of the surface layer is defined as the position where the distance
    from the surface is equal to -params.cascade.pmax_max.

    Parameters:
        proj (Projectile): Projectile to check (modified in-place)
        free_path (float): Free path length
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (bool): True if the projectile exits the target, False otherwise
    """
    pos_new = proj["pos"] + free_path * proj["dir"]
    dist_surf, is_on_beamside = get_distance_from_surface(pos_new, params)

    if dist_surf < -params.cascade.pmax_max:
        # edge of surface layer reached
        if is_on_beamside == proj["is_on_beamside"]:  # common case
            factor = ((-params.cascade.pmax_max - proj["dist_surf"]) 
                      / (dist_surf - proj["dist_surf"]))
        else:  # do it the long way
            pos_plane = (params.geometry.x_intf[0] - params.cascade.pmax_max 
                         if is_on_beamside 
                         else params.geometry.x_intf[params.geometry.nlayers] 
                              + params.cascade.pmax_max)
            factor = (pos_plane - proj["pos"][0]) / (pos_new[0] - proj["pos"][0])
        proj["pos"] += factor * free_path * proj["dir"]
        if proj["pos"][0] < (params.geometry.x_intf[0] 
                             - params.cascade.pmax_max - 0.01):
            print("x=", proj["pos"][0], ",dir=", proj["dir"])
        proj["dist_surf"] = -params.cascade.pmax_max
        proj["is_on_beamside"] = is_on_beamside
        proj["is_inside"] = False
        # Refraction or reflection?
        # alpha: angle wrt surface normal before refraction
        # beta: angle wrt surface normal after refraction
        cos_alpha = np.abs(proj["dir"][0])
        e_perp = proj["e"] * cos_alpha**2
        imat = 0 if is_on_beamside else params.geometry.nlayers - 1
        proj["ilayer"] = imat
        ielem_mat = params.materials[imat].ielem_mat[proj["ielem"]]
        e_surf = params.materials[imat].esurf[ielem_mat]
        if e_perp > e_surf:
            # refraction
            cos_beta = np.sqrt((e_perp - e_surf) / (proj["e"] - e_surf))
            sin_beta = np.sqrt(1.0 - cos_beta**2)
            proj["dir"][0] = np.sign(proj["dir"][0]) * cos_beta
            len_dirxy = np.linalg.norm(proj["dir"][1:])
            if len_dirxy > 0:
                proj["dir"][1:] *= sin_beta / len_dirxy
            proj["e"] -= e_surf
            #print("exiting check_exit_and_move (refract)...")
            return True
        else:
            # reflection
            proj["dir"][0] *= -1
            #print("exiting check_exit_and_move (reflect)...")
            return False
    else:
        # edge of surface layer not reached, just move the projectile
        proj["pos"] = pos_new
        proj["dist_surf"] = dist_surf
        proj["is_on_beamside"] = is_on_beamside
        proj["ilayer"] = get_layer_index(proj["pos"], params)
        proj["is_inside"] = is_inside_target(proj["pos"], params)
        #print("exiting check_exit_and_move (move)...")
        return False


@jit(inline = "always")
def get_layer_index(pos, params):
    """Get the layer index for a given position.

    For pos[0] < params.geometry.x_intf[0] return 0.
    For pos[0] >= params.geometry.x_intf[-1], return the last layer index.

    Note that material index = layer index.

    Parameters:
        pos (ndarray): Position to check (size 3)
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (int): Layer index
    """
    for i in range(1, params.geometry.nlayers):
        if pos[0] < params.geometry.x_intf[i]:
            return i - 1    # return layer index to the left of the interface
        
    return params.geometry.nlayers - 1  # If not found in any layer, 
                                        # return last layer index


@jit(inline = "always")
def get_element_index(proj, params):
    """Randomly select an element index for a given projectile's layer index.

    Parameters:
        proj (Projectile): Projectile to get and set the element index for
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (int): element index
    """
    imat = proj["ilayer"]
    nelem_mat = params.materials[imat].nelem_mat

    r = np.random.rand() * sum(
        params.materials[imat].atomic_fraction[:nelem_mat])
    for ielem_mat in range(nelem_mat - 1):
        if r < params.materials[imat].cumulative_fraction[ielem_mat]:
            return params.materials[imat].ielem[ielem_mat]

    return params.materials[imat].ielem[nelem_mat - 1]