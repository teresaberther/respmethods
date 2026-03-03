import numpy as np
import numba as nb
import tomllib as toml 
from itertools import product

def create_dict_from_toml(fname):
    
    cfg = {}
    ncombs_counter = 1
    with open(fname, 'rb') as f:
        fc = toml.load(f)        
        for cfg_name, cfg_values in fc.items():            
            if "fixed_range" in cfg_values.keys():            
                cfg.update({cfg_name : np.array(cfg_values["fixed_range"])})
            else:
                fr = np.linspace(cfg_values["strtval"],
                                 cfg_values["endval"],
                                 num=cfg_values["n_entries"])
                cfg.update({cfg_name : fr})
            ncombs_counter *= len(cfg[cfg_name])
        print(f"\n{ncombs_counter} combinations will be created. \n")
    return cfg

#%% interim: numba needs a preliminary dictionary example
fc = create_dict_from_toml("sim-cfg.toml")

# generation of all the parameters combination
combinations = product(*fc.values())
results = [dict(zip(fc.keys(), c)) for c in combinations]
tmp_dict = results[0]

# class definition
nb_d = nb.typed.Dict()
for k,v in tmp_dict.items():
    nb_d[k] = v


# key and value types
# kv_ty = (nb.types.unicode_type, nb.types.float64)
# nb_class_specs = [("cfg", nb.types.DictType(*kv_ty)),
#                   ]

nb_class_specs = [("cfg", nb.typeof(nb_d)),
                  ]


@nb.experimental.jitclass(nb_class_specs)
class SimUnity(object):
    def __init__(self, cfg):
        self.cfg = cfg


#%% testing pool

fc = create_dict_from_toml("sim-cfg.toml")

# generation of all the parameters combination
combinations = product(*fc.values())
results = [dict(zip(fc.keys(), c)) for c in combinations]

cnt = SimUnity(results[0])




