from ..containers.nucleic_acid import NucleicAcid
import numpy as np
from itertools import product
from typing import Union, Tuple, List



FRAGMENT_COUNT_CONFIG = {
    
    'helixes':{
        'get_fragment_func':lambda na: na.helixes, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 8), (7, 9), (8, np.inf)]
            }
        }
    }, 
    
    'hairpins':{
        'get_fragment_func':lambda na: na.hairpins, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(2, 6), (4, 7), (6, 8), (7, 9), (8, 10), (9, 11), 
                          (10, 13), (11, 17), (13, 21), (17, 25), (21, 30), (25, np.inf)]
            }
        }
    }, 
    
    'internal_loops':{
        'get_fragment_func':lambda na: na.internal_loops, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(6, 8), (7, 9), (8, 10), (9, 11), (10, 12), 
                          (11, 13), (12, 14), (13, 19), (14, 24), (19, np.inf)]
            }, 
            'asymmetry':{
                'get_property_func':lambda fragment: (fragment.short_side/fragment.long_side),
                'ranges':[(0.00325, 0.33333), (0.16829, 0.41667), 
                          (0.33333, 0.5), (0.41667, 0.625), (0.5, 0.75), 
                          (0.625, 0.875), (0.75, 1.0), (0.875, np.inf)]
            }
        }
    }, 
    
    'bulges':{
        'get_fragment_func':lambda na: na.bulges, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(5, 7), (6, 9), (7, 12), (9, 13), (12, 14), (13, np.inf)]
            }
        }
    }, 
    
    'junctions':{
        'get_fragment_func':lambda na: na.junctions, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(6, 10), (8, 11), (10, 12), (11, 13), (12, 14), 
                          (13, 15), (14, 16), (15, 17), (16, 19), (17, 21), 
                          (19, 23), (21, 28), (23, 33), (28, 43), (33, 54), (43, np.inf)]
            }, 
            'branches':{
                'get_property_func':lambda fragment: len(fragment.branches),
                'ranges':[(3, 5), (4, np.inf)]
            }
        }
    }, 
    
    'dangling_ends':{
        'get_fragment_func':lambda na: na.dangling_ends, 
        'properties':{
            'len':{
                'get_property_func':lambda fragment: len(fragment),
                'ranges':[(0, 2), (1, 3), (2, 4), (3, 5), (4, 7), (5, 9), 
                          (7, 12), (9, 17), (12, 23), (17, 25), (23, 28), (25, np.inf)]
            }
        }
    }, 
    
    # Knots
    'hairpins_knot_ratio':{
        'get_fragment_func':lambda na: [l for l in na.hairpins if l.knots is not None], 
        'properties':{
            'len_ratio':{
                'get_property_func':lambda fragment: sum([len(k) for k in fragment.knots])/len(fragment),
                'ranges':[(0.04762, 0.27778), (0.1627, 0.33333), (0.27778, 0.38889), 
                          (0.33333, 0.42778), (0.38889, 0.46667), (0.42778, 0.50256), 
                          (0.46667, 0.53846), (0.50256, 0.55495), (0.53846, 0.57143), 
                          (0.55495, 0.61905), (0.57143, 0.66667), (0.61905, 0.68333), 
                          (0.66667, 0.7), (0.68333, np.inf)]
            }
        }
    }, 
    
    'internal_loops_knot_ratio':{
        'get_fragment_func':lambda na: [l for l in na.internal_loops if l.knots is not None], 
        'properties':{
            'len_ratio':{
                'get_property_func':lambda fragment: sum([len(k) for k in fragment.knots])/len(fragment),
                'ranges':[(0.07965, 0.125), (0.10232, 0.25694), (0.125, 0.38889), 
                          (0.25694, 0.40873), (0.38889, 0.42857), (0.40873, np.inf)]
            }
        }
    }, 
    
    'bulges_knot_ratio':{
        'get_fragment_func':lambda na: [l for l in na.bulges if l.knots is not None], 
        'properties':{
            'len_ratio':{
                'get_property_func':lambda fragment: sum([len(k) for k in fragment.knots])/len(fragment),
                'ranges':[(0.01587, 0.2), (0.10794, 0.26667), (0.2, 0.33333), 
                          (0.26667, 0.36111), (0.33333, 0.38889), (0.36111, 0.46368), 
                          (0.38889, 0.53846), (0.46368, np.inf)]

            }
        }
    }, 
    
    'junctions_knot_ratio':{
        'get_fragment_func':lambda na: [l for l in na.junctions if l.knots is not None], 
        'properties':{
            'len_ratio':{
                'get_property_func':lambda fragment: sum([len(k) for k in fragment.knots])/len(fragment),
                'ranges':[(0.03704, 0.05263), (0.04483, 0.07587), (0.05263, 0.0991), 
                          (0.07587, 0.13411), (0.0991, 0.16912), (0.13411, 0.1773), 
                          (0.16912, 0.18548), (0.1773, 0.19444), (0.18548, 0.20339), 
                          (0.19444, 0.2174), (0.20339, 0.2314), (0.2174, 0.33309), 
                          (0.2314, 0.43478), (0.33309, np.inf)]
            }
        }
    }, 
    
}

FIRST_KNOT_FRAGMENT_INDEX = 6

FEATURE_INDEXES = {}
for i, fragment_dict in enumerate(FRAGMENT_COUNT_CONFIG.values()):
    all_features = product(*[p['ranges'] for p in fragment_dict['properties'].values()])
    for f in all_features:
        FEATURE_INDEXES[((i, ) + f)] = len(FEATURE_INDEXES)


def get_na_features(na):
    features = []
    
    # for each fragment type
    for i, fragment_dict in enumerate(FRAGMENT_COUNT_CONFIG.values()):
        fragments = fragment_dict['get_fragment_func'](na)
        if len(fragments)==0:
            continue
            
        properties = fragment_dict['properties']
        # for each fragment in NucleicAcid
        for f in fragments:
            
            # for each property of that fragment
            properties_range_match = {}
            for property_name, property_dict in properties.items():
                properties_range_match[property_name] = []
                p = property_dict['get_property_func'](f)
                
                for l, r in property_dict['ranges']:
                    if p>=l and p<r:
                        properties_range_match[property_name].append((l, r))
                        
            for f in product(*list(properties_range_match.values())):
                features.append((i,)+f)
                
    return features

def FragmentCount(na: NucleicAcid, 
                   with_knot_features: bool = True, 
                   return_features: bool = False
                  ) -> Union[np.ndarray, Tuple[np.ndarray, List]]:
    """
    Computes fragment count descriptor vector with 192 features.

    :param na: NucleicAcid object.
    :param with_knot_features: whether to count knot features. Default - True.
    :param return_features: return descriptor and features list.

    :return: descriptor numpy vector.
    """
    
    if not isinstance(na, NucleicAcid):
        raise ValueError("NA argument must be NucleicAcid object.")
        
    descriptor = np.zeros(len(FEATURE_INDEXES), dtype=np.int32)
    features = get_na_features(na)
    if not with_knot_features:
        features = list(filter(lambda x: x[0]<FIRST_KNOT_FRAGMENT_INDEX, features))
        
    # remove 3'end dangling end feature for linear structure
    if len(na.pairs)==0:
        dangling_end_first_range = FRAGMENT_COUNT_CONFIG["dangling_ends"]["properties"]["len"]["ranges"][0]
        dangling_end_fragment_idx = [i for i, f in enumerate(FRAGMENT_COUNT_CONFIG.keys()) if f=='dangling_ends'][0]
        dangling_end_feature = dangling_end_fragment_idx, dangling_end_first_range
        features = list(filter(lambda x: x!=dangling_end_feature, features))
        
        if len(na)==1: # add one end feature
            features.append(dangling_end_feature)

    for f in features:
        idx = FEATURE_INDEXES[f]
        descriptor[idx] += 1
        
    if return_features:
        return descriptor, features
    
    return descriptor


def FragmentFingerprint(na: NucleicAcid, 
                         with_knot_features: bool = True
                        ) -> np.ndarray:
    """
    Computes fragment descriptor with binary values instead of counts.
    """
    
    return (FragmentCount(na, with_knot_features) > 0).astype(np.int32)










