
def browse_object(obj, fun, name="", visited=None, **kw):
    if visited is None:
        visited = dict()
    obj_id = id(obj)
    if visited.has_key(obj_id):
        return
    else:
        fun(name, obj, **kw)
        visited[obj_id] = obj
    if isinstance(obj, dict):
        for k, v in obj.iteritems():
            fmt = "%s['%s']" if isinstance(k, str) else "%s[%s]"
            browse_object(v, fun, name=fmt % (name, k), visited=visited, **kw)
    elif hasattr(obj, '__dict__'):
        for k, v in obj.__dict__.iteritems():
            browse_object(v, fun, name="%s.%s" % (name, k), visited=visited, **kw)
    elif isinstance(obj, (list, tuple)) or \
            isinstance(obj, np.ndarray) and obj.dtype == object:
        for k, v in enumerate(obj):
            browse_object(v, fun, name="%s[%s]" % (name, k), visited=visited, **kw)
    return visited


import numpy as np

fact = 1. / (1024*1024)

def fun_ndarray_owning(name, obj, info):
    if isinstance(obj, np.ndarray):
        while obj.base is not None:
            obj = obj.base
        obj_id = id(obj)
        visited_owning = info['visited_owning']
        if visited_owning.has_key(obj_id):
            return
        else:
            print "%8.3f MB: %s" % (fact * obj.nbytes, name)
            info['nbytes_tot'] += obj.nbytes
            visited_owning[obj_id] = obj

def ndarray_owning(obj):
    info = dict(visited_owning=dict(), nbytes_tot=0.)
    browse_object(obj, fun_ndarray_owning, info=info)
    print "%8.3f MB total" % (fact * info['nbytes_tot'])
    return info


