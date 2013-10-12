from cffi import FFI
ffi = FFI()
ffi.cdef("""\
typedef size_t p_t;

int _H5Fopen(const char *name, char mode);
int _H5Fclose(int f);

int _H5Dopen(int f, const char *name);
int _H5Dclose(int d);

int _ls(int f, void *cb);

int _len(int d);
void _read(int d, p_t buf);

""")

H5 = ffi.verify("""\
typedef size_t p_t;

#include <hdf5.h>

int _H5Fopen(const char *name, char mode) {
    hid_t f;
    switch (mode) {
        case 'a': if ((f=H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT)) > 0) return f;
        case 'w': return H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        default: return H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    }
}

int _H5Fclose(int f) {
    return H5Fclose(f);
}

typedef void (callback)(const char *);

herr_t op(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data) {
    callback *cb = (callback*) op_data;
    cb(name);
    return 0;
}

int _ls(int f, void *cb) {
    H5Lvisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, op, cb);
    return 0;
}

int _H5Dopen(int f, const char *name) {
    return H5Dopen(f, name, H5P_DEFAULT);
}

int _H5Dclose(int d) {
    return H5Dclose(d);
}

int _len(int d) {
    hid_t s = H5Dget_space(d);
    int n = H5Sget_simple_extent_npoints(s);
    H5Sclose(s);
    return n;
}

void _read(int d, p_t buf) {
    hid_t t = H5Dget_type(d);
    H5Dread(d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*) buf);
    H5Tclose(t);
}

""", libraries=['hdf5'])


try:
    import numpy as np
except ImportError:
    import numpypy as np

def get_ptr_array(x):
    return x.__array_interface__['data'][0]

try:
    from accel import _get_ptr as get_ptr
except ImportError:
    get_ptr = get_ptr_array


class File:
    def __init__(self, name, mode):
        self.name, self.mode, self.id = name, mode, None

    def __enter__(self):
        self.id = H5._H5Fopen(self.name, self.mode)
        return self

    def __exit__(self, *args):
        if self.id is not None:
            H5._H5Fclose(self.id)

    def __getitem__(self, name):
        return Dataset(self.id, name)

    def __iter__(self):
        nodes = []
        @ffi.callback("void(const char*)")
        def python_callback(name):
            nodes.append(ffi.string(name))

        H5._ls(self.id, python_callback)

        for node in nodes:
            yield node

    def ls(self):
        return [node for node in self]


class Dataset:
    def __init__(self, fid, name):
        self.fid, self.name = fid, name

    def len(self):
        d = H5._H5Dopen(self.fid, self.name)
        n = H5._len(d)
        H5._H5Dclose(d)
        return n

    @property
    def value(self):
        d = H5._H5Dopen(self.fid, self.name)
        n = H5._len(d)
        a = np.empty(n, 'f')
        n = H5._read(d, get_ptr(a))
        H5._H5Dclose(d)
        return a


if __name__ == "__main__":
    with File('29864_XPR.h5', 'r') as f:
        d = f['/S4']
        print d.len()
        a = d.value
        print a[:20]

        print f.ls()



