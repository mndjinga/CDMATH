# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_medenum', [dirname(__file__)])
        except ImportError:
            import _medenum
            return _medenum
        if fp is not None:
            try:
                _mod = imp.load_module('_medenum', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _medenum = swig_import_helper()
    del swig_import_helper
else:
    import _medenum
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


ABSOLUTE_H5IPUBLIC_H = _medenum.ABSOLUTE_H5IPUBLIC_H
ABSOLUTE_H5PUBLIC_H = _medenum.ABSOLUTE_H5PUBLIC_H
HAVE_CC_C99 = _medenum.HAVE_CC_C99
HAVE_CUSERID = _medenum.HAVE_CUSERID
HAVE_DLFCN_H = _medenum.HAVE_DLFCN_H
HAVE_FTIME = _medenum.HAVE_FTIME
HAVE_GETEUID = _medenum.HAVE_GETEUID
HAVE_GETPWUID = _medenum.HAVE_GETPWUID
HAVE_GETTIMEOFDAY = _medenum.HAVE_GETTIMEOFDAY
HAVE_H5IPUBLIC_H = _medenum.HAVE_H5IPUBLIC_H
HAVE_H5PUBLIC_H = _medenum.HAVE_H5PUBLIC_H
HAVE_INTTYPES_H = _medenum.HAVE_INTTYPES_H
HAVE_LIBHDF5 = _medenum.HAVE_LIBHDF5
HAVE_MALLOC_H = _medenum.HAVE_MALLOC_H
HAVE_MEMORY_H = _medenum.HAVE_MEMORY_H
HAVE_PWD_H = _medenum.HAVE_PWD_H
HAVE_PYTHON = _medenum.HAVE_PYTHON
HAVE_STDBOOL_H = _medenum.HAVE_STDBOOL_H
HAVE_STDINT_H = _medenum.HAVE_STDINT_H
HAVE_STDLIB_H = _medenum.HAVE_STDLIB_H
HAVE_STRINGS_H = _medenum.HAVE_STRINGS_H
HAVE_STRING_H = _medenum.HAVE_STRING_H
HAVE_SYS_STAT_H = _medenum.HAVE_SYS_STAT_H
HAVE_SYS_TIMEB_H = _medenum.HAVE_SYS_TIMEB_H
HAVE_SYS_TIME_H = _medenum.HAVE_SYS_TIME_H
HAVE_SYS_TYPES_H = _medenum.HAVE_SYS_TYPES_H
HAVE_UNISTD_H = _medenum.HAVE_UNISTD_H
HAVE__BOOL = _medenum.HAVE__BOOL
LT_OBJDIR = _medenum.LT_OBJDIR
MED_API_23 = _medenum.MED_API_23
MED_CHECK_23FORMAT = _medenum.MED_CHECK_23FORMAT
MED_HAVE_FORTRAN = _medenum.MED_HAVE_FORTRAN
MED_HAVE_PYTHON = _medenum.MED_HAVE_PYTHON
MESGERR = _medenum.MESGERR
PACKAGE = _medenum.PACKAGE
PACKAGE_BUGREPORT = _medenum.PACKAGE_BUGREPORT
PACKAGE_NAME = _medenum.PACKAGE_NAME
PACKAGE_STRING = _medenum.PACKAGE_STRING
PACKAGE_TARNAME = _medenum.PACKAGE_TARNAME
PACKAGE_URL = _medenum.PACKAGE_URL
PACKAGE_VERSION = _medenum.PACKAGE_VERSION
SIZEOF_FORTRAN_INTEGER = _medenum.SIZEOF_FORTRAN_INTEGER
SIZEOF_INT = _medenum.SIZEOF_INT
SIZEOF_LONG = _medenum.SIZEOF_LONG
STDC_HEADERS = _medenum.STDC_HEADERS
TIME_WITH_SYS_TIME = _medenum.TIME_WITH_SYS_TIME
VERSION = _medenum.VERSION
H5F_LIBVER_18 = _medenum.H5F_LIBVER_18
MED_MAJOR_NUM = _medenum.MED_MAJOR_NUM
MED_MINOR_NUM = _medenum.MED_MINOR_NUM
MED_RELEASE_NUM = _medenum.MED_RELEASE_NUM
MED_NUM_MAJEUR = _medenum.MED_NUM_MAJEUR
MED_NUM_MINEUR = _medenum.MED_NUM_MINEUR
MED_NUM_RELEASE = _medenum.MED_NUM_RELEASE
MED_VERSION_STR = _medenum.MED_VERSION_STR
MED_MAX_PARA = _medenum.MED_MAX_PARA
MED_COMMENT_SIZE = _medenum.MED_COMMENT_SIZE
MED_IDENT_SIZE = _medenum.MED_IDENT_SIZE
MED_NAME_SIZE = _medenum.MED_NAME_SIZE
MED_SNAME_SIZE = _medenum.MED_SNAME_SIZE
MED_LNAME_SIZE = _medenum.MED_LNAME_SIZE
MED_SNAME_BLANK = _medenum.MED_SNAME_BLANK
MED_NAME_BLANK = _medenum.MED_NAME_BLANK
MED_PATHNAME_SIZE = _medenum.MED_PATHNAME_SIZE
MED_MAX_CHFID_PATH = _medenum.MED_MAX_CHFID_PATH
MED_FULL_INTERLACE = _medenum.MED_FULL_INTERLACE
MED_NO_INTERLACE = _medenum.MED_NO_INTERLACE
MED_UNDEF_INTERLACE = _medenum.MED_UNDEF_INTERLACE
MED_UNDEF_STMODE = _medenum.MED_UNDEF_STMODE
MED_GLOBAL_STMODE = _medenum.MED_GLOBAL_STMODE
MED_COMPACT_STMODE = _medenum.MED_COMPACT_STMODE
MED_GLOBAL_PFLMODE = _medenum.MED_GLOBAL_PFLMODE
MED_COMPACT_PFLMODE = _medenum.MED_COMPACT_PFLMODE
MED_UNDEF_PFLMODE = _medenum.MED_UNDEF_PFLMODE
MED_ACC_RDONLY = _medenum.MED_ACC_RDONLY
MED_ACC_RDWR = _medenum.MED_ACC_RDWR
MED_ACC_RDEXT = _medenum.MED_ACC_RDEXT
MED_ACC_CREAT = _medenum.MED_ACC_CREAT
MED_ACC_UNDEF = _medenum.MED_ACC_UNDEF
MED_UNSTRUCTURED_MESH = _medenum.MED_UNSTRUCTURED_MESH
MED_STRUCTURED_MESH = _medenum.MED_STRUCTURED_MESH
MED_UNDEF_MESH_TYPE = _medenum.MED_UNDEF_MESH_TYPE
MED_CARTESIAN_GRID = _medenum.MED_CARTESIAN_GRID
MED_POLAR_GRID = _medenum.MED_POLAR_GRID
MED_CURVILINEAR_GRID = _medenum.MED_CURVILINEAR_GRID
MED_UNDEF_GRID_TYPE = _medenum.MED_UNDEF_GRID_TYPE
MED_CELL = _medenum.MED_CELL
MED_DESCENDING_FACE = _medenum.MED_DESCENDING_FACE
MED_DESCENDING_EDGE = _medenum.MED_DESCENDING_EDGE
MED_NODE = _medenum.MED_NODE
MED_NODE_ELEMENT = _medenum.MED_NODE_ELEMENT
MED_STRUCT_ELEMENT = _medenum.MED_STRUCT_ELEMENT
MED_ALL_ENTITY_TYPE = _medenum.MED_ALL_ENTITY_TYPE
MED_UNDEF_ENTITY_TYPE = _medenum.MED_UNDEF_ENTITY_TYPE
MED_N_ENTITY_TYPES = _medenum.MED_N_ENTITY_TYPES
MED_COORDINATE = _medenum.MED_COORDINATE
MED_CONNECTIVITY = _medenum.MED_CONNECTIVITY
MED_NAME = _medenum.MED_NAME
MED_NUMBER = _medenum.MED_NUMBER
MED_FAMILY_NUMBER = _medenum.MED_FAMILY_NUMBER
MED_COORDINATE_AXIS1 = _medenum.MED_COORDINATE_AXIS1
MED_COORDINATE_AXIS2 = _medenum.MED_COORDINATE_AXIS2
MED_COORDINATE_AXIS3 = _medenum.MED_COORDINATE_AXIS3
MED_INDEX_FACE = _medenum.MED_INDEX_FACE
MED_INDEX_NODE = _medenum.MED_INDEX_NODE
MED_GLOBAL_NUMBER = _medenum.MED_GLOBAL_NUMBER
MED_VARIABLE_ATTRIBUTE = _medenum.MED_VARIABLE_ATTRIBUTE
MED_COORDINATE_TRSF = _medenum.MED_COORDINATE_TRSF
MED_UNDEF_DATATYPE = _medenum.MED_UNDEF_DATATYPE
MED_INTERNAL_FLOAT64 = _medenum.MED_INTERNAL_FLOAT64
MED_INTERNAL_INT32 = _medenum.MED_INTERNAL_INT32
MED_INTERNAL_INT64 = _medenum.MED_INTERNAL_INT64
MED_INTERNAL_INT = _medenum.MED_INTERNAL_INT
MED_INTERNAL_NAME = _medenum.MED_INTERNAL_NAME
MED_INTERNAL_SNAME = _medenum.MED_INTERNAL_SNAME
MED_INTERNAL_LNAME = _medenum.MED_INTERNAL_LNAME
MED_INTERNAL_IDENT = _medenum.MED_INTERNAL_IDENT
MED_INTERNAL_CHAR = _medenum.MED_INTERNAL_CHAR
MED_INTERNAL_UNDEF = _medenum.MED_INTERNAL_UNDEF
MED_FLOAT64 = _medenum.MED_FLOAT64
MED_INT32 = _medenum.MED_INT32
MED_INT64 = _medenum.MED_INT64
MED_INT = _medenum.MED_INT
MED_ATT_FLOAT64 = _medenum.MED_ATT_FLOAT64
MED_ATT_INT = _medenum.MED_ATT_INT
MED_ATT_NAME = _medenum.MED_ATT_NAME
MED_ATT_UNDEF = _medenum.MED_ATT_UNDEF
MED_MESH = _medenum.MED_MESH
MED_FIELD = _medenum.MED_FIELD
MED_LIBRARY = _medenum.MED_LIBRARY
MED_FILE = _medenum.MED_FILE
MED_MESH_SUPPORT = _medenum.MED_MESH_SUPPORT
MED_ELSTRUCT = _medenum.MED_ELSTRUCT
MED_FAMILY = _medenum.MED_FAMILY
MED_EQUIVALENCE = _medenum.MED_EQUIVALENCE
MED_GROUP = _medenum.MED_GROUP
MED_JOINT = _medenum.MED_JOINT
MED_LOCALIZATION = _medenum.MED_LOCALIZATION
MED_PROFILE = _medenum.MED_PROFILE
MED_FILTER = _medenum.MED_FILTER
MED_INTERPOLATION = _medenum.MED_INTERPOLATION
MED_NUMERICAL_DATA = _medenum.MED_NUMERICAL_DATA
MED_LINK = _medenum.MED_LINK
MED_CLASS_UNDEF = _medenum.MED_CLASS_UNDEF
MED_CLASS_ALL = _medenum.MED_CLASS_ALL
MED_POINT1 = _medenum.MED_POINT1
MED_SEG2 = _medenum.MED_SEG2
MED_SEG3 = _medenum.MED_SEG3
MED_SEG4 = _medenum.MED_SEG4
MED_TRIA3 = _medenum.MED_TRIA3
MED_QUAD4 = _medenum.MED_QUAD4
MED_TRIA6 = _medenum.MED_TRIA6
MED_TRIA7 = _medenum.MED_TRIA7
MED_QUAD8 = _medenum.MED_QUAD8
MED_QUAD9 = _medenum.MED_QUAD9
MED_TETRA4 = _medenum.MED_TETRA4
MED_PYRA5 = _medenum.MED_PYRA5
MED_PENTA6 = _medenum.MED_PENTA6
MED_HEXA8 = _medenum.MED_HEXA8
MED_TETRA10 = _medenum.MED_TETRA10
MED_OCTA12 = _medenum.MED_OCTA12
MED_PYRA13 = _medenum.MED_PYRA13
MED_PENTA15 = _medenum.MED_PENTA15
MED_HEXA20 = _medenum.MED_HEXA20
MED_HEXA27 = _medenum.MED_HEXA27
MED_POLYGON = _medenum.MED_POLYGON
MED_POLYGON2 = _medenum.MED_POLYGON2
MED_POLYHEDRON = _medenum.MED_POLYHEDRON
MED_STRUCT_GEO_INTERNAL = _medenum.MED_STRUCT_GEO_INTERNAL
MED_STRUCT_GEO_SUP_INTERNAL = _medenum.MED_STRUCT_GEO_SUP_INTERNAL
MED_NONE = _medenum.MED_NONE
MED_NO_GEOTYPE = _medenum.MED_NO_GEOTYPE
MED_UNDEF_GEOTYPE = _medenum.MED_UNDEF_GEOTYPE
MED_UNDEF_GEOMETRY_TYPE = _medenum.MED_UNDEF_GEOMETRY_TYPE
MED_ALL_GEOTYPE = _medenum.MED_ALL_GEOTYPE
MED_GEO_ALL = _medenum.MED_GEO_ALL
MED_N_CELL_GEO = _medenum.MED_N_CELL_GEO
MED_N_CELL_FIXED_GEO = _medenum.MED_N_CELL_FIXED_GEO
MED_N_CELL_GEO_FIXED_CON = _medenum.MED_N_CELL_GEO_FIXED_CON
MED_N_FACE_GEO = _medenum.MED_N_FACE_GEO
MED_N_FACE_FIXED_GEO = _medenum.MED_N_FACE_FIXED_GEO
MED_N_FACE_GEO_FIXED_CON = _medenum.MED_N_FACE_GEO_FIXED_CON
MED_N_EDGE_TYPES = _medenum.MED_N_EDGE_TYPES
MED_N_EDGE_FIXED_GEO = _medenum.MED_N_EDGE_FIXED_GEO
MED_N_EDGE_GEO_FIXED_CON = _medenum.MED_N_EDGE_GEO_FIXED_CON
MED_N_NODE_GEO = _medenum.MED_N_NODE_GEO
MED_N_NODE_FIXED_GEO = _medenum.MED_N_NODE_FIXED_GEO
MED_N_NODE_GEO_FIXED_CON = _medenum.MED_N_NODE_GEO_FIXED_CON
MED_NODAL = _medenum.MED_NODAL
MED_DESCENDING = _medenum.MED_DESCENDING
MED_UNDEF_CONNECTIVITY_MODE = _medenum.MED_UNDEF_CONNECTIVITY_MODE
MED_NO_CMODE = _medenum.MED_NO_CMODE
MED_CARTESIAN = _medenum.MED_CARTESIAN
MED_CYLINDRICAL = _medenum.MED_CYLINDRICAL
MED_SPHERICAL = _medenum.MED_SPHERICAL
MED_UNDEF_AXIS_TYPE = _medenum.MED_UNDEF_AXIS_TYPE
MED_FALSE = _medenum.MED_FALSE
MED_TRUE = _medenum.MED_TRUE
MED_GAUSS_ELNO = _medenum.MED_GAUSS_ELNO
MED_IPOINT_ELNO = _medenum.MED_IPOINT_ELNO
MED_NO_NAME = _medenum.MED_NO_NAME
MED_NO_MESHNAME = _medenum.MED_NO_MESHNAME
MED_NO_MESH = _medenum.MED_NO_MESH
MED_NO_MESH_SUPPORT = _medenum.MED_NO_MESH_SUPPORT
MED_NO_LOCALIZATION = _medenum.MED_NO_LOCALIZATION
MED_NO_INTERPOLATION = _medenum.MED_NO_INTERPOLATION
MED_NO_IPOINT_INTERNAL = _medenum.MED_NO_IPOINT_INTERNAL
MED_NO_PROFILE = _medenum.MED_NO_PROFILE
MED_NO_GROUP = _medenum.MED_NO_GROUP
MED_ALLENTITIES_PROFILE = _medenum.MED_ALLENTITIES_PROFILE
MED_NO_PROFILE_INTERNAL = _medenum.MED_NO_PROFILE_INTERNAL
MED_SAME_PROFILE_INTERNAL = _medenum.MED_SAME_PROFILE_INTERNAL
MED_ALL_CONSTITUENT = _medenum.MED_ALL_CONSTITUENT
MED_UNDEF_SIZE = _medenum.MED_UNDEF_SIZE
MED_NO_PROFILE_SIZE = _medenum.MED_NO_PROFILE_SIZE
MED_SORT_DTIT = _medenum.MED_SORT_DTIT
MED_SORT_ITDT = _medenum.MED_SORT_ITDT
MED_SORT_UNDEF = _medenum.MED_SORT_UNDEF
MED_NO_DT = _medenum.MED_NO_DT
MED_NO_IT = _medenum.MED_NO_IT
MED_UNDEF_DT = _medenum.MED_UNDEF_DT
MED_ATT_NOT_FILLED = _medenum.MED_ATT_NOT_FILLED
MED_MAX_FILTER_SPACES = _medenum.MED_MAX_FILTER_SPACES
class med_filter(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, med_filter, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, med_filter, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nspaces"] = _medenum.med_filter_nspaces_set
    __swig_getmethods__["nspaces"] = _medenum.med_filter_nspaces_get
    if _newclass:nspaces = _swig_property(_medenum.med_filter_nspaces_get, _medenum.med_filter_nspaces_set)
    __swig_setmethods__["memspace"] = _medenum.med_filter_memspace_set
    __swig_getmethods__["memspace"] = _medenum.med_filter_memspace_get
    if _newclass:memspace = _swig_property(_medenum.med_filter_memspace_get, _medenum.med_filter_memspace_set)
    __swig_setmethods__["diskspace"] = _medenum.med_filter_diskspace_set
    __swig_getmethods__["diskspace"] = _medenum.med_filter_diskspace_get
    if _newclass:diskspace = _swig_property(_medenum.med_filter_diskspace_get, _medenum.med_filter_diskspace_set)
    __swig_setmethods__["nentity"] = _medenum.med_filter_nentity_set
    __swig_getmethods__["nentity"] = _medenum.med_filter_nentity_get
    if _newclass:nentity = _swig_property(_medenum.med_filter_nentity_get, _medenum.med_filter_nentity_set)
    __swig_setmethods__["nvaluesperentity"] = _medenum.med_filter_nvaluesperentity_set
    __swig_getmethods__["nvaluesperentity"] = _medenum.med_filter_nvaluesperentity_get
    if _newclass:nvaluesperentity = _swig_property(_medenum.med_filter_nvaluesperentity_get, _medenum.med_filter_nvaluesperentity_set)
    __swig_setmethods__["nconstituentpervalue"] = _medenum.med_filter_nconstituentpervalue_set
    __swig_getmethods__["nconstituentpervalue"] = _medenum.med_filter_nconstituentpervalue_get
    if _newclass:nconstituentpervalue = _swig_property(_medenum.med_filter_nconstituentpervalue_get, _medenum.med_filter_nconstituentpervalue_set)
    __swig_setmethods__["constituentselect"] = _medenum.med_filter_constituentselect_set
    __swig_getmethods__["constituentselect"] = _medenum.med_filter_constituentselect_get
    if _newclass:constituentselect = _swig_property(_medenum.med_filter_constituentselect_get, _medenum.med_filter_constituentselect_set)
    __swig_setmethods__["switchmode"] = _medenum.med_filter_switchmode_set
    __swig_getmethods__["switchmode"] = _medenum.med_filter_switchmode_get
    if _newclass:switchmode = _swig_property(_medenum.med_filter_switchmode_get, _medenum.med_filter_switchmode_set)
    __swig_setmethods__["filterarraysize"] = _medenum.med_filter_filterarraysize_set
    __swig_getmethods__["filterarraysize"] = _medenum.med_filter_filterarraysize_get
    if _newclass:filterarraysize = _swig_property(_medenum.med_filter_filterarraysize_get, _medenum.med_filter_filterarraysize_set)
    __swig_setmethods__["filterarray23v30"] = _medenum.med_filter_filterarray23v30_set
    __swig_getmethods__["filterarray23v30"] = _medenum.med_filter_filterarray23v30_get
    if _newclass:filterarray23v30 = _swig_property(_medenum.med_filter_filterarray23v30_get, _medenum.med_filter_filterarray23v30_set)
    __swig_setmethods__["profilearraysize"] = _medenum.med_filter_profilearraysize_set
    __swig_getmethods__["profilearraysize"] = _medenum.med_filter_profilearraysize_get
    if _newclass:profilearraysize = _swig_property(_medenum.med_filter_profilearraysize_get, _medenum.med_filter_profilearraysize_set)
    __swig_setmethods__["storagemode"] = _medenum.med_filter_storagemode_set
    __swig_getmethods__["storagemode"] = _medenum.med_filter_storagemode_get
    if _newclass:storagemode = _swig_property(_medenum.med_filter_storagemode_get, _medenum.med_filter_storagemode_set)
    __swig_setmethods__["profilename"] = _medenum.med_filter_profilename_set
    __swig_getmethods__["profilename"] = _medenum.med_filter_profilename_get
    if _newclass:profilename = _swig_property(_medenum.med_filter_profilename_get, _medenum.med_filter_profilename_set)
    def __init__(self): 
        this = _medenum.new_med_filter()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _medenum.delete_med_filter
    __del__ = lambda self : None;
med_filter_swigregister = _medenum.med_filter_swigregister
med_filter_swigregister(med_filter)

MED_NO_FILTER_SIZE = _medenum.MED_NO_FILTER_SIZE
MED_NO_PROFILE_F = _medenum.MED_NO_PROFILE_F
class med_file_version(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, med_file_version, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, med_file_version, name)
    __repr__ = _swig_repr
    __swig_setmethods__["majeur"] = _medenum.med_file_version_majeur_set
    __swig_getmethods__["majeur"] = _medenum.med_file_version_majeur_get
    if _newclass:majeur = _swig_property(_medenum.med_file_version_majeur_get, _medenum.med_file_version_majeur_set)
    __swig_setmethods__["mineur"] = _medenum.med_file_version_mineur_set
    __swig_getmethods__["mineur"] = _medenum.med_file_version_mineur_get
    if _newclass:mineur = _swig_property(_medenum.med_file_version_mineur_get, _medenum.med_file_version_mineur_set)
    __swig_setmethods__["release"] = _medenum.med_file_version_release_set
    __swig_getmethods__["release"] = _medenum.med_file_version_release_get
    if _newclass:release = _swig_property(_medenum.med_file_version_release_get, _medenum.med_file_version_release_set)
    def __init__(self): 
        this = _medenum.new_med_file_version()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _medenum.delete_med_file_version
    __del__ = lambda self : None;
med_file_version_swigregister = _medenum.med_file_version_swigregister
med_file_version_swigregister(med_file_version)

MED_PARTICLE_NAME = _medenum.MED_PARTICLE_NAME
MED_BALL_NAME = _medenum.MED_BALL_NAME
MED_BEAM_NAME = _medenum.MED_BEAM_NAME
MED_PARTICLE_LABEL = _medenum.MED_PARTICLE_LABEL
MED_BALL_DIAMETER = _medenum.MED_BALL_DIAMETER
MED_BEAM_THICKNESS = _medenum.MED_BEAM_THICKNESS
class MED_MESH_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_MESH_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_MESH_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_MESH_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_MESH_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_MESH_TYPE_val_get, _medenum.MED_MESH_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_MESH_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_MESH_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_MESH_TYPE
    __del__ = lambda self : None;
MED_MESH_TYPE_swigregister = _medenum.MED_MESH_TYPE_swigregister
MED_MESH_TYPE_swigregister(MED_MESH_TYPE)
cvar = _medenum.cvar
MED_GET_ENTITY_TYPENAME = cvar.MED_GET_ENTITY_TYPENAME
MED_GET_CELL_GEOMETRY_TYPENAME = cvar.MED_GET_CELL_GEOMETRY_TYPENAME
MED_GET_FACE_GEOMETRY_TYPENAME = cvar.MED_GET_FACE_GEOMETRY_TYPENAME

MED_MESH_TYPE.__repr__= lambda self:"MED_MESH_TYPE" +"("+str(self.val)+")"

class MED_SORTING_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_SORTING_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_SORTING_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_SORTING_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_SORTING_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_SORTING_TYPE_val_get, _medenum.MED_SORTING_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_SORTING_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_SORTING_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_SORTING_TYPE
    __del__ = lambda self : None;
MED_SORTING_TYPE_swigregister = _medenum.MED_SORTING_TYPE_swigregister
MED_SORTING_TYPE_swigregister(MED_SORTING_TYPE)

MED_SORTING_TYPE.__repr__= lambda self:"MED_SORTING_TYPE" +"("+str(self.val)+")"

class MED_GRID_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_GRID_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_GRID_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_GRID_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_GRID_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_GRID_TYPE_val_get, _medenum.MED_GRID_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_GRID_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_GRID_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_GRID_TYPE
    __del__ = lambda self : None;
MED_GRID_TYPE_swigregister = _medenum.MED_GRID_TYPE_swigregister
MED_GRID_TYPE_swigregister(MED_GRID_TYPE)

MED_GRID_TYPE.__repr__= lambda self:"MED_GRID_TYPE" +"("+str(self.val)+")"

class MED_ENTITY_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_ENTITY_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_ENTITY_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_ENTITY_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_ENTITY_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_ENTITY_TYPE_val_get, _medenum.MED_ENTITY_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_ENTITY_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_ENTITY_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_ENTITY_TYPE
    __del__ = lambda self : None;
MED_ENTITY_TYPE_swigregister = _medenum.MED_ENTITY_TYPE_swigregister
MED_ENTITY_TYPE_swigregister(MED_ENTITY_TYPE)

MED_ENTITY_TYPE.__repr__= lambda self:"MED_ENTITY_TYPE" +"("+str(self.val)+")"

class MED_FIELD_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_FIELD_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_FIELD_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_FIELD_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_FIELD_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_FIELD_TYPE_val_get, _medenum.MED_FIELD_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_FIELD_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_FIELD_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_FIELD_TYPE
    __del__ = lambda self : None;
MED_FIELD_TYPE_swigregister = _medenum.MED_FIELD_TYPE_swigregister
MED_FIELD_TYPE_swigregister(MED_FIELD_TYPE)

MED_FIELD_TYPE.__repr__= lambda self:"MED_FIELD_TYPE" +"("+str(self.val)+")"

class MED_ATTRIBUTE_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_ATTRIBUTE_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_ATTRIBUTE_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_ATTRIBUTE_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_ATTRIBUTE_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_ATTRIBUTE_TYPE_val_get, _medenum.MED_ATTRIBUTE_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_ATTRIBUTE_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_ATTRIBUTE_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_ATTRIBUTE_TYPE
    __del__ = lambda self : None;
MED_ATTRIBUTE_TYPE_swigregister = _medenum.MED_ATTRIBUTE_TYPE_swigregister
MED_ATTRIBUTE_TYPE_swigregister(MED_ATTRIBUTE_TYPE)

MED_ATTRIBUTE_TYPE.__repr__= lambda self:"MED_ATTRIBUTE_TYPE" +"("+str(self.val)+")"

class MED_AXIS_TYPE(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MED_AXIS_TYPE, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MED_AXIS_TYPE, name)
    __repr__ = _swig_repr
    __swig_setmethods__["val"] = _medenum.MED_AXIS_TYPE_val_set
    __swig_getmethods__["val"] = _medenum.MED_AXIS_TYPE_val_get
    if _newclass:val = _swig_property(_medenum.MED_AXIS_TYPE_val_get, _medenum.MED_AXIS_TYPE_val_set)
    def __init__(self, *args): 
        this = _medenum.new_MED_AXIS_TYPE(*args)
        try: self.this.append(this)
        except: self.this = this
    def __str__(self): return _medenum.MED_AXIS_TYPE___str__(self)
    __swig_destroy__ = _medenum.delete_MED_AXIS_TYPE
    __del__ = lambda self : None;
MED_AXIS_TYPE_swigregister = _medenum.MED_AXIS_TYPE_swigregister
MED_AXIS_TYPE_swigregister(MED_AXIS_TYPE)

MED_AXIS_TYPE.__repr__= lambda self:"MED_AXIS_TYPE" +"("+str(self.val)+")"

MED_GET_ENTITY_TYPE=[
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE"),
    (MED_NODE,		   "MED_NODE"),
    (MED_CELL,		   "MED_CELL"),
    (MED_DESCENDING_FACE,  "MED_DESCENDING_FACE"),
    (MED_DESCENDING_EDGE,  "MED_DESCENDING_EDGE"),
    (MED_NODE_ELEMENT,	   "MED_NODE_ELEMENT"),
    (MED_STRUCT_ELEMENT,   "MED_STRUCT_ELEMENT")
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE")
]

MED_GET_NODAL_ENTITY_TYPE=[
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE"),
    (MED_NODE,		   "MED_NODE"),
    (MED_CELL,		   "MED_CELL"),
    (MED_NODE_ELEMENT,	   "MED_NODE_ELEMENT"),
    (MED_STRUCT_ELEMENT,   "MED_STRUCT_ELEMENT")
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE")
]

MED_GET_DESCENDING_ENTITY_TYPE=[
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE"),
    (MED_NODE,		   "MED_NODE"),
    (MED_CELL,		   "MED_CELL"),
    (MED_DESCENDING_FACE,  "MED_DESCENDING_FACE"),
    (MED_DESCENDING_EDGE,  "MED_DESCENDING_EDGE"),
    (MED_NODE_ELEMENT,	   "MED_NODE_ELEMENT"),
    (MED_STRUCT_ELEMENT,   "MED_STRUCT_ELEMENT")
#    (MED_UNDEF_ENTITY_TYPE,"MED_UNDEF_ENTITY_TYPE")
]

#Ordre iteration important pour la numerotation globale induite
MED_GET_CELL_GEOMETRY_TYPE=[
# (MED_NO_GEOTYPE, "MED_NO_GEOTYPE"),
 (MED_POINT1,     "MED_POINT1"),
 (MED_SEG2,       "MED_SEG2"),
 (MED_SEG3,       "MED_SEG3"),
 (MED_SEG4,       "MED_SEG4"),
 (MED_TRIA3,      "MED_TRIA3"),
 (MED_QUAD4,      "MED_QUAD4"),
 (MED_TRIA6,      "MED_TRIA6"),
 (MED_TRIA7,      "MED_TRIA7"),
 (MED_QUAD8,      "MED_QUAD8"),
 (MED_QUAD9,      "MED_QUAD9"),
 (MED_TETRA4,     "MED_TETRA4"),
 (MED_PYRA5,      "MED_PYRA5"),
 (MED_PENTA6,     "MED_PENTA6"),
 (MED_HEXA8,      "MED_HEXA8"),
 (MED_TETRA10,    "MED_TETRA10"),
 (MED_OCTA12,     "MED_OCTA12"),
 (MED_PYRA13,     "MED_PYRA13"),
 (MED_PENTA15,    "MED_PENTA15"),
 (MED_HEXA20,     "MED_HEXA20"),
 (MED_HEXA27,     "MED_HEXA27"),
 (MED_POLYGON,    "MED_POLYGON"),
 (MED_POLYGON2,    "MED_POLYGON2"),
 (MED_POLYHEDRON, "MED_POLYHEDRON")
# (MED_NO_GEOTYPE, "MED_NO_GEOTYPE")
]

#Ordre iteration important pour la numerotation globale induite
MED_GET_FACE_GEOMETRY_TYPE=[
#  (MED_NO_GEOTYPE, "MED_NO_GEOTYPE"),
  (MED_TRIA3,	   "MED_TRIA3"),
  (MED_QUAD4,	   "MED_QUAD4"),
  (MED_TRIA6,	   "MED_TRIA6"),
  (MED_TRIA7,	   "MED_TRIA7"),
  (MED_QUAD8,	   "MED_QUAD8"),
  (MED_QUAD9,	   "MED_QUAD9"),
  (MED_POLYGON,	   "MED_POLYGON")
#  (MED_NO_GEOTYPE, "MED_NO_GEOTYPE")
]

#Ordre iteration important pour la numerotation globale induite
MED_GET_EDGE_GEOMETRY_TYPE=[
#  (MED_NO_GEOTYPE,"MED_NO_GEOTYPE"),
  (MED_SEG2,	 "MED_SEG2"),
  (MED_SEG3,	 "MED_SEG3"),
  (MED_SEG4,	 "MED_SEG4")
#  (MED_NO_GEOTYPE,"MED_NO_GEOTYPE")
]

#Ordre iteration important pour la numerotation globale induite
MED_GET_NODE_GEOMETRY_TYPE=[
# (MED_NO_GEOTYPE,"MED_NO_GEOTYPE"),
 (MED_NO_GEOTYPE,"MED_NO_GEOTYPE")
# (MED_NO_GEOTYPE,"MED_NO_GEOTYPE")
]


MED_GET_GEOTYPE=dict([
    (MED_NODE,		   MED_GET_NODE_GEOMETRY_TYPE),
    (MED_DESCENDING_FACE,  MED_GET_FACE_GEOMETRY_TYPE),
    (MED_DESCENDING_EDGE,  MED_GET_EDGE_GEOMETRY_TYPE),
    (MED_CELL,		   MED_GET_CELL_GEOMETRY_TYPE),
    (MED_NODE_ELEMENT,	   MED_GET_CELL_GEOMETRY_TYPE),
    (MED_STRUCT_ELEMENT,   [(MED_NO_GEOTYPE,"MED_NO_GEOTYPE")])
])


# This file is compatible with both classic and new-style classes.


