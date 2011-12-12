# -*- mode: python; -*- 

import os
import os.path as op
import sys
import commands

import Configure 
import Options
from waflib import Context
from waflib.TaskGen import feature, before
from waflib.Tools import c_preproc
c_preproc.go_absolute = True

import frogsutils


def get_out_name():
    """
    Generate an (architecture dependant) 
    build directory name.
    """
    ret = os.uname()
    out_name = 'build.'+ ret[0]+'-'+ret[-1]
    return op.join('build', out_name)


def options(opt):
    """
    Default options for the frogs' packages.
    """
    opt.load('compiler_cc')
    opt.load('compiler_cxx')
    opt.load('compiler_fc')
    
    opt.add_option('--with-cfitsio', 
                   action='store', 
                   help='Path to the cfitsio root dir')
    opt.add_option('--with-sex', 
                   action='store', 
                   help='Path to the cfitsio root dir')    
    opt.add_option('--no-global-lapack', 
                   action='store_false',
                   default=True, 
                   dest='global_lapack',
                   help='use/no do use system lapack library')    
    opt.add_option("--no-cernlib", 
                   action='store_false', 
                   default=True, 
                   dest='cernlib', 
                   help='do not try to link with cernlib')
    # opt.add_option('--debug', 
    #                action='store', 
    #                default=1, 
    #                dest='debug', 
    #                help='turn on the -g option')
    # opt.add_option('--optimize', 
    #                action='store', 
    #                default=1, 
    #                dest='optimize', 
    #                help='turn on the -O3 option')    
    opt.add_option('--enable-lyon', 
                   action='store_true', 
                   default=False,
                   dest='enable_lyon',
                   help='We are at CCIN2P3. Things behave strangely there')


@Configure.conftest 
def configure(conf):
    """
    Loads all the tools usually needed to compile poloka and co 
    """
    
    print "*** WAF: loading FROGS specific stuff ***"
    
    # c compiler 
    conf.load( 'compiler_c' )
    conf.env['CCFLAGS'] = ['-fPIC', '-DPIC']
    
    # c++ compiler 
    conf.load( 'compiler_cxx' )
    conf.env['CXXFLAGS'] = ['-fPIC', '-DPIC']
    
    # fortran compilers 
    # apparently, we need to check for a c-compiler before.
    conf.load( 'compiler_fc')
    conf.check_fortran()
    conf.check_fortran_verbose_flag()
    conf.check_fortran_clib()
    conf.check_fortran_dummy_main()
    conf.check_fortran_mangling()
    
    # Debug option 
    debug_flag = False
    try:
        debug_flag = Context.g_module.debug
    except: pass
    if debug_flag:
        conf.env['CFLAGS'].append('-g')
        conf.env['CXXFLAGS'].append('-g')
    
    # Optimizer ? 
    opt_level = 0
    try:
        opt_level = Context.g_module.optimize
    except: pass
    if opt_level > 0:        
        conf.env['CFLAGS'].append('-O%d' % opt_level)
        conf.env['CXXFLAGS'].append('-O%d' % opt_level)
        conf.env['FCFLAGS'] = '-O%d' % opt_level
    

    # lapack
    if conf.options.global_lapack:
        if conf.check_cc( lib='lapack', msg='checking for lapack' ):
            conf.env.global_lapack = True
        else:
            conf.env.global_lapack = False

    # poloka-include
    pkg = Context.g_module.APPNAME
    ver = Context.g_module.VERSION
    conf.env['PKG_INCDIR'] = op.join('include', '%s-%s' % (pkg,ver))

    # poloka packages 
    pkgdir = os.getenv('POLOKA_PACKAGES')
    if pkgdir:
        conf.env['POLOKA_PACKAGES'] = pkgdir
    

def load_pkg_config(conf):
    """
    Load pkg-config (used to read the package metadata)
    """
    
    # first, we need pkg-config 
    if not conf.env['PKG_CONFIG'] or \
       not conf.env['PKG_CONFIG_PATH']:
        try:
            conf.find_program('pkg-config', var='PKG_CONFIG')
            pkgcpath = op.join( conf.env.PREFIX, 'lib', 'pkgconfig')
            if os.environ.has_key('PKG_CONFIG_PATH'):
                os.environ['PKG_CONFIG_PATH'] = pkgcpath + ":" + os.environ['PKG_CONFIG_PATH']
            else:
                os.environ['PKG_CONFIG_PATH'] = pkgcpath
        except Configure.ConfigurationError:
            conf.fatal('pkg-config not found')
    

@Configure.conftest
def check_cernlib(conf, mandatory=True):
    """
    Attempt to see whether the cernlib is around. 
    
    Using the cernlib is a bit tricky. As of today, the shared
    versions of cernlib are known not to work on 64-bit linux. So, we
    attempt to locate the static libs.
    """

    if not conf.options.cernlib:
        return
    
    def parse_new_cernlib(conf, l):
        
        s = ""
        static_libs = []
        static_libpaths = []
        dy_libs = []
        dy_libpaths = []
        
        current_libs  = None
        current_libpaths = None
        afs_sys_name = conf.env['AFS_SYS_NAME']
        if afs_sys_name == []: afs_sys_name = ""
        
        while True:
            try: 
                s = l.pop(0)
            except: 
                break
            
            if s.startswith('-Wl,-static'):
                current_libs = static_libs
                current_libpaths = static_libpaths
            elif s.startswith('-Wl,-dy'):
                current_libs = dy_libs
                current_libpaths = dy_libpaths
            else:
                if current_libs == None or \
                   current_libpaths == None:
                    conf.fatal('unable to parse cernlib output')
                elif s.startswith('-L'):
                    current_libpaths.append(s[2:].replace('@sys', afs_sys_name))
                elif s.startswith('-l'):
                    current_libs.append(s[2:])

        for stlib in static_libs:
            conf.check_cc(stlib=stlib, uselib_store='CERN')
        for dylib in dy_libs:
            conf.check_cc(lib=dylib, uselib_store='CERN')
            
                
    def parse_traditional_cernlib(conf, l):
        """
        """

        static_libs = []
        static_libpaths = []
        dy_libs = []
        dy_libpaths = []
        
        afs_sys_name = frogsutils.get_afs_sys_name()
        
        while True:
            try: s = l.pop(0)
            except: break
            if s.startswith('/'):
                # we should check that these things exist ...
                libname = op.basename(s)
                if libname.startswith('lib'):
                    libname = libname[3:-2]
                static_libs.append(libname)
                static_libpaths.append(op.dirname(s).replace('@sys', afs_sys_name))
            elif s.startswith('-l'):
                dy_libs.append(s[2:])
            elif s.startswith('-L'):
                dy_libpaths.append(s[2:].replace('@sys', afs_sys_name))

                
        conf.env.STLIB_CERN = static_libs
        conf.env.STLIBPATH = static_libpaths
        conf.env.LIB_CERN       = dy_libs
        conf.env.LIBPATH_CERN   = dy_libpaths

                
                    
    conf.find_program('cernlib', mandatory=mandatory)
    ret = commands.getstatusoutput('cernlib packlib kernlib')
    if ret[0] != 0:
        if mandatory:
            conf.fatal('unable to locate cernlib')
        else:
            conf.env['HAVE_CERN'] = 0
            return
        
    # if we got something from `cernlib`,
    # we attempt to extract the paths from it.
    l = ret[1].split()
    if l[0].startswith('-Wl,-static'):
        parse_new_cernlib(conf, l)
    else:
        parse_traditional_cernlib(conf, l)
    
    conf.env['HAVE_CERN'] = 1


@Configure.conftest
def check_packages(conf, pkg_list):
    """
    - Check whether the packages passed in argument 
      can be found by pkg-config. 
    - Parse the cflags and libs and store the information 
      in the uselib_store.
    """
    load_pkg_config(conf)
    # check for the packages passed in arguments 
    for pkg in pkg_list:
        try:
            pkg_name, pkg_version, pkg_mandatory = pkg
            conf.check_cfg(args='--cflags --libs', 
                           package = pkg_name + '-' + pkg_version, 
                           mandatory=pkg_mandatory, 
                           uselib_store=pkg_name.replace('-', '_').upper())
        except Configure.ConfigurationError:
            conf.fatal('unable to locate %s-%s (mandatory)' % (pkg_name, pkg_version))


@feature('cshlib')
@feature('cxxshlib')
@feature('cprogram')
@feature('cxxprogram')
@before('process_subst')
def hack_shlib_rpath(self):
    bld = self.bld 
    requirements = Context.g_module.requirements
    rpath = [self.bld.env.PREFIX + os.sep + 'lib']
    for pkg in requirements:
        pkg_name = pkg[0].replace('-', '_').upper()
        rp = self.bld.env['LIBPATH_' + pkg_name]
        if rp.__class__ is list:
            rpath.extend(rp)
        else:
            rpath.append(rp)
    rpath = set(rpath)
    for g in self.bld.groups:
        for tg in g:
            try:
                if 'cshlib' in tg.features or \
                   'cxxshlib' in tg.features or \
                   'cprogram' in tg.features or \
                   'cxxprogram' in tg.features:
                    ret = getattr(tg, 'rpath', [])
                    ret.extend(rpath)
                    tg.rpath = ret
            except:
                pass
    
# @feature('cprogram')
# @feature('cxxprogram')
# @before('process_subst')
# def hack_program_rpath(self):
#     rpath = self.bld.env.PREFIX + os.sep + 'lib'
#     for g in self.bld.groups:
#         for tg in g:
#             try:
#                 if 'cprogram' in tg.features or \
#                    'cxxprogram' in tg.features:
#                     ret = getattr(tg, 'rpath', [])
#                     if rpath in ret: continue
#                     ret.append(rpath)
#                     tg.rpath = ret
#                     print "hack_program_rpath !", tg.name, rpath
#             except:
#                 pass
    
def install_headers(bld, headers):
    """
    Just hide the install header commands 
    """
    name    = Context.g_module.APPNAME
    version = Context.g_module.VERSION
    install_dir = op.join('$PREFIX', 'include', '%s-%s' % (name, version))
    
    bld.install_files(install_dir, headers)


def gen_pkgconfig(bld):
    # solution directly from T. Nagy (see email in google-groups)
    from waflib.TaskGen import feature, before, after
    @feature('subst')
    @before('process_subst')
    def read_libs(self):
        self.env.append_value('ALL_LIBS', " ")
        for g in self.bld.groups:
            for tg in g:
                try:
                    if 'cshlib' in tg.features or \
                       'cxxshlib' in tg.features or \
                       'fcshlib' in tg.features:
                        # or 'cxxshlib'/'fcshlib'/'dshlib' in tg.features
                        self.env.append_value('ALL_LIBS', "-l" + tg.name)                    
                except:
                    pass
                
    @feature('subst')
    @after('process_subst')
    def link_pcfile_to_poloka_packages(self):
        for g in self.bld.groups:
            for tg in g:
                try:
                    if 'subst' in tg.features:
                        pcfile = tg.tasks[0].outputs[0]
                        if self.bld.env.POLOKA_PACKAGES:
                            self.bld.symlink_as(self.bld.env.POLOKA_PACKAGES + os.sep + obj.name, 
                                                pcfile.abspath())
                except:
                    pass
        
    appname = Context.g_module.APPNAME
    version = Context.g_module.VERSION    
    description = Context.g_module.description
    reqs = Context.g_module.requirements
    requirements = ""
    for r in reqs:
        requirements += r[0] + "-" + r[1] + " "
    obj = bld(features = 'subst', 
              target = '%s-%s.pc' % (appname, version),
              source = 'pc.in', 
              install_path =  '${PREFIX}/lib/pkgconfig/',
              PREFIX = bld.env['PREFIX'], 
              APPNAME   = appname,
              DESCRIPTION = description,
              VERSION = version,
              REQUIREMENTS = requirements)
    # LIBS = ["-L${libdir} "] + bld.env['ALL_LIBS']

# for the moment, we hook up this function 
# to the main module... May change in the future 
Context.g_module.__dict__['gen_pkgconfig'] = gen_pkgconfig
Context.g_module.__dict__['install_headers'] = install_headers
# Context.g_module.__dict__['get_required_package_libpaths'] =  get_required_package_libpaths

