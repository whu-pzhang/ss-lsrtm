import os
import sys

try:
    import bldutil
    glob_build = True  # scons command launched in RSFSRC
    srcroot = '../..'  # cwd is RSFSRC/build/user/pzhang
    Import('env bindir libdir pkgdir')
    env = env.Clone()
except:
    glob_build = False  # scons command launched in the local directory
    srcroot = os.environ.get('RSFSRC', '../..')
    sys.path.append(os.path.join(srcroot, 'framework'))
    import bldutil
    env = bldutil.Debug()  # Debugging flags for compilers
    # add -pg to enable gprof
    env.Prepend(CCFLAGS=['-pg'])
    env.Prepend(LINKFLAGS=['-pg'])
    ###############################
    bindir = libdir = pkgdir = None
    SConscript(os.path.join(srcroot, 'su/lib/SConstruct'))

targets = bldutil.UserSconsTargets()

# C mains
targets.c = '''
prertm2d_v03
lsprertm2d
lsrtmse
'''


# Python targets
targets.py = '''
setspk
'''

dynlib = env.get('DYNLIB', '')

env.Prepend(LIBS=['rsfpwd', 'rsf'])

targets.build_all(env, glob_build, srcroot, bindir, libdir, pkgdir)
