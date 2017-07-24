VERSION = '0.2.0'
APPNAME = 'redsvd'

srcdir = '.'
blddir = 'build'

def set_options(ctx):
  ctx.tool_options('compiler_cxx')
  ctx.tool_options('unittestt')

def configure(ctx):
  ctx.check_tool('compiler_cxx')
  ctx.check_tool('unittestt')
  ctx.env.CXXFLAGS += ['-O2', '-Wall', '-g']

import Scripting
Scripting.dist_exts += ['.sh']

def dist_hook():
    import os
    os.remove('googlecode_upload.py')

def build(bld):
  bld.recurse('vendor/redsvd-0.2.0/src vendor/redsvd-0.2.0/test vendor/redsvd-0.2.0/sample')
