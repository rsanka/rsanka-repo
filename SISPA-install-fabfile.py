#
# SISPA-install-fabfile.py
#

import os.path, re
from fabric.api import cd, env, hide, local, run, settings, sudo, task
from fabric.network import disconnect_all

@task
def install():
    try:
        _define_dirs()
        _print_env_variables()
        _initialize_dirs()
        _add_tools()
        _add_refs()
    finally:
        disconnect_all()

def _define_dirs():
    env.PROJECT_DIR = "%(ROOT_DIR)s/project" % env
    env.REFS_DIR = "%(ROOT_DIR)s/refs" % env
    env.TOOLS_DIR = "%(ROOT_DIR)s/tools" % env
    env.TOOLS_RUBYLIB_DIR = "%(TOOLS_DIR)s/RUBYLIB" % env
    env.TOOLS_BINARIES_DIR = "%(TOOLS_DIR)s/BINARIES" % env

def _print_env_variables():
    print("user: %(user)s" % env)
    print("host: %(host)s" % env)
    print("ROOT: %s(ROOT_DIR)s" % env)
    print("PROJECT: %(PROJECT_DIR)s" % env)
    print("REFS: %(REFS_DIR)s" % env)
    print("TOOLS: %(TOOLS_DIR)s" % env)
    print("TOOLS RUBYLIB: %(TOOLS_RUBYLIB_DIR)s" % env)
    print("TOOLS BINARIES: %(TOOLS_BINARIES_DIR)s" % env)

def _initialize_dirs():
    sudo("mkdir -p %(PROJECT_DIR)s" % env)
    sudo("mkdir -p %(REFS_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_RUBYLIB_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_BINARIES_DIR)s" % env)

def _add_tools():
    _add_tarball(env.RUBYLIB_URL,env.RUBYLIB_TARBALL,env.TOOLS_RUBYLIB_DIR)
    _add_tarball(env.BINARIES_URL,env.BINARIES_TARBALL,env.TOOLS_BINARIES_DIR)

def _add_refs():
    _add_tarball(env.REFS_URL,env.REFS_TARBALL,env.REFS_DIR)

def _add_tarball(download_url,tarball,install_dir):
    print("du: %s" % download_url)
    print("tf: %s" % tarball)
    print("id: %s" % install_dir)
    sudo("wget --no-check-certificate -O %s %s" % (tarball,download_url))
    sudo("mv %s %s" % (tarball,install_dir))
    _untar(install_dir,tarball)

def _untar(install_dir,tarball):
    sudo("tar xvf %s/%s -C %s" % (install_dir,tarball,install_dir))

@task
def clean_up():
    sudo("rm -fr %(ROOT_DIR)s" % env)
