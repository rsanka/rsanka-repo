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
        _chmod("%(ROOT_DIR)s" % env)
    finally:
        disconnect_all()

def _define_dirs():
    env.PROJECT_DIR = "%(ROOT_DIR)s/project" % env
    env.SAMPLE_DATA_DIR = "%(PROJECT_DIR)s/sample_data" % env
    env.REF_FILES = "%(REF_FILES)s" % env
    env.REF_DIR = "%(ROOT_DIR)s/refs" % env
    env.TOOLS_DIR = "%(ROOT_DIR)s/tools" % env
    env.TOOLS_RUBYLIB_DIR = "%(TOOLS_DIR)s/RUBYLIB" % env
    env.TOOLS_BINARIES_DIR = "%(TOOLS_DIR)s/BINARIES" % env

def _print_env_variables():
    print("user:                %(user)s" % env)
    print("host:                %(host)s" % env)
    print("ROOT DIR:            %(ROOT_DIR)s" % env)
    print("PROJECT DIR:         %(PROJECT_DIR)s" % env)
    print("SAMPLE DATA DIR:     %(SAMPLE_DATA_DIR)s" % env)
    print("REFS FILES:          %(REF_FILES)s" % env)
    print("REFS DIR:            %(REF_DIR)s" % env)
    print("TOOLS DIR:           %(TOOLS_DIR)s" % env)
    print("TOOLS RUBYLIB DIR:   %(TOOLS_RUBYLIB_DIR)s" % env)
    print("TOOLS BINARIES DIR:  %(TOOLS_BINARIES_DIR)s" % env)

def _initialize_dirs():
    sudo("mkdir -p %(PROJECT_DIR)s" % env)
    sudo("mkdir -p %(SAMPLE_DATA_DIR)s" % env)
    sudo("mkdir -p %(REF_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_RUBYLIB_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_BINARIES_DIR)s" % env)

def _add_tools():
    _add_tarball(env.RUBYLIB_URL,env.RUBYLIB_TARBALL,env.TOOLS_RUBYLIB_DIR)
    _add_tarball(env.BINARIES_URL,env.BINARIES_TARBALL,env.TOOLS_BINARIES_DIR)

def _add_refs():
    files = (env.REF_FILES).split(",")
    for file in files:
        sudo("mkdir -p %s/%s" % (env.SAMPLE_DATA_DIR,file))
        _add_tarball("%s/%s.tgz" % (env.REF_URL,file),file,env.REF_DIR)

def _chmod(path):
    sudo("chmod -R 755 %s" % path)

def _add_tarball(download_url,tarball,install_dir):
    print("du: %s" % download_url)
    print("tf: %s" % tarball)
    print("id: %s" % install_dir)
    sudo("wget --no-check-certificate -O %s %s" % (tarball,download_url))
    sudo("mv %s %s" % (tarball,install_dir))
    _untar(install_dir,tarball)

def _untar(install_dir,tarball):
    sudo("tar xfz %s/%s -C %s" % (install_dir,tarball,install_dir))

@task
def clean_up():
    sudo("rm -fr %(ROOT_DIR)s" % env)
