#
# SISPA-install-fabfile.py
#

import os.path, re
from fabric.api import cd, env, hide, local, run, settings, sudo, task
from fabric.network import disconnect_all

@task
def install():
    try:
        _print_env_variables()
        _initialize_dirs()
        _add_tools()
    finally:
        disconnect_all()

def _print_env_variables():
    print("user: %(user)s" % env)
    print("host: %(host)s" % env)
    print("PROJECT ROOT: %(PROJECT_ROOT_DIR)s" % env)
    print("SCRATCH ROOT: %(SCRATCH_ROOT_DIR)s" % env)
    print("TOOLS ROOT: %(TOOLS_ROOT_DIR)s" % env)

def _initialize_dirs():
    sudo("mkdir -p %(PROJECT_ROOT_DIR)s" % env)
    sudo("mkdir -p %(SCRATCH_ROOT_DIR)s" % env)
    sudo("mkdir -p %(TOOLS_ROOT_DIR)s" % env)
    sudo("mkdir -p %(RUBYLIB_DIR)s" % env)

def _add_tools():
    _add_rubylib(env.RUBYLIB_URL,env.RUBYLIB_DIR)

def _add_rubylib(download_url,install_dir):
    with cd(install_dir):
        print("du: %s" % download_url)
        print("id: %s" % install_dir)
        sudo("wget --no-check-certificate -O RUBYLIB.tar --directory-prefix=%s %s" % (install_dir,download_url))
        sudo("tar xvf %s/RUBYLIB.tar -C %s" % (install_dir,install_dir))

@task
def clean_up():
    sudo("rm -fr %(PROJECT_ROOT_DIR)s" % env)
    sudo("rm -fr %(SCRATCH_ROOT_DIR)s" % env)
    sudo("rm -fr %(TOOLS_ROOT_DIR)s" % env)
