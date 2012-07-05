#
# vir-assembly-install.py
#  - Configures the environment for running vir-assembly-pipeline.sh (creating directory structure and installs software). 
#

import os.path, re
from fabric.api import cd, env, hide, local, run, settings, sudo, task
from fabric.network import disconnect_all

directories = {}

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
    directories["PROJECT_DIR"] = "%(ROOT_DIR)s/project" % env
    directories["REF_DIR"] = "%(ROOT_DIR)s/references" % env
    directories["TOOLS_DIR"] = "%(ROOT_DIR)s/tools" % env
    directories["TOOLS_BINARIES_DIR"] = "%s/BINARIES" % directories["TOOLS_DIR"]
    directories["TOOLS_PERL_DIR"] = "%s/PERL" % directories["TOOLS_DIR"]
    directories["TOOLS_RUBYLIB_DIR"] = "%s/RUBYLIB" % directories["TOOLS_DIR"]

def _print_env_variables():
    print("user:                %(user)s" % env)
    print("host:                %(host)s" % env)
    print("ROOT DIR:            %(ROOT_DIR)s" % env)
    print("REFS FILES:          %(REF_FILES)s" % env)

def _initialize_dirs():
    for name in sorted(directories.keys()):
        sudo("mkdir -p %s" % directories[name])

def _add_tools():
    _add_tarball(env.BINARIES_URL,env.BINARIES_TARBALL,directories["TOOLS_BINARIES_DIR"])
    _add_tarball(env.PERL_URL,env.PERL_TARBALL,directories["TOOLS_PERL_DIR"])
    _add_tarball(env.RUBYLIB_URL,env.RUBYLIB_TARBALL,directories["TOOLS_RUBYLIB_DIR"])
    _initialize_bio_linux()

def _add_refs():
    files = (env.REF_FILES).split(",")
    for file in files:
        _add_tarball("%s/%s.tgz" % (env.REF_URL,file),file,directories["REF_DIR"])

def _chmod(path):
    sudo("chmod -R 755 %s" % path)

def _add_tarball(download_url,tarball,install_dir):
    print("du: %s" % download_url)
    print("tf: %s" % tarball)
    print("id: %s" % install_dir)
    sudo("wget --no-check-certificate -O %s %s" % (tarball,download_url))
    sudo("mv %s %s" % (tarball,install_dir))
    _untar(install_dir,tarball)

def _initialize_bio_linux():
    sudo("echo -e \"deb http://nebc.nerc.ac.uk/bio-linux/ unstable bio-linux\" >> /etc/apt/sources.list")
    sudo("sudo apt-get update")
    sudo("sudo apt-get install bio-linux-keyring")
    _apt_get_install("bowtie")
    _apt_get_install("samtools")
    _apt_get_install("bio-linux-cap3")
    _apt_get_install("emboss")

def _untar(install_dir,tarball):
    sudo("tar xfz %s/%s -C %s" % (install_dir,tarball,install_dir))

def _apt_get_install(tool):
    sudo("apt-get install %s" % tool)

@task
def clean_up():
    sudo("rm -fr %(ROOT_DIR)s" % env)
