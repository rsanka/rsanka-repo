#
# vir-assembly-install.py
#  - Configures the environment for running vir-assembly-pipeline.sh (creating directory structure and installs software). 
#

import os.path, re
from fabric.api import cd, env, hide, local, run, settings, sudo, task
from fabric.network import disconnect_all

directories = {}
urls = {}
tarballs = {}

@task
def install():
    try:
        _initialize_script()
        _add_tools()
        _add_refs()
        _chmod("%(ROOT_DIR)s" % env)
    finally:
        disconnect_all()

def _initialize_script():
    print("user: %(user)s" % env)
    print("host: %(host)s" % env)
    
    env.ROOT_DIR = "/usr/local/VHTNGS"
    print("ROOT DIR: %(ROOT_DIR)s" % env)
    
    directories["PROJECT_DIR"] = "%(ROOT_DIR)s/project" % env
    directories["REF_DIR"] = "%(ROOT_DIR)s/references" % env
    directories["TOOLS_DIR"] = "%(ROOT_DIR)s/tools" % env
    directories["TOOLS_BINARIES_DIR"] = "%s/BINARIES" % directories["TOOLS_DIR"]
    directories["TOOLS_PERL_DIR"] = "%s/PERL" % directories["TOOLS_DIR"]
    directories["TOOLS_FASTX_DIR"] = "%s/FASTX" % directories["TOOLS_DIR"]
    
    for name in sorted(directories.keys()):
        sudo("mkdir -p %s" % directories[name])
        print("%s: %s" % (name,directories[name]))
    
    env.VIR_ASSEMBLY_SCRIPT = "https://raw.github.com/rsanka/rsanka-repo/master/vir-assembly-pipeline.sh"
    print("VIR ASSEMBLY SCRIPT: %(VIR_ASSEMBLY_SCRIPT)s" % env)
    
    env.REF_FILES = "corona_virus,hadv,influenza_a_virus,jev,mpv,norv,rota_virus,rsv,veev,vzv,yfv"
    print("REFS FILES: %(REF_FILES)s" % env)
    
    urls["BINARIES_URL"] = "https://raw.github.com/rsanka/rsanka-repo/master/tools/BINARIES/BINARIES.tgz"
    urls["PERL_URL"] = "https://raw.github.com/rsanka/rsanka-repo/master/tools/PERL/PERL.tgz"
    urls["FASTX_URL"] = "http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2"
    urls["BIO_LINUX_URL"] = "http://nebc.nerc.ac.uk/bio-linux/"
    urls["REF_URL"] = "https://github.com/downloads/rsanka/rsanka-repo"
    
    for name in sorted(urls.keys()):
        sudo("mkdir -p %s" % urls[name])
        print("%s: %s" % (name,urls[name]))
    
    tarballs["BINARIES_TARBALL"] = "BINARIES.tgz"
    tarballs["PERL_TARBALL"] = "PERL.tgz"
    tarballs["FASTX_TARBALL"] = "fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2"
    
    for name in sorted(tarballs.keys()):
        print("%s: %s" % (name,tarballs[name]))

def _add_tools():
    sudo ("echo -e \"DEBIAN_FRONTEND=noninteractive\" >> /home/ubuntu/.bashrc")
    sudo("wget --no-check-certificate -O %s/vir-assembly-pipeline.sh %s" % (env.ROOT_DIR,env.VIR_ASSEMBLY_SCRIPT))
    _add_tarball(urls["BINARIES_URL"],tarballs["BINARIES_TARBALL"],directories["TOOLS_BINARIES_DIR"],"xfz")
    _add_tarball(urls["PERL_URL"],tarballs["PERL_TARBALL"],directories["TOOLS_PERL_DIR"],"xfz")
    _apt_get_install("csh")
    _apt_get_install("gawk")
    _add_fastx()
    _initialize_bio_linux()

def _add_refs():
    files = (env.REF_FILES).split(",")
    for file in files:
        _add_tarball("%s/%s.tgz" % (urls["REF_URL"],file),file,directories["REF_DIR"],"xfz")

def _chmod(path):
    sudo("chmod -R 755 %s" % path)

def _add_tarball(download_url,tarball,install_dir,options):
    sudo("wget --no-check-certificate -O %s %s" % (tarball,download_url))
    sudo("mv %s %s" % (tarball,install_dir))
    _untar(install_dir,tarball,options)

def _add_fastx():
    _add_tarball(urls["FASTX_URL"],tarballs["FASTX_TARBALL"],directories["TOOLS_FASTX_DIR"],"xfj")
    sudo("cp %s/bin/* %s" % (directories["TOOLS_FASTX_DIR"],directories["TOOLS_FASTX_DIR"]))
    sudo("rm -fr %s/bin %s/%s" % (directories["TOOLS_FASTX_DIR"],directories["TOOLS_FASTX_DIR"],tarballs["FASTX_TARBALL"]))

def _initialize_bio_linux():
    sudo("echo -e \"deb %s unstable bio-linux\" >> /etc/apt/sources.list" % urls["BIO_LINUX_URL"])
    sudo("sudo apt-get update")
    _apt_get_install("bio-linux-keyring")
    _apt_get_install("bwa")
    _apt_get_install("samtools")
    _apt_get_install("bio-linux-cap3")
    _apt_get_install("emboss")

def _untar(install_dir,tarball,options):
    sudo("tar %s %s/%s -C %s" % (options,install_dir,tarball,install_dir))

def _apt_get_install(tool):
    sudo("apt-get -q -y --force-yes install %s" % tool)

@task
def clean_up():
    sudo("rm -fr %(ROOT_DIR)s" % env)
