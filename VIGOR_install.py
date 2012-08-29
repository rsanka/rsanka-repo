import os.path, re
from fabric.api import cd, env, hide, local, run, settings, sudo, task
from fabric.network import disconnect_all

    
def install_cloudvigor():
    try:
        _initialize_script()
        _initialize_host()
        _install_vigor()
        _install_tools()
    finally:
        disconnect_all()

def install_testcloudvigor():
    try:
        _initialize_script()

        cmd = ("""%(VIGOR_RUNTIME_DIR)s/VIGOR3.pl \
                -D yfv \
                -i %(VIGOR_SAMPLE_DATA_DIR)s/westnile.fasta \
                -O %(VIGOR_TEST_OUTPUT_DIR)s/westnile \
                > %(VIGOR_SCRATCH_DIR)s/westnile_test_run.log 2>&1 \
                """) % env
        print("DEBUG: cmd[%s]" % cmd)
        run(cmd)

    finally:
        disconnect_all()

def install_validatevigor():
    try:
        _initialize_script()
        sudo("rm -f %s(VIGOR_TEST_OUTPUT_DIR)s/westnile.rpt" % env)
        with settings(hide("running","stdout")):
            results = run("""diff -Bwr %(VIGOR_VALIDATION_TEST_DATA_DIR)s \
                          %(VIGOR_TEST_OUTPUT_DIR)s \
                          || echo 'VALIDATION FAILED'""" % env)
        if results:
            print("\n\nValidation Failed:\n\n%s\n" % results)
    finally:
        disconnect_all()


def install_cleanall():
    try:
        _initialize_script()
        _initialize_host()
        _remove_vigor()
        _remove_tools()
        _remove_dir(env.SCRATCH_DIR)
    finally:
        disconnect_all()


def _create_scratch_dir():
    if not _path_is_dir(env.VIGOR_SCRATCH_DIR):
        sudo("mkdir -p %(VIGOR_SCRATCH_DIR)s" % env)
    sudo("chown -R %(user)s:%(user)s %(VIGOR_SCRATCH_DIR)s" % env)
    sudo("find %(VIGOR_SCRATCH_DIR)s -type d -exec chmod 775 {} \;" % env)

def _create_tools_dir():
    if not _path_is_dir(env.TOOLS_DIR):
        sudo("mkdir -p %(TOOLS_DIR)s" % env)
    sudo("chown -R %(user)s:%(user)s %(TOOLS_DIR)s" % env)

def _create_vigor_scratch_dir():
    if not _path_is_dir(env.VIGOR_SCRATCH_DIR):
        sudo("mkdir -p %(VIGOR_SCRATCH_DIR)s" % env)
    sudo("find %(VIGOR_SCRATCH_DIR)s -type f -exec chmod 666 {} \;" % env)
    sudo("find %(VIGOR_SCRATCH_DIR)s -type d -exec chmod 777 {} \;" % env)

def _create_vigor_tempspace_dir():
    if not _path_is_dir(env.VIGOR_TEMPSPACE_DIR):
        sudo("mkdir -p %(VIGOR_TEMPSPACE_DIR)s" % env)
    sudo("chown -R %(user)s:%(user)s %(VIGOR_TEMPSPACE_DIR)s" % env)
    sudo("find %(VIGOR_TEMPSPACE_DIR)s -type d -exec chmod 777 {} \;" % env)

def _fix_etc_hosts():
    internal_ip = sudo("hostname")
    print("internal_ip[%s]" % internal_ip)
    filespec = "/etc/hosts"
    sudo("echo '127.0.0.1 %s' >> %s" % (internal_ip, filespec))

def _initialize_host():
    local("ssh-keygen -R %(host)s" % env)
    _fix_etc_hosts()
    _create_scratch_dir()

def _initialize_script():

    machine = run("uname -m")

    if machine.find('64')>0:
        env.ARCH = 'x64-linux'
    else:
        env.ARCH = 'ia32-linux'

    env.ROOT_DIR = '/usr/local'
    env.SCRATCH_DIR = '/usr/local/scratch'
    env.GITHUB_URL = 'https://github.com/downloads/JCVI-Cloud/cloudbiolinux'
    env.BLAST_NAME = 'blast-2.2.15'
    env.CLUSTALW_NAME = 'clustalw-1.83'
    env.VIGOR_NAME = 'vigor-GSCcloud-release-20120828'

    print("user[%(user)s]" % env)
    print("host[%(host)s]" % env)
    print("ROOT_DIR[%(ROOT_DIR)s]" % env)
    print("SCRATCH_DIR[%(SCRATCH_DIR)s]" % env)
    print("GITHUB_URL[%(GITHUB_URL)s]" % env)
    print("BLAST_NAME[%(BLAST_NAME)s]" % env)
    print("CLUSTALW_NAME[%(CLUSTALW_NAME)s]" % env)
    print("VIGOR_NAME[%(VIGOR_NAME)s]" % env)
    print("ARCH[%(ARCH)s]" % env)
          
    env.HOME = os.path.join("/home", env.user)

    env.TOOLS_DIR = os.path.join(env.ROOT_DIR, "tools")
    env.VIGOR_STORED_DIR = os.path.join(env.TOOLS_DIR, "vigor")
    env.VIGOR_RUNTIME_DIR = os.path.join(env.VIGOR_STORED_DIR, "prod3")
    env.VIGOR_SCRATCH_DIR = os.path.join(env.SCRATCH_DIR, "vigor")
    env.VIGOR_TEMPSPACE_DIR = os.path.join(env.VIGOR_SCRATCH_DIR, "tempspace")
    env.VIGOR_SAMPLE_DATA_DIR = os.path.join(env.VIGOR_STORED_DIR, "samples")
    env.VIGOR_TAR_FILENAME = "%(VIGOR_NAME)s.tgz" % env
    env.VIGOR_TEST_OUTPUT_DIR = os.path.join(env.VIGOR_STORED_DIR, "test")
    env.VIGOR_VALIDATION_TEST_DATA_DIR = os.path.join(env.VIGOR_STORED_DIR, "samples_valid_output")
    env.BLAST_DIR = os.path.join(env.TOOLS_DIR, "blast")
    env.CLUSTALW_DIR = os.path.join(env.TOOLS_DIR, "clustalw")
    env.EXE_DIR = env.VIGOR_RUNTIME_DIR
    env.BLAST_TAR_FILENAME = "%(BLAST_NAME)s-%(ARCH)s.tar.gz" % env
    env.CLUSTALW_TAR_FILENAME = "%(CLUSTALW_NAME)s-%(ARCH)s.deb" % env

def _install_blast():
    print("    Installing blast...")
    _create_tools_dir()
    _install_package(env.GITHUB_URL, env.BLAST_TAR_FILENAME, env.BLAST_DIR, True)
    if not _path_exists(os.path.join(env.EXE_DIR, "blastall")):
        sudo("ln -s %s %s" % (os.path.join(env.BLAST_DIR, env.BLAST_NAME, "bin", "bl2seq"), env.EXE_DIR))
        sudo("ln -s %s %s" % (os.path.join(env.BLAST_DIR, env.BLAST_NAME, "bin", "blastall"), env.EXE_DIR))
        sudo("ln -s %s %s" % (os.path.join(env.BLAST_DIR, env.BLAST_NAME, "bin", "fastacmd"), env.EXE_DIR))
        sudo("ln -s %s %s" % (os.path.join(env.BLAST_DIR, env.BLAST_NAME, "bin", "formatdb"), env.EXE_DIR))

def _install_clustalw():
    print("    Installing clustalw...")
    _create_tools_dir()
    _install_package(env.GITHUB_URL, env.CLUSTALW_TAR_FILENAME,
                     env.CLUSTALW_DIR, False)
    if not _path_exists(os.path.join(env.EXE_DIR, "clustalw")):
            sudo("ln -s %s %s"
                 % (os.path.join(env.CLUSTALW_DIR, env.CLUSTALW_NAME, "clustalw"),
                    env.EXE_DIR))

def _install_tools():
    print("Install tools...")
    _create_tools_dir()
    _install_blast()
    _install_clustalw()
    _apt_get_install("libapache-dbi-perl")
    _apt_get_install("libclass-dbi-sqlite-perl")

def _install_vigor():
    print("Installing VIGOR...")
    _create_vigor_tempspace_dir()
    _create_vigor_scratch_dir()
    _install_package(env.GITHUB_URL, env.VIGOR_TAR_FILENAME,
                     env.VIGOR_STORED_DIR, True)
    sudo("chmod 755 %s" % os.path.join(env.VIGOR_RUNTIME_DIR, "*.pl"))
    if not _path_exists(os.path.join(env.EXE_DIR, "perl")):
        sudo("ln -s %s %s" % ("/usr/bin/perl", env.EXE_DIR))
        sudo("ln -s %s %s" % ("/usr/bin/perl", "/usr/local/bin"))
    if not _path_exists(os.path.join(env.EXE_DIR, "vigorscratch")):
        sudo("ln -s %s %s/vigorscratch" % (env.VIGOR_TEMPSPACE_DIR, env.EXE_DIR))

def _remove_blast():
    print("    Removing blast...")
    _remove_symlinks(os.path.join(env.BLAST_DIR, env.BLAST_NAME, "bin", "*"),
                     env.EXE_DIR)
    _remove_dir(env.BLAST_DIR)

def _remove_clustalw():
    print("    Removing clustalw...")
    _remove_symlinks(os.path.join(env.CLUSTALW_DIR, env.CLUSTALW_NAME, "*"),
                     env.EXE_DIR)
    _remove_dir(env.CLUSTALW_DIR)

def _remove_tools():
    print("Removing tools...")
    _remove_blast()
    _remove_clustalw()
    _remove_dir(env.TOOLS_DIR)

def _remove_vigor():
    print("Removing VIGOR...")
    _remove_dir(env.VIGOR_RUNTIME_DIR)
    _remove_dir(env.VIGOR_SAMPLE_DATA_DIR)
    _remove_dir(env.VIGOR_TEMPSPACE_DIR)
    _remove_dir(env.VIGOR_TEST_OUTPUT_DIR)
    _remove_dir(env.VIGOR_VALIDATION_TEST_DATA_DIR)
    _remove_dir(env.VIGOR_SCRATCH_DIR)

# Utility methods.

def _install_package(download_url, filename, install_dir, tar):
    if not _path_is_dir(install_dir):
        sudo("mkdir -p %s" % install_dir)
    with cd(install_dir):
        if not _path_exists(os.path.join(install_dir, filename)):
            sudo("""wget --no-host-directories --cut-dirs=1 \
                 --directory-prefix=%s %s/%s"""
                % (install_dir, download_url, filename))
            if tar:
                sudo("tar xvfz %s" % filename)
            else: 
                sudo("dpkg -x %s %s" % (filename,install_dir))
                sudo("mkdir %s/%s" % (install_dir, env.CLUSTALW_NAME))
                sudo("cp %s/usr/bin/* %s/%s" % (install_dir,install_dir,env.CLUSTALW_NAME))
    sudo("chown -R %s:%s %s" % (env.user, env.user, install_dir))
    sudo("find %s -type d -exec chmod 755 {} \;" % install_dir)

def _path_exists(path):
    found = False
    with settings(hide("running","stdout")):
        result = sudo("test -e '%s' || echo 'FALSE'" % path)
    if result != "FALSE": found = True
    return found

def _path_is_dir(path):
    found = False
    with settings(hide("running","stdout")):
        result = sudo("test -d '%s' || echo 'FALSE'" % path)
    if result != "FALSE": found = True
    return found

def _remove_dir(dirspec):
    if _path_is_dir(dirspec):
        _unlock_dir(dirspec)
        sudo("rm -rf %s" % dirspec)
    else:
        print("DEBUG: _remove_dir[%s] -- NOT FOUND" % dirspec)

def _remove_symlinks(link_from_filespec, link_to_dir):
    if _path_is_dir(link_to_dir):
        sudo("find %s -lname '%s' -delete" % (link_to_dir, link_from_filespec))

def _unlock_dir(dirspec):
    with settings(hide("running","stdout")):
        sudo("find %s -type d -exec chmod 755 {} \;" % dirspec)
        sudo("find %s -type d -exec chmod g+s {} \;" % dirspec)
        sudo("find %s -type f -exec chmod 644 {} \;" % dirspec)

def _apt_get_install(tool):
    sudo("apt-get install %s" % tool)
