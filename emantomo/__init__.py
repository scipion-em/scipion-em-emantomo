# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              David Herreros Calero (dherreros@cnb.csic.es)
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import subprocess

import pwem
import pyworkflow.utils as pwutils
from scipion.install.funcs import VOID_TGZ

from .constants import EMANTOMO_HOME, V2_39, COMMIT, CONDA_V2_39


_logo = "eman2_logo.png"
_references = ['Tang2007']
_url = "https://github.com/scipion-em/scipion-em-emantomo"


SCRATCHDIR = pwutils.getEnvVariable('EMANTOMOSCRATCHDIR', default='/tmp/')


class Plugin(pwem.Plugin):
    _homeVar = EMANTOMO_HOME
    _pathVars = [EMANTOMO_HOME]
    _supportedVersions = [V2_39]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(EMANTOMO_HOME, 'eman-' + V2_39)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Eman. """
        environ = pwutils.Environ(os.environ)
        pathList = [cls.getHome(d) for d in ['lib', 'bin']]

        # This environment variable is used to setup OpenGL (Mesa)
        # library in remote desktops
        if 'REMOTE_MESA_LIB' in os.environ:
            pathList.append(os.environ['REMOTE_MESA_LIB'])

        environ.update({'PATH': cls.getHome('bin')},
                       position=pwutils.Environ.BEGIN)

        environ.update({
            'LD_LIBRARY_PATH': os.pathsep.join(pathList),
            'PYTHONPATH': os.pathsep.join(pathList),
            'SCIPION_MPI_FLAGS': os.environ.get('EMANMPIOPTS', '')
        }, position=pwutils.Environ.REPLACE)

        return environ

    @classmethod
    def getActiveVersion(cls, home=None, versions=None):
        home = os.path.basename(home or cls.getHome())
        versions = versions or cls.getSupportedVersions()
        currentVersion = home.split('-')[-1]

        for v in versions:
            if v == currentVersion:
                return v

        return ''

    @classmethod
    def isVersion(cls, version=V2_39):
        return cls.getActiveVersion() == version

    @classmethod
    def getEmanActivation(cls, version=V2_39):
        return "conda activate emantomo-" + version

    @classmethod
    def getProgram(cls, program, python=False):
        """ Return the program binary that will be used. """
        # if cls.isVersion(V2_39):
        cmd = '%s %s && ' % (cls.getCondaActivationCmd(), cls.getEmanActivation())
        if python:
            python = subprocess.check_output(cmd + 'which python', shell=True).decode("utf-8")
            return '%(python)s %(program)s ' % locals()
        else:
            return cmd + '%(program)s ' % locals()
        # else:
        #     program = os.path.join(cls.getHome('bin'), program)
        #     if python:
        #         python = cls.getHome('bin/python')
        #         return '%(python)s %(program)s ' % locals()
        #     else:
        #         return '%(program)s ' % locals()

    @classmethod
    def getEmanCommand(cls, program, args, python=False):
        return cls.getProgram(program, python) + args

    @classmethod
    def getBoxerCommand(cls, boxerVersion='new'):
        cmd = 'e2boxer.py' if boxerVersion == 'new' else 'e2boxer_old.py'

        return os.path.join(cls.getHome('bin'), cmd)

    @classmethod
    def createEmanProcess(cls, script='e2converter.py', args=None, direc="."):
        """ Open a new Process with all EMAN environment (python...etc)
        that will server as an adaptor to use EMAN library
        """
        program = os.path.join(__path__[0], script)
        cmd = cls.getEmanCommand(program, args, python=True)

        print("** Running: '%s'" % cmd)
        cmd = cmd.split()
        proc = subprocess.Popen(cmd, env=cls.getEnviron(),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                cwd=direc,
                                universal_newlines=True)

        # Python 2 to 3 conversion: iterating over lines in subprocess stdout -> If universal_newlines is False the file
        # objects stdin, stdout and stderr will be opened as binary streams, and no line ending conversion is done.
        # If universal_newlines is True, these file objects will be opened as text streams in universal newlines mode
        # using the encoding returned by locale.getpreferredencoding(False). For stdin, line ending characters '\n' in
        # the input will be converted to the default line separator os.linesep. For stdout and stderr, all line endings
        # in the output will be converted to '\n'. For more information see the documentation of the io.TextIOWrapper
        # class when the newline argument to its constructor is None.

        return proc

    @classmethod
    def defineBinaries(cls, env):
        # SW_EM = env.getEmFolder()

        # shell = os.environ.get("SHELL", "bash")
        # eman23_commands = [
        #     (shell + ' ./emantomo.3.linux64.sh -b -p "%s/eman-2.3"' %
        #      SW_EM, '%s/eman-2.3/bin/python' % SW_EM)]
        # eman231_commands = [
        #     (shell + ' ./emantomo.31_sphire1.3.linux64.sh -b -p "%s/eman-2.31"' %
        #      SW_EM, '%s/eman-2.31/bin/python' % SW_EM)]

        # Eman Installation using Conda
        def getCondaInstallation(version=V2_39):
            installationCmd = cls.getCondaActivationCmd()
            installationCmd += 'conda create -y -n emantomo-' + version + ' --file ' + CONDA_V2_39 + ' && '
            installationCmd += 'conda activate emantomo-' + version + ' && '
            installationCmd += 'cd eman-build && '
            installationCmd += 'cmake ../eman-source/ -DENABLE_OPTIMIZE_MACHINE=ON && '
            installationCmd += 'make -j %d && make install' % env.getProcessors()
            return installationCmd

        # For Eman-2.39
        eman239_commands = []
        eman239_commands.append(('wget -c https://github.com/cryoem/eman2/archive/%s.tar.gz' % COMMIT, "2f7a976.tar.gz"))
        eman239_commands.append(("tar -xvf %s.tar.gz" % COMMIT, []))
        eman239_commands.append(("mv eman2*/ eman-source", "eman-source"))
        eman239_commands.append(('mkdir eman-build', 'eman-build'))
        installationCmd_239 = getCondaInstallation(V2_39)
        eman239_commands.append((installationCmd_239,
                                 "eman-build/libpyEM/CMakeFiles/pyPolarData2.dir/libpyPolarData2.cpp.o"))

        env.addPackage('eman', version=V2_39,
                       commands=eman239_commands,
                       tar=VOID_TGZ,
                       default=True)
