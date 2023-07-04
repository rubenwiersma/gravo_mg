import os, sys, re, subprocess, platform

import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

__version__ = '0.0.1'

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.16.0':
                raise RuntimeError("CMake >= 3.16.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)


    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j3']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

def main():
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
    
    setuptools.setup(
        name='gravomg',
        version=__version__,
        author='Ahmad Nasikun and Ruben Wiersma',
        author_email="rubenwiersma@gmail.com",
        description='Python bindings for fast intrinsic multigrid, letting one script multiple experiments',
        long_description=long_description,
        long_description_content_type="text/markdown",
        license="MIT",
        url="https://github.com/rubenwiersma/gravo_mg",
        project_urls={
            "Bug Tracker": "https://github.com/rubenwiersma/gravo_mg/issues",
        },
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: OS Independent",
        ],
        package_dir = {'': 'src'},
        packages=setuptools.find_packages(where="src"),
        ext_modules=[CMakeExtension('.')],
        install_requires=['numpy', 'scipy'],
        cmdclass=dict(build_ext=CMakeBuild),
        zip_safe=False,
    )

if __name__ == "__main__":
    main()