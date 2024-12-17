#!/bin/bash

# Run `setup-env.sh -h` for usage help

help() {
  echo "Setup an (LCG, compiler) view and, if on CERN servers, of the DIM DNS node name."
  echo "In case no LCG version or binary tag are specified, some defaults (with GCC) will be chosen for you."
  echo
  echo "Syntax: setup-env.sh [-h|l|b|c|r]"
  echo "options:"
  echo "h  Print this help."
  echo "l  LCG version to use (optional, e.g. 'LCG_107')."
  echo "b  Binary tag to use (optional, e.g. 'x86_64-el9-clang16-opt')."
  echo "c  Use Clang rather than GCC as compiler (used only if 'b' is not set)."
  echo "r  Use a local reduced setup of the LCG view rather than the standard one"
  echo "   (useful for CI tests, sometimes failing because of Herwig)"
  echo
  echo "Example: 'source setup-env.sh -l LCG_107 -b x86_64-el9-clang16-opt -r'"
  echo
}

# Read the options
clang=false
lcg_setup=false
binary_setup=false
use_custom_lcg_view=false
optstring=":hl:b:cr"
export OPTIND=1
while getopts ${optstring} opt; do
  case ${opt} in
    h) # Display help
      help
      return ;;
    c) # Setup Clang as compiler
      clang=true ;;
    l)
      lcg_setup=true
      export LCG_VERSION=${OPTARG} ;;
    b)
      binary_setup=true
      export BINARY_TAG=${OPTARG} ;;
    r) # Build LCG view locally with reduced number of libraries and dependencies
      use_custom_lcg_view=true ;;
    ?)
      echo "Unsupported option: ${OPTARG}"
      return ;;
  esac
done
if [[ $binary_setup = "true" && $clang = "true" ]] ; then
  echo "WARNING: The option -c is ignored as the binary tag ${BINARY_TAG} was set explicitly"
fi

# Function to setup a LCG view for a given compiler
setup_custom_lcg_view() {

  # Do not remove this - it allows for /cvmfs to be mounted elsewhere or to run this on views not installed on cvmfs.
  if [ -z $VIEWSDIR ]; then
      view_cvmfsPath=/cvmfs/sft.cern.ch/lcg/views
  else
      view_cvmfsPath=$VIEWSDIR
  fi

  #---Get the location this script (thisdir)
  SOURCE=${view_cvmfsPath}/${LCG_VERSION}/${BINARY_TAG}/setup.sh
  thisdir=$(cd "$(dirname "${SOURCE}")"; pwd)
  # Note: readlink -f is not working on OSX
  if [ "LINUX" = "OSX" ]; then
    thisdir=$(python3 -c "import os,sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" ${thisdir})
  else
    thisdir=$(readlink -f ${thisdir})
  fi

  # First the compiler
  if [ ${BINARY_TAG} == "x86_64-centos7-clang12-opt" ]; then
    folder="clang/12.0.0/x86_64-centos7"
  elif [ ${BINARY_TAG} == "x86_64-el9-clang16-opt" ]; then
    folder="clang/16.0.3/x86_64-el9"
  elif [ ${BINARY_TAG} == "x86_64-centos7-gcc8-opt" ]; then
    folder="gcc/8.3.0/x86_64-centos7"
  elif [ ${BINARY_TAG} == "x86_64-el9-gcc13-opt" ]; then
    folder="gcc/13.1.0/x86_64-el9"
  else
    echo "ERROR: The local configuration of binary tag ${BINARY_TAG} has not been set yet"
    return 1
  fi

  if [ "$COMPILER" != "native" ] && [ -e /cvmfs/sft.cern.ch/lcg/releases/${folder}/setup.sh ]; then
      source /cvmfs/sft.cern.ch/lcg/releases/${folder}/setup.sh
  fi

  #  then the rest...
  if [ -z "${PATH}" ]; then
      PATH=${thisdir}/bin; export PATH
  else
      PATH=${thisdir}/bin:$PATH; export PATH
  fi
  if [ -d ${thisdir}/scripts ]; then
      PATH=${thisdir}/scripts:$PATH; export PATH
  fi

  if [ -z "${LD_LIBRARY_PATH}" ]; then
      LD_LIBRARY_PATH=${thisdir}/lib; export LD_LIBRARY_PATH
  else
      LD_LIBRARY_PATH=${thisdir}/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
  fi
  if [ -d ${thisdir}/lib64 ]; then
      LD_LIBRARY_PATH=${thisdir}/lib64:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
  fi

  if [ -x ${thisdir}/bin/python ]; then
    PYTHON_VERSION=`expr $(readlink ${thisdir}/bin/python) : '.*/Python/\([0-9]\.[0-9]\+\).*'`
  else
    PYTHON_VERSION=`python3 -c "import sys; print('{0}.{1}'.format(*sys.version_info))"`
  fi
  PY_PATHS=${thisdir}/lib/python$PYTHON_VERSION/site-packages

  if [ -z "${PYTHONPATH}" ]; then
      PYTHONPATH=${thisdir}/lib:$PY_PATHS; export PYTHONPATH
  else
      PYTHONPATH=${thisdir}/lib:$PY_PATHS:$PYTHONPATH; export PYTHONPATH
  fi
  if [ -d ${thisdir}/python ]; then
      PYTHONPATH=${thisdir}/python:$PYTHONPATH; export PYTHONPATH
  fi

  if [ -z "${MANPATH}" ]; then
      MANPATH=${thisdir}/man:${thisdir}/share/man; export MANPATH
  else
      MANPATH=${thisdir}/man:${thisdir}/share/man:$MANPATH; export MANPATH
  fi
  if [ -z "${CMAKE_PREFIX_PATH}" ]; then
      CMAKE_PREFIX_PATH=${thisdir}; export CMAKE_PREFIX_PATH
  else
      CMAKE_PREFIX_PATH=${thisdir}:$CMAKE_PREFIX_PATH; export CMAKE_PREFIX_PATH
  fi
  if [ -z "${CPLUS_INCLUDE_PATH}" ]; then
      CPLUS_INCLUDE_PATH=${thisdir}/include; export CPLUS_INCLUDE_PATH
  else
      CPLUS_INCLUDE_PATH=${thisdir}/include:$CPLUS_INCLUDE_PATH; export CPLUS_INCLUDE_PATH
  fi
  if [ -z "${C_INCLUDE_PATH}" ]; then
      C_INCLUDE_PATH=${thisdir}/include; export C_INCLUDE_PATH
  else
      C_INCLUDE_PATH=${thisdir}/include:$C_INCLUDE_PATH; export C_INCLUDE_PATH
  fi

  #---check for compiler variables
  if [ -z "${CXX}" ]; then
      export FC=`command -v gfortran`
      export CC=`command -v gcc`
      export CXX=`command -v g++`
  fi

  #---Figure out the CMAKE_CXX_STANDARD (using Vc as a victim)
  if [ -f $thisdir/include/Vc/Vc ]; then
      vc_home=$(dirname $(dirname $(dirname $(readlink $thisdir/include/Vc/Vc))))
      std=$(cat $vc_home/logs/Vc*configure.cmake | \grep -Eo "CMAKE_CXX_STANDARD=[0-9]+" | \grep -Eo "[0-9]+")
      export CMAKE_CXX_STANDARD=$std
  fi

  #---then ROOT
  if [ -x $thisdir/bin/root ]; then
      if [ -x $thisdir/bin/python ]; then
          PYTHON_INCLUDE_PATH=$(dirname $(dirname $(readlink $thisdir/bin/python)))/include/$(\ls $(dirname $(dirname $(readlink $thisdir/bin/python)))/include)
      fi
      ROOTSYS=$(dirname $(dirname $(readlink $thisdir/bin/root))); export ROOTSYS
      if [ -z "${ROOT_INCLUDE_PATH}" ]; then
          ROOT_INCLUDE_PATH=${thisdir}/include:$PYTHON_INCLUDE_PATH; export ROOT_INCLUDE_PATH
      else
          ROOT_INCLUDE_PATH=${thisdir}/include:$PYTHON_INCLUDE_PATH:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
      fi
      if [ -d $thisdir/targets/x86_64-linux/include ]; then
          ROOT_INCLUDE_PATH=${thisdir}/targets/x86_64-linux/include:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
      fi
      if [ -z "${JUPYTER_PATH}" ]; then
          JUPYTER_PATH=${thisdir}/etc/notebook; export JUPYTER_PATH
      else
          JUPYTER_PATH=${thisdir}/etc/notebook:$JUPYTER_PATH; export JUPYTER_PATH
      fi
      export CPPYY_BACKEND_LIBRARY=$ROOTSYS/lib/libcppyy_backend${PYTHON_VERSION/./_}
      export CLING_STANDARD_PCH=none
  fi

  #---then Gaudi
  if [ -x $thisdir/include/Gaudi ]; then
      jsoninc=$(dirname $(dirname $(readlink ${thisdir}/include/nlohmann/json.hpp))) # see https://github.com/root-project/root/issues/7950
      ROOT_INCLUDE_PATH=${jsoninc}:${thisdir}/src/cpp:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
  fi
  if [ -x $thisdir/scripts/gaudirun.py ]; then
      Gaudi_DIR=$(dirname $(dirname $(readlink $thisdir/scripts/gaudirun.py)));
      export CMAKE_PREFIX_PATH=$Gaudi_DIR:$CMAKE_PREFIX_PATH
  fi

  #---then PYTHON
  # if [ -x $thisdir/bin/python ]; then
  #     PYTHONHOME=$(dirname $(dirname $(readlink $thisdir/bin/python))); export PYTHONHOME
  # elif [ -x $thisdir/bin/python3 ]; then
  #     PYTHONHOME=$(dirname $(dirname $(readlink $thisdir/bin/python3))); export PYTHONHOME
  # fi

  #---then Valgrind
  if [ -x $thisdir/libexec/valgrind/vgpreload_memcheck-amd64-linux.so ]; then
      VALGRIND_LIB=$thisdir/libexec/valgrind; export VALGRIND_LIB
  elif [ -x $thisdir/lib/valgrind/vgpreload_memcheck-amd64-linux.so ]; then
      VALGRIND_LIB=$thisdir/lib/valgrind; export VALGRIND_LIB
  fi

  #---then Graphviz
  if [ -f $thisdir/bin/dot ]; then
      GVBINDIR=$(dirname $(dirname $(readlink $thisdir/bin/dot)))/lib/graphviz; export GVBINDIR
  fi

  if [ -f $thisdir/lib/libQt5Gui.so ]; then
      export QT_PLUGIN_PATH=$(dirname $(dirname $(readlink $thisdir/lib/libQt5Gui.so )))/plugins
      export QT_XKB_CONFIG_ROOT=/usr/share/X11/xkb
  fi

  if [ -f $thisdir/etc/fonts/fonts.conf ]; then
      export FONTCONFIG_PATH=$thisdir/etc/fonts
  fi

  #---then tensorflow
  if [ -d $PY_PATHS/tensorflow_core ]; then   # version > 2.0
      export LD_LIBRARY_PATH=$PY_PATHS/tensorflow_core:$LD_LIBRARY_PATH
      TENSORFLOW=1
  elif [ -d $PY_PATHS/tensorflow ]; then
      export LD_LIBRARY_PATH=$PY_PATHS/tensorflow:$PY_PATHS/tensorflow/contrib/tensor_forest:$PY_PATHS/tensorflow/python/framework:$LD_LIBRARY_PATH
      TENSORFLOW=1
  fi

  #---then MySQL
  if [ -f $thisdir/bin/mariadb_config ]; then
      export MYSQL_HOME=$(dirname $(dirname $(readlink $thisdir/bin/mariadb_config )))
  fi

  #---then git
  if [ -f $thisdir/bin/git ]; then
      export GIT_EXEC_PATH=$(dirname $(dirname $(readlink $thisdir/bin/git )))/libexec/git-core
  fi

  #---then HTCondor
  if [ -x $thisdir/bin/condor_status ]; then
      PYTHONPATH=$(dirname $(dirname $(readlink $thisdir/bin/condor_status)))/lib/python3:$PYTHONPATH; export PYTHONPATH
  fi
}

# Setup the LCG version and binary tag ---------------------------------------------------------------------------------

echo "Will build Kepler with the following configuration:"
if grep -q "Red Hat Enterprise Linux 9\|AlmaLinux 9" /etc/os-release; then
  echo "  OS:           RHEL9/Alma9"
  if [ ${lcg_setup} = false ]; then
    export LCG_VERSION=LCG_107
  fi
  if [ ${binary_setup} = false ]; then
    if ${clang} ; then
      export BINARY_TAG=x86_64-el9-clang16-opt
    else
      export BINARY_TAG=x86_64-el9-gcc13-opt
    fi
  fi
else
  echo "The operating system of this computer is not supported by this script."
  return 1
fi
echo "  LCG version:  ${LCG_VERSION}"
echo "  Binary tag:   ${BINARY_TAG}"

# Build the local or standard LCG view ---------------------------------------------------------------------------------

if ${use_custom_lcg_view}; then
  echo "  LCG view:     local"
  setup_custom_lcg_view
else
  echo "  LCG view:     standard"
  source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh $LCG_VERSION $BINARY_TAG
fi
