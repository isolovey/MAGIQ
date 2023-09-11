ARG MATLAB_VERSION=R2022a
FROM mathworks/matlab:${MATLAB_VERSION}
ARG MATLAB_VERSION
USER 0

ENV DEBIAN_FRONTEND=noninteractive

ENV PYTHON_VERSION=3.9.17
ENV PYTHON2_VERSION=2.7.18
ENV FSL_VERSION=6.0.7.1

ENV PYENV_ROOT=/opt/pyenv
ENV FSLDIR=/opt/fsl
ENV PATH="$PYENV_ROOT/shims:$PYENV_ROOT/bin:$FSL:$PATH"
ENV PYTHONPATH /app

# get system packages
RUN apt-get update && apt install -y wget

# Install IDL
RUN cd /opt && wget https://ftp.lowell.edu/pub/temp/old/IDL_Clients/v8.7.3/idl873-linux.tar.gz
RUN cd /opt && tar -xzvf idl873-linux.tar.gz &&  \
    (echo 'y'; echo '/opt/idl'; echo 'n'; echo 'n'; echo 'y'; echo 'y'; echo 'n';) | MORE=-99999999 ./install.sh

RUN apt-get update && apt install -y build-essential libssl-dev zlib1g-dev libbz2-dev \
libreadline-dev libsqlite3-dev curl llvm libncurses5-dev libncursesw5-dev \
xz-utils tk-dev libffi-dev liblzma-dev git qtchooser libxcb-icccm4-dev libxcb-image0-dev libxcb-keysyms1-dev \
libxcb-render-util0-dev libxcb-xfixes0-dev libxcb-xinerama0-dev libxcb-xkb-dev \
libxkbcommon-x11-dev

# set up pyenv
RUN cd /opt && git clone https://github.com/pyenv/pyenv pyenv
RUN eval "$(pyenv init -)"
RUN echo 'eval "$(pyenv init -)"' >> ~/.profile
RUN pyenv install -v ${PYTHON2_VERSION}
RUN pyenv install -v ${PYTHON_VERSION} && pyenv global ${PYTHON_VERSION}

# Install Qt5
run apt install -y qtbase5-dev qt5-qmake qttools5-dev qttools5-dev-tools qtwebengine5-dev libqt5svg5-dev libqt5websockets5-dev

# Install Python dependencies
RUN pyenv exec pip install --upgrade pip && pyenv exec pip install PyQt5 scipy matplotlib pyfftw nibabel pygamma future

# Install FSL
RUN cd /tmp && wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py && \
    chmod 775 fslinstaller.py && \
    pyenv exec python fslinstaller.py -d ${FSLDIR} -V ${FSL_VERSION}

# Install Python Matlab module
RUN cd /opt/matlab/${MATLAB_VERSION}/extern/engines/python && pyenv exec python setup.py install

# Install dcm2niix
RUN apt install -y dcm2niix pigz

# Install python libraries needed by scripts
RUN pyenv exec pip install dicomweb_client python-keycloak ipython pyyaml dirsync

WORKDIR /app

CMD python ./main.py
