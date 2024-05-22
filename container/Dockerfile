FROM ubuntu:22.04 
# \ for seq of commands | for piping
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
     unzip xorg wget \
 && apt-get clean \ 
 && rm -rf \
     /tmp/hsperfdata* \
     /var/*/apt/*/partial \
     /var/lib/apt/lists/* \
     /var/log/apt/term*
# cleans file defs
# Install MATLAB MCR in /opt/mcr/
# check the below info with the matworks website for 2023a
ENV MATLAB_VERSION R2023a
ENV MCR_VERSION R2023a
ENV MCR_UPDATE 6
RUN mkdir /opt/mcr_install \
 && mkdir /opt/mcr \
 && wget --progress=bar:force -P /opt/mcr_install https://ssd.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/Release/${MCR_UPDATE}/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_${MATLAB_VERSION}_Update_${MCR_UPDATE}_glnxa64.zip \
 && unzip -q /opt/mcr_install/MATLAB_Runtime_${MATLAB_VERSION}_Update_${MCR_UPDATE}_glnxa64.zip -d /opt/mcr_install \
 && /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent \
 && rm -rf /opt/mcr_install /tmp/*

# Install hmri-standalone
ENV SPM_VERSION 12
ENV PATH="$PATH:/opt/standalone-hMRItoolboxv0.6.1"
ENV SPM_REVISION r7771
ENV LD_LIBRARY_PATH /opt/mcr/${MCR_VERSION}/runtime/glnxa64:/opt/mcr/${MCR_VERSION}/bin/glnxa64:/opt/mcr/${MCR_VERSION}/sys/os/glnxa64:/opt/mcr/${MCR_VERSION}/sys/opengl/lib/glnxa64:/opt/mcr/${MCR_VERSION}/extern/bin/glnxa64
ENV MCR_INHIBIT_CTF_LOCK 1
ENV SPM_HTML_BROWSER 0
# Running SPM once with "function exit" tests the succesfull installation *and*
# extracts the ctf archive which is necessary if singularity is going to be
# used later on, because singularity containers are read-only.
# Also, set +x on the entrypoint for non-root container invocations
RUN wget --no-check-certificate --progress=bar:force -P /opt https://github.com/hMRI-group/hMRI-toolbox/releases/download/v0.6.1/standalone-hMRItoolboxv0.6.1.zip \
 && unzip -q /opt/standalone-hMRItoolboxv0.6.1.zip -d /opt \
 && rm -f /opt/standalone-hMRItoolboxv0.6.1.zip \
 && /opt/standalone-hMRItoolboxv0.6.1/spm${SPM_VERSION} function exit \
 && chmod +x /opt/standalone-hMRItoolboxv0.6.1/spm${SPM_VERSION} \
 && chmod a+r /opt/standalone-hMRItoolboxv0.6.1 -R 


# Configure entry point
#ENTRYPOINT ["/opt/spm12/spm12"]
#CMD ["--help"]

# in root need to run: /opt/standalone-hMRItoolboxv0.6.1/spm12