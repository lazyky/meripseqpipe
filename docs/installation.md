## Installation
### Install Nextflow
Install [Nextflow](https://www.nextflow.io/) by using the following command:
```
curl -s https://get.nextflow.io | bash 
```
### Install m6APipe
```
git clone https://github.com/kingzhuky/m6APipe.git
```
### Build environment by conda or  docker
#### Build environment by docker
Using the following command for buildind envrionment:
```
docker pull kingzhuky/m6apipe
```
#### Build environment by conda
Using the following command for buildind envrionment:
```
cd m6APipe && \
conda env create -f environment.yml && \
conda activate nf-core-m6APipe-1.0dev
```
Waiting for the conda environment of m6APipe to be built.
##### Install MATK
[MATK](http://matk.renlab.org)  is one of the PeakCalling tools which arranged by m6APipe. ( options )
```
wget http://matk.renlab.org/download/MATK-1.0.jar
```
##### Install MeTPeak
[MeTPeak](https://github.com/compgenomics/MeTPeak) is one of the PeakCalling tools which arranged by m6APipe. ( options )
```
git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD INSTALL MeTPeak/ && \
    rm -rf MeTPeak*
```
##### Install QNB
[QNB](https://cran.r-project.org/src/contrib/Archive/QNB/) is one of the methylation analysis methods which arranged by m6APipe. ( options )
```
wget https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz && \
    R CMD INSTALL QNB_1.1.11.tar.gz && \
    rm QNB_1.1.11.tar.gz
```
##### Install MeTDiff
[MeTDiff](https://github.com/compgenomics/MeTDiff) is one of the methylation analysis methods which arranged by m6APipe. ( options )
```
git clone https://github.com/compgenomics/MeTDiff.git && \
    R CMD INSTALL MeTDiff/ && \
    rm -rf MeTDiff*
```

##### install MSPC
Download [MSPC](https://github.com/Genometric/MSPC) and Add it into enviroment variable $PATH
```
wget -O mspc.zip "https://github.com/Genometric/MSPC/releases/download/v4.0.0/linux-x64.zip" && \
    unzip mspc.zip -d mspc
chmod +x mspc/mspc
echo "export PATH=`pwd`/mspc:\$PATH" >> ~/.bashrc
```