## Conda
```
git clone https://github.com/kingzhuky/m6APipe.git && \
cd m6APipe && \
conda env create -f environment.yml && \
conda activate nf-core-m6APipe-1.0dev
```
Waiting for the conda environment of m6APipe to be built.
### Install MATK
```
wget http://matk.renlab.org/download/MATK-1.0.jar
```
### Install QNB
```
wget https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz && \
    R CMD INSTALL QNB_1.1.11.tar.gz && \
    rm QNB_1.1.11.tar.gz
```
### Install MeTDiff
```
git clone https://github.com/compgenomics/MeTDiff.git && \
    R CMD build MeTDiff/ && \
    R CMD INSTALL MeTDiff_1.0.tar.gz && \
    rm -rf MeTDiff*
```
### install MeTPeak
```
git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD build MeTPeak/ && \
    R CMD INSTALL MeTPeak_1.0.0.tar.gz && \
    rm -rf MeTPeak*
```
### install MSPC
Download [MSPC](https://github.com/Genometric/MSPC) and Add it into enviroment variable $PATH
```
wget -O mspc.zip "https://github.com/Genometric/MSPC/releases/latest/download/mspc.zip" && \
    unzip mspc.zip -d mspc
echo "export PATH=`pwd`/mspc:\$PATH" > ~/.bashrc
```