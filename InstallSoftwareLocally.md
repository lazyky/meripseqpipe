# Dependencies for run lncPipe Locally

Prerequisites install command (required when docker image is not favored, you should execute them via root)

* [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)


		aria2c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -q -o /opt/hisat2-2.1.0-Linux_x86_64.zip && \
		unzip -qq /opt/hisat2-2.1.0-Linux_x86_64.zip -d /opt/ && \
		rm /opt/hisat2-2.1.0-Linux_x86_64.zip && \
		cd /opt/hisat2-2.1.0 && \
		rm -rf doc example *debug MANUAL* NEWS TUTORIAL && \
		ln -s /opt/hisat2-2.1.0/hisat2* /usr/local/bin/ && \
		ln -sf /opt/hisat2-2.1.0/*.py /usr/local/bin/

* [fastp](https://github.com/OpenGene/fastp)

        RUN aria2c http://opengene.org/fastp/fastp -q -o /usr/local/bin/fastp && \
            chmod a+x /usr/local/bin/fastp


* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)


		aria2c https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -q -o /opt/fastqc_v0.11.5.zip && \
		unzip -qq /opt/fastqc_v0.11.5.zip -d /opt/ && \
		rm /opt/fastqc_v0.11.5.zip && \
		cd /opt/FastQC && \
		shopt -s extglob && \
		rm -rfv !\("fastqc"\|*.jar\) && \
		chmod 755 * && \
		ln -s /opt/FastQC/fastqc /usr/local/bin/


**Alternatively, when you are going to using STAR-Cufflinks in your analysis, the corresponding install cmd are as follows:**

* [STAR](https://github.com/alexdobin/STAR)


		aria2c https://raw.githubusercontent.com/alexdobin/STAR/master/bin/Linux_x86_64/STAR -q -o /opt/STAR && \
		chmod 755 /opt/STAR && \
		ln -s /opt/STAR /usr/local/bin


* [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks)


		aria2c https://github.com/bioinformatist/cufflinks/releases/download/v2.2.1/cufflinks-2.2.1.Linux_x86_64.tar.gz -q -o /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
		tar xf /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz --use-compress-prog=pigz -C /opt/ && \
		rm /opt/cufflinks-2.2.1.Linux_x86_64/README && \
		ln -s /opt/cufflinks-2.2.1.Linux_x86_64/* /usr/local/bin/ && \
		rm /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz


> The `gffcompare` utility share the same function as `cuffcompare`, therefore, in STAR-cufflinks analysis pipe, `gffcompare` is not required.
