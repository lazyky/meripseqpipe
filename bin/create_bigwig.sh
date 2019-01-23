#!/bin/bash
#$1 argv 1 : THREAD_NUM
THREAD_NUM=$1
#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done

for bam_file in *.bam
do
read -u 9
{
    bamCoverage -b $bam_file -o ${bam_file%.bam*}.bigwig
}& 
done
