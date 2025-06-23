library(stringr)
args = commandArgs(trailingOnly=TRUE)

#Parameter Check: 

if(length(args)<3) {
	cat("Program aborted - USAGE:\n window_SHAPE <Input alignment file in FASTA format> <window size> <window overlap> \n OPTIONAL: <SHAPE file for 1st sequence> <SHAPE file for 2nd sequence> <SHAPE file for ...>\n\n Parameters must be entered in correct order. Output is printed to file.")
	quit()
}


aln=read.table(args[1])
size=as.numeric(args[2])
overlap=as.numeric(args[3])
name1=gsub(".fasta","",args[1])

#aln=read.table("../0_raw/hiv_37_sampled_gaps.fasta")
#size=200
#overlap=100
#name1="../0_raw/hiv_37_sampled_gaps.fasta"




alnlen=str_count(aln$V1[2])
seqs=length(aln$V1)/2
windows=ceiling((alnlen-size)/overlap)+1

#Create windows of alignment
for (i in (1:windows)){
	name2=paste(name1,"_w_",i,".fasta",sep="")
	start=overlap*(i-1)+1
	stop=size+overlap*(i-1)
	if (stop>alnlen) 
		stop=alnlen
	for (i in (c(1:seqs)*2)) {
		write(file=name2,as.character(aln$V1[i-1]),append=TRUE)
		write(file=name2,substring(aln$V1[i],start,stop),append=TRUE)
	}
}


#Create SHAPE files 
a=4
while (a<=length(args)) {
	shape=read.table(args[a],skip=1)
	name3=gsub(".txt","",args[a])
	len=alnlen-str_count(aln$V1[2*(a-3)],"-")
	for (i in (1:windows)){
		name4=paste(name3,"_w_",i,".txt",sep="")
		start=overlap*(i-1)
		stop=size+overlap*(i-1)
		content=str_length(substring(aln$V1[2*(a-3)],start+1,stop))-str_count(substring(aln$V1[2*(a-3)],start+1,stop),"-")
		if (start==0)
			start=1
		else 
			start=str_length(substring(aln$V1[2*(a-3)],1,start))-str_count(substring(aln$V1[2*(a-3)],1,start),"-")+1
		
		stop=str_length(substring(aln$V1[2*(a-3)],1,stop))-str_count(substring(aln$V1[2*(a-3)],1,stop),"-")
		if ((stop!=0) && (content!=0)) {
			shape=read.table(args[a],fill=T)
			out=data.frame(shape$V1[2:(stop-start+2)],shape[(start+1):(stop+1),-1])
			write.table(file=name4,shape[1,],quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
			write.table(file=name4,out,quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE,sep="\t")
		}
		else 
			write(file=name4,"NO DATA")
	}
	a=a+1
}
