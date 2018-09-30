#module load R/3.4.0
#R

#batch1-2 v.s. CLIP

group_clip <- 175
group_no_clip <- 907
others_clip <- 568
others_no_clip <- 11549

overlap <- matrix (c(group_clip,group_no_clip,others_clip,others_no_clip),nrow=2, dimnames= list(c("CLIP","No CLIP"),c("DEG","others")))

print(paste("chisq:",chisq.test(overlap)$p.value,sep=""))
print(paste("fisher:",fisher.test(overlap, alternative = "greater")$p.value,sep=""))
print(paste("hypergeometric:",phyper(group_clip-1, (group_clip+others_clip), (group_clip+group_no_clip+others_clip+others_no_clip), (group_clip+group_no_clip), lower.tail=FALSE),sep=""))

##################
        DEG others
CLIP    175    568
No CLIP 907  11549

"chisq:4.05046889900191e-55"
"fisher:3.27841449115987e-40"
"hypergeometric:2.28019249097458e-43"

#batch3-4 v.s. CLIP

group_clip <- 13
group_no_clip <- 125
others_clip <- 730
others_no_clip <- 12376

####################
        DEG others
CLIP     13    730
No CLIP 125  12376

"chisq:0.0768319111552574"
"fisher:0.046180808395607"
"hypergeometric:0.0319634535734533"


