### Likelihood ratio accounting for IBD probabilities ###########

## p is the allele frequency of the allele A
## q is the allele frequency of the allele B

## AA-AA: k0*p^4 + k1*p^3 + k2*p^2
## AA-AB: k0*2*p^3*q + k1*p^2*q
## AA-BB: k0*p^2*q^2
## AB-AB: k0*4*p^2*q^2 + k1*p*q + k2*2*p*q

# 0 is major homozygote
# 1 is heterozogyte
# 2 is minor homozygote

## sample1: a numeric vector of 0, 1 and 2 for the first individual
## sample2: a numeric vector of 0, 1 and 2 for the second individual
## maf: a numeric vector with the minor allele frequencies
## rel.num: the IBD probabilities of the relationship from the numerator
## rel.den: the IBD probabilities of the relationship from the denominator
## e: the parameter for genotype error


LR_IBD <- function(sample1,sample2,maf,rel.num,rel.den,e){
  
  f1 = function(p,q,k0,k1,k2){return(k0*(p^4) + k1*(p^3) + k2*(p^2))}
  f2 = function(p,q,k0,k1,k2){return(k0*2*(p^3)*q + k1*(p^2)*q)}
  f3 = function(p,q,k0,k1,k2){return(k0*(p^2)*(q^2))}
  f4 = function(p,q,k0,k1,k2){return(k0*4*(p^2)*(q^2) + k1*p*q + k2*2*p*q)}
  
  p = 1-maf
  q = maf
  
  k0 = rel.num[1]
  k1 = rel.num[2]
  k2 = rel.num[3]
  
  lratio_num = NULL
  
  for(i in 1:length(sample1)){
    
    if(sample1[i]==0 & sample2[i]==0){
      
      lratio_num = c(lratio_num,f1(p[i],q[i],k0,k1,k2))}
    
    if((sample1[i]==0 & sample2[i]==1) | (sample1[i]==1 & sample2[i]==0)){
      
      lratio_num = c(lratio_num,f2(p[i],q[i],k0,k1,k2))}
    
    if((sample1[i]==0 & sample2[i]==2) | (sample1[i]==2 & sample2[i]==0)){
      
      likelihood = f3(p[i],q[i],k0,k1,k2)
      
      if(likelihood==0){lratio_num = c(lratio_num,e)}
      
      if(likelihood!=0){lratio_num = c(lratio_num,f3(p[i],q[i],k0,k1,k2))}}
    
    if(sample1[i]==1 & sample2[i]==1){
      
      lratio_num = c(lratio_num,f4(p[i],q[i],k0,k1,k2))}
    
    if(sample1[i]==2 & sample2[i]==2){
      
      lratio_num = c(lratio_num,f1(q[i],p[i],k0,k1,k2))}
    
    if((sample1[i]==1 & sample2[i]==2) | (sample1[i]==2 & sample2[i]==1)){
      
      lratio_num = c(lratio_num,f2(q[i],p[i],k0,k1,k2))}
  }
  
  
  k0 = rel.den[1]
  k1 = rel.den[2]
  k2 = rel.den[3]
  
  lratio_den = NULL
  
  for(i in 1:length(sample1)){
    
    if(sample1[i]==0 & sample2[i]==0){
      
      lratio_den = c(lratio_den,f1(p[i],q[i],k0,k1,k2))}
    
    if((sample1[i]==0 & sample2[i]==1) | (sample1[i]==1 & sample2[i]==0)){
      
      lratio_den = c(lratio_den,f2(p[i],q[i],k0,k1,k2))}
    
    if((sample1[i]==0 & sample2[i]==2) | (sample1[i]==2 & sample2[i]==0)){
      
      likelihood = f3(p[i],q[i],k0,k1,k2)
      
      if(likelihood==0){lratio_den = c(lratio_den,e)}
      
      if(likelihood!=0){lratio_den = c(lratio_den,f3(p[i],q[i],k0,k1,k2))}}
    
    if(sample1[i]==1 & sample2[i]==1){
      
      lratio_den = c(lratio_den,f4(p[i],q[i],k0,k1,k2))}
    
    if(sample1[i]==2 & sample2[i]==2){
      
      lratio_den = c(lratio_den,f1(q[i],p[i],k0,k1,k2))}
    
    if((sample1[i]==1 & sample2[i]==2) | (sample1[i]==2 & sample2[i]==1)){
      
      lratio_den = c(lratio_den,f2(q[i],p[i],k0,k1,k2))}
    
  }
  
  lratio = lratio_num/lratio_den
  
  return(sum(log10(lratio))/length(lratio)) 
}
