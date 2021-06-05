# crossover design
#1. set the data
Sequence <- c(rep("Left(AB)",7),rep("Right(BA)",7))
Period1 <- c(0.097, 0.07, 0.118, 0.067,0.1,0.018,0.021,0.121,0.003, 0.075,-0.013,0.026,0.087,0.027)
Period2 <- c(0.028,0.054,0.05,0.028,0.039,0.027,0.014,0.093,0.081,0.128,0.058,0.096,0.029,0.028)

Moisture <- data.frame(Sequence, Period1, Period2)
ML <- Moisture[Moisture$Sequence == "Left(AB)",-1]
MR <- Moisture[Moisture$Sequence == "Right(BA)",-1]

#Basic Statistics
mean(ML$Period1);sd(ML$Period1);
mean(ML$Period2);sd(ML$Period2);
mean(MR$Period1);sd(MR$Period1);
mean(MR$Period2);sd(MR$Period2);

Diff_ML <- ML$Period2-ML$Period1
Diff_MR <- MR$Period1-MR$Period2

mean(Diff_ML);sd(Diff_ML);
mean(Diff_MR);sd(Diff_MR);

Overall_ML <- c(ML$Period1,ML$Period2)
Overall_MR <- c(MR$Period1,MR$Period2)

mean(Overall_ML);sd(Overall_ML);
mean(Overall_MR);sd(Overall_MR);

#2. make the effect function for two sample T-test
#Treatment Effect
cross_treat <- function(A, B, sided = "two",conf = 0.05){
  A <- A
  B <- B
  side <- sided
  conf <- conf
  Diff1 <- A[,2]-A[,1]
  Diff2 <- B[,2]-B[,1]
  
  m1 <- mean(Diff1)
  m2 <- mean(Diff2)
  n1 <- dim(A)[1]
  n2 <- dim(B)[1]
  df <- n1+n2-2
  v_t <- sqrt(((n1-1)*var(Diff1) + (n2-1)*var(Diff2))/(df))
  
  
  T <- (m1-m2)/(v_t*sqrt(1/n1+1/n2))
  Treat <- (m1-m2)/2
  
  q_t <- qt(1-(conf/2),df)
  CI.lower <- Treat - (q_t * (v_t*sqrt(1/n1+1/n2)))
  CI.upper <- Treat + (q_t * (v_t*sqrt(1/n1+1/n2)))
  
  p_value = NULL
  if(side == "left"){
    p_value = pt(T,df)
  }else if(side == "right"){
    p_value = 1-pt(T,df)
  }else{
    a <- pt(T,df)
    if(a <= 0.5){
      p_value = 2*(pt(T,df))
    }else{
      p_value = 2*(1-pt(T,df))
    }
    
  }
  
  cat("Treatment effect","\n","Treatment effect : " , Treat , "\n" , "T_stat : " , T , "\n" 
      , "p_value : " , p_value , "\n" , "df : " , df, "\n"
      ,"CI.lower : ", CI.lower, "\n","CI.upper : ", CI.upper)
}

#Carryover Effect
cross_carry <- function(A,B,sided = "two"){
  A <- A
  B <- B
  side <- sided
  a1 <- A[,1]+A[,2]
  a2 <- B[,1]+B[,2]
  
  m1 <- mean(a1)
  m2 <- mean(a2)
  n1 <- dim(A)[1]
  n2 <- dim(B)[1]
  df <- n1+n2-2
  v_c <- sqrt(((n1-1)*var(a1) + (n2-1)*var(a2))/(df))
  
  T <- (m1-m2)/(v_c*sqrt(1/n1+1/n2))
  
  
  p_value = NULL
  if(side == "left"){
    p_value = pt(T,df)
  }else if(side == "right"){
    p_value = 1-pt(T,df)
  }else{
    a <- pt(T,df)
    if(a <= 0.5){
      p_value = 2*(pt(T,df))
    }else{
      p_value = 2*(1-pt(T,df))
    }
    
  }
  
  cat("Carryover effect","\n","T_stat : " , T , "\n" 
      , "p_value : " , p_value , "\n" , "df : " , df)
}

#Period Effect
cross_period <- function(A, B, sided = "two"){
  A <- A
  B <- B
  Diff1 <- A[,2]-A[,1]
  Diff2 <- B[,1]-B[,2]
  
  m1 <- mean(Diff1)
  m2 <- mean(Diff2)
  n1 <- dim(A)[1]
  n2 <- dim(B)[1]
  df <- n1+n2-2
  v_t <- sqrt(((n1-1)*var(Diff1) + (n2-1)*var(-Diff2))/(df))
  
  
  T <- (m1-m2)/(v_t*sqrt(1/n1+1/n2))
  Period <- (m1-m2)/2
  
  p_value = NULL
  if(side == "left"){
    p_value = pt(T,df)
  }else if(side == "right"){
    p_value = 1-pt(T,df)
  }else{
    a <- pt(T,df)
    if(a <= 0.5){
      p_value = 2*(pt(T,df))
    }else{
      p_value = 2*(1-pt(T,df))
    }
    
  }
  
  cat("Period effect","\n","Period effect : " , Period , "\n" , "T_stat : " , T , "\n" 
      , "p_value : " , p_value , "\n" , "df : " , df)
}

#3. get the result for Effect using T test
cross_treat(ML,MR)
cross_carry(ML,MR)
cross_period(ML,MR)

#4. Nonparametric test - Wilcoxon Rank sum Test
result <- sanon(cbind(Period1, Period2) ~ grp(Sequence) , data = Moisture)
result
summary(result)

#treatment * period interaction
contrast(result, C = cbind(1, 1))

