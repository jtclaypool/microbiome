cooccur_filter <- function(RA="relative abundance", co_per=0.5){
  #create logic (0/1) response for presence or absence of count data
  logic_RA=(RA>0)*1

  #check for exceeding co-occurence threshold (default: 50%)
  filter=(colSums(logic_RA)/nrow(logic_RA))>=co_per
  filter_RA=RA[,filter]
  return(filter_RA)
}
