# B0KstMuMu
Method of Moment:
1. change the directory of reading efficiency
   plugins/Moment.cc l382
2. change the directory of sample
   plugins/Moment.cc l213(gen) l442(reco)
3. compile the program: make Moment
4. run 
   if gen-level: ./Moment Paramater-type gen
   if reco-level: ./Moment Paramater-type reco 
