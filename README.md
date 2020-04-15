ROOT version:6.12.04

Method of Moment:
1. change the directory of reading efficiency
   plugins/Moment.cc L396
2. change the directory of sample
   plugins/Moment.cc L222(gen) L464(reco)
3. compile the program: make Moment
4. run 
   if gen-level: ./Moment Paramater-type gen
   if reco-level: ./Moment Paramater-type reco event-tag(mis or good)
	parameter-type: FlS,AFBS,P1S,P2S,P3S,P4pS,P5pS,P6pS,P8pS,S3S,S4S,S5S,S6S,S7S,S8S,S9S


