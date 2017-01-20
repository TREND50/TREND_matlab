function [ pos_ant ] = GetAntennaPosition( antID, GRANDproto)
% Retrieve antenna position from its ID
% OMH 04/02/2016

if GRANDproto == 1
   pos_ant = SetupCharacteristics_GRANDproto(antID,6800);
else  %TREND
   pos_ant = SetupCharacteristics_TREND50(antID,3003);
end
  